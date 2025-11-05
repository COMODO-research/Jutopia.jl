function get_material_matrix(model_type, E0, ν)
    if model_type == :dim2d
        C_voigt = E0 * [1.0 ν 0.0; ν 1.0 0.0; 0.0 0.0 (1-ν)/2] / (1 - ν^2)
        C_mat = fromvoigt(SymmetricTensor{4,2}, C_voigt)
    elseif model_type == :dim3d
        C_voigt = E0 / ((1 + ν) * (1 - 2 * ν)) * [
            1-ν ν ν 0 0 0;
            ν 1-ν ν 0 0 0;
            ν ν 1-ν 0 0 0;
            0 0 0 (1-2*ν)/2 0 0;
            0 0 0 0 (1-2*ν)/2 0;
            0 0 0 0 0 (1-2*ν)/2
        ]
        C_mat = fromvoigt(SymmetricTensor{4,3}, C_voigt)
    end
    return C_mat
end # end function

function assemble_cell!(model_type, ke, cell_values, ρ, penalty, E0, ν, Emin)
    C = get_material_matrix(model_type, E0, ν)
    for qp in 1:getnquadpoints(cell_values)
        dΩ = getdetJdV(cell_values, qp)
        for i in 1:getnbasefunctions(cell_values)
            ∇Ni = shape_gradient(cell_values, qp, i)
            for j in 1:getnbasefunctions(cell_values)
                ∇δNj = shape_symmetric_gradient(cell_values, qp, j)
                
                coeff = Emin .+ (ρ .^ penalty) .* (E0 .- Emin)
                ke[i, j] += coeff * (∇Ni ⊡ C ⊡ ∇δNj) * dΩ
            end
        end
    end
    return ke
end # end function

function assemble_global!(model_type, K, dh, cell_values, ρ, penalty, E0, ν, Emin)
    n_basefuncs = getnbasefunctions(cell_values)
    ke = zeros(n_basefuncs, n_basefuncs)
    assembler = start_assemble(K)

    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cell_values, cell)
        fill!(ke, 0.0)
        local_ρ = ρ[cell_index]
        assemble_cell!(model_type, ke, cell_values, local_ρ, penalty, E0, ν, Emin)
        assemble!(assembler, celldofs(cell), ke)
    end
    return K
end # end function


function assemble_external_forces!(f_ext, dh, facetset, facetvalues, traction_vector)
    # Create a temporary array for the facet's local contributions to the external force vector
    fe_ext = zeros(getnbasefunctions(facetvalues))
    for facet in FacetIterator(dh, facetset)
        # Update the facetvalues to the correct facet number
        reinit!(facetvalues, facet)
        # Reset the temporary array for the next facet
        fill!(fe_ext, 0.0)
        # Access the cell's coordinates
        for qp in 1:getnquadpoints(facetvalues)
            # Get the integration weight for the current quadrature point.
            dΓ = getdetJdV(facetvalues, qp)
            for i in 1:getnbasefunctions(facetvalues)
                Nᵢ = shape_value(facetvalues, qp, i)
                fe_ext[i] += traction_vector ⋅ Nᵢ * dΓ
            end
        end
        # Add the local contributions to the correct indices in the global external force vector
        assemble!(f_ext, celldofs(facet), fe_ext)
    end
    return f_ext
end # end function

#----------------------------------------------------
# Map node -> DOFs
Ferrite.get_n_copies(::Interpolation) = 1 # Extend internal function to work for non-vectorized interpolations
function node_to_dof!(dofmap::Matrix{Int}, sdh, fieldname)
    if geometric_interpolation(getcelltype(sdh)) != Ferrite.get_base_interpolation(Ferrite.getfieldinterpolation(sdh, fieldname))
        throw(ArgumentError("A node to dof map is only possible for isoparametric elements (same function and geometric interpolation)"))
    end
    #grid = Ferrite.get_grid(dh)
    ndofs_per_node = size(dofmap, 1)
    for cell in CellIterator(sdh)
        for (i, nodeidx) in enumerate(getnodes(cell))
            for j in 1:ndofs_per_node
                dofidx = dof_range(sdh, fieldname)[ndofs_per_node * (i - 1) + j]
                dofmap[j, nodeidx] = celldofs(cell)[dofidx]
            end
        end
    end
    return dofmap
end

#----------------------------------------------------
function node_to_dof(dh, fieldname)
    grid = Ferrite.get_grid(dh)
    ip = Ferrite.getfieldinterpolation(first(dh.subdofhandlers), fieldname);
    dofmap = zeros(Int, Ferrite.get_n_copies(ip)::Int, getnnodes(grid));
    for sdh in dh.subdofhandlers
        node_to_dof!(dofmap, sdh, fieldname)
    end
    return dofmap
end

#----------------------------------------------------

function apply_nodal_force!(f_ext, dh, nodeset, nodal_vector)
    # 1. Get the DOF mapping for displacement components (:u)
    node_dof_map = node_to_dof(dh, :u)
    
    # Check if nodal_vector size matches the number of components
    ncomp = size(node_dof_map, 1)
    if length(nodal_vector) != ncomp
        error("Nodal force vector size ($(length(nodal_vector))) must match DOFs per node ($ncomp).")
    end

    # 2. Loop over every node in the set
    for node_id in nodeset
        # 3. Loop over every component (DOF index) at that node
        for comp in 1:ncomp
            dof = node_dof_map[comp, node_id]

            # 4. Assign the component value
            f_ext[dof] = nodal_vector[comp]
        end
    end
    
    return f_ext
end

#----------------------------------------------------



#u = run_fem( :dim2d, :nodal,  grid, dh, ch, cell_values, ρ, penalty, E0, ν, Emin; nodeset, nodal_vector)
# u = run_fem(:dim2d, :traction, grid, dh, ch, cell_values, ρ, penalty, E0, ν, Emin; facetset ,facetvalues , traction_vector)

    
function run_fem(model_type, load_type, dh, ch, cell_values, ρ,  penalty, E0,ν, Emin;
    nodeset = nothing, nodal_vector = nothing, facetset = nothing, facetvalues = nothing, traction_vector = nothing)
    # --- Assemble stiffness matrix ---
    K = allocate_matrix(dh)
    K = assemble_global!(model_type, K, dh, cell_values, ρ, penalty, E0, ν, Emin)

    # --- Initialize external force vector ---
    f_ext = zeros(ndofs(dh))

    # --- Apply loads based on type ---
    if load_type == :nodal
        @assert nodeset !== nothing "nodeset must be provided for nodal load"
        @assert nodal_vector !== nothing "nodal_vector must be provided for nodal load"
        f_ext = apply_nodal_force!(f_ext, dh, nodeset, nodal_vector)

    elseif load_type == :traction
        @assert facetset !== nothing "facetset must be provided for traction load"
        @assert facetvalues !== nothing "facetvalues must be provided for traction load"
        @assert traction_vector !== nothing "traction_vector must be provided for traction load"
        f_ext = assemble_external_forces!(f_ext, dh, facetset, facetvalues, traction_vector)

    else
        error("Unknown load type: $load_type. Use :nodal or :traction.")
    end

    # --- Apply boundary conditions ---
    apply!(K, f_ext, ch)

    # --- Solve linear system ---
    u = K \ f_ext

    return u
end


"""
   calculate_cell_volume(cv)

computes and returns the local element volume (single element)
"""
function calculate_cell_volume(cv)
    cell_volume = 0.0
    for q_point in 1:getnquadpoints(cv)
        dΩ = getdetJdV(cv, q_point)
        cell_volume += dΩ
    end
    return cell_volume
end
"""
   calculate_all_cell_volumes(grid, dh, cell_values)

computes and returns the global element volume(all element)
"""
function calculate_all_cell_volumes(grid, dh, cell_values)
    volumes = zeros(getncells(grid))
    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cell_values, cell)
        volumes[cell_index] = calculate_cell_volume(cell_values)
    end
    return volumes
end
