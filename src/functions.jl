struct twoD end
struct threeD end 
struct Traction end
struct Nodal end
######################################################
"""
   get_material_matrix(twoD)
computes and returns the material stiffness matrix for 2D plane stress linear elasticity
"""
function get_material_matrix(::Type{twoD})
    E, ν = 1.0, 0.3
    C_voigt = E * [1.0 ν 0.0; ν 1.0 0.0; 0.0 0.0 (1-ν)/2] / (1 - ν^2)
    return fromvoigt(SymmetricTensor{4,2}, C_voigt)
end
######################################################
"""
    assemble_cell!(twoD, ke, cell_values, ρ, penalty)

computes and returns the local stiffness matrix ke for 2D plane stress linear elasticity
"""
function assemble_cell!(::Type{twoD}, ke, cell_values, ρ, penalty)
    C = get_material_matrix(twoD)
    for qp in 1:getnquadpoints(cell_values)
        dΩ = getdetJdV(cell_values, qp)
        for i in 1:getnbasefunctions(cell_values)
            ∇Ni = shape_gradient(cell_values, qp, i)
            for j in 1:getnbasefunctions(cell_values)
                ∇δNj = shape_symmetric_gradient(cell_values, qp, j)
                ke[i, j] +=  (ρ)^(penalty)* (∇Ni ⊡ C ⊡ ∇δNj) * dΩ
            end
        end
    end
    return ke
end
######################################################
"""
    assemble_global!(::Type{twoD}, K, dh, cell_values, ρ, penalty)

computes and returns the global stiffness matrix K for 2D plane stress linear elasticity
"""
function assemble_global!(::Type{twoD}, K, dh, cell_values, ρ, penalty)
    n_basefuncs = getnbasefunctions(cell_values)
    ke = zeros(n_basefuncs, n_basefuncs)
    assembler = start_assemble(K)

    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cell_values, cell)
        fill!(ke, 0.0)
        local_ρ = ρ[cell_index]
        assemble_cell!(twoD,ke, cell_values, local_ρ, penalty)
        assemble!(assembler, celldofs(cell), ke)
    end
    return K
end
######################################################
"""
    assemble_external_forces!(twoD, f_ext, dh, facetset, facetvalues, prescribed_traction)

computes and returns the global external forces for 2D plane stress linear elasticity (traction)
"""
function assemble_external_forces!(::Type{twoD}, f_ext, dh, facetset, facetvalues, prescribed_traction)
    # Create a temporary array for the facet's local contributions to the external force vector
    fe_ext = zeros(getnbasefunctions(facetvalues))
    for facet in FacetIterator(dh, facetset)
        # Update the facetvalues to the correct facet number
        reinit!(facetvalues, facet)
        # Reset the temporary array for the next facet
        fill!(fe_ext, 0.0)
        # Access the cell's coordinates
        for qp in 1:getnquadpoints(facetvalues)
            # Calculate the global coordinate of the quadrature point.
            # Get the integration weight for the current quadrature point.
            dΓ = getdetJdV(facetvalues, qp)
            for i in 1:getnbasefunctions(facetvalues)
                Nᵢ = shape_value(facetvalues, qp, i)
                fe_ext[i] += prescribed_traction ⋅ Nᵢ * dΓ
            end
        end
        # Add the local contributions to the correct indices in the global external force vector
        assemble!(f_ext, celldofs(facet), fe_ext)
    end
    return f_ext
end
######################################################
"""
    compliance_and_sensitivity(twoD, ρ, dh, u, cell_values, penalty)

computes and returns the compliance and differentiation of compliance with respect to density for 2D plane stress linear elasticity
"""
function compliance_and_sensitivity(::Type{twoD},ρ, dh, u, cell_values, penalty)
    dC_dρ = zeros(length(ρ))
    C = 0.0

    ke = zeros(getnbasefunctions(cell_values), getnbasefunctions(cell_values))

    for (cell_index, (cell)) in enumerate(CellIterator(dh))
        reinit!(cell_values, cell)
        fill!(ke, 0.0)

        ρ_local = ρ[cell_index]
        assemble_cell!(twoD,ke, cell_values, ρ_local, 0.0)  # K₀ only!

        eldofs = celldofs(cell)
        ue = u[eldofs]        
        C += 0.5 * ρ_local^penalty * ue' * ke * ue
        dC_dρ[cell_index] = -penalty * ρ_local^(penalty - 1) * ue' * ke * ue
    end
    return C, dC_dρ
end
######################################################
"""
   vertexdofs(twoD, dh, vertexid)
   nodeid_to_vertexindex(twoD, grid, nodeid)
   apply_nodal_force!(twoD, grid, nodeid, load_vector, f, dh)
computes and returns the global external forces for 2D plane stress linear elasticity (nodal force)
"""
function vertexdofs(::Type{twoD}, dh::DofHandler, vertexid::VertexIndex)
    cellid, lvidx = vertexid
    sdh = dh.subdofhandlers[dh.cell_to_subdofhandler[cellid]]
    local_vertex_dofs = Int[]

    for ifield in 1:length(sdh.field_names)
        offset = Ferrite.field_offset(sdh, ifield)
        field_dim = Ferrite.n_components(sdh, ifield)
        field_ip = isa(sdh.field_interpolations[ifield], Ferrite.VectorizedInterpolation) ?
                   sdh.field_interpolations[ifield].ip :
                   sdh.field_interpolations[ifield]

        vert = Ferrite.vertexdof_indices(field_ip)[lvidx]

        for vdof in vert, d in 1:field_dim
            push!(local_vertex_dofs, (vdof - 1) * field_dim + d + offset)
        end
    end

    dofs = zeros(Int, ndofs_per_cell(dh, cellid))
    celldofs!(dofs, dh, cellid)

    return dofs[local_vertex_dofs]
end
########
function nodeid_to_vertexindex(::Type{twoD}, grid::Grid, nodeid::Int)
    for (cellid, cell) in enumerate(grid.cells)
        for (i, nodeid2) in enumerate(cell.nodes)
            if nodeid == nodeid2
                return VertexIndex(cellid, i)
            end
        end
    end
    error("Node $(nodeid) does not belong to any cell")
end
########
function apply_nodal_force!(::Type{twoD}, grid, nodeid, load_vector, f, dh)
    # Get the coordinates of the nodes
    coords = [Ferrite.get_node_coordinate(grid, id) for id in nodeid]

    # Determine if the edge is vertical (x constant) or horizontal (y constant)
    is_vertical = all(abs(coords[1][1] - coord[1]) < 1e-8 for coord in coords)

    # Set primary direction based on edge orientation
    primary_direction = is_vertical ? 2 : 1  # 2 = y-direction, 1 = x-direction

    # Convert OrderedSet to Vector for indexing
    nodeid_vec = collect(nodeid)

    # Check if there is only one node
    if length(nodeid) == 1
        # Handle the single-node case
        single_node = nodeid[1]

        # Get the vertex index and DOFs for this node
        vertex = nodeid_to_vertexindex(twoD,grid, single_node)
        dofs = vertexdofs(twoD, dh, vertex)

        # Apply the entire load vector directly to this node's DOFs
        f[dofs[1:2]] .+= load_vector
    else
        # Handle the multi-node case
        # Sort nodes based on the primary direction
        sorted_indices = sortperm([coord[primary_direction] for coord in coords])
        sorted_nodeid = nodeid_vec[sorted_indices]

        # Compute segment lengths (dy or dx) along the primary direction
        sorted_coords = [coords[idx] for idx in sorted_indices]
        segment_lengths = diff([coord[primary_direction] for coord in sorted_coords])

        # Calculate the total length of the edge for proportional force distribution
        total_length = sum(segment_lengths)

        # Distribute forces among nodes proportionally to segment lengths
        for i in eachindex(segment_lengths)
            # Get the node IDs for the current segment
            node1 = sorted_nodeid[i]
            node2 = sorted_nodeid[i + 1]

            # Calculate the segment's proportional contribution to the load vector
            segment_force = load_vector * segment_lengths[i] / total_length

            # Get the vertex indices and DOFs for the nodes
            vertex1 = nodeid_to_vertexindex(twoD, grid, node1)
            vertex2 = nodeid_to_vertexindex(twoD, grid, node2)
            dofs1 = vertexdofs(twoD, dh, vertex1)
            dofs2 = vertexdofs(twoD, dh, vertex2)

            # Distribute half of the segment's force to each node
            f[dofs1[1:2]] .+= segment_force / 2
            f[dofs2[1:2]] .+= segment_force / 2
        end
    end
    return f
end
######################################################
"""
   get_centroid(grid, cell_idx)

computes and returns the centroid of an element
"""
function get_centroid(grid, cell_idx)
    node_coordinates = Ferrite.getcoordinates(grid, cell_idx)
    centroid = sum(node_coordinates) / length(node_coordinates)
    return centroid[1], centroid[2]
end

"""
   build_centroid_matrix(grid, n)

computes and returns a matrix for the centroid of an element
"""
function build_centroid_matrix(grid, n)
    centroids = Matrix{Float64}(undef, 2, n)
    for e in 1:n
        x, y = get_centroid(grid, e)
        centroids[1, e] = x
        centroids[2, e] = y
    end
    return centroids
end

"""
    create_optimized_sensitivity_filter(grid, n, r_min)

computes and returns neighbors, weights of an element 
"""
function create_optimized_sensitivity_filter(grid, n, r_min)
    # Build spatial index
    centroids = build_centroid_matrix(grid, n)
    tree = KDTree(centroids)
    
    # Find neighbors and weights using KD-tree
    neighbors = [Int[] for _ in 1:n]
    weights = [Float64[] for _ in 1:n]
    
    # Bulk search with correct indexing
    idxs_list = [Int[] for _ in 1:n]
    for e in 1:n
        idxs_list[e] = NearestNeighbors.inrange(tree, centroids[:, e], r_min)
    end
    
    # Now compute distances and weights
    for e in 1:n
        for i in idxs_list[e]
            dist = norm(centroids[:, e] - centroids[:, i])
            H_ei = max(0, r_min - dist)
            push!(neighbors[e], i)
            push!(weights[e], H_ei)
        end
    end
    
    return neighbors, weights
end
"""
    apply_sensitivity_filter(dc_raw, x, neighbors, weights, n)

computes and returns the sensitivity filter
"""
function apply_sensitivity_filter(dc_raw, x, neighbors, weights, n)
    dc_filtered = zeros(n)
    for e in 1:n
        denom = x[e] * sum(weights[e])
        numer = 0.0
        for k in 1:length(neighbors[e])
            i = neighbors[e][k]  # neighbor element index
            w = weights[e][k]    # weight for this neighbor
            numer += w * x[i] * dc_raw[i]
        end
        dc_filtered[e] = numer / denom
    end
    return dc_filtered
end
######################################################
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

"""
    update_density(x, volfrac, dc, element_volumes)

computes and returns the updated density 
"""
function update_density(x, volfrac, dc, element_volumes)
    l1 = 0.0
    l2 = 100_000.0
    move = 0.2
    
    total_volume = sum(element_volumes)
    # Initialize xnew to avoid undefined variable
    xnew = copy(x)
    
    while (l2 - l1 > 1e-4)
        lmid = 0.5 * (l2 + l1)
        xnew = @. max(0.001, max(x - move, min(1.0, min(x + move, x * sqrt(-dc / (lmid .* element_volumes)))))) ###
        
        # Use weighted volume constraint
        if sum(xnew .* element_volumes) - volfrac * total_volume > 0
            l1 = lmid
        else
            l2 = lmid
        end
    end
    return xnew
end
######################################################
"""
    compute_nodal_data(grid, element_data)

convert element data to node date
"""
function compute_nodal_data(grid, element_data)
    # Get the total number of nodes in the grid.
    nnodes = Ferrite.getnnodes(grid)
    nodal_data = zeros(Float64, nnodes)
    count_elements = zeros(Int, nnodes)
    cells = getcells(grid)
    for (element_id, cell) in enumerate(cells)
        node_ids = cell.nodes
        for node_id in node_ids
            nodal_data[node_id] += element_data[element_id] 
            count_elements[node_id] += 1  
        end
    end
    for node_id in 1:nnodes
        if count_elements[node_id] > 0
            nodal_data[node_id] /= count_elements[node_id]  
            # nodal_data[node_id] = nodal_data[node_id] / count_elements[node_id]
        end
    end
    return nodal_data
end
"""
    run_fem(::Type{twoD}, ::Type{Traction}, input::Dict)

finite element solver for 2D plane stress linear elasticity with Traction
"""
function run_fem(::Type{twoD}, ::Type{Traction}, input::Dict)

    ρ = input["ρ"]
    
    penalty = input["penalty"]
    prescribed_traction = input["prescribed_traction"]
    
    dh = input["dh"]
    ch = input["ch"]
    # # Create CellValues and FacetValues
    cell_values, facet_values = input["cell_values"], input["facet_values"]

    K = allocate_matrix(dh)
    K = assemble_global!(twoD, K, dh, cell_values, ρ, penalty)

    f_ext = zeros(ndofs(dh))
    facetset = input["facetset"]
    f_ext = assemble_external_forces!(twoD, f_ext, dh, facetset, facet_values, prescribed_traction)
    apply!(K, f_ext, ch)

    ## Solve linear system
    u = K \ f_ext

    C, dC_dρ  = compliance_and_sensitivity(twoD, ρ, dh, u, cell_values, penalty)    
    return C, dC_dρ
end
######################################################
"""
    run_fem(twoD, Nodal, input)

finite element solver for 2D plane stress linear elasticity with Nodal force
"""
function run_fem(::Type{twoD}, ::Type{Nodal}, input::Dict)

    ρ = input["ρ"]
    grid = input["grid"]
    penalty = input["penalty"]
    load_vector = input["load_vector"]
    
    dh = input["dh"]
    ch = input["ch"]
    # # Create CellValues and FacetValues
    cell_values = input["cell_values"]

    K = allocate_matrix(dh)
    K = assemble_global!(twoD, K, dh, cell_values, ρ, penalty)

    f_ext = zeros(ndofs(dh))
    nodeid = input["nodeid"]
    f_ext = apply_nodal_force!(twoD, grid, nodeid, load_vector, f_ext, dh)
    apply!(K, f_ext, ch)

    ## Solve linear system
    u = K \ f_ext
    C, dC_dρ  = compliance_and_sensitivity(twoD, ρ, dh, u, cell_values, penalty)    
    return C, dC_dρ
end
######################################################
struct top
    ρ_cells::Vector{Vector{Float64}} 
    ρ_nodes::Vector{Vector{Float64}} 
end

"""
    run_optimization(twoD, Nodal, input)

"""
function run_optimization(::Type{twoD}, ::Type{Nodal}, input::Dict)

    ##### Parameters
    ρ = input["ρ"]
    grid = input["grid"]
    n = getncells(grid)
    rmin = input["rmin"]
    max_iter = input["max_iter"]
    tol = input["tol"]
    volfrac = input["volfrac"]
    dh = input["dh"]
    cell_values = input["cell_values"]
    ####################
    loop = 0
    change = 1.0
    # Precompute filter
    neighbors, weights = create_optimized_sensitivity_filter(grid, n, rmin)

    ρ_cells =  Vector{Vector{Float64}}() 
    ρ_nodes =  Vector{Vector{Float64}}() 

    while change >tol  && loop < max_iter

        ρ_node = compute_nodal_data(grid, ρ)
        push!(ρ_cells, ρ)
        push!(ρ_nodes, ρ_node)
        loop += 1
        C, dC_dρ,  = run_fem(twoD, Nodal, input)
        dc_filtered = apply_sensitivity_filter(dC_dρ, ρ, neighbors, weights, n)
        element_volumes = calculate_all_cell_volumes(grid, dh, cell_values)
        ρnew = update_density(ρ, volfrac, dc_filtered, element_volumes)
        change = maximum(abs.(ρnew - ρ))
        ρ = copy(ρnew)  # Update for next iteration
        println("Iteration $loop, Compliance: $C, Change: $change")
    end 
    return top(ρ_cells, ρ_nodes)
end 
######################################################
function run_optimization(::Type{twoD}, ::Type{Traction}, input::Dict)

    ##### Parameters
    ρ = input["ρ"]
    grid = input["grid"]
    n = getncells(grid)
    rmin = input["rmin"]
    max_iter = input["max_iter"]
    tol = input["tol"]
    volfrac = input["volfrac"]
    dh = input["dh"]
    cell_values = input["cell_values"]
    ####################
    loop = 0
    change = 1.0
    # Precompute filter
    neighbors, weights = create_optimized_sensitivity_filter(grid, n, rmin)

    ρ_cells =  Vector{Vector{Float64}}() 
    ρ_nodes =  Vector{Vector{Float64}}() 

    while change >tol  && loop < max_iter

        ρ_node = compute_nodal_data(grid, ρ)
        push!(ρ_cells, ρ)
        push!(ρ_nodes, ρ_node)
        loop += 1
        C, dC_dρ,  = run_fem(twoD, Traction, input)
        dc_filtered = apply_sensitivity_filter(dC_dρ, ρ, neighbors, weights, n)
        element_volumes = calculate_all_cell_volumes(grid, dh, cell_values)
        ρnew = update_density(ρ, volfrac, dc_filtered, element_volumes)
        change = maximum(abs.(ρnew - ρ))
        ρ = copy(ρnew)  # Update for next iteration
        println("Iteration $loop, Compliance: $C, Change: $change")
    end 
    return top(ρ_cells, ρ_nodes)
end 