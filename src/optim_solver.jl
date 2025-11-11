function run_optimization(input::DotMap)
    #-------------------------------
    # input parameters
    grid = input.grid
    dh = input.dh
    ch = input.ch
    cv = input.cv
    penalty = input.penalty
    E0 = input.E0
    ν = input.ν
    Emin = input.Emin
    volfrac = input.volfrac
    rmin = input.rmin
    tol = input.tol
    max_iter = input.max_iter

    model_type = input.model_type
    load_type = input.load_type

    filter_type = input.filter_type
    # --- Optional load data ---
    nodeset = haskey(input, "nodeset") ? input.nodeset : nothing
    nodal_vector = haskey(input, "nodal_vector") ? input.nodal_vector : nothing
    facetset = haskey(input, "facetset") ? input.facetset : nothing
    facetvalues = haskey(input, :fv) ? input.fv : nothing
    traction_vector = haskey(input, "traction_vector") ? input.traction_vector : nothing


    n = getncells(grid)
    #neighbors, weights = create_optimized_sensitivity_filter(grid, n, rmin)
    centroid_data = get_all_element_centroids(grid, dh)
    neighbors, weights = create_neighbor_data(centroid_data, rmin)
    element_volumes = calculate_all_cell_volumes(grid, dh, cv)
    
    
    ####################
    loop = 0
    change = 1.0

   
    xPhys_cells = Vector{Vector{Float64}}()

    x = fill(volfrac, getncells(grid))
    xPhys = x
    loop = 0
    change = 1
    # Start with the initial density from input

    while change > tol && loop < max_iter

        push!(xPhys_cells, xPhys)
        loop += 1
        ## fem 
        u = run_fem(
            model_type, load_type, dh, ch, cv, xPhys, penalty, E0, ν, Emin;
            nodeset=nodeset,
            nodal_vector=nodal_vector,
            facetset=facetset,
            facetvalues=facetvalues,
            traction_vector=traction_vector
        )

        ## compliance and sensitivity
        c, dc = compliance_and_sensitivity(model_type, xPhys, dh, u, cv, penalty, E0, ν, Emin)
        
        dv = calculate_all_cell_volumes(grid, dh, cv)
        
        ## FILTERING/MODIFICATION OF SENSITIVITIES
        if filter_type == :sensitivity_filter
            #dc = sensitivity_filter(dc, x, neighbors, weights, n)
            dc = sensitivity_filter(n, x, dc,element_volumes, neighbors,weights)
        elseif filter_type == :density_filter
            #dc, dv = sensitivity_filter_chainrule(dc, dv, neighbors, weights, n)
            dc = calculate_all_filtered_compliance_sensitivities(n, element_volumes, dc ,neighbors,weights)
            dv = calculate_all_volume_constraint_sensitivities(n, element_volumes,neighbors,weights)
        else
            error("the filter type $filter_type is incorrect")
        end

        # xnew, xPhys = optimality_criteria(x, volfrac, dc, dv, element_volumes, filter_type, neighbors, weights, n)
        xnew, xPhys = optimality_criteria(x, volfrac, dc, dv, element_volumes, filter_type, neighbors, weights, n)

        change = maximum(abs.(xnew - x))

        x = xnew  # Update local variable for next iteration

        println("Iteration $loop, Compliance: $c, Change: $change")
        println("----------------------------------------")
    end
    FIGlet.render("This Is Jutopia", "Standard")
    FIGlet.render("Optimization Successful", "Standard")

    return xPhys_cells
end