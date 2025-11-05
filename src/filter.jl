"""
   get_centroid(grid, cell_idx)

computes and returns the centroid of an element
"""
# function get_centroid(grid, cell_idx)
#     node_coordinates = Ferrite.getcoordinates(grid, cell_idx)
#     centroid = sum(node_coordinates) / length(node_coordinates)
#     return centroid[1], centroid[2]
# end
function get_centroid(grid, cell_idx)
    # Get the coordinates of all nodes belonging to the cell
    node_coordinates = Ferrite.getcoordinates(grid, cell_idx)
    
    # Calculate the centroid by summing all node coordinates and dividing by the number of nodes.
    # The result will be a Vec{2} for 2D or Vec{3} for 3D.
    centroid = sum(node_coordinates) / length(node_coordinates)
    
    # Return the full vector/tuple of coordinates, not just the first two.
    return centroid
end
"""
   build_centroid_matrix(grid, n)

computes and returns a matrix for the centroid of an element
"""
# function build_centroid_matrix(grid, n)
#     centroids = Matrix{Float64}(undef, 2, n)
#     for e in 1:n
#         x, y = get_centroid(grid, e)
#         centroids[1, e] = x
#         centroids[2, e] = y
#     end
#     return centroids
# end
function build_centroid_matrix(grid, n)
    # Get the dimension of the grid from the first cell's centroid
    # This assumes all cells are the same dimension, which is standard.
    first_centroid = get_centroid(grid, 1)
    dim = length(first_centroid) # dim will be 2 or 3
    
    # Initialize the matrix with the correct number of rows (dim)
    centroids = Matrix{Float64}(undef, dim, n)
    
    for e in 1:n
        centroid = get_centroid(grid, e)
        # Assign all components of the centroid vector to the column of the matrix
        centroids[:, e] = centroid 
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
# function sensitivity_filter(dc_raw, x, neighbors, weights, n)
#     dc_filtered = zeros(n)
#     for e in 1:n
#         denom = x[e] * sum(weights[e])
#         numer = 0.0
#         for k in 1:length(neighbors[e])
#             i = neighbors[e][k]  # neighbor element index
#             w = weights[e][k]    # weight for this neighbor
#             numer += w * x[i] * dc_raw[i]
#         end
#         dc_filtered[e] = numer / (denom)
#     end
#     return dc_filtered
# end

function sensitivity_filter(dc_raw, x, neighbors, weights, n)
    dc_filtered = zeros(n)
    for e in 1:n
        denom = max(1e-3, x[e]) * sum(weights[e])
        numer = 0.0
        for k in 1:length(neighbors[e])
            i = neighbors[e][k]
            w = weights[e][k]
            numer += w * x[i] * dc_raw[i]
        end
        dc_filtered[e] = numer / denom
    end
    return dc_filtered
end

###############################################################
function rho_filter(x, neighbors, weights, n)
    rho_filtered = zeros(n)
    for e in 1:n
        denom = sum(weights[e])   
        numer = 0.0  
        for k in 1:length(neighbors[e])
            i = neighbors[e][k]
            w = weights[e][k]
            
            numer += w * x[i]
        end    
        rho_filtered[e] = numer / denom
    end
    return rho_filtered
end
###############################################################
function sensitivity_filter_chainrule(dc_raw, dV_raw, neighbors, weights, n)
    dc_filtered = zeros(n)
    dV_filtered = zeros(n)
    for j in 1:n
        numer_c = 0.0
        numer_v = 0.0
        denom = 0.0
        for k in 1:length(neighbors[j])
            e = neighbors[j][k]
            H_je = weights[j][k]

            numer_c += H_je * dc_raw[e]
            numer_v += H_je * dV_raw[e]
            denom += H_je
        end
        if denom > 0
            dc_filtered[j] = numer_c / denom
            dV_filtered[j] = numer_v / denom
        else
            dc_filtered[j] = 0.0
            dV_filtered[j] = 0.0
        end
    end
    return dc_filtered, dV_filtered
end
