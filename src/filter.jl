function calculate_raw_weight(distance_ij::Float64, R::Float64)
    return max(0.0, R - distance_ij)
end

function get_all_element_centroids(grid::Ferrite.AbstractGrid, dh::Ferrite.DofHandler)
    centroid_data = []
    for cell_idx in 1:getncells(grid)
        ferrite_coords = Ferrite.getcoordinates(grid, cell_idx)
        dim = length(ferrite_coords[1])  # 2 for 2D, 3 for 3D
        # Compute average for each dimension
        centroid = tuple([mean(c[d] for c in ferrite_coords) for d in 1:dim]...)
        push!(centroid_data, centroid)
    end
    return centroid_data
end

function centroids_to_matrix(centroid_data::Vector, dim::Int, n::Int)
    matrix = Matrix{Float64}(undef, dim, n)
    for i in 1:n
        # Assuming GeoInterface/GeometryOps results are indexable like a vector
        matrix[:, i] = [centroid_data[i][k] for k in 1:dim]
    end
    return matrix
end


function calculate_all_cell_volumes(dh::Ferrite.DofHandler, cv::Ferrite.CellValues)
    element_volumes = Vector{Float64}(undef, getncells(dh.grid))
    for (cell_idx, cell) in enumerate(Ferrite.CellIterator(dh))
        Ferrite.reinit!(cv, cell)
        cell_volume = 0.0
        for q_point in 1:Ferrite.getnquadpoints(cv)
            dΩ = Ferrite.getdetJdV(cv, q_point)
            cell_volume += dΩ
        end
        element_volumes[cell_idx] = cell_volume
    end
    return element_volumes
end


function create_neighbor_data(centroid_data::Vector, R::Float64)
    N = length(centroid_data)
    dim = length(centroid_data[1])
    
    centroids_matrix = centroids_to_matrix(centroid_data, dim, N)
    tree = KDTree(centroids_matrix)
    neighbor_idxs_list = NearestNeighbors.inrange(tree, centroids_matrix, R)
    
    neighbors = [Int[] for _ in 1:N]
    weights = [Float64[] for _ in 1:N]
    
    for e in 1:N
        centroid_e = centroid_data[e]
        for i in neighbor_idxs_list[e]
            centroid_i = centroid_data[i]
            dist = GeometryOps.distance(centroid_e, centroid_i)
            H_ei = calculate_raw_weight(dist, R) 
            
            push!(neighbors[e], i)
            push!(weights[e], H_ei)
        end
    end
    return neighbors, weights
end

# --- 1. THE SENSITIVITY FILTER FUNCTION (Optimized) ---

function sensitivity_filter(
    n::Int,
    x::Vector{Float64},
    dc_raw::Vector{Float64},
    v::Vector{Float64},
    neighbors::Vector{Vector{Int}},
    weights::Vector{Vector{Float64}}
)
    EPSILON = 1e-3 
    dc_filtered = Vector{Float64}(undef, n)

    for e in 1:n
        rho_e = x[e]
        v_e = v[e]

        numerator_sum = 0.0
        denominator_sum_w = 0.0

        for k in 1:length(neighbors[e])
            i = neighbors[e][k]
            w_ei = weights[e][k]
            
            rho_i = x[i]
            raw_sensitivity_i = dc_raw[i]
            v_i = v[i]

            numerator_sum += w_ei * rho_i * (raw_sensitivity_i / v_i)
            denominator_sum_w += w_ei
        end

        safe_rho_e = max(rho_e, EPSILON)
        denominator = (safe_rho_e / v_e) * denominator_sum_w

        dc_filtered[e] = denominator == 0.0 ? 0.0 : numerator_sum / denominator
    end

    return dc_filtered
end

function density_filter(
    n::Int,
    x::Vector{Float64},
    v::Vector{Float64},
    neighbors::Vector{Vector{Int}},
    weights::Vector{Vector{Float64}}
)
    filtered_densities = Vector{Float64}(undef, n)

    for e in 1:n
        
        numerator_sum = 0.0
        denominator_sum = 0.0

        for k in 1:length(neighbors[e])
            i = neighbors[e][k]
            w_ei = weights[e][k]
            
            rho_i = x[i]
            v_i = v[i]
            
            weighted_volume = w_ei * v_i
            
            numerator_sum += weighted_volume * rho_i
            
            denominator_sum += weighted_volume
        end
        
        if denominator_sum == 0.0
            filtered_densities[e] = x[e]
        else
            filtered_densities[e] = numerator_sum / denominator_sum
        end
    end

    return filtered_densities
end

# Use precomputed neighbors::Vector{Vector{Int}} and weights::Vector{Vector{Float64}}
function calculate_all_volume_constraint_sensitivities(
    n::Int,
    element_volumes::Vector{Float64},
    neighbors::Vector{Vector{Int}},
    weights::Vector{Vector{Float64}}
)
    sensitivities = zeros(Float64, n)

    for e in 1:n
        v_e = element_volumes[e]
        total = 0.0

        for (k, i) in enumerate(neighbors[e])
            H_ei = weights[e][k]
            v_i = element_volumes[i]

            denom = 0.0
            for (m, j) in enumerate(neighbors[i])
                H_ij = weights[i][m]
                denom += H_ij * element_volumes[j]
            end

            if denom > 0.0
                total += v_i * (H_ei * v_e) / denom
            end
        end

        sensitivities[e] = total
    end

    return sensitivities
end

function calculate_all_filtered_compliance_sensitivities(
    n::Int,
    element_volumes::Vector{Float64},
    raw_compliance_sensitivities::Vector{Float64},
    neighbors::Vector{Vector{Int}},
    weights::Vector{Vector{Float64}}
)
    sensitivities = zeros(Float64, n)

    for e in 1:n
        v_e = element_volumes[e]
        total = 0.0

        for (k, i) in enumerate(neighbors[e])
            H_ei = weights[e][k]
            denom = 0.0

            for (m, j) in enumerate(neighbors[i])
                H_ij = weights[i][m]
                denom += H_ij * element_volumes[j]
            end

            if denom > 0.0
                raw_i = raw_compliance_sensitivities[i]
                total += raw_i * (H_ei * v_e) / denom
            end
        end

        sensitivities[e] = total
    end

    return sensitivities
end
