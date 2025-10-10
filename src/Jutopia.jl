module Jutopia
using Ferrite
using NearestNeighbors
using Printf

export twoD, Traction, Nodal, threeD
############################################################
###### fem
export get_material_matrix
export assemble_cell!
export assemble_global!
export assemble_external_forces!
export compliance_and_sensitivity
export vertexdofs , nodeid_to_vertexindex, apply_nodal_force!

########## top 
export get_centroid, build_centroid_matrix, create_optimized_sensitivity_filter, apply_sensitivity_filter
export calculate_cell_volume, calculate_all_cell_volumes, update_density

#### solver
export run_fem
export run_optimization

include("functions.jl")
end # module Jutopia
