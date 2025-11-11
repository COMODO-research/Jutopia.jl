module Jutopia

using Ferrite
using NearestNeighbors
using Printf
using FIGlet
using DotMaps
using GeoInterface
using GeometryOps

function Jutopia_dir()
    joinpath(@__DIR__, "..")
end

export Jutopia_dir
#------------------------------------------------
#------stiffness matrix
export get_material_matrix
export assemble_cell!, assemble_global!
#------load traction
export assemble_external_forces!
#------nodal load
export node_to_dof!
export node_to_dof
export apply_nodal_force!
#------fem solver
export run_fem
#------check volume
export calculate_cell_volume
export calculate_all_cell_volumes
include("fem_functions.jl")
#------------------------------------------------
export compliance_and_sensitivity
include("sensitivity.jl")
#------------------------------------------------
export calculate_raw_weight
export get_all_element_centroids
export calculate_all_cell_volumes
export centroids_to_matrix
export create_neighbor_data
export sensitivity_filter
export density_filter
export calculate_all_volume_constraint_sensitivities
export calculate_all_filtered_compliance_sensitivities

include("filter.jl")

#------------------------------------------------

export optimality_criteria
include("optimality_criteria.jl")

export run_optimization
include("optim_solver.jl")
end # module Jutopia
