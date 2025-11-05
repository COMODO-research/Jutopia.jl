module Jutopia

using Ferrite
using NearestNeighbors
using Printf
using FIGlet
using DotMaps
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

export get_centroid
export build_centroid_matrix
export create_optimized_sensitivity_filter
export sensitivity_filter
export rho_filter
export sensitivity_filter_chainrule
include("filter.jl")

#------------------------------------------------

export optimality_criteria
include("optimality_criteria.jl")

export run_optimization
include("optim_solver.jl")
end # module Jutopia
