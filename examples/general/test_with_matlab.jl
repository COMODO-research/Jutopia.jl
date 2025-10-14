using Jutopia
using Jutopia.Ferrite
using ComodoFerrite
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo

input = Dict()
function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly))
    ]
    grid = Ferrite.generate_grid(Ferrite.Triangle, (nx, ny), corners)
    addnodeset!(grid, "support_1", x -> x[1] ≈ 0.0) #fixed in x-direction
    addnodeset!(grid, "support_2", x -> x[1] ≈ Lx && x[2] ≈ 0.0) # fixed in y direction
    addnodeset!(grid, "nodal_force", x -> x[1] ≈ 0.0 && x[2] ≈ Ly) # nodal_force
    return grid
end
# Function to create CellValues and FacetValues
function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{RefTriangle,order}()^dim
    qr = Ferrite.QuadratureRule{RefTriangle}(2)
    qr_face = Ferrite.FacetQuadratureRule{RefTriangle}(1)
    cell_values = Ferrite.CellValues(qr, ip)
    facet_values = Ferrite.FacetValues(qr_face, ip)
    return cell_values, facet_values
end
# Function to create DofHandler
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{RefTriangle,1}()^2)
    Ferrite.close!(dh)
    return dh
end
# Function to create Dirichlet boundary conditions
function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getnodeset(dh.grid, "support_1"), (x, t) -> [0.0], [1]))
    add!(ch, Dirichlet(:u, getnodeset(dh.grid, "support_2"), (x, t) -> [0.0], [2]))
    Ferrite.close!(ch)
    return ch
end
####### 
Lx, Ly = 2.0, 1.0  # Plate dimensions
nx, ny = 40, 40   # Number of elements along x and y
grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid

input["grid"] = grid
input["dh"] = create_dofhandler(grid)
input["ch"] = create_bc(input["dh"])

input["cell_values"], _ = create_values()

input["volfrac"] = 0.5
input["penalty"] = 3.0


input["rmin"] = 0.1 * min(Lx, Ly)

input["load_vector"] = (0.0, -1.0)
input["nodeid"] = getnodeset(grid, "nodal_force")
input["ρ"] = fill(0.5, getncells(grid))
input["max_iter"] = 1000
input["tol"] = 0.01

C, dC_dρ = run_fem(twoD,Nodal, input)

@info "sum dc " sum(dC_dρ)

n = getncells(grid)
rmin = input["rmin"]
neighbors, weights = create_optimized_sensitivity_filter(grid, n, rmin)
dc_filtered = apply_sensitivity_filter(dC_dρ, input["ρ"], neighbors, weights, n)

@info "sum dc filter" sum(dc_filtered)

element_volumes = calculate_all_cell_volumes(grid, input["dh"], input["cell_values"])
ρnew = update_density(input["ρ"], input["volfrac"], dc_filtered, element_volumes)

@info "sum ρnew" sum(ρnew)