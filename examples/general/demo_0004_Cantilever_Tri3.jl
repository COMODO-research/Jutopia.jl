using Jutopia
using Jutopia.Ferrite
using ComodoFerrite
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo

## GLMakie setting 
GLMakie.closeall()

# Create empty dict
input = Dict()

function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly))
    ]
    grid = Ferrite.generate_grid(Ferrite.Triangle, (nx, ny), corners)
    addnodeset!(grid, "clamped", x -> x[1] ≈ 0.0)
    addfacetset!(grid, "traction", x -> x[1] ≈ Lx && norm(x[2] - 0.5) <= 0.05)
    return grid
end
function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{RefTriangle, order}()^dim
    qr = Ferrite.QuadratureRule{RefTriangle}(2)
    qr_face = Ferrite.FacetQuadratureRule{RefTriangle}(1)
    cell_values = Ferrite.CellValues(qr, ip)
    facet_values = Ferrite.FacetValues(qr_face, ip)
    return cell_values, facet_values
end
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{RefTriangle, 1}()^2)
    Ferrite.close!(dh)
    return dh
end
# Function to create Dirichlet boundary conditions
function create_bc(dh)
    ch = Ferrite.ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, getnodeset(dh.grid, "clamped"), (x, t) -> [0.0, 0.0], [1, 2]))
    Ferrite.close!(ch)
    return ch
end
Lx, Ly = 2.0, 1.0  # Plate dimensions
nx, ny = 120, 60   # Number of elements along x and y
grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid
F, V = FerriteToComodo(grid, Ferrite.Triangle)

M = GeometryBasics.Mesh(V, F)
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], aspect=DataAspect(), xlabel="X", ylabel="Y", title="Mesh with Boundary Conditions")
xlims!(ax, -0.5, 2.5)
ylims!(ax, -0.5, 1.5)
poly!(ax, M, color=(Gray(0.95), 0.3), strokecolor=:black, strokewidth=1, shading=true, transparency=false)

faceset1 = get_boundary_points(grid, getfacetset(grid, "traction"), Faces, Ferrite.Triangle)
scatter!(ax, faceset1, color=:blue, markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label="traction")

nodeset2 = get_boundary_points(grid, getnodeset(grid, "clamped"), Nodes, Ferrite.Triangle)
scatter!(ax, nodeset2, color=:red, markersize=8.0, marker=:hexagon, strokecolor=:black, strokewidth=2, label="Fixed XY")


axislegend(ax, position=:rb, backgroundcolor=(:white, 0.7), framecolor=:gray)
display(GLMakie.Screen(), fig)

grid = grid
dh = create_dofhandler(grid)
ch = create_bc(dh)

cv, fv = create_values()

volfrac = 0.5
penalty = 3.0
rmin = 0.1 * min(Lx, Ly)

traction = (0.0, -1.0)
facetset = getfacetset(grid, "traction")
ρ = fill(0.5, getncells(grid))
max_iter = 1000
tol = 0.01

# === Run optimization and extract densities ===
top = run_optimization(twoD, Traction, grid, dh, ch , cv, fv, ρ, penalty,  facetset, traction, volfrac, rmin, tol, max_iter )

ρ_cells = top.ρ_cells
ρ_nodes = top.ρ_nodes  # List of nodal density fields for each timestep

nothing