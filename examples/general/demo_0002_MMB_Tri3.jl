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
nx, ny = 30, 30   # Number of elements along x and y
grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid
F, V = FerriteToComodo(grid, Ferrite.Triangle)

M = GeometryBasics.Mesh(V, F)
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], aspect=DataAspect(), xlabel="X", ylabel="Y", title="Mesh with Boundary Conditions")
xlims!(ax, -0.5, 2.5)
ylims!(ax, -0.5, 1.5)
poly!(ax, M, color=(Gray(0.95), 0.3), strokecolor=:black, strokewidth=1, shading=true, transparency=false)

nodeset1 = get_boundary_points(grid, getnodeset(grid, "nodal_force"), Nodes, Ferrite.Triangle)
scatter!(ax, nodeset1, color=:blue, markersize=20.0, marker=:circle, strokecolor=:black, strokewidth=2, label="nodal force")

nodeset2 = get_boundary_points(grid, getnodeset(grid, "support_1"), Nodes, Ferrite.Triangle)
scatter!(ax, nodeset2, color=:red, markersize=8.0, marker=:hexagon, strokecolor=:black, strokewidth=2, label="Fixed X")

nodeset3 = get_boundary_points(grid, getnodeset(grid, "support_2"), Nodes, Ferrite.Triangle)
scatter!(ax, nodeset3, color=:green, markersize=15.0, marker=:diamond, strokecolor=:black, strokewidth=2, label="Fixed Y")


axislegend(ax, position=:rb, backgroundcolor=(:white, 0.7), framecolor=:gray)
display(GLMakie.Screen(), fig)

dh = create_dofhandler(grid)
ch = create_bc(dh)

cv, _ = create_values() # cv: cell_values

volfrac = 0.5
penalty = 3.0

rmin = 0.1 * min(Lx, Ly)

load_vector = (0.0, -1.0)
nodeid = getnodeset(grid, "nodal_force")
ρ = fill(0.5, getncells(grid))
max_iter = 1000
tol = 0.01

# === Run optimization and extract densities ===
top = run_optimization(twoD, Nodal, grid, dh, ch , cv, ρ, penalty, nodeid, load_vector, volfrac, rmin, tol, max_iter)

ρ_cells = top.ρ_cells
ρ_nodes = top.ρ_nodes  # List of nodal density fields for each timestep

nothing