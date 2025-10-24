using Jutopia
using Jutopia.Ferrite
using ComodoFerrite
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo

## GLMakie setting 
GLMakie.closeall()
function create_grid(Lx, Ly, nx, ny)
    corners = [
        Ferrite.Vec{2}((0.0, 0.0)), Ferrite.Vec{2}((Lx, 0.0)),
        Ferrite.Vec{2}((Lx, Ly)), Ferrite.Vec{2}((0.0, Ly))
    ]
    grid = Ferrite.generate_grid(Ferrite.Quadrilateral, (nx, ny), corners)
    addnodeset!(grid, "clamped", x -> x[1] ≈ 0.0)
    addfacetset!(grid, "traction", x -> x[1] ≈ Lx && norm(x[2] - 0.5) <= 0.05)
    return grid
end
function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{Ferrite.RefQuadrilateral, order}()^dim
    qr = Ferrite.QuadratureRule{Ferrite.RefQuadrilateral}(2)
    qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefQuadrilateral}(1)
    cell_values = Ferrite.CellValues(qr, ip)
    facet_values = Ferrite.FacetValues(qr_face, ip)
    return cell_values, facet_values
end
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()^2)
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
Lx, Ly = 3.0, 1.0  # Plate dimensions
nx, ny = 200, 100   # Number of elements along x and y
grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid
F, V = FerriteToComodo(grid, Ferrite.Quadrilateral)

M = GeometryBasics.Mesh(V, F)
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], aspect=DataAspect(), xlabel="X", ylabel="Y", title="Mesh with Boundary Conditions")
xlims!(ax, -0.5, 2.5)
ylims!(ax, -0.5, 1.5)
poly!(ax, M, color=(Gray(0.95), 0.3), strokecolor=:black, strokewidth=1, shading=true, transparency=false)

faceset1 = get_boundary_points(grid, getfacetset(grid, "traction"), Faces, Ferrite.Quadrilateral)
scatter!(ax, faceset1, color=:blue, markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label="traction")

nodeset2 = get_boundary_points(grid, getnodeset(grid, "clamped"), Nodes, Ferrite.Quadrilateral)
scatter!(ax, nodeset2, color=:red, markersize=8.0, marker=:hexagon, strokecolor=:black, strokewidth=2, label="Fixed XY")


axislegend(ax, position=:rb, backgroundcolor=(:white, 0.7), framecolor=:gray)
display(GLMakie.Screen(), fig)

grid = grid
dh = create_dofhandler(grid)
ch = create_bc(dh)

cv, fv = create_values()

volfrac = 0.5
penalty = 3.0
rmin = 0.05 * min(Lx, Ly)

traction = (0.0, -1.0)
facetset = getfacetset(grid, "traction")
ρ = fill(0.5, getncells(grid))
max_iter = 1000
tol = 0.01

# === Run optimization and extract densities ===
top = run_optimization(twoD, Traction, grid, dh, ch , cv, fv, ρ, penalty,  facetset, traction, volfrac, rmin, tol, max_iter )

ρ_cells = top.ρ_cells
ρ_nodes = top.ρ_nodes 

function final_figure(data, loop, nx, ny, Lx, Ly)
    f = Figure(size = (600, 600))
    ax = Axis(f[1,1], aspect = Lx/Ly)
    

    x = range(0, Lx, length = nx + 1)
    y = range(Ly, 0, length = ny + 1)  
    
    data_obs = Observable(reverse(reshape(data[1], nx, ny); dims=2))
    hm = heatmap!(ax, x, y, data_obs, colormap=Reverse(:grays), colorrange=(0.0,1.0))
    Colorbar(f[1,2], hm)
    hSlider = Slider(f[2,1], range = 1:loop, startvalue = 1, linewidth=30)
   
    on(hSlider.value) do i 
        data_obs[] = reverse(reshape(data[i], nx, ny); dims=2)
        ax.title = "Step: $i"
    end 
    display(f)
    return f 
end

loop = length(ρ_cells)
f = final_figure(ρ_cells, loop, nx, ny, Lx, Ly)

# For node data
# loop_nodes = length(ρ_nodes)
# nx_nodes = nx + 1  
# ny_nodes = ny + 1

# f_nodes = final_figure(ρ_nodes, loop_nodes, nx_nodes, ny_nodes, Lx, Ly)


