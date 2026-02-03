using Jutopia
using Jutopia.Ferrite
using ComodoFerrite
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo
using DotMaps
#### tested and verified with MATLAB code
dict = Dict()
input = DotMap(dict)
## GLMakie setting 
GLMakie.closeall()

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
nx, ny = 200, 100   # Number of elements along x and y
grid = create_grid(Lx, Ly, nx, ny)  # Generate the grid
F, V = FerriteToComodo(grid)

M = GeometryBasics.Mesh(V, F)
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], aspect=DataAspect(), xlabel="X", ylabel="Y", title="Mesh with Boundary Conditions")
xlims!(ax, -0.5, 2.5)
ylims!(ax, -0.5, 1.5)
poly!(ax, M, color=(Gray(0.95), 0.3), strokecolor=:black, strokewidth=1, shading=true, transparency=false)

faceset1 = get_boundary_points(grid, getfacetset(grid, "traction"))
scatter!(ax, faceset1, color=:blue, markersize=15.0, marker=:circle,label="traction")

nodeset2 = get_boundary_points(grid, getnodeset(grid, "clamped"))
scatter!(ax, nodeset2, color=:red, markersize=8.0, marker=:hexagon,label="Fixed XY")


axislegend(ax, position=:rb, backgroundcolor=(:white, 0.7), framecolor=:gray)
display(GLMakie.Screen(), fig)

input.grid  = grid
input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh)

input.cv, input.fv = create_values()

input.volfrac = 0.5
input.penalty = 3.0

input.rmin = 0.06

input.traction_vector = (0.0, -1.0)
input.facetset = getfacetset(grid, "traction")
input.max_iter = 300
input.tol = 0.01
input.ν= 0.3
input.Emin = 1e-9
input.E0 = 1.0

input.filter_type = :sensitivity_filter
input.model_type = :dim2d
input.load_type=:traction 
ρ_cells = run_optimization(input)


function final_figure(F, V, ρ_cells; plot_type=:elements, colormap = :Spectral, strokewidth=0.0, strokecolor=:black)           
    V = [Point{3,Float64}(v[1], v[2], 0.0) for v in V]
    C_F = ρ_cells[1]

    if plot_type == :elements 
        CF = FaceView(C_F, eltype(F).(eachindex(F)))
        NF = FaceView(facenormal(F,V), eltype(F).(eachindex(F)))
        M = GeometryBasics.Mesh(V, F, normal=NF, color=CF)
    elseif plot_type == :nodes
        CV = simplex2vertexdata(F, C_F,  V; weighting=:size)
    end

    f = Figure(size = (800, 600))
    ax1 = Axis(f[1,1], aspect= DataAspect())        
    if plot_type == :elements
        hm1 = meshplot!(ax1, M; strokewidth=strokewidth, strokecolor=strokecolor, colormap=colormap, colorrange=(0.0, 1.0))
    elseif plot_type == :nodes    
        hm1 = meshplot!(ax1, F, V; strokewidth=strokewidth, strokecolor=strokecolor, color=CV, colormap=colormap, colorrange=(0.0, 1.0))
    end
    Colorbar(f[1,2], hm1)    
    
    hSlider = Slider(f[2,:], range = 1:length(ρ_cells), startvalue = 1, linewidth=30)
   
    on(hSlider.value) do i   
        C_F = ρ_cells[i]
        if plot_type == :elements                          
            CF = FaceView(C_F, eltype(F).(eachindex(F)))
            M = GeometryBasics.Mesh(V, F, normal=NF, color=CF)
            hm1[1] = M
            ax1.title = "Iteration: $i"
        elseif plot_type == :nodes  
            CV = simplex2vertexdata(F, C_F,  V; weighting=:size)
            hm1.color = CV
            ax1.title = "Iteration: $i"
        end
    end 
    # slider2anim(f, hSlider, joinpath(Jutopia_dir(),"assets","temp","cantilever_tri3.gif"); backforth=true, duration=3)
    display(f)
    return f 
end

fig = final_figure(F, V, ρ_cells; plot_type=:elements, colormap = (:turbo), strokewidth=0.0, strokecolor=:black) 