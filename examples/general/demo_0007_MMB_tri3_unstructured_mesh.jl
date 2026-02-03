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

w = 2.0
h = 1.0
pointSpacing = 0.03
V = rectanglepoints(w, h, pointSpacing; dir=:acw)

VT = (V,)
R = ([1],)
P = (pointSpacing)

F, V, C = regiontrimesh(VT, R, P)

grid = ComodoToFerrite(F, V)

addnodeset!(grid, "support_1", x -> x[1] ≈ -1.0) #fixed in x-direction
addnodeset!(grid, "support_2", x -> x[1] ≈ 1.0 && x[2] ≈ -0.5) # fixed in y direction
addnodeset!(grid, "nodal_force", x -> x[1] ≈ -1.0 && x[2] ≈ 0.5) # nodal_force

M = GeometryBasics.Mesh(V, F)
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], aspect=DataAspect(), xlabel="X", ylabel="Y", title="Mesh with Boundary Conditions")

poly!(ax, M, color=(Gray(0.95), 0.3), strokecolor=:black, strokewidth=1, shading=true, transparency=false)

nodeset1 = get_boundary_points(grid, getnodeset(grid, "nodal_force"))
scatter!(ax, nodeset1, color=:blue, markersize=20.0, marker=:circle,  label="nodal force")

nodeset2 = get_boundary_points(grid, getnodeset(grid, "support_1"))
scatter!(ax, nodeset2, color=:red, markersize= 15.0, marker=:hexagon, label="Fixed X")

nodeset3 = get_boundary_points(grid, getnodeset(grid, "support_2"))
scatter!(ax, nodeset3, color=:green, markersize=15.0, marker=:diamond,  label="Fixed Y")

axislegend(ax, position=:rb, backgroundcolor=(:white, 0.7), framecolor=:gray)
display(GLMakie.Screen(), fig)

# Function to create CellValues and FacetValues
function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{Ferrite.RefTriangle,order}()^dim
    qr = Ferrite.QuadratureRule{Ferrite.RefTriangle}(2)
    qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefTriangle}(1)
    cell_values = Ferrite.CellValues(qr, ip)
    facet_values = Ferrite.FacetValues(qr_face, ip)
    return cell_values, facet_values
end
# Function to create DofHandler
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefTriangle,1}()^2)
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
input.grid = grid
input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh)


input.cv, _ = create_values() # cv: cell_values
input.volfrac = 0.5
input.penalty = 3.0

input.rmin = 0.1
input.nodal_vector =  (0.0, -1.0)
input.nodeset = getnodeset(grid, "nodal_force")
input.max_iter = 300
input.tol = 0.01
input.ν= 0.3
input.Emin = 1e-9
input.E0 = 1.0

input.filter_type = :sensitivity_filter # :sensitivity_filter # :density_filter
input.model_type = :dim2d
input.load_type=:nodal 
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
    # slider2anim(f, hSlider, joinpath(Jutopia_dir(),"assets","temp","tri3_unstructured.gif"); backforth=true, duration=3)
    display(f)
    return f 
end

## colormap = (:turbo)
fig = final_figure(F, V, ρ_cells; plot_type=:nodes, colormap = Reverse(:grays), strokewidth=0.0, strokecolor=:black)