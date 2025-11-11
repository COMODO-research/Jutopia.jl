using Jutopia
using Jutopia.Ferrite
using ComodoFerrite
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo
using DotMaps
using FerriteGmsh
dict = Dict()
input = DotMap(dict)

fileName_mesh = joinpath(Jutopia_dir(),"assets","fine_rect.msh")
grid = togrid(fileName_mesh)
addnodeset!(grid, "support_1", x -> x[1] ≈ 0.0) #fixed in x-direction
addnodeset!(grid, "support_2", x -> x[1] ≈ 2.0 && x[2] ≈ 0.0) # fixed in y direction
addnodeset!(grid, "nodal_force", x -> x[1] ≈ 0.0 && x[2] ≈ 1.0) # nodal_force
F, V = FerriteToComodo(grid, Ferrite.Triangle)

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

input.filter_type = :sensitivity_filter # :sensitivity_filter
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
    display(f)
    return f 
end

## colormap = (:turbo)
f_nodes = final_figure(F, V, ρ_cells; plot_type=:nodes, colormap = Reverse(:grays), strokewidth=0.0, strokecolor=:black)