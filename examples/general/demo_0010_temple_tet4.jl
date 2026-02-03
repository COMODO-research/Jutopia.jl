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

# Function to check if a point is inside a rectangle on a given plane
function in_rectangle(x, cx, cy, cz, width, depth)
    abs(x[1]-cx) <= width/2 && abs(x[2]-cy) <= depth/2 && abs(x[3]-cz) <= 1e-8
end

# Function to create a 3D grid
function create_grid(Lx, Ly, Lz, nx, ny, nz)
    left = Ferrite.Vec(0.0, 0.0, 0.0)
    right = Ferrite.Vec(Lx, Ly, Lz)
    generate_grid(Ferrite.Tetrahedron, (nx, ny, nz), left, right)
end

function create_boundary(grid,Lx,Ly,Lz)
    # Top rectangle parameters
    rect_width, rect_depth = 0.2, 0.2  # be careful about mesh size to capture the boundary
    cx_top, cy_top = Lx/2, Ly/2
    addfacetset!(grid, "top_rectangle", x -> in_rectangle(x, cx_top, cy_top, Lz, rect_width, rect_depth))

    # Bottom small lines at center of each edge
    edge_centers = [
        (Lx/2, 0.0, :y),
        (Lx/2, Ly, :y),
        (0.0, Ly/2, :x),
        (Lx, Ly/2, :x)
    ]
    line_length = 0.1   # half-length along the line   # be careful about mesh size to capture the boundary
    tol = 0.05          # thickness perpendicular to line   # be careful about mesh size to capture the boundary

    for (i,(cx,cy,dir)) in enumerate(edge_centers)
        if dir == :x
            addnodeset!(grid, "bottom_edge_center_$i", x -> (x[3] ≈ 0.0) && abs(x[1]-cx) <= line_length/2 && abs(x[2]-cy) <= tol)
        else
            addnodeset!(grid, "bottom_edge_center_$i", x -> (x[3] ≈ 0.0) && abs(x[1]-cx) <= tol && abs(x[2]-cy) <= line_length/2)
        end
    end
end

# Function to create cell and facet values
function create_values()
    order = 1
    dim = 3
    ip = Lagrange{RefTetrahedron, order}()^dim
    qr = QuadratureRule{RefTetrahedron}(2)
    qr_face = FacetQuadratureRule{RefTetrahedron}(1)
    cell_values = CellValues(qr, ip)
    facet_values = FacetValues(qr_face, ip)
    return cell_values, facet_values
end

# Function to create a DOF handler
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefTetrahedron, 1}()^3)
    Ferrite.close!(dh)
    return dh
end

# Function to create boundary conditions
function create_bc(dh, grid)
    ch = Ferrite.ConstraintHandler(dh)
    for i in 1:4
        dbc = Dirichlet(:u, getnodeset(grid, "bottom_edge_center_$i"), (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3])
        add!(ch, dbc)
    end
    Ferrite.close!(ch)
    return ch
end

Lx, Ly, Lz = 1.0, 1.0, 1.0
nx, ny, nz = 15, 15, 15

grid = create_grid(Lx, Ly, Lz, nx, ny, nz)
create_boundary(grid, Lx, Ly, Lz)
# Convert to Comodo mesh
E, V, F, Fb, CFb_type = FerriteToComodo(grid)
M = GeometryBasics.Mesh(V, F, normal=face_normals(V, F))
Mb = GeometryBasics.Mesh(V, Fb, normal=face_normals(V, Fb))

# Plot
fig_mesh = Figure(size=(1200, 800))

ax1 = AxisGeom(fig_mesh[1, 1], title="Tet4 mesh")
meshplot!(ax1, M, color=:gray, strokecolor=:black, strokewidth=3.0, shading=false, transparency=false)

ax2 = AxisGeom(fig_mesh[1, 2], title="Boundary condition")
meshplot!(ax2, Mb, color=(Gray(0.95), 0.3), strokecolor=:black, strokewidth=2.0, shading=true, transparency=true)

# Top rectangle nodes
facesset_top = get_boundary_points(grid, getfacetset(grid, "top_rectangle"))
scatter!(ax2, facesset_top, color=:red, markersize=15.0, marker=:circle,label="Top rectangle")

# Bottom edge-center nodes
for i in 1:4
    bottom_nodes = get_boundary_points(grid, getnodeset(grid, "bottom_edge_center_$i"))
    scatter!(ax2, bottom_nodes, color=:blue, markersize=15.0, marker=:circle,label="Bottom edge center $i")
end

axislegend(ax2, position=:rb, backgroundcolor=(:white, 0.7), framecolor=:gray)
display(fig_mesh)

input.grid = grid
input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh, grid)

input.cv, input.fv = create_values()

input.volfrac = 0.5
input.penalty = 3.0

# Calculate element size and rmin
element_size_x = Lx / nx
element_size_y = Ly / ny 
element_size_z = Lz / nz
avg_element_size = (element_size_x + element_size_y + element_size_z) / 3

# Choose rmin based on element size (typical range: 1.5-3.0 times element size)
input.rmin = 2.0 * avg_element_size


input.traction_vector = (0.0, 0.0, 1.0)
input.facetset = getfacetset(grid, "top_rectangle")
input.max_iter = 200
input.tol = 0.01
input.ν = 0.3
input.Emin = 1e-9
input.E0 = 1.0

input.filter_type = :sensitivity_filter
input.model_type = :dim3d
input.load_type = :traction
ρ_cells = run_optimization(input)

# VTKGridFile("/Users/aminalibakhshi/Desktop/vtu_results/amin.vtu", input.dh) do vtk
#     write_cell_data(vtk, ρ_cells[end], "density")
#     Ferrite.write_cellset(vtk, grid)
# end


function final_figure(E, V, ρ_cells, plotThreshold = 0.1; plot_type=:elements, colormap = :Spectral, strokewidth=0.0, strokecolor=:black)           
    C_E = ρ_cells[length(ρ_cells)]
    
    B_E = C_E .> plotThreshold
    F = element2faces(E[B_E])
    C_F = repeat(C_E[B_E],inner=4)

    if plot_type == :elements 
        CF = FaceView(C_F, eltype(F).(eachindex(F)))
        NF = FaceView(facenormal(F,V), eltype(F).(eachindex(F)))
        M = GeometryBasics.Mesh(V, F, normal=NF, color=CF)
    elseif plot_type == :nodes
        CV = simplex2vertexdata(F, C_F,  V; weighting=:size)
    end

    f = Figure(size = (800, 600))
    ax1 = AxisGeom(f[1,1])        
    if plot_type == :elements
        hm1 = meshplot!(ax1, M; strokewidth=strokewidth, strokecolor=strokecolor, colormap=colormap, colorrange=(0.0, 1.0))
    elseif plot_type == :nodes    
        hm1 = meshplot!(ax1, F, V; strokewidth=strokewidth, strokecolor=strokecolor, color=CV, colormap=colormap, colorrange=(0.0, 1.0))
    end
    Colorbar(f[1,2], hm1)    
    
    hSlider = Slider(f[2,:], range = 1:length(ρ_cells), startvalue = length(ρ_cells), linewidth=30)
   
    on(hSlider.value) do i   
        C_E = ρ_cells[i]
        B_E = C_E .> plotThreshold
        F = element2faces(E[B_E])
        C_F = repeat(C_E[B_E],inner=4)
        if plot_type == :elements                          
            CF = FaceView(C_F, eltype(F).(eachindex(F)))
            NF = FaceView(facenormal(F,V), eltype(F).(eachindex(F)))
            M = GeometryBasics.Mesh(V, F, normal=NF, color=CF)
            hm1[1] = M
            ax1.title = "Iteration: $i"
        elseif plot_type == :nodes  
            hm1[1]=GeometryBasics.Mesh(V, F)
            CV = simplex2vertexdata(F, C_F,  V; weighting=:size)
            hm1.color = CV
            ax1.title = "Iteration: $i"
        end
    end 
    # slider2anim(f, hSlider, joinpath(Jutopia_dir(),"assets","temp","temp.gif"); backforth=true, duration=3)
    display(f)    
    return f 
end

plotThreshold = 0.1
fig = final_figure(E, V, ρ_cells, plotThreshold; plot_type=:nodes, colormap = Reverse(:Spectral), strokewidth=0.0, strokecolor=:black)    