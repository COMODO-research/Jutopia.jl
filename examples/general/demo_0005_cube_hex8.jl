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

# Function to check if a point is inside a circle on a given plane
function in_circle(x, r, cx, cy, cz)
    return (x[3] ≈ cz) && ((x[1] - cx)^2 + (x[2] - cy)^2 <= r^2)
end

# Function to create a 3D grid
function create_grid(Lx, Ly, Lz, nx, ny, nz)
    # Define the domain boundaries using Vec{3}
    left = Ferrite.Vec(0.0, 0.0, 0.0)
    right = Ferrite.Vec(Lx, Ly, Lz)

    # Generate the 3D hexahedral grid
    grid = generate_grid(Hexahedron, (nx, ny, nz), left, right)
    return grid
end

# Function to define boundary conditions and nodesets
function create_boundary(grid, Lx, Ly, Lz)
    # Parameters for the circle on the top surface
    r_top = 25.0
    cx_top, cy_top = Lx / 2, Ly / 2

    # Parameters for the circles on the bottom surface
    r_bottom = 10.0

    # Define the bottom corner centers
    corner_centers = [
        (0.0, 0.0),
        (Lx, 0.0),
        (0.0, Ly),
        (Lx, Ly)
    ]

    # Add nodeset for the top surface inside the circle
    # addnodeset!(grid, "top_circle", x -> in_circle(x, r_top, cx_top, cy_top, Lz))
    addfacetset!(grid, "top_circle", x -> in_circle(x, r_top, cx_top, cy_top, Lz))

    # Add nodesets for each bottom corner circle
    for (i, (cx, cy)) in enumerate(corner_centers)
        addnodeset!(grid, "bottom_corner_circle_$i", x -> in_circle(x, r_bottom, cx, cy, 0.0))
    end
end

# Function to create cell and facet values
function create_values()
    order = 1
    dim = 3
    ip = Lagrange{RefHexahedron, order}()^dim
    qr = QuadratureRule{RefHexahedron}(2)
    qr_face = FacetQuadratureRule{RefHexahedron}(1)
    cell_values = CellValues(qr, ip)
    facet_values = FacetValues(qr_face, ip)
    return cell_values, facet_values
end

# Function to create a DOF handler
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefHexahedron, 1}()^3)
    Ferrite.close!(dh)
    return dh
end

# Function to create boundary conditions
function create_bc(dh, grid)
    ch = Ferrite.ConstraintHandler(dh)
    for i in 1:4
        dbc = Dirichlet(:u, getnodeset(grid, "bottom_corner_circle_$i"), (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3])
        add!(ch, dbc)
    end
    Ferrite.close!(ch)
    return ch
end

Lx, Ly, Lz = 100.0, 100.0, 100.0
nx, ny, nz = 15, 15, 15

grid = create_grid(Lx, Ly, Lz, nx, ny, nz)
F, V = FerriteToComodo(grid, Ferrite.Hexahedron)

create_boundary(grid, Lx, Ly, Lz)

input.grid = grid
input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh, grid)

input.cv, input.fv = create_values()

input.volfrac = 0.5
input.penalty = 3.0

input.rmin = 0.05

input.traction_vector = (0.0, 0.0, 1.0)
input.facetset = getfacetset(grid, "top")
input.max_iter = 200
input.tol = 0.01
input.ν = 0.3
input.Emin = 1e-9
input.E0 = 1.0

input.filter_type = :sensitivity_filter
input.model_type = :dim3d
input.load_type = :traction
ρ_cells = run_optimization(input)

VTKGridFile("linear_elasticity", input.dh) do vtk
    write_cell_data(vtk, ρ_cells[end], "density")
    Ferrite.write_cellset(vtk, grid)
end
