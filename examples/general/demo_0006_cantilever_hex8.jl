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

# Function to create a 3D grid
function create_grid(Lx,  Ly, Lz, nx, ny, nz)
    # Define the domain boundaries using Vec{3}
    left = Ferrite.Vec(0.0, 0.0, 0.0)
    right = Ferrite.Vec(Lx, Ly, Lz)
    # Generate the 3D hexahedral grid
    grid = generate_grid(Hexahedron, (nx, ny, nz), left, right)
    # Define a small load region in the center of the y-z plane at x = L
    addnodeset!(grid, "load", x -> x[1] ≈Lx && x[3]<= 0.0)
    return grid
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
    dbcs = ConstraintHandler(dh)
    # Clamped on the left side
    dofs = [1, 2, 3]
    dbc = Dirichlet(:u, getfacetset(grid, "left"), (x,t) -> [0.0, 0.0, 0.0], dofs)
    add!(dbcs, dbc)
    close!(dbcs)
    return dbcs
end;


# Main script
Lx, Ly, Lz = 30.0, 2.0, 10.0
nx, ny, nz = 40, 5, 10


grid = create_grid(Lx, Ly, Lz, nx, ny, nz)
F, V = FerriteToComodo(grid, Ferrite.Hexahedron)


input.grid = grid
input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh, grid)

input.cv, input.fv = create_values()

input.volfrac = 0.5
input.penalty = 3.0

input.rmin = 1.5

input.nodal_vector = (0.0, 0.0, -1.0)
input.nodeset = getnodeset(grid, "load")
input.max_iter = 200
input.tol = 0.01
input.ν = 0.3
input.Emin = 1e-9
input.E0 = 1.0

input.filter_type = :sensitivity_filter
input.model_type = :dim3d
input.load_type = :nodal
ρ_cells = run_optimization(input)

VTKGridFile("linear_elasticity", input.dh) do vtk
    write_cell_data(vtk, ρ_cells[end], "density")
    Ferrite.write_cellset(vtk, grid)
end
