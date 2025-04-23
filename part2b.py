import matplotlib.pyplot as plt
import numpy as np
from petsc4py import PETSc
from mpi4py import MPI
from dolfinx import fem, mesh
from dolfinx.fem.petsc import assemble_matrix, assemble_vector, create_vector, apply_lifting, set_bc
import ufl

# Define the domain and analytical solution
def analytical_solution(x):
    return np.sin(np.pi * x[0]) * np.sin(np.pi * x[1])

# Compute the L^2 norm of the error
def compute_l2_error(u_numerical, V, domain):
    u_analytical = fem.Function(V)
    u_analytical.interpolate(analytical_solution)
    error_form = ufl.inner(u_numerical - u_analytical, u_numerical - u_analytical) * ufl.dx
    error_l2 = np.sqrt(fem.assemble_scalar(fem.form(error_form)))
    norm_analytical = np.sqrt(fem.assemble_scalar(fem.form(ufl.inner(u_analytical, u_analytical) * ufl.dx)))
    relative_error = error_l2 / norm_analytical if norm_analytical > 0 else error_l2
    return relative_error

# Solve the Poisson equation for a given mesh and polynomial degree
def solve_poisson(nx, ny, degree, cell_type, case_name):
    # Create mesh
    domain = mesh.create_rectangle(MPI.COMM_WORLD, [np.array([0, 0]), np.array([1, 1])], [nx, ny], cell_type)
    
    # Define function space
    V = fem.functionspace(domain, ("Lagrange", degree))
    
    # Boundary condition (u = 0 on the boundary)
    fdim = domain.topology.dim - 1
    boundary_facets = mesh.locate_entities_boundary(
        domain, fdim, lambda x: np.full(x.shape[1], True, dtype=bool))
    bc = fem.dirichletbc(PETSc.ScalarType(0), fem.locate_dofs_topological(V, fdim, boundary_facets), V)
    
    # Define variational problem: -Δu = f
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    x = ufl.SpatialCoordinate(domain)  # Define spatial coordinates
    f = 2 * np.pi**2 * ufl.sin(np.pi * x[0]) * ufl.sin(np.pi * x[1])  # Right-hand side
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx  # Bilinear form
    L = f * v * ufl.dx  # Linear form
    
    # Assemble system
    bilinear_form = fem.form(a)
    linear_form = fem.form(L)
    A = assemble_matrix(bilinear_form, bcs=[bc])
    A.assemble()
    b = create_vector(linear_form)
    
    # Assemble right-hand side
    with b.localForm() as loc_b:
        loc_b.set(0)
    assemble_vector(b, linear_form)
    apply_lifting(b, [bilinear_form], [[bc]])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    set_bc(b, [bc])
    
    # Solve
    uh = fem.Function(V)
    solver = PETSc.KSP().create(domain.comm)
    solver.setOperators(A)
    solver.setType(PETSc.KSP.Type.PREONLY)
    solver.getPC().setType(PETSc.PC.Type.LU)
    solver.solve(b, uh.x.petsc_vec)
    uh.x.scatter_forward()
    
    # Compute L^2 error
    relative_error = compute_l2_error(uh, V, domain)
    return relative_error

# Define cases
h_refinement_cases = [
    (5, 5, 1, mesh.CellType.quadrilateral, "h-ref: 5x5, p=1"),
    (20, 20, 1, mesh.CellType.quadrilateral, "h-ref: 20x20, p=1"),
    (50, 50, 1, mesh.CellType.quadrilateral, "h-ref: 50x50, p=1"),
    (100, 100, 1, mesh.CellType.quadrilateral, "h-ref: 100x100, p=1")
]

p_refinement_cases = [
    (20, 20, 1, mesh.CellType.triangle, "p-ref: degree=1, tri"),
    (20, 20, 3, mesh.CellType.triangle, "p-ref: degree=3, tri"),
    (20, 20, 4, mesh.CellType.quadrilateral, "p-ref: degree=4, quad"),
    (20, 20, 6, mesh.CellType.quadrilateral, "p-ref: degree=6, quad")
]

all_cases = h_refinement_cases + p_refinement_cases

# Store results
relative_errors = []
case_names = []

# Run simulations for all cases
for nx, ny, degree, cell_type, case_name in all_cases:
    print(f"Running case: {case_name}")
    rel_err = solve_poisson(nx, ny, degree, cell_type, case_name)
    relative_errors.append(rel_err)
    case_names.append(case_name)

# Plot comparison of relative errors
plt.figure(figsize=(10, 6))
bars = plt.bar(range(len(case_names)), relative_errors, tick_label=case_names)
plt.xlabel('Case')
plt.ylabel('Relative L^2 Error')
plt.title('Relative L^2 Error Comparison for h- and p-Refinement')
plt.yscale('log')  # Use logarithmic scale for better visualization
plt.xticks(rotation=45, ha='right')
plt.grid(True, which="both", ls="--")

# Add error values on top of bars
for bar, error in zip(bars, relative_errors):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height(), f'{error:.2e}',
             ha='center', va='bottom')

plt.tight_layout()
plt.savefig('relative_l2_error_comparison_poisson.png')
plt.show()

# Print relative errors
print("\nRelative L^2 Errors:")
for name, error in zip(case_names, relative_errors):
    print(f"{name}: {error:.2e}")
