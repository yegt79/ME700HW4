# Assignment 4: Finite Element Analysis with FEniCSx

This repository contains tutorials for **Assignment 4**, focusing on solving partial differential equations (PDEs) using FEniCSx, with comparisons to analytical solutions and studies on mesh refinement and boundary condition effects. The tutorials are divided into three parts: **Part A**, **Part B**, and **Part C**, each addressing different aspects of finite element analysis (FEA).

---

## Part A: 2D Heat Equation on a Square with Analytical Comparison

**Tutorial_Part_A** demonstrates the solution of the **2D Heat Equation** on a square domain, comparing numerical results obtained with FEniCSx to analytical solutions.

- **Methodology**:
  - The domain is visualized over time using **PyVista**.
  - Maximum temperature is computed for both the analytical solution and the numerical solution using FEniCSx.
  
- **Results**:
  - The numerical results from FEniCSx closely match the analytical solution, confirming the accuracy of FEniCSx in capturing the heat equation dynamics.

---

## Part B: 2D Poisson Equation with h- and p-Refinement Study

**Tutorial_Part_B** explores the **2D Poisson Equation** with a focus on **h-refinement** (mesh refinement) and **p-refinement** (polynomial degree refinement) to study their impact on solution accuracy.

- **Refinement Cases**:
  - **h-refinement** (increasing mesh resolution, fixed polynomial degree):
    - Mesh sizes: `(5x5, p=1)`, `(20x20, p=1)`, `(50x50, p=1)`, `(100x100, p=1)` with quadrilateral elements.
  - **p-refinement** (increasing polynomial degree, fixed mesh):
    - Cases: `(20x20, p=3, triangle)`, `(20x20, p=6, triangle)`, `(20x20, p=4, quadrilateral)`.

- **Analysis**:
  - The **L² error** is calculated for analytical and numerical solutions across all refinement cases.
  - Results are visualized in a **GIF**, showcasing the domain for each case.

- **Findings**:
  - Both h- and p-refinements improve the accuracy of FEniCSx solutions, with error comparisons demonstrating the effectiveness of refinement strategies.

---

## Part C: Failures in FEniCSx Due to Mesh and Boundary Conditions

**Tutorial_Part_C** examines two scenarios where FEniCSx fails to produce accurate results due to poor mesh quality or inappropriate boundary conditions.

### Tutorial_Part_C_Mesh: Poor Mesh Quality
- **Issue**:
  - A very low number of mesh elements leads to inadequate resolution of the domain.
- **Result**:
  - The GIF visualization of the solution clearly shows that the model fails to capture the physical behavior accurately due to insufficient mesh refinement.

### Tutorial_Part_C_Boundary: Inappropriate Boundary Conditions
- **Issue**:
  - A Gaussian hill is modeled with **Neumann boundary conditions** that trap the hill, preventing it from flattening naturally.
- **Result**:
  - The solver converges, but the solution is **unphysical** because the boundary conditions do not allow proper dissipation of the Gaussian hill.

---

## Conclusion

This assignment demonstrates the capabilities and limitations of FEniCSx in solving PDEs:
- **Part A** validates FEniCSx against analytical solutions for the heat equation.
- **Part B** highlights the importance of mesh and polynomial refinements in improving solution accuracy.
- **Part C** underscores the critical role of proper mesh quality and boundary conditions in achieving physically meaningful results.

For detailed implementations, refer to the respective tutorial files in the repository.

---

# Download FEniCSx

**Note**: Additional installs since last time!

**Recommendation**: I strongly recommend working on the **SCC** (Shared Computing Cluster).

### Create a Conda Environment

See: [FEniCSx Download Page](https://fenicsproject.org/download/)

On the SCC, run the following commands to set up the environment:

```bash
module load miniconda
mamba create -n fenicsx-env
mamba activate fenicsx-env
mamba install -c conda-forge fenics-dolfinx mpich pyvista
pip install imageio
pip install gmsh
pip install PyYAML
```

**Note**: This installation may take a couple of minutes.

Once installed, you can launch a **VSCode server** and select `fenicsx-env` as your Conda environment.

## Automation Example

1. Run `script.py`.
2. Consider how you would perform a mesh refinement study using this script.
3. Refer to `updated_script_*.py` for an example of implementing a mesh refinement study.
4. Use `run_simulations.sh` for a bash script example to run an array job for multiple simulations simultaneously. Execute it with:

```bash
qsub run_simulations.sh
```

To monitor jobs, use commands like:

```bash
qstat -u <username>
```

For more examples tailored to your needs, see: [SCC Batch Script Examples](https://www.bu.edu/tech/support/research/system-usage/running-jobs/batch-script-examples/)

**Note**: For the job scheduler, avoid submitting many tiny jobs. Aim for job durations between **10 minutes and 12 hours** (for this class). You can use the bash script to run subsets of array jobs in series as needed.