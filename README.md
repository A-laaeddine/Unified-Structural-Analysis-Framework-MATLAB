# Unified 2D Structural Analysis Framework (MATLAB)

This repository contains a comprehensive, modular MATLAB implementation of the **Direct Stiffness Method (DSM)** for the linear static analysis of 2D structural systems, including **trusses, beams, and frames**.

Developed as a unified computational tool, this solver bridges the gap between theoretical matrix structural analysis and practical numerical execution. It automates the transition from local member properties to global system equilibrium, providing an end-to-end pipeline from input modeling to post-processing visualization.



---

## 🚀 Key Features

* **Unified Solver:** A single algorithmic core capable of handling three distinct structural types (Truss, Beam, Frame) by dynamically adjusting Degrees of Freedom (DOF).
* **Automated Load Treatment:** Built-in logic to convert member-level loads (Uniformly Distributed Loads and mid-span point loads) into equivalent nodal force vectors using fixed-end force transformations.
* **Numerical Efficiency:** Utilizes vectorized assembly and MATLAB’s `sparse` matrix storage to optimize memory usage and computational speed for larger systems.
* **Robust Partitioning:** Implements automated system partitioning to handle varied boundary conditions (pinned, roller, and fixed supports).
* **Visual Post-Processing:** Includes dedicated routines for plotting structural geometry and support conditions.

---

## 📂 Repository Structure

### Core Script
* **`main_analysis.m`**: The primary entry point. It orchestrates the workflow: initializing data, calling the assembly functions, solving the system of linear equations ($KU = F$), and outputting results for displacements and reactions.

### Helper Functions
To maintain modularity, the framework is supported by several specialized functions:
* **`manualInput.m`**: Configures the structural topology, material properties ($E, A, I$), and loading conditions.
* **`stiffness_matrix.m`**: Computes the local stiffness matrix for elements based on their type.
* **`transformation_matrix.m`**: Handles the rotation of local matrices into the global coordinate system.
* **`equivalent_loads.m`**: Processes member loads and maps them to the global load vector.
* **`supportCenterPoint.m`**: A graphical utility that ensures support icons and reaction vectors are aligned correctly in post-processing plots.



---

## 🛠️ Usage

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/YourUsername/Unified-Structural-Analysis.git](https://github.com/YourUsername/Unified-Structural-Analysis.git)
    ```
2.  Open **MATLAB** and navigate to the repository folder.
3.  Define your structural parameters (nodes, members, loads) within the `manualInput.m` function.
4.  Execute `main_analysis.m` to run the solver.

---

## 📊 Validation

The solver has been rigorously benchmarked against analytical solutions from standard engineering texts (e.g., *Kassimali*). Validation cases included:
* Multi-span continuous beams.
* Symmetric and unsymmetric portal frames.
* Indeterminate planar trusses.

All test cases demonstrated a high degree of numerical accuracy, typically maintaining an error margin of **< 1%**.

---

## 📝 License

This project is licensed under the **MIT License** - see the `LICENSE` file for details.

## 🤝 Contact

**Mohamad Alaaeddine** - [Your Email/LinkedIn]  
Project Link: [https://github.com/YourUsername/Unified-Structural-Analysis](https://github.com/YourUsername/Unified-Structural-Analysis)
