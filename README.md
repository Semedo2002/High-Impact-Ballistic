* Solver: 3D explicit FEA with Newmark-beta time integration.
* Contact Model: Penalty-based force calculation and detection.
* Visualization: Dynamic 3D rendering with Von Mises stress-gradient coloring.
* Analysis: Quantitative assessment of pressure distribution and energy conservation.

Execution 

1. ballistic_main.m: Primary engine. Run first to compute simulation data .
2. ballistic_visualize.m: Post tool for 3D animation and deformation review.
3. ballistic_plots.m: Generates analytical figures for stress evolution and damage.
Default Parameters

- Projectile: 5 kg, 800 m/s, Steel (E=200 GPa, Yield=350 MPa).
- Target Plate: 1m x 1m x 5cm, Hardened Steel (E=210 GPa, Yield=500 MPa).
- Computational: 1 microsecond time step, 19,200 nodes.

Analytical Interpretation

- Stress Levels: Blue/Green (<200 MPa) is safe; Red (>500 MPa) indicates failure.
- Results: The console will output "PLATE INTACT" or 'PLATE FAILURE' based on yield criteria.
- Energy: Evaluates kinetic energy loss and dissipation through plastic deformation.

Theoretical Basis

- Governing Equations: Newton’s Second Law and Hooke’s Law.
- Yield Criterion: Von Mises stress formulation.
- Contact: Spring-damper penalty method for interface dynamics.

System Requirements
- Preferred Minimum 8GB RAM.
- Results are saved to 'ballistic_results.mat'.
