# Gravo MG Python Bindings
[[Paper]](https://graphics.tudelft.nl/~klaus/papers/Gravo_MG.pdf) [[Project page]](https://rubenwiersma.nl/gravomg)

![](https://rubenwiersma.nl/assets/img/publications/gravomg/teaser_gravomg.png)

Python binding for Gravo MG. Gravo MG is a Geometric Multigrid Method for solving linear systems on curved surfaces. For more information, check out our [project page](https://rubenwiersma.nl/gravomg).

## How to use
Given positions of points on a mesh or point clouds and a neighborhood graph (e.g., the edges on a mesh), you can use Gravo MG to solve linear systems as follows:
```python
import mgsolver

V, F = # Mesh with vertices and faces
M = # Mass matrix
S = # Stiffness matrix (Laplacian)

# Create the hierarchy
neigh = gravomg.util.neighbors_from_stiffness(S)
solver = gravomg.MultigridSolver(V, neigh, M)

# Solve a linear system
lhs = # Some left hand side, e.g. M + 0.01 * S, make sure it is in csr format.
rhs = # Some right hand side, e.g. M @ V
solution = solver.solve(lhs, rhs)
```

## Installation
The Gravo MG package is available on pip
```bash
pip install gravomg
```

If you would like to recompile yourself, you can clone this repository
```bash
git clone --recurse https://github.com/rubenwiersma/gravo_mg_python.git
```

And install from the folder
```bash
cd gravo_mg_python
pip install ./
```

## Citations
Please cite our paper if this code contributes to an academic publication:

```bib
@Article{WiersmaNasikun2023GravoMG,
author = {Ruben Wiersma, Ahmad Nasikun, Elmar Eisemann, Klaus Hildebrandt},
journal = {SIGGRAPH 2023},
title = {A Fast Geometric Multigrid Method for Curved Surfaces},
year = {2023},
month = jul,
number = {4},
volume = {41},
doi = {10.1145/3588432.3591502},
publisher = {ACM}
}
```
