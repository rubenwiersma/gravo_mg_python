# Gravo MG: Graph Voronoi Multigrid
[[Paper]](https://graphics.tudelft.nl/~klaus/papers/Gravo_MG.pdf) [[Project page]](https://rubenwiersma.nl/gravomg)

![](https://rubenwiersma.nl/assets/img/publications/gravomg/teaser_gravomg.png)

Repository for **"A Fast Geometric Multigrid Method for Curved Surfaces"**, published at **SIGGRAPH 2023**
<br />
by Ruben Wiersma, [Ahmad Nasikun](https://github.com/a-nasikun) (equal contribution); Elmar Eisemann; Klaus Hildebrandt.

If you need to solve linear systems on meshes or point clouds with >50.000 vertices, Gravo MG is what you need. Gravo MG is a fast geometric multigrid method that quickly computes a hierarchy used in an iterative multigrid solver. The use of graph Voronoi cells enables fast convergence, while allowing for fast construction.

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
