import numpy as np
from  scipy.sparse import csr_matrix
import gravomg_bindings

from gravomg_bindings import Sampling, Weighting

class MultigridSolver(object):
    def __init__(
        self, pos, neigh, mass,
        ratio=8.0, lower_bound=1000, cycle_type=0, tolerance=1e-4, stopping_criteria=2, pre_iters=2, post_iters=2, max_iter=100, 
        nested=False, sampling_strategy=Sampling.FASTDISK, weighting=Weighting.BARYCENTRIC,
        verbose=False 
        ):
        """Creates the Gravo MG solver for linear systems on curved surfaces (meshes and point clouds).
        
        Args:
            pos (np.ndarray): The positions of the points in the mesh or point cloud.
            neigh (np.ndarray): The neighbors of each point in the mesh or point cloud.
                This should be given as a homogeneous array of size (n_points, max_neighbors),
                padded with -1.
            mass (scipy.sparse.csr_matrix): The mass matrix of the mesh or point cloud.
            ratio (float, optional): The coarsening ratio. Defaults to 8.0.
            lower_bound (int, optional): The lower bound on the number of points in the coarsest level.
                Defaults to 1000.
            cycle_type (int, optional): The type of cycle to use. Defaults to 0 (V-cycle).
            tolerance (float, optional): The tolerance for the stopping criteria. Defaults to 1e-4.
            stopping_criteria (int, optional): The stopping criteria to use. Defaults to 2 (relative residual).
            pre_iters (int, optional): The number of pre-smoothing iterations. Defaults to 2.
            post_iters (int, optional): The number of post-smoothing iterations. Defaults to 2.
            max_iter (int, optional): The maximum number of iterations. Defaults to 100.
            nested (bool, optional): When set to True, does not shift the coarser points to the barycenter.
                Defaults to False.
            sampling_strategy (int, optional): The sampling strategy to use. Defaults to Sampling.FASTDISK.
            weighting (int, optional): The weighting scheme to use. Defaults to Weighting.BARYCENTRIC.
            verbose (bool, optional): Whether to print verbose output. Defaults to False.
        """
        super().__init__()
        if not mass.getformat() == 'csr':
            mass = mass.tocsr()
        self.solver = gravomg_bindings.MultigridSolver(
            pos, neigh, mass,
            ratio, lower_bound,
            cycle_type, tolerance, stopping_criteria, pre_iters, post_iters, max_iter,
            nested, sampling_strategy, weighting,
            verbose
            )

    def solve(self, lhs, rhs):
        """Solves a linear system Ax = b, where lhs is A and rhs is b.

        Args:
            lhs (scipy.sparse.csr_matrix): The left-hand side of the linear system.
            rhs (np.ndarray): The right-hand side of the linear system.
        """
        if not lhs.getformat() == 'csr':
            print('LHS is not in CSR format, converting to CSR')
            lhs = lhs.tocsr()
        return self.solver.solve(lhs, rhs)

    def residual(self, lhs, rhs, solution, type=2):
        """Computes the relative residual of Ax=b,
        where lhs is A, rhs is b and solution is x.
        """
        return self.solver.residual(lhs, rhs, solution, type)

    # Getters and setters

    @property
    def prolongation_matrices(self):
        return self.solver.prolongation_matrices()

    def set_prolongation_matrices(self, U):
        self.solver.set_prolongation_matrices(U)

    @property
    def sampling_indices(self):
        return self.solver.sampling_indices()

    @property
    def all_triangles(self):
        return self.solver.all_triangles()

    @property
    def nearest_source(self):
        return self.solver.nearest_source()