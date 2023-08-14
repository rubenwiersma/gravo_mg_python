#include "gravomg/multigrid_solver.h"
#include "gravomg/utility.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Eigen/Sparse"

namespace py = pybind11;

class MultigridSolver {

public:
    // Expects positions as N x 3 matrix and neighbors as N x K,
    // where K is the maximum number of neighbors.
    // Row i contains the indices of the neighbors of node i;
    // the neighbors should be padded with -1 to get K entries per row.
    MultigridSolver(
        Eigen::MatrixXd positions, Eigen::MatrixXi neighbors, Eigen::SparseMatrix<double> mass,
        double ratio, int low_bound, // Hierarchy settings
        int cycle_type, double tolerance, int stopping_criteria, int pre_iters, int post_iters, int max_iter, // Solver settings
        bool nested, Sampling sampling_strategy, Weighting weighting, bool verbose // Further customization
        ) {
        solver.reset(new GravoMG::MultigridSolver(positions, neighbors, mass));
        solver->nested = nested;
        solver->samplingStrategy = sampling_strategy;
        solver->weightingScheme = weighting;
        solver->maxIter = max_iter;

        solver->verbose = verbose;

        // Building hierarchy
        solver->ratio = ratio;
        solver->lowBound = low_bound;
        solver->buildHierarchy();

        // Set solver settings
        solver->cycleType = cycle_type;
        solver->accuracy = tolerance;
        solver->stoppingCriteria = stopping_criteria;
        solver->preIters = pre_iters;
        solver->postIters = post_iters;
        solver->isSmootherGaussSeidel = true;
    }

    Eigen::MatrixXd solve(Eigen::SparseMatrix<double> lhs, Eigen::MatrixXd rhs) {
        Eigen::MatrixXd x = rhs;
        solver->solve(lhs, rhs, x);
        return x;
    }
    
    // -- Metrics

    double residual(Eigen::SparseMatrix<double> lhs, Eigen::MatrixXd rhs, Eigen::MatrixXd solution, int type=2) {
        return solver->residualCheck(lhs, rhs, solution, type);
    }

    //-- Data access

    std::vector<Eigen::SparseMatrix<double>> prolongation_matrices() {
        return solver->U;
    }

    void set_prolongation_matrices(std::vector<Eigen::SparseMatrix<double>> U) {
        solver->U = U;
    }

    std::vector<std::vector<int>> sampling_indices() {
        return solver->samples;
    }

    std::vector<std::vector<std::vector<int>>> all_triangles() {
        return solver->allTriangles;
    }

    std::vector<std::vector<size_t>> nearest_source() {
        return solver->nearestSource;
    }

private:
    std::unique_ptr<GravoMG::MultigridSolver> solver; 
};


PYBIND11_MODULE(gravomg_bindings, m) {
    m.doc() = "Multigrid solver bindings";

    py::class_<MultigridSolver>(m, "MultigridSolver")
        .def(py::init<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::SparseMatrix<double>, double, int, int, double, int, int, int, int, bool, Sampling, Weighting, bool>())
        .def("solve", &MultigridSolver::solve, py::arg("lhs"), py::arg("rhs"))
        .def("residual", &MultigridSolver::residual, py::arg("lhs"), py::arg("rhs"), py::arg("solution"), py::arg("type") = 2)
        .def("prolongation_matrices", &MultigridSolver::prolongation_matrices)
        .def("set_prolongation_matrices", &MultigridSolver::set_prolongation_matrices, py::arg("U"))
        .def("sampling_indices", &MultigridSolver::sampling_indices)
        .def("all_triangles", &MultigridSolver::all_triangles)
        .def("nearest_source", &MultigridSolver::nearest_source);

    py::enum_<Sampling>(m, "Sampling")
        .value("FASTDISK", FASTDISK)
        .value("RANDOM", RANDOM)
        .value("MIS", MIS);

    py::enum_<Weighting>(m, "Weighting")
        .value("BARYCENTRIC", BARYCENTRIC)
        .value("UNIFORM", UNIFORM)
        .value("INVDIST", INVDIST);
}
