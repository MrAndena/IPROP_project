#ifndef INSIMEX_HPP
#define INSIMEX_HPP

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_point_data.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>


#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/solver_bicgstab.h>  //Se vogliamo usare bicgstab solver

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <filesystem>

#include "BlockSchurPreconditioner.hpp"
#include "Time.hpp"
#include "BoundaryValues.hpp"

// The incompressible Navier-Stokes solver
template <int dim>
class InsIMEX
{
public:
  InsIMEX(parallel::distributed::Triangulation<dim> &tria);              //We deserve a parallel::DistributedTriangulation so our mesh must be in this class or similar to parallelize it
  void run();
  ~InsIMEX() { timer.print_summary(); }

private:
  void setup_dofs();

  void make_constraints_init();
  void make_constraints_update();

  void initialize_system();
  void assemble(bool use_nonzero_constraints, bool assemble_system);

  std::pair<unsigned int, double> solve(bool use_nonzero_constraints, bool assemble_system, double time_step);
  
  void output_results(const unsigned int) const;

  double viscosity;
  double gamma;
  const unsigned int degree;
  std::vector<types::global_dof_index> dofs_per_block;

  parallel::distributed::Triangulation<dim> &triangulation;
  FESystem<dim> fe;
  DoFHandler<dim> dof_handler;                             //The DoFHandler object that describes which degrees of freedom live on which cells.
  QGauss<dim> volume_quad_formula;
  QGauss<dim - 1> face_quad_formula;

  AffineConstraints<double> initial_NS_constraints;              //Each "line" in objects of this class corresponds to one constrained degree of freedom
  AffineConstraints<double> update_NS_constraints;

  BlockSparsityPattern sparsity_pattern;
  
  PETScWrappers::MPI::BlockSparseMatrix system_matrix;     // System matrix to be solved
  
  PETScWrappers::MPI::BlockSparseMatrix mass_matrix;       // Block matrix which includes both velocity mass matrix and pressure mass matrix.
  
  PETScWrappers::MPI::BlockSparseMatrix mass_schur;        // The schur complement of mass matrix is not a block matrix. It is defined as a block matrix where only one block is actually used to reuse the partition we created for the system matrix

  PETScWrappers::MPI::BlockVector present_solution;        // The latest known solution.
  
  PETScWrappers::MPI::BlockVector solution_increment;      // The increment at a certain time step.
  
  PETScWrappers::MPI::BlockVector system_rhs;              // System RHS

  MPI_Comm mpi_communicator;

  ConditionalOStream pcout;                                //A class allows you to print an output stream, useful in parallel computations

  std::vector<IndexSet> owned_partitioning;                // The IndexSets of owned velocity and pressure respectively.

  std::vector<IndexSet> relevant_partitioning;             // The IndexSets of relevant velocity and pressure respectively.

  IndexSet locally_relevant_dofs;                          //IndexSet: a class that represents a subset of indices among a larger set. 

  std::shared_ptr<BlockSchurPreconditioner> preconditioner;   // The BlockSchurPreconditioner for the entire system.

  Time time;
  mutable TimerOutput timer;
};

#include "InsIMEX_impl.hpp"

#endif