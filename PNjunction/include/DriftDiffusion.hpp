#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>                 // for the timer
#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/conditional_ostream.h>  // serve per usare ConditionalOStream pcout
#include <deal.II/base/index_set.h>            // serve per la classe indexset
#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h> // For UMFPACK
#include <deal.II/lac/vector.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_tools.h> // per distribute_sparsity_pattern
//#include <deal.II/lac/matrix_out.h> // For matrix output
  
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h> // For Laplace Matrix
#include <deal.II/numerics/error_estimator.h>

//To parallelize
#include <deal.II/distributed/grid_refinement.h>     //tools to operate on parallel distributed triangulations
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>                //parallel distributed triangulation


#include <fstream>
#include <iostream>
#include <limits>

#include "Electrical_Constants.hpp"   
#include "Electrical_Values.hpp"

namespace {

  using namespace dealii;
  
  template <int dim>
    class DriftDiffusion
    {
    public:

      DriftDiffusion(parallel::distributed::Triangulation<dim> &tria);

      void run(const unsigned int max_iter, const double dd_toll); // we pass to run the tolerance and the max number of iterations for newton

    private:
  
      void setup_system();   //for both Poisson and Drift problems

      void initialization(); //for potential, electron_density and hole density

      void assemble_laplace_matrix();
      void assemble_mass_matrix();
      void assemble_nonlinear_poisson();
      void assemble_drift_diffusion_matrix();
      
      void compute_densities();
      void cycle_newton_poisson(const unsigned int max_iter_newton, const double toll_newton);

      void solve_poisson();
      void apply_drift_diffusion_boundary_conditions(); 
      void solve_drift_diffusion();

      void output_results(const unsigned int cycle);


      MPI_Comm mpi_communicator;

      parallel::distributed::Triangulation<dim> &triangulation;

      FE_Q<dim>       fe;
      DoFHandler<dim> dof_handler;
      
      //Whole Problem Index Sets 
      IndexSet locally_owned_dofs;
      IndexSet locally_relevant_dofs;
      
      //Poisson Problem AffineConstraints
      AffineConstraints<double> zero_constraints;
   
      
      //Drift Diffusion AffineConstraints
      //AffineConstraints<double> density_constraints;
      AffineConstraints<double> elec_constraints; 
      AffineConstraints<double> hole_constraints; 

      
      //Poisson Problem Matrices
      PETScWrappers::MPI::SparseMatrix poisson_system_matrix;
      PETScWrappers::MPI::SparseMatrix laplace_matrix;
      PETScWrappers::MPI::SparseMatrix mass_matrix;
      PETScWrappers::MPI::SparseMatrix density_matrix;

      //Drift Diffusion Matrices
      PETScWrappers::MPI::SparseMatrix electron_matrix;
      PETScWrappers::MPI::SparseMatrix hole_matrix;



      
      //Poisson Problem Vectors
      PETScWrappers::MPI::Vector       current_solution;  // potential
      PETScWrappers::MPI::Vector       newton_update;     
      PETScWrappers::MPI::Vector       poisson_system_rhs;
      
      //Drift Diffusion Vectors
      PETScWrappers::MPI::Vector       electron_density;
      PETScWrappers::MPI::Vector       hole_density;

      PETScWrappers::MPI::Vector       old_electron_density;
      PETScWrappers::MPI::Vector       old_hole_density;

      PETScWrappers::MPI::Vector       rhs_hole_density;
      PETScWrappers::MPI::Vector       rhs_electron_density;

      ConditionalOStream pcout;
      MappingQ1<dim> mapping;
      
      Timer timer;
    };

    // HELPER FUNCTION
    void bernoulli (double x, double &bp, double &bn);
    double side_length (const Point<2> a, const Point<2> b);
    double triangle_denom(const Point<2> a, const Point<2> b, const Point<2> c);
    Tensor<1,2> face_normal(const Point<2> a, const Point<2> b);
    FullMatrix<double> compute_triangle_matrix(const Point<2> a, const Point<2> b, const Point<2> c, const double alpha12, const double alpha23, const double alpha31, const double D);


} // namespace
#include "DriftDiffusion_impl.hpp"
