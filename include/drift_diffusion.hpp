#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h> // for the timer

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h> // For neighbor renumbering

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_interface_values.h> // For gradient evaluator

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h> // To use manifolds
#include <deal.II/grid/grid_in.h> // For GMSH
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_direct.h> // For UMFPACK
#include <deal.II/lac/solver_gmres.h> // For GMRES
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/sparse_ilu.h> // ILU preconditioning

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h> // For Laplace Matrix
#include <deal.II/numerics/fe_field_function.h> // For boundary values
#include <deal.II/numerics/solution_transfer.h> // For the solution transfer
#include <deal.II/numerics/error_estimator.h> // Kelly error estimator

#include <fstream>
#include <cmath>
#include "data_struct.hpp"
#include "CollectorGeometryNS.hpp"



using namespace dealii;

template <int dim>
class drift_diffusion
{

public:

    drift_diffusion(parallel::distributed::Triangulation<dim> &tria, const data_struct &d);

    void run();

private:

    void setup_poisson();
    void assemble_initial_system();
    void assemble_poisson_laplace_matrix();
    void assemble_poisson_mass_matrix();
    void assemble_nonlinear_poisson();
    double solve_poisson(); // update poisson and eta. it return the residual (L_INF norm of newton_update)
    void solve_homogeneous_poisson(); 
    void solve_nonlinear_poisson(const unsigned int max_iterations,const double tol); // update poisson and eta

    void setup_drift_diffusion(); //setup DD
    void assemble_drift_diffusion_mass_matrix();
    void assemble_drift_diffusion_matrix(); // build DD matrix
    void solve_drift_diffusion();  // used inside "perform_dd.." to update ion_density
    void perform_drift_diffusion_fixed_point_iteration_step();  // this method update ion_density

    void evaluate_electric_field(); // usato sia in assemble_DD che in solve_NS 
    void output_results(const unsigned int step); // preso dal nostro DD dovrebbe funzionare dovrebbere essere const method no?
    
    // Data for the simulation
    data_struct m_data;
    
    // Parallel data
    MPI_Comm mpi_communicator;
    ConditionalOStream pcout;
    
    // Object for the mesh
    parallel::distributed::Triangulation<dim> &triangulation; 
    
    // DRIFT-DIFFUSION PART

    // Indexsets
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    // FE - DofHandler and Mapping
    FE_Q<dim>       fe;
    DoFHandler<dim> dof_handler;
    MappingQ1<dim>  mapping;

    // Constraints (messi come l'ultimo dd, non ci sono gli elettorni)
    AffineConstraints<double> ion_constraints;
    AffineConstraints<double> zero_constraints_poisson;
    AffineConstraints<double> constraints_poisson;
    
    // Poisson Matrices
    PETScWrappers::MPI::SparseMatrix laplace_matrix_poisson;
    PETScWrappers::MPI::SparseMatrix mass_matrix_poisson;
    PETScWrappers::MPI::SparseMatrix system_matrix_poisson;

    PETScWrappers::MPI::SparseMatrix initial_matrix_poisson;
    
    // Drift-Diffusion Matrices
    PETScWrappers::MPI::SparseMatrix ion_system_matrix;
    PETScWrappers::MPI::SparseMatrix ion_mass_matrix;
    
    // Poisson Vectors
    PETScWrappers::MPI::Vector poisson_newton_update;
    PETScWrappers::MPI::Vector potential;
    PETScWrappers::MPI::Vector poisson_rhs;

    PETScWrappers::MPI::Vector initial_poisson_rhs;

    PETScWrappers::MPI::Vector Field_X; //electric field
    PETScWrappers::MPI::Vector Field_Y;
    
    // Drift-Diffusion Vectors
    PETScWrappers::MPI::Vector old_ion_density;
    PETScWrappers::MPI::Vector ion_density;
    PETScWrappers::MPI::Vector ion_rhs;
    PETScWrappers::MPI::Vector eta; 
    
    //mutable TimerOutput timer;
    double timestep = 0;
    SparsityPattern      sparsity_pattern_poisson;
};

// HELPER FUNCTIONS FOR LOCAL TRIANGLE ASSEMBLE (Drift-Diffusion)
void bernoulli (double x, double &bp, double &bn);
double side_length (const Point<2> a, const Point<2> b);
double triangle_denom(const Point<2> a, const Point<2> b, const Point<2> c);
Tensor<1,2> face_normal(const Point<2> a, const Point<2> b);
FullMatrix<double> compute_triangle_matrix(const Point<2> a, const Point<2> b, const Point<2> c, const double alpha12, const double alpha23, const double alpha31, const double D);



#include "drift_diffusion_impl.hpp"
