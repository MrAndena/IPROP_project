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
#include "BlockSchurPreconditioner.hpp"
#include "CollectorGeometryNS.hpp"



using namespace dealii;

template <int dim>
class CompleateProblem
{

public:

    CompleateProblem(parallel::distributed::Triangulation<dim> &tria, const data_struct &d,unsigned short int i);

    void run(); // ha senso mettere toleranze qua ?

private:

    void setup_poisson();
    void assemble_laplace_matrix();
    void assemble_mass_matrix();
    void assemble_nonlinear_poisson();
    void solve_poisson();
    void solve_homogeneous_poisson();
    void solve_nonlinear_poisson(const double tol, const unsigned int max_iterations);

    void setup_drift_diffusion(const bool reinitialize_densities);
    void assemble_drift_diffusion_matrix();
    void apply_drift_diffusion_boundary_conditions(Vector<double> &solution);
    void solve_drift_diffusion();
    void perform_drift_diffusion_fixed_point_iteration_step();

    void setup_navier_stokes();
    void assemble_navier_stokes(const bool nonzero_constraints);
    void solve_nonlinear_navier_stokes_step(const bool nonzero_constraints);
    void navier_stokes_newton_iteration( const double tolerance,const unsigned int max_n_line_searches);
    void solve_navier_stokes();
    void estimate_thrust();

    void evaluate_electric_field();
    void refine_mesh();
    void output_results(const unsigned int step);
    
    // Data for the simulation
    data_struct m_data;
    unsigned short int simulation_tag;
    
    // Parallel data
    MPI_Comm mpi_communicator;
    ConditionalOStream pcout;
    
    // Indexsets
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    // Object for the mesh
    parallel::distributed::Triangulation<dim> &triangulation; 
    
    // DRIFT-DIFFUSION PART

    // FE - DofHandler and Mapping
    FE_Q<dim>       fe;
    DoFHandler<dim> dof_handler;
    MappingQ1<dim>  mapping;

    // Constraints (messi come l'ultimo dd, non ci sono gli elettorni!?)
    AffineConstraints<double> ion_constraints;
    AffineConstraints<double> zero_constraints_poisson;

    // non capisco sti due a che servono
    PETScWrappers::MPI::Vector Field_X; 
    PETScWrappers::MPI::Vector Field_Y;
    
    // Poisson Matrices
    PETScWrappers::MPI::SparseMatrix laplace_matrix_poisson;
    PETScWrappers::MPI::SparseMatrix mass_matrix_poisson;
    PETScWrappers::MPI::SparseMatrix system_matrix_poisson; //manca la density matrix
    
    // Drift-Diffusion Matrices
    PETScWrappers::MPI::SparseMatrix ion_system_matrix;
    PETScWrappers::MPI::SparseMatrix ion_mass_matrix;
    PETScWrappers::MPI::SparseMatrix drift_diffusion_matrix; //NB eliminati entrambi gli sparsity pattern 
    
    // Poisson Vectors
    PETScWrappers::MPI::Vector poisson_newton_update;
    PETScWrappers::MPI::Vector potential;
    PETScWrappers::MPI::Vector poisson_rhs;
    
    // Drift-Diffusion Vectors
    PETScWrappers::MPI::Vector old_ion_density;
    PETScWrappers::MPI::Vector ion_density;
    PETScWrappers::MPI::Vector ion_rhs;
    PETScWrappers::MPI::Vector eta; // chi sei tu ?
    

    // NAVIER-STOKES PART

    // Vectors
    PETScWrappers::MPI::Vector Vel_X;  // in NS originale viene usato Block Vector, qua solo vector, perchè spezza le componenti
    PETScWrappers::MPI::Vector Vel_Y;  // poi però mette anche quello a blocchi
    PETScWrappers::MPI::Vector pressure;

    PETScWrappers::MPI::Vector current_values; // a che serve ?
    
    // NS Parameters
    double viscosity;
    double gamma;
    const unsigned int degree;
    
    // FE - DofHandler and Mapping
    std::vector<types::global_dof_index> dofs_per_block;

    MappingQ1<dim> NS_mapping; // non c'è nel nostro INSIEMEX

    FESystem<dim>    NS_fe;
    DoFHandler<dim>  NS_dof_handler;
    QGauss<dim> volume_quad_formula;   // non cera nell'originale compleate problem
    QGauss<dim - 1> face_quad_formula; // non cera nell'originale complrate problem
    
    // Constraints
    AffineConstraints<double> zero_NS_constraints;
    AffineConstraints<double> nonzero_NS_constraints;
    
    // Matrices
    BlockSparsityPattern      NS_sparsity_pattern;
    PETScWrappers::MPI::BlockSparseMatrix NS_system_matrix;
    PETScWrappers::MPI::SparseMatrix      pressure_mass_matrix;
    
    // BlockVectors
    PETScWrappers::MPI::BlockVector NS_solution;
    PETScWrappers::MPI::BlockVector NS_newton_update; // perchè c'è newton qui?
    PETScWrappers::MPI::BlockVector NS_system_rhs;
    
    std::shared_ptr<BlockSchurPreconditioner> preconditioner; // non cera nell'originale compleate problem
    std::vector<IndexSet> owned_partitioning;                //  non cera nell'originale compleate problem
    std::vector<IndexSet> relevant_partitioning;             // non cera nell'originale compleate problem
    
    // Step and Timestep
    unsigned int step_number = 0;
    double timestep = 0;
    
    // Timer
    Timer timer;
};

// HELPER FUNCTIONS FOR LOCAL TRIANGLE ASSEMBLE (Drift-Diffusion)
void bernoulli (double x, double &bp, double &bn);
double side_length (const Point<2> a, const Point<2> b);
double triangle_denom(const Point<2> a, const Point<2> b, const Point<2> c);
Tensor<1,2> face_normal(const Point<2> a, const Point<2> b);
FullMatrix<double> compute_triangle_matrix(const Point<2> a, const Point<2> b, const Point<2> c, const double alpha12, const double alpha23, const double alpha31, const double D);
// l'originale era un po diversa perche non dava come input D
#include "Compleate_problem_impl.hpp"

//NB: la struttura del preconditioner blockshur che abbiamo usato nel nostro NS è diversa da quella del problema completo originale
//    verificare che siano compatibili

//NB: bisogna riscrivere le funzionio Electrical Values perchè sono settate per una geometria rettangolare

//NB: il nostro CollectoGeometry è diverso ma penso sia ok, menessini passava il segno per descrivere i due profili della naca
//    noi passiamo tutti i punti

//NB: ATTENZIONE ai tag delle bcs ! noi in teoria abbiamo 1-2-3-4 inlet outlet emi coll!!