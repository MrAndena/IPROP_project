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

#include <deal.II/lac/petsc_solver.h>


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
    void update_ion_boundary_condition();
    void assemble_drift_diffusion_mass_matrix();
    void assemble_drift_diffusion_matrix(); // build DD matrix
    void solve_drift_diffusion();  // used inside "perform_dd.." to update ion_density
    void perform_drift_diffusion_fixed_point_iteration_step();  // this method update ion_density

    void setup_NS();
    void assemble_NS(bool use_nonzero_constraints, bool assemble_system);
    std::pair<unsigned int, double> solver_NS(bool use_nonzero_constraints, bool assemble_system, double time_step);
    void solve_navier_stokes();
    // void estimate_thrust();

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
    PETScWrappers::MPI::SparseMatrix density_matrix;

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


    FESystem<dim>    NS_fe;
    DoFHandler<dim>  NS_dof_handler;
    QGauss<dim> volume_quad_formula;   // non cera nell'originale complete problem
    QGauss<dim - 1> face_quad_formula; // non cera nell'originale complete problem

    MappingQ1<dim> NS_mapping; // non c'è nel nostro INSIEMEX

    // Constraints
    AffineConstraints<double> zero_NS_constraints;
    AffineConstraints<double> nonzero_NS_constraints;

    // Matrices
    BlockSparsityPattern      NS_sparsity_pattern;
    PETScWrappers::MPI::BlockSparseMatrix  NS_system_matrix;
    PETScWrappers::MPI::BlockSparseMatrix  pressure_mass_matrix;
    PETScWrappers::MPI::BlockSparseMatrix  NS_mass_matrix;

    // BlockVectors
    PETScWrappers::MPI::BlockVector NS_solution;
    PETScWrappers::MPI::BlockVector NS_newton_update; // perchè c'è newton qui?
    PETScWrappers::MPI::BlockVector NS_solution_update;
    PETScWrappers::MPI::BlockVector NS_system_rhs;

    std::shared_ptr<BlockSchurPreconditioner> preconditioner; // non c'era nell'originale complete problem
    std::vector<IndexSet> owned_partitioning;                //  non c'era nell'originale complete problem
    std::vector<IndexSet> relevant_partitioning;             // non c'era nell'originale complete problem

    IndexSet owned_partitioning_p;

    IndexSet NS_locally_relevant_dofs; //nuovo


    double time_NS = 0;
    double timestep_NS;
    mutable TimerOutput timer;
    double timestep = 0;
    SparsityPattern      sparsity_pattern_poisson;
};

// HELPER FUNCTIONS FOR LOCAL TRIANGLE ASSEMBLE (Drift-Diffusion)
void bernoulli (double x, double &bp, double &bn);
double side_length (const Point<2> a, const Point<2> b);
double triangle_denom(const Point<2> a, const Point<2> b, const Point<2> c);
Tensor<1,2> face_normal(const Point<2> a, const Point<2> b);
FullMatrix<double> compute_triangle_matrix(const Point<2> a, const Point<2> b, const Point<2> c, const double alpha12, const double alpha23, const double alpha31, const double D);
Tensor<1,2> get_emitter_normal(const Point<2> a);



// ESPERIMENTO CON TEMPLATE PER LE BCS
template <int dim>
class DynamicBoundaryValues : public Function<dim>
{
public:

    DynamicBoundaryValues(const PETScWrappers::MPI::Vector &vec_1, 
                          const PETScWrappers::MPI::Vector &vec_2,
                          const PETScWrappers::MPI::Vector &vec_3,
                          const DoFHandler<dim> &dh,
                          const MappingQ1<dim>  &map,
                          const data_struct &d): 
    Function<dim>(), pot(vec_1),ion(vec_2),old_ion(vec_3), dof_handler(dh), mapping(map), m_data(d){}

    virtual double value(const dealii::Point<dim> &p, const unsigned int component = 0) const override;
    
private:
    const PETScWrappers::MPI::Vector &pot;
    const PETScWrappers::MPI::Vector &ion;
    const PETScWrappers::MPI::Vector &old_ion;
    const DoFHandler<dim> &dof_handler;
    const MappingQ1<dim>  &mapping;
    const data_struct &m_data;
};


template <int dim>
double DynamicBoundaryValues<dim>::value(const Point<dim> & p,const unsigned int component) const 
{ 

const double E_ON  = m_data.electrical_parameters.E_ON;
const double E_ref = m_data.electrical_parameters.E_ref;
const double N_ref = m_data.electrical_parameters.N_ref;
const double N_min = m_data.electrical_parameters.N_min;

double theta;
const double k_min = 0.5;
const double k_max = 2;

const double ion_norm = ion.l2_norm();
const double old_ion_norm = old_ion.l2_norm();
const double condition = ion_norm / old_ion_norm;

if(condition <= k_max && condition >= k_min){
    theta = 1;
}

if(condition > k_max){
   theta = (k_max -1)/(condition -1);
}

if(condition < k_min){
   theta = (k_min -1)/(condition -1);
}

std::cout<<"theta is: "<<theta<<std::endl;

/*
Functions::FEFieldFunction<dim, dealii::PETScWrappers::MPI::Vector> solution_as_function_object_1(dof_handler, pot, mapping);
Functions::FEFieldFunction<dim, dealii::PETScWrappers::MPI::Vector> solution_as_function_object_2(dof_handler, ion, mapping);

Tensor<1,dim> grad_pot = solution_as_function_object_1.gradient(p);      //gradiente del potenziale in p
const double ion_density_value = solution_as_function_object_2.value(p); //valore della densità nel punto p
Tensor<1,dim> normal = get_emitter_normal(p);                            //versore normale nel punto p

const double En = grad_pot*normal; // campo elettrico normale in p: gradiente del potenziale * normale

const double kappa_over_alpha = N_ref*std::exp((En-E_ON)/E_ref);

const double n = std::max(N_min, kappa_over_alpha);

const double scalar = En*ion_density_value;

const double returned_value = theta*n*scalar + comp_theta*ion_density_value;

return returned_value;
*/
Functions::FEFieldFunction<dim, dealii::PETScWrappers::MPI::Vector> solution_as_function_object_1(dof_handler, ion, mapping);
Functions::FEFieldFunction<dim, dealii::PETScWrappers::MPI::Vector> solution_as_function_object_2(dof_handler, old_ion, mapping);

const double ion_value = solution_as_function_object_1.value(p);
const double old_ion_value = solution_as_function_object_2.value(p);

const double value = theta*old_ion_value +(1-theta)*ion_value;

return value;

}

#include "drift_diffusion_impl.hpp"
