template <int dim>
drift_diffusion<dim>::drift_diffusion(parallel::distributed::Triangulation<dim> &tria,  // triangulation for the mesh
                                                             const data_struct &d)  // data struct from the user
  : m_data(d)
  , mpi_communicator(MPI_COMM_WORLD)
  , pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  , triangulation(tria)
  , fe(1) //fe for poisson / DD
  , dof_handler(tria) //dof h for poisson / DD
  , mapping() // mapping for poisson / DD
  , timestep(1e-2) // timestep for the DD algorithm
{}

//---------------------------------------------------------------------------------------------------------------------------------------------------------

template <int dim>
void drift_diffusion<dim>::setup_poisson()
{

//fix the constant
double const Ve = m_data.electrical_parameters.Ve;

const double Re = 30e-5; // ??? DA AUTOMATIZZARE CON JSON
const double E_ON = m_data.electrical_parameters.E_ON;
const double p_amb = m_data.electrical_parameters.stratosphere ? 5474 : 101325; 
const double T = m_data.electrical_parameters.stratosphere ? 217. : 303.; // [K] fluid temperature
const double delta = p_amb/101325*298/T;                                       
const double eps = 1.; // wire surface roughness correction coefficient
const double Ep = E_ON*delta*eps*(1+0.308/std::sqrt(Re*1.e+2*delta));

const double Vi = Ve - Ep*Re*log(Ep/E_ON); // [V] voltage on ionization region boundary


dof_handler.distribute_dofs(fe);

// INDEX SETS INITIALIZATION
locally_owned_dofs = dof_handler.locally_owned_dofs();                           //local dofs
//estrae gli indici globali dei gradi di libertà (DoFs) che sono posseduti localmente dal processore
locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);    //local dofs + ghost dofs


// PETSC VECTORS DECLARATIONS 
potential.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);   //ghosted
poisson_newton_update.reinit(locally_owned_dofs, mpi_communicator);              //non-ghosted
poisson_rhs.reinit(locally_owned_dofs, mpi_communicator);                        //non-ghosted

Field_X.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);     //ghosted Electric field X
Field_Y.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);     //ghosted Electric field Y

std::unordered_map<std::string, int> stringToCase{
  {"NACA", 1},
  {"WW", 2},
  {"CYL", 3}
  };

const std::string input = m_data.simulation_specification.mesh_TAG;
auto iter = stringToCase.find(input);

// ZERO_CONSTRAINTS FOR NEWTON POISSON PROBLEM (coerente con il nostro modo di procedere)
zero_constraints_poisson.clear();
zero_constraints_poisson.reinit(locally_relevant_dofs);

//VectorTools::interpolate_boundary_values(dof_handler, 3, Functions::ZeroFunction<dim>(), zero_constraints_poisson); //emitter
//VectorTools::interpolate_boundary_values(dof_handler, 4, Functions::ZeroFunction<dim>(), zero_constraints_poisson); // collector
zero_constraints_poisson.close(); 


// NON ZERO CONSTRAINTS FOR THE INITIAL SYSTEM
constraints_poisson.clear();
constraints_poisson.reinit(locally_relevant_dofs);
                
if (iter != stringToCase.end()) {
    switch (iter->second) {

        case 3:{
            VectorTools::interpolate_boundary_values(dof_handler, 3, Functions::ConstantFunction<dim>(Ve), constraints_poisson); //emitter
            VectorTools::interpolate_boundary_values(dof_handler, 4, Functions::ZeroFunction<dim>(), constraints_poisson); // collector

            break;
        }

        case 2:{

            // VectorTools::interpolate_boundary_values(dof_handler, 0, Functions::ZeroFunction<dim>(), constraints_poisson);   // Up and down
            VectorTools::interpolate_boundary_values(dof_handler, 3, Functions::ConstantFunction<dim>(Vi), constraints_poisson);    // Emitter
            VectorTools::interpolate_boundary_values(dof_handler, 4, Functions::ZeroFunction<dim>(), constraints_poisson);   // Collector
            // VectorTools::interpolate_boundary_values(dof_handler, 1, Functions::ZeroFunction<dim>(), constraints_poisson);   // Inlet
            // VectorTools::interpolate_boundary_values(dof_handler, 2, Functions::ZeroFunction<dim>(), constraints_poisson);   // Outlet

            break;
        }


        case 1:{

            break;
        }


        default:{
            std::cout << "   This TAG does not exists\n";
            break;
        }
    }
  }

constraints_poisson.close();


// DYNAMIC SPARSITY PATTERN AND POISSON MATRICES 
DynamicSparsityPattern dsp(locally_relevant_dofs);
DoFTools::make_sparsity_pattern(dof_handler, dsp, zero_constraints_poisson, false); 
SparsityTools::distribute_sparsity_pattern(dsp,
                                            dof_handler.locally_owned_dofs(),
                                            mpi_communicator,
                                            locally_relevant_dofs);


system_matrix_poisson.clear(); //store the matrix used to solve nlpoisson
system_matrix_poisson.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

laplace_matrix_poisson.clear(); //store laplace matrix
laplace_matrix_poisson.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

mass_matrix_poisson.clear();  //store mass matrix
mass_matrix_poisson.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

density_matrix.clear();  //store density matrix in the non linear poisson system
density_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);



// INITIAL SYSTEM SETUP
DynamicSparsityPattern dsp_init(locally_relevant_dofs);
DoFTools::make_sparsity_pattern(dof_handler, dsp_init, constraints_poisson, false); 
SparsityTools::distribute_sparsity_pattern(dsp_init,
                                            dof_handler.locally_owned_dofs(),
                                            mpi_communicator,
                                            locally_relevant_dofs);

initial_matrix_poisson.clear();  //store the initial matrix to initialize a good potential guess in order to make newton converge
initial_matrix_poisson.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

initial_poisson_rhs.reinit(locally_owned_dofs, mpi_communicator);                       

}
//----------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void drift_diffusion<dim>::assemble_initial_system()
{

const QTrapezoid<dim> quadrature_formula;

initial_matrix_poisson = 0;
initial_poisson_rhs = 0;


FEValues<dim> fe_values(fe,
            quadrature_formula,
            update_values | update_gradients |
            update_quadrature_points | update_JxW_values);

const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
const unsigned int n_q_points    = quadrature_formula.size();

FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
Vector<double>     cell_rhs(dofs_per_cell);

std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

unsigned int q_point = 0, idof = 0, jdof = 0;

for (const auto &cell : dof_handler.active_cell_iterators()){

    if (cell->is_locally_owned()){
      cell_matrix = 0.;
      cell_rhs = 0.;
      fe_values.reinit(cell);

        for (q_point = 0; q_point < n_q_points; ++q_point) {

            for (idof = 0; idof < dofs_per_cell; ++idof) {

                for (jdof = 0; jdof < dofs_per_cell; ++jdof){
                    cell_matrix(idof, jdof) += fe_values.shape_grad(idof, q_point) *
                    fe_values.shape_grad(jdof, q_point) * fe_values.JxW(q_point);
                }
            }
        }

        cell->get_dof_indices(local_dof_indices);
        constraints_poisson.distribute_local_to_global(cell_matrix,
                                                       cell_rhs,       
                                                       local_dof_indices,
                                                       initial_matrix_poisson,
                                                       initial_poisson_rhs);
    }
 }

initial_matrix_poisson.compress(VectorOperation::add);
initial_poisson_rhs.compress(VectorOperation::add);

}
//--------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>  
void drift_diffusion<dim>::solve_homogeneous_poisson() // find the initial value for the potential
{
  
  PETScWrappers::MPI::Vector temp;
  temp.reinit(locally_owned_dofs, mpi_communicator);
  
  //Solve homo poisson system problem
  SolverControl sc_p(dof_handler.n_dofs(), 1e-10);     
  PETScWrappers::SparseDirectMUMPS solverMUMPS(sc_p); 
  solverMUMPS.solve(initial_matrix_poisson, temp, initial_poisson_rhs);

  constraints_poisson.distribute(temp);
  temp.compress(VectorOperation::insert);
  
  // initialize the potential with the solution of this first system
  potential = temp ;

}
//--------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void drift_diffusion<dim>::assemble_poisson_laplace_matrix()
{
const QTrapezoid<dim> quadrature_formula;

laplace_matrix_poisson = 0;

FEValues<dim> fe_values(fe,
            quadrature_formula,
            update_values | update_gradients |
            update_quadrature_points | update_JxW_values);

const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
const unsigned int n_q_points    = quadrature_formula.size();

FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

unsigned int q_point = 0, idof = 0, jdof = 0;
for (const auto &cell : dof_handler.active_cell_iterators()){
    if (cell->is_locally_owned()){
        cell_matrix = 0.;
    fe_values.reinit(cell);

    for (q_point = 0; q_point < n_q_points; ++q_point) {
    for (idof = 0; idof < dofs_per_cell; ++idof) {
    for (jdof = 0; jdof < dofs_per_cell; ++jdof)
        cell_matrix(idof, jdof) += fe_values.shape_grad(idof, q_point) *
    fe_values.shape_grad(jdof, q_point) * fe_values.JxW(q_point);
    }
}

    cell->get_dof_indices(local_dof_indices);
    zero_constraints_poisson.distribute_local_to_global(cell_matrix,       
                                            local_dof_indices,
                                            laplace_matrix_poisson);
    }
}

laplace_matrix_poisson.compress(VectorOperation::add);

}

//------------------------------------------------------------------------------------------------------------------------------

template <int dim>
void drift_diffusion<dim>::assemble_poisson_mass_matrix()
{
const QTrapezoid<dim> quadrature_formula;

mass_matrix_poisson = 0;

FEValues<dim> fe_values(fe,
            quadrature_formula,
            update_values | update_gradients |
            update_quadrature_points | update_JxW_values);

const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
const unsigned int n_q_points    = quadrature_formula.size();

FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

unsigned int q_point = 0, idof = 0, jdof = 0;
for (const auto &cell : dof_handler.active_cell_iterators()){
    if (cell->is_locally_owned())
{
    cell_matrix = 0.;

    fe_values.reinit(cell);

    for (q_point = 0; q_point < n_q_points; ++q_point) {
    for (idof = 0; idof < dofs_per_cell; ++idof) {
        for (jdof = 0; jdof < dofs_per_cell; ++jdof)
    cell_matrix(idof, jdof) += fe_values.shape_value(idof, q_point) *
        fe_values.shape_value(jdof, q_point) * fe_values.JxW(q_point);
    }
    }

    cell->get_dof_indices(local_dof_indices);
    zero_constraints_poisson.distribute_local_to_global(cell_matrix,   
                                                    local_dof_indices,
                                                    mass_matrix_poisson );
}
}

mass_matrix_poisson.compress(VectorOperation::add);

}
  
//----------------------------------------------------------------------------------------------------------------------------  
  
template <int dim>
void drift_diffusion<dim>::setup_drift_diffusion()
{ 

//fix the constants
const double E_ON  = m_data.electrical_parameters.E_ON;
const double E_ref = m_data.electrical_parameters.E_ref;
const double N_ref = m_data.electrical_parameters.N_ref;
const double N_min = m_data.electrical_parameters.N_min;

// Corona inception condition: boundary condition for ion density
Functions::FEFieldFunction<dim, dealii::PETScWrappers::MPI::Vector> solution_as_function_object(dof_handler, potential, mapping);


auto boundary_evaluator = [&] (const Point<dim> &p) //lambda function
{

  Tensor<1,dim> grad_U = solution_as_function_object.gradient(p);
  Tensor<1,dim> normal = get_emitter_normal(p);

  const double En = grad_U*normal;

  const double EXP = std::exp((-En-E_ON)/E_ref); 

  const double value =  N_ref * EXP; 

  const double n = std::max(N_min, value);

  return n;

};


std::unordered_map<std::string, int> stringToCase{
  {"NACA", 1},
  {"WW", 2},
  {"CYL", 3}
  };

const std::string input = m_data.simulation_specification.mesh_TAG;
auto iter = stringToCase.find(input);

const double N_0 = m_data.electrical_parameters.stratosphere ? 2.2e-3 : 0.5e-3; // [m^-3] ambient ion density 

// DENSITY IONS CONSTRAINTS
// BCS 1 : DANNO LO STESSO RISULTATO LE DUE BCS
/*
ion_constraints.clear();
ion_constraints.reinit(locally_relevant_dofs);
VectorTools::interpolate_boundary_values(dof_handler,3, ScalarFunctionFromFunctionObject<2>(boundary_evaluator), ion_constraints); //emitter
ion_constraints.close();
*/
// BCS 2
/*
ion_constraints.clear();
ion_constraints.reinit(locally_relevant_dofs);
VectorTools::interpolate_boundary_values(dof_handler,3, ScalarFunctionFromFunctionObject<2>(boundary_evaluator), ion_constraints); //emitter
VectorTools::interpolate_boundary_values(dof_handler,4, Functions::ConstantFunction<dim>(N_min), ion_constraints); //collector
ion_constraints.close();
*/
// BCS 3

ion_constraints.clear();
ion_constraints.reinit(locally_relevant_dofs);


if (iter != stringToCase.end()) {
    switch (iter->second) {

        case 3:{
            VectorTools::interpolate_boundary_values(dof_handler,3, ScalarFunctionFromFunctionObject<dim>(boundary_evaluator), ion_constraints); //emitter
            VectorTools::interpolate_boundary_values(dof_handler,4, ScalarFunctionFromFunctionObject<dim>(boundary_evaluator), ion_constraints); //collector

            break;
        }

        case 2:{

            // VectorTools::interpolate_boundary_values(dof_handler, 0, Functions::ZeroFunction<dim>(), ion_constraints);   // Up and down
            VectorTools::interpolate_boundary_values(dof_handler,3, ScalarFunctionFromFunctionObject<dim>(boundary_evaluator), ion_constraints);    // Emitter
            // VectorTools::interpolate_boundary_values(dof_handler,4, ScalarFunctionFromFunctionObject<dim>(boundary_evaluator), ion_constraints);  // Collector
            VectorTools::interpolate_boundary_values(dof_handler, 1, Functions::ConstantFunction<dim>(N_0), ion_constraints);   // Inlet
            // VectorTools::interpolate_boundary_values(dof_handler, 2, Functions::ZeroFunction<dim>(), ion_constraints);   // Outlet

            break;
        }


        case 1:{

            break;
        }


        default:{
            std::cout << "   This TAG does not exists\n";
            break;
        }
    }
  }


ion_constraints.close();



DynamicSparsityPattern ion_dsp(locally_relevant_dofs);
DoFTools::make_sparsity_pattern(dof_handler, ion_dsp, ion_constraints, false);
SparsityTools::distribute_sparsity_pattern(ion_dsp,
                                           dof_handler.locally_owned_dofs(),
                                           mpi_communicator,
                                           locally_relevant_dofs);


// DENSITY VECTORS AND MATRICES
eta.reinit(locally_owned_dofs, mpi_communicator);                                // non-ghosted
ion_density.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator); // ghosted
old_ion_density.reinit(locally_owned_dofs, mpi_communicator);                    // non-ghosted 
ion_rhs.reinit(locally_owned_dofs, mpi_communicator);                            // non-ghosted


ion_mass_matrix.clear();
ion_mass_matrix.reinit(locally_owned_dofs, locally_owned_dofs, ion_dsp,  mpi_communicator); 

ion_system_matrix.clear();
ion_system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, ion_dsp,  mpi_communicator); 

}
//--------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void drift_diffusion<dim>::update_ion_boundary_condition(){


double theta;
const double k_min = 0.9;
const double k_max = 1.1;

const double ion_norm = ion_density.linfty_norm();
const double old_ion_norm = old_ion_density.linfty_norm();
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

pcout<<" theta is: "<<theta<<std::endl;


PETScWrappers::MPI::Vector temp;
temp.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
temp = old_ion_density;

// Corona inception condition: boundary condition for ion density
Functions::FEFieldFunction<dim, dealii::PETScWrappers::MPI::Vector> solution_as_function_object_1(dof_handler, ion_density, mapping);
Functions::FEFieldFunction<dim, dealii::PETScWrappers::MPI::Vector> solution_as_function_object_2(dof_handler, temp, mapping);


auto boundary_evaluator = [&] (const Point<dim> &p) 
{

const double ion_value = solution_as_function_object_1.value(p);
const double old_ion_value = solution_as_function_object_2.value(p);

const double value = theta*old_ion_value +(1-theta)*ion_value;

return value;

};


std::unordered_map<std::string, int> stringToCase{
  {"NACA", 1},
  {"WW", 2},
  {"CYL", 3}
  };

const std::string input = m_data.simulation_specification.mesh_TAG;
auto iter = stringToCase.find(input);

const double N_0 = m_data.electrical_parameters.stratosphere ? 2.2e-3 : 0.5e-3; // [m^-3] ambient ion density 

// DENSITY IONS CONSTRAINTS
/*
ion_constraints.clear();
VectorTools::interpolate_boundary_values(dof_handler,
                                         3,
                                         DynamicBoundaryValues<dim>(potential,ion_density,old_ion_density,dof_handler,mapping,m_data), 
                                         ion_constraints); //emitter
ion_constraints.close();
*/

ion_constraints.clear();

if (iter != stringToCase.end()) {
    switch (iter->second) {

        case 3:{
            VectorTools::interpolate_boundary_values(dof_handler,3, ScalarFunctionFromFunctionObject<dim>(boundary_evaluator), ion_constraints); //emitter
            VectorTools::interpolate_boundary_values(dof_handler,4, ScalarFunctionFromFunctionObject<dim>(boundary_evaluator), ion_constraints); //collector

            break;
        }

        case 2:{

            // VectorTools::interpolate_boundary_values(dof_handler, 0, Functions::ZeroFunction<dim>(), ion_constraints);   // Up and down
            VectorTools::interpolate_boundary_values(dof_handler,3, ScalarFunctionFromFunctionObject<dim>(boundary_evaluator), ion_constraints);    // Emitter
            VectorTools::interpolate_boundary_values(dof_handler,4, ScalarFunctionFromFunctionObject<dim>(boundary_evaluator), ion_constraints);  // Collector
            VectorTools::interpolate_boundary_values(dof_handler, 1, Functions::ConstantFunction<dim>(N_0), ion_constraints);   // Inlet
            // VectorTools::interpolate_boundary_values(dof_handler, 2, Functions::ZeroFunction<dim>(), ion_constraints);   // Outlet

            break;
        }


        case 1:{

            break;
        }


        default:{
            std::cout << "   This TAG does not exists\n";
            break;
        }
    }
  }

ion_constraints.close();


DynamicSparsityPattern ion_dsp(locally_relevant_dofs);
DoFTools::make_sparsity_pattern(dof_handler, ion_dsp, ion_constraints, false);
SparsityTools::distribute_sparsity_pattern(ion_dsp,
                                           dof_handler.locally_owned_dofs(),
                                           mpi_communicator,
                                           locally_relevant_dofs);

ion_mass_matrix.clear();
ion_mass_matrix.reinit(locally_owned_dofs, locally_owned_dofs, ion_dsp,  mpi_communicator); 

ion_system_matrix.clear();
ion_system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, ion_dsp,  mpi_communicator); 



}
//------------------------------------------------------------------------------------------------------------------------------

template <int dim>
void drift_diffusion<dim>::assemble_drift_diffusion_mass_matrix()
{
  const QTrapezoid<dim> quadrature_formula;

  ion_mass_matrix = 0;

  FEValues<dim> fe_values(fe,
        quadrature_formula,
        update_values | update_gradients |
        update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  unsigned int q_point = 0, idof = 0, jdof = 0;
  for (const auto &cell : dof_handler.active_cell_iterators()){
    if (cell->is_locally_owned())
{
  cell_matrix = 0.;

  fe_values.reinit(cell);

  for (q_point = 0; q_point < n_q_points; ++q_point) {
    for (idof = 0; idof < dofs_per_cell; ++idof) {
      for (jdof = 0; jdof < dofs_per_cell; ++jdof){
         cell_matrix(idof, jdof) += fe_values.shape_value(idof, q_point) *
         fe_values.shape_value(jdof, q_point) * fe_values.JxW(q_point);
      }
    }
  }

  cell->get_dof_indices(local_dof_indices);
  ion_constraints.distribute_local_to_global(cell_matrix,    
                                             local_dof_indices,
                                             ion_mass_matrix);
}
  }

  ion_mass_matrix.compress(VectorOperation::add);

}
//---------------------------------------------------------------------------------------------------------------------------------------------------

template <int dim>
void drift_diffusion<dim>::assemble_nonlinear_poisson()
  { 
    // Fix the constants
    const double q0 = m_data.electrical_parameters.q0;
    const double eps_r = m_data.electrical_parameters.eps_r;
    const double eps_0 = m_data.electrical_parameters.eps_0;
    const double kB = m_data.electrical_parameters.kB;
    const double T = m_data.electrical_parameters.stratosphere ? 217. : 303.; 
    const double V_TH = kB*T/q0;

    //BUILDING POISSON SYSTEM MATRIX (for newton)

    system_matrix_poisson = 0;
    //ion_mass_matrix = 0;
    density_matrix = 0; // we store in a new matrix the  ion density

    double new_value = 0;
  
    // Generate the term:  (eta)*MASS_MAT   lumped version stored in ion_mass_matrix
    for (auto iter = locally_owned_dofs.begin(); iter != locally_owned_dofs.end(); ++iter){ 

      new_value = mass_matrix_poisson(*iter, *iter) * eta(*iter);
      //ion_mass_matrix.set(*iter,*iter,new_value);
      density_matrix.set(*iter,*iter,new_value);

    }

    //ion_mass_matrix.compress(VectorOperation::insert);
    density_matrix.compress(VectorOperation::insert);
    

    system_matrix_poisson.add(eps_r * eps_0, laplace_matrix_poisson); // SYS_MAT = SYS_MAT +  eps*A
    system_matrix_poisson.add(q0 / V_TH, density_matrix);             // SYS_MAT = SYS_MAT + q0/V_TH * (eta)*MASS_MAT

  
    // BUILDING SYSTEM RHS
    poisson_rhs = 0;

    PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator); 
    mass_matrix_poisson.vmult(temp,eta);
    poisson_rhs.add(q0,temp);
    laplace_matrix_poisson.vmult(temp,potential); 
    poisson_rhs.add(- eps_r * eps_0, temp);

  
    poisson_rhs.compress(VectorOperation::insert);

    //coerentemente con il nostro approccio non ci sono apply bcs

    //pcout << "   The L_INF norm of the poisson system RHS is: "<<poisson_system_rhs.linfty_norm() <<std::endl;
    //pcout << "   The L2 norm of the poisson system RHS is: " << poisson_system_rhs.l2_norm() << std::endl;
    //pcout << "   End of solve assemble non-linear poisson"<< std::endl<<std::endl;

  }
  //----------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  double drift_diffusion<dim>::solve_poisson()
  { 

    // Fixing costants
    const double q0 = m_data.electrical_parameters.q0;
    const double kB = m_data.electrical_parameters.kB;
    const double T = m_data.electrical_parameters.stratosphere ? 217. : 303.; 
    const double V_TH = kB*T/q0;

    //Apply zero boundary conditions to the whole linear newton poisson system
    //We apply the BCs on tags 3 (emitter) and 4 (collector)


    std::unordered_map<std::string, int> stringToCase{
      {"NACA", 1},
      {"WW", 2},
      {"CYL", 3}
      };

    const std::string input = m_data.simulation_specification.mesh_TAG;
    auto iter = stringToCase.find(input);


    if (iter != stringToCase.end()) {
    switch (iter->second) {

        case 3:{
            std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

            VectorTools::interpolate_boundary_values(mapping, dof_handler,3, Functions::ZeroFunction<dim>(), emitter_boundary_values);
            MatrixTools::apply_boundary_values(emitter_boundary_values, system_matrix_poisson, poisson_newton_update, poisson_rhs);

            VectorTools::interpolate_boundary_values(mapping, dof_handler,4, Functions::ZeroFunction<dim>(), collector_boundary_values);
            MatrixTools::apply_boundary_values(collector_boundary_values, system_matrix_poisson, poisson_newton_update, poisson_rhs);

            break;
        }

        case 2:{
            std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values, up_down_boundary_values, inlet_boundary_values, outlet_boundary_values;

            VectorTools::interpolate_boundary_values(mapping, dof_handler,3, Functions::ZeroFunction<dim>(), emitter_boundary_values);
            MatrixTools::apply_boundary_values(emitter_boundary_values, system_matrix_poisson, poisson_newton_update, poisson_rhs);

            VectorTools::interpolate_boundary_values(mapping, dof_handler,4, Functions::ZeroFunction<dim>(), collector_boundary_values);
            MatrixTools::apply_boundary_values(collector_boundary_values, system_matrix_poisson, poisson_newton_update, poisson_rhs);

            // VectorTools::interpolate_boundary_values(mapping, dof_handler,0, Functions::ZeroFunction<dim>(), up_down_boundary_values);
            // MatrixTools::apply_boundary_values(up_down_boundary_values, system_matrix_poisson, poisson_newton_update, poisson_rhs);

            // VectorTools::interpolate_boundary_values(mapping, dof_handler,1, Functions::ZeroFunction<dim>(), inlet_boundary_values);
            // MatrixTools::apply_boundary_values(inlet_boundary_values, system_matrix_poisson, poisson_newton_update, poisson_rhs);

            // VectorTools::interpolate_boundary_values(mapping, dof_handler,2, Functions::ZeroFunction<dim>(), outlet_boundary_values);
            // MatrixTools::apply_boundary_values(outlet_boundary_values, system_matrix_poisson, poisson_newton_update, poisson_rhs);

            break;
        }


        case 1:{

            break;
        }


        default:{
            std::cout << "   This TAG does not exists\n";
            break;
        }
      }
    }

    //Solve poisson system problem
    const double coeff = 1e-3;
    SolverControl sc_p(dof_handler.n_dofs(), /*1e-10*/coeff *poisson_rhs.l2_norm());     
    PETScWrappers::SparseDirectMUMPS solverMUMPS(sc_p); 
    solverMUMPS.solve(system_matrix_poisson, poisson_newton_update, poisson_rhs); // sulla matrice ci sono

    //Compute the residual before the clamping
    double residual = poisson_newton_update.linfty_norm();

    //Clamping 
    for (auto iter = locally_owned_dofs.begin(); iter != locally_owned_dofs.end(); ++iter){
  
      if (poisson_newton_update[*iter] < -V_TH) { poisson_newton_update[*iter] = -V_TH; }
      else if (poisson_newton_update[*iter] > V_TH) { poisson_newton_update[*iter] = V_TH; }

      eta[*iter] *= std::exp(-poisson_newton_update[*iter]/V_TH); // aggiorno qua le cariche, non serve update charge
      
    }
   
    
    poisson_newton_update.compress(VectorOperation::insert);
    eta.compress(VectorOperation::insert);
    
    ion_constraints.distribute(eta);  
    eta.compress(VectorOperation::insert); 

    //Update current solution
    PETScWrappers::MPI::Vector temp;
    temp.reinit(locally_owned_dofs, mpi_communicator);

    temp = potential;
    temp.add(1.0, poisson_newton_update);

    potential = temp;

    return residual;

    // pcout << "L2 norm of the current solution: " << current_solution.l2_norm() << std::endl;
    // pcout << "L_INF norm of the current solution: " << current_solution.linfty_norm() << std::endl;
    // pcout << "   End of solve poisson problem"<< std::endl<<std::endl;
  
}

//-----------------------------------------------------------------------------------------------------------------------------
template <int dim>
void drift_diffusion<dim>::solve_nonlinear_poisson(const unsigned int max_iter_newton, 
                                                   const double toll_newton){

  unsigned int counter = 0; // it keeps track of newton iteration

  double increment_norm = std::numeric_limits<double>::max(); // the increment norm is + inf 

  while(counter < max_iter_newton && increment_norm > toll_newton){

    //NB: Mass and Laplace matrices are already build
    
    assemble_nonlinear_poisson();
    increment_norm = solve_poisson();  //residual computation, 
                                       //clamping on newton update, 
                                       //BCs and update of the charges are inside this method  

    counter ++;

    if(counter == max_iter_newton){
      pcout<< "   MAX NUMBER OF NEWTON ITERATIONS REACHED!"<<std::endl;
    }

  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------

template <int dim>
void drift_diffusion<dim>::assemble_drift_diffusion_matrix()
{ 
  
  //fix the constants
  const double kB = m_data.electrical_parameters.kB;
  const double q0 = m_data.electrical_parameters.q0;
  const double mu0 = m_data.electrical_parameters.mu0;
  const double T = m_data.electrical_parameters.stratosphere ? 217. : 303.; // [K] fluid temperature
             
  const double V_TH = kB * T / q0; // [V] ion temperature                           
  const double D = mu0 * V_TH; //                                                    
  
  // ion_mass_matrix already built
  // start to build the flux matrix
  

  const unsigned int vertices_per_cell = 4;
  //FullMatrix<double> Robin(vertices_per_cell, vertices_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(vertices_per_cell);

  FullMatrix<double> cell_matrix(vertices_per_cell, vertices_per_cell);
  Vector<double>     cell_rhs(vertices_per_cell);

  const unsigned int t_size = 3;

  Vector<double> A_cell_rhs(t_size), B_cell_rhs(t_size);
  FullMatrix<double> A(t_size,t_size), B(t_size,t_size);

  std::vector<types::global_dof_index> A_local_dof_indices(t_size);
  std::vector<types::global_dof_index> B_local_dof_indices(t_size);

  
  evaluate_electric_field(); // compute Field_X and Field_Y

  //const double Vh = std::sqrt(8.*numbers::PI*kB * T / Mm * Avo / 2. / numbers::PI); // Hopf velocity
  //const double Vh = std::sqrt(kB * T / Mm * Avo / 2. / numbers::PI); // Thermal velocity


  QTrapezoid<dim-1>	face_quadrature;
  const unsigned int n_q_points = face_quadrature.size();
  FEFaceValues<dim> face_values(fe, face_quadrature, update_values | update_quadrature_points /*| update_normal_vectors*/ | update_JxW_values);


  for (const auto &cell : dof_handler.active_cell_iterators())
    {
        if (cell->is_locally_owned()){

            A = 0;
            B = 0;
            cell_matrix = 0;
            cell_rhs = 0;
            //Robin = 0;
            A_cell_rhs = 0;
            B_cell_rhs = 0;


            cell->get_dof_indices(local_dof_indices); //Recupera gli indici globali dei gradi di libertà associati alla cella corrente

            /*
            // Robin conditions at outlet and (optional) at collector
            if (cell->at_boundary()) {

                for (const auto &face : cell->face_iterators()) {

                    if (face->at_boundary() && face->boundary_id() == 2) { // Outlet

                        face_values.reinit(cell, face);

                        for (unsigned int i = 0; i < vertices_per_cell; ++i) { // "i" gira su i global dof della cella

                            const double vel_f = Vel_X(local_dof_indices[i]); // creato in setup_NS - riempito in solve_NS
                            
                            for (unsigned int q = 0; q < n_q_points; ++q) {

                                for (unsigned int j = 0; j < vertices_per_cell; ++j) {

                                    Robin(i,j) += face_values.JxW(q) * face_values.shape_value(i,q) * face_values.shape_value(j,q) * vel_f;

                                }
                            }
                        }

                    } else if (face->at_boundary() && face->boundary_id() == 4) { // Collector

                        face_values.reinit(cell,face);

                        for (unsigned int i = 0; i < vertices_per_cell; ++i) {

                            Tensor<1,dim> Evec;
                            Evec[0] = Field_X(local_dof_indices[i]);
                            Evec[1] = Field_Y(local_dof_indices[i]);

                            for (unsigned int q = 0; q < n_q_points; ++q) {

                                //const Tensor<1,dim> n = face_values.normal_vector(q); // Normal  pointing into the electrode
                                //const double En = Evec * n;

                                // Alternatively, as the field should be normal to the electrode (equipotential surface):
                                const double En = std::sqrt(Evec[0]*Evec[0] + Evec[1]*Evec[1]);

                                // CHECK
                                if (std::isnan(En) || En <0.) pcout << "WARNING! Scalar product is " << En << " in " << face->center() << std::endl;

                                for (unsigned int j = 0; j < vertices_per_cell; ++j) {

                                        Robin(i,j) += face_values.JxW(q) * face_values.shape_value(i,q) * face_values.shape_value(j,q) * (Vh + mu * En);

                                }
                            }
                        }
                    }
                }
            }
            // End Robin conditions
            */

            // Lexicographic ordering
            const Point<dim> v1 = cell->vertex(2); // top left
            const Point<dim> v2 = cell->vertex(3); // top right
            const Point<dim> v3 = cell->vertex(0); // bottom left
            const Point<dim> v4 = cell->vertex(1); // bottom right

            const double u1 = -potential[local_dof_indices[2]]/V_TH;
            const double u2 = -potential[local_dof_indices[3]]/V_TH;
            const double u3 = -potential[local_dof_indices[0]]/V_TH;
            const double u4 = -potential[local_dof_indices[1]]/V_TH;

            const double l_12 = side_length(v1,v2);
            const double l_31 = side_length(v1,v3);
            const double l_24 = side_length(v4,v2);
            const double l_43 = side_length(v3,v4);

            const double l_alpha = std::sqrt(l_12*l_12 + l_24*l_24 - 2*((v1 - v2) * (v4 - v2)) );
            const double l_beta = std::sqrt(l_43*l_43 + l_24*l_24 - 2*((v2 - v4) * (v3 - v4)) );
            
            /*
            Tensor<1,dim> u_f_1, u_f_2, u_f_3, u_f_4;
            u_f_1[0] = Vel_X(local_dof_indices[2]);
            u_f_1[1] = Vel_Y(local_dof_indices[2]);
            u_f_2[0] = Vel_X(local_dof_indices[3]);
            u_f_2[1] = Vel_Y(local_dof_indices[3]);
            u_f_3[0] = Vel_X(local_dof_indices[0]);
            u_f_3[1] = Vel_Y(local_dof_indices[0]);
            u_f_4[0] = Vel_X(local_dof_indices[1]);
            u_f_4[1] = Vel_Y(local_dof_indices[1]);*/

            const Tensor<1,dim> dir_21 = (v1 - v2)/l_12;
            const Tensor<1,dim> dir_42 = (v2 - v4)/l_24;
            const Tensor<1,dim> dir_34 = (v4 - v3)/l_43;
            const Tensor<1,dim> dir_13 = (v3 - v1)/l_31;

            const double alpha21 = /*(u_f_2 * dir_21)/D*l_12*/ + (u1 - u2);
            const double alpha42 = /*(u_f_4 * dir_42)/D*l_24*/ + (u2 - u4);
            const double alpha34 = /*(u_f_3 * dir_34)/D*l_43*/ + (u4 - u3);
            const double alpha13 = /*(u_f_1 * dir_13)/D*l_31*/ + (u3 - u1);

            if (l_alpha >= l_beta) { // l_alpha is the longest diagonal: split by beta

                        const double l_23 = side_length(v2,v3);
                        const Tensor<1,dim> dir_23 = (v3 - v2)/l_beta;

                        const double alpha23 = /*(u_f_2 * dir_23)/D*l_23*/ + (u3 - u2);

                        // Triangle A:
                        A = compute_triangle_matrix(v2,v1,v3, alpha21, alpha13, -alpha23, D);

                        // Triangle B:
                        B = compute_triangle_matrix(v3,v4,v2, alpha34, alpha42, alpha23, D);

                        // Matrix assemble
                        A_local_dof_indices[0] = local_dof_indices[3];
                        A_local_dof_indices[1] = local_dof_indices[2];
                        A_local_dof_indices[2] = local_dof_indices[0];

                        B_local_dof_indices[0] = local_dof_indices[0];
                        B_local_dof_indices[1] = local_dof_indices[1];
                        B_local_dof_indices[2] = local_dof_indices[3];


                    } else { // l_beta is the longest diagonal: split by alpha
                        const double l_14 = side_length(v1,v4);
                        const Tensor<1,dim> dir_14 = (v4 - v1)/l_alpha;

                        const double alpha14 = /*(u_f_1 * dir_14)/D*l_14*/ + (u4 - u1);

                        // Triangle A:
                        A = compute_triangle_matrix(v4,v2,v1, alpha42, alpha21, alpha14, D);

                        // Triangle B:
                        B = compute_triangle_matrix(v1,v3,v4, alpha13, alpha34, -alpha14, D);

                        A_local_dof_indices[0] = local_dof_indices[1];
                        A_local_dof_indices[1] = local_dof_indices[3];
                        A_local_dof_indices[2] = local_dof_indices[2];

                        B_local_dof_indices[0] = local_dof_indices[2];
                        B_local_dof_indices[1] = local_dof_indices[0];
                        B_local_dof_indices[2] = local_dof_indices[1];
                    }

                    // As the ion system matrix is M + delta t DD, the contributions are multiplied by the timestep. M is already there!
                    for (unsigned int i = 0; i < t_size; ++i) {
                        for (unsigned int j = 0; j < t_size; ++j) {
                            A(i,j) = A(i,j)*timestep;
                            B(i,j) = B(i,j)*timestep;
                        }
                    }

                    /*
                    for (unsigned int i = 0; i < vertices_per_cell; ++i) {
                        for (unsigned int j = 0; j < vertices_per_cell; ++j) {
                            Robin(i,j) = Robin(i,j)*timestep;
                        }
                    }*/
                    
                    
                    ion_constraints.distribute_local_to_global(A, A_cell_rhs,  A_local_dof_indices, ion_system_matrix, ion_rhs);
                    ion_constraints.distribute_local_to_global(B, B_cell_rhs,  B_local_dof_indices, ion_system_matrix, ion_rhs);
                    //ion_constraints.distribute_local_to_global(Robin,  cell_rhs, local_dof_indices, ion_system_matrix, ion_rhs);

        }
	 }
   
   ion_system_matrix.compress(VectorOperation::add);
   ion_rhs.compress(VectorOperation::add);
   
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void drift_diffusion<dim>::solve_drift_diffusion()
{ 
  const double coeff = 1e-3;
  PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator);
  SolverControl sc_ion(dof_handler.n_dofs(), /*1e-10*/coeff *ion_rhs.l2_norm());
  PETScWrappers::SparseDirectMUMPS solverMUMPS_ion(sc_ion); 
  solverMUMPS_ion.solve(ion_system_matrix, temp, ion_rhs);

  ion_constraints.distribute(temp); // come il nostro DD ultimo
  temp.compress(VectorOperation::insert);

  ion_density = temp;

}
//--------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void drift_diffusion<dim>::perform_drift_diffusion_fixed_point_iteration_step() // method used to update ion_density
{  
    
  //sistema da risolvere: (MASS  +  DELTA_T * DD_MATRIX)*ION_DENS = MASS*OLD_ION_DENS   (q0 del paper???) ???


  // Fix the constants
  const double q0 = m_data.electrical_parameters.q0;

  PETScWrappers::MPI::Vector temp;
  temp.reinit(locally_owned_dofs, mpi_communicator);

  ion_rhs = 0;            // mettiamo a zero il rhs
  ion_system_matrix = 0;  // rimettiamo a zero la matrice DD ( flussi + massa )
    
  // generiamo RHS: MASS_ION * OLD_DENSITY 
  ion_mass_matrix.vmult(temp, old_ion_density);
  ion_rhs.add(1.0,temp);
  
  // generiamo LHS: (MASS + DELTA_T * DD_MATRIX)
  
  //ion_system_matrix.copy_from(ion_mass_matrix);  
  ion_system_matrix.add(1.0, ion_mass_matrix); // inizio mettendo il contributo di massa
  assemble_drift_diffusion_matrix();           // poi aggiungiamo quello di flusso

    // Uncomment to add a non-zero forcing term to DD equations ...
	/*
    Vector<double> forcing_terms(old_ion_density.size());
    VectorTools::create_right_hand_side(dof_handler, QTrapezoid<dim>(), Functions::ZeroFunction<dim>(), tmp); // ... by changing the ZeroFunction to an appropriate one
	forcing_terms = tmp;
	forcing_terms *= timestep;
	ion_rhs += forcing_terms;
	*/

  solve_drift_diffusion();

}
//------------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void drift_diffusion<dim>::evaluate_electric_field()
{
  // Inizializzazione dei vettori per i campi elettrici
  // Field_X.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);      
  // Field_Y.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  // Uso di un vettore MPI per global_dof_hits
  PETScWrappers::MPI::Vector global_dof_hits(locally_owned_dofs, mpi_communicator); 

  PETScWrappers::MPI::Vector el_field_X(locally_owned_dofs, mpi_communicator);
  PETScWrappers::MPI::Vector el_field_Y(locally_owned_dofs, mpi_communicator);

  QTrapezoid<dim-1> iv_quadrature;
  FEInterfaceValues<dim> fe_iv(fe, iv_quadrature, update_gradients);

  const unsigned int n_q_points = iv_quadrature.size();
  std::vector<Tensor<1,dim>> iv_gradients(n_q_points);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // Loop sulle celle attive
  for (auto &cell : dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      for (const auto face_index : GeometryInfo<dim>::face_indices())
      {
        fe_iv.reinit(cell, face_index);
        local_dof_indices = fe_iv.get_interface_dof_indices();
        fe_iv.get_average_of_function_gradients(potential, iv_gradients);

        for (const auto q : fe_iv.quadrature_point_indices()) 
        {
          for (const auto i : fe_iv.dof_indices()) 
          {
            // Incrementiamo global_dof_hits per tenere traccia del numero di contributi
            global_dof_hits[local_dof_indices[i]] += 1.0;

            for (unsigned int d = 0; d < dim; ++d) 
            {
              if (d == 0)
                el_field_X(local_dof_indices[i]) += -iv_gradients[q][d]; // Campo elettrico in X
              else if (d == 1)
                el_field_Y(local_dof_indices[i]) += -iv_gradients[q][d]; // Campo elettrico in Y
              else
                Assert(false, ExcNotImplemented());
            }
          }
        }
      }
    }
  }

  // Compressione parallela (somma dei valori)
  el_field_X.compress(VectorOperation::add);
  el_field_Y.compress(VectorOperation::add);
  global_dof_hits.compress(VectorOperation::add);

  // Aggiornamento delle ghost cells (sincronizzazione tra processori)
  el_field_X.update_ghost_values();
  el_field_Y.update_ghost_values();
  global_dof_hits.update_ghost_values();

  // Divisione del campo elettrico per il numero di contributi in ogni DOF
  for (auto iter = locally_owned_dofs.begin(); iter != locally_owned_dofs.end(); ++iter)
  {
    // Estrai il valore dal vettore MPI prima di utilizzarlo in std::max
    const double hit_count = global_dof_hits[*iter];

    el_field_X[*iter] /= std::max(1.0, hit_count);  // Dividi per hit_count, con protezione su 1.0
    el_field_Y[*iter] /= std::max(1.0, hit_count);  // Dividi per hit_count, con protezione su 1.0
  }

  // Compressione con inserimento per garantire che i valori siano consistenti sui processori
  el_field_X.compress(VectorOperation::insert);
  el_field_Y.compress(VectorOperation::insert);

  // Assegno i valori finali ai campi elettrici
  Field_X = el_field_X;
  Field_Y = el_field_Y;
}



//------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void drift_diffusion<dim>::run()
{

  // fix the constants
  const double N_0 = m_data.electrical_parameters.stratosphere ? 2.2e-3 : 0.5e-3; // [m^-3] ambient ion density 
  
  pcout << "   SETUP POISSON PROBLEM ... ";
  setup_poisson();
  pcout << "   Done !"<<std::endl;


  pcout << "   ASSEMBLE INITIAL SYSTEM ... ";
  assemble_initial_system();
  pcout << "   Done !"<<std::endl;


  pcout << "   SOLVE INITIAL SYSTEM (INITIALIZE POTENTIAL) ... ";
  solve_homogeneous_poisson();
  pcout << "   Done !"<<std::endl;
  

  pcout << "   ASSEMBLE MASS AND LAPLACE POISSON MATRICES ... ";
  assemble_poisson_laplace_matrix();
  assemble_poisson_mass_matrix();
  pcout << "   Done !"<<std::endl;


  pcout << "   SETUP DRIFT-DIFFUSION PROBLEM ... ";
  setup_drift_diffusion(); 
  pcout << "   Done !"<<std::endl;
  

  pcout << "   ASSEMBLE MASS DRIFT DIFFUSION MATRIX ... ";
  assemble_drift_diffusion_mass_matrix();
  pcout << "   Done !"<<std::endl;


  pcout << "   INITIALIZE OLD_ION_DENSITY ... ";
  VectorTools::interpolate(mapping, dof_handler , Functions::ConstantFunction<dim>(N_0), old_ion_density); 
  pcout << "   Done !"<<std::endl;
  

  pcout << "   STORE INITIAL CONDITIONS ... ";
  output_results(0);
  pcout << "   Done !"<<std::endl<<std::endl;


  // SET ERRORS AND TOLERANCES
  int step_number = 0;   // time evolution step

  const unsigned int max_it = 1e+3;    // riferito al gummel algorithm
  const unsigned int max_steps = 500;  // riferito al time loop

  const double time_tol = 1.0e-3;    // riferito al time loop
  const double tol = 1.e-9;         // tolleranza poisson newton

  double time_err = 1. + time_tol; // error time loop

  eta = old_ion_density; // basically eta = N_0 function

  // START THE ALGORITHM
pcout << "   START Time dependent drift_diffusion Loop ... "<< std::endl<<std::endl;

while (step_number < max_steps && time_err > time_tol){
    

  ++step_number;

  // Faster time-stepping (adaptive time-stepping would be MUCH better!)
  if (step_number % 40 == 1 && step_number > 1 && timestep < 1.e-3)
    timestep*= 10.;

  int it = 0;            // internal gummel iterator
  const double gummel_tol = 6.e-4;  // riferito al gummel algorithm
  double err = gummel_tol + 1.;    // errore gummel

  PETScWrappers::MPI::Vector previous_density;
  previous_density.reinit(locally_owned_dofs, mpi_communicator);

  pcout << "   GUMMEL ALGORITHM N. "<<step_number;

  while (err > gummel_tol && it < max_it) { // in questo ciclo NON si aggiorna old_ion_density

    solve_nonlinear_poisson(max_it,tol); // UPDATE potential AND eta (loop over k)
   
    previous_density = ion_density; // save the previously computed ion_density

    perform_drift_diffusion_fixed_point_iteration_step(); // UPDATE ion_density    

    previous_density -= ion_density;

    err = previous_density.linfty_norm()/ion_density.linfty_norm();

    eta = ion_density;

    it++; 
      
  }

  if (it >= max_it){
  pcout << "WARNING! DD achieved a relative error " << err << " after " << it << " iterations" << std::endl;
  }

  pcout << "   Done !"<< std::endl;

  pcout << "   UPDATE ION BOUNDARY CONDITIONS ... ";
  update_ion_boundary_condition();
  assemble_drift_diffusion_mass_matrix();
  pcout<<"   Done !"<<std::endl; 

  previous_density = old_ion_density;
  previous_density -= ion_density;
  
  time_err = previous_density.linfty_norm()/old_ion_density.linfty_norm();
  pcout << "   Density change from previous time-step is: " << time_err*100. << " %" << std::endl;
  pcout << "   time error is:  "<<time_err << std::endl<<std::endl;

  old_ion_density = ion_density;

  output_results(step_number);

 }
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void drift_diffusion<dim>::output_results(const unsigned int cycle)
{

DataOut<dim> data_out;
data_out.attach_dof_handler(dof_handler);


// Base directory for output
std::string base_directory = "../output";

// Directory to store the results of this simulation
std::string output_directory = base_directory + "/DD_Simulation/";

// Ensure the output directory is created (if it doesn't exist)
if (!std::filesystem::exists(output_directory))
{
    std::filesystem::create_directory(output_directory);
}

data_out.add_data_vector(potential, "potential");
data_out.add_data_vector(ion_density, "Ion_Density");
data_out.add_data_vector(Field_X, "Field_X");
data_out.add_data_vector(Field_Y, "Field_Y");

Vector<float> subdomain(triangulation.n_active_cells());
for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();

data_out.add_data_vector(subdomain, "subdomain");
data_out.build_patches();

data_out.write_vtu_with_pvtu_record(output_directory, "solution", cycle, mpi_communicator, 2, 1);


}
//##################### - HELPER FUNCTION IMPLEMENTATION - #####################################################################################

  void bernoulli (double x, double &bp, double &bn)
  {
    const double xlim = 1.0e-2;
    double ax  = fabs(x);

    bp  = 0.0;
    bn  = 0.0;

    //  X=0
    if (x == 0.0)
      {
        bp = 1.0;
        bn = 1.0;
        return;
      }

    // ASYMPTOTICS
    if (ax > 80.0)
      {
        if (x > 0.0)
          {
            bp = 0.0;
            bn = x;
          }
        else
          {
            bp = -x;
            bn = 0.0;
          }
        return;
      }

    // INTERMEDIATE VALUES
    if (ax <= 80 &&  ax > xlim)
      {
        bp = x / (exp (x) - 1.0);
        bn = x + bp;
        return;
      }

    // SMALL VALUES
    if (ax <= xlim &&  ax != 0.0)
      {
        double jj = 1.0;
        double fp = 1.0;
        double fn = 1.0;
        double df = 1.0;
        double segno = 1.0;
        while (fabs (df) > 1.0e-16)
          {
            jj += 1.0;
            segno = -segno;
            df = df * x / jj;
            fp = fp + df;
            fn = fn + segno * df;
          }
        bp = 1 / fp;
        bn = 1 / fn;
        return;
      }

  };

//-------------------------------------------------------------------------------------------------------------------------------------------------

double side_length (const Point<2> a, const Point<2> b)
{
	double length = 0.;

	if (a[0] == b[0])
		length = std::abs(a[1] - b[1]);
	else if (a[1] == b[1])
		length = std::abs(a[0] - b[0]);
	else
		length = std::sqrt(a[0]*a[0] + b[0]*b[0] - 2.*a[0]*b[0] + a[1]*a[1] + b[1]*b[1] - 2.*a[1]*b[1]);

	return length;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------

double triangle_denom(const Point<2> a, const Point<2> b, const Point<2> c)
{
	const double x1 = a[0];
	const double y1 = a[1];

	const double x2 = b[0];
	const double y2 = b[1];

	const double x3 = c[0];
	const double y3 = c[1];

	const double denom = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2);

	return denom;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------

Tensor<1,2> face_normal(const Point<2> a, const Point<2> b) {

	Tensor<1,2> tangent, normal;

	tangent[0] = b[0] - a[0];
	tangent[1] = b[1] - a[1];

	normal[0] = -tangent[1];
	normal[1] = tangent[0];

	return normal;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------

FullMatrix<double> compute_triangle_matrix(const Point<2> a, const Point<2> b, const Point<2> c, const double alpha12, const double alpha23, const double alpha31, const double D)
{
	const unsigned int size = 3;
	FullMatrix<double> tria_matrix(size,size);

	tria_matrix = 0;

	const double denom = triangle_denom(a,b,c);
	const double area = 0.5*std::abs(denom);

	const Tensor<1,2> grad_psi_1 = face_normal(b,c)/denom;
	const Tensor<1,2> grad_psi_2 = face_normal(c,a)/denom;
	const Tensor<1,2> grad_psi_3 = face_normal(a,b)/denom;

	const double l_12 = grad_psi_1 * grad_psi_2;
	const double l_23 = grad_psi_2 * grad_psi_3;
	const double l_31 = grad_psi_1 * grad_psi_3;

	double bp12, bn12, bp23, bn23, bp31, bn31;

	bernoulli(alpha12,bp12,bn12);
	bernoulli(alpha23,bp23,bn23);
	bernoulli(alpha31,bp31,bn31);

	tria_matrix(0,1) = D * area * bp12 * l_12;
	tria_matrix(0,2) = D * area * bn31 * l_31;

	tria_matrix(1,0) = D * area * bn12 * l_12;
	tria_matrix(1,2) = D * area * bp23 * l_23;

	tria_matrix(2,0) = D * area * bp31 * l_31;
	tria_matrix(2,1) = D * area * bn23 * l_23;
	
	tria_matrix(0,0) = - (tria_matrix(1,0)+tria_matrix(2,0));
	tria_matrix(1,1) = - (tria_matrix(0,1)+tria_matrix(2,1));
	tria_matrix(2,2) = - (tria_matrix(0,2)+tria_matrix(1,2));

	return tria_matrix;
}
//---------------------------------------------------------------------------------------------------------------------------------------------------
Tensor<1,2> get_emitter_normal(const Point<2> a) {

	Tensor<1,2> normal;

	normal[0] = a[0];
	normal[1] = a[1];

	const double norm = std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]);

	return normal/norm;
}