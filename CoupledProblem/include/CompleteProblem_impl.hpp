template <int dim>
CompleteProblem<dim>::CompleteProblem(parallel::distributed::Triangulation<dim> &tria,  // triangulation for the mesh
                                                             const data_struct &d)  // data struct from the user
  : m_data(d)
  , mpi_communicator(MPI_COMM_WORLD)
  , pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  , triangulation(tria)
  , fe(1) //fe for poisson / DD
  , dof_handler(tria) //dof h for poisson / DD
  , mapping() // mapping for poisson / DD
  ,	viscosity(d.fluid_parameters.viscosity)
  , gamma(d.fluid_parameters.gamma) 
  , degree(1) 
  , NS_fe(FE_Q<dim>(degree+1), dim, FE_Q<dim>(degree), 1)
  , NS_dof_handler(tria)
  , volume_quad_formula(degree + 2) 
  , face_quad_formula(degree + 2)   
  , NS_mapping() 
  , timestep(1e-5) // timestep for the DD algorithm
  , timestep_NS(1e-4) // timestep for the NS algorithm
  , timer(mpi_communicator, pcout, TimerOutput::never, TimerOutput::wall_times)
{}

//---------------------------------------------------------------------------------------------------------------------------------------------------------

template <int dim>
void CompleteProblem<dim>::setup_poisson()
{

  //fix the constant
  double const Ve = m_data.electrical_parameters.Ve;

  dof_handler.distribute_dofs(fe);

  // INDEX SETS INITIALIZATION
  locally_owned_dofs = dof_handler.locally_owned_dofs();                           //local dofs
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);    //local dofs + ghost dofs


  // PETSC VECTORS DECLARATIONS 
  potential.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);   //ghosted
  poisson_newton_update.reinit(locally_owned_dofs, mpi_communicator);              //non-ghosted
  poisson_rhs.reinit(locally_owned_dofs, mpi_communicator);                        //non-ghosted

  Field_X.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);     //ghosted Electric field X
  Field_Y.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);     //ghosted Electric field Y


  // ZERO_CONSTRAINTS FOR NEWTON POISSON PROBLEM 
  zero_constraints_poisson.clear();
  zero_constraints_poisson.reinit(locally_relevant_dofs);
  zero_constraints_poisson.close(); 


  // NON ZERO CONSTRAINTS FOR THE INITIAL SYSTEM
  constraints_poisson.clear();
  constraints_poisson.reinit(locally_relevant_dofs);

  VectorTools::interpolate_boundary_values(dof_handler, 3, Functions::ConstantFunction<dim>(Ve), constraints_poisson); //emitter
  VectorTools::interpolate_boundary_values(dof_handler, 4, Functions::ZeroFunction<dim>(), constraints_poisson); // collector

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


  Vel_X.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  Vel_Y.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  pressure.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator); 

}
//----------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::assemble_initial_system()
{

  //Set the quadrature formula
  const QTrapezoid<dim> quadrature_formula;

  // Initialize the global matrix and RHS vector to zero
  initial_matrix_poisson = 0;
  initial_poisson_rhs = 0;


  // Set up FEValues to compute necessary values (gradients, JxW, etc.) for the quadrature points
  FEValues<dim> fe_values(fe,
              quadrature_formula,
              update_values | update_gradients |
              update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int n_q_points    = quadrature_formula.size();

  // Initialize local matrix and RHS vector for each cell
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  // Vector to store global DoF indices for local cell
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  unsigned int q_point = 0, idof = 0, jdof = 0;

  // Loop over all active cells in the domain
  for (const auto &cell : dof_handler.active_cell_iterators()){

      // Only process cells that are locally owned
      if (cell->is_locally_owned()){
        cell_matrix = 0.;
        cell_rhs = 0.;
        fe_values.reinit(cell);

          // Loop over quadrature points in the cell
          for (q_point = 0; q_point < n_q_points; ++q_point) {

              // Loop over degrees of freedom for the current cell
              for (idof = 0; idof < dofs_per_cell; ++idof) {

                  // Compute contributions to the local cell matrix
                  for (jdof = 0; jdof < dofs_per_cell; ++jdof){
                      cell_matrix(idof, jdof) += fe_values.shape_grad(idof, q_point) *
                      fe_values.shape_grad(jdof, q_point) * fe_values.JxW(q_point);
                  }
              }
          }

          cell->get_dof_indices(local_dof_indices);

          // Transfer local contributions to global system with constraints
          constraints_poisson.distribute_local_to_global(cell_matrix,
                                                        cell_rhs,       
                                                        local_dof_indices,
                                                        initial_matrix_poisson,
                                                        initial_poisson_rhs);
      }
  }

  // Finalize the assembly by compressing the global matrix and RHS vector
  initial_matrix_poisson.compress(VectorOperation::add);
  initial_poisson_rhs.compress(VectorOperation::add);

}
//--------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>  
void CompleteProblem<dim>::solve_homogeneous_poisson() // find the initial value for the potential
{
  
  PETScWrappers::MPI::Vector temp;
  temp.reinit(locally_owned_dofs, mpi_communicator);
  
  //Solve homogeneous poisson system problem
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
void CompleteProblem<dim>::assemble_poisson_laplace_matrix()
{
  //Set the quadrature formula
  const QTrapezoid<dim> quadrature_formula;

  laplace_matrix_poisson = 0;

  // Set up FEValues to compute necessary values (gradients, JxW, etc.) for the quadrature points
  FEValues<dim> fe_values(fe,
              quadrature_formula,
              update_values | update_gradients |
              update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int n_q_points    = quadrature_formula.size();

  // Initialize local matrix for each cell
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

      // Transfer local contributions to global system with constraints
      zero_constraints_poisson.distribute_local_to_global(cell_matrix,       
                                              local_dof_indices,
                                              laplace_matrix_poisson);
      }
  }

  laplace_matrix_poisson.compress(VectorOperation::add);
}

//------------------------------------------------------------------------------------------------------------------------------

template <int dim>
void CompleteProblem<dim>::assemble_poisson_mass_matrix()
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
void CompleteProblem<dim>::setup_drift_diffusion()
{ 

  //fix the constants
  const double E_ON  = m_data.electrical_parameters.E_ON;
  const double E_ref = m_data.electrical_parameters.E_ref;
  const double N_ref = m_data.electrical_parameters.N_ref;
  const double N_min = m_data.electrical_parameters.N_min;

  // Set up boundary condition for ion density based on potential field gradient
  Functions::FEFieldFunction<dim, dealii::PETScWrappers::MPI::Vector> solution_as_function_object(dof_handler, potential, mapping);

  // Map strings to cases for handling different geometries (NACA, WW, CYL)
  std::unordered_map<std::string, int> stringToCase{
    {"NACA", 1},
    {"WW", 2},
    {"CYL", 3}
    };

  const std::string input = m_data.simulation_specification.mesh_TAG;
  auto iter = stringToCase.find(input);

  // Lambda function to evaluate boundary ion density based on the electric field normal component for CYL case
  auto boundary_evaluator = [&] (const Point<dim> &p) //lambda function
  {

    Point<dim> center(0.0,0.0);
    Tensor<1,dim> grad_U = solution_as_function_object.gradient(p);
    Tensor<1,dim> normal = get_emitter_normal(p,center);

    const double En = grad_U*normal;

    const double EXP = std::exp((-En-E_ON)/E_ref); 

    const double value =  N_ref * EXP; 

    // Ensure ion density is above minimum threshold
    const double n = std::max(N_min, value);

    return n;

  };

//?????????????????????????????????????????????????????????????????????????????????????????????????

  // Lambda function for emitter boundary ion density with a different center position
  auto boundary_evaluator_emitter = [&] (const Point<dim> &p) //lambda function
  {

    Point<dim> center(-0.0253,0.0);
    Tensor<1,dim> grad_U = solution_as_function_object.gradient(p);
    Tensor<1,dim> normal = get_emitter_normal(p,center);

    const double En = grad_U*normal;

    const double EXP = std::exp((-En-E_ON)/E_ref); 

    const double value =  N_ref * EXP; 

    const double n = std::max(N_min, value);

    return n;
  };

  // Set the initial ambient ion density based on stratosphere condition
  const double N_0 = m_data.electrical_parameters.stratosphere ? 2.2e-3 : 0.5e-3; // [m^-3] ambient ion density 

  ion_constraints.clear();
  ion_constraints.reinit(locally_relevant_dofs);


  // Set boundary conditions based on the specified mesh_TAG
  if (iter != stringToCase.end()) {
      switch (iter->second) {

          case 3:{// CYL case
              VectorTools::interpolate_boundary_values(dof_handler,3, ScalarFunctionFromFunctionObject<dim>(boundary_evaluator), ion_constraints); // Emitter
              VectorTools::interpolate_boundary_values(dof_handler,4, ScalarFunctionFromFunctionObject<dim>(boundary_evaluator), ion_constraints); // Collector

              break;
          }

          case 2:{// WW case

              VectorTools::interpolate_boundary_values(dof_handler,3, ScalarFunctionFromFunctionObject<dim>(boundary_evaluator_emitter), ion_constraints);    // Emitter
              VectorTools::interpolate_boundary_values(dof_handler,4, Functions::ConstantFunction<dim>(1e9), ion_constraints);  // Collector
              VectorTools::interpolate_boundary_values(dof_handler, 1, Functions::ConstantFunction<dim>(N_0), ion_constraints);   // Inlet

              break;
          }


          case 1:{// NACA case (currently no specific conditions)

              break;
          }


          default:{
              std::cout << "   This TAG does not exists\n";
              break;
          }
      }
    }


  ion_constraints.close();


  // Set up sparsity pattern for the ion density matrices
  DynamicSparsityPattern ion_dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, ion_dsp, ion_constraints, false);
  SparsityTools::distribute_sparsity_pattern(ion_dsp,
                                            dof_handler.locally_owned_dofs(),
                                            mpi_communicator,
                                            locally_relevant_dofs);


  // Initialize ion density vectors and matrices
  eta.reinit(locally_owned_dofs, mpi_communicator);                                // non-ghosted
  ion_density.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator); // ghosted
  ion_density = N_0; // Imposta tutti i valori iniziali di ion_density a N_0
  old_ion_density.reinit(locally_owned_dofs, mpi_communicator); // non-ghosted      
  old_ion_density = N_0;              
  ion_rhs.reinit(locally_owned_dofs, mpi_communicator);                            // non-ghosted

  // Clear and reinitialize matrices for ion mass and system matrices
  ion_mass_matrix.clear();
  ion_mass_matrix.reinit(locally_owned_dofs, locally_owned_dofs, ion_dsp,  mpi_communicator); 

  ion_system_matrix.clear();
  ion_system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, ion_dsp,  mpi_communicator); 

}

//--------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::update_ion_boundary_condition(){

  // Initialize variables to relaxed ion density boundary condition
  double theta;
  const double k_min = 0.995;
  const double k_max = 1.005;

  // Compute L2 norms of ion densities for update condition check
  const double ion_norm = ion_density.l2_norm();
  const double old_ion_norm = old_ion_density.l2_norm();
  const double condition = ion_norm / old_ion_norm;

  // Update `theta` based on relative change in ion densities
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

  // Temporary vector for storing previous ion density state
  PETScWrappers::MPI::Vector temp;
  temp.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  temp = old_ion_density;

  // Functions to evaluate ion density based on boundary conditions
  Functions::FEFieldFunction<dim, dealii::PETScWrappers::MPI::Vector> solution_as_function_object_1(dof_handler, ion_density, mapping);
  Functions::FEFieldFunction<dim, dealii::PETScWrappers::MPI::Vector> solution_as_function_object_2(dof_handler, temp, mapping);


  // Lambda function to evaluate boundary ion density
  auto boundary_evaluator = [&] (const Point<dim> &p) 
  {
    
    const double ion_value = solution_as_function_object_1.value(p);
    const double old_ion_value = solution_as_function_object_2.value(p);

    // Use current ion density value as boundary condition
    const double value = ion_value;

    return value;

  };

  // Map string identifiers to case numbers for different mesh types
  std::unordered_map<std::string, int> stringToCase{
    {"NACA", 1},
    {"WW", 2},
    {"CYL", 3}
    };

  const std::string input = m_data.simulation_specification.mesh_TAG;
  auto iter = stringToCase.find(input);

  const double N_0 = m_data.electrical_parameters.stratosphere ? 2.2e-3 : 0.5e-3; // [m^-3] ambient ion density 

  ion_constraints.clear();

  // Apply boundary conditions based on mesh type (emitter, collector, inlet)
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

  // Generate sparsity pattern for ion density matrices with updated constraints
  DynamicSparsityPattern ion_dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, ion_dsp, ion_constraints, false);
  SparsityTools::distribute_sparsity_pattern(ion_dsp,
                                            dof_handler.locally_owned_dofs(),
                                            mpi_communicator,
                                            locally_relevant_dofs);

  // Clear and reinitialize ion mass and system matrices based on the new constraints
  ion_mass_matrix.clear();
  ion_mass_matrix.reinit(locally_owned_dofs, locally_owned_dofs, ion_dsp,  mpi_communicator); 

  ion_system_matrix.clear();
  ion_system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, ion_dsp,  mpi_communicator); 

}
//------------------------------------------------------------------------------------------------------------------------------

template <int dim>
void CompleteProblem<dim>::assemble_drift_diffusion_mass_matrix()
{
  // Define quadrature rule for integration
  const QTrapezoid<dim> quadrature_formula;

  ion_mass_matrix = 0;

  // Set up FEValues object to compute shape functions and gradients at quadrature points
  FEValues<dim> fe_values(fe,
        quadrature_formula,
        update_values | update_gradients |
        update_quadrature_points | update_JxW_values);

  // Get the number of degrees of freedom (DoFs) per cell and quadrature points
  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int n_q_points    = quadrature_formula.size();

  // Local cell matrix for storing intermediate computations
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  unsigned int q_point = 0, idof = 0, jdof = 0;
  for (const auto &cell : dof_handler.active_cell_iterators()){
    if (cell->is_locally_owned())
    {
      cell_matrix = 0.;

      // Update FEValues for the current cell
      fe_values.reinit(cell);

      // Loop over quadrature points and DoFs to compute the local mass matrix
      for (q_point = 0; q_point < n_q_points; ++q_point) {
        for (idof = 0; idof < dofs_per_cell; ++idof) {
          for (jdof = 0; jdof < dofs_per_cell; ++jdof){
            // Compute cell matrix entry using shape values at quadrature points
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

  // Compress the matrix to finalize the assembly for parallel environments
  ion_mass_matrix.compress(VectorOperation::add);

}
//---------------------------------------------------------------------------------------------------------------------------------------------------

template <int dim>
void CompleteProblem<dim>::assemble_nonlinear_poisson()
  { 
    // Fix the constants
    const double q0 = m_data.electrical_parameters.q0;   // Elementary charge
    const double eps_r = m_data.electrical_parameters.eps_r;
    const double eps_0 = m_data.electrical_parameters.eps_0;
    const double kB = m_data.electrical_parameters.kB;  // Boltzmann constant
    const double T = m_data.electrical_parameters.stratosphere ? 217. : 303.;    //Temperature
    const double V_TH = kB*T/q0;       // Thermal voltage

    // Initialize system matrix and density matrix
    system_matrix_poisson = 0;
    density_matrix = 0;   // we store in a new matrix the ion density

    double new_value = 0;
  
    // Populate density matrix with values based on the product of eta and the mass matrix's diagonal terms
    for (auto iter = locally_owned_dofs.begin(); iter != locally_owned_dofs.end(); ++iter){ 

      new_value = mass_matrix_poisson(*iter, *iter) * eta(*iter);
      //ion_mass_matrix.set(*iter,*iter,new_value);
      density_matrix.set(*iter,*iter,new_value);

    }

    density_matrix.compress(VectorOperation::insert);
    
    // Assemble the Poisson system matrix: add permittivity term and ion density term
    system_matrix_poisson.add(eps_r * eps_0, laplace_matrix_poisson); // Add the permittivity-weighted Laplace matrix
    system_matrix_poisson.add(q0 / V_TH, density_matrix);             // Add the density matrix weighted by q0/V_TH

  
    // BUILDING SYSTEM RHS
    poisson_rhs = 0;

    // Temporary vector for storing intermediate calculations
    PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator); 

    // Calculate the ion density contribution to RHS: q0 * (mass_matrix_poisson * eta)
    mass_matrix_poisson.vmult(temp,eta);
    poisson_rhs.add(q0,temp);

    // Calculate the potential contribution to RHS: -eps_r * eps_0 * (laplace_matrix_poisson * potential)
    laplace_matrix_poisson.vmult(temp,potential); 
    poisson_rhs.add(- eps_r * eps_0, temp);

  
    poisson_rhs.compress(VectorOperation::insert);

  }
  //----------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  double CompleteProblem<dim>::solve_poisson()
  { 

  // Fixing costants
  const double q0 = m_data.electrical_parameters.q0;
  const double kB = m_data.electrical_parameters.kB;
  const double T = m_data.electrical_parameters.stratosphere ? 217. : 303.; 
  const double V_TH = kB*T/q0;

  //Apply zero boundary conditions to the whole linear newton poisson system
  //We apply the BCs on tags 3 (emitter) and 4 (collector)

  // Define maps to store boundary conditions for emitter (tag 3) and collector (tag 4)
  std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;


  // Apply zero boundary condition on emitter and collector
  VectorTools::interpolate_boundary_values(mapping, dof_handler,3, Functions::ZeroFunction<dim>(), emitter_boundary_values);
  MatrixTools::apply_boundary_values(emitter_boundary_values, system_matrix_poisson, poisson_newton_update, poisson_rhs);

  VectorTools::interpolate_boundary_values(mapping, dof_handler,4, Functions::ZeroFunction<dim>(), collector_boundary_values);
  MatrixTools::apply_boundary_values(collector_boundary_values, system_matrix_poisson, poisson_newton_update, poisson_rhs);

  //Set up solver for the Poisson system
  const double coeff = 1e-2;
  SolverControl sc_p(dof_handler.n_dofs(), /*1e-10*/coeff *poisson_rhs.l2_norm());     
  PETScWrappers::SparseDirectMUMPS solverMUMPS(sc_p); 
  solverMUMPS.solve(system_matrix_poisson, poisson_newton_update, poisson_rhs); // sulla matrice ci sono

  // Calculate maximum residual (infinity norm) as a measure of convergence
  double residual = poisson_newton_update.linfty_norm();

  // Clamp the solution to stay within the range [-V_TH, V_TH]
  for (auto iter = locally_owned_dofs.begin(); iter != locally_owned_dofs.end(); ++iter){

    if (poisson_newton_update[*iter] < -V_TH) { poisson_newton_update[*iter] = -V_TH; }
    else if (poisson_newton_update[*iter] > V_TH) { poisson_newton_update[*iter] = V_TH; }

    eta[*iter] *= std::exp(-poisson_newton_update[*iter]/V_TH); // aggiorno qua le cariche, non serve update charge
    
  }
  
  // Ensure consistency of vectors across distributed MPI processes
  poisson_newton_update.compress(VectorOperation::insert);
  eta.compress(VectorOperation::insert);
  
  // Distribute constraints across `eta` to ensure charge density continuity
  ion_constraints.distribute(eta);   
  eta.compress(VectorOperation::insert); 

  // Update the potential field with the current solution from Newton's update
  PETScWrappers::MPI::Vector temp;
  temp.reinit(locally_owned_dofs, mpi_communicator);

  temp = potential;
  temp.add(1.0, poisson_newton_update);

  potential = temp;

  return residual;   // Return the residual for convergence check
}

//-----------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::solve_nonlinear_poisson(const unsigned int max_iter_newton, 
                                                   const double toll_newton){

  unsigned int counter = 0;       // Track the Newton iteration count

  double increment_norm = std::numeric_limits<double>::max();    // Initialize increment norm to a high value

  while(counter < max_iter_newton && increment_norm > toll_newton){

    assemble_nonlinear_poisson();  // Assemble the system matrix and right-hand side vector for the nonlinear Poisson problem

    // Solve the assembled Poisson system, returning the residual (norm of the Newton update)
    increment_norm = solve_poisson();  

    counter ++;

    if(counter == max_iter_newton){
      pcout<< "   MAX NUMBER OF NEWTON ITERATIONS REACHED!"<<std::endl;
    }

  }

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------

template <int dim>
void CompleteProblem<dim>::assemble_drift_diffusion_matrix()
{ 
  
  //fix the constants
  const double kB = m_data.electrical_parameters.kB;
  const double q0 = m_data.electrical_parameters.q0;
  const double mu0 = m_data.electrical_parameters.mu0;
  const double T = m_data.electrical_parameters.stratosphere ? 217. : 303.; 
             
  const double V_TH = kB * T / q0; // Ion temperature                           
  const double D = mu0 * V_TH;    // diffusion coefficient


  const unsigned int vertices_per_cell = 4;
  FullMatrix<double> Robin(vertices_per_cell, vertices_per_cell);   // Robin boundary matrix

  std::vector<types::global_dof_index> local_dof_indices(vertices_per_cell);

  FullMatrix<double> cell_matrix(vertices_per_cell, vertices_per_cell);
  Vector<double>     cell_rhs(vertices_per_cell);

  const unsigned int t_size = 3;  // number of vertices per triangle


  Vector<double> A_cell_rhs(t_size), B_cell_rhs(t_size);   // RHS vectors for triangles A and B
  FullMatrix<double> A(t_size,t_size), B(t_size,t_size);   // local matrices for triangles A and B

  std::vector<types::global_dof_index> A_local_dof_indices(t_size);
  std::vector<types::global_dof_index> B_local_dof_indices(t_size);

  
  evaluate_electric_field();   // calculate Field_X and Field_Y, electric field components

  const double p_amb = m_data.electrical_parameters.stratosphere ? 5474 : 101325; 
  const double delta = p_amb/101325*298/T;   
  const double Mm = 29.e-3; // kg m^-3,average air molar mass
  const double Avo = 6.022e+23; // Avogadro's number
  const double mu = mu0 * delta; // scaled mobility from Moseley

  const double Vh = std::sqrt(8.*numbers::PI*kB * T / Mm * Avo / 2. / numbers::PI); // Hopf velocity

  QTrapezoid<dim-1>	face_quadrature;
  const unsigned int n_q_points = face_quadrature.size();
  FEFaceValues<dim> face_values(fe, face_quadrature, update_values | update_quadrature_points /*| update_normal_vectors*/ | update_JxW_values);


  for (const auto &cell : dof_handler.active_cell_iterators())
    {
        if (cell->is_locally_owned()){

            // Initialize matrices and vectors for the current cell
            A = 0;
            B = 0;
            cell_matrix = 0;
            cell_rhs = 0;
            //Robin = 0;
            A_cell_rhs = 0;
            B_cell_rhs = 0;


            cell->get_dof_indices(local_dof_indices); //Recupera gli indici globali dei gradi di libertà associati alla cella corrente

            
            // Robin boundary conditions (outlet and collector)
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

                    // Robin condition at the collector (boundary_id = 4)
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
            

            // Lexicographic ordering of the cell vertices
            const Point<dim> v1 = cell->vertex(2); // top left
            const Point<dim> v2 = cell->vertex(3); // top right
            const Point<dim> v3 = cell->vertex(0); // bottom left
            const Point<dim> v4 = cell->vertex(1); // bottom right

            // Normalized potential at each vertex
            const double u1 = -potential[local_dof_indices[2]]/V_TH;
            const double u2 = -potential[local_dof_indices[3]]/V_TH;
            const double u3 = -potential[local_dof_indices[0]]/V_TH;
            const double u4 = -potential[local_dof_indices[1]]/V_TH;

            // Side lengths
            const double l_12 = side_length(v1,v2);
            const double l_31 = side_length(v1,v3);
            const double l_24 = side_length(v4,v2);
            const double l_43 = side_length(v3,v4);

            // Determine the longest diagonal for triangle splitting
            const double l_alpha = std::sqrt(l_12*l_12 + l_24*l_24 - 2*((v1 - v2) * (v4 - v2)) );
            const double l_beta = std::sqrt(l_43*l_43 + l_24*l_24 - 2*((v2 - v4) * (v3 - v4)) );
            
            
            Tensor<1,dim> u_f_1, u_f_2, u_f_3, u_f_4;
            u_f_1[0] = Vel_X(local_dof_indices[2]);
            u_f_1[1] = Vel_Y(local_dof_indices[2]);
            u_f_2[0] = Vel_X(local_dof_indices[3]);
            u_f_2[1] = Vel_Y(local_dof_indices[3]);
            u_f_3[0] = Vel_X(local_dof_indices[0]);
            u_f_3[1] = Vel_Y(local_dof_indices[0]);
            u_f_4[0] = Vel_X(local_dof_indices[1]);
            u_f_4[1] = Vel_Y(local_dof_indices[1]);

            const Tensor<1,dim> dir_21 = (v1 - v2)/l_12;
            const Tensor<1,dim> dir_42 = (v2 - v4)/l_24;
            const Tensor<1,dim> dir_34 = (v4 - v3)/l_43;
            const Tensor<1,dim> dir_13 = (v3 - v1)/l_31;

            const double alpha21 = (u_f_2 * dir_21)/D*l_12 + (u1 - u2);
            const double alpha42 = (u_f_4 * dir_42)/D*l_24 + (u2 - u4);
            const double alpha34 = (u_f_3 * dir_34)/D*l_43 + (u4 - u3);
            const double alpha13 = (u_f_1 * dir_13)/D*l_31 + (u3 - u1);

            // Calculate drift-diffusion coefficients for triangles A and B
            if (l_alpha >= l_beta) { // l_alpha is the longest diagonal: split by beta

              const double l_23 = side_length(v2,v3);
              const Tensor<1,dim> dir_23 = (v3 - v2)/l_beta;

              const double alpha23 = (u_f_2 * dir_23)/D*l_23 + (u3 - u2);

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

              const double alpha14 = (u_f_1 * dir_14)/D*l_14 + (u4 - u1);

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

          
            for (unsigned int i = 0; i < vertices_per_cell; ++i) {
                for (unsigned int j = 0; j < vertices_per_cell; ++j) {
                    Robin(i,j) = Robin(i,j)*timestep;
                }
            }
          
            ion_constraints.distribute_local_to_global(A, A_cell_rhs,  A_local_dof_indices, ion_system_matrix, ion_rhs);
            ion_constraints.distribute_local_to_global(B, B_cell_rhs,  B_local_dof_indices, ion_system_matrix, ion_rhs);
            ion_constraints.distribute_local_to_global(Robin,  cell_rhs, local_dof_indices, ion_system_matrix, ion_rhs);

        }
	 }
   
   ion_system_matrix.compress(VectorOperation::add);
   ion_rhs.compress(VectorOperation::add);
   
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::solve_drift_diffusion()
{ 
  // Set a coefficient for the solver control
  const double coeff = 1e-2;

  PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator);

  SolverControl sc_ion(dof_handler.n_dofs(), /*1e-10*/coeff *ion_rhs.l2_norm());
  PETScWrappers::SparseDirectMUMPS solverMUMPS_ion(sc_ion); 

  // Solve the linear system of equations: ion_system_matrix * temp = ion_rhs
  solverMUMPS_ion.solve(ion_system_matrix, temp, ion_rhs);

  ion_constraints.distribute(temp); // come il nostro DD ultimo
  temp.compress(VectorOperation::insert);

  ion_density = temp;

}
//--------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::perform_drift_diffusion_fixed_point_iteration_step() // method used to update ion_density
{  

  // Fix the constants
  const double q0 = m_data.electrical_parameters.q0;

  // Create a temporary vector to hold the right-hand side values for the linear system
  PETScWrappers::MPI::Vector temp;
  temp.reinit(locally_owned_dofs, mpi_communicator);

  ion_rhs = 0;           
  ion_system_matrix = 0;  // Reset the drift-diffusion system matrix to zero (including mass and flux contributions)
    
  // Generate rhs: mass_ion * old_density
  ion_mass_matrix.vmult(temp, old_ion_density);   // Multiply the ion mass matrix by the old ion density
  ion_rhs.add(1.0,temp);                          // Add the result to the RHS vector, effectively storing the contribution from the mass term
  
  // Generate lhs
  ion_system_matrix.add(1.0, ion_mass_matrix); // Start by adding the mass contribution to the system matrix
  assemble_drift_diffusion_matrix();           // Assemble the drift-diffusion matrix to include flux contributions

  solve_drift_diffusion();     // Solve the drift-diffusion equation with the updated system matrix and RHS

}
//------------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::evaluate_electric_field()
{

  // Get the number of degrees of freedom (DOFs) per cell for the finite element
  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  // Use an MPI vector for global_dof_hits to keep track of contributions from each DOF
  PETScWrappers::MPI::Vector global_dof_hits(locally_owned_dofs, mpi_communicator); 

  // Vectors to hold the electric field components, initialized for MPI parallelism
  PETScWrappers::MPI::Vector el_field_X(locally_owned_dofs, mpi_communicator);
  PETScWrappers::MPI::Vector el_field_Y(locally_owned_dofs, mpi_communicator);

  // Set up quadrature for evaluating interface values
  QTrapezoid<dim-1> iv_quadrature;
  FEInterfaceValues<dim> fe_iv(fe, iv_quadrature, update_gradients);

  const unsigned int n_q_points = iv_quadrature.size();   // Number of quadrature points for integration
  std::vector<Tensor<1,dim>> iv_gradients(n_q_points);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);   // Local DOF indices for the current cell

  // Loop over active cells in the mesh
  for (auto &cell : dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      for (const auto face_index : GeometryInfo<dim>::face_indices())
      {
        fe_iv.reinit(cell, face_index);    // Reinitialize interface values for the current cell and face
        local_dof_indices = fe_iv.get_interface_dof_indices();   // Get DOF indices at the interface
        fe_iv.get_average_of_function_gradients(potential, iv_gradients);    // Calculate average gradients of the potential function

        // Loop over quadrature points to compute the electric field
        for (const auto q : fe_iv.quadrature_point_indices()) 
        {
          for (const auto i : fe_iv.dof_indices()) 
          {
            // Increment global_dof_hits to track the number of contributions for each DOF
            global_dof_hits[local_dof_indices[i]] += 1.0;

            for (unsigned int d = 0; d < dim; ++d) 
            {
              if (d == 0)
                el_field_X(local_dof_indices[i]) += -iv_gradients[q][d]; // Electric field in the X direction
              else if (d == 1)
                el_field_Y(local_dof_indices[i]) += -iv_gradients[q][d]; // Electric field in the Y direction
              else
                Assert(false, ExcNotImplemented());    // Ensure that dimensionality is correctly handled
            }
          }
        }
      }
    }
  }

  
  el_field_X.compress(VectorOperation::add);
  el_field_Y.compress(VectorOperation::add);
  global_dof_hits.compress(VectorOperation::add);

  // Update ghost cells to synchronize values across different processes
  el_field_X.update_ghost_values();
  el_field_Y.update_ghost_values();
  global_dof_hits.update_ghost_values();

  // Divide the electric field components by the number of contributions for each DOF
  for (auto iter = locally_owned_dofs.begin(); iter != locally_owned_dofs.end(); ++iter)
  {
    // Extract hit count from the MPI vector before using it in std::max
    const double hit_count = global_dof_hits[*iter];

    el_field_X[*iter] /= std::max(1.0, hit_count);  // Protect against division by zero
    el_field_Y[*iter] /= std::max(1.0, hit_count);  
  }

  // Final compression with insertion to ensure consistency of values across processors
  el_field_X.compress(VectorOperation::insert);
  el_field_Y.compress(VectorOperation::insert);

  // Assign the final computed electric field values to class members
  Field_X = el_field_X;
  Field_Y = el_field_Y;
}

//##################### - Navier-Stokes - #####################################################################################

template <int dim>
void CompleteProblem<dim>::setup_NS()
  {
   
    // SET UP DOFS
    NS_dof_handler.distribute_dofs(NS_fe);

    // Create a mask for reordering DOFs, where the last component corresponds to pressure
    std::vector<unsigned int> block_component(dim + 1, 0);    
    block_component[dim] = 1;      // Set the last block for pressure
    DoFRenumbering::component_wise(NS_dof_handler, block_component);    // Renumber DOFs component-wise

    // Count the degrees of freedom for each block (velocity and pressure)
    dofs_per_block = DoFTools::count_dofs_per_fe_block(NS_dof_handler, block_component);

    // Partitioning for owned DOFs
    unsigned int dof_u = dofs_per_block[0];
    unsigned int dof_p = dofs_per_block[1];

    // Resize partitioning vectors for owned DOFs
	  owned_partitioning.resize(2);
    owned_partitioning[0] = NS_dof_handler.locally_owned_dofs().get_view(0, dof_u);      //Extract the set of locally owned DoF indices for each component within the mask that are owned by the current processor.
    owned_partitioning[1] = NS_dof_handler.locally_owned_dofs().get_view(dof_u, dof_u + dof_p);

    // Extract locally relevant DOFs, including ghost cells
    DoFTools::extract_locally_relevant_dofs(NS_dof_handler,NS_locally_relevant_dofs );     //Extract the set of global DoF indices that are active on the current DoFHandler. This is the union of DoFHandler::locally_owned_dofs() and the DoF indices on all ghost cells.

    // Resize relevant partitioning vectors
    relevant_partitioning.resize(2);
    relevant_partitioning[0] = NS_locally_relevant_dofs.get_view(0, dof_u);
    relevant_partitioning[1] = NS_locally_relevant_dofs.get_view(dof_u, dof_u + dof_p);

   
    //MAKE CONSTRAINTS
    const FEValuesExtractors::Vector velocities(0);           // Extractor for velocity components
    const FEValuesExtractors::Scalar vertical_velocity(1);    // Extractor for vertical velocity
    const FEValuesExtractors::Vector vertical_velocity_and_pressure(1);    // Extractor for both vertical velocity and pressure

    // Map to link simulation types with their identifiers
    std::unordered_map<std::string, int> stringToCase{
      {"NACA", 1},
      {"WW", 2},
      {"CYL", 3}
      };

    const std::string input = m_data.simulation_specification.mesh_TAG;   // Get mesh type from simulation specification
    auto iter = stringToCase.find(input);      // Find corresponding case in the map

    nonzero_NS_constraints.clear();      // Clear previous constraints
    nonzero_NS_constraints.reinit(NS_locally_relevant_dofs);     // Reinitialize constraints based on relevant DOFs

    // Create hanging node constraints to ensure continuity
    DoFTools::make_hanging_node_constraints(NS_dof_handler, nonzero_NS_constraints);  //Lo lasciamo nel caso si faccia refine mesh adattivo

    // Set boundary conditions based on the selected mesh type
    if (iter != stringToCase.end()) {
    switch (iter->second) {

        case 3:{

            std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

            // Interpolate boundary values for the cylindrical case
            VectorTools::interpolate_boundary_values(NS_dof_handler, 3, Functions::ZeroFunction<dim>(dim+1), nonzero_NS_constraints, NS_fe.component_mask(velocities));
            VectorTools::interpolate_boundary_values(NS_dof_handler, 4, Functions::ZeroFunction<dim>(dim+1), nonzero_NS_constraints, NS_fe.component_mask(velocities));
            VectorTools::interpolate_boundary_values(NS_dof_handler, 4, Functions::ZeroFunction<dim>(dim+1), nonzero_NS_constraints, NS_fe.component_mask(vertical_velocity_and_pressure)); 
            break;
        }

        case 2:{
              VectorTools::interpolate_boundary_values(NS_dof_handler,                           
                                            0,                                      // Up and down 
                                            Functions::ZeroFunction<dim>(dim+1),    // For each degree of freedom at the boundary, its boundary value will be overwritten if its index already exists in boundary_values. Otherwise, a new entry with proper index and boundary value for this degree of freedom will be inserted into boundary_values.
                                            nonzero_NS_constraints,
                                            NS_fe.component_mask(vertical_velocity));


              VectorTools::interpolate_boundary_values(NS_dof_handler,
                                                      3,                                      //Emitter                  
                                                      Functions::ZeroFunction<dim>(dim+1),
                                                      nonzero_NS_constraints,
                                                      NS_fe.component_mask(velocities));
                                                      
              VectorTools::interpolate_boundary_values(NS_dof_handler,
                                                      4,                                      //Collector                                     
                                                      Functions::ZeroFunction<dim>(dim+1),
                                                      nonzero_NS_constraints,
                                                      NS_fe.component_mask(velocities));
                                                      
              VectorTools::interpolate_boundary_values(NS_dof_handler, 
                                                      1,                    // Inlet
                                                      BoundaryValues<dim>(1.0), // Functions::ZeroFunction<dim>(dim+1), 
                                                      nonzero_NS_constraints, 
                                                      NS_fe.component_mask(velocities));
                                                      
            
              VectorTools::interpolate_boundary_values(NS_dof_handler,
                                                      2,                    // Outlet
                                                      Functions::ZeroFunction<dim>(dim+1),//BoundaryValues<dim>()
                                                      nonzero_NS_constraints,
                                                      NS_fe.component_mask(vertical_velocity_and_pressure));  //Vertical velocity and pressure at outlet equal to 0

            break;
        }


        case 1:{
              // Additional conditions can be added here for the NACA case
            break;
        }


        default:{
            std::cout << "   This TAG does not exists\n";
            break;
        }
      }
    }                         
    
    nonzero_NS_constraints.close();
    
    // Prepare for zero constraints (e.g., for Dirichlet boundary conditions)
    zero_NS_constraints.clear();               
	  zero_NS_constraints.reinit(NS_locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(NS_dof_handler, zero_NS_constraints);

    // Set zero boundary conditions based on the selected mesh type
    if (iter != stringToCase.end()) {
      switch (iter->second) {

        case 3:{
            std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

            VectorTools::interpolate_boundary_values(NS_dof_handler, 3, Functions::ZeroFunction<dim>(dim+1), zero_NS_constraints, NS_fe.component_mask(velocities));
            VectorTools::interpolate_boundary_values(NS_dof_handler, 4, Functions::ZeroFunction<dim>(dim+1), zero_NS_constraints, NS_fe.component_mask(velocities));
            VectorTools::interpolate_boundary_values(NS_dof_handler, 4, Functions::ZeroFunction<dim>(dim+1), zero_NS_constraints, NS_fe.component_mask(vertical_velocity_and_pressure)); 

            break;
        }

        case 2:{
            VectorTools::interpolate_boundary_values(NS_dof_handler,
                                        0,
                                        Functions::ZeroFunction<dim>(dim+1),
                                        zero_NS_constraints,
                                        NS_fe.component_mask(vertical_velocity));
            VectorTools::interpolate_boundary_values(NS_dof_handler,
                                                    3,
                                                    Functions::ZeroFunction<dim>(dim+1),
                                                    zero_NS_constraints,
                                                    NS_fe.component_mask(velocities));
            VectorTools::interpolate_boundary_values(NS_dof_handler,
                                                    4,
                                                    Functions::ZeroFunction<dim>(dim+1),
                                                    zero_NS_constraints,
                                                    NS_fe.component_mask(velocities));
            VectorTools::interpolate_boundary_values(NS_dof_handler,
                                                    1,
                                                    Functions::ZeroFunction<dim>(dim+1),
                                                    zero_NS_constraints,
                                                    NS_fe.component_mask(velocities));
            VectorTools::interpolate_boundary_values(NS_dof_handler,
                                                    2, 
                                                    Functions::ZeroFunction<dim>(dim+1),
                                                    zero_NS_constraints,
                                                    NS_fe.component_mask(vertical_velocity_and_pressure));
            break;
        }


        case 1:{
              // Additional conditions can be added here for the NACA case
            break;
        }


        default:{
            std::cout << "   This TAG does not exists\n";
            break;
        }
      }
    }   
    
    zero_NS_constraints.close();

    // Output information about the problem setup
    pcout << "   Number of active fluid cells: "
      << triangulation.n_global_active_cells() << std::endl
      << "   Number of degrees of freedom: " << NS_dof_handler.n_dofs() << " ("
      << dof_u << '+' << dof_p << ')' << std::endl;

    
    //INITIALIZE SYSTEM
    BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
    DoFTools::make_sparsity_pattern(NS_dof_handler, dsp, nonzero_NS_constraints);
    NS_sparsity_pattern.copy_from(dsp);
    SparsityTools::distribute_sparsity_pattern(dsp, NS_dof_handler.locally_owned_dofs(), mpi_communicator, NS_locally_relevant_dofs); 

   
    preconditioner.reset();
    NS_system_matrix.clear();
    NS_mass_matrix.clear();
    pressure_mass_matrix.clear(); 

   
    NS_system_matrix.reinit(owned_partitioning, dsp, mpi_communicator);
    NS_solution.reinit(owned_partitioning, relevant_partitioning, mpi_communicator);
    NS_solution_update.reinit(owned_partitioning, mpi_communicator);
    NS_system_rhs.reinit(owned_partitioning, mpi_communicator);


    NS_mass_matrix.reinit(owned_partitioning, dsp, mpi_communicator);

    // Only the $(1, 1)$ block in the mass schur matrix is used. (B M_u B.T)
    // Compute the sparsity pattern for mass schur in advance.
    // The only nonzero block has the same sparsity pattern as $BB^T$.
    BlockDynamicSparsityPattern schur_dsp(dofs_per_block, dofs_per_block);
    schur_dsp.block(1, 1).compute_mmult_pattern(NS_sparsity_pattern.block(1, 0), NS_sparsity_pattern.block(0, 1));
    pressure_mass_matrix.reinit(owned_partitioning, schur_dsp, mpi_communicator);       //We want to reinitialize our matrix with new partitions of data

    //Initialize the NS_solution 
    PETScWrappers::MPI::BlockVector tmp1;
    tmp1.reinit(owned_partitioning, mpi_communicator);
    tmp1 = NS_solution;
    VectorTools::interpolate(NS_mapping, NS_dof_handler, Functions::ZeroFunction<dim>(dim+1), tmp1);   //Zero initial solution
    NS_solution = tmp1;
      
    owned_partitioning_p = NS_dof_handler.locally_owned_dofs().get_view(dof_u, dof_u + dof_p);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::assemble_NS(bool use_nonzero_constraints,
                                       bool assemble_system)
{ 
  
  TimerOutput::Scope timer_section(timer, "Assemble system");        //Enter the given section in the timer

  if (assemble_system)           
      NS_system_matrix = 0;
      NS_mass_matrix = 0;
  }

  NS_system_rhs = 0;

  // Initialize finite element values for volume quadrature
  FEValues<dim> fe_values(NS_fe,                                                         
                          volume_quad_formula,             
                          update_values | update_quadrature_points |
                          update_JxW_values | update_gradients);
  
  // Initialize finite element values for face quadrature
  FEFaceValues<dim> fe_face_values(NS_fe,
                                  face_quad_formula,
                                  update_values | update_normal_vectors |
                                  update_quadrature_points |
                                  update_JxW_values); 

  // Define degrees of freedom and number of quadrature points
  const unsigned int dofs_per_cell = NS_fe.dofs_per_cell;     // Number of DoFs per cell (22 for NS)
  const unsigned int n_q_points = volume_quad_formula.size(); // Number of quadrature points

  const FEValuesExtractors::Vector velocities(0);   //Extractor cells velocities that takes the vector in position 0
  const FEValuesExtractors::Scalar pressure(dim);   //Extractor cells pressure that takes the scalar in position dim

  FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> local_mass_matrix(dofs_per_cell, dofs_per_cell);

  Vector<double> local_rhs(dofs_per_cell);

  // Local DoF indices and temporary storage for velocity and pressure values
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  std::vector<Tensor<1, dim>> current_velocity_values(n_q_points);                  
  std::vector<Tensor<2, dim>> current_velocity_gradients(n_q_points);
  std::vector<double> current_velocity_divergences(n_q_points);
  std::vector<double> current_pressure_values(n_q_points);

  // Prepare vectors for shape function values
  std::vector<double> div_phi_u(dofs_per_cell);
  std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
  std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
  std::vector<double> phi_p(dofs_per_cell);

  Tensor<1,dim> f;   // Force vector
  
  // Prepare for handling electrical variables
  auto ion_cell = dof_handler.begin_active();  
  const auto ion_endc = dof_handler.end();     
  
  const unsigned int ion_dofs_per_cell = 4;    // Number of DoFs for ion cell
  std::vector<types::global_dof_index> ion_local_dof_indices(ion_dofs_per_cell);   // Local DoF indices for ion cell

  const double q0 = m_data.electrical_parameters.q0;    // Electric charge
  const double rho = m_data.electrical_parameters.stratosphere ? 0.089 : 1.225; // kg m^-3
  
  // Set up quadrature for ion problem
  const QGauss<dim> quadrature_formula_ion(degree + 2);  
  FEValues<dim> fe_values_ion(fe, 
                              quadrature_formula_ion, 
                              update_values | update_quadrature_points | update_JxW_values);
  
  // Iterate over NS and ion cells simultaneously
  for (auto cell = NS_dof_handler.begin_active(), ion_cell = dof_handler.begin_active(); 
       cell != NS_dof_handler.end() && ion_cell != dof_handler.end(); 
       ++cell, ++ion_cell) 
  {
      Assert(cell->index() == ion_cell->index(), ExcMessage("Mismatch between NS and ion cells!"));

      if (cell->is_locally_owned()) {   // Proceed only if the cell is locally owned

          fe_values.reinit(cell);          // Initialize FEValues for current NS cell
          fe_values_ion.reinit(ion_cell);  // Initialize for corresponding ion cell

          if (assemble_system) {
              local_matrix = 0;
              local_mass_matrix = 0;
          } 

          local_rhs = 0;

          // Fetch velocity and pressure values from the solution vector
          fe_values[velocities].get_function_values(NS_solution, current_velocity_values);
          fe_values[velocities].get_function_gradients(NS_solution, current_velocity_gradients);
          fe_values[velocities].get_function_divergences(NS_solution, current_velocity_divergences);
          fe_values[pressure].get_function_values(NS_solution, current_pressure_values);

          // Obtain field values at quadrature points
          std::vector<double> field_x_values(n_q_points);
          std::vector<double> field_y_values(n_q_points);
          std::vector<double> ion_density_values(n_q_points);

          fe_values_ion.get_function_values(Field_X, field_x_values);
          fe_values_ion.get_function_values(Field_Y, field_y_values);
          fe_values_ion.get_function_values(ion_density, ion_density_values);

          if (ion_cell->is_locally_owned()) {
             // Get the DoF for the corresponding ion cell
            ion_cell->get_dof_indices(ion_local_dof_indices);
          }

          // Compute electric force based on fields and ion density
          for (unsigned int q = 0; q < n_q_points; ++q) {
              // Prepare electric force at quadrature points
              Tensor<1, dim> f;
              f[0] = q0 * field_x_values[q] / rho * ion_density_values[q];  // Forza in X
              f[1] = q0 * field_y_values[q] / rho * ion_density_values[q];  // Forza in Y

              // Loop over the DoFs of the NS system
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                div_phi_u[k] = fe_values[velocities].divergence(k, q);           //Returns the value of the k-th shape function in q quadrature point
                grad_phi_u[k] = fe_values[velocities].gradient(k, q);
                phi_u[k] = fe_values[velocities].value(k, q);
                phi_p[k] = fe_values[pressure].value(k, q); 
              }

              // Assemble contributions to the local matrix and right-hand side
              for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                  if (assemble_system) {
                      for (unsigned int j = 0; j < dofs_per_cell; ++j) {

                          local_matrix(i, j) +=
                              (viscosity * scalar_product(grad_phi_u[j], grad_phi_u[i]) -
                              div_phi_u[i] * phi_p[j] -
                              phi_p[i] * div_phi_u[j] +
                              gamma * div_phi_u[j] * div_phi_u[i] +
                              phi_u[i] * phi_u[j] / timestep_NS) *
                              fe_values.JxW(q);
                              
                          local_mass_matrix(i, j) +=
                              (phi_u[i] * phi_u[j] + phi_p[i] * phi_p[j]) *
                              fe_values.JxW(q);
                      }
                  }

                  local_rhs(i) -=
                      (viscosity * scalar_product(current_velocity_gradients[q], grad_phi_u[i]) -
                      current_velocity_divergences[q] * phi_p[i] -
                      current_pressure_values[q] * div_phi_u[i] +
                      gamma * current_velocity_divergences[q] * div_phi_u[i] +
                      current_velocity_gradients[q] * current_velocity_values[q] * phi_u[i] -
                      scalar_product(phi_u[i], f)) * fe_values.JxW(q);
                }
          }

          cell->get_dof_indices(local_dof_indices);

          // Select the appropriate constraints
          const AffineConstraints<double> &constraints_used = use_nonzero_constraints ? nonzero_NS_constraints : zero_NS_constraints;

          // Distribute local contributions to global matrices and RHS
          constraints_used.distribute_local_to_global(local_matrix, local_rhs, local_dof_indices, NS_system_matrix, NS_system_rhs);
          
          if (assemble_system) {
            constraints_used.distribute_local_to_global(local_matrix,             //In practice this function implements a scatter operation
                                                          local_rhs,
                                                          local_dof_indices,        //Contains the corresponding global indexes
                                                          NS_system_matrix,
                                                          NS_system_rhs);

            constraints_used.distribute_local_to_global(local_mass_matrix, 
                                                        local_dof_indices, 
                                                        NS_mass_matrix);
          }

          constraints_used.distribute_local_to_global(local_rhs, local_dof_indices, NS_system_rhs);
      }

  }

  if (assemble_system) {
      NS_system_matrix.compress(VectorOperation::add);
      NS_mass_matrix.compress(VectorOperation::add);
  }

  NS_system_rhs.compress(VectorOperation::add);

}


//-----------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
std::pair<unsigned int, double>
CompleteProblem<dim>::solver_NS(bool use_nonzero_constraints, bool assemble_system, double time_step)
{
  // Reset the preconditioner if assembling the system
  if (assemble_system)
  {
    // Create a new BlockSchurPreconditioner for iterative solver
    preconditioner.reset(new BlockSchurPreconditioner(timer,
                                                      gamma,
                                                      viscosity,
                                                      timestep_NS,
                                                      owned_partitioning,
                                                      NS_system_matrix,
                                                      NS_mass_matrix,
                                                      pressure_mass_matrix));
  }
  
  // Initialize a coefficient for the stopping criteria of the iterative solver
  double coeff = 0.0 ;   
  if (time_step < 2) {
      coeff = 1e-2;
  } else {
      coeff = 1e-1;
  }

  // Initialize the solver control to manage the convergence of the solver
  SolverControl solver_control(10*NS_system_matrix.m(), coeff * NS_system_rhs.l2_norm(), true);

  // Instantiate the BiCGStab solver with the control settings
  SolverBicgstab<PETScWrappers::MPI::BlockVector> bicg(solver_control);

  pcout << "   NS_system_matrix frob norm is " << NS_system_matrix.frobenius_norm() << std::endl;
  pcout << "   NS_system_rhs l2 norm is " << NS_system_rhs.l2_norm() << std::endl;
  // The solution vector must be non-ghosted

  // Solve the linear system using the BiCGStab method
  // The solution vector (NS_solution_update) must be non-ghosted
  bicg.solve(NS_system_matrix, NS_solution_update, NS_system_rhs, *preconditioner);

  pcout<<"   Navier-Stokes solution computed, L2 norm: "<< NS_solution.l2_norm() << std::endl; 

  // Determine which constraints to apply based on the input flag
  const AffineConstraints<double> &constraints_used =
  use_nonzero_constraints ? nonzero_NS_constraints : zero_NS_constraints;
  constraints_used.distribute(NS_solution_update);

  // Return the last step taken by the solver and the last value of the residual
  return {solver_control.last_step(), solver_control.last_value()};
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::solve_navier_stokes()
{
  // Evaluate the electric field before solving the Navier-Stokes equations
	evaluate_electric_field();  

  NS_solution_update = 0;

  // Determine whether to apply non-zero constraints based on the time step
  bool apply_nonzero_constraints = (time_NS == 1);  

  bool assemble_system = true;

  pcout << "   ASSEMBLE NAVIER STOKES SYSTEM ..." ;
  assemble_NS(apply_nonzero_constraints, assemble_system);
  pcout << "   Done !"<<std::endl;


  pcout << "   SOLVE NAVIER STOKES SYSTEM ..."  << std::endl;
  auto state = solver_NS(apply_nonzero_constraints, assemble_system, time_NS);


  // Note we have to use a non-ghosted vector to do the addition.
  PETScWrappers::MPI::BlockVector tmp;
  tmp.reinit(owned_partitioning, mpi_communicator);   
  // Copy the current solution (which is a ghost vector and read-only) to the temporary vector
  tmp = NS_solution;
  tmp += NS_solution_update;
  // Update the current solution with the new values
  NS_solution = tmp;

  pcout << std::scientific << std::left << "   GMRES_ITR = " << std::setw(3)
      << state.first << "   GMRES_RES = " << state.second << std::endl;

  pcout << "   L2 norm of the present solution: " << NS_solution.l2_norm() << std::endl;
  pcout << "   Linf norm of the present solution: " << NS_solution.linfty_norm() << std::endl;

	pcout << "   Recovering velocity and pressure values for output... " << std::endl;

  // Temporary vectors for velocity components and pressure
  PETScWrappers::MPI::Vector temp_X;  
  PETScWrappers::MPI::Vector temp_Y;
  PETScWrappers::MPI::Vector temp_pressure;

  temp_X.reinit(locally_owned_dofs,  mpi_communicator);
  temp_Y.reinit(locally_owned_dofs,  mpi_communicator);
  temp_pressure.reinit(owned_partitioning[1], mpi_communicator); // non-ghosted pressure in order to look elements

  // Get the number of degrees of freedom (DoFs) per cell
	const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
	std::vector<types::global_dof_index> ion_local_dof_indices(dofs_per_cell);

  // Get DoFs for Navier-Stokes cells
	const unsigned int dofs_per_NS_cell = NS_fe.n_dofs_per_cell();
	std::vector<types::global_dof_index> NS_local_dof_indices(dofs_per_NS_cell);

	auto cell = dof_handler.begin_active();    
	auto NS_cell = NS_dof_handler.begin_active();

  // Define end iterators for both dof handlers
	const auto endc = dof_handler.end();
	const auto NS_endc = NS_dof_handler.end();

	double vel_max = 0.;  


  // Initialize min and max DoF values for X, Y, and pressure
  unsigned int min_dof_x = std::numeric_limits<unsigned int>::max();
  unsigned int max_dof_x = 0;
  unsigned int min_dof_y = std::numeric_limits<unsigned int>::max();
  unsigned int max_dof_y = 0;
  unsigned int min_dof_p = std::numeric_limits<unsigned int>::max();
  unsigned int max_dof_p = 0;

  // Index vector for accessing specific velocity degrees of freedom
  std::vector<int> index_x = {0, 3, 6, 9};

  // Set up the block component mapping for DoF renumbering
  std::vector<unsigned int> block_component(dim + 1, 0);    
  block_component[dim] = 1;
  DoFRenumbering::component_wise(NS_dof_handler, block_component);
  // Count the DoFs per block for the Navier-Stokes handler
  dofs_per_block = DoFTools::count_dofs_per_fe_block(NS_dof_handler, block_component);
  unsigned int dof_u = dofs_per_block[0];


  for (auto NS_cell = NS_dof_handler.begin_active(), ion_cell = dof_handler.begin_active();
      NS_cell != NS_dof_handler.end() && ion_cell != dof_handler.end();
      ++NS_cell, ++ion_cell)
   {
      // Ensure indices match between NS and ion cells
      Assert(NS_cell->index() == ion_cell->index(), ExcMessage("Mismatch between NS and ion cells!"));

      // Only proceed if both cells are locally owned
      if (NS_cell->is_locally_owned() && ion_cell->is_locally_owned()) {
          
          // Obtain local DoF indices for both NS and ion cells
          NS_cell->get_dof_indices(NS_local_dof_indices);
          ion_cell->get_dof_indices(ion_local_dof_indices);

          // Verify indices and block sizes
          for (unsigned int i = 0; i < 4; ++i) {
              const unsigned int dof_x = NS_local_dof_indices[index_x[i]];        // X velocity DoF index
              const unsigned int dof_y = NS_local_dof_indices[index_x[i] + 1];    // Y velocity DoF index
              const unsigned int dof_p = NS_local_dof_indices[3 * i + 2];         // Pressure DoF index

              // Check that DoF indices are within valid range before accessing
              if (dof_x < NS_solution.block(0).size() && dof_y < NS_solution.block(0).size() && dof_p < NS_solution.block(1).size() + dof_u) {
                  // Assign values to the temporary vectors if indices are valid
                  temp_X(ion_local_dof_indices[i]) = NS_solution.block(0)[dof_x];  
                  temp_Y(ion_local_dof_indices[i]) = NS_solution.block(0)[dof_y];  
                  temp_pressure(ion_local_dof_indices[i]) = NS_solution.block(1)[dof_p-dof_u];  
              } else {
                  // Log an error if any DoF index is out of bounds
                  std::cerr << "Errore: Indice DoF fuori dai limiti. DoF X: " << dof_x << ", DoF Y: " << dof_y << ", DoF P: " << dof_p << std::endl;
              }
          }

      }
  }

  // Compress the temporary vectors to finalize the assignments
  temp_X.compress(VectorOperation::insert);
  temp_Y.compress(VectorOperation::insert);
  temp_pressure.compress(VectorOperation::insert);

  // Assign the final values to class members for velocity and pressure
  Vel_X = temp_X;
  Vel_Y = temp_Y;
  pressure = temp_pressure;

}


//############ - RUN AND OUTPUT RESULTS - #######################################################################################################

template <int dim>
void CompleteProblem<dim>::run()
{

  // Fix the constant
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
  

  pcout << "   SETUP NAVIER STOKES PROBLEM ... "<< std::endl;
	setup_NS();

  pcout << "   STORE INITIAL CONDITIONS ... ";
  output_results(0);
  pcout << "   Done !"<<std::endl<<std::endl;


  // SET ERRORS AND TOLERANCES
  int step_number = 0;   // time evolution step

  const unsigned int max_it = 1e+3;    // Max iterations for Gummel
  const unsigned int max_steps = 1000;  // Max steps for the time loop

  const double tol = 1.e-9;         // Tolerance for the Newton 

  eta = old_ion_density;    // eta = N_0 function

  // START THE ALGORITHM
  pcout << "   START COMPLETE PROBLEM ... "<< std::endl;

  while (step_number < max_steps){
      
    ++step_number;

    if (step_number == 200 && timestep < 1.e-1){
    timestep*= 5.;
    timestep_NS*=5.;
    }

    int it = 0;                       // internal gummel iterator
    const double gummel_tol = 6.e-4;  // riferito al gummel algorithm
    double err = gummel_tol + 1.;     // errore gummel

    PETScWrappers::MPI::Vector previous_density;
    previous_density.reinit(locally_owned_dofs, mpi_communicator);

    pcout << "   GUMMEL ALGORITHM N. "<<step_number;

    while (err > gummel_tol && it < max_it) { // in questo ciclo NON si aggiorna old_ion_density

      solve_nonlinear_poisson(max_it,tol); // UPDATE potential AND eta (loop over k)
    
      previous_density = ion_density; // save the previously computed ion_density

      perform_drift_diffusion_fixed_point_iteration_step(); // UPDATE ion_density    

      previous_density -= ion_density;

      err = previous_density.linfty_norm()/ion_density.linfty_norm(); // compute error

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
    pcout << "   Density change from previous time-step is: " << time_err << std::endl;

    old_ion_density = ion_density;


    if (step_number % 10 == 1){ // NS solution update every 10 timesteps
      time_NS += 1;
      pcout << "   START NAVIER STOKES PROBLEM ... "<< std::endl; 
      solve_navier_stokes();

    }

    output_results(step_number);

  }
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::output_results(const unsigned int cycle)
{

// Base directory for output
std::string base_directory = "../output";

// Directory to store the results of this simulation
std::string output_directory = base_directory + "/NS_DD_Simulation/";

// Ensure the output directory is created (if it doesn't exist)
if (!std::filesystem::exists(output_directory))
{
    std::filesystem::create_directory(output_directory);
}

DataOut<dim> data_out;

data_out.attach_dof_handler(dof_handler);

data_out.add_data_vector(potential, "potential");
data_out.add_data_vector(pressure, "pressure");
data_out.add_data_vector(Vel_X, "Vel_X");
data_out.add_data_vector(Vel_Y, "Vel_Y");
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
// Tensor<1,2> get_emitter_normal(const Point<2> a) {

// 	Tensor<1,2> normal;

// 	normal[0] = a[0];
// 	normal[1] = a[1];

// 	const double norm = std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]);

// 	return normal/norm;
// }


Tensor<1,2> get_emitter_normal(const Point<2> &a, const Point<2> &emitter_center) {

    Tensor<1,2> normal;

    // Calcola la differenza tra il punto e il centro dell'emettitore
    normal[0] = a[0] - emitter_center[0];
    normal[1] = a[1] - emitter_center[1];

    // Calcola la norma del vettore
    const double norm = std::sqrt(normal[0]*normal[0] + normal[1]*normal[1]);

    // Evita la divisione per zero
    if (norm == 0) {
        throw std::runtime_error("Punto coincide con il centro dell'emettitore. La normale non è definita.");
    }

    // Restituisce il vettore normalizzato
    return normal / norm;
}