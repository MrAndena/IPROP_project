// Constructor of the class

template <int dim>
CompleteProblem<dim>::CompleteProblem(parallel::distributed::Triangulation<dim> &tria,  // triangulation for the mesh
                                      const data_struct &d,  // data struct from the user
                                      unsigned short int i)  // tag of the simulation
  : m_data(d)
  , simulation_tag(i)
  , mpi_communicator(MPI_COMM_WORLD)
  , pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  , triangulation(tria)
  , fe(1) //fe for poisson / DD
  , dof_handler(tria) //dof h for poisson / DD
  , mapping() // mapping for poisson / DD
  ,	viscosity(d.navier_stokes.physical_parameters.viscosity)
  , gamma(d.navier_stokes.physical_parameters.gamma) 
  , degree(1) 
  , NS_fe(FE_Q<dim>(degree+1), dim, FE_Q<dim>(degree), 1)
  , NS_dof_handler(tria)
  , volume_quad_formula(degree + 2) //non cerano nell'originale
  , face_quad_formula(degree + 2)   //non cerano nell'originale
  , NS_mapping() // serve ? nel nostro NS non cè
  , timestep(1.e-5) // vogliamo passare uno dalla struct ? è diverso da NS 
  , time_NS(d.navier_stokes) // come in INSIMEX
  , timer(mpi_communicator, pcout, TimerOutput::never, TimerOutput::wall_times)
{}

//---------------------------------------------------------------------------------------------------------------------------------------------------------

  template <int dim>
  void CompleteProblem<dim>::setup_poisson()
  {
          // fix the constants
    const double Re = m_data.geometrical_parameters.emitter_radius[simulation_tag]; // per adesso un singolo valore
    const double E_ON = m_data.drift_diffusion.physical_parameters.E_ON;
    const double Ve = m_data.drift_diffusion.physical_parameters.Ve;

    // const double N_0 = m_data.drift_diffusion.physical_parameters.stratosphere ? 2.2e-3 : 0.5e-3; // [m^-3] ambient ion density  NON SERVE IN QUESTA FUNZIONE
    const double p_amb = m_data.drift_diffusion.physical_parameters.stratosphere ? 5474 : 101325; 
    const double T = m_data.drift_diffusion.physical_parameters.stratosphere ? 217. : 303.; // [K] fluid temperature
    // const double rho = m_data.drift_diffusion.physical_parameters.stratosphere ? 0.089 : 1.225; // kg m^-3     !!GUARDA CHE DOPO LA DEFINISCO A MANO

    const double delta = p_amb/101325*298/T; //                                       
        
    const double eps = 1.; // wire surface roughness correction coefficient
    const double Ep = E_ON*delta*eps*(1+0.308/std::sqrt(Re*1.e+2*delta));
    // const double Ri = Ep/E_ON*Re; // [m] ionization radius
    const double Vi = Ve - Ep*Re*log(Ep/E_ON); // [V] voltage on ionization region boundary

    dof_handler.distribute_dofs(fe);

    // INDEX SETS INITIALIZATION
    locally_owned_dofs = dof_handler.locally_owned_dofs();                           //local dofs
    locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);    //local dofs + ghost dofs
    
    // PETSC VECTORS DECLARATIONS 
    potential.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);   //ghosted
    poisson_newton_update.reinit(locally_owned_dofs, mpi_communicator);              //non-ghosted
    poisson_rhs.reinit(locally_owned_dofs, mpi_communicator);                        //non-ghosted

    Field_X.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);     //ghosted, servono in output, campo elettrico
    Field_Y.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);     //ghosted

    // ZERO_CONSTRAINTS FOR NEWTON POISSON PROBLEM (coerente con il nostro modo di procedere)
    zero_constraints_poisson.clear();
    zero_constraints_poisson.reinit(locally_relevant_dofs);
    zero_constraints_poisson.close(); 
    
    //POTREBBERO SERVIRE ADDED !! nuovo
    constraints_poisson.clear();
    constraints_poisson.reinit(locally_relevant_dofs);
    VectorTools::interpolate_boundary_values(dof_handler,3, Functions::ConstantFunction<dim>(Vi), constraints_poisson); //emitter
    VectorTools::interpolate_boundary_values(dof_handler, 4,Functions::ZeroFunction<dim>(), constraints_poisson); // inlet
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

}
//----------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::initialize_potential() //non viene più usato!!
{ 
  // fix the constants
  const double Re = m_data.geometrical_parameters.emitter_radius[simulation_tag]; // per adesso un singolo valore
  const double E_ON = m_data.drift_diffusion.physical_parameters.E_ON;
  const double Ve = m_data.drift_diffusion.physical_parameters.Ve;

  // const double N_0 = m_data.drift_diffusion.physical_parameters.stratosphere ? 2.2e-3 : 0.5e-3; // [m^-3] ambient ion density  NON SERVE IN QUESTA FUNZIONE
  const double p_amb = m_data.drift_diffusion.physical_parameters.stratosphere ? 5474 : 101325; 
  const double T = m_data.drift_diffusion.physical_parameters.stratosphere ? 217. : 303.; // [K] fluid temperature
  // const double rho = m_data.drift_diffusion.physical_parameters.stratosphere ? 0.089 : 1.225; // kg m^-3     !!GUARDA CHE DOPO LA DEFINISCO A MANO

  const double delta = p_amb/101325*298/T; //                                       
     
  const double eps = 1.; // wire surface roughness correction coefficient
  const double Ep = E_ON*delta*eps*(1+0.308/std::sqrt(Re*1.e+2*delta));
  // const double Ri = Ep/E_ON*Re; // [m] ionization radius
  const double Vi = Ve - Ep*Re*log(Ep/E_ON); // [V] voltage on ionization region boundary

  PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator); 

  temp = potential;
  
  // IC potential: everywhere Vi expect on the collector: zero, PROVARE AL CONTRARIO MA NON CONVERGE
  //VectorTools::interpolate(mapping, dof_handler, Functions::ConstantFunction<dim>(Vi) , temp);

  //IC potential: everywhere null except on the emitter: Vi
  VectorTools::interpolate(mapping, dof_handler, Functions::ZeroFunction<dim>(), temp);
  
  std::map<types::global_dof_index, double> boundary_values;

  VectorTools::interpolate_boundary_values(dof_handler, 3, Functions::ConstantFunction<dim>(Vi), boundary_values);
  VectorTools::interpolate_boundary_values(dof_handler, 4, Functions::ZeroFunction<dim>(), boundary_values);

  for (auto &boundary_value : boundary_values){
    temp(boundary_value.first) = boundary_value.second;
  }

  temp.compress(VectorOperation::insert);
  potential = temp;

}
//--------------------------------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  void CompleteProblem<dim>::assemble_poisson_laplace_matrix() // ADESSO CREA ANCHE IL RHS POTENTIAL
  {
    const QTrapezoid<dim> quadrature_formula;

    laplace_matrix_poisson = 0;
    poisson_rhs = 0; // ADDED

    FEValues<dim> fe_values(fe,
			    quadrature_formula,
			    update_values | update_gradients |
			    update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell); // ADDED

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  
    unsigned int q_point = 0, idof = 0, jdof = 0;
    for (const auto &cell : dof_handler.active_cell_iterators()){
      if (cell->is_locally_owned()){
	      cell_matrix = 0.;
        cell_rhs = 0.;
        fe_values.reinit(cell);
	
        for (q_point = 0; q_point < n_q_points; ++q_point) {
	  for (idof = 0; idof < dofs_per_cell; ++idof) {
	    for (jdof = 0; jdof < dofs_per_cell; ++jdof)
	      cell_matrix(idof, jdof) += fe_values.shape_grad(idof, q_point) *
		fe_values.shape_grad(jdof, q_point) * fe_values.JxW(q_point);
	  }
	}

        cell->get_dof_indices(local_dof_indices);
        constraints_poisson.distribute_local_to_global(cell_matrix,         // PRIMA ERA zero 
                                                            cell_rhs,       //ADDED  
                                                            local_dof_indices,
                                                            laplace_matrix_poisson,
                                                            poisson_rhs);
      }
    }

    laplace_matrix_poisson.compress(VectorOperation::add);
    poisson_rhs.compress(VectorOperation::add); // ADDED

    // pcout << " The L_INF norm of the laplace matrix is "<<laplace_matrix.linfty_norm() <<std::endl;
    // pcout << " The L_FROB norm of the laplace matrix is "<<laplace_matrix.frobenius_norm() <<std::endl<<std::endl;
    // pcout << "   End of Assembling Laplce matrix "<< std::endl<<std::endl;
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

    // pcout << " The L_INF norm of the mass matrix is "<<mass_matrix.linfty_norm() <<std::endl;
    // pcout << " The L_FROB norm of the mass matrix is "<<mass_matrix.frobenius_norm() <<std::endl<<std::endl;
    // pcout << "   End of Assembling Mass matrix "<< std::endl<<std::endl;

  }
  
//----------------------------------------------------------------------------------------------------------------------------  
  
template <int dim>
void CompleteProblem<dim>::setup_drift_diffusion(const bool reinitialize_densities)
{ 
  //fix the constants
  const double E_ON = m_data.drift_diffusion.physical_parameters.E_ON;
  const double E_ref = m_data.drift_diffusion.physical_parameters.E_ref;
  const double N_ref = m_data.drift_diffusion.physical_parameters.N_ref;
  const double N_0 = m_data.drift_diffusion.physical_parameters.stratosphere ? 2.2e-3 : 0.5e-3; // [m^-3] ambient ion density 

	if (reinitialize_densities) {  // questa condizione  non penso a noi serva
		eta.reinit(locally_owned_dofs, mpi_communicator); //non-ghosted, non serve come output
		ion_density.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator); // ghosted, lo vogliamo in output
		old_ion_density.reinit(locally_owned_dofs, mpi_communicator); // non-ghosted non serve come output (magari però utile per confronti)
  }

	// Corona inception condition
  // Functions::FEFieldFunction<dim> solution_as_function_object(dof_handler, potential, mapping);
  Functions::FEFieldFunction<dim, dealii::PETScWrappers::MPI::Vector> solution_as_function_object(dof_handler, potential, mapping);


	auto boundary_evaluator = [&] (const Point<dim> &p)
		{
			Tensor<1,dim> grad_U = solution_as_function_object.gradient(p);

			const double EXP = std::exp((grad_U.norm()-E_ON)/E_ref);

			return N_ref * EXP;
		};

    
  //vedere commento tesi pag. 38 per capire queste bcs 

  // DENSITY IONS CONSTRAINTS
	ion_constraints.clear();
	ion_constraints.reinit(locally_relevant_dofs);
	VectorTools::interpolate_boundary_values(dof_handler,3, ScalarFunctionFromFunctionObject<2>(boundary_evaluator), ion_constraints); //emitter
  VectorTools::interpolate_boundary_values(dof_handler, 1, Functions::ConstantFunction<dim>(N_0), ion_constraints); // inlet
	ion_constraints.close();
  

  DynamicSparsityPattern ion_dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, ion_dsp, ion_constraints, false);
  SparsityTools::distribute_sparsity_pattern(ion_dsp,
                                             dof_handler.locally_owned_dofs(),
                                             mpi_communicator,
                                             locally_relevant_dofs);
  
  // DENSITY VECTORS AND MATRICES
	ion_mass_matrix.clear();
  ion_mass_matrix.reinit(locally_owned_dofs, locally_owned_dofs, ion_dsp,  mpi_communicator); 

  // questa matrice viene costruita con un metodo nuovo e parallelo
	//MatrixCreator::create_mass_matrix(mapping, dof_handler, QTrapezoid<dim>(),ion_mass_matrix, (const Function<dim> *const)nullptr, ion_constraints);

	ion_rhs.reinit(locally_owned_dofs, mpi_communicator); 

	ion_system_matrix.clear();
  ion_system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, ion_dsp,  mpi_communicator); 
	
  drift_diffusion_matrix.clear();
  drift_diffusion_matrix.reinit(locally_owned_dofs, locally_owned_dofs, ion_dsp,  mpi_communicator); 
}
//------------------------------------------------------------------------------------------------------------------------------

template <int dim>
void CompleteProblem<dim>::assemble_drift_diffusion_mass_matrix()
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
      for (jdof = 0; jdof < dofs_per_cell; ++jdof)
  cell_matrix(idof, jdof) += fe_values.shape_value(idof, q_point) *
    fe_values.shape_value(jdof, q_point) * fe_values.JxW(q_point);
    }
  }

  cell->get_dof_indices(local_dof_indices);
  ion_constraints.distribute_local_to_global(cell_matrix,   
                                             local_dof_indices,
                                             ion_mass_matrix);
}
  }

  ion_mass_matrix.compress(VectorOperation::add);

  // pcout << " The L_INF norm of the mass matrix is "<<mass_matrix.linfty_norm() <<std::endl;
  // pcout << " The L_FROB norm of the mass matrix is "<<mass_matrix.frobenius_norm() <<std::endl<<std::endl;
  // pcout << "   End of Assembling Mass matrix "<< std::endl<<std::endl;

}
//---------------------------------------------------------------------------------------------------------------------------------------------------

template <int dim>
void CompleteProblem<dim>::assemble_nonlinear_poisson()
  { 
    // Fix the constants
    const double q0 = m_data.drift_diffusion.physical_parameters.q0;
    const double eps_r = m_data.drift_diffusion.physical_parameters.eps_r;
    const double eps_0 = m_data.drift_diffusion.physical_parameters.eps_0;
    const double kB = m_data.drift_diffusion.physical_parameters.kB;
    const double T = m_data.drift_diffusion.physical_parameters.stratosphere ? 217. : 303.; 
    const double V_TH = kB*T/q0;

    //BUILDING POISSON SYSTEM MATRIX (for newton)
    system_matrix_poisson = 0;
    ion_mass_matrix = 0;

    double new_value = 0;
  
    // Generate the term:  (eta)*MASS_MAT   lumped version stored in ion_mass_matrix
    for (auto iter = locally_owned_dofs.begin(); iter != locally_owned_dofs.end(); ++iter){ 

      new_value = mass_matrix_poisson(*iter, *iter) * eta(*iter);
      ion_mass_matrix.set(*iter,*iter,new_value);

    }

    ion_mass_matrix.compress(VectorOperation::insert);
    //bisognerà mettere le bcs per la ion mass????????????? non pesno, c'è distribute local to global quando viene creata

    system_matrix_poisson.add(eps_r * eps_0, laplace_matrix_poisson); // SYS_MAT = SYS_MAT +  eps*A
    system_matrix_poisson.add(q0 / V_TH, ion_mass_matrix);     // SYS_MAT = SYS_MAT + q0/V_TH * (eta)*MASS_MAT

  
    // BUILDING SYSTEM RHS
    poisson_rhs = 0;

    PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator); //temporary non-ghosted vector number 1
    mass_matrix_poisson.vmult(temp,eta);
    poisson_rhs.add(q0,temp);
    laplace_matrix_poisson.vmult(temp,potential); //se con potential non funziona, creare il secondo temp
    poisson_rhs.add(- eps_r * eps_0, temp);

    //coerentemente con il nostro approccio non ci sono apply bcs

    //pcout << "   The L_INF norm of the poisson system RHS is: "<<poisson_system_rhs.linfty_norm() <<std::endl;
    //pcout << "   The L2 norm of the poisson system RHS is: " << poisson_system_rhs.l2_norm() << std::endl;
    //pcout << "   End of solve assemble non-linear poisson"<< std::endl<<std::endl;

  }
  //----------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  double CompleteProblem<dim>::solve_poisson()
  { 

    // Fixing costants
    const double q0 = m_data.drift_diffusion.physical_parameters.q0;
    const double kB = m_data.drift_diffusion.physical_parameters.kB;
    const double T = m_data.drift_diffusion.physical_parameters.stratosphere ? 217. : 303.; 
    const double V_TH = kB*T/q0;

    //Apply zero boundary conditions to the whole newton poisson system
    //We apply the BCs on tags 3 (emitter) and 4 (collector)

    std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

    VectorTools::interpolate_boundary_values(mapping, dof_handler,3, Functions::ZeroFunction<dim>(), emitter_boundary_values);
    MatrixTools::apply_boundary_values(emitter_boundary_values, system_matrix_poisson, poisson_newton_update, poisson_rhs);

    VectorTools::interpolate_boundary_values(mapping, dof_handler,4, Functions::ZeroFunction<dim>(), collector_boundary_values);
    MatrixTools::apply_boundary_values(collector_boundary_values, system_matrix_poisson, poisson_newton_update, poisson_rhs);
    

    //Solve poisson system problem
    SolverControl sc_p(dof_handler.n_dofs(), 1e-10);     
    PETScWrappers::SparseDirectMUMPS solverMUMPS(sc_p); 
    solverMUMPS.solve(system_matrix_poisson, poisson_newton_update, poisson_rhs);

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
    temp.add(0.9, poisson_newton_update);

    potential = temp;

    return residual;

    // pcout << "L2 norm of the current solution: " << current_solution.l2_norm() << std::endl;
    // pcout << "L_INF norm of the current solution: " << current_solution.linfty_norm() << std::endl;
    // pcout << "   End of solve poisson problem"<< std::endl<<std::endl;
  
}
//--------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::solve_homogeneous_poisson()
{

	system_matrix_poisson.copy_from(laplace_matrix_poisson);

  PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator); 

  //Solve homo poisson system problem
  SolverControl sc_p(dof_handler.n_dofs(), 1e-10);     
  PETScWrappers::SparseDirectMUMPS solverMUMPS(sc_p); 
  solverMUMPS.solve(system_matrix_poisson, temp, poisson_rhs);

  
	constraints_poisson.distribute(temp);
  temp.compress(VectorOperation::insert);

  potential = temp;


}

//-----------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::solve_nonlinear_poisson(const unsigned int max_iter_newton, 
                                                   const double toll_newton){

  unsigned int counter = 0; // it keeps track of newton iteration

  pcout << "   - START POISSON NEWTON - "<< std::endl;

  double increment_norm = std::numeric_limits<double>::max(); // the increment norm is + inf

  pcout << "   Initial Newton Increment Norm dphi: " << increment_norm << std::endl; 
  

  while(counter < max_iter_newton && increment_norm > toll_newton){

    pcout << "   NEWTON ITERATION NUMBER: "<< counter +1<<std::endl;
    pcout << "   Assemble System Poisson Matrix"<< std::endl; 

    //NB: Mass and Laplace matrices are already build
    
    assemble_nonlinear_poisson();
    increment_norm = solve_poisson();  //residual computation, clamping on newton update, BCs and update of the charges are inside this method  
    
    pcout << "   Update Increment: "<<increment_norm<<std::endl<<std::endl;

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
  const double kB = m_data.drift_diffusion.physical_parameters.kB;
  const double Mm = m_data.drift_diffusion.physical_parameters.Mm;
  const double Avo = m_data.drift_diffusion.physical_parameters.Avo;
  const double q0 = m_data.drift_diffusion.physical_parameters.q0;
  const double mu0 = m_data.drift_diffusion.physical_parameters.mu0;
  const double T = m_data.drift_diffusion.physical_parameters.stratosphere ? 217. : 303.; // [K] fluid temperature
  const double p_amb = m_data.drift_diffusion.physical_parameters.stratosphere ? 5474 : 101325; 
  const double delta = p_amb/101325*298/T; // 

  const double mu = mu0 * delta; // scaled mobility from Moseley                   
  const double V_TH = kB * T / q0; // [V] ion temperature                           
  const double D = mu * V_TH; //                                                    

  // start to build the matrix
  drift_diffusion_matrix = 0;

  const unsigned int vertices_per_cell = 4;
  FullMatrix<double> Robin(vertices_per_cell,vertices_per_cell);

  Vector<double> cell_rhs(vertices_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(vertices_per_cell);
  FullMatrix<double> cell_matrix(vertices_per_cell, vertices_per_cell);

  const unsigned int t_size = 3;

  Vector<double> A_cell_rhs(t_size), B_cell_rhs(t_size);
  FullMatrix<double> A(t_size,t_size), B(t_size,t_size);

  std::vector<types::global_dof_index> A_local_dof_indices(t_size);
  std::vector<types::global_dof_index> B_local_dof_indices(t_size);

  pcout << "   FATTO ASS 1 "<< std::endl;
  evaluate_electric_field();
  pcout << "   FATTO ASS 2 "<< std::endl;

  const double Vh = std::sqrt(8.*numbers::PI*kB * T / Mm * Avo / 2. / numbers::PI); // Hopf velocity
  //const double Vh = std::sqrt(kB * T / Mm * Avo / 2. / numbers::PI); // Thermal velocity


  QTrapezoid<dim-1>	face_quadrature;
  const unsigned int n_q_points = face_quadrature.size();
  FEFaceValues<dim> face_values(fe, face_quadrature, update_values | update_quadrature_points /*| update_normal_vectors*/ | update_JxW_values);

  pcout << "   FATTO ASS 3 "<< std::endl;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
        if (cell->is_locally_owned()){

            A = 0;
            B = 0;
            cell_matrix = 0;
            cell_rhs = 0;
            Robin = 0;
            A_cell_rhs = 0;
            B_cell_rhs = 0;
            cell->get_dof_indices(local_dof_indices);

            // Robin conditions at outlet and (optional) at collector
            if (cell->at_boundary()) {
                for (const auto &face : cell->face_iterators()) {
                    if (face->at_boundary() && face->boundary_id() == 2) { // Outlet

                        face_values.reinit(cell,face);

                        for (unsigned int i = 0; i < vertices_per_cell; ++i) {

                            const double vel_f = Vel_X(local_dof_indices[i]);
                            
                            for (unsigned int q = 0; q < n_q_points; ++q) {
                                for (unsigned int j = 0; j < vertices_per_cell; ++j) {
                                    Robin(i,j) += face_values.JxW(q) * face_values.shape_value(i,q) * face_values.shape_value(j,q) * vel_f;
                                }
                            }
                        }

                    } else if (face->at_boundary() && (face->boundary_id() == 4)) { // Collector
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

            if (l_alpha >= l_beta) { // l_alpha is the longest diagonal: split by beta
                        const double l_23 = side_length(v2,v3);
                        const Tensor<1,dim> dir_23 = (v3 - v2)/l_beta;

                        const double alpha23 = (u_f_2 * dir_23)/D*l_23 + (u3 - u2);

                        // Triangle A:
                        A= compute_triangle_matrix(v2,v1,v3, alpha21, alpha13, -alpha23, this->m_data);

                        // Triangle B:
                        B = compute_triangle_matrix(v3,v4,v2, alpha34, alpha42, alpha23, this->m_data);

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
                        A = compute_triangle_matrix(v4,v2,v1, alpha42, alpha21, alpha14, this->m_data);

                        // Triangle B:
                        B = compute_triangle_matrix(v1,v3,v4, alpha13, alpha34, -alpha14, this->m_data);

                        A_local_dof_indices[0] = local_dof_indices[1];
                        A_local_dof_indices[1] = local_dof_indices[3];
                        A_local_dof_indices[2] = local_dof_indices[2];

                        B_local_dof_indices[0] = local_dof_indices[2];
                        B_local_dof_indices[1] = local_dof_indices[0];
                        B_local_dof_indices[2] = local_dof_indices[1];
                    }

                    // As the ion system matrix is M + delta t DD, the contributions are multiplied by the timestep
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
   pcout << "   FATTO ass 4 "<< std::endl;
   //i compress prima non c'erano
   ion_system_matrix.compress(VectorOperation::add);
   ion_rhs.compress(VectorOperation::add);
   
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::solve_drift_diffusion()
{ 
  // aggiornare e mettere come il nostro
  PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator);
  SolverControl sc_ion(dof_handler.n_dofs(), 1e-10);
  PETScWrappers::SparseDirectMUMPS solverMUMPS_ion(sc_ion); 
  solverMUMPS_ion.solve(ion_system_matrix, temp, ion_rhs);

  ion_constraints.distribute(temp); // qua constraints per le ions (paragona a nostro DD)
  temp.compress(VectorOperation::insert);

  ion_density = temp;

  /*
  pcout << "   L_INF norm of the hole_density: " << hole_density.linfty_norm() << std::endl;
  pcout << "   L_INF norm of the electron_density: " << electron_density.linfty_norm() << std::endl;
  */

}
//--------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::perform_drift_diffusion_fixed_point_iteration_step() // method used to update ion_density
{  
  
  PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator);

	ion_rhs = 0;
	ion_system_matrix = 0;

	ion_mass_matrix.vmult(temp, old_ion_density);
	ion_rhs += temp;

  
	// Integration in time with BE:
	ion_system_matrix.copy_from(ion_mass_matrix);


  assemble_drift_diffusion_matrix();  //INDAGARA QUA DENTRO

  pcout << "   FATTO D2 "<< std::endl;

    // Uncomment to add a non-zero forcing term to DD equations ...
	/*
    Vector<double> forcing_terms(old_ion_density.size());
    VectorTools::create_right_hand_side(dof_handler, QTrapezoid<dim>(), Functions::ZeroFunction<dim>(), tmp); // ... by changing the ZeroFunction to an appropriate one
	forcing_terms = tmp;
	forcing_terms *= timestep;
	ion_rhs += forcing_terms;
	*/

  solve_drift_diffusion(); // INDAGARE QUA DENTRO
  pcout << "   FATTO D3 "<< std::endl;

}
//------------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::evaluate_electric_field()
{

  Field_X.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);     
  Field_Y.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);     

  const unsigned int 		dofs_per_cell = fe.n_dofs_per_cell(); // GRADI DI LIBERTA PER OGNI CELLA

  //std::vector<double>		global_dof_hits(dof_handler.n_dofs()); //sicuro funzioni ? noi dobbiamo mettere quelli del processore

  std::vector<double> global_dof_hits(locally_owned_dofs.n_elements());


  PETScWrappers::MPI::Vector		el_field_X(locally_owned_dofs, mpi_communicator);
  PETScWrappers::MPI::Vector		el_field_Y(locally_owned_dofs, mpi_communicator);

	QTrapezoid<dim-1>			    iv_quadrature;
	FEInterfaceValues<dim> 		fe_iv(fe, iv_quadrature, update_gradients);

	const unsigned int 		      n_q_points = iv_quadrature.size();
  std::vector<Tensor<1,dim>> 	iv_gradients(n_q_points);

	std::vector<types::global_dof_index> 	local_dof_indices(dofs_per_cell);
    

	for (auto &cell : dof_handler.active_cell_iterators()){

    if (cell->is_locally_owned()){
      
      for (const auto face_index : GeometryInfo<dim>::face_indices())
          {

            fe_iv.reinit(cell, face_index);
            local_dof_indices = fe_iv.get_interface_dof_indices();

            fe_iv.get_average_of_function_gradients(potential, iv_gradients);

            for (const auto q : fe_iv.quadrature_point_indices()) {
              for (const auto i : fe_iv.dof_indices()) {

                    global_dof_hits[local_dof_indices[i]] += 1.;

                    for (unsigned int d = 0; d < dim; ++d) {

                      if (d == 0)
                        el_field_X(local_dof_indices[i]) += - iv_gradients[q][d]; //-grad_phi_i[d] / n_q_points;
                      else if (d == 1)
                        el_field_Y(local_dof_indices[i]) +=  - iv_gradients[q][d]; // * (J_inverse[d] * shape_gradient) * dx;
                      else
                        Assert(false, ExcNotImplemented());
                    }

              }
            }
          }
      
    }
   }

   el_field_X.compress(VectorOperation::add);
   el_field_Y.compress(VectorOperation::add);

   // Take the average of all dof values
   /*
   for (unsigned int k = 0; k < dof_handler.n_dofs(); k++) {
	   el_field_X(k) /= std::max(1., global_dof_hits[k]);
	   el_field_Y(k) /= std::max(1., global_dof_hits[k]);
   }*/

  // è giusta scritta così? il pezzo sotto?

  for (auto iter = locally_owned_dofs.begin(); iter != locally_owned_dofs.end(); ++iter){ 

    el_field_X[*iter] /= std::max(1., global_dof_hits[*iter]);
    el_field_Y[*iter] /= std::max(1., global_dof_hits[*iter]);

  }


   el_field_X.compress(VectorOperation::insert);
   el_field_Y.compress(VectorOperation::insert);

   Field_X = el_field_X;
   Field_Y = el_field_Y;

   
}
  
//##################### - Navier-Stokes - #####################################################################################

template <int dim>
void CompleteProblem<dim>::setup_NS()
  {
    
    //SET UP DOFS
    NS_dof_handler.distribute_dofs(NS_fe);

	  DoFRenumbering::Cuthill_McKee(NS_dof_handler);  //Renumber the degrees of freedom according to the Cuthill-McKee method.

    std::vector<unsigned int> block_component(dim + 1, 0);   //Musk for the reording 
    block_component[dim] = 1;
    DoFRenumbering::component_wise(NS_dof_handler, block_component);

    dofs_per_block = DoFTools::count_dofs_per_fe_block(NS_dof_handler, block_component);

    // Partitioning.
    unsigned int dof_u = dofs_per_block[0];
    unsigned int dof_p = dofs_per_block[1];

	  owned_partitioning.resize(2);
    owned_partitioning[0] = NS_dof_handler.locally_owned_dofs().get_view(0, dof_u);      //Extract the set of locally owned DoF indices for each component within the mask that are owned by the current processor.
    owned_partitioning[1] = NS_dof_handler.locally_owned_dofs().get_view(dof_u, dof_u + dof_p);

    DoFTools::extract_locally_relevant_dofs(NS_dof_handler,NS_locally_relevant_dofs );     //Extract the set of global DoF indices that are active on the current DoFHandler. This is the union of DoFHandler::locally_owned_dofs() and the DoF indices on all ghost cells.

    relevant_partitioning.resize(2);
    relevant_partitioning[0] = NS_locally_relevant_dofs.get_view(0, dof_u);
    relevant_partitioning[1] = NS_locally_relevant_dofs.get_view(dof_u, dof_u + dof_p);

   
    //MAKE CONSTRAINTS
    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar vertical_velocity(1);
    const FEValuesExtractors::Vector vertical_velocity_and_pressure(1);

    nonzero_NS_constraints.clear();
    nonzero_NS_constraints.reinit(NS_locally_relevant_dofs); 

    
    DoFTools::make_hanging_node_constraints(NS_dof_handler, nonzero_NS_constraints);  //Lo lasciamo nel caso si faccia refine mesh adattivo
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
                                            BoundaryValues<dim>(), // Functions::ZeroFunction<dim>(dim+1), 
                                            nonzero_NS_constraints, 
                                            NS_fe.component_mask(velocities));
                                            
  
    VectorTools::interpolate_boundary_values(NS_dof_handler,
                                            2,                    // Outlet
                                            Functions::ZeroFunction<dim>(dim+1),//BoundaryValues<dim>()
                                            nonzero_NS_constraints,
                                            NS_fe.component_mask(vertical_velocity_and_pressure));  //Vertical velocity and pressure at outlet equal to 0
                              
    
    nonzero_NS_constraints.close();
    
    
    zero_NS_constraints.clear();                // The IndexSets of owned velocity and pressure respectively.
	  zero_NS_constraints.reinit(NS_locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(NS_dof_handler, zero_NS_constraints);

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
    
    zero_NS_constraints.close();


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


    //??Farlo ereditare anche a Vel_X e Vel_Y ??
   
    owned_partitioning_u = NS_dof_handler.locally_owned_dofs().get_view(0, dof_u);      
    owned_partitioning_p = NS_dof_handler.locally_owned_dofs().get_view(dof_u, dof_u + dof_p);
    

    Vel_X.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);  //NB dof handler di DD !!! va così
    Vel_Y.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

    pressure.reinit(owned_partitioning_p, mpi_communicator); // lui non da problemi in output, inizializzartlo come VelX e Y??
  
  }



//-----------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::assemble_NS(bool use_nonzero_constraints,
                                       bool assemble_system)
{ 
  
  TimerOutput::Scope timer_section(timer, "Assemble system");        //Enter the given section in the timer

  if (assemble_system)           //If assemble_system is false i assembly only rhs
  {
      NS_system_matrix = 0;
      NS_mass_matrix = 0;
  }

  NS_system_rhs = 0;

  FEValues<dim> fe_values(NS_fe,                              //FEValues contains finite element evaluated in quadrature points of a cell.                            
                          volume_quad_formula,             //It implicitely uses a Q1 mapping
                          update_values | update_quadrature_points |
                          update_JxW_values | update_gradients);
  /*
  FEFaceValues<dim> fe_face_values(NS_fe,
                                  face_quad_formula,
                                  update_values | update_normal_vectors |
                                  update_quadrature_points |
                                  update_JxW_values); INUTILE NON VIENE USATO */

  const unsigned int dofs_per_cell = NS_fe.dofs_per_cell;
  const unsigned int n_q_points = volume_quad_formula.size();

  const FEValuesExtractors::Vector velocities(0);   //Extractor calls velocities that takes the vector in position 0
  const FEValuesExtractors::Scalar pressure(dim);   //Extractor calls pressure that takes the scalar in position dim

  FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> local_mass_matrix(dofs_per_cell, dofs_per_cell);

  Vector<double> local_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  std::vector<Tensor<1, dim>> current_velocity_values(n_q_points);     //vector of values of velocity in quadrature points              
  std::vector<Tensor<2, dim>> current_velocity_gradients(n_q_points);
  std::vector<double> current_velocity_divergences(n_q_points);
  std::vector<double> current_pressure_values(n_q_points);


  std::vector<double> div_phi_u(dofs_per_cell);
  std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
  std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
  std::vector<double> phi_p(dofs_per_cell);

  Tensor<1,dim> f;

  auto ion_cell = dof_handler.begin_active();  
  const auto ion_endc = dof_handler.end();

  std::vector<types::global_dof_index> ion_local_dof_indices(4);


  for (auto cell = NS_dof_handler.begin_active(); cell != NS_dof_handler.end(); ++cell)     //Iterator from the first active cell to the last one
  {   
      
      if (cell->is_locally_owned())
      {
          fe_values.reinit(cell);    //Reinitialize the gradients, Jacobi determinants, etc for the given cell of type "iterator into a Triangulation object", and the given finite element.

          if (assemble_system)
          {
              local_matrix = 0;
              local_mass_matrix = 0;
          } 

          local_rhs = 0;
          
          fe_values[velocities].get_function_values(NS_solution,
                                                    current_velocity_values);     //in current_velocity_values i stored values of NS_solution in quadrature points

          fe_values[velocities].get_function_gradients(NS_solution, 
                                                       current_velocity_gradients);

          fe_values[velocities].get_function_divergences(NS_solution, 
                                                         current_velocity_divergences);

          fe_values[pressure].get_function_values(NS_solution,
                                                  current_pressure_values);
          if (ion_cell->is_locally_owned())
          {
          if (ion_cell != ion_endc){
        	  ion_cell->get_dof_indices(ion_local_dof_indices);
          }  else{
        	  pcout << "Warning! Reached end of ion cells at NS cell " << cell->index() << std::endl;
          }
          }
          // Assemble the system matrix and mass matrix simultaneouly.
          // The mass matrix only uses the $(0, 0)$ and $(1, 1)$ blocks.
          
          for (unsigned int q = 0; q < n_q_points; ++q)
          { 
            
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
            {
              div_phi_u[k] = fe_values[velocities].divergence(k, q);           //Returns the value of the k-th shape function in q quadrature point
              grad_phi_u[k] = fe_values[velocities].gradient(k, q);
              phi_u[k] = fe_values[velocities].value(k, q);
              phi_p[k] = fe_values[pressure].value(k, q); 
            }
        
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            { 
              if (assemble_system)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                { 
                  
                  local_matrix(i, j) +=
                  (viscosity *
                   scalar_product(grad_phi_u[j], grad_phi_u[i]) -
                   div_phi_u[i] * phi_p[j] -
                   phi_p[i] * div_phi_u[j] +
                   gamma * div_phi_u[j] * div_phi_u[i] +
                   phi_u[i] * phi_u[j] / time_NS.get_delta_t() //+

                  // phi_u[i] * (current_velocity_gradients[q] * phi_u[j]) +
                  // phi_u[i] * (grad_phi_u[j] * current_velocity_values[q]) 
                  ) *
                  fe_values.JxW(q);
             
                  local_mass_matrix(i, j) +=
                  (phi_u[i] * phi_u[j] + phi_p[i] * phi_p[j]) *
                  fe_values.JxW(q);
                  
                }
              }

              const double q0 = m_data.drift_diffusion.physical_parameters.q0;

              const double rho = m_data.drift_diffusion.physical_parameters.stratosphere ? 0.089 : 1.225; // kg m^-3

              //const double rho = 1.225; // kg m^-3  //Attenzione: buono solo se non nella stratosphere
              
              /*
              PETScWrappers::MPI::Vector temp_X(locally_owned_dofs, mpi_communicator);
              PETScWrappers::MPI::Vector temp_Y(locally_owned_dofs, mpi_communicator);
              PETScWrappers::MPI::Vector temp_ion(locally_owned_dofs, mpi_communicator);

              temp_X = Field_X;
              temp_Y = Field_Y;
              temp_ion = ion_density;
               */
              
              if (ion_cell->is_locally_owned())
              {
              if (ion_cell != ion_endc && i < 12) {
                
                double E_x,E_y,ions;
                E_x = Field_X(ion_local_dof_indices[i % 3]);  //Mi servono prese da evaluate electric field
                E_y = Field_Y(ion_local_dof_indices[i % 3]);
                ions = ion_density(ion_local_dof_indices[i % 3]);
                f[0] = q0 * E_x / rho * ions;
                f[1] = q0 * E_y / rho * ions;
              }
              }
              
              
              local_rhs(i) -=
              (viscosity * scalar_product(current_velocity_gradients[q], grad_phi_u[i]) -
              current_velocity_divergences[q] * phi_p[i] -
              current_pressure_values[q] * div_phi_u[i] +
              gamma * current_velocity_divergences[q] * div_phi_u[i] +
              current_velocity_gradients[q] * current_velocity_values[q] * phi_u[i] -
              scalar_product(phi_u[i], f))*     // ADDED
              fe_values.JxW(q);
              
            }
          }


          cell->get_dof_indices(local_dof_indices);
        
          const AffineConstraints<double> &constraints_used = use_nonzero_constraints ? nonzero_NS_constraints : zero_NS_constraints;
          
          if (cell->is_locally_owned())
          {
          if (ion_cell != ion_endc){
        	   ion_cell++;
          }
          }

          if (assemble_system)
          {   
          
              constraints_used.distribute_local_to_global(local_matrix,             //In practice this function implements a scatter operation
                                                          local_rhs,
                                                          local_dof_indices,        //Contains the corresponding global indexes
                                                          NS_system_matrix,
                                                          NS_system_rhs);
              constraints_used.distribute_local_to_global(
              local_mass_matrix, local_dof_indices, NS_mass_matrix);
          }
          else
          {
              constraints_used.distribute_local_to_global(
              local_rhs, local_dof_indices, NS_system_rhs);

          
          }
      }
  }
   
  if (assemble_system)
  {
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
    if (assemble_system)
    {
      preconditioner.reset(new BlockSchurPreconditioner(timer,
                                                      gamma,
                                                      viscosity,
                                                      time_NS.get_delta_t(),
                                                      owned_partitioning,
                                                      NS_system_matrix,
                                                      NS_mass_matrix,
                                                      pressure_mass_matrix));
    }
    
    double coeff = 0.0 ;   // to avoid to have a tolerance too small
    if (time_step < 4) {
        coeff = 1e-5;
    } else if (time_step >= 4 && time_step < 10) {
        coeff = 1e-3;
    } else {
        coeff = 1e-2;
    }

    SolverControl solver_control(                                         //Used by iterative methods to determine whether the iteration should be continued
    NS_system_matrix.m(), coeff * NS_system_rhs.l2_norm(), true);

    SolverBicgstab<PETScWrappers::MPI::BlockVector> bicg(solver_control);

    pcout << "   NS_system_matrix frob norm is " << NS_system_matrix.frobenius_norm() << std::endl;
    pcout << "   NS_system_rhs l2 norm is " << NS_system_rhs.l2_norm() << std::endl;
    // The solution vector must be non-ghosted
    bicg.solve(NS_system_matrix, NS_solution_update, NS_system_rhs, *preconditioner);

    pcout << "   solver di NS fatto: " << NS_solution.l2_norm() << std::endl; 

    const AffineConstraints<double> &constraints_used =
    use_nonzero_constraints ? nonzero_NS_constraints : zero_NS_constraints;
    constraints_used.distribute(NS_solution_update);

    return {solver_control.last_step(), solver_control.last_value()};
}





//-----------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::solve_navier_stokes()
{
	evaluate_electric_field();   //Serve ad assembly_NS, crea Field_X e Field_Y

	

	// if (step_number == 1) {   
	// 	assemble_navier_stokes(true);
	// 	solve_nonlinear_navier_stokes_step(true);
	// 	NS_solution = NS_newton_update;
	// 	nonzero_NS_constraints.distribute(NS_solution);
	// 	double current_res = NS_newton_update.block(0).linfty_norm();
	// 	std::cout << "The residual of  the initial guess is " << current_res << std::endl;
	// }

  NS_solution_update = 0;
  // Only use nonzero constraints at the very first time step (Così i valori nei constraints non vengono modificati)
  bool apply_nonzero_constraints = (time_NS.get_timestep() == 1);  //Capire se usare step_number o la classe Time

  // We have to assemble the LHS for the initial two time steps:
  // once using nonzero_constraints, once using zero_constraints.
  bool assemble_system = (time_NS.get_timestep() < 3);

  pcout << "   ASSEMBLE NAVIER STOKES SYSTEM ..."  << std::endl;
  assemble_NS(apply_nonzero_constraints, assemble_system);


  pcout << "   SOLVE NAVIER STOKES SYSTEM ..."  << std::endl;
  auto state = solver_NS(apply_nonzero_constraints, assemble_system, time_NS.get_timestep());


  // Note we have to use a non-ghosted vector to do the addition.
  PETScWrappers::MPI::BlockVector tmp;
  tmp.reinit(owned_partitioning, mpi_communicator);   //We do this since present solution is a ghost vector so is read-only
  tmp = NS_solution;
  tmp += NS_solution_update;
  NS_solution = tmp;

  pcout << std::scientific << std::left << "   GMRES_ITR = " << std::setw(3)
      << state.first << "   GMRES_RES = " << state.second << std::endl;

  pcout << "   L2 norm of the present solution: " << NS_solution.l2_norm() << std::endl;

  //CI VIENE UNA SOL NULLA !!! La norma L2 in solver_NS viene zero, qua fuori e-20


  //Ho ottenuto la nuova soluzione di NS, spezziamola ora per essere usata dalla parte elettrica


  //Questa parte con Vel_X e Y capire se effettivamnete crearli, per pressure possiamo accedervi in modo più intelligente essendo il 
  //secondo blocco come avevamo fatto anche noi in passato. Il problema è che vengono usati Vel_X e Y in DD. Sicurametne c'è modo più intelligente per accedervi

	pcout << "   Recovering velocity and pressure values for output... " << std::endl;

  //????
  PETScWrappers::MPI::Vector temp_X;  //NB dof handler di DD !!! va così, come al primo giro
  PETScWrappers::MPI::Vector temp_Y;

  temp_X.reinit(locally_owned_dofs,  mpi_communicator);
  temp_Y.reinit(locally_owned_dofs,  mpi_communicator);

  pressure.reinit(owned_partitioning_p, mpi_communicator);

	const unsigned int dofs_per_cell = 4;
	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

	const unsigned int dofs_per_NS_cell = 22;
	std::vector<types::global_dof_index> NS_local_dof_indices(dofs_per_NS_cell);

	auto cell = dof_handler.begin_active();    //Numero dof tra DD (ordine 1) e NS (ordine velocità e pressione diversi) sono diversi
	auto NS_cell = NS_dof_handler.begin_active();

	const auto endc = dof_handler.end();
	const auto NS_endc = NS_dof_handler.end();

	double vel_max = 0.;  //vel_max non serve da nessuna parte credo, capire se serve calcolarla

	while (cell != endc && NS_cell != NS_endc) {
    if (cell->is_locally_owned() && NS_cell->is_locally_owned())
    {

        cell->get_dof_indices(local_dof_indices);
        NS_cell->get_dof_indices(NS_local_dof_indices);

        for (unsigned int k = 0; k < dofs_per_cell; ++k) {

          const unsigned int ind = local_dof_indices[k];

          temp_X(ind) = NS_solution[NS_local_dof_indices[3*k]]; // Not so sure about this...
          temp_Y(ind) = NS_solution[NS_local_dof_indices[3*k+1]]; // ... or this...
          pressure(ind) = NS_solution[NS_local_dof_indices[3*k+2]]; // ... or even this
          // But they all seem to work in the output!

          // vel_max = std::max(vel_max,Vel_X(ind));
          vel_max = std::max(vel_max, static_cast<double>(Vel_X(ind)));    //static_cast altrimenti Vel_X(ind) è un tipo dealii::PETScWrappers::internal::VectorReference

        }
    }

		++cell;
		++NS_cell;

	}

  temp_X.compress(VectorOperation::insert);
  temp_Y.compress(VectorOperation::insert);
  
  pressure.compress(VectorOperation::insert);
  
  Vel_X = temp_X;
  Vel_Y = temp_Y;

	// cout << "Estimating thrust..." << endl; estimate_thrust();
}


//############ - RUN AND OUTPUT RESULTS - #######################################################################################################

template <int dim>
void CompleteProblem<dim>::run()
{

  // fix the constants
  const double Re = m_data.geometrical_parameters.emitter_radius[simulation_tag]; // per adesso un singolo valore
  const double E_ON = m_data.drift_diffusion.physical_parameters.E_ON;
  const double Ve = m_data.drift_diffusion.physical_parameters.Ve;

  const double N_0 = m_data.drift_diffusion.physical_parameters.stratosphere ? 2.2e-3 : 0.5e-3; // [m^-3] ambient ion density 
  const double p_amb = m_data.drift_diffusion.physical_parameters.stratosphere ? 5474 : 101325; 
  const double T = m_data.drift_diffusion.physical_parameters.stratosphere ? 217. : 303.; // [K] fluid temperature

  const double delta = p_amb/101325*298/T;                                       
     
  const double eps = 1.; // wire surface roughness correction coefficient
  const double Ep = E_ON*delta*eps*(1+0.308/std::sqrt(Re*1.e+2*delta));
  const double Ri = Ep/E_ON*Re; // [m] ionization radius
  const double Vi = Ve - Ep*Re*log(Ep/E_ON); // [V] voltage on ionization region boundary


  pcout << std::endl;
	pcout << "From Peek's law:  " << std::endl<<std::endl;
	pcout << "	the field at emitter surface is " << Ep*1e-6 << " [MV/m] " <<std::endl;
	pcout << "	the ionization radius is "<< Ri*1e+3 << " [mm]" << std::endl;
	pcout << " 	the ionization potential is " << Vi*1e-3  << " [kV]" << std::endl;
	pcout << std::endl;


  
  pcout << "   SETUP POISSON PROBLEM ... "<< std::endl;
	setup_poisson();

  
  pcout << "   ASSEMBLE MASS AND LAPLACE POISSON MATRICES ... "<< std::endl;
  assemble_poisson_laplace_matrix();
  assemble_poisson_mass_matrix();


  pcout << "   INITIALIZE POTENTIAL ... "<< std::endl;
  //initialize_potential();
  solve_homogeneous_poisson(); // per avere una buona guess per newton
	
  
  pcout << "   SETUP DRIFT-DIFFUSION PROBLEM ... "<< std::endl;
	setup_drift_diffusion(/*re-initialize densities = */ true); // sta cosa del bool dentro è strana
                                                              // nel codice originale, false, non viene mai messo
  
  pcout << "   ASSEMBLE MASS DRIFT DIFFUSION MATRIX ... "<< std::endl;
  assemble_drift_diffusion_mass_matrix();


  pcout << "   INITIALIZE OLD_ION_DENSITY ... "<< std::endl;
	VectorTools::interpolate(mapping, dof_handler , Functions::ConstantFunction<dim>(N_0), old_ion_density); //N_0 deriva da stratosphere bool
  

  pcout << "   SETUP NAVIER STOKES PROBLEM ... "<< std::endl;
	setup_NS();
  

  pcout << "   STORE INITIAL CONDITIONS ... "<< std::endl<<std::endl;
	output_results(0);


  // SET ERRORS AND TOLERANCES
  step_number = 0; 

  const unsigned int max_steps = 500;
	const double tol = 5.e-9;

	const unsigned int max_it = 1e+3; // riferito al gummel algorithm

	const double time_tol = 5.e-3; //1.e-3;
	double time_err = 1. + time_tol;


  pcout << "   START COMPLETE PROBLEM ... "<< std::endl;

	while (step_number < max_steps && time_err > time_tol)
	{

              ++step_number;

              pcout << "   DD-NS ITERATION:  "<<step_number<< std::endl;
              // Faster time-stepping (adaptive time-stepping would be MUCH better!)
              // timestep varia, diventando sempre più grande, col procedere delle iterazioni temporali
              if (step_number % 40 == 1 && step_number > 1 && timestep < 1.e-3)
                timestep*= 10.;


              const double gummel_tol = 5.e-3; // il nostro codice arriva a 3.7e-3, prima era settata a 1.e-4!!!
              double err = gummel_tol + 1.;
              unsigned int it = 0;

              eta = old_ion_density; // basically eta = N_0 function

              PETScWrappers::MPI::Vector previous_density(locally_owned_dofs, mpi_communicator); //non-ghosted

              pcout << "   GUMMEL ALGORITHM ... "<< std::endl; 
            
              while (err > gummel_tol && it < max_it) {

                    pcout << "   GUMMEL ITERATION: "<< it+1<<std::endl<<std::endl;

                    solve_nonlinear_poisson(max_it,tol); // UPDATE potential AND eta
                    
                    previous_density = ion_density; // save the previously computed ion_density

                    perform_drift_diffusion_fixed_point_iteration_step(); // UPDATE ion_density
                    pcout << "   FATTO 3 "<< std::endl;    

                    previous_density -= ion_density;

                    err = previous_density.linfty_norm()/ion_density.linfty_norm();
                    pcout << "   FATTO 4 "<< std::endl;

                    eta = ion_density;

                    pcout <<"   ERROR: " << err <<std::endl<<std::endl;
                    
                    it++; // questo ciclo lo fa ma l'errore si abbassa poco alla volta arriva a 3.7e-3, in circa 25 iter
              }

              if (it >= max_it){
                pcout << "WARNING! DD achieved a relative error " << err << " after " << it << " iterations" << std::endl;
              }

            
              previous_density = old_ion_density;
              previous_density -= ion_density;
              time_err = previous_density.linfty_norm()/old_ion_density.linfty_norm();

              pcout << "   Density change from previous time-step is: " << time_err*100. << " %" << std::endl<<std::endl;

              old_ion_density = ion_density;

              if (step_number % 40 == 1){ // NS solution update every 40 timesteps

                pcout << "   START NAVIER STOKES PROBLEM ... "<< std::endl; // FINO A QUA FUNZIONA 
                solve_navier_stokes();

                }

              output_results(step_number);

  }

    //pcout << " 	Elapsed CPU time: " << timer.cpu_time()/60. << " minutes.\n" << std::endl << std::endl;

}


//-----------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void CompleteProblem<dim>::output_results(const unsigned int cycle)
  {

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);

    
    // Base directory for output
    std::string base_directory = "../output/results";

    // Directory to store the results of this simulation
    std::string output_directory = base_directory + "/NS_DD_Simulation_"+std::to_string(simulation_tag)+"/";

    // Ensure the output directory is created (if it doesn't exist)
    if (!std::filesystem::exists(output_directory))
    {
        std::filesystem::create_directory(output_directory);
    }

    data_out.add_data_vector(potential, "phi");
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

FullMatrix<double> compute_triangle_matrix(const Point<2> a, const Point<2> b, const Point<2> c, const double alpha12, const double alpha23, const double alpha31, const data_struct& m_data)
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

    //fix the constants
  const double kB = m_data.drift_diffusion.physical_parameters.kB;
  const double q0 = m_data.drift_diffusion.physical_parameters.q0;
  const double mu0 = m_data.drift_diffusion.physical_parameters.mu0;
  const double T = m_data.drift_diffusion.physical_parameters.stratosphere ? 217. : 303.; // [K] fluid temperature
  const double p_amb = m_data.drift_diffusion.physical_parameters.stratosphere ? 5474 : 101325; 
  const double delta = p_amb/101325*298/T; // 

  const double mu = mu0 * delta; // scaled mobility from Moseley                   
  const double V_E = kB * T / q0; // [V] ion temperature                           
  const double D = mu * V_E; //

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
//----------------------------------------------------------------------------------------------------------------------------------------------------
/*Tensor<1,2> get_emitter_normal(const Point<2> a) { //serve per calcolare il thrust

	Tensor<1,2> normal;

	normal[0] = a[0] - X;
	normal[1] = a[1];

	const double norm = std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]);

	return normal/norm;
}*/