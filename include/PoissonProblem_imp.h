namespace {

  using namespace dealii;
  
  //------------------------------------------------------------------------------------------------------------------------------

  template <int dim>
  void PoissonProblem<dim>::run(const unsigned int max_iter, const double toll) 
  {
  
    unsigned int counter = 0; 

    pcout << " -- START NEWTON METHOD --" << std::endl;
    pcout << "    SETUP SYSTEM ... ";
    setup_system();
    pcout << " done! "  << std::endl;
  
    pcout << "    INITIALIZE POTENTIAL AND DENSITIES ... ";
    initialize_current_solution();
    //compute_densities();
    output_results(counter);
    pcout << " done! "  << std::endl;
  
    pcout << "    BUILD LAPLACE - MASS MATRICES ... ";
    assemble_laplace_matrix();
    assemble_mass_matrix();
    pcout << " done! "  << std::endl;

    double increment_norm = std::numeric_limits<double>::max(); 

    pcout << "    Initial Increment Norm: " << increment_norm << std::endl;
  

    while(counter < max_iter && increment_norm > toll){
      counter ++; 
      compute_densities();
      
      pcout << "    NEWTON ITERATION NUMBER: "<< counter  <<std::endl;
      pcout << "    Assemble System Matrix ... ";
      assemble_system_matrix();


      std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;
      
      PETScWrappers::MPI::Vector temp00(locally_owned_dofs, mpi_communicator); //non-ghosted auxiliary vector
      temp00=newton_update;

      VectorTools::interpolate_boundary_values(mapping, dof_handler,1, Functions::ZeroFunction<dim>(), emitter_boundary_values);
      MatrixTools::apply_boundary_values(emitter_boundary_values, system_matrix, temp00, system_rhs);
   
      VectorTools::interpolate_boundary_values(mapping, dof_handler,2, Functions::ZeroFunction<dim>(), collector_boundary_values);
      MatrixTools::apply_boundary_values(collector_boundary_values, system_matrix, temp00, system_rhs);
      
      newton_update=temp00;

      pcout << " done! "  << std::endl;
      
      pcout << "    Solve System ... ";
      solve(); // dentro c'e anche il clamping
      pcout << " done! "  << std::endl;
      
      increment_norm = newton_update.l2_norm();
      pcout << "    Update Increment: " << increment_norm << std::endl << std::endl;

      pcout << "    OUTPUT RESULT "<< std::endl;
      output_results(counter);

      
    }
  

    if(counter == max_iter){
      pcout << "    WARNING! MAX NUMBER OF ITERATIONS REACHED!" << std::endl;
    }

    pcout << " -- END NEWTON METHOD -- "<< std::endl;
    
    //stop the timer to see the elapsed time
    timer.stop();
 
    pcout << "   Elapsed CPU time: " << timer.cpu_time() << " seconds."<<std::endl;
    pcout << "   Elapsed wall time: " << timer.wall_time() << " seconds."<<std::endl;

    // reset timer for the next thing it shall do
    timer.reset();
  
  }

  
  //-----------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  void PoissonProblem<dim>::setup_system()
  {

    dof_handler.distribute_dofs(fe);

    // pcout << "   Number of active cells:       "
    // 	  << triangulation.n_global_active_cells() << std::endl
    // 	  << "   Number of degrees of freedom: " << dof_handler.n_dofs()
    // 	  << std::endl;


    // INDEX SETS INITIALIZATION
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

    // PETSC VECTORS DECLARATIONS 
    current_solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);  //ghosted
    newton_update.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);     //ghosted
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);                               //non ghosted

    electron_density.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);  // ghosted
    hole_density.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);      // ghosted
 
    // ZERO CONSTRAINTS FOR NEWTON
    zero_constraints.clear();
    zero_constraints.reinit(locally_relevant_dofs);
    // VectorTools::interpolate_boundary_values(dof_handler, 1, Functions::ZeroFunction<dim>(), zero_constraints);
    // VectorTools::interpolate_boundary_values(dof_handler, 2, Functions::ZeroFunction<dim>(), zero_constraints);
    zero_constraints.close();

    // DYNAMIC SPARSITY PATTERN    
    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, zero_constraints, false); 
    SparsityTools::distribute_sparsity_pattern(dsp,
					       dof_handler.locally_owned_dofs(),
					       mpi_communicator,
					       locally_relevant_dofs);

  
    system_matrix.clear(); // store the matrix that will be solved in the newton iterations
    system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

    laplace_matrix.clear(); //store laplace matrix
    laplace_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

    mass_matrix.clear();  //store mass matrix
    mass_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

    density_matrix.clear(); // store the term: M(n+p)q0/V_TH
    density_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);
  
    // pcout << " End of setup_system "<< std::endl<<std::endl;

  }

  
  //------------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  PoissonProblem<dim>::PoissonProblem(parallel::distributed::Triangulation<dim> &tria, const data_struct& d)
    : mpi_communicator(MPI_COMM_WORLD)
    , m_data(d)
    , triangulation(tria)
    , fe(m_data.drift_diffusion.numerical_parameters.FE_degree)
    , dof_handler(tria)
    , pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
  {}


  //------------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  void PoissonProblem<dim>:: initialize_current_solution(){
    

   double V_TH = m_data.drift_diffusion.physical_parameters.V_TH;
   double D = m_data.drift_diffusion.physical_parameters.D;
   double ni = m_data.drift_diffusion.physical_parameters.ni;
   double A = m_data.drift_diffusion.physical_parameters.A;

    PETScWrappers::MPI::Vector temp;
    temp.reinit(locally_owned_dofs, mpi_communicator);   //non ghosted, serve per imporre i valori delle BCs
    temp = current_solution; //current solution here is zero by default constructor
	  
    std::map<types::global_dof_index, double> boundary_values;
  
    VectorTools::interpolate_boundary_values(dof_handler,
					     1,
					     Functions::ConstantFunction<dim>(V_TH*std::log(D/ni)),
					     boundary_values);

    VectorTools::interpolate_boundary_values(dof_handler,
					     2,
					     Functions::ConstantFunction<dim>(-V_TH*std::log(A/ni)),
					     boundary_values);
  
    for (auto &boundary_value : boundary_values){
      temp(boundary_value.first) = boundary_value.second;
    }
  
    temp.compress(VectorOperation::insert); //giusto insert, add non funziona
    current_solution = temp;
    // pcout << " End of initialization_current_solution "<< std::endl<<std::endl;


  }
  
  //-----------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  void PoissonProblem<dim>:: compute_densities(){   
    
    double V_TH = m_data.drift_diffusion.physical_parameters.V_TH;
    double ni = m_data.drift_diffusion.physical_parameters.ni;
    // in this function we compute the densities of electrons and holes starting from current solution that is the potential
    
    PETScWrappers::MPI::Vector temp_elec; //for elec density
    PETScWrappers::MPI::Vector temp_hole; //for hole density

    temp_elec.reinit(locally_owned_dofs, mpi_communicator); //non-ghosted
    temp_hole.reinit(locally_owned_dofs, mpi_communicator); //non-ghosted

    temp_elec = current_solution;
    temp_hole = current_solution;
  

    //double check = 0;
    for (unsigned int i = temp_elec.local_range().first; i < temp_elec.local_range().second; ++i){ 

      temp_elec[i] = ni*std::exp( temp_elec[i]/V_TH); // electrons
      temp_hole[i] = ni*std::exp(-temp_hole[i]/V_TH); // holes
    
      //check = temp_elec[i]*temp_hole[i]; //deve essere 10^20 --> lo stampa giusto
      //pcout << "check density: " <<check << std::endl;

    }
  
    temp_elec.compress(VectorOperation::insert);
    temp_hole.compress(VectorOperation::insert);
  
    electron_density = temp_elec;
    hole_density = temp_hole;
    
    // pcout << " The L2 norm of electorn density is: "<< electron_density.l2_norm()<< std::endl;
    // pcout << " The L_INF norm of electron density is: "<< electron_density.linfty_norm()<< std::endl;

    // pcout << " The L2 norm of hole density is: "<< hole_density.l2_norm()<< std::endl;
    // pcout << " The L_INF norm of hole density is: "<< hole_density.linfty_norm()<< std::endl;

    // pcout << " End of compute densities "<< std::endl<<std::endl;

  }

  //---------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  void PoissonProblem<dim>::assemble_laplace_matrix()
  {
    const QTrapezoid<dim> quadrature_formula;

    laplace_matrix = 0;

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
        zero_constraints.distribute_local_to_global(cell_matrix,
                                                    local_dof_indices,
                                                    laplace_matrix);
      }
    }

    laplace_matrix.compress(VectorOperation::add);

    // pcout << " The L_INF norm of the laplace matrix is "<<laplace_matrix.linfty_norm() <<std::endl;
    // pcout << " The L_FROB norm of the laplace matrix is "<<laplace_matrix.frobenius_norm() <<std::endl<<std::endl;
    // pcout << " End of Assembling Laplce matrix "<< std::endl<<std::endl;
  }

  //------------------------------------------------------------------------------------------------------------------------------

  template <int dim>
  void PoissonProblem<dim>::assemble_mass_matrix()
  {
    const QTrapezoid<dim> quadrature_formula;

    mass_matrix = 0;

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
	  zero_constraints.distribute_local_to_global(cell_matrix,
						      local_dof_indices,
						      mass_matrix );
	}
    }

    mass_matrix.compress(VectorOperation::add);

    // pcout << " The L_INF norm of the mass matrix is "<<mass_matrix.linfty_norm() <<std::endl;
    // pcout << " The L_FROB norm of the mass matrix is "<<mass_matrix.frobenius_norm() <<std::endl<<std::endl;
    // pcout << " End of Assembling Mass matrix "<< std::endl<<std::endl;

  }
  
  //----------------------------------------------------------------------------------------------------------------------------  
  template <int dim>
  void PoissonProblem<dim>::assemble_system_matrix()
  {

    double eps_r = m_data.drift_diffusion.physical_parameters.eps_r;
    double eps_0 = m_data.drift_diffusion.physical_parameters.eps_0;
    double q0 = m_data.drift_diffusion.physical_parameters.q0;
    double V_TH = m_data.drift_diffusion.physical_parameters.V_TH;

    //BUILDING SYSTEM MATRIX
    system_matrix = 0;
    density_matrix = 0;
  
    PETScWrappers::MPI::Vector temp_1; //temporary non-ghosted vector number 1
    PETScWrappers::MPI::Vector temp_2; //temporary non-ghosted vector number 2
  
    temp_1.reinit(locally_owned_dofs, mpi_communicator);
    temp_2.reinit(locally_owned_dofs, mpi_communicator);

    temp_1 = electron_density; //temp_1 store electron_density "n"
    temp_2 = hole_density;     //temp_2 store hole_density "p"

    double new_value = 0;
  
    // generate the term:  (n+p)*MASS_MAT   lumped version stored in density_matrix

    for (auto iter = locally_owned_dofs.begin(); iter != locally_owned_dofs.end(); ++iter){ 

      new_value = mass_matrix(*iter, *iter) * (temp_1[*iter] + temp_2[*iter]);
      density_matrix.set(*iter,*iter,new_value);

    }
  
    density_matrix.compress(VectorOperation::insert);
  
    // pcout << " The L_INF norm of the density matrix is "<<density_matrix.linfty_norm() <<std::endl;
    // pcout << " The L_FROB norm of the density matrix is "<<density_matrix.frobenius_norm() <<std::endl<<std::endl;

    system_matrix.add(eps_r * eps_0, laplace_matrix); //ho checkato che passi giusto.  This term is: SYS_MAT = SYS_MAT +  eps*A
    system_matrix.add(q0 / V_TH, density_matrix);   // SYS_MAT = SYS_MAT + q0/V_TH * (n+p)*MASS_MAT

  
    // pcout << " The L_INF norm of the assembled matrix is "<<system_matrix.linfty_norm() <<std::endl;
    // pcout << " The L_FROB norm of the assembled matrix is "<<system_matrix.frobenius_norm() <<std::endl<<std::endl;

  
    // BUILDING SYSTEM RHS
    system_rhs = 0;

    temp_1.add(-1., temp_2);   // temp_1 = temp_1 - temp_2 -->  temp_1 = n - p

    VectorTools::interpolate(mapping, dof_handler, DopingValues<dim>(m_data), temp_2);  //temp_2 = N ; where N is 1e+22 on the left side and 1e-22 on the right

    temp_1.add(-1., temp_2);      // temp_1 = n -p -N

    // basically: temp_1 = (n -p -N)
    mass_matrix.vmult(temp_2,temp_1);  // temp_2 = MASS*(n-p-N)
    system_rhs.add(-q0, temp_2);       // SYS_RHS = -q0*MASS*(n-p-N)

    temp_1 = current_solution;                    // temp_1 = phi 
    laplace_matrix.vmult(temp_2, temp_1);         // temp_2 = A*phi
    system_rhs.add(- eps_r * eps_0, temp_2);       //SYS_RHS = SYS_RHS - eps*A*phi

    // zero_constraints.distribute(system_rhs);
  }
  //----------------------------------------------------------------------------------------------------------------------------

  template <int dim>
  void PoissonProblem<dim>::solve()
  {
    double V_TH = m_data.drift_diffusion.physical_parameters.V_TH;

    PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator);
    SolverControl sc_p(dof_handler.n_dofs(), 1e-10);     
    PETScWrappers::SparseDirectMUMPS solverMUMPS(sc_p);     
    solverMUMPS.solve(system_matrix, temp, system_rhs);

    //CLAMPING  
    for (auto iter = locally_owned_dofs.begin(); iter != locally_owned_dofs.end(); ++iter){ 
  
      if (temp[*iter] < -V_TH) { temp[*iter] = -V_TH; }
      else if (temp[*iter] > V_TH) { temp[*iter] = V_TH; }
      
    }

    temp.compress(VectorOperation::insert);
    //zero_constraints.distribute(temp);
    newton_update = temp;
    temp = current_solution;
    temp.add(1.0, newton_update);
    current_solution = temp;

  
  }
  
  //----------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  void PoissonProblem<dim>::output_results(const unsigned int cycle)
  {

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(current_solution, "phi");
    data_out.add_data_vector(electron_density, "n");
    data_out.add_data_vector(hole_density,     "p");

    PETScWrappers::MPI::Vector doping;     //store the term N, that is 1e+22 on the left side and 1e-22 on the right
    doping.reinit(locally_owned_dofs, mpi_communicator);
    doping = 0;
    VectorTools::interpolate(mapping, dof_handler, DopingValues<dim>(m_data), doping); // We interpolate the previusly
                                                                                 // created vector with the values
                                                                                 // of Doping provided by DopingValues

    data_out.add_data_vector(doping, "doping");

    // Base directory for output
    std::string base_directory = "../output/results";

    // Directory to store the results of this simulation
    std::string output_directory = base_directory + "/NLP_Simulation/";

    // Ensure the output directory is created (if it doesn't exist)
    if (!std::filesystem::exists(output_directory))
    {
        std::filesystem::create_directory(output_directory);
    }
    

    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");
    data_out.build_patches();
    data_out.write_vtu_with_pvtu_record(output_directory, "solution", cycle, mpi_communicator, 2, 1);



  }

}
