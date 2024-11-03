namespace {

  using namespace dealii;
  
  //----------------------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  void DriftDiffusion<dim>::run(const unsigned int max_iter, const double dd_toll) 
  {
  
    unsigned int cycle_drift_diffusion = 0;       
    double hole_tol = dd_toll;
    double electron_tol = dd_toll;

    double hole_err = hole_tol + 1.;
    double electron_err = electron_tol + 1.;

    unsigned int max_newton_iterations = 200;
    double tol_newton = 1e-10;

    pcout << " -- START OF DRIFT DIFFUSION PROBLEM --" <<std::endl;
    pcout << "    SETUP SYSTEM ... ";
    setup_system();
    PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator);
    pcout << " done! "  << std::endl;
    

    pcout << "    INITIALIZE POTENTIAL AND DENSITIES ... ";
    initialization();
    pcout << " done! "  << std::endl;
    
    output_results(cycle_drift_diffusion); //for initial condition
    

    pcout << "    BUILD LAPLACE - MASS MATRICES ... ";
    assemble_laplace_matrix();
    assemble_mass_matrix();
    pcout << " done! "  << std::endl;
    

    while ( (hole_err > hole_tol || electron_err > electron_tol) && cycle_drift_diffusion < max_iter)
    {

    ++cycle_drift_diffusion;

    pcout<< "   CYCLE NUMBER " <<cycle_drift_diffusion<<"  OF DRIFT DIFFUSION"<<std::endl;
    
    pcout<< "   NEWTON POISSON" <<std::endl;
    cycle_newton_poisson(max_newton_iterations, tol_newton);



    pcout<< "   Assemble drift diffusion matrix" <<std::endl;
    assemble_drift_diffusion_matrix();
    /*
    if(cycle_drift_diffusion==1){
      std::ofstream outFile("nuove_bc_assemble.dat");
      hole_matrix.print(outFile);
      outFile.close();
    }*/
/*
    PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator);
    temp = hole_density;
    
    for (unsigned int i = temp.local_range().first; i < temp.local_range().second; ++i){ 

    std::cout<< temp[i]<<std::endl;
 

    }*/
   
 
    //pcout<<"   Apply BCs"<<std::endl;
    //apply_drift_diffusion_boundary_conditions();

 
  
  pcout<< "   Solve drift diffusion"<<std::endl;
  solve_drift_diffusion();
/*
    
    temp = hole_density;
    
    for (unsigned int i = temp.local_range().first; i < temp.local_range().second; ++i){ 

    std::cout<< temp[i]<<std::endl;
 

    }*/

      
/*
       for (unsigned int i = rhs_hole_density.local_range().first; i < rhs_hole_density.local_range().second; ++i){ 

    std::cout<< rhs_hole_density[i]<<std::endl;
 

    } */

 
    //stampa matrici
    
    // if(cycle_drift_diffusion==1 ){
    //   std::ofstream outFile("our_toy_hole_1.dat");
    //   hole_matrix.print(outFile);
    //   outFile.close();
    //   std::ofstream outFile2("our_toy_elec_1.dat");
    //   electron_matrix.print(outFile2);
    //   outFile.close();
    // }

    // if(cycle_drift_diffusion==2 ){
    //   std::ofstream outFile("our_toy_hole_2.dat");
    //   hole_matrix.print(outFile);
    //   outFile.close();
    //   std::ofstream outFile2("our_toy_elec_2.dat");
    //   electron_matrix.print(outFile2);
    //   outFile.close();
    // }
    


    
    pcout<< "   Update error for convergence"<<std::endl;

    electron_tol = 1.e-10*old_electron_density.linfty_norm();
    hole_tol = 1.e-10*old_hole_density.linfty_norm();

    temp = hole_density;
    temp.add(-1.,old_hole_density);
    hole_err = temp.linfty_norm();

    temp = electron_density;
    temp.add(-1.,old_electron_density);
    electron_err = temp.linfty_norm();
    
    pcout<<"   Hole density error: "<< hole_err<<std::endl;
    pcout<<"   Electron density error: "<< electron_err<<std::endl<<std::endl;

    output_results(cycle_drift_diffusion);
    

    pcout << "   The L2 norm old hole density: " << old_hole_density.l2_norm() << std::endl;
    pcout << "   The L2 norm old electron density: " << old_electron_density.l2_norm() << std::endl;

    old_hole_density = hole_density;
    old_electron_density = electron_density;
    
    pcout << "   The L2 norm new hole density: " << old_hole_density.l2_norm() << std::endl;
    pcout << "   The L2 norm new electron density: " << old_electron_density.l2_norm() << std::endl;
    
    }


    if(cycle_drift_diffusion == max_iter){
      pcout << "    WARNING! MAX NUMBER OF ITERATIONS REACHED!" << std::endl;
    }

    pcout << " -- END DRIFT DIFFUSION METHOD -- "<< std::endl;
    
    //stop the timer to see the elapsed time
    timer.stop();
 
    pcout << "   Elapsed CPU time: " << timer.cpu_time() << " seconds."<<std::endl;
    pcout << "   Elapsed wall time: " << timer.wall_time() << " seconds."<<std::endl;

    // reset timer for the next thing it shall do
    timer.reset();
  
  }


  //-----------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  void DriftDiffusion<dim>::setup_system()
  {

    dof_handler.distribute_dofs(fe);

    // pcout << "   Number of active cells:       "
    // 	  << triangulation.n_global_active_cells() << std::endl
    // 	  << "   Number of degrees of freedom: " << dof_handler.n_dofs()
    // 	  << std::endl;


    // INDEX SETS INITIALIZATION
    locally_owned_dofs = dof_handler.locally_owned_dofs();                           //local dofs
    locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);    //local dofs + ghost dofs
    

    // PETSC VECTORS DECLARATIONS 
    current_solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);  //ghosted
    newton_update.reinit(locally_owned_dofs, mpi_communicator);                            //non-ghosted
    poisson_system_rhs.reinit(locally_owned_dofs, mpi_communicator);                       //non-ghosted

    electron_density.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);  //ghosted
    hole_density.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);      //ghosted

    old_electron_density.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);  // ghosted
    old_hole_density.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);      // ghosted

    rhs_electron_density.reinit(locally_owned_dofs, mpi_communicator);  // non-ghosted
    rhs_hole_density.reinit(locally_owned_dofs, mpi_communicator);      // non-ghosted
    

    // ZERO_CONSTRAINTS FOR NEWTON POISSON PROBLEM
    zero_constraints.clear();
    zero_constraints.reinit(locally_relevant_dofs);
    //VectorTools::interpolate_boundary_values(dof_handler, 1, Functions::ZeroFunction<dim>(), zero_constraints); 
    //VectorTools::interpolate_boundary_values(dof_handler, 2, Functions::ZeroFunction<dim>(), zero_constraints); 
    zero_constraints.close(); 


    
    //DENSITY CONSTRAINTS FOR DRIFT DIFFUSION SYSTEM
   
   // density_constraints.clear();
    //density_constraints.reinit(locally_relevant_dofs);   
    //density_constraints.close();
    
    // tutto nuovo @
    elec_constraints.clear();
    elec_constraints.reinit(locally_relevant_dofs);
    VectorTools::interpolate_boundary_values(dof_handler, 1, Functions::ConstantFunction<dim>(N1), elec_constraints); 
    VectorTools::interpolate_boundary_values(dof_handler, 2, Functions::ConstantFunction<dim>(N2), elec_constraints); 
    elec_constraints.close();

    hole_constraints.clear();
    hole_constraints.reinit(locally_relevant_dofs);
    VectorTools::interpolate_boundary_values(dof_handler, 1, Functions::ConstantFunction<dim>(P1), hole_constraints); 
    VectorTools::interpolate_boundary_values(dof_handler, 2, Functions::ConstantFunction<dim>(P2), hole_constraints); 
    hole_constraints.close();


    // DYNAMIC SPARSITY PATTERN AND POISSON MATRICES 
    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, zero_constraints, false); 
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               dof_handler.locally_owned_dofs(),
                                               mpi_communicator,
                                               locally_relevant_dofs);

  
    poisson_system_matrix.clear(); //store the matrix used to solve nlpoisson
    poisson_system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

    laplace_matrix.clear(); //store laplace matrix
    laplace_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

    mass_matrix.clear();  //store mass matrix
    mass_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

    density_matrix.clear(); // store the term: M(n+p)q0/V_TH
    density_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);
    

    // DYNAMIC SPARSITY PATTERN AND DRIFT DIFFUSION MATRICES
    /*
    DynamicSparsityPattern dsp_dd(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp_dd, density_constraints, false);  

    // SparsityTools::distribute_sparsity_pattern(dsp_dd,
    //                                           dof_handler.locally_owned_dofs(),
    //                                           mpi_communicator,
    //                                           locally_relevant_dofs);

    // hole_matrix.clear(); //store holes density matrix
    // hole_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp_dd,  mpi_communicator);

    // electron_matrix.clear();// store electron density matrix
    // electron_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp_dd,  mpi_communicator);

    // //pcout << "   End of setup_system "<< std::endl<<std::endl;

    DynamicSparsityPattern elec_dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, elec_dsp, elec_constraints, false);
    SparsityTools::distribute_sparsity_pattern(elec_dsp,
                                            dof_handler.locally_owned_dofs(),
                                            mpi_communicator,
                                            locally_relevant_dofs);

    DynamicSparsityPattern hole_dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, hole_dsp, hole_constraints, false);
    SparsityTools::distribute_sparsity_pattern(hole_dsp,
                                            dof_handler.locally_owned_dofs(),
                                            mpi_communicator,
                                            locally_relevant_dofs);

    hole_matrix.clear(); //store holes density matrix
    hole_matrix.reinit(locally_owned_dofs, locally_owned_dofs, hole_dsp,  mpi_communicator);

    electron_matrix.clear();// store electron density matrix
    electron_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp_dd,  mpi_communicator);
*/
    //pcout << "   End of setup_system "<< std::endl<<std::endl;

    DynamicSparsityPattern elec_dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, elec_dsp, elec_constraints, false);
    SparsityTools::distribute_sparsity_pattern(elec_dsp,
                                            dof_handler.locally_owned_dofs(),
                                            mpi_communicator,
                                            locally_relevant_dofs);

    DynamicSparsityPattern hole_dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, hole_dsp, hole_constraints, false);
    SparsityTools::distribute_sparsity_pattern(hole_dsp,
                                            dof_handler.locally_owned_dofs(),
                                            mpi_communicator,
                                            locally_relevant_dofs);

    hole_matrix.clear(); //store holes density matrix
    hole_matrix.reinit(locally_owned_dofs, locally_owned_dofs, hole_dsp,  mpi_communicator);

    electron_matrix.clear();// store electron density matrix
    electron_matrix.reinit(locally_owned_dofs, locally_owned_dofs, elec_dsp,  mpi_communicator);                                     

    
  }

  //------------------------------------------------------------------------------------------------------------------------------
  //CONSTRUCTOR
  template <int dim>
  DriftDiffusion<dim>::DriftDiffusion(parallel::distributed::Triangulation<dim> &tria)
    : mpi_communicator(MPI_COMM_WORLD)
    , triangulation(tria)
    , fe(1)
    , dof_handler(tria)
    , pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , mapping()
  {}


  //------------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  void DriftDiffusion<dim>:: initialization(){
  
  // POTENTIAL INITIAL CONDITION
  // null with the right BCs
  /*
  PETScWrappers::MPI::Vector temp(locally_owned_dofs, mpi_communicator); //non ghosted vector, needed for imposing BCs

  temp = current_solution; //current solution here is zero by default constructor
  
  std::map<types::global_dof_index, double> boundary_values;

  VectorTools::interpolate_boundary_values(dof_handler, 1, Functions::ConstantFunction<dim>(V_TH*std::log(D/ni)), boundary_values);
  VectorTools::interpolate_boundary_values(dof_handler, 2, Functions::ConstantFunction<dim>(-V_TH*std::log(A/ni)),boundary_values);

  for (auto &boundary_value : boundary_values){
    temp(boundary_value.first) = boundary_value.second;
  }

  temp.compress(VectorOperation::insert);
  current_solution = temp;

  // DENSITIES INITIAL CONDITIONS
  //electrons
  temp = old_electron_density;
  VectorTools::interpolate_boundary_values(dof_handler,1, Functions::ConstantFunction<dim>(N1), boundary_values);
  VectorTools::interpolate_boundary_values(dof_handler,2, Functions::ConstantFunction<dim>(N2), boundary_values);
  
  for (auto &boundary_value : boundary_values){
    temp(boundary_value.first) = boundary_value.second;
  }

  temp.compress(VectorOperation::insert);
  old_electron_density = temp;


  //holes
  temp = old_hole_density;
  VectorTools::interpolate_boundary_values(dof_handler,1, Functions::ConstantFunction<dim>(P1), boundary_values);
  VectorTools::interpolate_boundary_values(dof_handler,2, Functions::ConstantFunction<dim>(P2), boundary_values);
  
  for (auto &boundary_value : boundary_values){
    temp(boundary_value.first) = boundary_value.second;
  }
  temp.compress(VectorOperation::insert);
  old_hole_density = temp;
  */
  //pcout << "   End of initialization_current_solution "<< std::endl;

  //initialization with full interpolation of the right bcs values
  PETScWrappers::MPI::Vector temp_pot(locally_owned_dofs, mpi_communicator);
  PETScWrappers::MPI::Vector temp_hole(locally_owned_dofs, mpi_communicator);
  PETScWrappers::MPI::Vector temp_elec(locally_owned_dofs, mpi_communicator);

  temp_pot = current_solution;
  temp_hole = old_hole_density;
  temp_elec = old_electron_density;

	VectorTools::interpolate(mapping, dof_handler, PotentialValues<dim>(), temp_pot);
	VectorTools::interpolate(mapping, dof_handler, HoleInitialValues<dim>(), temp_hole);
	VectorTools::interpolate(mapping, dof_handler, ElectronInitialValues<dim>(), temp_elec);
  
  temp_pot.compress(VectorOperation::insert);
  temp_hole.compress(VectorOperation::insert);
  temp_elec.compress(VectorOperation::insert);

  current_solution = temp_pot;
  old_hole_density = temp_hole;
  old_electron_density = temp_elec;
  }
  
  //-----------------------------------------------------------------------------------------------------------------------------
  
  template <int dim>
  void DriftDiffusion<dim>::assemble_laplace_matrix()
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
    // pcout << "   End of Assembling Laplce matrix "<< std::endl<<std::endl;
  }

  //------------------------------------------------------------------------------------------------------------------------------

  template <int dim>
  void DriftDiffusion<dim>::assemble_mass_matrix()
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
    // pcout << "   End of Assembling Mass matrix "<< std::endl<<std::endl;

  }
  
  //----------------------------------------------------------------------------------------------------------------------------  
  template <int dim>
  void DriftDiffusion<dim>::assemble_nonlinear_poisson()
  {

    //BUILDING POISSON SYSTEM MATRIX (for newton)
    poisson_system_matrix = 0;
    density_matrix = 0;
  
    PETScWrappers::MPI::Vector temp_1(locally_owned_dofs, mpi_communicator); //temporary non-ghosted vector number 1
    PETScWrappers::MPI::Vector temp_2(locally_owned_dofs, mpi_communicator); //temporary non-ghosted vector number 2

    temp_1 = old_electron_density; //temp_1 store electron_density "n"
    temp_2 = old_hole_density;     //temp_2 store hole_density "p"

    double new_value = 0;
  
    // generate the term:  (n+p)*MASS_MAT   lumped version stored in density_matrix

    for (auto iter = locally_owned_dofs.begin(); iter != locally_owned_dofs.end(); ++iter){ 

      new_value = mass_matrix(*iter, *iter) * (temp_1[*iter] + temp_2[*iter]);
      density_matrix.set(*iter,*iter,new_value);

    }

    density_matrix.compress(VectorOperation::insert);
  
    //pcout << "   The L_INF norm of the density matrix is: "<<density_matrix.linfty_norm() <<std::endl;
    //pcout << "   The L_FROB norm of the density matrix is: "<<density_matrix.frobenius_norm() <<std::endl<<std::endl;

    poisson_system_matrix.add(eps_r * eps_0, laplace_matrix); // SYS_MAT = SYS_MAT +  eps*A
    poisson_system_matrix.add(q0 / V_TH, density_matrix);     // SYS_MAT = SYS_MAT + q0/V_TH * (n+p)*MASS_MAT


    //pcout << "   The L_INF norm of the poisson system matrix is: "<<poisson_system_matrix.linfty_norm() <<std::endl;
    //pcout << "   The L_FROB norm of the poisson system matrix is: "<<poisson_system_matrix.frobenius_norm() <<std::endl<<std::endl;

  
    // BUILDING SYSTEM RHS
    poisson_system_rhs = 0;

    temp_1.add(-1., temp_2);   // temp_1 = temp_1 - temp_2 -->  temp_1 = n - p

    VectorTools::interpolate(mapping, dof_handler, DopingValues<dim>(), temp_2);  //temp_2 = N ; where N is 1e+22 on the left side and 1e-22 on the right

    temp_1.add(-1., temp_2);      // temp_1 = n -p -N

    // basically: temp_1 = (n -p -N)
    mass_matrix.vmult(temp_2,temp_1);  // temp_2 = MASS*(n-p-N)
    poisson_system_rhs.add(-q0, temp_2);       // SYS_RHS = -q0*MASS*(n-p-N)

    temp_1 = current_solution;                    // temp_1 = phi 
    laplace_matrix.vmult(temp_2, temp_1);         // temp_2 = A*phi
    poisson_system_rhs.add(- eps_r * eps_0, temp_2);       //SYS_RHS = SYS_RHS - eps*A*phi



    //pcout << "   The L_INF norm of the poisson system RHS is: "<<poisson_system_rhs.linfty_norm() <<std::endl;
    //pcout << "   The L2 norm of the poisson system RHS is: " << poisson_system_rhs.l2_norm() << std::endl;
    //pcout << "   End of solve assemble non-linear poisson"<< std::endl<<std::endl;

  }
  //----------------------------------------------------------------------------------------------------------------------------

  template <int dim>
  void DriftDiffusion<dim>::solve_poisson()
  {

    //Apply zero boundary conditions to the whole newton poisson system
    std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

    VectorTools::interpolate_boundary_values(mapping, dof_handler,1, Functions::ZeroFunction<dim>(), emitter_boundary_values);
    MatrixTools::apply_boundary_values(emitter_boundary_values, poisson_system_matrix, newton_update, poisson_system_rhs);

    VectorTools::interpolate_boundary_values(mapping, dof_handler,2, Functions::ZeroFunction<dim>(), collector_boundary_values);
    MatrixTools::apply_boundary_values(collector_boundary_values, poisson_system_matrix, newton_update, poisson_system_rhs);
    

    //Solve poisson system problem
    SolverControl sc_p(dof_handler.n_dofs(), 1e-10);     
    PETScWrappers::SparseDirectMUMPS solverMUMPS(sc_p); 

    solverMUMPS.solve(poisson_system_matrix, newton_update, poisson_system_rhs);

    //Clamping 
    for (auto iter = locally_owned_dofs.begin(); iter != locally_owned_dofs.end(); ++iter){ 
  
      if (newton_update[*iter] < -V_TH) { newton_update[*iter] = -V_TH; }
      else if (newton_update[*iter] > V_TH) { newton_update[*iter] = V_TH; }
      
    }

    newton_update.compress(VectorOperation::insert); 
    
    //Update current solution
    PETScWrappers::MPI::Vector temp;
    temp.reinit(locally_owned_dofs, mpi_communicator);

    temp = current_solution;
    temp.add(1.0, newton_update);


    current_solution = temp;

    // pcout << "L2 norm of the current solution: " << current_solution.l2_norm() << std::endl;
    // pcout << "L_INF norm of the current solution: " << current_solution.linfty_norm() << std::endl;
    // pcout << "   End of solve poisson problem"<< std::endl<<std::endl;
  
  }
  //-----------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  void DriftDiffusion<dim>:: compute_densities(){   
  
    // in this function we compute the densities of electrons and holes starting from current solution that is the potential
    // this function is used inside the newton cycle
    
    //update the densities
    PETScWrappers::MPI::Vector old_temp_elec(locally_owned_dofs, mpi_communicator);
    PETScWrappers::MPI::Vector old_temp_hole(locally_owned_dofs, mpi_communicator);

    old_temp_elec  = current_solution;
    old_temp_hole  = current_solution;

    //double check = 0;
    for (unsigned int i = old_temp_elec.local_range().first; i < old_temp_elec.local_range().second; ++i){ 

    old_temp_elec[i] =ni*std::exp(old_temp_elec[i]/V_TH);
    old_temp_hole[i] =ni*std::exp(-old_temp_hole[i]/V_TH);

    //check = old_temp_elec[i]*old_temp_hole[i]; //deve essere ni^2 = 10^32 --> lo stampa giusto
    //pcout << "check density: " <<check << std::endl;

    }

    old_temp_elec.compress(VectorOperation::insert);
    old_temp_hole.compress(VectorOperation::insert);

    old_electron_density = old_temp_elec;
    old_hole_density = old_temp_hole;
    
    // pcout << " The L2 norm of electorn density is: "<< electron_density.l2_norm()<< std::endl;
    // pcout << " The L_INF norm of electron density is: "<< electron_density.linfty_norm()<< std::endl;

    // pcout << " The L2 norm of hole density is: "<< hole_density.l2_norm()<< std::endl;
    // pcout << " The L_INF norm of hole density is: "<< hole_density.linfty_norm()<< std::endl;

    // pcout << " End of compute densities "<< std::endl<<std::endl;

  }

//--------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void DriftDiffusion<dim>::cycle_newton_poisson(const unsigned int max_iter_newton, 
                                               const double toll_newton){

  unsigned int counter = 0; // it keeps track of newton iteration

  pcout << "   - START NEWTON METHOD - "<< std::endl;

  double increment_norm = std::numeric_limits<double>::max(); // the increment norm is + inf

  pcout << "   Initial Newton Increment Norm dphi: " << increment_norm << std::endl<<std::endl;
  

  while(counter < max_iter_newton && increment_norm > toll_newton){

    pcout << "   NEWTON ITERATION NUMBER: "<< counter +1<<std::endl;
    pcout << "   Assemble System Poisson Matrix"<< std::endl;

    //NB: Mass and Laplace matrices are already build (see run method)
    
    assemble_nonlinear_poisson();     
    solve_poisson();  //clamping on newton update and BCs are inside this method  
    compute_densities();
    increment_norm = newton_update.l2_norm();
    pcout << "   Update Increment: "<<increment_norm<<std::endl<<std::endl;

    counter ++;

    if(counter == max_iter_newton){
      pcout<< "   MAX NUMBER OF NEWTON ITERATIONS REACHED!"<<std::endl;
    }

  }

  // pcout << "   end one cycle of newton method "<< std::endl<<std::endl;

}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
/*
template <int dim>
void DriftDiffusion<dim>::assemble_drift_diffusion_matrix()
{

  //initialization
	rhs_electron_density = 0;
	rhs_hole_density = 0;
  hole_matrix = 0;
	electron_matrix = 0;

  
  const unsigned int vertices_per_cell = 4; // 4 number of dofs per cell, 4 is dofs per cell

  std::vector<types::global_dof_index> local_dof_indices(vertices_per_cell);

  const unsigned int t_size = 3;

  Vector<double> cell_rhs(t_size);
  FullMatrix<double> A(t_size,t_size), B(t_size,t_size), neg_A(t_size,t_size), neg_B(t_size,t_size);

  std::vector<types::global_dof_index> A_local_dof_indices(t_size);
  std::vector<types::global_dof_index> B_local_dof_indices(t_size);


  int cell_index = 0;
  std::cout << "#cell\t\t#nvertex(local)\t\t#nvertex(global)\t\tcoords\t\tvalues\n"<<std::endl;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned()){ // tieni calcolo coordinate e potenziale da 644 puntini 

        A = 0;
        B = 0;
        neg_A = 0;
        neg_B = 0;

        cell_rhs = 0;
        
        cell->get_dof_indices(local_dof_indices); //controlla valori

        // Lexicographic ordering
        const Point<dim> v1 = cell->vertex(0); // "top left"      2            // in realtà top right
        const Point<dim> v2 = cell->vertex(1); // "top right"     3            // in realtà bottom right
        const Point<dim> v3 = cell->vertex(2); // "bottom left"   0           // in realtà top left
        const Point<dim> v4 = cell->vertex(3); // "bottom right"  1          // in realtà bottom left
        
        std::cout << cell_index   << "\t\t   " << 0 << "\t\t   " << local_dof_indices[0] << "\t\t   " << v1[0] << ", " << v1[1] <<"\t\t   "<<current_solution[local_dof_indices[0]]<< "\n";
        std::cout << cell_index   << "\t\t   " << 1 << "\t\t   " << local_dof_indices[1] << "\t\t   " << v2[0] << ", " << v2[1] <<"\t\t   "<<current_solution[local_dof_indices[1]]<< "\n";
        std::cout << cell_index   << "\t\t   " << 2 << "\t\t   " << local_dof_indices[2] << "\t\t   " << v3[0] << ", " << v3[1] <<"\t\t   "<<current_solution[local_dof_indices[2]]<< "\n";
        std::cout << cell_index++ << "\t\t   " << 3 << "\t\t   " << local_dof_indices[3] << "\t\t   " << v4[0] << ", " << v4[1] <<"\t\t   "<<current_solution[local_dof_indices[3]]<< "\n";

        const double u1 = -current_solution[local_dof_indices[0]]/V_TH;
        const double u2 = -current_solution[local_dof_indices[1]]/V_TH;
        const double u3 = -current_solution[local_dof_indices[2]]/V_TH;
        const double u4 = -current_solution[local_dof_indices[3]]/V_TH;

        const double l_alpha = side_length(v1,v4);
        const double l_beta = side_length(v2,v3);
        
        const double alpha21 = (u1 - u2);
        const double alpha42 = (u2 - u4);
        const double alpha34 = (u4 - u3);
        const double alpha13 = (u3 - u1);
        
        const double neg_alpha21 =  - (u1 - u2);
        const double neg_alpha42 =  - (u2 - u4);
        const double neg_alpha34 = - (u4 - u3);
        const double neg_alpha13 = - (u3 - u1);  

        if (l_alpha >= l_beta) { // l_alpha is the longest diagonal: split by beta
              const double alpha23 =  (u3 - u2);
              const double neg_alpha23 = - (u3 - u2);
              //cout << "Alpha 23 is: " << alpha23 << endl;
              
              // Triangle A:
              A= compute_triangle_matrix(v2,v1,v3, alpha21, alpha13, -alpha23, Dp);
              neg_A= compute_triangle_matrix(v2,v1,v3, neg_alpha21, neg_alpha13, -neg_alpha23, Dn);
              
              // Triangle B:
              B = compute_triangle_matrix(v3,v4,v2, alpha34, alpha42, alpha23, Dp);
              neg_B = compute_triangle_matrix(v3,v4,v2, neg_alpha34, neg_alpha42, neg_alpha23, Dn);
              
              // Matrix assemble
              A_local_dof_indices[0] = local_dof_indices[1];
              A_local_dof_indices[1] = local_dof_indices[0];
              A_local_dof_indices[2] = local_dof_indices[2];
              
              B_local_dof_indices[0] = local_dof_indices[2];
              B_local_dof_indices[1] = local_dof_indices[3];
              B_local_dof_indices[2] = local_dof_indices[1];
          
            } else { // l_beta is the longest diagonal: split by alpha
              const double alpha14 = (u4 - u1);
              const double neg_alpha14 = - (u4 - u1);
              //cout << "Alpha 14 is: " << alpha14 << endl;
              
              // Triangle A:
              A = compute_triangle_matrix(v4,v2,v1, alpha42, alpha21, alpha14, Dp);
              neg_A = compute_triangle_matrix(v4,v2,v1, neg_alpha42, neg_alpha21, neg_alpha14, Dn);
              
              // Triangle B:
              B = compute_triangle_matrix(v1,v3,v4, alpha13, alpha34, -alpha14, Dp);
              neg_B = compute_triangle_matrix(v1,v3,v4, neg_alpha13, neg_alpha34, -neg_alpha14, Dn);
              
              A_local_dof_indices[0] = local_dof_indices[3];
              A_local_dof_indices[1] = local_dof_indices[1];
              A_local_dof_indices[2] = local_dof_indices[0];
            
              B_local_dof_indices[0] = local_dof_indices[0];
              B_local_dof_indices[1] = local_dof_indices[2];
              B_local_dof_indices[2] = local_dof_indices[3];
              
				     }
        
				density_constraints.distribute_local_to_global(A, cell_rhs,  A_local_dof_indices, hole_matrix, rhs_hole_density);
				density_constraints.distribute_local_to_global(B, cell_rhs,  B_local_dof_indices, hole_matrix, rhs_hole_density);

				density_constraints.distribute_local_to_global(neg_A, cell_rhs,  A_local_dof_indices, electron_matrix, rhs_electron_density);
				density_constraints.distribute_local_to_global(neg_B, cell_rhs,  B_local_dof_indices, electron_matrix, rhs_electron_density);

        }
		  }
    
    
    hole_matrix.compress(VectorOperation::add);
    electron_matrix.compress(VectorOperation::add);
    
    rhs_hole_density.compress(VectorOperation::add);
    rhs_electron_density.compress(VectorOperation::add);*/

//  INDICI ORIGINALI

template <int dim>
void DriftDiffusion<dim>::assemble_drift_diffusion_matrix()
{

  //initialization
	rhs_electron_density = 0;
	rhs_hole_density = 0;
  hole_matrix = 0;
	electron_matrix = 0;

  
  const unsigned int vertices_per_cell = 4; // 4 number of dofs per cell, 4 is dofs per cell

  std::vector<types::global_dof_index> local_dof_indices(vertices_per_cell);

  const unsigned int t_size = 3;

  Vector<double> cell_rhs(t_size);
  FullMatrix<double> A(t_size,t_size), B(t_size,t_size), neg_A(t_size,t_size), neg_B(t_size,t_size);

  std::vector<types::global_dof_index> A_local_dof_indices(t_size);
  std::vector<types::global_dof_index> B_local_dof_indices(t_size);


  int cell_index = 0;
  // std::cout << "#cell\t\t    #nvertex(local)   \t\t    #nvertex(global)\t\t    coords\n";

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned()){ // tieni calcolo coordinate e potenziale da 644 puntini 

        A = 0;
        B = 0;
        neg_A = 0;
        neg_B = 0;

        cell_rhs = 0;
        
        cell->get_dof_indices(local_dof_indices); //controlla valori

        // Lexicographic ordering
        const Point<dim> v1 = cell->vertex(2); // top left
        const Point<dim> v2 = cell->vertex(3); // top right
        const Point<dim> v3 = cell->vertex(0); // bottom left
        const Point<dim> v4 = cell->vertex(1); // bottom right
        
        /*
        std::cout << cell_index   << "\t\t   " << 2 << "\t\t   " << local_dof_indices[2] << "\t\t   " << v1[0] << ", " << v1[1] << "\n";
        std::cout << cell_index   << "\t\t   " << 3 << "\t\t   " << local_dof_indices[3] << "\t\t   " << v2[0] << ", " << v2[1] << "\n";
        std::cout << cell_index   << "\t\t   " << 0 << "\t\t   " << local_dof_indices[0] << "\t\t   " << v3[0] << ", " << v3[1] << "\n";
        std::cout << cell_index++ << "\t\t   " << 1 << "\t\t   " << local_dof_indices[1] << "\t\t   " << v4[0] << ", " << v4[1] << "\n";
*/
        const double u1 = -current_solution[local_dof_indices[2]]/V_TH;
        const double u2 = -current_solution[local_dof_indices[3]]/V_TH;
        const double u3 = -current_solution[local_dof_indices[0]]/V_TH;
        const double u4 = -current_solution[local_dof_indices[1]]/V_TH;

        const double l_alpha = side_length(v1,v4);
        const double l_beta = side_length(v2,v3);
        
        const double alpha21 = (u1 - u2);
        const double alpha42 = (u2 - u4);
        const double alpha34 = (u4 - u3);
        const double alpha13 = (u3 - u1);
        
        const double neg_alpha21 =  - (u1 - u2);
        const double neg_alpha42 =  - (u2 - u4);
        const double neg_alpha34 = - (u4 - u3);
        const double neg_alpha13 = - (u3 - u1);

        if (l_alpha >= l_beta) { // l_alpha is the longest diagonal: split by beta
              const double alpha23 =  (u3 - u2);
              const double neg_alpha23 = - (u3 - u2);
              //cout << "Alpha 23 is: " << alpha23 << endl;
              
              // Triangle A:
              A= compute_triangle_matrix(v2,v1,v3, alpha21, alpha13, -alpha23, Dp);
              neg_A= compute_triangle_matrix(v2,v1,v3, neg_alpha21, neg_alpha13, -neg_alpha23, Dn);
              
              // Triangle B:
              B = compute_triangle_matrix(v3,v4,v2, alpha34, alpha42, alpha23, Dp);
              neg_B = compute_triangle_matrix(v3,v4,v2, neg_alpha34, neg_alpha42, neg_alpha23, Dn);
              
              // Matrix assemble
              A_local_dof_indices[0] = local_dof_indices[3];
              A_local_dof_indices[1] = local_dof_indices[2];
              A_local_dof_indices[2] = local_dof_indices[0];
              
              B_local_dof_indices[0] = local_dof_indices[0];
              B_local_dof_indices[1] = local_dof_indices[1];
              B_local_dof_indices[2] = local_dof_indices[3];
          
            } else { // l_beta is the longest diagonal: split by alpha
              const double alpha14 = (u4 - u1);
              const double neg_alpha14 = - (u4 - u1);
              //cout << "Alpha 14 is: " << alpha14 << endl;
              
              // Triangle A:
              A = compute_triangle_matrix(v4,v2,v1, alpha42, alpha21, alpha14, Dp);
              neg_A = compute_triangle_matrix(v4,v2,v1, neg_alpha42, neg_alpha21, neg_alpha14, Dn);
              
              // Triangle B:
              B = compute_triangle_matrix(v1,v3,v4, alpha13, alpha34, -alpha14, Dp);
              neg_B = compute_triangle_matrix(v1,v3,v4, neg_alpha13, neg_alpha34, -neg_alpha14, Dn);
              
              A_local_dof_indices[0] = local_dof_indices[1];
              A_local_dof_indices[1] = local_dof_indices[3];
              A_local_dof_indices[2] = local_dof_indices[2];
            
              B_local_dof_indices[0] = local_dof_indices[2];
              B_local_dof_indices[1] = local_dof_indices[0];
              B_local_dof_indices[2] = local_dof_indices[1];
              
				     }
        
				hole_constraints.distribute_local_to_global(A, cell_rhs,  A_local_dof_indices, hole_matrix, rhs_hole_density);
				hole_constraints.distribute_local_to_global(B, cell_rhs,  B_local_dof_indices, hole_matrix, rhs_hole_density);

				elec_constraints.distribute_local_to_global(neg_A, cell_rhs,  A_local_dof_indices, electron_matrix, rhs_electron_density);
				elec_constraints.distribute_local_to_global(neg_B, cell_rhs,  B_local_dof_indices, electron_matrix, rhs_electron_density);

        }
		  }
    
    
    hole_matrix.compress(VectorOperation::add);
    electron_matrix.compress(VectorOperation::add);
    
    rhs_hole_density.compress(VectorOperation::add);
    rhs_electron_density.compress(VectorOperation::add);
   
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void DriftDiffusion<dim>::apply_drift_diffusion_boundary_conditions() 
{
  
  //appy boundary condition for the densities systems to solve
  PETScWrappers::MPI::Vector temp_elec(locally_owned_dofs, mpi_communicator);
  PETScWrappers::MPI::Vector temp_hole(locally_owned_dofs, mpi_communicator);     

  temp_elec = electron_density;
  temp_hole = hole_density;
  
  std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

  VectorTools::interpolate_boundary_values(mapping, dof_handler,1, Functions::ConstantFunction<dim>(N1), emitter_boundary_values);
  MatrixTools::apply_boundary_values(emitter_boundary_values, electron_matrix, temp_elec, rhs_electron_density);

  VectorTools::interpolate_boundary_values(mapping, dof_handler,2, Functions::ConstantFunction<dim>(N2), collector_boundary_values);
  MatrixTools::apply_boundary_values(collector_boundary_values, electron_matrix, temp_elec, rhs_electron_density);

  VectorTools::interpolate_boundary_values(mapping, dof_handler,1, Functions::ConstantFunction<dim>(P1), emitter_boundary_values);
  MatrixTools::apply_boundary_values(emitter_boundary_values, hole_matrix, temp_hole, rhs_hole_density);

  VectorTools::interpolate_boundary_values(mapping, dof_handler,2, Functions::ConstantFunction<dim>(P2), collector_boundary_values);
  MatrixTools::apply_boundary_values(collector_boundary_values, hole_matrix, temp_hole, rhs_hole_density);
  
  electron_density = temp_elec;
  hole_density = temp_hole;
 /* 
  pcout << "   (apply bcs) L_INF norm of the hole matrix:   "<<hole_matrix.linfty_norm() <<std::endl;           
  pcout << "   (apply bcs) L_INF norm of the electron matrix:  "<<electron_matrix.linfty_norm() <<std::endl;

  pcout << "   L_INF norm of the hole RHS:  " << rhs_hole_density.linfty_norm() << std::endl;        
  pcout << "   L_INF norm of the electron RHS:  " << rhs_electron_density.linfty_norm() << std::endl;
  */
  //pcout<<" end of apply boundary conditions for drift diffusion"<<std::endl<<std::endl;

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void DriftDiffusion<dim>::solve_drift_diffusion()
{ 

  PETScWrappers::MPI::Vector temp_elec(locally_owned_dofs, mpi_communicator);
  PETScWrappers::MPI::Vector temp_hole(locally_owned_dofs, mpi_communicator);
  
  SolverControl sc_hole(dof_handler.n_dofs(), 1e-10);
  SolverControl sc_elec(dof_handler.n_dofs(), 1e-10); 
  
  PETScWrappers::SparseDirectMUMPS solverMUMPS_hole(sc_hole); 
  PETScWrappers::SparseDirectMUMPS solverMUMPS_elec(sc_elec); 
  
  solverMUMPS_hole.solve(hole_matrix, temp_hole, rhs_hole_density);
  solverMUMPS_elec.solve(electron_matrix, temp_elec, rhs_electron_density);

  hole_constraints.distribute(temp_hole);
  elec_constraints.distribute(temp_elec);




  temp_elec.compress(VectorOperation::insert);
  temp_hole.compress(VectorOperation::insert);

  hole_density = temp_hole;
  electron_density = temp_elec;

  /*
  pcout << "   L_INF norm of the hole_density: " << hole_density.linfty_norm() << std::endl;
  pcout << "   L_INF norm of the electron_density: " << electron_density.linfty_norm() << std::endl;
  */

  //pcout << "   End of solve drift diffusion"<< std::endl<<std::endl;
    

}
//----------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  void DriftDiffusion<dim>::output_results(const unsigned int cycle)
  {
    // Base directory for output
    std::string base_directory = "../output";

    // Directory to store the results of this simulation
    std::string output_directory = base_directory + "/PN_Simulation/";

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(current_solution, "phi");
    data_out.add_data_vector(electron_density, "n");
    data_out.add_data_vector(hole_density,     "p");
    data_out.add_data_vector(old_electron_density, "old_n");
    data_out.add_data_vector(old_hole_density,     "old_p");
    

    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");
    data_out.build_patches();
    data_out.write_vtu_with_pvtu_record(output_directory, "solution", cycle, mpi_communicator, 2, 1);

    //pcout << " End of output_results"<< std::endl<<std::endl;

  }
  //###################################################################################################################################################
  
  //HELPER FUNCTION IMPLEMENTATION

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


}
