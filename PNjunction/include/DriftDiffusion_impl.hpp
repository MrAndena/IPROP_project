namespace {

  using namespace dealii;

  //------------------------------------------------------------------------------------------------------------------------------
  //CONSTRUCTOR
  template <int dim>
  DriftDiffusion<dim>::DriftDiffusion(parallel::distributed::Triangulation<dim> &tria)
    : mpi_communicator(MPI_COMM_WORLD)   // Initialize MPI communicator to handle parallel processes
    , triangulation(tria)                // Assign the passed triangulation object to the class variable
    , fe(1)                              // Set up a finite element object with polynomial degree 1
    , dof_handler(tria)                  // Initialize DoF handler with the triangulation
    , pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))   // Set up a conditional output stream (only rank 0 outputs to console)
    , mapping()                          // Initialize the default mapping object (identity mapping in this case)
  {}
  
  //-----------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  void DriftDiffusion<dim>::setup_system()
  {
    // Distribute DoFs (degrees of freedom) across all processors for the finite element (FE) system.
    dof_handler.distribute_dofs(fe);

    // INDEX SETS INITIALIZATION
    locally_owned_dofs = dof_handler.locally_owned_dofs();                           // Locally owned DoFs (non-ghosted).
    locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);    // Locally relevant DoFs (includes ghost DoFs).
    

    // PETSC VECTORS DECLARATIONS 
    current_solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);  // Ghosted vector for the current solution.
    newton_update.reinit(locally_owned_dofs, mpi_communicator);                            // Non-ghosted vector for Newton updates.
    poisson_system_rhs.reinit(locally_owned_dofs, mpi_communicator);                       // Non-ghosted vector for Poisson right-hand side.

    electron_density.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);  // Ghosted vector for electron density.
    hole_density.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);      // Ghosted vector for hole density.

    old_electron_density.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);  
    old_hole_density.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);      

    rhs_electron_density.reinit(locally_owned_dofs, mpi_communicator);  
    rhs_hole_density.reinit(locally_owned_dofs, mpi_communicator);      
    

    // zero_constraints for newton poisson problem
    zero_constraints.clear();
    zero_constraints.reinit(locally_relevant_dofs);   // Define empty zero constraints for Poisson problem.
    VectorTools::interpolate_boundary_values(dof_handler, 1, Functions::ZeroFunction<dim>(), elec_constraints); 
    VectorTools::interpolate_boundary_values(dof_handler, 2, Functions::ZeroFunction<dim>(), elec_constraints); 
    zero_constraints.close(); 
    
    // Boundary constraints for electrons
    elec_constraints.clear();
    elec_constraints.reinit(locally_relevant_dofs);
    VectorTools::interpolate_boundary_values(dof_handler, 1, Functions::ConstantFunction<dim>(N1), elec_constraints); 
    VectorTools::interpolate_boundary_values(dof_handler, 2, Functions::ConstantFunction<dim>(N2), elec_constraints); 
    elec_constraints.close();

    // Boundary constraints for holes
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


    // Initialize the Poisson system matrix
    poisson_system_matrix.clear(); 
    poisson_system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

    // Initialize Laplace, mass, and density matrices
    laplace_matrix.clear(); 
    laplace_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

    mass_matrix.clear();  
    mass_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

    density_matrix.clear(); 
    density_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp,  mpi_communicator);

  // Sparsity pattern for electron constraints
    DynamicSparsityPattern elec_dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, elec_dsp, elec_constraints, false);
    SparsityTools::distribute_sparsity_pattern(elec_dsp,
                                            dof_handler.locally_owned_dofs(),
                                            mpi_communicator,
                                            locally_relevant_dofs);

    // Sparsity pattern for hole constraints
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
  template <int dim>
  void DriftDiffusion<dim>:: initialization(){

  // Create temporary vectors for potential, hole density, and electron density
  PETScWrappers::MPI::Vector temp_pot(locally_owned_dofs, mpi_communicator);
  PETScWrappers::MPI::Vector temp_hole(locally_owned_dofs, mpi_communicator);
  PETScWrappers::MPI::Vector temp_elec(locally_owned_dofs, mpi_communicator);

  // Initialize the temporary vectors with current solution values
  temp_pot = current_solution;
  temp_hole = old_hole_density;
  temp_elec = old_electron_density;

  // Interpolate initial or boundary values into each vector:
  // - Potential values
  // - Initial values for hole density
  // - Initial values for electron density
	VectorTools::interpolate(mapping, dof_handler, PotentialValues<dim>(), temp_pot);
	VectorTools::interpolate(mapping, dof_handler, HoleInitialValues<dim>(), temp_hole);
	VectorTools::interpolate(mapping, dof_handler, ElectronInitialValues<dim>(), temp_elec);
  
  // Compress the vectors after insertion to prepare them for parallel use
  temp_pot.compress(VectorOperation::insert);
  temp_hole.compress(VectorOperation::insert);
  temp_elec.compress(VectorOperation::insert);

  // Assign the initialized values to the main solution vectors
  current_solution = temp_pot;
  old_hole_density = temp_hole;
  old_electron_density = temp_elec;
  }
  
  //-----------------------------------------------------------------------------------------------------------------------------
  
  template <int dim>
  void DriftDiffusion<dim>::assemble_laplace_matrix()
  {
    // Define quadrature formula for numerical integration
    const QTrapezoid<dim> quadrature_formula;

    // Initialize the global Laplace matrix to zero
    laplace_matrix = 0;

    // Set up FEValues to compute gradients, values, points, and Jacobians for the finite element
    FEValues<dim> fe_values(fe,
			    quadrature_formula,
			    update_values | update_gradients |
			    update_quadrature_points | update_JxW_values);

    // Get number of degrees of freedom (DoFs) per cell and quadrature points
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    // Temporary variables for the local cell matrix and DoF indices
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  
    // Loop over all cells in the domain
    unsigned int q_point = 0, idof = 0, jdof = 0;
    for (const auto &cell : dof_handler.active_cell_iterators()){
      if (cell->is_locally_owned()){
	      cell_matrix = 0.;
        fe_values.reinit(cell);   // Reinitialize fe_values for the current cell
	
        // Loop over all quadrature points
        for (q_point = 0; q_point < n_q_points; ++q_point) {
          // Loop over all pairs of test functions to assemble local matrix
          for (idof = 0; idof < dofs_per_cell; ++idof) {
            for (jdof = 0; jdof < dofs_per_cell; ++jdof)
              // Compute entry for Laplace matrix using shape function gradients and JxW
              cell_matrix(idof, jdof) += fe_values.shape_grad(idof, q_point) *
                fe_values.shape_grad(jdof, q_point) * fe_values.JxW(q_point);
          }
	      }

        // Get global indices of the DoFs for the current cell
        cell->get_dof_indices(local_dof_indices);

        // Distribute local cell matrix values to the global matrix
        zero_constraints.distribute_local_to_global(cell_matrix,       
                                               local_dof_indices,
                                               laplace_matrix);
      }
    }

     // Compress the global matrix to finalize assembly for parallel processing
    laplace_matrix.compress(VectorOperation::add);
  }

  //------------------------------------------------------------------------------------------------------------------------------

  template <int dim>
  void DriftDiffusion<dim>::assemble_mass_matrix()
  {
    const QTrapezoid<dim> quadrature_formula;

    mass_matrix = 0;

    // Set up FEValues to compute shape function values, gradients, points, and Jacobian weights for the finite element
    FEValues<dim> fe_values(fe,
			    quadrature_formula,
			    update_values | update_gradients |
			    update_quadrature_points | update_JxW_values);

    // Number of degrees of freedom (DoFs) per cell and quadrature points
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    // Local cell matrix and vector for DoF indices
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    unsigned int q_point = 0, idof = 0, jdof = 0;

    // Loop over all cells in the domain
    for (const auto &cell : dof_handler.active_cell_iterators()){
      if (cell->is_locally_owned())
      {
        cell_matrix = 0.;

        fe_values.reinit(cell);

        for (q_point = 0; q_point < n_q_points; ++q_point) {

          // Assemble local matrix by looping over pairs of test functions
          for (idof = 0; idof < dofs_per_cell; ++idof) {
            for (jdof = 0; jdof < dofs_per_cell; ++jdof)

              // Calculate entry for the mass matrix using shape function values and JxW
              cell_matrix(idof, jdof) += fe_values.shape_value(idof, q_point) *
                fe_values.shape_value(jdof, q_point) * fe_values.JxW(q_point);
          }

        }

        cell->get_dof_indices(local_dof_indices);

        // Distribute local cell matrix values to the global mass matrix
        zero_constraints.distribute_local_to_global(cell_matrix,   
                                              local_dof_indices,
                                              mass_matrix );
      }
    }

    // Finalize assembly for the parallel global matrix
    mass_matrix.compress(VectorOperation::add);

  }
  
  //----------------------------------------------------------------------------------------------------------------------------  
  template <int dim>
  void DriftDiffusion<dim>::assemble_nonlinear_poisson()
  {

    // Initialize Poisson system matrix and density matrix for the Newton method
    poisson_system_matrix = 0;
    density_matrix = 0;
  
    // Temporary non-ghosted vectors for storing electron and hole densities
    PETScWrappers::MPI::Vector temp_1(locally_owned_dofs, mpi_communicator); 
    PETScWrappers::MPI::Vector temp_2(locally_owned_dofs, mpi_communicator); 

    // Copy old electron and hole densities into temporary vectors
    temp_1 = old_electron_density; //temp_1 store electron_density "n"
    temp_2 = old_hole_density;     //temp_2 store hole_density "p"

    double new_value = 0;

    // Compute the term (n + p) * MASS_MATRIX in a lumped version, storing in density_matrix
    for (auto iter = locally_owned_dofs.begin(); iter != locally_owned_dofs.end(); ++iter){ 

      // Compute and store the diagonal entries as (n + p) * mass_matrix(*iter, *iter)
      new_value = mass_matrix(*iter, *iter) * (temp_1[*iter] + temp_2[*iter]);
      density_matrix.set(*iter,*iter,new_value);
    }

    density_matrix.compress(VectorOperation::insert);

    // Build Poisson system matrix for Newton's method:
    // SYS_MAT = SYS_MAT + eps * Laplace matrix
    poisson_system_matrix.add(eps_r * eps_0, laplace_matrix);

    // Add the density matrix term: SYS_MAT = SYS_MAT + (q0 / V_TH) * density_matrix
    poisson_system_matrix.add(q0 / V_TH, density_matrix);    
  
    // Initialize the system right-hand side (RHS)
    poisson_system_rhs = 0;

    // Compute (n - p - N) and store it in temp_1
    temp_1.add(-1., temp_2);   

    // Interpolate doping profile values (N) into temp_2
    VectorTools::interpolate(mapping, dof_handler, DopingValues<dim>(), temp_2);  

    // Compute (n - p - N) and store in temp_1
    temp_1.add(-1., temp_2);     

    // Calculate mass_matrix * (n - p - N) and store result in temp_2
    mass_matrix.vmult(temp_2,temp_1);

    // Compute RHS term: -q0 * MASS * (n - p - N)
    poisson_system_rhs.add(-q0, temp_2);       

    // Incorporate potential term in the RHS: -eps * Laplace_matrix * phi
    temp_1 = current_solution;                           // temp_1 = phi (current potential solution)
    laplace_matrix.vmult(temp_2, temp_1);                // temp_2 = Laplace matrix * phi
    poisson_system_rhs.add(- eps_r * eps_0, temp_2);     // Update RHS with -eps * Laplace * phi

  }
  //----------------------------------------------------------------------------------------------------------------------------

  template <int dim>
  void DriftDiffusion<dim>::solve_poisson()
  {

    // Step 1: Apply Zero Boundary Conditions
    // This sets zero boundary conditions on both emitter and collector for the Poisson system.
    std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

    // Set zero potential boundary conditions on emitter (boundary id = 1)
    VectorTools::interpolate_boundary_values(mapping, dof_handler,1, Functions::ZeroFunction<dim>(), emitter_boundary_values);
    MatrixTools::apply_boundary_values(emitter_boundary_values, poisson_system_matrix, newton_update, poisson_system_rhs);

    // Set zero potential boundary conditions on collector (boundary id = 2)
    VectorTools::interpolate_boundary_values(mapping, dof_handler,2, Functions::ZeroFunction<dim>(), collector_boundary_values);
    MatrixTools::apply_boundary_values(collector_boundary_values, poisson_system_matrix, newton_update, poisson_system_rhs);
    

    // Step 2: Solve the Poisson System
    // Set up and use the MUMPS direct solver to solve the Poisson system for newton_update
    SolverControl sc_p(dof_handler.n_dofs(), 1e-10);      // Solver control with tolerance 1e-10
    PETScWrappers::SparseDirectMUMPS solverMUMPS(sc_p); 

    // Solve the system: poisson_system_matrix * newton_update = poisson_system_rhs
    solverMUMPS.solve(poisson_system_matrix, newton_update, poisson_system_rhs);

    // Step 3: Clamp the Solution Update
    // Ensure that the values in newton_update remain within [-V_TH, V_TH] to avoid divergence
    for (auto iter = locally_owned_dofs.begin(); iter != locally_owned_dofs.end(); ++iter){ 
  
      if (newton_update[*iter] < -V_TH) { newton_update[*iter] = -V_TH; }
      else if (newton_update[*iter] > V_TH) { newton_update[*iter] = V_TH; }
      
    }

    newton_update.compress(VectorOperation::insert); 
    
    // Step 4: Update the Current Solution
    // Add the computed Newton update to the current potential solution
    PETScWrappers::MPI::Vector temp;
    temp.reinit(locally_owned_dofs, mpi_communicator);

    temp = current_solution;
    temp.add(1.0, newton_update);   // Add the Newton update to the existing solution

    current_solution = temp;
  
  }

  //-----------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  void DriftDiffusion<dim>:: compute_densities(){   
  
    // This function calculates the electron and hole densities based on the current potential solution.
    // It is typically called within each Newton cycle iteration to update density values.
    
    // Step 1: Initialize Temporary Vectors for Electron and Hole Densities
    PETScWrappers::MPI::Vector old_temp_elec(locally_owned_dofs, mpi_communicator);
    PETScWrappers::MPI::Vector old_temp_hole(locally_owned_dofs, mpi_communicator);

    old_temp_elec  = current_solution;
    old_temp_hole  = current_solution;

    // Step 2: Calculate Electron and Hole Densities
    // Loop over all locally owned degrees of freedom to update densities based on the potential
    for (unsigned int i = old_temp_elec.local_range().first; i < old_temp_elec.local_range().second; ++i){ 

      // Calculate electron density using n = ni * exp(ϕ / V_TH)
      old_temp_elec[i] =ni*std::exp(old_temp_elec[i]/V_TH);

      // Calculate hole density using p = ni * exp(-ϕ / V_TH)
      old_temp_hole[i] =ni*std::exp(-old_temp_hole[i]/V_TH);

    }

    // Step 3: Finalize and Update Densities
    old_temp_elec.compress(VectorOperation::insert);
    old_temp_hole.compress(VectorOperation::insert);

    // Update class members to store the computed electron and hole densities
    old_electron_density = old_temp_elec;
    old_hole_density = old_temp_hole;

  }

//--------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void DriftDiffusion<dim>::cycle_newton_poisson(const unsigned int max_iter_newton, 
                                               const double toll_newton){

  unsigned int counter = 0; // To keep track of newton iteration

  pcout << "   - START NEWTON METHOD - "<< std::endl;

  // Initialize increment norm to a very large value to ensure the loop begins
  double increment_norm = std::numeric_limits<double>::max(); // the increment norm is + inf

  pcout << "   Initial Newton Increment Norm dphi: " << increment_norm << std::endl<<std::endl;
  
  // Start Newton iteration loop
  while(counter < max_iter_newton && increment_norm > toll_newton){

    pcout << "   NEWTON ITERATION NUMBER: "<< counter +1<<std::endl;
    pcout << "   Assemble System Poisson Matrix"<< std::endl;

    // The Mass and Laplace matrices are already built in the run() method
    // Step 1: Assemble the nonlinear Poisson matrix for the current iteration
    assemble_nonlinear_poisson();     

    // Step 2: Solve the Poisson equation to get the potential increment (dphi)
    // Boundary conditions and clamping are applied within the solve_poisson() method
    solve_poisson(); 

    // Step 3: Update electron and hole densities based on the new potential solution
    compute_densities();

    // Calculate the L2 norm of the Newton update (dphi) to monitor convergence
    increment_norm = newton_update.l2_norm();
    pcout << "   Update Increment: "<<increment_norm<<std::endl<<std::endl;

    counter ++;

    // Check if the maximum number of iterations has been reached
    if(counter == max_iter_newton){
      pcout<< "   MAX NUMBER OF NEWTON ITERATIONS REACHED!"<<std::endl;
    }

  }

}

//---------------------------------------------------------------------------------------------------------------------------
template <int dim>
void DriftDiffusion<dim>::assemble_drift_diffusion_matrix()
{

  // Initialization of right-hand side vectors and matrices
	rhs_electron_density = 0;
	rhs_hole_density = 0;
  hole_matrix = 0;
	electron_matrix = 0;

  
  const unsigned int vertices_per_cell = 4;   // Assuming a quadrilateral element with 4 vertices

  std::vector<types::global_dof_index> local_dof_indices(vertices_per_cell);    // Store local dof indices

  const unsigned int t_size = 3;    // Size for the local matrices and vectors

  Vector<double> cell_rhs(t_size);  // Right-hand side vector for the current cell
  FullMatrix<double> A(t_size,t_size), B(t_size,t_size), neg_A(t_size,t_size), neg_B(t_size,t_size);    // Local matrices

  std::vector<types::global_dof_index> A_local_dof_indices(t_size);    // Local dof indices for matrix A
  std::vector<types::global_dof_index> B_local_dof_indices(t_size);    // Local dof indices for matrix B

  // Loop over all active cells in the mesh
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned()){   // Ensure that only locally owned cells are processed

        A = 0;
        B = 0;
        neg_A = 0;
        neg_B = 0;
        cell_rhs = 0;
        
        // Get the local dof indices for the current cell
        cell->get_dof_indices(local_dof_indices); 

        // // Get the coordinates of the vertices, Lexicographic ordering
        const Point<dim> v1 = cell->vertex(2); // top left
        const Point<dim> v2 = cell->vertex(3); // top right
        const Point<dim> v3 = cell->vertex(0); // bottom left
        const Point<dim> v4 = cell->vertex(1); // bottom right
        
        // Calculate potentials at the vertices (normalized by V_TH)
        const double u1 = -current_solution[local_dof_indices[2]]/V_TH;
        const double u2 = -current_solution[local_dof_indices[3]]/V_TH;
        const double u3 = -current_solution[local_dof_indices[0]]/V_TH;
        const double u4 = -current_solution[local_dof_indices[1]]/V_TH;

        // Calculate lengths of the sides
        const double l_alpha = side_length(v1,v4);
        const double l_beta = side_length(v2,v3);
        
        // Calculate differences in potentials
        const double alpha21 = (u1 - u2);
        const double alpha42 = (u2 - u4);
        const double alpha34 = (u4 - u3);
        const double alpha13 = (u3 - u1);
        
        // Negate the alpha values for matrix B
        const double neg_alpha21 =  - (u1 - u2);
        const double neg_alpha42 =  - (u2 - u4);
        const double neg_alpha34 = - (u4 - u3);
        const double neg_alpha13 = - (u3 - u1);

        // Determine how to split the cell based on the lengths
        if (l_alpha >= l_beta) {        // l_alpha is the longest diagonal: split by beta

              const double alpha23 =  (u3 - u2);
              const double neg_alpha23 = - (u3 - u2);
              
              // Triangle A:
              A= compute_triangle_matrix(v2,v1,v3, alpha21, alpha13, -alpha23, Dp);
              neg_A= compute_triangle_matrix(v2,v1,v3, neg_alpha21, neg_alpha13, -neg_alpha23, Dn);
              
              // Triangle B:
              B = compute_triangle_matrix(v3,v4,v2, alpha34, alpha42, alpha23, Dp);
              neg_B = compute_triangle_matrix(v3,v4,v2, neg_alpha34, neg_alpha42, neg_alpha23, Dn);
              
              // Assign local dof indices for the matrices
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
        
          // Distribute local matrices and rhs contributions to global matrices and vectors
          hole_constraints.distribute_local_to_global(A, cell_rhs,  A_local_dof_indices, hole_matrix, rhs_hole_density);
          hole_constraints.distribute_local_to_global(B, cell_rhs,  B_local_dof_indices, hole_matrix, rhs_hole_density);

          elec_constraints.distribute_local_to_global(neg_A, cell_rhs,  A_local_dof_indices, electron_matrix, rhs_electron_density);
          elec_constraints.distribute_local_to_global(neg_B, cell_rhs,  B_local_dof_indices, electron_matrix, rhs_electron_density);

        }
		  }
    
    // Compress the matrices and vectors to finalize assembly
    hole_matrix.compress(VectorOperation::add);
    electron_matrix.compress(VectorOperation::add);
    
    rhs_hole_density.compress(VectorOperation::add);
    rhs_electron_density.compress(VectorOperation::add);
   
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void DriftDiffusion<dim>::apply_drift_diffusion_boundary_conditions() 
{
  
  // Initialize temporary vectors for electron and hole densities
  PETScWrappers::MPI::Vector temp_elec(locally_owned_dofs, mpi_communicator);
  PETScWrappers::MPI::Vector temp_hole(locally_owned_dofs, mpi_communicator);     

  // Copy current densities to temporary vectors
  temp_elec = electron_density;
  temp_hole = hole_density;
  
  // Maps to store boundary values for emitter and collector boundaries
  std::map<types::global_dof_index, double> emitter_boundary_values, collector_boundary_values;

  // Apply boundary conditions for electron density
  VectorTools::interpolate_boundary_values(mapping, dof_handler,1, Functions::ConstantFunction<dim>(N1), emitter_boundary_values);
  MatrixTools::apply_boundary_values(emitter_boundary_values, electron_matrix, temp_elec, rhs_electron_density);

  // Apply boundary conditions for hole density
  VectorTools::interpolate_boundary_values(mapping, dof_handler,2, Functions::ConstantFunction<dim>(N2), collector_boundary_values);
  MatrixTools::apply_boundary_values(collector_boundary_values, electron_matrix, temp_elec, rhs_electron_density);

  VectorTools::interpolate_boundary_values(mapping, dof_handler,1, Functions::ConstantFunction<dim>(P1), emitter_boundary_values);
  MatrixTools::apply_boundary_values(emitter_boundary_values, hole_matrix, temp_hole, rhs_hole_density);

  VectorTools::interpolate_boundary_values(mapping, dof_handler,2, Functions::ConstantFunction<dim>(P2), collector_boundary_values);
  MatrixTools::apply_boundary_values(collector_boundary_values, hole_matrix, temp_hole, rhs_hole_density);
  
  // Update the electron and hole densities with the values from the temporary vectors
  electron_density = temp_elec;
  hole_density = temp_hole;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------
template <int dim>
void DriftDiffusion<dim>::solve_drift_diffusion()
{ 

  PETScWrappers::MPI::Vector temp_elec(locally_owned_dofs, mpi_communicator);
  PETScWrappers::MPI::Vector temp_hole(locally_owned_dofs, mpi_communicator);
  
  // Set up solver controls for holes and electrons
  SolverControl sc_hole(dof_handler.n_dofs(), 1e-10);
  SolverControl sc_elec(dof_handler.n_dofs(), 1e-10); 
  
  // Initialize solvers for holes and electrons using MUMPS
  PETScWrappers::SparseDirectMUMPS solverMUMPS_hole(sc_hole); 
  PETScWrappers::SparseDirectMUMPS solverMUMPS_elec(sc_elec); 
  
  // Solve the hole density system
  solverMUMPS_hole.solve(hole_matrix, temp_hole, rhs_hole_density);
  // Solve the electron density system
  solverMUMPS_elec.solve(electron_matrix, temp_elec, rhs_electron_density);

  hole_constraints.distribute(temp_hole);
  elec_constraints.distribute(temp_elec);

  temp_elec.compress(VectorOperation::insert);
  temp_hole.compress(VectorOperation::insert);

  // Update the main density vectors with the new solutions
  hole_density = temp_hole;
  electron_density = temp_elec;
}

//----------------------------------------------------------------------------------------------------------------------------------------
  template <int dim>
  void DriftDiffusion<dim>::run(const unsigned int max_iter, const double dd_toll) 
  {
  
    // Initialize cycle counter and error tolerances
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
    
    // Main iterative solving loop
    while ( (hole_err > hole_tol || electron_err > electron_tol) && cycle_drift_diffusion < max_iter)
    {

      ++cycle_drift_diffusion;

      pcout<< "   CYCLE NUMBER " <<cycle_drift_diffusion<<"  OF DRIFT DIFFUSION"<<std::endl;

      pcout<< "   NEWTON POISSON" <<std::endl;
      cycle_newton_poisson(max_newton_iterations, tol_newton);


      pcout<< "   Assemble drift diffusion matrix" <<std::endl;
      assemble_drift_diffusion_matrix();

      pcout<< "   Solve drift diffusion"<<std::endl;
      solve_drift_diffusion();

      // Update error for convergence
      temp = hole_density;

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

      // Output results for current cycle
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
