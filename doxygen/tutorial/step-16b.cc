

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2003 - 2021 by the deal.II authors 
 * 
 * This file is part of the deal.II library. 
 * 
 * The deal.II library is free software; you can use it, redistribute 
 * it, and/or modify it under the terms of the GNU Lesser General 
 * Public License as published by the Free Software Foundation; either 
 * version 2.1 of the License, or (at your option) any later version. 
 * The full text of the license can be found in the file LICENSE.md at 
 * the top level directory of deal.II. 
 * 
 * --------------------------------------------------------------------- 
 * 
 * Authors: Guido Kanschat, University of Heidelberg, 2003 
 *          Baerbel Janssen, University of Heidelberg, 2010 
 *          Wolfgang Bangerth, Texas A&M University, 2010 
 */ 




#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 



#include <deal.II/multigrid/mg_constrained_dofs.h> 
#include <deal.II/multigrid/multigrid.h> 
#include <deal.II/multigrid/mg_transfer.h> 
#include <deal.II/multigrid/mg_tools.h> 
#include <deal.II/multigrid/mg_coarse.h> 
#include <deal.II/multigrid/mg_smoother.h> 
#include <deal.II/multigrid/mg_matrix.h> 


#include <deal.II/meshworker/dof_info.h> 
#include <deal.II/meshworker/integration_info.h> 
#include <deal.II/meshworker/simple.h> 
#include <deal.II/meshworker/output.h> 
#include <deal.II/meshworker/loop.h> 


#include <deal.II/integrators/laplace.h> 
#include <deal.II/integrators/l2.h> 


#include <iostream> 
#include <fstream> 

using namespace dealii; 

namespace Step16 
{ 


  template <int dim> 
  class LaplaceIntegrator : public MeshWorker::LocalIntegrator<dim> 
  { 
  public: 
    LaplaceIntegrator(); 
    virtual void cell(MeshWorker::DoFInfo<dim> &        dinfo, 
                      MeshWorker::IntegrationInfo<dim> &info) const override; 
  }; 

  template <int dim> 
  LaplaceIntegrator<dim>::LaplaceIntegrator() 
    : MeshWorker::LocalIntegrator<dim>(true, false, false) 
  {} 





  template <int dim> 
  void 
  LaplaceIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &        dinfo, 
                               MeshWorker::IntegrationInfo<dim> &info) const 
  { 
    AssertDimension(dinfo.n_matrices(), 1); 
    const double coefficient = (dinfo.cell->center()(0) > 0.) ? .1 : 1.; 

    LocalIntegrators::Laplace::cell_matrix(dinfo.matrix(0, false).matrix, 
                                           info.fe_values(0), 
                                           coefficient); 

    if (dinfo.n_vectors() > 0) 
      { 
        std::vector<double> rhs(info.fe_values(0).n_quadrature_points, 1.); 
        LocalIntegrators::L2::L2(dinfo.vector(0).block(0), 
                                 info.fe_values(0), 
                                 rhs); 
      } 
  } 


  template <int dim> 
  class LaplaceProblem 
  { 
  public: 
    LaplaceProblem(const unsigned int degree); 
    void run(); 

  private: 
    void setup_system(); 
    void assemble_system(); 
    void assemble_multigrid(); 
    void solve(); 
    void refine_grid(); 
    void output_results(const unsigned int cycle) const; 

    Triangulation<dim> triangulation; 
    FE_Q<dim>          fe; 
    DoFHandler<dim>    dof_handler; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    AffineConstraints<double> constraints; 

    Vector<double> solution; 
    Vector<double> system_rhs; 

    const unsigned int degree; 




    MGLevelObject<SparsityPattern>      mg_sparsity_patterns; 
    MGLevelObject<SparseMatrix<double>> mg_matrices; 
    MGLevelObject<SparseMatrix<double>> mg_interface_in; 
    MGLevelObject<SparseMatrix<double>> mg_interface_out; 
    MGConstrainedDoFs                   mg_constrained_dofs; 
  }; 




  template <int dim> 
  LaplaceProblem<dim>::LaplaceProblem(const unsigned int degree) 
    : triangulation(Triangulation<dim>::limit_level_difference_at_vertices) 
    , fe(degree) 
    , dof_handler(triangulation) 
    , degree(degree) 
  {} 



  template <int dim> 
  void LaplaceProblem<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 
    dof_handler.distribute_mg_dofs(); 

    deallog << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
            << " (by level: "; 
    for (unsigned int level = 0; level < triangulation.n_levels(); ++level) 
      deallog << dof_handler.n_dofs(level) 
              << (level == triangulation.n_levels() - 1 ? ")" : ", "); 
    deallog << std::endl; 

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 

    constraints.clear(); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 

    std::set<types::boundary_id> dirichlet_boundary_ids = {0}; 
    Functions::ZeroFunction<dim> homogeneous_dirichlet_bc; 
    const std::map<types::boundary_id, const Function<dim> *> 
      dirichlet_boundary_functions = { 
        {types::boundary_id(0), &homogeneous_dirichlet_bc}}; 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             dirichlet_boundary_functions, 
                                             constraints); 
    constraints.close(); 
    constraints.condense(dsp); 
    sparsity_pattern.copy_from(dsp); 
    system_matrix.reinit(sparsity_pattern); 


    mg_constrained_dofs.clear(); 
    mg_constrained_dofs.initialize(dof_handler); 
    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, 
                                                       dirichlet_boundary_ids); 


    const unsigned int n_levels = triangulation.n_levels(); 

    mg_interface_in.resize(0, n_levels - 1); 
    mg_interface_in.clear_elements(); 
    mg_interface_out.resize(0, n_levels - 1); 
    mg_interface_out.clear_elements(); 
    mg_matrices.resize(0, n_levels - 1); 
    mg_matrices.clear_elements(); 
    mg_sparsity_patterns.resize(0, n_levels - 1); 



    for (unsigned int level = 0; level < n_levels; ++level) 
      { 
        DynamicSparsityPattern dsp(dof_handler.n_dofs(level), 
                                   dof_handler.n_dofs(level)); 
        MGTools::make_sparsity_pattern(dof_handler, dsp, level); 

        mg_sparsity_patterns[level].copy_from(dsp); 

        mg_matrices[level].reinit(mg_sparsity_patterns[level]); 
        mg_interface_in[level].reinit(mg_sparsity_patterns[level]); 
        mg_interface_out[level].reinit(mg_sparsity_patterns[level]); 
      } 
  } 







  template <int dim> 
  void LaplaceProblem<dim>::assemble_system() 
  { 
    MappingQ1<dim>                      mapping; 
    MeshWorker::IntegrationInfoBox<dim> info_box; 
    UpdateFlags                         update_flags = 
      update_values | update_gradients | update_hessians; 
    info_box.add_update_flags_all(update_flags); 
    info_box.initialize(fe, mapping); 

    MeshWorker::DoFInfo<dim> dof_info(dof_handler); 

    MeshWorker::Assembler::SystemSimple<SparseMatrix<double>, Vector<double>> 
      assembler; 
    assembler.initialize(constraints); 
    assembler.initialize(system_matrix, system_rhs); 

    LaplaceIntegrator<dim> matrix_integrator; 
    MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(), 
                                           dof_handler.end(), 
                                           dof_info, 
                                           info_box, 
                                           matrix_integrator, 
                                           assembler); 

    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) 
      if (constraints.is_constrained(i)) 
        system_matrix.set(i, i, 1.); 
  } 


  template <int dim> 
  void LaplaceProblem<dim>::assemble_multigrid() 
  { 
    MappingQ1<dim>                      mapping; 
    MeshWorker::IntegrationInfoBox<dim> info_box; 
    UpdateFlags                         update_flags = 
      update_values | update_gradients | update_hessians; 
    info_box.add_update_flags_all(update_flags); 
    info_box.initialize(fe, mapping); 

    MeshWorker::DoFInfo<dim> dof_info(dof_handler); 

    MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double>> assembler; 
    assembler.initialize(mg_constrained_dofs); 
    assembler.initialize(mg_matrices); 
    assembler.initialize_interfaces(mg_interface_in, mg_interface_out); 

    LaplaceIntegrator<dim> matrix_integrator; 
    MeshWorker::integration_loop<dim, dim>(dof_handler.begin_mg(), 
                                           dof_handler.end_mg(), 
                                           dof_info, 
                                           info_box, 
                                           matrix_integrator, 
                                           assembler); 

    const unsigned int nlevels = triangulation.n_levels(); 
    for (unsigned int level = 0; level < nlevels; ++level) 
      { 
        for (unsigned int i = 0; i < dof_handler.n_dofs(level); ++i) 
          if (mg_constrained_dofs.is_boundary_index(level, i) || 
              mg_constrained_dofs.at_refinement_edge(level, i)) 
            mg_matrices[level].set(i, i, 1.); 
      } 
  } 





  template <int dim> 
  void LaplaceProblem<dim>::solve() 
  { 
    MGTransferPrebuilt<Vector<double>> mg_transfer(mg_constrained_dofs); 
    mg_transfer.build(dof_handler); 

    FullMatrix<double> coarse_matrix; 
    coarse_matrix.copy_from(mg_matrices[0]); 
    MGCoarseGridHouseholder<double, Vector<double>> coarse_grid_solver; 
    coarse_grid_solver.initialize(coarse_matrix); 






    using Smoother = PreconditionSOR<SparseMatrix<double>>; 
    mg::SmootherRelaxation<Smoother, Vector<double>> mg_smoother; 
    mg_smoother.initialize(mg_matrices); 
    mg_smoother.set_steps(2); 
    mg_smoother.set_symmetric(true); 


    mg::Matrix<Vector<double>> mg_matrix(mg_matrices); 
    mg::Matrix<Vector<double>> mg_interface_up(mg_interface_in); 
    mg::Matrix<Vector<double>> mg_interface_down(mg_interface_out); 


    Multigrid<Vector<double>> mg( 
      mg_matrix, coarse_grid_solver, mg_transfer, mg_smoother, mg_smoother); 
    mg.set_edge_matrices(mg_interface_down, mg_interface_up); 

    PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>> 
      preconditioner(dof_handler, mg, mg_transfer); 


    SolverControl            solver_control(1000, 1e-12); 
    SolverCG<Vector<double>> solver(solver_control); 

    solution = 0; 

    solver.solve(system_matrix, solution, system_rhs, preconditioner); 
    constraints.distribute(solution); 
  } 



  template <int dim> 
  void LaplaceProblem<dim>::refine_grid() 
  { 
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      QGauss<dim - 1>(fe.degree + 1), 
      std::map<types::boundary_id, const Function<dim> *>(), 
      solution, 
      estimated_error_per_cell); 
    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    estimated_error_per_cell, 
                                                    0.3, 
                                                    0.03); 
    triangulation.execute_coarsening_and_refinement(); 
  } 

  template <int dim> 
  void LaplaceProblem<dim>::output_results(const unsigned int cycle) const 
  { 
    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "solution"); 
    data_out.build_patches(); 

    std::ofstream output("solution-" + std::to_string(cycle) + ".vtk"); 
    data_out.write_vtk(output); 
  } 


  template <int dim> 
  void LaplaceProblem<dim>::run() 
  { 
    for (unsigned int cycle = 0; cycle < 8; ++cycle) 
      { 
        deallog << "Cycle " << cycle << std::endl; 

        if (cycle == 0) 
          { 
            GridGenerator::hyper_ball(triangulation); 
            triangulation.refine_global(1); 
          } 
        else 
          refine_grid(); 

        deallog << "   Number of active cells:       " 
                << triangulation.n_active_cells() << std::endl; 

        setup_system(); 

        assemble_system(); 
        assemble_multigrid(); 

        solve(); 
        output_results(cycle); 
      } 
  } 
} // namespace Step16 


int main() 
{ 
  try 
    { 
      using namespace Step16; 

      deallog.depth_console(2); 

      LaplaceProblem<2> laplace_problem(1); 
      laplace_problem.run(); 
    } 
  catch (std::exception &exc) 
    { 
      std::cerr << std::endl 
                << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      std::cerr << "Exception on processing: " << std::endl 
                << exc.what() << std::endl 
                << "Aborting!" << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 

      return 1; 
    } 
  catch (...) 
    { 
      std::cerr << std::endl 
                << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      std::cerr << "Unknown exception!" << std::endl 
                << "Aborting!" << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      return 1; 
    } 

  return 0; 
} 

