

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
 *          Timo Heister, Clemson University, 2018 
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


#include <deal.II/meshworker/mesh_loop.h> 


#include <iostream> 
#include <fstream> 

using namespace dealii; 

namespace Step16 
{ 


  template <int dim> 
  struct ScratchData 
  { 
    ScratchData(const Mapping<dim> &      mapping, 
                const FiniteElement<dim> &fe, 
                const unsigned int        quadrature_degree, 
                const UpdateFlags         update_flags) 
      : fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags) 
    {} 

    ScratchData(const ScratchData<dim> &scratch_data) 
      : fe_values(scratch_data.fe_values.get_mapping(), 
                  scratch_data.fe_values.get_fe(), 
                  scratch_data.fe_values.get_quadrature(), 
                  scratch_data.fe_values.get_update_flags()) 
    {} 

    FEValues<dim> fe_values; 
  }; 

  struct CopyData 
  { 
    unsigned int                         level; 
    FullMatrix<double>                   cell_matrix; 
    Vector<double>                       cell_rhs; 
    std::vector<types::global_dof_index> local_dof_indices; 

    template <class Iterator> 
    void reinit(const Iterator &cell, unsigned int dofs_per_cell) 
    { 
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell); 
      cell_rhs.reinit(dofs_per_cell); 

      local_dof_indices.resize(dofs_per_cell); 
      cell->get_active_or_mg_dof_indices(local_dof_indices); 
      level = cell->level(); 
    } 
  }; 




  template <int dim> 
  class LaplaceProblem 
  { 
  public: 
    LaplaceProblem(const unsigned int degree); 
    void run(); 

  private: 
    template <class Iterator> 
    void cell_worker(const Iterator &  cell, 
                     ScratchData<dim> &scratch_data, 
                     CopyData &        copy_data); 

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




    MGLevelObject<SparsityPattern> mg_sparsity_patterns; 
    MGLevelObject<SparsityPattern> mg_interface_sparsity_patterns; 

    MGLevelObject<SparseMatrix<double>> mg_matrices; 
    MGLevelObject<SparseMatrix<double>> mg_interface_matrices; 
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

    std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
              << " (by level: "; 
    for (unsigned int level = 0; level < triangulation.n_levels(); ++level) 
      std::cout << dof_handler.n_dofs(level) 
                << (level == triangulation.n_levels() - 1 ? ")" : ", "); 
    std::cout << std::endl; 

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

    { 
      DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints); 
      sparsity_pattern.copy_from(dsp); 
    } 
    system_matrix.reinit(sparsity_pattern); 


    mg_constrained_dofs.clear(); 
    mg_constrained_dofs.initialize(dof_handler); 
    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, 
                                                       dirichlet_boundary_ids); 


    const unsigned int n_levels = triangulation.n_levels(); 

    mg_interface_matrices.resize(0, n_levels - 1); 
    mg_matrices.resize(0, n_levels - 1); 
    mg_sparsity_patterns.resize(0, n_levels - 1); 
    mg_interface_sparsity_patterns.resize(0, n_levels - 1); 



    for (unsigned int level = 0; level < n_levels; ++level) 
      { 
        { 
          DynamicSparsityPattern dsp(dof_handler.n_dofs(level), 
                                     dof_handler.n_dofs(level)); 
          MGTools::make_sparsity_pattern(dof_handler, dsp, level); 

          mg_sparsity_patterns[level].copy_from(dsp); 
          mg_matrices[level].reinit(mg_sparsity_patterns[level]); 
        } 
        { 
          DynamicSparsityPattern dsp(dof_handler.n_dofs(level), 
                                     dof_handler.n_dofs(level)); 
          MGTools::make_interface_sparsity_pattern(dof_handler, 
                                                   mg_constrained_dofs, 
                                                   dsp, 
                                                   level); 
          mg_interface_sparsity_patterns[level].copy_from(dsp); 
          mg_interface_matrices[level].reinit( 
            mg_interface_sparsity_patterns[level]); 
        } 
      } 
  } 



  template <int dim> 
  template <class Iterator> 
  void LaplaceProblem<dim>::cell_worker(const Iterator &  cell, 
                                        ScratchData<dim> &scratch_data, 
                                        CopyData &        copy_data) 
  { 
    FEValues<dim> &fe_values = scratch_data.fe_values; 
    fe_values.reinit(cell); 

    const unsigned int dofs_per_cell = fe_values.get_fe().n_dofs_per_cell(); 
    const unsigned int n_q_points    = fe_values.get_quadrature().size(); 

    copy_data.reinit(cell, dofs_per_cell); 

    const std::vector<double> &JxW = fe_values.get_JxW_values(); 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      { 
        const double coefficient = 
          (fe_values.get_quadrature_points()[q][0] < 0.0) ? 1.0 : 0.1; 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          { 
            for (unsigned int j = 0; j < dofs_per_cell; ++j) 
              { 
                copy_data.cell_matrix(i, j) += 
                  coefficient * 
                  (fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q)) * 
                  JxW[q]; 
              } 
            copy_data.cell_rhs(i) += 1.0 * fe_values.shape_value(i, q) * JxW[q]; 
          } 
      } 
  } 



  template <int dim> 
  void LaplaceProblem<dim>::assemble_system() 
  { 
    MappingQ1<dim> mapping; 

    auto cell_worker = 
      [&](const typename DoFHandler<dim>::active_cell_iterator &cell, 
          ScratchData<dim> &                                    scratch_data, 
          CopyData &                                            copy_data) { 
        this->cell_worker(cell, scratch_data, copy_data); 
      }; 

    auto copier = [&](const CopyData &cd) { 
      this->constraints.distribute_local_to_global(cd.cell_matrix, 
                                                   cd.cell_rhs, 
                                                   cd.local_dof_indices, 
                                                   system_matrix, 
                                                   system_rhs); 
    }; 

    const unsigned int n_gauss_points = degree + 1; 

    ScratchData<dim> scratch_data(mapping, 
                                  fe, 
                                  n_gauss_points, 
                                  update_values | update_gradients | 
                                    update_JxW_values | 
                                    update_quadrature_points); 

    MeshWorker::mesh_loop(dof_handler.begin_active(), 
                          dof_handler.end(), 
                          cell_worker, 
                          copier, 
                          scratch_data, 
                          CopyData(), 
                          MeshWorker::assemble_own_cells); 
  } 



  template <int dim> 
  void LaplaceProblem<dim>::assemble_multigrid() 
  { 
    MappingQ1<dim>     mapping; 
    const unsigned int n_levels = triangulation.n_levels(); 

    std::vector<AffineConstraints<double>> boundary_constraints(n_levels); 
    for (unsigned int level = 0; level < n_levels; ++level) 
      { 
        IndexSet dofset; 
        DoFTools::extract_locally_relevant_level_dofs(dof_handler, 
                                                      level, 
                                                      dofset); 
        boundary_constraints[level].reinit(dofset); 
        boundary_constraints[level].add_lines( 
          mg_constrained_dofs.get_refinement_edge_indices(level)); 
        boundary_constraints[level].add_lines( 
          mg_constrained_dofs.get_boundary_indices(level)); 
        boundary_constraints[level].close(); 
      } 

    auto cell_worker = 
      [&](const typename DoFHandler<dim>::level_cell_iterator &cell, 
          ScratchData<dim> &                                   scratch_data, 
          CopyData &                                           copy_data) { 
        this->cell_worker(cell, scratch_data, copy_data); 
      }; 

    auto copier = [&](const CopyData &cd) { 
      boundary_constraints[cd.level].distribute_local_to_global( 
        cd.cell_matrix, cd.local_dof_indices, mg_matrices[cd.level]); 

      const unsigned int dofs_per_cell = cd.local_dof_indices.size(); 


      for (unsigned int i = 0; i < dofs_per_cell; ++i) 
        for (unsigned int j = 0; j < dofs_per_cell; ++j) 
          if (mg_constrained_dofs.is_interface_matrix_entry( 
                cd.level, cd.local_dof_indices[i], cd.local_dof_indices[j])) 
            { 
              mg_interface_matrices[cd.level].add(cd.local_dof_indices[i], 
                                                  cd.local_dof_indices[j], 
                                                  cd.cell_matrix(i, j)); 
            } 
    }; 

    const unsigned int n_gauss_points = degree + 1; 

    ScratchData<dim> scratch_data(mapping, 
                                  fe, 
                                  n_gauss_points, 
                                  update_values | update_gradients | 
                                    update_JxW_values | 
                                    update_quadrature_points); 

    MeshWorker::mesh_loop(dof_handler.begin_mg(), 
                          dof_handler.end_mg(), 
                          cell_worker, 
                          copier, 
                          scratch_data, 
                          CopyData(), 
                          MeshWorker::assemble_own_cells); 
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
    mg::Matrix<Vector<double>> mg_interface_up(mg_interface_matrices); 
    mg::Matrix<Vector<double>> mg_interface_down(mg_interface_matrices); 


    Multigrid<Vector<double>> mg( 
      mg_matrix, coarse_grid_solver, mg_transfer, mg_smoother, mg_smoother); 
    mg.set_edge_matrices(mg_interface_down, mg_interface_up); 

    PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>> 
      preconditioner(dof_handler, mg, mg_transfer); 


    SolverControl            solver_control(1000, 1e-12); 
    SolverCG<Vector<double>> solver(solver_control); 

    solution = 0; 

    solver.solve(system_matrix, solution, system_rhs, preconditioner); 
    std::cout << "   Number of CG iterations: " << solver_control.last_step() 
              << "\n" 
              << std::endl; 
    constraints.distribute(solution); 
  } 



  template <int dim> 
  void LaplaceProblem<dim>::refine_grid() 
  { 
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      QGauss<dim - 1>(degree + 2), 
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
        std::cout << "Cycle " << cycle << std::endl; 

        if (cycle == 0) 
          { 
            GridGenerator::hyper_ball(triangulation); 
            triangulation.refine_global(2); 
          } 
        else 
          refine_grid(); 

        std::cout << "   Number of active cells:       " 
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



