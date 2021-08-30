

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2021 by the deal.II authors 
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
 * Author: Justin O'Connor, Colorado State University, 2021. 
 */ 


#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/tensor.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/base/signaling_nan.h> 

#include <deal.II/lac/block_vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/block_sparse_matrix.h> 
#include <deal.II/lac/linear_operator.h> 
#include <deal.II/lac/packaged_operation.h> 
#include <deal.II/lac/sparse_direct.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/fe_dgq.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_q.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 

#include <iostream> 
#include <fstream> 
#include <algorithm> 






namespace SAND 
{ 
  using namespace dealii; 


  namespace SolutionComponents 
  { 
    template <int dim> 
    constexpr unsigned int density = 0; 
    template <int dim> 
    constexpr unsigned int displacement = 1; 
    template <int dim> 
    constexpr unsigned int unfiltered_density = 1 + dim; 
    template <int dim> 
    constexpr unsigned int displacement_multiplier = 2 + dim; 
    template <int dim> 
    constexpr unsigned int unfiltered_density_multiplier = 2 + 2 * dim; 
    template <int dim> 
    constexpr unsigned int density_lower_slack = 3 + 2 * dim; 
    template <int dim> 
    constexpr unsigned int density_lower_slack_multiplier = 4 + 2 * dim; 
    template <int dim> 
    constexpr unsigned int density_upper_slack = 5 + 2 * dim; 
    template <int dim> 
    constexpr unsigned int density_upper_slack_multiplier = 6 + 2 * dim; 
  } // namespace SolutionComponents 


  namespace SolutionBlocks 
  { 
    constexpr unsigned int density                        = 0; 
    constexpr unsigned int displacement                   = 1; 
    constexpr unsigned int unfiltered_density             = 2; 
    constexpr unsigned int displacement_multiplier        = 3; 
    constexpr unsigned int unfiltered_density_multiplier  = 4; 
    constexpr unsigned int density_lower_slack            = 5; 
    constexpr unsigned int density_lower_slack_multiplier = 6; 
    constexpr unsigned int density_upper_slack            = 7; 
    constexpr unsigned int density_upper_slack_multiplier = 8; 
  } // namespace SolutionBlocks 

  namespace BoundaryIds 
  { 
    constexpr types::boundary_id down_force = 101; 
    constexpr types::boundary_id no_force   = 102; 
  } // namespace BoundaryIds 

  namespace ValueExtractors 
  { 
    template <int dim> 
    const FEValuesExtractors::Scalar 
      densities(SolutionComponents::density<dim>); 
    template <int dim> 
    const FEValuesExtractors::Vector 
      displacements(SolutionComponents::displacement<dim>); 
    template <int dim> 
    const FEValuesExtractors::Scalar 
      unfiltered_densities(SolutionComponents::unfiltered_density<dim>); 
    template <int dim> 
    const FEValuesExtractors::Vector displacement_multipliers( 
      SolutionComponents::displacement_multiplier<dim>); 
    template <int dim> 
    const FEValuesExtractors::Scalar unfiltered_density_multipliers( 
      SolutionComponents::unfiltered_density_multiplier<dim>); 
    template <int dim> 
    const FEValuesExtractors::Scalar 
      density_lower_slacks(SolutionComponents::density_lower_slack<dim>); 
    template <int dim> 
    const FEValuesExtractors::Scalar density_lower_slack_multipliers( 
      SolutionComponents::density_lower_slack_multiplier<dim>); 
    template <int dim> 
    const FEValuesExtractors::Scalar 
      density_upper_slacks(SolutionComponents::density_upper_slack<dim>); 
    template <int dim> 
    const FEValuesExtractors::Scalar density_upper_slack_multipliers( 
      SolutionComponents::density_upper_slack_multiplier<dim>); 
  } // namespace ValueExtractors 



  template <int dim> 
  class SANDTopOpt 
  { 
  public: 
    SANDTopOpt(); 

    void run(); 

  private: 
    void create_triangulation(); 

    void setup_boundary_values(); 

    void setup_block_system(); 

    void setup_filter_matrix(); 

    void assemble_system(); 

    BlockVector<double> solve(); 

    std::pair<double, double> 
    calculate_max_step_size(const BlockVector<double> &state, 
                            const BlockVector<double> &step) const; 

    BlockVector<double> 
    calculate_test_rhs(const BlockVector<double> &test_solution) const; 

    double calculate_exact_merit(const BlockVector<double> &test_solution); 

    BlockVector<double> find_max_step(); 

    BlockVector<double> compute_scaled_step(const BlockVector<double> &state, 
                                            const BlockVector<double> &step, 
                                            const double descent_requirement); 

    bool check_convergence(const BlockVector<double> &state); 

    void output_results(const unsigned int j) const; 

    void write_as_stl(); 

    std::set<typename Triangulation<dim>::cell_iterator> 
    find_relevant_neighbors( 
      typename Triangulation<dim>::cell_iterator cell) const; 


    Triangulation<dim>        triangulation; 
    FESystem<dim>             fe; 
    DoFHandler<dim>           dof_handler; 
    AffineConstraints<double> constraints; 

    std::map<types::global_dof_index, double> boundary_values; 

    BlockSparsityPattern      sparsity_pattern; 
    BlockSparseMatrix<double> system_matrix; 

    SparsityPattern      filter_sparsity_pattern; 
    SparseMatrix<double> filter_matrix; 

    BlockVector<double> system_rhs; 
    BlockVector<double> nonlinear_solution; 

    const double density_ratio; 
    const double density_penalty_exponent; 
    const double filter_r; 
    double       penalty_multiplier; 
    double       barrier_size; 

    TimerOutput timer; 
  }; 



  template <int dim> 
  SANDTopOpt<dim>::SANDTopOpt() 
    : fe(FE_DGQ<dim>(0), 
         1, 
         (FESystem<dim>(FE_Q<dim>(1) ^ dim)), 
         1, 
         FE_DGQ<dim>(0), 
         1, 
         (FESystem<dim>(FE_Q<dim>(1) ^ dim)), 
         1, 
         FE_DGQ<dim>(0), 
         5) 
    , dof_handler(triangulation) 
    , density_ratio(.5) 
    , density_penalty_exponent(3) 
    , filter_r(.251) 
    , penalty_multiplier(1) 
    , timer(std::cout, TimerOutput::summary, TimerOutput::wall_times) 
  { 
    Assert(dim > 1, ExcNotImplemented()); 
  } 



  template <int dim> 
  void SANDTopOpt<dim>::create_triangulation() 
  { 
    Assert(dim == 2, ExcNotImplemented()); 
    GridGenerator::subdivided_hyper_rectangle(triangulation, 
                                              {6, 1}, 
                                              Point<dim>(0, 0), 
                                              Point<dim>(6, 1)); 

    triangulation.refine_global(3); 


    for (const auto &cell : triangulation.active_cell_iterators()) 
      { 
        for (const auto &face : cell->face_iterators()) 
          { 
            if (face->at_boundary()) 
              { 
                const auto center = face->center(); 
                if (std::fabs(center(1) - 1) < 1e-12) 
                  { 
                    if ((std::fabs(center(0) - 3) < .3)) 
                      face->set_boundary_id(BoundaryIds::down_force); 
                    else 
                      face->set_boundary_id(BoundaryIds::no_force); 
                  } 
                else 
                  face->set_boundary_id(BoundaryIds::no_force); 
              } 
          } 
      } 
  } 


  template <int dim> 
  void SANDTopOpt<dim>::setup_boundary_values() 
  { 
    boundary_values.clear(); 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        for (const auto &face : cell->face_iterators()) 
          { 
            if (face->at_boundary()) 
              { 
                const auto center = face->center(); 


                if (std::fabs(center(1) - 0) < 1e-12) 
                  { 
                    for (const auto vertex_number : cell->vertex_indices()) 
                      { 
                        const auto vert = cell->vertex(vertex_number); 

                        if (std::fabs(vert(0) - 0) < 1e-12 && 
                            std::fabs(vert(1) - 0) < 1e-12) 
                          { 
                            types::global_dof_index x_displacement = 
                              cell->vertex_dof_index(vertex_number, 0); 
                            types::global_dof_index y_displacement = 
                              cell->vertex_dof_index(vertex_number, 1); 
                            types::global_dof_index x_displacement_multiplier = 
                              cell->vertex_dof_index(vertex_number, 2); 
                            types::global_dof_index y_displacement_multiplier = 
                              cell->vertex_dof_index(vertex_number, 3); 

                            boundary_values[x_displacement]            = 0; 
                            boundary_values[y_displacement]            = 0; 
                            boundary_values[x_displacement_multiplier] = 0; 
                            boundary_values[y_displacement_multiplier] = 0; 
                          } 

                        else if (std::fabs(vert(0) - 6) < 1e-12 && 
                                 std::fabs(vert(1) - 0) < 1e-12) 
                          { 
                            types::global_dof_index y_displacement = 
                              cell->vertex_dof_index(vertex_number, 1); 
                            types::global_dof_index y_displacement_multiplier = 
                              cell->vertex_dof_index(vertex_number, 3); 

                            boundary_values[y_displacement]            = 0; 
                            boundary_values[y_displacement_multiplier] = 0; 
                          } 
                      } 
                  } 
              } 
          } 
      } 
  } 



  template <int dim> 
  void SANDTopOpt<dim>::setup_block_system() 
  { 
    std::vector<unsigned int> block_component(9, 2); 
    block_component[0] = 0; 
    block_component[1] = 1; 
    const std::vector<types::global_dof_index> dofs_per_block = 
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component); 

    const types::global_dof_index                     n_p = dofs_per_block[0]; 
    const types::global_dof_index                     n_u = dofs_per_block[1]; 
    const std::vector<BlockVector<double>::size_type> block_sizes = { 
      n_p, n_u, n_p, n_u, n_p, n_p, n_p, n_p, n_p}; 

    BlockDynamicSparsityPattern dsp(9, 9); 
    for (unsigned int k = 0; k < 9; ++k) 
      for (unsigned int j = 0; j < 9; ++j) 
        dsp.block(j, k).reinit(block_sizes[j], block_sizes[k]); 
    dsp.collect_sizes(); 




    Table<2, DoFTools::Coupling> coupling(2 * dim + 7, 2 * dim + 7); 
    { 
      using namespace SolutionComponents; 

      coupling[density<dim>][density<dim>] = DoFTools::always; 

      for (unsigned int i = 0; i < dim; ++i) 
        { 
          coupling[density<dim>][displacement<dim> + i] = DoFTools::always; 
          coupling[displacement<dim> + i][density<dim>] = DoFTools::always; 
        } 

      for (unsigned int i = 0; i < dim; ++i) 
        { 
          coupling[density<dim>][displacement_multiplier<dim> + i] = 
            DoFTools::always; 
          coupling[displacement_multiplier<dim> + i][density<dim>] = 
            DoFTools::always; 
        } 

      coupling[density<dim>][unfiltered_density_multiplier<dim>] = 
        DoFTools::always; 
      coupling[unfiltered_density_multiplier<dim>][density<dim>] = 
        DoFTools::always; 
      /*位移的联结  */ 
      for (unsigned int i = 0; i < dim; ++i) 
        { 
          for (unsigned int k = 0; k < dim; ++k) 
            { 
              coupling[displacement<dim> + i] 
                      [displacement_multiplier<dim> + k] = DoFTools::always; 
              coupling[displacement_multiplier<dim> + k] 
                      [displacement<dim> + i] = DoFTools::always; 
            } 
        } 
      /*松弛变量的耦合 */ 
      coupling[density_lower_slack<dim>][density_lower_slack<dim>] = 
        DoFTools::always; 
      coupling[density_lower_slack<dim>][density_upper_slack<dim>] = 
        DoFTools::always; 
      coupling[density_upper_slack<dim>][density_lower_slack<dim>] = 
        DoFTools::always; 

      coupling[density_lower_slack_multiplier<dim>] 
              [density_lower_slack_multiplier<dim>] = DoFTools::always; 
      coupling[density_lower_slack_multiplier<dim>] 
              [density_upper_slack_multiplier<dim>] = DoFTools::always; 
      coupling[density_upper_slack_multiplier<dim>] 
              [density_lower_slack_multiplier<dim>] = DoFTools::always; 
    } 


    const ComponentMask density_mask = 
      fe.component_mask(ValueExtractors::densities<dim>); 
    const IndexSet density_dofs = 
      DoFTools::extract_dofs(dof_handler, density_mask); 

    types::global_dof_index last_density_dof = 
      density_dofs.nth_index_in_set(density_dofs.n_elements() - 1); 
    constraints.clear(); 
    constraints.add_line(last_density_dof); 
    for (unsigned int i = 0; i < density_dofs.n_elements() - 1; ++i) 
      constraints.add_entry(last_density_dof, 
                            density_dofs.nth_index_in_set(i), 
                            -1); 
    constraints.set_inhomogeneity(last_density_dof, 0); 

    constraints.close(); 


    DoFTools::make_sparsity_pattern(dof_handler, coupling, dsp, constraints); 


    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        const unsigned int i = cell->active_cell_index(); 
        for (const auto &check_cell : find_relevant_neighbors(cell)) 
          { 
            const double distance = 
              cell->center().distance(check_cell->center()); 
            if (distance < filter_r) 
              { 
                dsp 
                  .block(SolutionBlocks::unfiltered_density, 
                         SolutionBlocks::unfiltered_density_multiplier) 
                  .add(i, check_cell->active_cell_index()); 
                dsp 
                  .block(SolutionBlocks::unfiltered_density_multiplier, 
                         SolutionBlocks::unfiltered_density) 
                  .add(i, check_cell->active_cell_index()); 
              } 
          } 
      } 


    sparsity_pattern.copy_from(dsp); 

    std::ofstream out("sparsity.plt"); 
    sparsity_pattern.print_gnuplot(out); 

    system_matrix.reinit(sparsity_pattern); 


    nonlinear_solution.reinit(block_sizes); 
    system_rhs.reinit(block_sizes); 

    { 
      using namespace SolutionBlocks; 
      nonlinear_solution.block(density).add(density_ratio); 
      nonlinear_solution.block(unfiltered_density).add(density_ratio); 
      nonlinear_solution.block(unfiltered_density_multiplier) 
        .add(density_ratio); 
      nonlinear_solution.block(density_lower_slack).add(density_ratio); 
      nonlinear_solution.block(density_lower_slack_multiplier).add(50); 
      nonlinear_solution.block(density_upper_slack).add(1 - density_ratio); 
      nonlinear_solution.block(density_upper_slack_multiplier).add(50); 
    } 
  } 



  template <int dim> 
  void SANDTopOpt<dim>::setup_filter_matrix() 
  { 


    filter_sparsity_pattern.copy_from( 
      sparsity_pattern.block(SolutionBlocks::unfiltered_density, 
                             SolutionBlocks::unfiltered_density_multiplier)); 
    filter_matrix.reinit(filter_sparsity_pattern); 


    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        const unsigned int i = cell->active_cell_index(); 
        for (const auto &check_cell : find_relevant_neighbors(cell)) 
          { 
            const double distance = 
              cell->center().distance(check_cell->center()); 
            if (distance < filter_r) 
              { 
                filter_matrix.add(i, 
                                  check_cell->active_cell_index(), 
                                  filter_r - distance); 


              } 
          } 
      } 


    for (unsigned int i = 0; i < filter_matrix.m(); ++i) 
      { 
        double denominator = 0; 
        for (SparseMatrix<double>::iterator iter = filter_matrix.begin(i); 
             iter != filter_matrix.end(i); 
             iter++) 
          denominator = denominator + iter->value(); 
        for (SparseMatrix<double>::iterator iter = filter_matrix.begin(i); 
             iter != filter_matrix.end(i); 
             iter++) 
          iter->value() = iter->value() / denominator; 
      } 
  } 


  template <int dim> 
  std::set<typename Triangulation<dim>::cell_iterator> 
  SANDTopOpt<dim>::find_relevant_neighbors( 
    typename Triangulation<dim>::cell_iterator cell) const 
  { 
    std::set<unsigned int>                               neighbor_ids; 
    std::set<typename Triangulation<dim>::cell_iterator> cells_to_check; 

    neighbor_ids.insert(cell->active_cell_index()); 
    cells_to_check.insert(cell); 

    bool new_neighbors_found; 
    do 
      { 
        new_neighbors_found = false; 
        for (const auto &check_cell : 
             std::vector<typename Triangulation<dim>::cell_iterator>( 
               cells_to_check.begin(), cells_to_check.end())) 
          { 
            for (const auto n : check_cell->face_indices()) 
              { 
                if (!(check_cell->face(n)->at_boundary())) 
                  { 
                    const auto & neighbor = check_cell->neighbor(n); 
                    const double distance = 
                      cell->center().distance(neighbor->center()); 
                    if ((distance < filter_r) && 
                        !(neighbor_ids.count(neighbor->active_cell_index()))) 
                      { 
                        cells_to_check.insert(neighbor); 
                        neighbor_ids.insert(neighbor->active_cell_index()); 
                        new_neighbors_found = true; 
                      } 
                  } 
              } 
          } 
      } 
    while (new_neighbors_found); 
    return cells_to_check; 
  } 



  template <int dim> 
  void SANDTopOpt<dim>::assemble_system() 
  { 
    TimerOutput::Scope t(timer, "assembly"); 

    system_matrix = 0; 
    system_rhs    = 0; 

    MappingQGeneric<dim> mapping(1); 
    QGauss<dim>          quadrature_formula(fe.degree + 1); 
    QGauss<dim - 1>      face_quadrature_formula(fe.degree + 1); 
    FEValues<dim>        fe_values(mapping, 
                            fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 
    FEFaceValues<dim>    fe_face_values(mapping, 
                                     fe, 
                                     face_quadrature_formula, 
                                     update_values | update_quadrature_points | 
                                       update_normal_vectors | 
                                       update_JxW_values); 

    const unsigned int dofs_per_cell = fe.dofs_per_cell; 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     dummy_cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    std::vector<double>                    lambda_values(n_q_points); 
    std::vector<double>                    mu_values(n_q_points); 
    const Functions::ConstantFunction<dim> lambda(1.); 
    const Functions::ConstantFunction<dim> mu(1.); 
    std::vector<Tensor<1, dim>>            rhs_values(n_q_points); 


    BlockVector<double> filtered_unfiltered_density_solution = 
      nonlinear_solution; 
    BlockVector<double> filter_adjoint_unfiltered_density_multiplier_solution = 
      nonlinear_solution; 

    filter_matrix.vmult(filtered_unfiltered_density_solution.block( 
                          SolutionBlocks::unfiltered_density), 
                        nonlinear_solution.block( 
                          SolutionBlocks::unfiltered_density)); 
    filter_matrix.Tvmult( 
      filter_adjoint_unfiltered_density_multiplier_solution.block( 
        SolutionBlocks::unfiltered_density_multiplier), 
      nonlinear_solution.block(SolutionBlocks::unfiltered_density_multiplier)); 

    std::vector<double>                  old_density_values(n_q_points); 
    std::vector<Tensor<1, dim>>          old_displacement_values(n_q_points); 
    std::vector<double>                  old_displacement_divs(n_q_points); 
    std::vector<SymmetricTensor<2, dim>> old_displacement_symmgrads(n_q_points); 
    std::vector<Tensor<1, dim>> old_displacement_multiplier_values(n_q_points); 
    std::vector<double>         old_displacement_multiplier_divs(n_q_points); 
    std::vector<SymmetricTensor<2, dim>> old_displacement_multiplier_symmgrads( 
      n_q_points); 
    std::vector<double> old_lower_slack_multiplier_values(n_q_points); 
    std::vector<double> old_upper_slack_multiplier_values(n_q_points); 
    std::vector<double> old_lower_slack_values(n_q_points); 
    std::vector<double> old_upper_slack_values(n_q_points); 
    std::vector<double> old_unfiltered_density_values(n_q_points); 
    std::vector<double> old_unfiltered_density_multiplier_values(n_q_points); 
    std::vector<double> filtered_unfiltered_density_values(n_q_points); 
    std::vector<double> filter_adjoint_unfiltered_density_multiplier_values( 
      n_q_points); 

    using namespace ValueExtractors; 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_matrix = 0; 

        cell->get_dof_indices(local_dof_indices); 

        fe_values.reinit(cell); 

        lambda.value_list(fe_values.get_quadrature_points(), lambda_values); 
        mu.value_list(fe_values.get_quadrature_points(), mu_values); 


        fe_values[densities<dim>].get_function_values(nonlinear_solution, 
                                                      old_density_values); 
        fe_values[displacements<dim>].get_function_values( 
          nonlinear_solution, old_displacement_values); 
        fe_values[displacements<dim>].get_function_divergences( 
          nonlinear_solution, old_displacement_divs); 
        fe_values[displacements<dim>].get_function_symmetric_gradients( 
          nonlinear_solution, old_displacement_symmgrads); 
        fe_values[displacement_multipliers<dim>].get_function_values( 
          nonlinear_solution, old_displacement_multiplier_values); 
        fe_values[displacement_multipliers<dim>].get_function_divergences( 
          nonlinear_solution, old_displacement_multiplier_divs); 
        fe_values[displacement_multipliers<dim>] 
          .get_function_symmetric_gradients( 
            nonlinear_solution, old_displacement_multiplier_symmgrads); 
        fe_values[density_lower_slacks<dim>].get_function_values( 
          nonlinear_solution, old_lower_slack_values); 
        fe_values[density_lower_slack_multipliers<dim>].get_function_values( 
          nonlinear_solution, old_lower_slack_multiplier_values); 
        fe_values[density_upper_slacks<dim>].get_function_values( 
          nonlinear_solution, old_upper_slack_values); 
        fe_values[density_upper_slack_multipliers<dim>].get_function_values( 
          nonlinear_solution, old_upper_slack_multiplier_values); 
        fe_values[unfiltered_densities<dim>].get_function_values( 
          nonlinear_solution, old_unfiltered_density_values); 
        fe_values[unfiltered_density_multipliers<dim>].get_function_values( 
          nonlinear_solution, old_unfiltered_density_multiplier_values); 
        fe_values[unfiltered_densities<dim>].get_function_values( 
          filtered_unfiltered_density_solution, 
          filtered_unfiltered_density_values); 
        fe_values[unfiltered_density_multipliers<dim>].get_function_values( 
          filter_adjoint_unfiltered_density_multiplier_solution, 
          filter_adjoint_unfiltered_density_multiplier_values); 

        for (const auto q_point : fe_values.quadrature_point_indices()) 
          { 


            for (const auto i : fe_values.dof_indices()) 
              { 
                const SymmetricTensor<2, dim> displacement_phi_i_symmgrad = 
                  fe_values[displacements<dim>].symmetric_gradient(i, q_point); 
                const double displacement_phi_i_div = 
                  fe_values[displacements<dim>].divergence(i, q_point); 

                const SymmetricTensor<2, dim> 
                  displacement_multiplier_phi_i_symmgrad = 
                    fe_values[displacement_multipliers<dim>].symmetric_gradient( 
                      i, q_point); 
                const double displacement_multiplier_phi_i_div = 
                  fe_values[displacement_multipliers<dim>].divergence(i, 
                                                                      q_point); 

                const double density_phi_i = 
                  fe_values[densities<dim>].value(i, q_point); 
                const double unfiltered_density_phi_i = 
                  fe_values[unfiltered_densities<dim>].value(i, q_point); 
                const double unfiltered_density_multiplier_phi_i = 
                  fe_values[unfiltered_density_multipliers<dim>].value(i, 
                                                                       q_point); 

                const double lower_slack_multiplier_phi_i = 
                  fe_values[density_lower_slack_multipliers<dim>].value( 
                    i, q_point); 

                const double lower_slack_phi_i = 
                  fe_values[density_lower_slacks<dim>].value(i, q_point); 

                const double upper_slack_phi_i = 
                  fe_values[density_upper_slacks<dim>].value(i, q_point); 

                const double upper_slack_multiplier_phi_i = 
                  fe_values[density_upper_slack_multipliers<dim>].value( 
                    i, q_point); 

                for (const auto j : fe_values.dof_indices()) 
                  { 


                    const SymmetricTensor<2, dim> displacement_phi_j_symmgrad = 
                      fe_values[displacements<dim>].symmetric_gradient(j, 
                                                                       q_point); 
                    const double displacement_phi_j_div = 
                      fe_values[displacements<dim>].divergence(j, q_point); 

                    const SymmetricTensor<2, dim> 
                      displacement_multiplier_phi_j_symmgrad = 
                        fe_values[displacement_multipliers<dim>] 
                          .symmetric_gradient(j, q_point); 
                    const double displacement_multiplier_phi_j_div = 
                      fe_values[displacement_multipliers<dim>].divergence( 
                        j, q_point); 

                    const double density_phi_j = 
                      fe_values[densities<dim>].value(j, q_point); 

                    const double unfiltered_density_phi_j = 
                      fe_values[unfiltered_densities<dim>].value(j, q_point); 
                    const double unfiltered_density_multiplier_phi_j = 
                      fe_values[unfiltered_density_multipliers<dim>].value( 
                        j, q_point); 

                    const double lower_slack_phi_j = 
                      fe_values[density_lower_slacks<dim>].value(j, q_point); 

                    const double upper_slack_phi_j = 
                      fe_values[density_upper_slacks<dim>].value(j, q_point); 

                    const double lower_slack_multiplier_phi_j = 
                      fe_values[density_lower_slack_multipliers<dim>].value( 
                        j, q_point); 

                    const double upper_slack_multiplier_phi_j = 
                      fe_values[density_upper_slack_multipliers<dim>].value( 
                        j, q_point); 


                    /* 方程1  */ 
                    cell_matrix(i, j) += 
                      fe_values.JxW(q_point) * 
                      ( 
                        -density_phi_i * unfiltered_density_multiplier_phi_j 
                        + density_penalty_exponent * 
                            (density_penalty_exponent - 1) * 
                            std::pow(old_density_values[q_point], 
                                     density_penalty_exponent - 2) * 
                            density_phi_i * density_phi_j * 
                            (old_displacement_multiplier_divs[q_point] * 
                               old_displacement_divs[q_point] * 
                               lambda_values[q_point] + 
                             2 * mu_values[q_point] * 
                               (old_displacement_symmgrads[q_point] * 
                                old_displacement_multiplier_symmgrads[q_point])) 
                        + density_penalty_exponent * 
                            std::pow(old_density_values[q_point], 
                                     density_penalty_exponent - 1) * 
                            density_phi_i * 
                            (displacement_multiplier_phi_j_div * 
                               old_displacement_divs[q_point] * 
                               lambda_values[q_point] + 
                             2 * mu_values[q_point] * 
                               (old_displacement_symmgrads[q_point] * 
                                displacement_multiplier_phi_j_symmgrad)) 
                        + density_penalty_exponent * 
                            std::pow(old_density_values[q_point], 
                                     density_penalty_exponent - 1) * 
                            density_phi_i * 
                            (displacement_phi_j_div * 
                               old_displacement_multiplier_divs[q_point] * 
                               lambda_values[q_point] + 
                             2 * mu_values[q_point] * 
                               (old_displacement_multiplier_symmgrads[q_point] * 
                                displacement_phi_j_symmgrad))); 
                   
                    /* 方程2  */ 
                    cell_matrix(i, j) += 
                      fe_values.JxW(q_point) * 
                      (density_penalty_exponent * 
                         std::pow(old_density_values[q_point], 
                                  density_penalty_exponent - 1) * 
                         density_phi_j * 
                         (old_displacement_multiplier_divs[q_point] * 
                            displacement_phi_i_div * lambda_values[q_point] + 
                          2 * mu_values[q_point] * 
                            (old_displacement_multiplier_symmgrads[q_point] * 
                             displacement_phi_i_symmgrad)) 
                       + std::pow(old_density_values[q_point], 
                                  density_penalty_exponent) * 
                           (displacement_multiplier_phi_j_div * 
                              displacement_phi_i_div * lambda_values[q_point] + 
                            2 * mu_values[q_point] * 
                              (displacement_multiplier_phi_j_symmgrad * 
                               displacement_phi_i_symmgrad)) 
                      ); 

                   /*方程3，这与过滤器有关 */ 
                    cell_matrix(i, j) += 
                      fe_values.JxW(q_point) * 
                      (-1 * unfiltered_density_phi_i * 
                         lower_slack_multiplier_phi_j + 
                       unfiltered_density_phi_i * upper_slack_multiplier_phi_j); 

                     /* 方程4：原始可行性  */ 
                    cell_matrix(i, j) += 
                      fe_values.JxW(q_point) * 
                      ( 
                        density_penalty_exponent * 
                          std::pow(old_density_values[q_point], 
                                   density_penalty_exponent - 1) * 
                          density_phi_j * 
                          (old_displacement_divs[q_point] * 
                             displacement_multiplier_phi_i_div * 
                             lambda_values[q_point] + 
                           2 * mu_values[q_point] * 
                             (old_displacement_symmgrads[q_point] * 
                              displacement_multiplier_phi_i_symmgrad)) 

                        + std::pow(old_density_values[q_point], 
                                   density_penalty_exponent) * 
                            (displacement_phi_j_div * 
                               displacement_multiplier_phi_i_div * 
                               lambda_values[q_point] + 
                             2 * mu_values[q_point] * 
                               (displacement_phi_j_symmgrad * 
                                displacement_multiplier_phi_i_symmgrad))); 

                   /*等式5：原始可行性  */ 
                    cell_matrix(i, j) += 
                      -1 * fe_values.JxW(q_point) * 
                      lower_slack_multiplier_phi_i * 
                      (unfiltered_density_phi_j - lower_slack_phi_j); 
                  /* 等式6：原始可行性  */ 
                    cell_matrix(i, j) += 
                      -1 * fe_values.JxW(q_point) * 
                      upper_slack_multiplier_phi_i * 
                      (-1 * unfiltered_density_phi_j - upper_slack_phi_j); 
                    /* Equation 7: Primal feasibility - the part with the filter
                     * is added later */
                    cell_matrix(i, j) += -1 * fe_values.JxW(q_point) * 
                                         unfiltered_density_multiplier_phi_i * 
                                         (density_phi_j); 
                    /* Equation 8: Complementary slackness */
                    cell_matrix(i, j) += 
                      fe_values.JxW(q_point) * 
                      (lower_slack_phi_i * lower_slack_multiplier_phi_j 

                       + lower_slack_phi_i * lower_slack_phi_j * 
                           old_lower_slack_multiplier_values[q_point] / 
                           old_lower_slack_values[q_point]); 
                    /* Equation 9: Complementary slackness */
                    cell_matrix(i, j) += 
                      fe_values.JxW(q_point) * 
                      (upper_slack_phi_i * upper_slack_multiplier_phi_j 

                       + upper_slack_phi_i * upper_slack_phi_j * 
                           old_upper_slack_multiplier_values[q_point] / 
                           old_upper_slack_values[q_point]); 
                  } 
              } 
          } 


        MatrixTools::local_apply_boundary_values(boundary_values, 
                                                 local_dof_indices, 
                                                 cell_matrix, 
                                                 dummy_cell_rhs, 
                                                 true); 

        constraints.distribute_local_to_global(cell_matrix, 
                                               local_dof_indices, 
                                               system_matrix); 
      } 


    system_rhs = calculate_test_rhs(nonlinear_solution); 


    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        const unsigned int i = cell->active_cell_index(); 
        for (typename SparseMatrix<double>::iterator iter = 
               filter_matrix.begin(i); 
             iter != filter_matrix.end(i); 
             ++iter) 
          { 
            const unsigned int j     = iter->column(); 
            const double       value = iter->value() * cell->measure(); 

            system_matrix 
              .block(SolutionBlocks::unfiltered_density_multiplier, 
                     SolutionBlocks::unfiltered_density) 
              .add(i, j, value); 
            system_matrix 
              .block(SolutionBlocks::unfiltered_density, 
                     SolutionBlocks::unfiltered_density_multiplier) 
              .add(j, i, value); 
          } 
      } 
  } 


  template <int dim> 
  BlockVector<double> SANDTopOpt<dim>::solve() 
  { 
    TimerOutput::Scope t(timer, "solver"); 

    BlockVector<double> linear_solution; 
    linear_solution.reinit(nonlinear_solution); 

    SparseDirectUMFPACK A_direct; 
    A_direct.initialize(system_matrix); 
    A_direct.vmult(linear_solution, system_rhs); 

    constraints.distribute(linear_solution); 

    return linear_solution; 
  } 




  template <int dim> 
  std::pair<double, double> SANDTopOpt<dim>::calculate_max_step_size( 
    const BlockVector<double> &state, 
    const BlockVector<double> &step) const 
  { 
    double       fraction_to_boundary; 
    const double min_fraction_to_boundary = .8; 
    const double max_fraction_to_boundary = 1. - 1e-5; 

    if (min_fraction_to_boundary < 1 - barrier_size) 
      { 
        if (1 - barrier_size < max_fraction_to_boundary) 
          fraction_to_boundary = 1 - barrier_size; 
        else 
          fraction_to_boundary = max_fraction_to_boundary; 
      } 
    else 
      fraction_to_boundary = min_fraction_to_boundary; 

    double step_size_s_low  = 0; 
    double step_size_z_low  = 0; 
    double step_size_s_high = 1; 
    double step_size_z_high = 1; 
    double step_size_s, step_size_z; 

    const int max_bisection_method_steps = 50; 
    for (unsigned int k = 0; k < max_bisection_method_steps; ++k) 
      { 
        step_size_s = (step_size_s_low + step_size_s_high) / 2; 
        step_size_z = (step_size_z_low + step_size_z_high) / 2; 

        const BlockVector<double> state_test_s = 
          (fraction_to_boundary * state) + (step_size_s * step); 

        const BlockVector<double> state_test_z = 
          (fraction_to_boundary * state) + (step_size_z * step); 

        const bool accept_s = 
          (state_test_s.block(SolutionBlocks::density_lower_slack) 
             .is_non_negative()) && 
          (state_test_s.block(SolutionBlocks::density_upper_slack) 
             .is_non_negative()); 
        const bool accept_z = 
          (state_test_z.block(SolutionBlocks::density_lower_slack_multiplier) 
             .is_non_negative()) && 
          (state_test_z.block(SolutionBlocks::density_upper_slack_multiplier) 
             .is_non_negative()); 

        if (accept_s) 
          step_size_s_low = step_size_s; 
        else 
          step_size_s_high = step_size_s; 

        if (accept_z) 
          step_size_z_low = step_size_z; 
        else 
          step_size_z_high = step_size_z; 
      } 

    return {step_size_s_low, step_size_z_low}; 
  } 



  template <int dim> 
  BlockVector<double> SANDTopOpt<dim>::calculate_test_rhs( 
    const BlockVector<double> &test_solution) const 
  { 


    BlockVector<double> test_rhs; 
    test_rhs.reinit(system_rhs); 

    MappingQGeneric<dim>  mapping(1); 
    const QGauss<dim>     quadrature_formula(fe.degree + 1); 
    const QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 
    FEValues<dim>         fe_values(mapping, 
                            fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 
    FEFaceValues<dim>     fe_face_values(mapping, 
                                     fe, 
                                     face_quadrature_formula, 
                                     update_values | update_quadrature_points | 
                                       update_normal_vectors | 
                                       update_JxW_values); 

    const unsigned int dofs_per_cell = fe.dofs_per_cell; 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    Vector<double>     cell_rhs(dofs_per_cell); 
    FullMatrix<double> dummy_cell_matrix(dofs_per_cell, dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    std::vector<double> lambda_values(n_q_points); 
    std::vector<double> mu_values(n_q_points); 

    const Functions::ConstantFunction<dim> lambda(1.), mu(1.); 
    std::vector<Tensor<1, dim>>            rhs_values(n_q_points); 

    BlockVector<double> filtered_unfiltered_density_solution = test_solution; 
    BlockVector<double> filter_adjoint_unfiltered_density_multiplier_solution = 
      test_solution; 
    filtered_unfiltered_density_solution.block( 
      SolutionBlocks::unfiltered_density) = 0; 
    filter_adjoint_unfiltered_density_multiplier_solution.block( 
      SolutionBlocks::unfiltered_density_multiplier) = 0; 

    filter_matrix.vmult(filtered_unfiltered_density_solution.block( 
                          SolutionBlocks::unfiltered_density), 
                        test_solution.block( 
                          SolutionBlocks::unfiltered_density)); 
    filter_matrix.Tvmult( 
      filter_adjoint_unfiltered_density_multiplier_solution.block( 
        SolutionBlocks::unfiltered_density_multiplier), 
      test_solution.block(SolutionBlocks::unfiltered_density_multiplier)); 

    std::vector<double>                  old_density_values(n_q_points); 
    std::vector<Tensor<1, dim>>          old_displacement_values(n_q_points); 
    std::vector<double>                  old_displacement_divs(n_q_points); 
    std::vector<SymmetricTensor<2, dim>> old_displacement_symmgrads(n_q_points); 
    std::vector<Tensor<1, dim>> old_displacement_multiplier_values(n_q_points); 
    std::vector<double>         old_displacement_multiplier_divs(n_q_points); 
    std::vector<SymmetricTensor<2, dim>> old_displacement_multiplier_symmgrads( 
      n_q_points); 
    std::vector<double> old_lower_slack_multiplier_values(n_q_points); 
    std::vector<double> old_upper_slack_multiplier_values(n_q_points); 
    std::vector<double> old_lower_slack_values(n_q_points); 
    std::vector<double> old_upper_slack_values(n_q_points); 
    std::vector<double> old_unfiltered_density_values(n_q_points); 
    std::vector<double> old_unfiltered_density_multiplier_values(n_q_points); 
    std::vector<double> filtered_unfiltered_density_values(n_q_points); 
    std::vector<double> filter_adjoint_unfiltered_density_multiplier_values( 
      n_q_points); 

    using namespace ValueExtractors; 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_rhs = 0; 

        cell->get_dof_indices(local_dof_indices); 

        fe_values.reinit(cell); 

 
        mu.value_list(fe_values.get_quadrature_points(), mu_values); 

        fe_values[densities<dim>].get_function_values(test_solution, 
                                                      old_density_values); 
        fe_values[displacements<dim>].get_function_values( 
          test_solution, old_displacement_values); 
        fe_values[displacements<dim>].get_function_divergences( 
          test_solution, old_displacement_divs); 
        fe_values[displacements<dim>].get_function_symmetric_gradients( 
          test_solution, old_displacement_symmgrads); 
        fe_values[displacement_multipliers<dim>].get_function_values( 
          test_solution, old_displacement_multiplier_values); 
        fe_values[displacement_multipliers<dim>].get_function_divergences( 
          test_solution, old_displacement_multiplier_divs); 
        fe_values[displacement_multipliers<dim>] 
          .get_function_symmetric_gradients( 
            test_solution, old_displacement_multiplier_symmgrads); 
        fe_values[density_lower_slacks<dim>].get_function_values( 
          test_solution, old_lower_slack_values); 
        fe_values[density_lower_slack_multipliers<dim>].get_function_values( 
          test_solution, old_lower_slack_multiplier_values); 
        fe_values[density_upper_slacks<dim>].get_function_values( 
          test_solution, old_upper_slack_values); 
        fe_values[density_upper_slack_multipliers<dim>].get_function_values( 
          test_solution, old_upper_slack_multiplier_values); 
        fe_values[unfiltered_densities<dim>].get_function_values( 
          test_solution, old_unfiltered_density_values); 
        fe_values[unfiltered_density_multipliers<dim>].get_function_values( 
          test_solution, old_unfiltered_density_multiplier_values); 
        fe_values[unfiltered_densities<dim>].get_function_values( 
          filtered_unfiltered_density_solution, 
          filtered_unfiltered_density_values); 
        fe_values[unfiltered_density_multipliers<dim>].get_function_values( 
          filter_adjoint_unfiltered_density_multiplier_solution, 
          filter_adjoint_unfiltered_density_multiplier_values); 

        for (const auto q_point : fe_values.quadrature_point_indices()) 
          { 
            for (const auto i : fe_values.dof_indices()) 
              { 
                const SymmetricTensor<2, dim> displacement_phi_i_symmgrad = 
                  fe_values[displacements<dim>].symmetric_gradient(i, q_point); 
                const double displacement_phi_i_div = 
                  fe_values[displacements<dim>].divergence(i, q_point); 

                const SymmetricTensor<2, dim> 
                  displacement_multiplier_phi_i_symmgrad = 
                    fe_values[displacement_multipliers<dim>].symmetric_gradient( 
                      i, q_point); 
                const double displacement_multiplier_phi_i_div = 
                  fe_values[displacement_multipliers<dim>].divergence(i, 
                                                                      q_point); 

                const double density_phi_i = 
                  fe_values[densities<dim>].value(i, q_point); 
                const double unfiltered_density_phi_i = 
                  fe_values[unfiltered_densities<dim>].value(i, q_point); 
                const double unfiltered_density_multiplier_phi_i = 
                  fe_values[unfiltered_density_multipliers<dim>].value(i, 
                                                                       q_point); 

                const double lower_slack_multiplier_phi_i = 
                  fe_values[density_lower_slack_multipliers<dim>].value( 
                    i, q_point); 

                const double lower_slack_phi_i = 
                  fe_values[density_lower_slacks<dim>].value(i, q_point); 

                const double upper_slack_phi_i = 
                  fe_values[density_upper_slacks<dim>].value(i, q_point); 

                const double upper_slack_multiplier_phi_i = 
                  fe_values[density_upper_slack_multipliers<dim>].value( 
                    i, q_point); 

                /* 方程1：这个方程以及方程
                 * 2 and 3, are the variational derivatives of the 
                 * Lagrangian with respect to the decision 
                 * variables - the density, displacement, and 
                 * unfiltered density. */ 


                cell_rhs(i) += 
                  -1 * fe_values.JxW(q_point) * 
                  (density_penalty_exponent * 
                     std::pow(old_density_values[q_point], 
                              density_penalty_exponent - 1) * 
                     density_phi_i * 
                     (old_displacement_multiplier_divs[q_point] * 
                        old_displacement_divs[q_point] * 
                        lambda_values[q_point] + 
                      2 * mu_values[q_point] * 
                        (old_displacement_symmgrads[q_point] * 
                         old_displacement_multiplier_symmgrads[q_point])) - 
                   density_phi_i * 
                     old_unfiltered_density_multiplier_values[q_point]); 

                /*方程2；边界项将被进一步添加。
                 * below. */ 


                cell_rhs(i) += 
                  -1 * fe_values.JxW(q_point) * 
                  (std::pow(old_density_values[q_point], 
                            density_penalty_exponent) * 
                   (old_displacement_multiplier_divs[q_point] * 
                      displacement_phi_i_div * lambda_values[q_point] + 
                    2 * mu_values[q_point] * 
                      (old_displacement_multiplier_symmgrads[q_point] * 
                       displacement_phi_i_symmgrad))); 
               /* 方程3  */ 
                cell_rhs(i) += 
                  -1 * fe_values.JxW(q_point) * 
                  (unfiltered_density_phi_i * 
                     filter_adjoint_unfiltered_density_multiplier_values 
                       [q_point] + 
                   unfiltered_density_phi_i * 
                     old_upper_slack_multiplier_values[q_point] + 
                   -1 * unfiltered_density_phi_i * 
                     old_lower_slack_multiplier_values[q_point]); 

               /* 方程4；边界项将再次被处理。with below. 
                * This equation being driven to 0 ensures that the elasticity 
                * equation is met as a constraint. */ 
                cell_rhs(i) += -1 * fe_values.JxW(q_point) * 
                               (std::pow(old_density_values[q_point], 
                                         density_penalty_exponent) * 
                                (old_displacement_divs[q_point] * 
                                   displacement_multiplier_phi_i_div * 
                                   lambda_values[q_point] + 
                                 2 * mu_values[q_point] * 
                                   (displacement_multiplier_phi_i_symmgrad * 
                                    old_displacement_symmgrads[q_point]))); 

                /* 方程5：该方程设定了下限的松弛量， giving a minimum density of 0. */ 
                cell_rhs(i) += fe_values.JxW(q_point) * 
                               (lower_slack_multiplier_phi_i * 
                                (old_unfiltered_density_values[q_point] - 
                                 old_lower_slack_values[q_point])); 

                /* 方程6：该方程设定了上层松弛量variable equal to one minus the unfiltered density. */ 
                cell_rhs(i) += fe_values.JxW(q_point) * 
                               (upper_slack_multiplier_phi_i * 
                                (1 - old_unfiltered_density_values[q_point] - 
                                 old_upper_slack_values[q_point])); 

                /*等式7：这是在
                 * density and the filter applied to the 
                 * unfiltered density. This being driven to 0 by 
                 * the Newton steps ensures that the filter is 
                 * applied correctly. */ 
                cell_rhs(i) += fe_values.JxW(q_point) * 
                               (unfiltered_density_multiplier_phi_i * 
                                (old_density_values[q_point] - 
                                 filtered_unfiltered_density_values[q_point])); 

                /*方程8：这与方程9一起给出了
                 * requirement that $s*z = \alpha$ for the barrier 
                 * size alpha, and gives complementary slackness 
                 * from KKT conditions when $\alpha$ goes to 0. */ 
                cell_rhs(i) += 
                  -1 * fe_values.JxW(q_point) * 
                  (lower_slack_phi_i * 
                   (old_lower_slack_multiplier_values[q_point] - 
                    barrier_size / old_lower_slack_values[q_point])); 

                /*方程9  */ 
                cell_rhs(i) += 
                  -1 * fe_values.JxW(q_point) * 
                  (upper_slack_phi_i * 
                   (old_upper_slack_multiplier_values[q_point] - 
                    barrier_size / old_upper_slack_values[q_point])); 
              } 
          } 

        for (const auto &face : cell->face_iterators()) 
          { 
            if (face->at_boundary() && 
                face->boundary_id() == BoundaryIds::down_force) 
              { 
                fe_face_values.reinit(cell, face); 

                for (const auto face_q_point : 
                     fe_face_values.quadrature_point_indices()) 
                  { 
                    for (const auto i : fe_face_values.dof_indices()) 
                      { 
                        Tensor<1, dim> traction; 
                        traction[1] = -1.; 

                        cell_rhs(i) += 
                          -1 * 
                          (traction * fe_face_values[displacements<dim>].value( 
                                        i, face_q_point)) * 
                          fe_face_values.JxW(face_q_point); 

                        cell_rhs(i) += 
                          (traction * 
                           fe_face_values[displacement_multipliers<dim>].value( 
                             i, face_q_point)) * 
                          fe_face_values.JxW(face_q_point); 
                      } 
                  } 
              } 
          } 

        MatrixTools::local_apply_boundary_values(boundary_values, 
                                                 local_dof_indices, 
                                                 dummy_cell_matrix, 
                                                 cell_rhs, 
                                                 true); 

        constraints.distribute_local_to_global(cell_rhs, 
                                               local_dof_indices, 
                                               test_rhs); 
      } 

    return test_rhs; 
  } 



  template <int dim> 
  double SANDTopOpt<dim>::calculate_exact_merit( 
    const BlockVector<double> &test_solution) 
  { 
    TimerOutput::Scope t(timer, "merit function"); 

    double objective_function_merit = 0; 
    { 
      MappingQGeneric<dim>  mapping(1); 
      const QGauss<dim>     quadrature_formula(fe.degree + 1); 
      const QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 
      FEValues<dim>         fe_values(mapping, 
                              fe, 
                              quadrature_formula, 
                              update_values | update_gradients | 
                                update_quadrature_points | update_JxW_values); 
      FEFaceValues<dim>     fe_face_values(mapping, 
                                       fe, 
                                       face_quadrature_formula, 
                                       update_values | 
                                         update_quadrature_points | 
                                         update_normal_vectors | 
                                         update_JxW_values); 

      const unsigned int n_face_q_points = face_quadrature_formula.size(); 

      std::vector<Tensor<1, dim>> displacement_face_values(n_face_q_points); 

      for (const auto &cell : dof_handler.active_cell_iterators()) 
        { 
          for (const auto &face : cell->face_iterators()) 
            { 
              if (face->at_boundary() && 
                  face->boundary_id() == BoundaryIds::down_force) 
                { 
                  fe_face_values.reinit(cell, face); 
                  fe_face_values[ValueExtractors::displacements<dim>] 
                    .get_function_values(test_solution, 
                                         displacement_face_values); 
                  for (unsigned int face_q_point = 0; 
                       face_q_point < n_face_q_points; 
                       ++face_q_point) 
                    { 
                      Tensor<1, dim> traction; 
                      traction[1] = -1.; 

                      objective_function_merit += 
                        (traction * displacement_face_values[face_q_point]) * 
                        fe_face_values.JxW(face_q_point); 
                    } 
                } 
            } 
        } 
    } 

    for (const auto &cell : triangulation.active_cell_iterators()) 
      { 
        objective_function_merit = 
          objective_function_merit - 
          barrier_size * cell->measure() * 
            std::log(test_solution.block( 
              SolutionBlocks::density_lower_slack)[cell->active_cell_index()]); 
        objective_function_merit = 
          objective_function_merit - 
          barrier_size * cell->measure() * 
            std::log(test_solution.block( 
              SolutionBlocks::density_upper_slack)[cell->active_cell_index()]); 
      } 
    const BlockVector<double> test_rhs = calculate_test_rhs(test_solution); 

    const double elasticity_constraint_merit = 
      penalty_multiplier * 
      test_rhs.block(SolutionBlocks::displacement_multiplier).l1_norm(); 
    const double filter_constraint_merit = 
      penalty_multiplier * 
      test_rhs.block(SolutionBlocks::unfiltered_density_multiplier).l1_norm(); 
    const double lower_slack_merit = 
      penalty_multiplier * 
      test_rhs.block(SolutionBlocks::density_lower_slack_multiplier).l1_norm(); 
    const double upper_slack_merit = 
      penalty_multiplier * 
      test_rhs.block(SolutionBlocks::density_upper_slack_multiplier).l1_norm(); 

    const double total_merit = 
      objective_function_merit + elasticity_constraint_merit + 
      filter_constraint_merit + lower_slack_merit + upper_slack_merit; 
    return total_merit; 
  } 




  template <int dim> 
  BlockVector<double> SANDTopOpt<dim>::find_max_step() 
  { 
    assemble_system(); 
    BlockVector<double> step = solve(); 



    const std::vector<unsigned int> decision_variables = { 
      SolutionBlocks::density, 
      SolutionBlocks::displacement, 
      SolutionBlocks::unfiltered_density, 
      SolutionBlocks::density_upper_slack, 
      SolutionBlocks::density_lower_slack}; 
    double hess_part = 0; 
    double grad_part = 0; 
    for (const unsigned int decision_variable_i : decision_variables) 
      { 
        for (const unsigned int decision_variable_j : decision_variables) 
          { 
            Vector<double> temp_vector(step.block(decision_variable_i).size()); 
            system_matrix.block(decision_variable_i, decision_variable_j) 
              .vmult(temp_vector, step.block(decision_variable_j)); 
            hess_part += step.block(decision_variable_i) * temp_vector; 
          } 
        grad_part -= system_rhs.block(decision_variable_i) * 
                     step.block(decision_variable_i); 
      } 

    const std::vector<unsigned int> equality_constraint_multipliers = { 
      SolutionBlocks::displacement_multiplier, 
      SolutionBlocks::unfiltered_density_multiplier, 
      SolutionBlocks::density_lower_slack_multiplier, 
      SolutionBlocks::density_upper_slack_multiplier}; 
    double constraint_norm = 0; 
    for (unsigned int multiplier_i : equality_constraint_multipliers) 
      constraint_norm += system_rhs.block(multiplier_i).linfty_norm(); 

    double test_penalty_multiplier; 
    if (hess_part > 0) 
      test_penalty_multiplier = 
        (grad_part + .5 * hess_part) / (.05 * constraint_norm); 
    else 
      test_penalty_multiplier = (grad_part) / (.05 * constraint_norm); 

    penalty_multiplier = std::max(penalty_multiplier, test_penalty_multiplier); 


    const std::pair<double, double> max_step_sizes = 
      calculate_max_step_size(nonlinear_solution, step); 
    const double step_size_s = max_step_sizes.first; 
    const double step_size_z = max_step_sizes.second; 

    step.block(SolutionBlocks::density) *= step_size_s; 
    step.block(SolutionBlocks::displacement) *= step_size_s; 
    step.block(SolutionBlocks::unfiltered_density) *= step_size_s; 
    step.block(SolutionBlocks::displacement_multiplier) *= step_size_z; 
    step.block(SolutionBlocks::unfiltered_density_multiplier) *= step_size_z; 
    step.block(SolutionBlocks::density_lower_slack) *= step_size_s; 
    step.block(SolutionBlocks::density_lower_slack_multiplier) *= step_size_z; 
    step.block(SolutionBlocks::density_upper_slack) *= step_size_s; 
    step.block(SolutionBlocks::density_upper_slack_multiplier) *= step_size_z; 

    return step; 
  } 



  template <int dim> 
  BlockVector<double> 
  SANDTopOpt<dim>::compute_scaled_step(const BlockVector<double> &state, 
                                       const BlockVector<double> &max_step, 
                                       const double descent_requirement) 
  { 
    const double merit_derivative = 
      (calculate_exact_merit(state + 1e-4 * max_step) - 
       calculate_exact_merit(state)) / 
      1e-4; 
    double       step_size                 = 1; 
    unsigned int max_linesearch_iterations = 10; 
    for (unsigned int k = 0; k < max_linesearch_iterations; ++k) 
      { 
        if (calculate_exact_merit(state + step_size * max_step) < 
            calculate_exact_merit(state) + 
              step_size * descent_requirement * merit_derivative) 
          break; 
        else 
          step_size = step_size / 2; 
      } 
    return state + (step_size * max_step); 
  } 



  template <int dim> 
  bool SANDTopOpt<dim>::check_convergence(const BlockVector<double> &state) 
  { 
    const BlockVector<double> test_rhs      = calculate_test_rhs(state); 
    const double              test_rhs_norm = test_rhs.l1_norm(); 

    const double convergence_condition = 1e-2; 
    const double target_norm           = convergence_condition * barrier_size; 

    std::cout << "    Checking convergence. Current rhs norm is " 
              << test_rhs_norm << ", target is " << target_norm << std::endl; 

    return (test_rhs_norm < target_norm); 
  } 



  template <int dim> 
  void SANDTopOpt<dim>::output_results(const unsigned int iteration) const 
  { 
    std::vector<std::string> solution_names(1, "density"); 
    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      data_component_interpretation( 
        1, DataComponentInterpretation::component_is_scalar); 
    for (unsigned int i = 0; i < dim; ++i) 
      { 
        solution_names.emplace_back("displacement"); 
        data_component_interpretation.push_back( 
          DataComponentInterpretation::component_is_part_of_vector); 
      } 
    solution_names.emplace_back("unfiltered_density"); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 
    for (unsigned int i = 0; i < dim; ++i) 
      { 
        solution_names.emplace_back("displacement_multiplier"); 
        data_component_interpretation.push_back( 
          DataComponentInterpretation::component_is_part_of_vector); 
      } 
    solution_names.emplace_back("unfiltered_density_multiplier"); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 
    solution_names.emplace_back("low_slack"); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 
    solution_names.emplace_back("low_slack_multiplier"); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 
    solution_names.emplace_back("high_slack"); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 
    solution_names.emplace_back("high_slack_multiplier"); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 

    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(nonlinear_solution, 
                             solution_names, 
                             DataOut<dim>::type_dof_data, 
                             data_component_interpretation); 
    data_out.build_patches(); 

    std::ofstream output("solution" + std::to_string(iteration) + ".vtu"); 
    data_out.write_vtu(output); 
  } 

  template <int dim> 
  void SANDTopOpt<dim>::write_as_stl() 
  { 
    static_assert(dim == 2, 
                  "This function is not implemented for anything " 
                  "other than the 2d case."); 

    std::ofstream stlfile; 
    stlfile.open("bridge.stl"); 

    stlfile << "solid bridge\n" << std::scientific; 
    double height = .25; 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        if (nonlinear_solution.block( 
              SolutionBlocks::density)[cell->active_cell_index()] > 0.5) 
          { 
            const Tensor<1, dim> edge_directions[2] = {cell->vertex(1) - 
                                                         cell->vertex(0), 
                                                       cell->vertex(2) - 
                                                         cell->vertex(0)}; 
            const Tensor<2, dim> edge_tensor( 
              {{edge_directions[0][0], edge_directions[0][1]}, 
               {edge_directions[1][0], edge_directions[1][1]}}); 
            const bool is_right_handed_cell = (determinant(edge_tensor) > 0); 

            if (is_right_handed_cell) 
              { 

               /*在z=0处写出一个边。  */ 
                stlfile << "   facet normal " << 0.000000e+00 << " " 
                        << 0.000000e+00 << " " << -1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(0)[0] << " " 
                        << cell->vertex(0)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 
                stlfile << "   facet normal " << 0.000000e+00 << " " 
                        << 0.000000e+00 << " " << -1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(3)[0] << " " 
                        << cell->vertex(3)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 

               /*在z=高度处写下一个边。  */  
                stlfile << "   facet normal " << 0.000000e+00 << " " 
                        << 0.000000e+00 << " " << 1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(0)[0] << " " 
                        << cell->vertex(0)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << height << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 
                stlfile << "   facet normal " << 0.000000e+00 << " " 
                        << 0.000000e+00 << " " << 1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(3)[0] << " " 
                        << cell->vertex(3)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << height << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 
              } 
            else /* The cell has a left-handed set up */ 
              { 
               /* 在z=0处写出一边。  */ 
                stlfile << "   facet normal " << 0.000000e+00 << " "
                        << 0.000000e+00 << " " << -1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(0)[0] << " " 
                        << cell->vertex(0)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 
                stlfile << "   facet normal " << 0.000000e+00 << " " 
                        << 0.000000e+00 << " " << -1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(3)[0] << " " 
                        << cell->vertex(3)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 

               /*在z=高度处写出一个边。  */ 
                stlfile << "   facet normal " << 0.000000e+00 << " " 
                        << 0.000000e+00 << " " << 1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(0)[0] << " " 
                        << cell->vertex(0)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << height << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 
                stlfile << "   facet normal " << 0.000000e+00 << " " 
                        << 0.000000e+00 << " " << 1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(3)[0] << " " 
                        << cell->vertex(3)[1] << " " << height << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 
              } 


            for (unsigned int face_number = 0; 
                 face_number < GeometryInfo<dim>::faces_per_cell; 
                 ++face_number) 
              { 
                const typename DoFHandler<dim>::face_iterator face = 
                  cell->face(face_number); 

                if ((face->at_boundary()) || 
                    (!face->at_boundary() && 
                     (nonlinear_solution.block( 
                        0)[cell->neighbor(face_number)->active_cell_index()] < 
                      0.5))) 
                  { 
                    const Tensor<1, dim> normal_vector = 
                      (face->center() - cell->center()); 
                    const double normal_norm = normal_vector.norm(); 
                    if ((face->vertex(0)[0] - face->vertex(0)[0]) * 
                            (face->vertex(1)[1] - face->vertex(0)[1]) * 
                            0.000000e+00 + 
                          (face->vertex(0)[1] - face->vertex(0)[1]) * (0 - 0) * 
                            normal_vector[0] + 
                          (height - 0) * 
                            (face->vertex(1)[0] - face->vertex(0)[0]) * 
                            normal_vector[1] - 
                          (face->vertex(0)[0] - face->vertex(0)[0]) * (0 - 0) * 
                            normal_vector[1] - 
                          (face->vertex(0)[1] - face->vertex(0)[1]) * 
                            (face->vertex(1)[0] - face->vertex(0)[0]) * 
                            normal_vector[0] - 
                          (height - 0) * 
                            (face->vertex(1)[1] - face->vertex(0)[1]) * 0 > 
                        0) 
                      { 
                        stlfile << "   facet normal " 
                                << normal_vector[0] / normal_norm << " " 
                                << normal_vector[1] / normal_norm << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "      outer loop\n"; 
                        stlfile << "         vertex " << face->vertex(0)[0] 
                                << " " << face->vertex(0)[1] << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "         vertex " << face->vertex(0)[0] 
                                << " " << face->vertex(0)[1] << " " << height 
                                << "\n"; 
                        stlfile << "         vertex " << face->vertex(1)[0] 
                                << " " << face->vertex(1)[1] << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "      endloop\n"; 
                        stlfile << "   endfacet\n"; 
                        stlfile << "   facet normal " 
                                << normal_vector[0] / normal_norm << " " 
                                << normal_vector[1] / normal_norm << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "      outer loop\n"; 
                        stlfile << "         vertex " << face->vertex(0)[0] 
                                << " " << face->vertex(0)[1] << " " << height 
                                << "\n"; 
                        stlfile << "         vertex " << face->vertex(1)[0] 
                                << " " << face->vertex(1)[1] << " " << height 
                                << "\n"; 
                        stlfile << "         vertex " << face->vertex(1)[0] 
                                << " " << face->vertex(1)[1] << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "      endloop\n"; 
                        stlfile << "   endfacet\n"; 
                      } 
                    else 
                      { 
                        stlfile << "   facet normal " 
                                << normal_vector[0] / normal_norm << " " 
                                << normal_vector[1] / normal_norm << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "      outer loop\n"; 
                        stlfile << "         vertex " << face->vertex(0)[0] 
                                << " " << face->vertex(0)[1] << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "         vertex " << face->vertex(1)[0] 
                                << " " << face->vertex(1)[1] << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "         vertex " << face->vertex(0)[0] 
                                << " " << face->vertex(0)[1] << " " << height 
                                << "\n"; 
                        stlfile << "      endloop\n"; 
                        stlfile << "   endfacet\n"; 
                        stlfile << "   facet normal " 
                                << normal_vector[0] / normal_norm << " " 
                                << normal_vector[1] / normal_norm << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "      outer loop\n"; 
                        stlfile << "         vertex " << face->vertex(0)[0] 
                                << " " << face->vertex(0)[1] << " " << height 
                                << "\n"; 
                        stlfile << "         vertex " << face->vertex(1)[0] 
                                << " " << face->vertex(1)[1] << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "         vertex " << face->vertex(1)[0] 
                                << " " << face->vertex(1)[1] << " " << height 
                                << "\n"; 
                        stlfile << "      endloop\n"; 
                        stlfile << "   endfacet\n"; 
                      } 
                  } 
              } 
          } 
      } 
    stlfile << "endsolid bridge"; 
  } 




  template <int dim> 
  void SANDTopOpt<dim>::run() 
  { 
    std::cout << "filter r is: " << filter_r << std::endl; 

    { 
      TimerOutput::Scope t(timer, "setup"); 

      create_triangulation(); 

      dof_handler.distribute_dofs(fe); 
      DoFRenumbering::component_wise(dof_handler); 

      setup_boundary_values(); 
      setup_block_system(); 
      setup_filter_matrix(); 
    } 


    barrier_size                  = 25; 
    const double min_barrier_size = .0005; 

    const unsigned int max_uphill_steps    = 8; 
    const double       descent_requirement = .0001; 


    unsigned int       iteration_number = 0; 
    const unsigned int max_iterations   = 10000; 

    do 
      { 
        std::cout << "Starting outer step in iteration " << iteration_number 
                  << " with barrier parameter " << barrier_size << std::endl; 



        do 
          { 
            std::cout << "  Starting inner step in iteration " 
                      << iteration_number 
                      << " with merit function penalty multiplier " 
                      << penalty_multiplier << std::endl; 

            bool watchdog_step_found = false; 

            const BlockVector<double> watchdog_state = nonlinear_solution; 
            BlockVector<double>       first_step; 
            double target_merit     = numbers::signaling_nan<double>(); 
            double merit_derivative = numbers::signaling_nan<double>(); 

            for (unsigned int k = 0; k < max_uphill_steps; ++k) 
              { 
                ++iteration_number; 
                const BlockVector<double> update_step = find_max_step(); 

                if (k == 0) 
                  { 
                    first_step = update_step; 
                    merit_derivative = 
                      ((calculate_exact_merit(watchdog_state + 
                                              .0001 * first_step) - 
                        calculate_exact_merit(watchdog_state)) / 
                       .0001); 
                    target_merit = calculate_exact_merit(watchdog_state) + 
                                   descent_requirement * merit_derivative; 
                  } 

                nonlinear_solution += update_step; 
                const double current_merit = 
                  calculate_exact_merit(nonlinear_solution); 

                std::cout << "    current watchdog state merit is: " 
                          << current_merit << "; target merit is " 
                          << target_merit << std::endl; 

                if (current_merit < target_merit) 
                  { 
                    watchdog_step_found = true; 
                    std::cout << "    found workable step after " << k + 1 
                              << " iterations" << std::endl; 
                    break; 
                  } 
              } 

            if (watchdog_step_found == false) 
              { 
                ++iteration_number; 
                const BlockVector<double> update_step = find_max_step(); 
                const BlockVector<double> stretch_state = 
                  compute_scaled_step(nonlinear_solution, 
                                      update_step, 
                                      descent_requirement); 


                if ((calculate_exact_merit(nonlinear_solution) < 
                     calculate_exact_merit(watchdog_state)) || 
                    (calculate_exact_merit(stretch_state) < target_merit)) 
                  { 
                    std::cout << "    Taking scaled step from end of watchdog" 
                              << std::endl; 
                    nonlinear_solution = stretch_state; 
                  } 
                else 
                  { 
                    std::cout 
                      << "    Taking scaled step from beginning of watchdog" 
                      << std::endl; 
                    if (calculate_exact_merit(stretch_state) > 
                        calculate_exact_merit(watchdog_state)) 
                      { 
                        nonlinear_solution = 
                          compute_scaled_step(watchdog_state, 
                                              first_step, 
                                              descent_requirement); 
                      } 
                    else 
                      { 
                        ++iteration_number; 
                        nonlinear_solution = stretch_state; 
                        const BlockVector<double> stretch_step = 
                          find_max_step(); 
                        nonlinear_solution = 
                          compute_scaled_step(nonlinear_solution, 
                                              stretch_step, 
                                              descent_requirement); 
                      } 
                  } 
              } 

            output_results(iteration_number); 
          } 
        while ((iteration_number < max_iterations) && 
               (check_convergence(nonlinear_solution) == false)); 


        const double barrier_size_multiplier = .8; 
        const double barrier_size_exponent   = 1.2; 

        barrier_size = 
          std::max(std::min(barrier_size * barrier_size_multiplier, 
                            std::pow(barrier_size, barrier_size_exponent)), 
                   min_barrier_size); 

        std::cout << std::endl; 
      } 
    while (((barrier_size > min_barrier_size) || 
            (check_convergence(nonlinear_solution) == false)) && 
           (iteration_number < max_iterations)); 

    write_as_stl(); 
    timer.print_summary(); 
  } 
} // namespace SAND 


int main() 
{ 
  try 
    { 
      SAND::SANDTopOpt<2> elastic_problem_2d; 
      elastic_problem_2d.run(); 
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
                << "Aborting!" << std::endl;
      return 1;
    }
  return 0;
}
