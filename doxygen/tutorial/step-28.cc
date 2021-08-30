

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2009 - 2021 by the deal.II authors 
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
 * Author: Yaqi Wang, Texas A&M University, 2009, 2010 
 */ 




#include <deal.II/base/timer.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/parameter_handler.h> 
#include <deal.II/base/thread_management.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparsity_pattern.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/grid_generator.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

#include <fstream> 
#include <iostream> 


#include <deal.II/lac/block_vector.h> 


#include <deal.II/numerics/solution_transfer.h> 


#include <deal.II/grid/grid_tools.h> 


#include <boost/io/ios_state.hpp> 


#include <list> 
#include <iomanip> 


namespace Step28 
{ 
  using namespace dealii; 





  class MaterialData 
  { 
  public: 
    MaterialData(const unsigned int n_groups); 

    double get_diffusion_coefficient(const unsigned int group, 
                                     const unsigned int material_id) const; 
    double get_removal_XS(const unsigned int group, 
                          const unsigned int material_id) const; 
    double get_fission_XS(const unsigned int group, 
                          const unsigned int material_id) const; 
    double get_fission_dist_XS(const unsigned int group_1, 
                               const unsigned int group_2, 
                               const unsigned int material_id) const; 
    double get_scattering_XS(const unsigned int group_1, 
                             const unsigned int group_2, 
                             const unsigned int material_id) const; 
    double get_fission_spectrum(const unsigned int group, 
                                const unsigned int material_id) const; 

  private: 
    const unsigned int n_groups; 
    const unsigned int n_materials; 

    Table<2, double> diffusion; 
    Table<2, double> sigma_r; 
    Table<2, double> nu_sigma_f; 
    Table<3, double> sigma_s; 
    Table<2, double> chi; 
  }; 



  MaterialData::MaterialData(const unsigned int n_groups) 
    : n_groups(n_groups) 
    , n_materials(8) 
    , diffusion(n_materials, n_groups) 
    , sigma_r(n_materials, n_groups) 
    , nu_sigma_f(n_materials, n_groups) 
    , sigma_s(n_materials, n_groups, n_groups) 
    , chi(n_materials, n_groups) 
  { 
    switch (this->n_groups) 
      { 
        case 2: 
          { 
            for (unsigned int m = 0; m < n_materials; ++m) 
              { 
                diffusion[m][0] = 1.2; 
                diffusion[m][1] = 0.4; 
                chi[m][0]       = 1.0; 
                chi[m][1]       = 0.0; 
                sigma_r[m][0]   = 0.03; 
                for (unsigned int group_1 = 0; group_1 < n_groups; ++group_1) 
                  for (unsigned int group_2 = 0; group_2 < n_groups; ++group_2) 
                    sigma_s[m][group_1][group_2] = 0.0; 
              } 

            diffusion[5][1] = 0.2; 

            sigma_r[4][0] = 0.026; 
            sigma_r[5][0] = 0.051; 
            sigma_r[6][0] = 0.026; 
            sigma_r[7][0] = 0.050; 

            sigma_r[0][1] = 0.100; 
            sigma_r[1][1] = 0.200; 
            sigma_r[2][1] = 0.250; 
            sigma_r[3][1] = 0.300; 
            sigma_r[4][1] = 0.020; 
            sigma_r[5][1] = 0.040; 
            sigma_r[6][1] = 0.020; 
            sigma_r[7][1] = 0.800; 

            nu_sigma_f[0][0] = 0.0050; 
            nu_sigma_f[1][0] = 0.0075; 
            nu_sigma_f[2][0] = 0.0075; 
            nu_sigma_f[3][0] = 0.0075; 
            nu_sigma_f[4][0] = 0.000; 
            nu_sigma_f[5][0] = 0.000; 
            nu_sigma_f[6][0] = 1e-7; 
            nu_sigma_f[7][0] = 0.00; 

            nu_sigma_f[0][1] = 0.125; 
            nu_sigma_f[1][1] = 0.300; 
            nu_sigma_f[2][1] = 0.375; 
            nu_sigma_f[3][1] = 0.450; 
            nu_sigma_f[4][1] = 0.000; 
            nu_sigma_f[5][1] = 0.000; 
            nu_sigma_f[6][1] = 3e-6; 
            nu_sigma_f[7][1] = 0.00; 

            sigma_s[0][0][1] = 0.020; 
            sigma_s[1][0][1] = 0.015; 
            sigma_s[2][0][1] = 0.015; 
            sigma_s[3][0][1] = 0.015; 
            sigma_s[4][0][1] = 0.025; 
            sigma_s[5][0][1] = 0.050; 
            sigma_s[6][0][1] = 0.025; 
            sigma_s[7][0][1] = 0.010; 

            break; 
          } 

        default: 
          Assert(false, 
                 ExcMessage( 
                   "Presently, only data for 2 groups is implemented")); 
      } 
  } 


  double 
  MaterialData::get_diffusion_coefficient(const unsigned int group, 
                                          const unsigned int material_id) const 
  { 
    Assert(group < n_groups, ExcIndexRange(group, 0, n_groups)); 
    Assert(material_id < n_materials, 
           ExcIndexRange(material_id, 0, n_materials)); 

    return diffusion[material_id][group]; 
  } 

  double MaterialData::get_removal_XS(const unsigned int group, 
                                      const unsigned int material_id) const 
  { 
    Assert(group < n_groups, ExcIndexRange(group, 0, n_groups)); 
    Assert(material_id < n_materials, 
           ExcIndexRange(material_id, 0, n_materials)); 

    return sigma_r[material_id][group]; 
  } 

  double MaterialData::get_fission_XS(const unsigned int group, 
                                      const unsigned int material_id) const 
  { 
    Assert(group < n_groups, ExcIndexRange(group, 0, n_groups)); 
    Assert(material_id < n_materials, 
           ExcIndexRange(material_id, 0, n_materials)); 

    return nu_sigma_f[material_id][group]; 
  } 

  double MaterialData::get_scattering_XS(const unsigned int group_1, 
                                         const unsigned int group_2, 
                                         const unsigned int material_id) const 
  { 
    Assert(group_1 < n_groups, ExcIndexRange(group_1, 0, n_groups)); 
    Assert(group_2 < n_groups, ExcIndexRange(group_2, 0, n_groups)); 
    Assert(material_id < n_materials, 
           ExcIndexRange(material_id, 0, n_materials)); 

    return sigma_s[material_id][group_1][group_2]; 
  } 

  double 
  MaterialData::get_fission_spectrum(const unsigned int group, 
                                     const unsigned int material_id) const 
  { 
    Assert(group < n_groups, ExcIndexRange(group, 0, n_groups)); 
    Assert(material_id < n_materials, 
           ExcIndexRange(material_id, 0, n_materials)); 

    return chi[material_id][group]; 
  } 


  double MaterialData::get_fission_dist_XS(const unsigned int group_1, 
                                           const unsigned int group_2, 
                                           const unsigned int material_id) const 
  { 
    return (get_fission_spectrum(group_1, material_id) * 
            get_fission_XS(group_2, material_id)); 
  } 










  template <int dim> 
  class EnergyGroup 
  { 
  public: 



    EnergyGroup(const unsigned int        group, 
                const MaterialData &      material_data, 
                const Triangulation<dim> &coarse_grid, 
                const FiniteElement<dim> &fe); 

    void setup_linear_system(); 

    unsigned int n_active_cells() const; 
    unsigned int n_dofs() const; 


    void assemble_system_matrix(); 
    void assemble_ingroup_rhs(const Function<dim> &extraneous_source); 
    void assemble_cross_group_rhs(const EnergyGroup<dim> &g_prime); 


    void solve(); 

    double get_fission_source() const; 

    void output_results(const unsigned int cycle) const; 

    void estimate_errors(Vector<float> &error_indicators) const; 

    void refine_grid(const Vector<float> &error_indicators, 
                     const double         refine_threshold, 
                     const double         coarsen_threshold); 


  public: 
    Vector<double> solution; 
    Vector<double> solution_old; 



  private: 
    const unsigned int  group; 
    const MaterialData &material_data; 

    Triangulation<dim>        triangulation; 
    const FiniteElement<dim> &fe; 
    DoFHandler<dim>           dof_handler; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> system_rhs; 

    std::map<types::global_dof_index, double> boundary_values; 
    AffineConstraints<double>                 hanging_node_constraints; 


  private: 
    void assemble_cross_group_rhs_recursive( 
      const EnergyGroup<dim> &                       g_prime, 
      const typename DoFHandler<dim>::cell_iterator &cell_g, 
      const typename DoFHandler<dim>::cell_iterator &cell_g_prime, 
      const FullMatrix<double> &                     prolongation_matrix); 
  }; 


  template <int dim> 
  EnergyGroup<dim>::EnergyGroup(const unsigned int        group, 
                                const MaterialData &      material_data, 
                                const Triangulation<dim> &coarse_grid, 
                                const FiniteElement<dim> &fe) 
    : group(group) 
    , material_data(material_data) 
    , fe(fe) 
    , dof_handler(triangulation) 
  { 
    triangulation.copy_triangulation(coarse_grid); 
    dof_handler.distribute_dofs(fe); 
  } 

  template <int dim> 
  unsigned int EnergyGroup<dim>::n_active_cells() const 
  { 
    return triangulation.n_active_cells(); 
  } 

  template <int dim> 
  unsigned int EnergyGroup<dim>::n_dofs() const 
  { 
    return dof_handler.n_dofs(); 
  } 



  template <int dim> 
  void EnergyGroup<dim>::setup_linear_system() 
  { 
    const unsigned int n_dofs = dof_handler.n_dofs(); 

    hanging_node_constraints.clear(); 
    DoFTools::make_hanging_node_constraints(dof_handler, 
                                            hanging_node_constraints); 
    hanging_node_constraints.close(); 

    system_matrix.clear(); 

    DynamicSparsityPattern dsp(n_dofs, n_dofs); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 
    hanging_node_constraints.condense(dsp); 
    sparsity_pattern.copy_from(dsp); 

    system_matrix.reinit(sparsity_pattern); 

    system_rhs.reinit(n_dofs); 

    if (solution.size() == 0) 
      { 
        solution.reinit(n_dofs); 
        solution_old.reinit(n_dofs); 
        solution_old = 1.0; 
        solution     = solution_old; 
      } 




    boundary_values.clear(); 

    for (unsigned int i = 0; i < dim; ++i) 
      VectorTools::interpolate_boundary_values(dof_handler, 
                                               2 * i + 1, 
                                               Functions::ZeroFunction<dim>(), 
                                               boundary_values); 
  } 




  template <int dim> 
  void EnergyGroup<dim>::assemble_system_matrix() 
  { 
    const QGauss<dim> quadrature_formula(fe.degree + 1); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_matrix = 0; 

        fe_values.reinit(cell); 

        const double diffusion_coefficient = 
          material_data.get_diffusion_coefficient(group, cell->material_id()); 
        const double removal_XS = 
          material_data.get_removal_XS(group, cell->material_id()); 

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            for (unsigned int j = 0; j < dofs_per_cell; ++j) 
              cell_matrix(i, j) += 
                ((diffusion_coefficient * fe_values.shape_grad(i, q_point) * 
                    fe_values.shape_grad(j, q_point) + 
                  removal_XS * fe_values.shape_value(i, q_point) * 
                    fe_values.shape_value(j, q_point)) * 
                 fe_values.JxW(q_point)); 

        cell->get_dof_indices(local_dof_indices); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
            system_matrix.add(local_dof_indices[i], 
                              local_dof_indices[j], 
                              cell_matrix(i, j)); 
      } 

    hanging_node_constraints.condense(system_matrix); 
  } 



  template <int dim> 
  void 
  EnergyGroup<dim>::assemble_ingroup_rhs(const Function<dim> &extraneous_source) 
  { 
    system_rhs.reinit(dof_handler.n_dofs()); 

    const QGauss<dim> quadrature_formula(fe.degree + 1); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_quadrature_points | 
                              update_JxW_values); 

    Vector<double>      cell_rhs(dofs_per_cell); 
    std::vector<double> extraneous_source_values(n_q_points); 
    std::vector<double> solution_old_values(n_q_points); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_rhs = 0; 

        fe_values.reinit(cell); 

        const double fission_dist_XS = 
          material_data.get_fission_dist_XS(group, group, cell->material_id()); 

        extraneous_source.value_list(fe_values.get_quadrature_points(), 
                                     extraneous_source_values); 

        fe_values.get_function_values(solution_old, solution_old_values); 

        cell->get_dof_indices(local_dof_indices); 

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            cell_rhs(i) += 
              ((extraneous_source_values[q_point] + 
                fission_dist_XS * solution_old_values[q_point]) * 
               fe_values.shape_value(i, q_point) * fe_values.JxW(q_point)); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          system_rhs(local_dof_indices[i]) += cell_rhs(i); 
      } 
  } 




  template <int dim> 
  void 
  EnergyGroup<dim>::assemble_cross_group_rhs(const EnergyGroup<dim> &g_prime) 
  { 
    if (group == g_prime.group) 
      return; 

    const std::list<std::pair<typename DoFHandler<dim>::cell_iterator, 
                              typename DoFHandler<dim>::cell_iterator>> 
      cell_list = 
        GridTools::get_finest_common_cells(dof_handler, g_prime.dof_handler); 

    for (const auto &cell_pair : cell_list) 
      { 
        FullMatrix<double> unit_matrix(fe.n_dofs_per_cell()); 
        for (unsigned int i = 0; i < unit_matrix.m(); ++i) 
          unit_matrix(i, i) = 1; 
        assemble_cross_group_rhs_recursive(g_prime, 
                                           cell_pair.first, 
                                           cell_pair.second, 
                                           unit_matrix); 
      } 
  } 





  template <int dim> 
  void EnergyGroup<dim>::assemble_cross_group_rhs_recursive( 
    const EnergyGroup<dim> &                       g_prime, 
    const typename DoFHandler<dim>::cell_iterator &cell_g, 
    const typename DoFHandler<dim>::cell_iterator &cell_g_prime, 
    const FullMatrix<double> &                     prolongation_matrix) 
  { 


    if (!cell_g->has_children() && !cell_g_prime->has_children()) 
      { 
        const QGauss<dim>  quadrature_formula(fe.degree + 1); 
        const unsigned int n_q_points = quadrature_formula.size(); 

        FEValues<dim> fe_values(fe, 
                                quadrature_formula, 
                                update_values | update_JxW_values); 

        if (cell_g->level() > cell_g_prime->level()) 
          fe_values.reinit(cell_g); 
        else 
          fe_values.reinit(cell_g_prime); 

        const double fission_dist_XS = 
          material_data.get_fission_dist_XS(group, 
                                            g_prime.group, 
                                            cell_g_prime->material_id()); 

        const double scattering_XS = 
          material_data.get_scattering_XS(g_prime.group, 
                                          group, 
                                          cell_g_prime->material_id()); 

        FullMatrix<double> local_mass_matrix_f(fe.n_dofs_per_cell(), 
                                               fe.n_dofs_per_cell()); 
        FullMatrix<double> local_mass_matrix_g(fe.n_dofs_per_cell(), 
                                               fe.n_dofs_per_cell()); 

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i) 
            for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j) 
              { 
                local_mass_matrix_f(i, j) += 
                  (fission_dist_XS * fe_values.shape_value(i, q_point) * 
                   fe_values.shape_value(j, q_point) * fe_values.JxW(q_point)); 
                local_mass_matrix_g(i, j) += 
                  (scattering_XS * fe_values.shape_value(i, q_point) * 
                   fe_values.shape_value(j, q_point) * fe_values.JxW(q_point)); 
              } 





        Vector<double> g_prime_new_values(fe.n_dofs_per_cell()); 
        Vector<double> g_prime_old_values(fe.n_dofs_per_cell()); 
        cell_g_prime->get_dof_values(g_prime.solution_old, g_prime_old_values); 
        cell_g_prime->get_dof_values(g_prime.solution, g_prime_new_values); 

        Vector<double> cell_rhs(fe.n_dofs_per_cell()); 
        Vector<double> tmp(fe.n_dofs_per_cell()); 

        if (cell_g->level() > cell_g_prime->level()) 
          { 
            prolongation_matrix.vmult(tmp, g_prime_old_values); 
            local_mass_matrix_f.vmult(cell_rhs, tmp); 

            prolongation_matrix.vmult(tmp, g_prime_new_values); 
            local_mass_matrix_g.vmult_add(cell_rhs, tmp); 
          } 
        else 
          { 
            local_mass_matrix_f.vmult(tmp, g_prime_old_values); 
            prolongation_matrix.Tvmult(cell_rhs, tmp); 

            local_mass_matrix_g.vmult(tmp, g_prime_new_values); 
            prolongation_matrix.Tvmult_add(cell_rhs, tmp); 
          } 

        std::vector<types::global_dof_index> local_dof_indices( 
          fe.n_dofs_per_cell()); 
        cell_g->get_dof_indices(local_dof_indices); 

        for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i) 
          system_rhs(local_dof_indices[i]) += cell_rhs(i); 
      } 


    else 
      for (unsigned int child = 0; 
           child < GeometryInfo<dim>::max_children_per_cell; 
           ++child) 
        { 
          FullMatrix<double> new_matrix(fe.n_dofs_per_cell(), 
                                        fe.n_dofs_per_cell()); 
          fe.get_prolongation_matrix(child).mmult(new_matrix, 
                                                  prolongation_matrix); 

          if (cell_g->has_children()) 
            assemble_cross_group_rhs_recursive(g_prime, 
                                               cell_g->child(child), 
                                               cell_g_prime, 
                                               new_matrix); 
          else 
            assemble_cross_group_rhs_recursive(g_prime, 
                                               cell_g, 
                                               cell_g_prime->child(child), 
                                               new_matrix); 
        } 
  } 


  template <int dim> 
  double EnergyGroup<dim>::get_fission_source() const 
  { 
    const QGauss<dim>  quadrature_formula(fe.degree + 1); 
    const unsigned int n_q_points = quadrature_formula.size(); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_JxW_values); 

    std::vector<double> solution_values(n_q_points); 

    double fission_source = 0; 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_values.reinit(cell); 

        const double fission_XS = 
          material_data.get_fission_XS(group, cell->material_id()); 

        fe_values.get_function_values(solution, solution_values); 

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          fission_source += 
            (fission_XS * solution_values[q_point] * fe_values.JxW(q_point)); 
      } 

    return fission_source; 
  } 


  template <int dim> 
  void EnergyGroup<dim>::solve() 
  { 
    hanging_node_constraints.condense(system_rhs); 
    MatrixTools::apply_boundary_values(boundary_values, 
                                       system_matrix, 
                                       solution, 
                                       system_rhs); 

    SolverControl            solver_control(system_matrix.m(), 
                                 1e-12 * system_rhs.l2_norm()); 
    SolverCG<Vector<double>> cg(solver_control); 

    PreconditionSSOR<SparseMatrix<double>> preconditioner; 
    preconditioner.initialize(system_matrix, 1.2); 

    cg.solve(system_matrix, solution, system_rhs, preconditioner); 

    hanging_node_constraints.distribute(solution); 
  } 



  template <int dim> 
  void EnergyGroup<dim>::estimate_errors(Vector<float> &error_indicators) const 
  { 
    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      QGauss<dim - 1>(fe.degree + 1), 
      std::map<types::boundary_id, const Function<dim> *>(), 
      solution, 
      error_indicators); 
    error_indicators /= solution.linfty_norm(); 
  } 




  template <int dim> 
  void EnergyGroup<dim>::refine_grid(const Vector<float> &error_indicators, 
                                     const double         refine_threshold, 
                                     const double         coarsen_threshold) 
  { 
    for (const auto &cell : triangulation.active_cell_iterators()) 
      if (error_indicators(cell->active_cell_index()) > refine_threshold) 
        cell->set_refine_flag(); 
      else if (error_indicators(cell->active_cell_index()) < coarsen_threshold) 
        cell->set_coarsen_flag(); 

    SolutionTransfer<dim> soltrans(dof_handler); 

    triangulation.prepare_coarsening_and_refinement(); 
    soltrans.prepare_for_coarsening_and_refinement(solution); 

    triangulation.execute_coarsening_and_refinement(); 
    dof_handler.distribute_dofs(fe); 
    setup_linear_system(); 

    solution.reinit(dof_handler.n_dofs()); 
    soltrans.interpolate(solution_old, solution); 


    hanging_node_constraints.distribute(solution); 

    solution_old.reinit(dof_handler.n_dofs()); 
    solution_old = solution; 
  } 


  template <int dim> 
  void EnergyGroup<dim>::output_results(const unsigned int cycle) const 
  { 
    const std::string filename = std::string("solution-") + 
                                 Utilities::int_to_string(group, 2) + "." + 
                                 Utilities::int_to_string(cycle, 2) + ".vtu"; 

    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "solution"); 
    data_out.build_patches(); 

    std::ofstream output(filename); 
    data_out.write_vtu(output); 
  } 




  template <int dim> 
  class NeutronDiffusionProblem 
  { 
  public: 
    class Parameters 
    { 
    public: 
      Parameters(); 

      static void declare_parameters(ParameterHandler &prm); 
      void        get_parameters(ParameterHandler &prm); 

      unsigned int n_groups; 
      unsigned int n_refinement_cycles; 

      unsigned int fe_degree; 

      double convergence_tolerance; 
    }; 

    NeutronDiffusionProblem(const Parameters &parameters); 

    void run(); 

  private: 


    void initialize_problem(); 

    void refine_grid(); 

    double get_total_fission_source() const; 


    const Parameters & parameters; 
    const MaterialData material_data; 
    FE_Q<dim>          fe; 


    double k_eff; 


    std::vector<std::unique_ptr<EnergyGroup<dim>>> energy_groups; 


    std::ofstream convergence_table_stream; 
  }; 


  template <int dim> 
  NeutronDiffusionProblem<dim>::Parameters::Parameters() 
    : n_groups(2) 
    , n_refinement_cycles(5) 
    , fe_degree(2) 
    , convergence_tolerance(1e-12) 
  {} 

  template <int dim> 
  void NeutronDiffusionProblem<dim>::Parameters::declare_parameters( 
    ParameterHandler &prm) 
  { 
    prm.declare_entry("Number of energy groups", 
                      "2", 
                      Patterns::Integer(), 
                      "The number of energy different groups considered"); 
    prm.declare_entry("Refinement cycles", 
                      "5", 
                      Patterns::Integer(), 
                      "Number of refinement cycles to be performed"); 
    prm.declare_entry("Finite element degree", 
                      "2", 
                      Patterns::Integer(), 
                      "Polynomial degree of the finite element to be used"); 
    prm.declare_entry( 
      "Power iteration tolerance", 
      "1e-12", 
      Patterns::Double(), 
      "Inner power iterations are stopped when the change in k_eff falls " 
      "below this tolerance"); 
  } 

  template <int dim> 
  void NeutronDiffusionProblem<dim>::Parameters::get_parameters( 
    ParameterHandler &prm) 
  { 
    n_groups              = prm.get_integer("Number of energy groups"); 
    n_refinement_cycles   = prm.get_integer("Refinement cycles"); 
    fe_degree             = prm.get_integer("Finite element degree"); 
    convergence_tolerance = prm.get_double("Power iteration tolerance"); 
  } 



  template <int dim> 
  NeutronDiffusionProblem<dim>::NeutronDiffusionProblem( 
    const Parameters &parameters) 
    : parameters(parameters) 
    , material_data(parameters.n_groups) 
    , fe(parameters.fe_degree) 
    , k_eff(std::numeric_limits<double>::quiet_NaN()) 
  {} 




  template <int dim> 
  void NeutronDiffusionProblem<dim>::initialize_problem() 
  { 
    const unsigned int rods_per_assembly_x = 17, rods_per_assembly_y = 17; 
    const double       pin_pitch_x = 1.26, pin_pitch_y = 1.26; 
    const double       assembly_height = 200; 

    const unsigned int assemblies_x = 2, assemblies_y = 2, assemblies_z = 1; 

    const Point<dim> bottom_left = Point<dim>(); 
    const Point<dim> upper_right = 
      (dim == 2 ? Point<dim>(assemblies_x * rods_per_assembly_x * pin_pitch_x, 
                             assemblies_y * rods_per_assembly_y * pin_pitch_y) : 
                  Point<dim>(assemblies_x * rods_per_assembly_x * pin_pitch_x, 
                             assemblies_y * rods_per_assembly_y * pin_pitch_y, 
                             assemblies_z * assembly_height)); 

    std::vector<unsigned int> n_subdivisions; 
    n_subdivisions.push_back(assemblies_x * rods_per_assembly_x); 
    if (dim >= 2) 
      n_subdivisions.push_back(assemblies_y * rods_per_assembly_y); 
    if (dim >= 3) 
      n_subdivisions.push_back(assemblies_z); 

    Triangulation<dim> coarse_grid; 
    GridGenerator::subdivided_hyper_rectangle( 
      coarse_grid, n_subdivisions, bottom_left, upper_right, true); 




    const unsigned int n_assemblies = 4; 
    const unsigned int assembly_materials 
      [n_assemblies][rods_per_assembly_x][rods_per_assembly_y] = { 
        {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 5, 1, 1, 5, 1, 1, 7, 1, 1, 5, 1, 1, 5, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}}, 
        {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 8, 1, 1, 8, 1, 1, 7, 1, 1, 8, 1, 1, 8, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}}, 
        {{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}, 
         {2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2}, 
         {2, 3, 3, 3, 3, 5, 3, 3, 5, 3, 3, 5, 3, 3, 3, 3, 2}, 
         {2, 3, 3, 5, 3, 4, 4, 4, 4, 4, 4, 4, 3, 5, 3, 3, 2}, 
         {2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2}, 
         {2, 3, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 3, 2}, 
         {2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2}, 
         {2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2}, 
         {2, 3, 5, 4, 4, 5, 4, 4, 7, 4, 4, 5, 4, 4, 5, 3, 2}, 
         {2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2}, 
         {2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2}, 
         {2, 3, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 3, 2}, 
         {2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2}, 
         {2, 3, 3, 5, 3, 4, 4, 4, 4, 4, 4, 4, 3, 5, 3, 3, 2}, 
         {2, 3, 3, 3, 3, 5, 3, 3, 5, 3, 3, 5, 3, 3, 3, 3, 2}, 
         {2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2}, 
         {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}}, 
        {{6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}}}; 


    const unsigned int core[assemblies_x][assemblies_y][assemblies_z] = { 
      {{0}, {2}}, {{2}, {0}}}; 


    for (auto &cell : coarse_grid.active_cell_iterators()) 
      { 
        const Point<dim> cell_center = cell->center(); 

        const unsigned int tmp_x = int(cell_center[0] / pin_pitch_x); 
        const unsigned int ax    = tmp_x / rods_per_assembly_x; 
        const unsigned int cx    = tmp_x - ax * rods_per_assembly_x; 

        const unsigned     tmp_y = int(cell_center[1] / pin_pitch_y); 
        const unsigned int ay    = tmp_y / rods_per_assembly_y; 
        const unsigned int cy    = tmp_y - ay * rods_per_assembly_y; 

        const unsigned int az = 
          (dim == 2 ? 0 : int(cell_center[dim - 1] / assembly_height)); 

        Assert(ax < assemblies_x, ExcInternalError()); 
        Assert(ay < assemblies_y, ExcInternalError()); 
        Assert(az < assemblies_z, ExcInternalError()); 

        Assert(core[ax][ay][az] < n_assemblies, ExcInternalError()); 

        Assert(cx < rods_per_assembly_x, ExcInternalError()); 
        Assert(cy < rods_per_assembly_y, ExcInternalError()); 

        cell->set_material_id(assembly_materials[core[ax][ay][az]][cx][cy] - 1); 
      } 


    for (unsigned int group = 0; group < parameters.n_groups; ++group) 
      energy_groups.emplace_back(std::make_unique<EnergyGroup<dim>>( 
        group, material_data, coarse_grid, fe)); 
    convergence_table_stream.open("convergence_table"); 
    convergence_table_stream.precision(12); 
  } 




  template <int dim> 
  double NeutronDiffusionProblem<dim>::get_total_fission_source() const 
  { 
    std::vector<double>  fission_sources(parameters.n_groups); 
    Threads::TaskGroup<> tasks; 
    for (unsigned int group = 0; group < parameters.n_groups; ++group) 
      tasks += Threads::new_task<>([&, group]() { 
        fission_sources[group] = energy_groups[group]->get_fission_source(); 
      }); 
    tasks.join_all(); 

    return std::accumulate(fission_sources.begin(), fission_sources.end(), 0.0); 
  } 



  template <int dim> 
  void NeutronDiffusionProblem<dim>::refine_grid() 
  { 
    std::vector<types::global_dof_index> n_cells(parameters.n_groups); 
    for (unsigned int group = 0; group < parameters.n_groups; ++group) 
      n_cells[group] = energy_groups[group]->n_active_cells(); 

    BlockVector<float> group_error_indicators(n_cells); 

    { 
      Threads::TaskGroup<> tasks; 
      for (unsigned int group = 0; group < parameters.n_groups; ++group) 
        tasks += Threads::new_task([&, group]() { 
          energy_groups[group]->estimate_errors( 
            group_error_indicators.block(group)); 
        }); 
    } 


    const float max_error         = group_error_indicators.linfty_norm(); 
    const float refine_threshold  = 0.3 * max_error; 
    const float coarsen_threshold = 0.01 * max_error; 

    { 
      Threads::TaskGroup<void> tasks; 
      for (unsigned int group = 0; group < parameters.n_groups; ++group) 
        tasks += Threads::new_task([&, group]() { 
          energy_groups[group]->refine_grid(group_error_indicators.block(group), 
                                            refine_threshold, 
                                            coarsen_threshold); 
        }); 
    } 
  } 



  template <int dim> 
  void NeutronDiffusionProblem<dim>::run() 
  { 


    boost::io::ios_flags_saver restore_flags(std::cout); 
    std::cout << std::setprecision(12) << std::fixed; 


    double k_eff_old = 0.0; 

    for (unsigned int cycle = 0; cycle < parameters.n_refinement_cycles; 
         ++cycle) 
      { 


        Timer timer; 

        std::cout << "Cycle " << cycle << ':' << std::endl; 

        if (cycle == 0) 
          { 
            initialize_problem(); 
            for (unsigned int group = 0; group < parameters.n_groups; ++group) 
              energy_groups[group]->setup_linear_system(); 
          } 

        else 
          { 
            refine_grid(); 
            for (unsigned int group = 0; group < parameters.n_groups; ++group) 
              energy_groups[group]->solution *= k_eff; 
          } 

        std::cout << "   Numbers of active cells:       "; 
        for (unsigned int group = 0; group < parameters.n_groups; ++group) 
          std::cout << energy_groups[group]->n_active_cells() << ' '; 
        std::cout << std::endl; 
        std::cout << "   Numbers of degrees of freedom: "; 
        for (unsigned int group = 0; group < parameters.n_groups; ++group) 
          std::cout << energy_groups[group]->n_dofs() << ' '; 
        std::cout << std::endl << std::endl; 

        Threads::TaskGroup<> tasks; 
        for (unsigned int group = 0; group < parameters.n_groups; ++group) 
          tasks += Threads::new_task( 
            [&, group]() { energy_groups[group]->assemble_system_matrix(); }); 
        tasks.join_all(); 

        double       error; 
        unsigned int iteration = 1; 
        do 
          { 
            for (unsigned int group = 0; group < parameters.n_groups; ++group) 
              { 
                energy_groups[group]->assemble_ingroup_rhs( 
                  Functions::ZeroFunction<dim>()); 

                for (unsigned int bgroup = 0; bgroup < parameters.n_groups; 
                     ++bgroup) 
                  energy_groups[group]->assemble_cross_group_rhs( 
                    *energy_groups[bgroup]); 

                energy_groups[group]->solve(); 
              } 

            k_eff = get_total_fission_source(); 
            error = std::abs(k_eff - k_eff_old) / std::abs(k_eff); 
            const double flux_ratio = energy_groups[0]->solution.linfty_norm() / 
                                      energy_groups[1]->solution.linfty_norm(); 
            const double max_thermal = energy_groups[1]->solution.linfty_norm(); 
            std::cout << "Iter number:" << std::setw(2) << std::right 
                      << iteration << " k_eff=" << k_eff 
                      << " flux_ratio=" << flux_ratio 
                      << " max_thermal=" << max_thermal << std::endl; 
            k_eff_old = k_eff; 

            for (unsigned int group = 0; group < parameters.n_groups; ++group) 
              { 
                energy_groups[group]->solution_old = 
                  energy_groups[group]->solution; 
                energy_groups[group]->solution_old /= k_eff; 
              } 

            ++iteration; 
          } 
        while ((error > parameters.convergence_tolerance) && (iteration < 500)); 
        convergence_table_stream << cycle << " " << energy_groups[0]->n_dofs() 
                                 << " " << energy_groups[1]->n_dofs() << " " 
                                 << k_eff << " " 
                                 << energy_groups[0]->solution.linfty_norm() / 
                                      energy_groups[1]->solution.linfty_norm() 
                                 << '\n'; 

        for (unsigned int group = 0; group < parameters.n_groups; ++group) 
          energy_groups[group]->output_results(cycle); 


        std::cout << std::endl; 
        std::cout << "   Cycle=" << cycle << ", n_dofs=" 
                  << energy_groups[0]->n_dofs() + energy_groups[1]->n_dofs() 
                  << ",  k_eff=" << k_eff << ", time=" << timer.cpu_time() 
                  << std::endl; 

        std::cout << std::endl << std::endl; 
      } 
  } 
} // namespace Step28 




int main(int argc, char **argv) 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step28; 

      std::string filename; 
      if (argc < 2) 
        filename = "project.prm"; 
      else 
        filename = argv[1]; 

      const unsigned int dim = 2; 

      ParameterHandler parameter_handler; 

      NeutronDiffusionProblem<dim>::Parameters parameters; 
      parameters.declare_parameters(parameter_handler); 

      parameter_handler.parse_input(filename); 

      parameters.get_parameters(parameter_handler); 

      NeutronDiffusionProblem<dim> neutron_diffusion_problem(parameters); 
      neutron_diffusion_problem.run(); 
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

