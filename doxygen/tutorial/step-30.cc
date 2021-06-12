

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2007 - 2021 by the deal.II authors 
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
 * Author: Tobias Leicht, 2007 
 */ 




#include <deal.II/base/function.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/lac/precondition_block.h> 
#include <deal.II/lac/solver_richardson.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/mapping_q1.h> 
#include <deal.II/fe/fe_dgq.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/derivative_approximation.h> 


#include <array> 
#include <iostream> 
#include <fstream> 


namespace Step30 
{ 
  using namespace dealii; 


  template <int dim> 
  class RHS : public Function<dim> 
  { 
  public: 
    virtual void value_list(const std::vector<Point<dim>> &points, 
                            std::vector<double> &          values, 
                            const unsigned int /*component*/ = 0) const override 
    { 
      (void)points; 
      Assert(values.size() == points.size(), 
             ExcDimensionMismatch(values.size(), points.size())); 

      std::fill(values.begin(), values.end(), 0.); 
    } 
  }; 

  template <int dim> 
  class BoundaryValues : public Function<dim> 
  { 
  public: 
    virtual void value_list(const std::vector<Point<dim>> &points, 
                            std::vector<double> &          values, 
                            const unsigned int /*component*/ = 0) const override 
    { 
      Assert(values.size() == points.size(), 
             ExcDimensionMismatch(values.size(), points.size())); 

      for (unsigned int i = 0; i < values.size(); ++i) 
        { 
          if (points[i](0) < 0.5) 
            values[i] = 1.; 
          else 
            values[i] = 0.; 
        } 
    } 
  }; 

  template <int dim> 
  class Beta 
  { 
  public: 


    void value_list(const std::vector<Point<dim>> &points, 
                    std::vector<Point<dim>> &      values) const 
    { 
      Assert(values.size() == points.size(), 
             ExcDimensionMismatch(values.size(), points.size())); 

      for (unsigned int i = 0; i < points.size(); ++i) 
        { 
          if (points[i](0) > 0) 
            { 
              values[i](0) = -points[i](1); 
              values[i](1) = points[i](0); 
            } 
          else 
            { 
              values[i]    = Point<dim>(); 
              values[i](0) = -points[i](1); 
            } 
        } 
    } 
  }; 



  template <int dim> 
  class DGTransportEquation 
  { 
  public: 
    DGTransportEquation(); 

    void assemble_cell_term(const FEValues<dim> &fe_v, 
                            FullMatrix<double> & ui_vi_matrix, 
                            Vector<double> &     cell_vector) const; 

    void assemble_boundary_term(const FEFaceValues<dim> &fe_v, 
                                FullMatrix<double> &     ui_vi_matrix, 
                                Vector<double> &         cell_vector) const; 

    void assemble_face_term(const FEFaceValuesBase<dim> &fe_v, 
                            const FEFaceValuesBase<dim> &fe_v_neighbor, 
                            FullMatrix<double> &         ui_vi_matrix, 
                            FullMatrix<double> &         ue_vi_matrix, 
                            FullMatrix<double> &         ui_ve_matrix, 
                            FullMatrix<double> &         ue_ve_matrix) const; 

  private: 
    const Beta<dim>           beta_function; 
    const RHS<dim>            rhs_function; 
    const BoundaryValues<dim> boundary_function; 
  }; 


  template <int dim> 
  DGTransportEquation<dim>::DGTransportEquation() 
    : beta_function() 
    , rhs_function() 
    , boundary_function() 
  {} 

  template <int dim> 
  void DGTransportEquation<dim>::assemble_cell_term( 
    const FEValues<dim> &fe_v, 
    FullMatrix<double> & ui_vi_matrix, 
    Vector<double> &     cell_vector) const 
  { 
    const std::vector<double> &JxW = fe_v.get_JxW_values(); 

    std::vector<Point<dim>> beta(fe_v.n_quadrature_points); 
    std::vector<double>     rhs(fe_v.n_quadrature_points); 

    beta_function.value_list(fe_v.get_quadrature_points(), beta); 
    rhs_function.value_list(fe_v.get_quadrature_points(), rhs); 

    for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point) 
      for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
        { 
          for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j) 
            ui_vi_matrix(i, j) -= beta[point] * fe_v.shape_grad(i, point) * 
                                  fe_v.shape_value(j, point) * JxW[point]; 

          cell_vector(i) += 
            rhs[point] * fe_v.shape_value(i, point) * JxW[point]; 
        } 
  } 

  template <int dim> 
  void DGTransportEquation<dim>::assemble_boundary_term( 
    const FEFaceValues<dim> &fe_v, 
    FullMatrix<double> &     ui_vi_matrix, 
    Vector<double> &         cell_vector) const 
  { 
    const std::vector<double> &        JxW     = fe_v.get_JxW_values(); 
    const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors(); 

    std::vector<Point<dim>> beta(fe_v.n_quadrature_points); 
    std::vector<double>     g(fe_v.n_quadrature_points); 

    beta_function.value_list(fe_v.get_quadrature_points(), beta); 
    boundary_function.value_list(fe_v.get_quadrature_points(), g); 

    for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point) 
      { 
        const double beta_n = beta[point] * normals[point]; 
        if (beta_n > 0) 
          for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
            for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j) 
              ui_vi_matrix(i, j) += beta_n * fe_v.shape_value(j, point) * 
                                    fe_v.shape_value(i, point) * JxW[point]; 
        else 
          for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
            cell_vector(i) -= 
              beta_n * g[point] * fe_v.shape_value(i, point) * JxW[point]; 
      } 
  } 

  template <int dim> 
  void DGTransportEquation<dim>::assemble_face_term( 
    const FEFaceValuesBase<dim> &fe_v, 
    const FEFaceValuesBase<dim> &fe_v_neighbor, 
    FullMatrix<double> &         ui_vi_matrix, 
    FullMatrix<double> &         ue_vi_matrix, 
    FullMatrix<double> &         ui_ve_matrix, 
    FullMatrix<double> &         ue_ve_matrix) const 
  { 
    const std::vector<double> &        JxW     = fe_v.get_JxW_values(); 
    const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors(); 

    std::vector<Point<dim>> beta(fe_v.n_quadrature_points); 

    beta_function.value_list(fe_v.get_quadrature_points(), beta); 

    for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point) 
      { 
        const double beta_n = beta[point] * normals[point]; 
        if (beta_n > 0) 
          { 
            for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
              for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j) 
                ui_vi_matrix(i, j) += beta_n * fe_v.shape_value(j, point) * 
                                      fe_v.shape_value(i, point) * JxW[point]; 

            for (unsigned int k = 0; k < fe_v_neighbor.dofs_per_cell; ++k) 
              for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j) 
                ui_ve_matrix(k, j) -= beta_n * fe_v.shape_value(j, point) * 
                                      fe_v_neighbor.shape_value(k, point) * 
                                      JxW[point]; 
          } 
        else 
          { 
            for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
              for (unsigned int l = 0; l < fe_v_neighbor.dofs_per_cell; ++l) 
                ue_vi_matrix(i, l) += beta_n * 
                                      fe_v_neighbor.shape_value(l, point) * 
                                      fe_v.shape_value(i, point) * JxW[point]; 

            for (unsigned int k = 0; k < fe_v_neighbor.dofs_per_cell; ++k) 
              for (unsigned int l = 0; l < fe_v_neighbor.dofs_per_cell; ++l) 
                ue_ve_matrix(k, l) -= 
                  beta_n * fe_v_neighbor.shape_value(l, point) * 
                  fe_v_neighbor.shape_value(k, point) * JxW[point]; 
          } 
      } 
  } 


  template <int dim> 
  class DGMethod 
  { 
  public: 
    DGMethod(const bool anisotropic); 

    void run(); 

  private: 
    void setup_system(); 
    void assemble_system(); 
    void solve(Vector<double> &solution); 
    void refine_grid(); 
    void set_anisotropic_flags(); 
    void output_results(const unsigned int cycle) const; 

    Triangulation<dim>   triangulation; 
    const MappingQ1<dim> mapping; 


    const unsigned int degree; 
    FE_DGQ<dim>        fe; 
    DoFHandler<dim>    dof_handler; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 


    const double anisotropic_threshold_ratio; 


    const bool anisotropic; 

    const QGauss<dim>     quadrature; 
    const QGauss<dim - 1> face_quadrature; 

    Vector<double> solution2; 
    Vector<double> right_hand_side; 

    const DGTransportEquation<dim> dg; 
  }; 

  template <int dim> 
  DGMethod<dim>::DGMethod(const bool anisotropic) 
    : mapping() 
    , 


    degree(1) 
    , fe(degree) 
    , dof_handler(triangulation) 
    , anisotropic_threshold_ratio(3.) 
    , anisotropic(anisotropic) 
    , 


    quadrature(degree + 1) 
    , face_quadrature(degree + 1) 
    , dg() 
  {} 

  template <int dim> 
  void DGMethod<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 
    sparsity_pattern.reinit(dof_handler.n_dofs(), 
                            dof_handler.n_dofs(), 
                            (GeometryInfo<dim>::faces_per_cell * 
                               GeometryInfo<dim>::max_children_per_face + 
                             1) * 
                              fe.n_dofs_per_cell()); 

    DoFTools::make_flux_sparsity_pattern(dof_handler, sparsity_pattern); 

    sparsity_pattern.compress(); 

    system_matrix.reinit(sparsity_pattern); 

    solution2.reinit(dof_handler.n_dofs()); 
    right_hand_side.reinit(dof_handler.n_dofs()); 
  } 


  template <int dim> 
  void DGMethod<dim>::assemble_system() 
  { 
    const unsigned int dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell(); 
    std::vector<types::global_dof_index> dofs(dofs_per_cell); 
    std::vector<types::global_dof_index> dofs_neighbor(dofs_per_cell); 

    const UpdateFlags update_flags = update_values | update_gradients | 
                                     update_quadrature_points | 
                                     update_JxW_values; 

    const UpdateFlags face_update_flags = 
      update_values | update_quadrature_points | update_JxW_values | 
      update_normal_vectors; 

    const UpdateFlags neighbor_face_update_flags = update_values; 

    FEValues<dim>        fe_v(mapping, fe, quadrature, update_flags); 
    FEFaceValues<dim>    fe_v_face(mapping, 
                                fe, 
                                face_quadrature, 
                                face_update_flags); 
    FESubfaceValues<dim> fe_v_subface(mapping, 
                                      fe, 
                                      face_quadrature, 
                                      face_update_flags); 
    FEFaceValues<dim>    fe_v_face_neighbor(mapping, 
                                         fe, 
                                         face_quadrature, 
                                         neighbor_face_update_flags); 

    FullMatrix<double> ui_vi_matrix(dofs_per_cell, dofs_per_cell); 
    FullMatrix<double> ue_vi_matrix(dofs_per_cell, dofs_per_cell); 

    FullMatrix<double> ui_ve_matrix(dofs_per_cell, dofs_per_cell); 
    FullMatrix<double> ue_ve_matrix(dofs_per_cell, dofs_per_cell); 

    Vector<double> cell_vector(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        ui_vi_matrix = 0; 
        cell_vector  = 0; 

        fe_v.reinit(cell); 

        dg.assemble_cell_term(fe_v, ui_vi_matrix, cell_vector); 

        cell->get_dof_indices(dofs); 

        for (const auto face_no : cell->face_indices()) 
          { 
            const auto face = cell->face(face_no); 


            if (face->at_boundary()) 
              { 
                fe_v_face.reinit(cell, face_no); 

                dg.assemble_boundary_term(fe_v_face, ui_vi_matrix, cell_vector); 
              } 
            else 
              { 
                Assert(cell->neighbor(face_no).state() == IteratorState::valid, 
                       ExcInternalError()); 
                const auto neighbor = cell->neighbor(face_no); 


                if (face->has_children()) 
                  { 


                    const unsigned int neighbor2 = 
                      cell->neighbor_face_no(face_no); 


                    for (unsigned int subface_no = 0; 
                         subface_no < face->n_active_descendants(); 
                         ++subface_no) 
                      { 


                        const auto neighbor_child = 
                          cell->neighbor_child_on_subface(face_no, subface_no); 
                        Assert(!neighbor_child->has_children(), 
                               ExcInternalError()); 


                        ue_vi_matrix = 0; 
                        ui_ve_matrix = 0; 
                        ue_ve_matrix = 0; 

                        fe_v_subface.reinit(cell, face_no, subface_no); 
                        fe_v_face_neighbor.reinit(neighbor_child, neighbor2); 

                        dg.assemble_face_term(fe_v_subface, 
                                              fe_v_face_neighbor, 
                                              ui_vi_matrix, 
                                              ue_vi_matrix, 
                                              ui_ve_matrix, 
                                              ue_ve_matrix); 

                        neighbor_child->get_dof_indices(dofs_neighbor); 

                        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                            { 
                              system_matrix.add(dofs[i], 
                                                dofs_neighbor[j], 
                                                ue_vi_matrix(i, j)); 
                              system_matrix.add(dofs_neighbor[i], 
                                                dofs[j], 
                                                ui_ve_matrix(i, j)); 
                              system_matrix.add(dofs_neighbor[i], 
                                                dofs_neighbor[j], 
                                                ue_ve_matrix(i, j)); 
                            } 
                      } 
                  } 
                else 
                  { 


                    if (dim > 1 && cell->neighbor_is_coarser(face_no)) 
                      continue; 


                    if (((dim > 1) && (cell->index() < neighbor->index())) || 
                        ((dim == 1) && ((cell->level() < neighbor->level()) || 
                                        ((cell->level() == neighbor->level()) && 
                                         (cell->index() < neighbor->index()))))) 
                      { 


                        const unsigned int neighbor2 = 
                          cell->neighbor_of_neighbor(face_no); 

                        ue_vi_matrix = 0; 
                        ui_ve_matrix = 0; 
                        ue_ve_matrix = 0; 

                        fe_v_face.reinit(cell, face_no); 
                        fe_v_face_neighbor.reinit(neighbor, neighbor2); 

                        dg.assemble_face_term(fe_v_face, 
                                              fe_v_face_neighbor, 
                                              ui_vi_matrix, 
                                              ue_vi_matrix, 
                                              ui_ve_matrix, 
                                              ue_ve_matrix); 

                        neighbor->get_dof_indices(dofs_neighbor); 

                        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                            { 
                              system_matrix.add(dofs[i], 
                                                dofs_neighbor[j], 
                                                ue_vi_matrix(i, j)); 
                              system_matrix.add(dofs_neighbor[i], 
                                                dofs[j], 
                                                ui_ve_matrix(i, j)); 
                              system_matrix.add(dofs_neighbor[i], 
                                                dofs_neighbor[j], 
                                                ue_ve_matrix(i, j)); 
                            } 
                      } 


                  } 
              } 
          } 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
            system_matrix.add(dofs[i], dofs[j], ui_vi_matrix(i, j)); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          right_hand_side(dofs[i]) += cell_vector(i); 
      } 
  } 


  template <int dim> 
  void DGMethod<dim>::solve(Vector<double> &solution) 
  { 
    SolverControl                    solver_control(1000, 1e-12, false, false); 
    SolverRichardson<Vector<double>> solver(solver_control); 

    PreconditionBlockSSOR<SparseMatrix<double>> preconditioner; 

    preconditioner.initialize(system_matrix, fe.n_dofs_per_cell()); 

    solver.solve(system_matrix, solution, right_hand_side, preconditioner); 
  } 


  template <int dim> 
  void DGMethod<dim>::refine_grid() 
  { 
    Vector<float> gradient_indicator(triangulation.n_active_cells()); 


    DerivativeApproximation::approximate_gradient(mapping, 
                                                  dof_handler, 
                                                  solution2, 
                                                  gradient_indicator); 


    for (const auto &cell : triangulation.active_cell_iterators()) 
      gradient_indicator[cell->active_cell_index()] *= 
        std::pow(cell->diameter(), 1 + 1.0 * dim / 2); 


    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    gradient_indicator, 
                                                    0.3, 
                                                    0.1); 


    if (anisotropic) 
      set_anisotropic_flags(); 


    triangulation.execute_coarsening_and_refinement(); 
  } 


  template <int dim> 
  void DGMethod<dim>::set_anisotropic_flags() 
  { 


    UpdateFlags face_update_flags = 
      UpdateFlags(update_values | update_JxW_values); 

    FEFaceValues<dim>    fe_v_face(mapping, 
                                fe, 
                                face_quadrature, 
                                face_update_flags); 
    FESubfaceValues<dim> fe_v_subface(mapping, 
                                      fe, 
                                      face_quadrature, 
                                      face_update_flags); 
    FEFaceValues<dim>    fe_v_face_neighbor(mapping, 
                                         fe, 
                                         face_quadrature, 
                                         update_values); 


    for (const auto &cell : dof_handler.active_cell_iterators()) 


      if (cell->refine_flag_set()) 
        { 
          Point<dim> jump; 
          Point<dim> area; 

          for (const auto face_no : cell->face_indices()) 
            { 
              const auto face = cell->face(face_no); 

              if (!face->at_boundary()) 
                { 
                  Assert(cell->neighbor(face_no).state() == 
                           IteratorState::valid, 
                         ExcInternalError()); 
                  const auto neighbor = cell->neighbor(face_no); 

                  std::vector<double> u(fe_v_face.n_quadrature_points); 
                  std::vector<double> u_neighbor(fe_v_face.n_quadrature_points); 


                  if (face->has_children()) 
                    { 


                      unsigned int neighbor2 = cell->neighbor_face_no(face_no); 


                      for (unsigned int subface_no = 0; 
                           subface_no < face->n_active_descendants(); 
                           ++subface_no) 
                        { 


                          const auto neighbor_child = 
                            cell->neighbor_child_on_subface(face_no, 
                                                            subface_no); 
                          Assert(!neighbor_child->has_children(), 
                                 ExcInternalError()); 


                          fe_v_subface.reinit(cell, face_no, subface_no); 
                          fe_v_face_neighbor.reinit(neighbor_child, neighbor2); 


                          fe_v_subface.get_function_values(solution2, u); 
                          fe_v_face_neighbor.get_function_values(solution2, 
                                                                 u_neighbor); 


                          const std::vector<double> &JxW = 
                            fe_v_subface.get_JxW_values(); 


                          for (unsigned int x = 0; 
                               x < fe_v_subface.n_quadrature_points; 
                               ++x) 
                            { 


                              jump[face_no / 2] += 
                                std::abs(u[x] - u_neighbor[x]) * JxW[x]; 


                              area[face_no / 2] += JxW[x]; 
                            } 
                        } 
                    } 
                  else 
                    { 
                      if (!cell->neighbor_is_coarser(face_no)) 
                        { 


                          unsigned int neighbor2 = 
                            cell->neighbor_of_neighbor(face_no); 

                          fe_v_face.reinit(cell, face_no); 
                          fe_v_face_neighbor.reinit(neighbor, neighbor2); 

                          fe_v_face.get_function_values(solution2, u); 
                          fe_v_face_neighbor.get_function_values(solution2, 
                                                                 u_neighbor); 

                          const std::vector<double> &JxW = 
                            fe_v_face.get_JxW_values(); 

                          for (unsigned int x = 0; 
                               x < fe_v_face.n_quadrature_points; 
                               ++x) 
                            { 
                              jump[face_no / 2] += 
                                std::abs(u[x] - u_neighbor[x]) * JxW[x]; 
                              area[face_no / 2] += JxW[x]; 
                            } 
                        } 
                      else // i.e. neighbor is coarser than cell 
                        { 


                          std::pair<unsigned int, unsigned int> 
                            neighbor_face_subface = 
                              cell->neighbor_of_coarser_neighbor(face_no); 
                          Assert(neighbor_face_subface.first < cell->n_faces(), 
                                 ExcInternalError()); 
                          Assert(neighbor_face_subface.second < 
                                   neighbor->face(neighbor_face_subface.first) 
                                     ->n_active_descendants(), 
                                 ExcInternalError()); 
                          Assert(neighbor->neighbor_child_on_subface( 
                                   neighbor_face_subface.first, 
                                   neighbor_face_subface.second) == cell, 
                                 ExcInternalError()); 

                          fe_v_face.reinit(cell, face_no); 
                          fe_v_subface.reinit(neighbor, 
                                              neighbor_face_subface.first, 
                                              neighbor_face_subface.second); 

                          fe_v_face.get_function_values(solution2, u); 
                          fe_v_subface.get_function_values(solution2, 
                                                           u_neighbor); 

                          const std::vector<double> &JxW = 
                            fe_v_face.get_JxW_values(); 

                          for (unsigned int x = 0; 
                               x < fe_v_face.n_quadrature_points; 
                               ++x) 
                            { 
                              jump[face_no / 2] += 
                                std::abs(u[x] - u_neighbor[x]) * JxW[x]; 
                              area[face_no / 2] += JxW[x]; 
                            } 
                        } 
                    } 
                } 
            } 


          std::array<double, dim> average_jumps; 
          double                  sum_of_average_jumps = 0.; 
          for (unsigned int i = 0; i < dim; ++i) 
            { 
              average_jumps[i] = jump(i) / area(i); 
              sum_of_average_jumps += average_jumps[i]; 
            } 


          for (unsigned int i = 0; i < dim; ++i) 
            if (average_jumps[i] > anisotropic_threshold_ratio * 
                                     (sum_of_average_jumps - average_jumps[i])) 
              cell->set_refine_flag(RefinementCase<dim>::cut_axis(i)); 
        } 
  } 


  template <int dim> 
  void DGMethod<dim>::output_results(const unsigned int cycle) const 
  { 
    std::string refine_type; 
    if (anisotropic) 
      refine_type = ".aniso"; 
    else 
      refine_type = ".iso"; 

    { 
      const std::string filename = 
        "grid-" + std::to_string(cycle) + refine_type + ".svg"; 
      std::cout << "   Writing grid to <" << filename << ">..." << std::endl; 
      std::ofstream svg_output(filename); 

      GridOut grid_out; 
      grid_out.write_svg(triangulation, svg_output); 
    } 

    { 
      const std::string filename = 
        "sol-" + std::to_string(cycle) + refine_type + ".vtu"; 
      std::cout << "   Writing solution to <" << filename << ">..." 
                << std::endl; 
      std::ofstream gnuplot_output(filename); 

      DataOut<dim> data_out; 
      data_out.attach_dof_handler(dof_handler); 
      data_out.add_data_vector(solution2, "u"); 

      data_out.build_patches(degree); 

      data_out.write_vtu(gnuplot_output); 
    } 
  } 

  template <int dim> 
  void DGMethod<dim>::run() 
  { 
    for (unsigned int cycle = 0; cycle < 6; ++cycle) 
      { 
        std::cout << "Cycle " << cycle << ':' << std::endl; 

        if (cycle == 0) 
          { 


            Point<dim> p1, p2; 
            p1(0) = 0; 
            p1(0) = -1; 
            for (unsigned int i = 0; i < dim; ++i) 
              p2(i) = 1.; 


            std::vector<unsigned int> repetitions(dim, 1); 
            repetitions[0] = 2; 
            GridGenerator::subdivided_hyper_rectangle(triangulation, 
                                                      repetitions, 
                                                      p1, 
                                                      p2); 

            triangulation.refine_global(5 - dim); 
          } 
        else 
          refine_grid(); 

        std::cout << "   Number of active cells:       " 
                  << triangulation.n_active_cells() << std::endl; 

        setup_system(); 

        std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
                  << std::endl; 

        Timer assemble_timer; 
        assemble_system(); 
        std::cout << "   Time of assemble_system: " << assemble_timer.cpu_time() 
                  << std::endl; 
        solve(solution2); 

        output_results(cycle); 

        std::cout << std::endl; 
      } 
  } 
} // namespace Step30 

int main() 
{ 
  try 
    { 
      using namespace Step30; 


      const unsigned int dim = 2; 

      { 


        std::cout << "Performing a " << dim 
                  << "D run with isotropic refinement..." << std::endl 
                  << "------------------------------------------------" 
                  << std::endl; 
        DGMethod<dim> dgmethod_iso(false); 
        dgmethod_iso.run(); 
      } 

      { 


        std::cout << std::endl 
                  << "Performing a " << dim 
                  << "D run with anisotropic refinement..." << std::endl 
                  << "--------------------------------------------------" 
                  << std::endl; 
        DGMethod<dim> dgmethod_aniso(true); 
        dgmethod_aniso.run(); 
      } 
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
    }; 

  return 0; 
} 

