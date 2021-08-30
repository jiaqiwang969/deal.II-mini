

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2018 - 2021 by the deal.II authors 
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
 *      Author: Zhuoran Wang, Colorado State University, 2018 
 */ 



#include <deal.II/base/quadrature.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/tensor_function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/point.h> 
#include <deal.II/lac/block_vector.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/block_sparse_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_dgq.h> 
#include <deal.II/fe/fe_raviart_thomas.h> 
#include <deal.II/fe/fe_dg_vector.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/fe_face.h> 
#include <deal.II/fe/component_mask.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/data_out_faces.h> 

#include <fstream> 
#include <iostream> 


namespace Step61 
{ 
  using namespace dealii; 



  template <int dim> 
  class WGDarcyEquation 
  { 
  public: 
    WGDarcyEquation(const unsigned int degree); 
    void run(); 

  private: 
    void make_grid(); 
    void setup_system(); 
    void assemble_system(); 
    void solve(); 
    void compute_postprocessed_velocity(); 
    void compute_velocity_errors(); 
    void compute_pressure_error(); 
    void output_results() const; 

    Triangulation<dim> triangulation; 

    FESystem<dim>   fe; 
    DoFHandler<dim> dof_handler; 

    AffineConstraints<double> constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> solution; 
    Vector<double> system_rhs; 

    FE_DGRaviartThomas<dim> fe_dgrt; 
    DoFHandler<dim>         dof_handler_dgrt; 
    Vector<double>          darcy_velocity; 
  }; 



  template <int dim> 
  class Coefficient : public TensorFunction<2, dim> 
  { 
  public: 
    Coefficient() 
      : TensorFunction<2, dim>() 
    {} 

    virtual void value_list(const std::vector<Point<dim>> &points, 
                            std::vector<Tensor<2, dim>> &values) const override; 
  }; 

  template <int dim> 
  void Coefficient<dim>::value_list(const std::vector<Point<dim>> &points, 
                                    std::vector<Tensor<2, dim>> &  values) const 
  { 
    Assert(points.size() == values.size(), 
           ExcDimensionMismatch(points.size(), values.size())); 
    for (unsigned int p = 0; p < points.size(); ++p) 
      values[p] = unit_symmetric_tensor<dim>(); 
  } 

  template <int dim> 
  class BoundaryValues : public Function<dim> 
  { 
  public: 
    BoundaryValues() 
      : Function<dim>(2) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 
  }; 

  template <int dim> 
  double BoundaryValues<dim>::value(const Point<dim> & /*p*/, 
                                    const unsigned int /*component*/) const 
  { 
    return 0; 
  } 

  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 
  };

  template <int dim> 
  double RightHandSide<dim>::value(const Point<dim> &p, 
                                   const unsigned int /*component*/) const 
  { 
    return (2 * numbers::PI * numbers::PI * std::sin(numbers::PI * p[0]) * 
            std::sin(numbers::PI * p[1])); 
  } 


  template <int dim> 
  class ExactPressure : public Function<dim> 
  { 
  public: 
    ExactPressure() 
      : Function<dim>(2) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component) const override; 
  }; 

  template <int dim> 
  double ExactPressure<dim>::value(const Point<dim> &p, 
                                   const unsigned int /*component*/) const 
  { 
    return std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]); 
  } 

  template <int dim> 
  class ExactVelocity : public TensorFunction<1, dim> 
  { 
  public: 
    ExactVelocity() 
      : TensorFunction<1, dim>() 
    {} 

    virtual Tensor<1, dim> value(const Point<dim> &p) const override; 
  }; 

  template <int dim> 
  Tensor<1, dim> ExactVelocity<dim>::value(const Point<dim> &p) const 
  { 
    Tensor<1, dim> return_value; 
    return_value[0] = -numbers::PI * std::cos(numbers::PI * p[0]) * 
                      std::sin(numbers::PI * p[1]); 
    return_value[1] = -numbers::PI * std::sin(numbers::PI * p[0]) * 
                      std::cos(numbers::PI * p[1]); 
    return return_value; 
  } 



  template <int dim> 
  WGDarcyEquation<dim>::WGDarcyEquation(const unsigned int degree) 
    : fe(FE_DGQ<dim>(degree), 1, FE_FaceQ<dim>(degree), 1) 
    , dof_handler(triangulation) 
    , fe_dgrt(degree) 
    , dof_handler_dgrt(triangulation) 
  {} 



  template <int dim> 
  void WGDarcyEquation<dim>::make_grid() 
  { 
    GridGenerator::hyper_cube(triangulation, 0, 1); 
    triangulation.refine_global(5); 

    std::cout << "   Number of active cells: " << triangulation.n_active_cells() 
              << std::endl 
              << "   Total number of cells: " << triangulation.n_cells() 
              << std::endl; 
  } 



  template <int dim> 
  void WGDarcyEquation<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 
    dof_handler_dgrt.distribute_dofs(fe_dgrt); 

    std::cout << "   Number of pressure degrees of freedom: " 
              << dof_handler.n_dofs() << std::endl; 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 

    { 
      constraints.clear(); 
      const FEValuesExtractors::Scalar interface_pressure(1); 
      const ComponentMask              interface_pressure_mask = 
        fe.component_mask(interface_pressure); 
      VectorTools::interpolate_boundary_values(dof_handler, 
                                               0, 
                                               BoundaryValues<dim>(), 
                                               constraints, 
                                               interface_pressure_mask); 
      constraints.close(); 
    } 


    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints); 
    sparsity_pattern.copy_from(dsp); 

    system_matrix.reinit(sparsity_pattern); 
  } 






  template <int dim> 
  void WGDarcyEquation<dim>::assemble_system() 
  { 
    const QGauss<dim>     quadrature_formula(fe_dgrt.degree + 1); 
    const QGauss<dim - 1> face_quadrature_formula(fe_dgrt.degree + 1); 

    FEValues<dim>     fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_quadrature_points | 
                              update_JxW_values); 
    FEFaceValues<dim> fe_face_values(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_normal_vectors | 
                                       update_quadrature_points | 
                                       update_JxW_values); 

    FEValues<dim>     fe_values_dgrt(fe_dgrt, 
                                 quadrature_formula, 
                                 update_values | update_gradients | 
                                   update_quadrature_points | 
                                   update_JxW_values); 
    FEFaceValues<dim> fe_face_values_dgrt(fe_dgrt, 
                                          face_quadrature_formula, 
                                          update_values | 
                                            update_normal_vectors | 
                                            update_quadrature_points | 
                                            update_JxW_values); 

    const unsigned int dofs_per_cell      = fe.n_dofs_per_cell(); 
    const unsigned int dofs_per_cell_dgrt = fe_dgrt.n_dofs_per_cell(); 

    const unsigned int n_q_points      = fe_values.get_quadrature().size(); 
    const unsigned int n_q_points_dgrt = fe_values_dgrt.get_quadrature().size(); 

    const unsigned int n_face_q_points = fe_face_values.get_quadrature().size(); 

    RightHandSide<dim>  right_hand_side; 
    std::vector<double> right_hand_side_values(n_q_points); 

    const Coefficient<dim>      coefficient; 
    std::vector<Tensor<2, dim>> coefficient_values(n_q_points); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 


    FullMatrix<double> cell_matrix_M(dofs_per_cell_dgrt, dofs_per_cell_dgrt); 
    FullMatrix<double> cell_matrix_G(dofs_per_cell_dgrt, dofs_per_cell); 
    FullMatrix<double> cell_matrix_C(dofs_per_cell, dofs_per_cell_dgrt); 
    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 
    Vector<double>     cell_solution(dofs_per_cell); 


    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure_interior(0); 
    const FEValuesExtractors::Scalar pressure_face(1); 


    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_values.reinit(cell); 

        const typename Triangulation<dim>::active_cell_iterator cell_dgrt = 
          cell; 
        fe_values_dgrt.reinit(cell_dgrt); 

        right_hand_side.value_list(fe_values.get_quadrature_points(), 
                                   right_hand_side_values); 
        coefficient.value_list(fe_values.get_quadrature_points(), 
                               coefficient_values); 


        cell_matrix_M = 0; 
        for (unsigned int q = 0; q < n_q_points_dgrt; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
            { 
              const Tensor<1, dim> v_i = fe_values_dgrt[velocities].value(i, q); 
              for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k) 
                { 
                  const Tensor<1, dim> v_k = 
                    fe_values_dgrt[velocities].value(k, q); 
                  cell_matrix_M(i, k) += (v_i * v_k * fe_values_dgrt.JxW(q)); 
                } 
            } 


        cell_matrix_M.gauss_jordan(); 


        cell_matrix_G = 0; 
        for (unsigned int q = 0; q < n_q_points; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
            { 
              const double div_v_i = 
                fe_values_dgrt[velocities].divergence(i, q); 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                { 
                  const double phi_j_interior = 
                    fe_values[pressure_interior].value(j, q); 

                  cell_matrix_G(i, j) -= 
                    (div_v_i * phi_j_interior * fe_values.JxW(q)); 
                } 
            } 


        for (const auto &face : cell->face_iterators()) 
          { 
            fe_face_values.reinit(cell, face); 
            fe_face_values_dgrt.reinit(cell_dgrt, face); 

            for (unsigned int q = 0; q < n_face_q_points; ++q) 
              { 
                const Tensor<1, dim> &normal = fe_face_values.normal_vector(q); 

                for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
                  { 
                    const Tensor<1, dim> v_i = 
                      fe_face_values_dgrt[velocities].value(i, q); 
                    for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                      { 
                        const double phi_j_face = 
                          fe_face_values[pressure_face].value(j, q); 

                        cell_matrix_G(i, j) += 
                          ((v_i * normal) * phi_j_face * fe_face_values.JxW(q)); 
                      } 
                  } 
              } 
          } 
        cell_matrix_G.Tmmult(cell_matrix_C, cell_matrix_M); 


        local_matrix = 0; 
        for (unsigned int q = 0; q < n_q_points_dgrt; ++q) 
          { 
            for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k) 
              { 
                const Tensor<1, dim> v_k = 
                  fe_values_dgrt[velocities].value(k, q); 
                for (unsigned int l = 0; l < dofs_per_cell_dgrt; ++l) 
                  { 
                    const Tensor<1, dim> v_l = 
                      fe_values_dgrt[velocities].value(l, q); 

                    for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                      for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                        local_matrix(i, j) += 
                          (coefficient_values[q] * cell_matrix_C[i][k] * v_k) * 
                          cell_matrix_C[j][l] * v_l * fe_values_dgrt.JxW(q); 
                  } 
              } 
          } 


        cell_rhs = 0; 
        for (unsigned int q = 0; q < n_q_points; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              cell_rhs(i) += (fe_values[pressure_interior].value(i, q) * 
                              right_hand_side_values[q] * fe_values.JxW(q)); 
            } 


        cell->get_dof_indices(local_dof_indices); 
        constraints.distribute_local_to_global( 
          local_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs); 
      } 
  } 



  template <int dim> 
  void WGDarcyEquation<dim>::solve() 
  { 
    SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm()); 
    SolverCG<Vector<double>> solver(solver_control); 
    solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity()); 
    constraints.distribute(solution); 
  } 



  template <int dim> 
  void WGDarcyEquation<dim>::compute_postprocessed_velocity() 
  { 
    darcy_velocity.reinit(dof_handler_dgrt.n_dofs()); 

    const QGauss<dim>     quadrature_formula(fe_dgrt.degree + 1); 
    const QGauss<dim - 1> face_quadrature_formula(fe_dgrt.degree + 1); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_quadrature_points | 
                              update_JxW_values); 

    FEFaceValues<dim> fe_face_values(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_normal_vectors | 
                                       update_quadrature_points | 
                                       update_JxW_values); 

    FEValues<dim> fe_values_dgrt(fe_dgrt, 
                                 quadrature_formula, 
                                 update_values | update_gradients | 
                                   update_quadrature_points | 
                                   update_JxW_values); 

    FEFaceValues<dim> fe_face_values_dgrt(fe_dgrt, 
                                          face_quadrature_formula, 
                                          update_values | 
                                            update_normal_vectors | 
                                            update_quadrature_points | 
                                            update_JxW_values); 

    const unsigned int dofs_per_cell      = fe.n_dofs_per_cell(); 
    const unsigned int dofs_per_cell_dgrt = fe_dgrt.n_dofs_per_cell(); 

    const unsigned int n_q_points      = fe_values.get_quadrature().size(); 
    const unsigned int n_q_points_dgrt = fe_values_dgrt.get_quadrature().size(); 

    const unsigned int n_face_q_points = fe_face_values.get_quadrature().size(); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
    std::vector<types::global_dof_index> local_dof_indices_dgrt( 
      dofs_per_cell_dgrt); 

    FullMatrix<double> cell_matrix_M(dofs_per_cell_dgrt, dofs_per_cell_dgrt); 
    FullMatrix<double> cell_matrix_G(dofs_per_cell_dgrt, dofs_per_cell); 
    FullMatrix<double> cell_matrix_C(dofs_per_cell, dofs_per_cell_dgrt); 
    FullMatrix<double> cell_matrix_D(dofs_per_cell_dgrt, dofs_per_cell_dgrt); 
    FullMatrix<double> cell_matrix_E(dofs_per_cell_dgrt, dofs_per_cell_dgrt); 

    Vector<double> cell_solution(dofs_per_cell); 
    Vector<double> cell_velocity(dofs_per_cell_dgrt); 

    const Coefficient<dim>      coefficient; 
    std::vector<Tensor<2, dim>> coefficient_values(n_q_points_dgrt); 

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure_interior(0); 
    const FEValuesExtractors::Scalar pressure_face(1); 


    typename DoFHandler<dim>::active_cell_iterator 
      cell = dof_handler.begin_active(), 
      endc = dof_handler.end(), cell_dgrt = dof_handler_dgrt.begin_active(); 
    for (; cell != endc; ++cell, ++cell_dgrt) 
      { 
        fe_values.reinit(cell); 
        fe_values_dgrt.reinit(cell_dgrt); 

        coefficient.value_list(fe_values_dgrt.get_quadrature_points(), 
                               coefficient_values); 


        cell_matrix_M = 0; 
        cell_matrix_E = 0; 
        for (unsigned int q = 0; q < n_q_points_dgrt; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
            { 
              const Tensor<1, dim> v_i = fe_values_dgrt[velocities].value(i, q); 
              for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k) 
                { 
                  const Tensor<1, dim> v_k = 
                    fe_values_dgrt[velocities].value(k, q); 

                  cell_matrix_E(i, k) += 
                    (coefficient_values[q] * v_i * v_k * fe_values_dgrt.JxW(q)); 

                  cell_matrix_M(i, k) += (v_i * v_k * fe_values_dgrt.JxW(q)); 
                } 
            } 


        cell_matrix_M.gauss_jordan(); 
        cell_matrix_M.mmult(cell_matrix_D, cell_matrix_E); 


        cell_matrix_G = 0; 
        for (unsigned int q = 0; q < n_q_points; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
            { 
              const double div_v_i = 
                fe_values_dgrt[velocities].divergence(i, q); 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                { 
                  const double phi_j_interior = 
                    fe_values[pressure_interior].value(j, q); 

                  cell_matrix_G(i, j) -= 
                    (div_v_i * phi_j_interior * fe_values.JxW(q)); 
                } 
            } 

        for (const auto &face : cell->face_iterators()) 
          { 
            fe_face_values.reinit(cell, face); 
            fe_face_values_dgrt.reinit(cell_dgrt, face); 

            for (unsigned int q = 0; q < n_face_q_points; ++q) 
              { 
                const Tensor<1, dim> &normal = fe_face_values.normal_vector(q); 

                for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
                  { 
                    const Tensor<1, dim> v_i = 
                      fe_face_values_dgrt[velocities].value(i, q); 
                    for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                      { 
                        const double phi_j_face = 
                          fe_face_values[pressure_face].value(j, q); 

                        cell_matrix_G(i, j) += 
                          ((v_i * normal) * phi_j_face * fe_face_values.JxW(q)); 
                      } 
                  } 
              } 
          } 
        cell_matrix_G.Tmmult(cell_matrix_C, cell_matrix_M); 


        cell->get_dof_values(solution, cell_solution); 


        cell_velocity = 0; 
        for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k) 
          for (unsigned int j = 0; j < dofs_per_cell_dgrt; ++j) 
            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              cell_velocity(k) += 
                -(cell_solution(i) * cell_matrix_C(i, j) * cell_matrix_D(k, j)); 


        cell_dgrt->get_dof_indices(local_dof_indices_dgrt); 
        for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k) 
          for (unsigned int j = 0; j < dofs_per_cell_dgrt; ++j) 
            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              darcy_velocity(local_dof_indices_dgrt[k]) += 
                -(cell_solution(i) * cell_matrix_C(i, j) * cell_matrix_D(k, j)); 
      } 
  } 



  template <int dim> 
  void WGDarcyEquation<dim>::compute_pressure_error() 
  { 
    Vector<float> difference_per_cell(triangulation.n_active_cells()); 
    const ComponentSelectFunction<dim> select_interior_pressure(0, 2); 
    VectorTools::integrate_difference(dof_handler, 
                                      solution, 
                                      ExactPressure<dim>(), 
                                      difference_per_cell, 
                                      QGauss<dim>(fe.degree + 2), 
                                      VectorTools::L2_norm, 
                                      &select_interior_pressure); 

    const double L2_error = difference_per_cell.l2_norm(); 
    std::cout << "L2_error_pressure " << L2_error << std::endl; 
  } 




  template <int dim> 
  void WGDarcyEquation<dim>::compute_velocity_errors() 
  { 
    const QGauss<dim>     quadrature_formula(fe_dgrt.degree + 1); 
    const QGauss<dim - 1> face_quadrature_formula(fe_dgrt.degree + 1); 

    FEValues<dim> fe_values_dgrt(fe_dgrt, 
                                 quadrature_formula, 
                                 update_values | update_gradients | 
                                   update_quadrature_points | 
                                   update_JxW_values); 

    FEFaceValues<dim> fe_face_values_dgrt(fe_dgrt, 
                                          face_quadrature_formula, 
                                          update_values | 
                                            update_normal_vectors | 
                                            update_quadrature_points | 
                                            update_JxW_values); 

    const unsigned int n_q_points_dgrt = fe_values_dgrt.get_quadrature().size(); 
    const unsigned int n_face_q_points_dgrt = 
      fe_face_values_dgrt.get_quadrature().size(); 

    std::vector<Tensor<1, dim>> velocity_values(n_q_points_dgrt); 
    std::vector<Tensor<1, dim>> velocity_face_values(n_face_q_points_dgrt); 

    const FEValuesExtractors::Vector velocities(0); 

    const ExactVelocity<dim> exact_velocity; 

    double L2_err_velocity_cell_sqr_global = 0; 
    double L2_err_flux_sqr                 = 0; 


    for (const auto &cell_dgrt : dof_handler_dgrt.active_cell_iterators()) 
      { 
        fe_values_dgrt.reinit(cell_dgrt); 


        fe_values_dgrt[velocities].get_function_values(darcy_velocity, 
                                                       velocity_values); 
        double L2_err_velocity_cell_sqr_local = 0; 
        for (unsigned int q = 0; q < n_q_points_dgrt; ++q) 
          { 
            const Tensor<1, dim> velocity = velocity_values[q]; 
            const Tensor<1, dim> true_velocity = 
              exact_velocity.value(fe_values_dgrt.quadrature_point(q)); 

            L2_err_velocity_cell_sqr_local += 
              ((velocity - true_velocity) * (velocity - true_velocity) * 
               fe_values_dgrt.JxW(q)); 
          } 
        L2_err_velocity_cell_sqr_global += L2_err_velocity_cell_sqr_local; 


        const double cell_area = cell_dgrt->measure(); 
        for (const auto &face_dgrt : cell_dgrt->face_iterators()) 
          { 
            const double face_length = face_dgrt->measure(); 
            fe_face_values_dgrt.reinit(cell_dgrt, face_dgrt); 
            fe_face_values_dgrt[velocities].get_function_values( 
              darcy_velocity, velocity_face_values); 

            double L2_err_flux_face_sqr_local = 0; 
            for (unsigned int q = 0; q < n_face_q_points_dgrt; ++q) 
              { 
                const Tensor<1, dim> velocity = velocity_face_values[q]; 
                const Tensor<1, dim> true_velocity = 
                  exact_velocity.value(fe_face_values_dgrt.quadrature_point(q)); 

                const Tensor<1, dim> &normal = 
                  fe_face_values_dgrt.normal_vector(q); 

 
                  ((velocity * normal - true_velocity * normal) * 
                   (velocity * normal - true_velocity * normal) * 
                   fe_face_values_dgrt.JxW(q)); 
              } 
            const double err_flux_each_face = 
              L2_err_flux_face_sqr_local / face_length * cell_area; 
            L2_err_flux_sqr += err_flux_each_face; 
          } 
      } 


    const double L2_err_velocity_cell = 
      std::sqrt(L2_err_velocity_cell_sqr_global); 
    const double L2_err_flux_face = std::sqrt(L2_err_flux_sqr); 

    std::cout << "L2_error_vel:  " << L2_err_velocity_cell << std::endl 
              << "L2_error_flux: " << L2_err_flux_face << std::endl; 
  } 




  template <int dim> 
  void WGDarcyEquation<dim>::output_results() const 
  { 
    { 
      DataOut<dim> data_out; 


      const std::vector<std::string> solution_names = {"interior_pressure", 
                                                       "interface_pressure"}; 
      data_out.add_data_vector(dof_handler, solution, solution_names); 


      const std::vector<std::string> velocity_names(dim, "velocity"); 
      const std::vector< 
        DataComponentInterpretation::DataComponentInterpretation> 
        velocity_component_interpretation( 
          dim, DataComponentInterpretation::component_is_part_of_vector); 
      data_out.add_data_vector(dof_handler_dgrt, 
                               darcy_velocity, 
                               velocity_names, 
                               velocity_component_interpretation); 

      data_out.build_patches(fe.degree); 
      std::ofstream output("solution_interior.vtu"); 
      data_out.write_vtu(output); 
    } 

    { 
      DataOutFaces<dim> data_out_faces(false); 
      data_out_faces.attach_dof_handler(dof_handler); 
      data_out_faces.add_data_vector(solution, "Pressure_Face"); 
      data_out_faces.build_patches(fe.degree); 
      std::ofstream face_output("solution_interface.vtu"); 
      data_out_faces.write_vtu(face_output); 
    } 
  } 


  template <int dim> 
  void WGDarcyEquation<dim>::run() 
  { 
    std::cout << "Solving problem in " << dim << " space dimensions." 
              << std::endl; 
    make_grid(); 
    setup_system(); 
    assemble_system(); 
    solve(); 
    compute_postprocessed_velocity(); 
    compute_pressure_error(); 
    compute_velocity_errors(); 
    output_results(); 
  } 

} // namespace Step61 


int main() 
{ 
  try 
    { 
      Step61::WGDarcyEquation<2> wg_darcy(0); 
      wg_darcy.run(); 
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

