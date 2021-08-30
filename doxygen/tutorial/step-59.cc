

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
 * 
 * Authors: Katharina Kormann, Martin Kronbichler, 2018 
 */ 




#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/timer.h> 

#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/la_parallel_vector.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/tensor_product_matrix.h> 

#include <deal.II/fe/fe_dgq.h> 
#include <deal.II/fe/fe_tools.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_tools.h> 

#include <deal.II/multigrid/multigrid.h> 
#include <deal.II/multigrid/mg_transfer_matrix_free.h> 
#include <deal.II/multigrid/mg_tools.h> 
#include <deal.II/multigrid/mg_coarse.h> 
#include <deal.II/multigrid/mg_smoother.h> 
#include <deal.II/multigrid/mg_matrix.h> 

#include <deal.II/numerics/vector_tools.h> 

#include <deal.II/matrix_free/matrix_free.h> 
#include <deal.II/matrix_free/fe_evaluation.h> 

#include <iostream> 
#include <fstream> 

namespace Step59 
{ 
  using namespace dealii; 


  const unsigned int degree_finite_element = 8; 
  const unsigned int dimension             = 3; 


  template <int dim> 
  class Solution : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> &p, 
                         const unsigned int = 0) const override final 
    { 
      double val = 1.; 
      for (unsigned int d = 0; d < dim; ++d) 
        val *= std::cos(numbers::PI * 2.4 * p[d]); 
      return val; 
    } 

    virtual Tensor<1, dim> gradient(const Point<dim> &p, 
                                    const unsigned int = 0) const override final 
    { 
      const double   arg = numbers::PI * 2.4; 
      Tensor<1, dim> grad; 
      for (unsigned int d = 0; d < dim; ++d) 
        { 
          grad[d] = 1.; 
          for (unsigned int e = 0; e < dim; ++e) 
            if (d == e) 
              grad[d] *= -arg * std::sin(arg * p[e]); 
            else 
              grad[d] *= std::cos(arg * p[e]); 
        } 
      return grad; 
    } 
  }; 

  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> &p, 
                         const unsigned int = 0) const override final 
    { 
      const double arg = numbers::PI * 2.4; 
      double       val = 1.; 
      for (unsigned int d = 0; d < dim; ++d) 
        val *= std::cos(arg * p[d]); 
      return dim * arg * arg * val; 
    } 
  }; 




  template <int dim, int fe_degree, typename number> 
  class LaplaceOperator : public Subscriptor 
  { 
  public: 
    using value_type = number; 

    LaplaceOperator() = default; 

    void initialize(std::shared_ptr<const MatrixFree<dim, number>> data); 

    void clear(); 

    types::global_dof_index m() const; 

    void initialize_dof_vector( 
      LinearAlgebra::distributed::Vector<number> &vec) const; 

    std::shared_ptr<const MatrixFree<dim, number>> get_matrix_free() const; 

    void vmult(LinearAlgebra::distributed::Vector<number> &      dst, 
               const LinearAlgebra::distributed::Vector<number> &src) const; 

    void Tvmult(LinearAlgebra::distributed::Vector<number> &      dst, 
                const LinearAlgebra::distributed::Vector<number> &src) const; 

    number get_penalty_factor() const 
    { 
      return 1.0 * fe_degree * (fe_degree + 1); 
    } 

  private: 
    void 
    apply_cell(const MatrixFree<dim, number> &                   data, 
               LinearAlgebra::distributed::Vector<number> &      dst, 
               const LinearAlgebra::distributed::Vector<number> &src, 
               const std::pair<unsigned int, unsigned int> &cell_range) const; 

    void 
    apply_face(const MatrixFree<dim, number> &                   data, 
               LinearAlgebra::distributed::Vector<number> &      dst, 
               const LinearAlgebra::distributed::Vector<number> &src, 
               const std::pair<unsigned int, unsigned int> &face_range) const; 

    void apply_boundary( 
      const MatrixFree<dim, number> &                   data, 
      LinearAlgebra::distributed::Vector<number> &      dst, 
      const LinearAlgebra::distributed::Vector<number> &src, 
      const std::pair<unsigned int, unsigned int> &     face_range) const; 

    std::shared_ptr<const MatrixFree<dim, number>> data; 
  }; 


  template <int dim, int fe_degree, typename number> 
  class PreconditionBlockJacobi 
  { 
  public: 
    using value_type = number; 

    void clear() 
    { 
      cell_matrices.clear(); 
    } 

    void initialize(const LaplaceOperator<dim, fe_degree, number> &op); 

    void vmult(LinearAlgebra::distributed::Vector<number> &      dst, 
               const LinearAlgebra::distributed::Vector<number> &src) const; 

    void Tvmult(LinearAlgebra::distributed::Vector<number> &      dst, 
                const LinearAlgebra::distributed::Vector<number> &src) const 
    { 
      vmult(dst, src); 
    } 

  private: 
    std::shared_ptr<const MatrixFree<dim, number>> data; 
    std::vector<TensorProductMatrixSymmetricSum<dim, 
                                                VectorizedArray<number>, 
                                                fe_degree + 1>> 
      cell_matrices; 
  }; 


  template <int dim, typename number> 
  void adjust_ghost_range_if_necessary( 
    const MatrixFree<dim, number> &                   data, 
    const LinearAlgebra::distributed::Vector<number> &vec) 
  { 
    if (vec.get_partitioner().get() == 
        data.get_dof_info(0).vector_partitioner.get()) 
      return; 

    LinearAlgebra::distributed::Vector<number> copy_vec(vec); 
    const_cast<LinearAlgebra::distributed::Vector<number> &>(vec).reinit( 
      data.get_dof_info(0).vector_partitioner); 
    const_cast<LinearAlgebra::distributed::Vector<number> &>(vec) 
      .copy_locally_owned_data_from(copy_vec); 
  } 

  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::clear() 
  { 
    data.reset(); 
  } 

  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::initialize( 
    std::shared_ptr<const MatrixFree<dim, number>> data) 
  { 
    this->data = data; 
  } 

  template <int dim, int fe_degree, typename number> 
  std::shared_ptr<const MatrixFree<dim, number>> 
  LaplaceOperator<dim, fe_degree, number>::get_matrix_free() const 
  { 
    return data; 
  } 

  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::initialize_dof_vector( 
    LinearAlgebra::distributed::Vector<number> &vec) const 
  { 
    data->initialize_dof_vector(vec); 
  } 

  template <int dim, int fe_degree, typename number> 
  types::global_dof_index LaplaceOperator<dim, fe_degree, number>::m() const 
  { 
    Assert(data.get() != nullptr, ExcNotInitialized()); 
    return data->get_dof_handler().n_dofs(); 
  } 







  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::vmult( 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src) const 
  { 
    adjust_ghost_range_if_necessary(*data, dst); 
    adjust_ghost_range_if_necessary(*data, src); 
    data->loop(&LaplaceOperator::apply_cell, 
               &LaplaceOperator::apply_face, 
               &LaplaceOperator::apply_boundary, 
               this, 
               dst, 
               src, 
               /* zero_dst =  */ true,
               MatrixFree<dim, number>::DataAccessOnFaces::gradients,
               MatrixFree<dim, number>::DataAccessOnFaces::gradients);
  } 


  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::Tvmult( 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src) const 
  { 
    vmult(dst, src); 
  } 


  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::apply_cell( 
    const MatrixFree<dim, number> &                   data, 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src, 
    const std::pair<unsigned int, unsigned int> &     cell_range) const 
  { 
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(data); 
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) 
      { 
        phi.reinit(cell); 
        phi.gather_evaluate(src, EvaluationFlags::gradients); 
        for (unsigned int q = 0; q < phi.n_q_points; ++q) 
          phi.submit_gradient(phi.get_gradient(q), q); 
        phi.integrate_scatter(EvaluationFlags::gradients, dst); 
      } 
  } 



  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::apply_face( 
    const MatrixFree<dim, number> &                   data, 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src, 
    const std::pair<unsigned int, unsigned int> &     face_range) const 
  { 
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi_inner(data, 
                                                                         true); 
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi_outer(data, 
                                                                         false); 
    for (unsigned int face = face_range.first; face < face_range.second; ++face) 
      { 


        phi_inner.reinit(face); 
        phi_inner.gather_evaluate(src, 
                                  EvaluationFlags::values | 
                                    EvaluationFlags::gradients); 
        phi_outer.reinit(face); 
        phi_outer.gather_evaluate(src, 
                                  EvaluationFlags::values | 
                                    EvaluationFlags::gradients); 


        const VectorizedArray<number> inverse_length_normal_to_face = 
          0.5 * (std::abs((phi_inner.get_normal_vector(0) * 
                           phi_inner.inverse_jacobian(0))[dim - 1]) + 
                 std::abs((phi_outer.get_normal_vector(0) * 
                           phi_outer.inverse_jacobian(0))[dim - 1])); 
        const VectorizedArray<number> sigma = 
          inverse_length_normal_to_face * get_penalty_factor(); 


        for (unsigned int q = 0; q < phi_inner.n_q_points; ++q) 
          { 
            const VectorizedArray<number> solution_jump = 
              (phi_inner.get_value(q) - phi_outer.get_value(q)); 
            const VectorizedArray<number> average_normal_derivative = 
              (phi_inner.get_normal_derivative(q) + 
               phi_outer.get_normal_derivative(q)) * 
              number(0.5); 
            const VectorizedArray<number> test_by_value = 
              solution_jump * sigma - average_normal_derivative; 

            phi_inner.submit_value(test_by_value, q); 
            phi_outer.submit_value(-test_by_value, q); 

            phi_inner.submit_normal_derivative(-solution_jump * number(0.5), q); 
            phi_outer.submit_normal_derivative(-solution_jump * number(0.5), q); 
          } 


        phi_inner.integrate_scatter(EvaluationFlags::values | 
                                      EvaluationFlags::gradients, 
                                    dst); 
        phi_outer.integrate_scatter(EvaluationFlags::values | 
                                      EvaluationFlags::gradients, 
                                    dst); 
      } 
  } 



  template <int dim, int fe_degree, typename number> 
  void LaplaceOperator<dim, fe_degree, number>::apply_boundary( 
    const MatrixFree<dim, number> &                   data, 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src, 
    const std::pair<unsigned int, unsigned int> &     face_range) const 
  { 
    FEFaceEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi_inner(data, 
                                                                         true); 
    for (unsigned int face = face_range.first; face < face_range.second; ++face) 
      { 
        phi_inner.reinit(face); 
        phi_inner.gather_evaluate(src, 
                                  EvaluationFlags::values | 
                                    EvaluationFlags::gradients); 

        const VectorizedArray<number> inverse_length_normal_to_face = 
          std::abs((phi_inner.get_normal_vector(0) * 
                    phi_inner.inverse_jacobian(0))[dim - 1]); 
        const VectorizedArray<number> sigma = 
          inverse_length_normal_to_face * get_penalty_factor(); 

        const bool is_dirichlet = (data.get_boundary_id(face) == 0); 

        for (unsigned int q = 0; q < phi_inner.n_q_points; ++q) 
          { 
            const VectorizedArray<number> u_inner = phi_inner.get_value(q); 
            const VectorizedArray<number> u_outer = 
              is_dirichlet ? -u_inner : u_inner; 
            const VectorizedArray<number> normal_derivative_inner = 
              phi_inner.get_normal_derivative(q); 
            const VectorizedArray<number> normal_derivative_outer = 
              is_dirichlet ? normal_derivative_inner : -normal_derivative_inner; 
            const VectorizedArray<number> solution_jump = (u_inner - u_outer); 
            const VectorizedArray<number> average_normal_derivative = 
              (normal_derivative_inner + normal_derivative_outer) * number(0.5); 
            const VectorizedArray<number> test_by_value = 
              solution_jump * sigma - average_normal_derivative; 
            phi_inner.submit_normal_derivative(-solution_jump * number(0.5), q); 
            phi_inner.submit_value(test_by_value, q); 
          } 
        phi_inner.integrate_scatter(EvaluationFlags::values | 
                                      EvaluationFlags::gradients, 
                                    dst); 
      } 
  } 


  template <int dim, int fe_degree, typename number> 
  void PreconditionBlockJacobi<dim, fe_degree, number>::initialize( 
    const LaplaceOperator<dim, fe_degree, number> &op) 
  { 
    data = op.get_matrix_free(); 

    std::string name = data->get_dof_handler().get_fe().get_name(); 
    name.replace(name.find('<') + 1, 1, "1"); 
    std::unique_ptr<FiniteElement<1>> fe_1d = FETools::get_fe_by_name<1>(name); 


    const unsigned int                                 N = fe_degree + 1; 
    FullMatrix<double>                                 laplace_unscaled(N, N); 
    std::array<Table<2, VectorizedArray<number>>, dim> mass_matrices; 
    std::array<Table<2, VectorizedArray<number>>, dim> laplace_matrices; 
    for (unsigned int d = 0; d < dim; ++d) 
      { 
        mass_matrices[d].reinit(N, N); 
        laplace_matrices[d].reinit(N, N); 
      } 

    QGauss<1> quadrature(N); 
    for (unsigned int i = 0; i < N; ++i) 
      for (unsigned int j = 0; j < N; ++j) 
        { 
          double sum_mass = 0, sum_laplace = 0; 
          for (unsigned int q = 0; q < quadrature.size(); ++q) 
            { 
              sum_mass += (fe_1d->shape_value(i, quadrature.point(q)) * 
                           fe_1d->shape_value(j, quadrature.point(q))) * 
                          quadrature.weight(q); 
              sum_laplace += (fe_1d->shape_grad(i, quadrature.point(q))[0] * 
                              fe_1d->shape_grad(j, quadrature.point(q))[0]) * 
                             quadrature.weight(q); 
            } 
          for (unsigned int d = 0; d < dim; ++d) 
            mass_matrices[d](i, j) = sum_mass; 


          sum_laplace += 
            (1. * fe_1d->shape_value(i, Point<1>()) * 
               fe_1d->shape_value(j, Point<1>()) * op.get_penalty_factor() + 
             0.5 * fe_1d->shape_grad(i, Point<1>())[0] * 
               fe_1d->shape_value(j, Point<1>()) + 
             0.5 * fe_1d->shape_grad(j, Point<1>())[0] * 
               fe_1d->shape_value(i, Point<1>())); 

          sum_laplace += 
            (1. * fe_1d->shape_value(i, Point<1>(1.0)) * 
               fe_1d->shape_value(j, Point<1>(1.0)) * op.get_penalty_factor() - 
             0.5 * fe_1d->shape_grad(i, Point<1>(1.0))[0] * 
               fe_1d->shape_value(j, Point<1>(1.0)) - 
             0.5 * fe_1d->shape_grad(j, Point<1>(1.0))[0] * 
               fe_1d->shape_value(i, Point<1>(1.0))); 

          laplace_unscaled(i, j) = sum_laplace; 
        } 



    cell_matrices.clear(); 
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(*data); 
    unsigned int old_mapping_data_index = numbers::invalid_unsigned_int; 
    for (unsigned int cell = 0; cell < data->n_cell_batches(); ++cell) 
      { 
        phi.reinit(cell); 

        if (phi.get_mapping_data_index_offset() == old_mapping_data_index) 
          continue; 

        Tensor<2, dim, VectorizedArray<number>> inverse_jacobian = 
          phi.inverse_jacobian(0); 

        for (unsigned int d = 0; d < dim; ++d) 
          for (unsigned int e = 0; e < dim; ++e) 
            if (d != e) 
              for (unsigned int v = 0; v < VectorizedArray<number>::size(); ++v) 
                AssertThrow(inverse_jacobian[d][e][v] == 0., 
                            ExcNotImplemented()); 

        VectorizedArray<number> jacobian_determinant = inverse_jacobian[0][0]; 
        for (unsigned int e = 1; e < dim; ++e) 
          jacobian_determinant *= inverse_jacobian[e][e]; 
        jacobian_determinant = 1. / jacobian_determinant; 

        for (unsigned int d = 0; d < dim; ++d) 
          { 
            const VectorizedArray<number> scaling_factor = 
              inverse_jacobian[d][d] * inverse_jacobian[d][d] * 
              jacobian_determinant; 


            for (unsigned int i = 0; i < N; ++i) 
              for (unsigned int j = 0; j < N; ++j) 
                laplace_matrices[d](i, j) = 
                  scaling_factor * laplace_unscaled(i, j); 
          } 
        if (cell_matrices.size() <= phi.get_mapping_data_index_offset()) 
          cell_matrices.resize(phi.get_mapping_data_index_offset() + 1); 
        cell_matrices[phi.get_mapping_data_index_offset()].reinit( 
          mass_matrices, laplace_matrices); 
      } 
  } 


  template <int dim, int fe_degree, typename number> 
  void PreconditionBlockJacobi<dim, fe_degree, number>::vmult( 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src) const 
  { 
    adjust_ghost_range_if_necessary(*data, dst); 
    adjust_ghost_range_if_necessary(*data, src); 

    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number> phi(*data); 
    for (unsigned int cell = 0; cell < data->n_cell_batches(); ++cell) 
      { 
        phi.reinit(cell); 
        phi.read_dof_values(src); 
        cell_matrices[phi.get_mapping_data_index_offset()].apply_inverse( 
          ArrayView<VectorizedArray<number>>(phi.begin_dof_values(), 
                                             phi.dofs_per_cell), 
          ArrayView<const VectorizedArray<number>>(phi.begin_dof_values(), 
                                                   phi.dofs_per_cell)); 
        phi.set_dof_values(dst); 
      } 
  } 


  template <int dim, int fe_degree> 
  class LaplaceProblem 
  { 
  public: 
    LaplaceProblem(); 
    void run(); 

  private: 
    void setup_system(); 
    void compute_rhs(); 
    void solve(); 
    void analyze_results() const; 

#ifdef DEAL_II_WITH_P4EST 
    parallel::distributed::Triangulation<dim> triangulation; 
#else 
    Triangulation<dim> triangulation; 
#endif 

    FE_DGQHermite<dim> fe; 
    DoFHandler<dim>    dof_handler; 

    MappingQ1<dim> mapping; 

    using SystemMatrixType = LaplaceOperator<dim, fe_degree, double>; 
    SystemMatrixType system_matrix; 

    using LevelMatrixType = LaplaceOperator<dim, fe_degree, float>; 
    MGLevelObject<LevelMatrixType> mg_matrices; 

    LinearAlgebra::distributed::Vector<double> solution; 
    LinearAlgebra::distributed::Vector<double> system_rhs; 

    double             setup_time; 
    ConditionalOStream pcout; 
    ConditionalOStream time_details; 
  }; 

  template <int dim, int fe_degree> 
  LaplaceProblem<dim, fe_degree>::LaplaceProblem() 
    : 
#ifdef DEAL_II_WITH_P4EST 
    triangulation( 
      MPI_COMM_WORLD, 
      Triangulation<dim>::limit_level_difference_at_vertices, 
      parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy) 
    , 
#else 
    triangulation(Triangulation<dim>::limit_level_difference_at_vertices) 
    , 
#endif 
    fe(fe_degree) 
    , dof_handler(triangulation) 
    , setup_time(0.) 
    , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
    , time_details(std::cout, 
                   false && 
                     Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
  {} 


  template <int dim, int fe_degree> 
  void LaplaceProblem<dim, fe_degree>::setup_system() 
  { 
    Timer time; 
    setup_time = 0; 

    system_matrix.clear(); 
    mg_matrices.clear_elements(); 

    dof_handler.distribute_dofs(fe); 
    dof_handler.distribute_mg_dofs(); 

    pcout << "Number of degrees of freedom: " << dof_handler.n_dofs() 
          << std::endl; 

    setup_time += time.wall_time(); 
    time_details << "Distribute DoFs               " << time.wall_time() << " s" 
                 << std::endl; 
    time.restart(); 

    AffineConstraints<double> dummy; 
    dummy.close(); 

    { 
      typename MatrixFree<dim, double>::AdditionalData additional_data; 
      additional_data.tasks_parallel_scheme = 
        MatrixFree<dim, double>::AdditionalData::none; 
      additional_data.mapping_update_flags = 
        (update_gradients | update_JxW_values | update_quadrature_points); 
      additional_data.mapping_update_flags_inner_faces = 
        (update_gradients | update_JxW_values | update_normal_vectors); 
      additional_data.mapping_update_flags_boundary_faces = 
        (update_gradients | update_JxW_values | update_normal_vectors | 
         update_quadrature_points); 
      const auto system_mf_storage = 
        std::make_shared<MatrixFree<dim, double>>(); 
      system_mf_storage->reinit( 
        mapping, dof_handler, dummy, QGauss<1>(fe.degree + 1), additional_data); 
      system_matrix.initialize(system_mf_storage); 
    } 

    system_matrix.initialize_dof_vector(solution); 
    system_matrix.initialize_dof_vector(system_rhs); 

    setup_time += time.wall_time(); 
    time_details << "Setup matrix-free system      " << time.wall_time() << " s" 
                 << std::endl; 
    time.restart(); 

    const unsigned int nlevels = triangulation.n_global_levels(); 
    mg_matrices.resize(0, nlevels - 1); 

    for (unsigned int level = 0; level < nlevels; ++level) 
      { 
        typename MatrixFree<dim, float>::AdditionalData additional_data; 
        additional_data.tasks_parallel_scheme = 
          MatrixFree<dim, float>::AdditionalData::none; 
        additional_data.mapping_update_flags = 
          (update_gradients | update_JxW_values); 
        additional_data.mapping_update_flags_inner_faces = 
          (update_gradients | update_JxW_values); 
        additional_data.mapping_update_flags_boundary_faces = 
          (update_gradients | update_JxW_values); 
        additional_data.mg_level = level; 
        const auto mg_mf_storage_level = 
          std::make_shared<MatrixFree<dim, float>>(); 
        mg_mf_storage_level->reinit(mapping, 
                                    dof_handler, 
                                    dummy, 
                                    QGauss<1>(fe.degree + 1), 
                                    additional_data); 

        mg_matrices[level].initialize(mg_mf_storage_level); 
      } 
    setup_time += time.wall_time(); 
    time_details << "Setup matrix-free levels      " << time.wall_time() << " s" 
                 << std::endl; 
  } 


  template <int dim, int fe_degree> 
  void LaplaceProblem<dim, fe_degree>::compute_rhs() 
  { 
    Timer time; 
    system_rhs                          = 0; 
    const MatrixFree<dim, double> &data = *system_matrix.get_matrix_free(); 
    FEEvaluation<dim, fe_degree>   phi(data); 
    RightHandSide<dim>             rhs_func; 
    Solution<dim>                  exact_solution; 
    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell) 
      { 
        phi.reinit(cell); 
        for (unsigned int q = 0; q < phi.n_q_points; ++q) 
          { 
            VectorizedArray<double> rhs_val = VectorizedArray<double>(); 
            Point<dim, VectorizedArray<double>> point_batch = 
              phi.quadrature_point(q); 
            for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v) 
              { 
                Point<dim> single_point; 
                for (unsigned int d = 0; d < dim; ++d) 
                  single_point[d] = point_batch[d][v]; 
                rhs_val[v] = rhs_func.value(single_point); 
              } 
            phi.submit_value(rhs_val, q); 
          } 
        phi.integrate_scatter(EvaluationFlags::values, system_rhs); 
      } 


    FEFaceEvaluation<dim, fe_degree> phi_face(data, true); 
    for (unsigned int face = data.n_inner_face_batches(); 
         face < data.n_inner_face_batches() + data.n_boundary_face_batches(); 
         ++face) 
      { 
        phi_face.reinit(face); 

        const VectorizedArray<double> inverse_length_normal_to_face = 
          std::abs((phi_face.get_normal_vector(0) * 
                    phi_face.inverse_jacobian(0))[dim - 1]); 
        const VectorizedArray<double> sigma = 
          inverse_length_normal_to_face * system_matrix.get_penalty_factor(); 

        for (unsigned int q = 0; q < phi_face.n_q_points; ++q) 
          { 
            VectorizedArray<double> test_value = VectorizedArray<double>(), 
                                    test_normal_derivative = 
                                      VectorizedArray<double>(); 
            Point<dim, VectorizedArray<double>> point_batch = 
              phi_face.quadrature_point(q); 

            for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v) 
              { 
                Point<dim> single_point; 
                for (unsigned int d = 0; d < dim; ++d) 
                  single_point[d] = point_batch[d][v]; 


                if (data.get_boundary_id(face) == 0) 
                  test_value[v] = 2.0 * exact_solution.value(single_point); 
                else 
                  { 
                    Tensor<1, dim> normal; 
                    for (unsigned int d = 0; d < dim; ++d) 
                      normal[d] = phi_face.get_normal_vector(q)[d][v]; 
                    test_normal_derivative[v] = 
                      -normal * exact_solution.gradient(single_point); 
                  } 
              } 
            phi_face.submit_value(test_value * sigma - test_normal_derivative, 
                                  q); 
            phi_face.submit_normal_derivative(-0.5 * test_value, q); 
          } 
        phi_face.integrate_scatter(EvaluationFlags::values | 
                                     EvaluationFlags::gradients, 
                                   system_rhs); 
      } 


    system_rhs.compress(VectorOperation::add); 
    setup_time += time.wall_time(); 
    time_details << "Compute right hand side       " << time.wall_time() 
                 << " s\n"; 
  } 


  template <int dim, int fe_degree> 
  void LaplaceProblem<dim, fe_degree>::solve() 
  { 
    Timer                            time; 
    MGTransferMatrixFree<dim, float> mg_transfer; 
    mg_transfer.build(dof_handler); 
    setup_time += time.wall_time(); 
    time_details << "MG build transfer time        " << time.wall_time() 
                 << " s\n"; 
    time.restart(); 

    using SmootherType = 
      PreconditionChebyshev<LevelMatrixType, 
                            LinearAlgebra::distributed::Vector<float>, 
                            PreconditionBlockJacobi<dim, fe_degree, float>>; 
    mg::SmootherRelaxation<SmootherType, 
                           LinearAlgebra::distributed::Vector<float>> 
                                                         mg_smoother; 
    MGLevelObject<typename SmootherType::AdditionalData> smoother_data; 
    smoother_data.resize(0, triangulation.n_global_levels() - 1); 
    for (unsigned int level = 0; level < triangulation.n_global_levels(); 
         ++level) 
      { 
        if (level > 0) 
          { 
            smoother_data[level].smoothing_range     = 15.; 
            smoother_data[level].degree              = 3; 
            smoother_data[level].eig_cg_n_iterations = 10; 
          } 
        else 
          { 
            smoother_data[0].smoothing_range = 2e-2; 
            smoother_data[0].degree          = numbers::invalid_unsigned_int; 
            smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m(); 
          } 
        smoother_data[level].preconditioner = 
          std::make_shared<PreconditionBlockJacobi<dim, fe_degree, float>>(); 
        smoother_data[level].preconditioner->initialize(mg_matrices[level]); 
      } 
    mg_smoother.initialize(mg_matrices, smoother_data); 

    MGCoarseGridApplySmoother<LinearAlgebra::distributed::Vector<float>> 
      mg_coarse; 
    mg_coarse.initialize(mg_smoother); 

    mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_matrix( 
      mg_matrices); 

    Multigrid<LinearAlgebra::distributed::Vector<float>> mg( 
      mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother); 

    PreconditionMG<dim, 
                   LinearAlgebra::distributed::Vector<float>, 
                   MGTransferMatrixFree<dim, float>> 
      preconditioner(dof_handler, mg, mg_transfer); 

    SolverControl solver_control(10000, 1e-12 * system_rhs.l2_norm()); 
    SolverCG<LinearAlgebra::distributed::Vector<double>> cg(solver_control); 
    setup_time += time.wall_time(); 
    time_details << "MG build smoother time        " << time.wall_time() 
                 << "s\n"; 
    pcout << "Total setup time              " << setup_time << " s\n"; 

    time.reset(); 
    time.start(); 
    cg.solve(system_matrix, solution, system_rhs, preconditioner); 

    pcout << "Time solve (" << solver_control.last_step() << " iterations)    " 
          << time.wall_time() << " s" << std::endl; 
  } 


  template <int dim, int fe_degree> 
  void LaplaceProblem<dim, fe_degree>::analyze_results() const 
  { 
    Vector<float> error_per_cell(triangulation.n_active_cells()); 
    VectorTools::integrate_difference(mapping, 
                                      dof_handler, 
                                      solution, 
                                      Solution<dim>(), 
                                      error_per_cell, 
                                      QGauss<dim>(fe.degree + 2), 
                                      VectorTools::L2_norm); 
    pcout << "Verification via L2 error:    " 
          << std::sqrt( 
               Utilities::MPI::sum(error_per_cell.norm_sqr(), MPI_COMM_WORLD)) 
          << std::endl; 
  } 


  template <int dim, int fe_degree> 
  void LaplaceProblem<dim, fe_degree>::run() 
  { 
    const unsigned int n_ranks = 
      Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); 
    pcout << "Running with " << n_ranks << " MPI process" 
          << (n_ranks > 1 ? "es" : "") << ", element " << fe.get_name() 
          << std::endl 
          << std::endl; 
    for (unsigned int cycle = 0; cycle < 9 - dim; ++cycle) 
      { 
        pcout << "Cycle " << cycle << std::endl; 

        if (cycle == 0) 
          { 
            Point<dim> upper_right; 
            upper_right[0] = 2.5; 
            for (unsigned int d = 1; d < dim; ++d) 
              upper_right[d] = 2.8; 
            GridGenerator::hyper_rectangle(triangulation, 
                                           Point<dim>(), 
                                           upper_right); 
            triangulation.begin_active()->face(0)->set_boundary_id(10); 
            triangulation.begin_active()->face(1)->set_boundary_id(11); 
            triangulation.begin_active()->face(2)->set_boundary_id(0); 
            for (unsigned int f = 3; 
                 f < triangulation.begin_active()->n_faces(); 
                 ++f) 
              triangulation.begin_active()->face(f)->set_boundary_id(1); 

            std::vector<GridTools::PeriodicFacePair< 
              typename Triangulation<dim>::cell_iterator>> 
              periodic_faces; 
            GridTools::collect_periodic_faces( 
              triangulation, 10, 11, 0, periodic_faces); 
            triangulation.add_periodicity(periodic_faces); 

            triangulation.refine_global(6 - 2 * dim); 
          } 
        triangulation.refine_global(1); 
        setup_system(); 
        compute_rhs(); 
        solve(); 
        analyze_results(); 
        pcout << std::endl; 
      }; 
  } 
} // namespace Step59 


int main(int argc, char *argv[]) 
{ 
  try 
    { 
      using namespace Step59; 

      Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1); 

      LaplaceProblem<dimension, degree_finite_element> laplace_problem; 
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


