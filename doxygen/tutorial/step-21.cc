

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2006 - 2021 by the deal.II authors 
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
 * Authors: Yan Li, Wolfgang Bangerth, Texas A&M University, 2006 
 */ 






#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/function.h> 

#include <deal.II/lac/block_vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/block_sparse_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_tools.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_raviart_thomas.h> 
#include <deal.II/fe/fe_dgq.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 

#include <iostream> 
#include <fstream> 


#include <deal.II/base/tensor_function.h> 


#include <deal.II/base/discrete_time.h> 


namespace Step21 
{ 
  using namespace dealii; 






  template <int dim> 
  class TwoPhaseFlowProblem 
  { 
  public: 
    TwoPhaseFlowProblem(const unsigned int degree); 
    void run(); 

  private: 
    void   make_grid_and_dofs(); 
    void   assemble_system(); 
    void   assemble_rhs_S(); 
    double get_maximal_velocity() const; 
    void   solve(); 
    void   project_back_saturation(); 
    void   output_results() const; 

    const unsigned int degree; 

    Triangulation<dim> triangulation; 
    FESystem<dim>      fe; 
    DoFHandler<dim>    dof_handler; 

    BlockSparsityPattern      sparsity_pattern; 
    BlockSparseMatrix<double> system_matrix; 

    const unsigned int n_refinement_steps; 

    DiscreteTime time; 
    double       viscosity; 

    BlockVector<double> solution; 
    BlockVector<double> old_solution; 
    BlockVector<double> system_rhs; 
  }; 


  template <int dim> 
  class PressureRightHandSide : public Function<dim> 
  { 
  public: 
    PressureRightHandSide() 
      : Function<dim>(1) 
    {} 

    virtual double value(const Point<dim> & /*p*/, 
                         const unsigned int /*component*/ = 0) const override 
    { 
      return 0; 
    } 
  }; 



  template <int dim> 
  class PressureBoundaryValues : public Function<dim> 
  { 
  public: 
    PressureBoundaryValues() 
      : Function<dim>(1) 
    {} 

    virtual double value(const Point<dim> &p, 
                         const unsigned int /*component*/ = 0) const override 
    { 
      return 1 - p[0]; 
    } 
  }; 



  template <int dim> 
  class SaturationBoundaryValues : public Function<dim> 
  { 
  public: 
    SaturationBoundaryValues() 
      : Function<dim>(1) 
    {} 

    virtual double value(const Point<dim> &p, 
                         const unsigned int /*component*/ = 0) const override 
    { 
      if (p[0] == 0) 
        return 1; 
      else 
        return 0; 
    } 
  }; 



  template <int dim> 
  class InitialValues : public Function<dim> 
  { 
  public: 
    InitialValues() 
      : Function<dim>(dim + 2) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override 
    { 
      return Functions::ZeroFunction<dim>(dim + 2).value(p, component); 
    } 

    virtual void vector_value(const Point<dim> &p, 
                              Vector<double> &  values) const override 
    { 
      Functions::ZeroFunction<dim>(dim + 2).vector_value(p, values); 
    } 
  }; 





  namespace SingleCurvingCrack 
  { 
    template <int dim> 
    class KInverse : public TensorFunction<2, dim> 
    { 
    public: 
      KInverse() 
        : TensorFunction<2, dim>() 
      {} 

      virtual void 
      value_list(const std::vector<Point<dim>> &points, 
                 std::vector<Tensor<2, dim>> &  values) const override 
      { 
        Assert(points.size() == values.size(), 
               ExcDimensionMismatch(points.size(), values.size())); 

        for (unsigned int p = 0; p < points.size(); ++p) 
          { 
            values[p].clear(); 

            const double distance_to_flowline = 
              std::fabs(points[p][1] - 0.5 - 0.1 * std::sin(10 * points[p][0])); 

            const double permeability = 
              std::max(std::exp(-(distance_to_flowline * distance_to_flowline) / 
                                (0.1 * 0.1)), 
                       0.01); 

            for (unsigned int d = 0; d < dim; ++d) 
              values[p][d][d] = 1. / permeability; 
          } 
      } 
    }; 
  } // namespace SingleCurvingCrack 





  namespace RandomMedium 
  { 
    template <int dim> 
    class KInverse : public TensorFunction<2, dim> 
    { 
    public: 
      KInverse() 
        : TensorFunction<2, dim>() 
      {} 

      virtual void 
      value_list(const std::vector<Point<dim>> &points, 
                 std::vector<Tensor<2, dim>> &  values) const override 
      { 
        Assert(points.size() == values.size(), 
               ExcDimensionMismatch(points.size(), values.size())); 

        for (unsigned int p = 0; p < points.size(); ++p) 
          { 
            values[p].clear(); 

            double permeability = 0; 
            for (unsigned int i = 0; i < centers.size(); ++i) 
              permeability += std::exp(-(points[p] - centers[i]).norm_square() / 
                                       (0.05 * 0.05)); 

            const double normalized_permeability = 
              std::min(std::max(permeability, 0.01), 4.); 

            for (unsigned int d = 0; d < dim; ++d) 
              values[p][d][d] = 1. / normalized_permeability; 
          } 
      } 

    private: 
      static std::vector<Point<dim>> centers; 

      static std::vector<Point<dim>> get_centers() 
      { 
        const unsigned int N = 
          (dim == 2 ? 40 : (dim == 3 ? 100 : throw ExcNotImplemented())); 

        std::vector<Point<dim>> centers_list(N); 
        for (unsigned int i = 0; i < N; ++i) 
          for (unsigned int d = 0; d < dim; ++d) 
            centers_list[i][d] = static_cast<double>(rand()) / RAND_MAX; 

        return centers_list; 
      } 
    }; 

    template <int dim> 
    std::vector<Point<dim>> 
      KInverse<dim>::centers = KInverse<dim>::get_centers(); 
  } // namespace RandomMedium 



  double mobility_inverse(const double S, const double viscosity) 
  { 
    return 1.0 / (1.0 / viscosity * S * S + (1 - S) * (1 - S)); 
  } 

  double fractional_flow(const double S, const double viscosity) 
  { 
    return S * S / (S * S + viscosity * (1 - S) * (1 - S)); 
  } 



  template <class MatrixType> 
  class InverseMatrix : public Subscriptor 
  { 
  public: 
    InverseMatrix(const MatrixType &m) 
      : matrix(&m) 
    {} 

    void vmult(Vector<double> &dst, const Vector<double> &src) const 
    { 
      SolverControl solver_control(std::max<unsigned int>(src.size(), 200), 
                                   1e-8 * src.l2_norm()); 
      SolverCG<Vector<double>> cg(solver_control); 

      dst = 0; 

      cg.solve(*matrix, dst, src, PreconditionIdentity()); 
    } 

  private: 
    const SmartPointer<const MatrixType> matrix; 
  }; 

  class SchurComplement : public Subscriptor 
  { 
  public: 
    SchurComplement(const BlockSparseMatrix<double> &          A, 
                    const InverseMatrix<SparseMatrix<double>> &Minv) 
      : system_matrix(&A) 
      , m_inverse(&Minv) 
      , tmp1(A.block(0, 0).m()) 
      , tmp2(A.block(0, 0).m()) 
    {} 

    void vmult(Vector<double> &dst, const Vector<double> &src) const 
    { 
      system_matrix->block(0, 1).vmult(tmp1, src); 
      m_inverse->vmult(tmp2, tmp1); 
      system_matrix->block(1, 0).vmult(dst, tmp2); 
    } 

  private: 
    const SmartPointer<const BlockSparseMatrix<double>>           system_matrix; 
    const SmartPointer<const InverseMatrix<SparseMatrix<double>>> m_inverse; 

    mutable Vector<double> tmp1, tmp2; 
  }; 

  class ApproximateSchurComplement : public Subscriptor 
  { 
  public: 
    ApproximateSchurComplement(const BlockSparseMatrix<double> &A) 
      : system_matrix(&A) 
      , tmp1(A.block(0, 0).m()) 
      , tmp2(A.block(0, 0).m()) 
    {} 

    void vmult(Vector<double> &dst, const Vector<double> &src) const 
    { 
      system_matrix->block(0, 1).vmult(tmp1, src); 
      system_matrix->block(0, 0).precondition_Jacobi(tmp2, tmp1); 
      system_matrix->block(1, 0).vmult(dst, tmp2); 
    } 

  private: 
    const SmartPointer<const BlockSparseMatrix<double>> system_matrix; 

    mutable Vector<double> tmp1, tmp2; 
  }; 





  template <int dim> 
  TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem(const unsigned int degree) 
    : degree(degree) 
    , fe(FE_RaviartThomas<dim>(degree), 
         1, 
         FE_DGQ<dim>(degree), 
         1, 
         FE_DGQ<dim>(degree), 
         1) 
    , dof_handler(triangulation) 
    , n_refinement_steps(5) 
    , time(/*start time*/ 0., /*end time*/ 1.) 
    , viscosity(0.2) 
  {} 



  template <int dim> 
  void TwoPhaseFlowProblem<dim>::make_grid_and_dofs() 
  { 
    GridGenerator::hyper_cube(triangulation, 0, 1); 
    triangulation.refine_global(n_refinement_steps); 

    dof_handler.distribute_dofs(fe); 
    DoFRenumbering::component_wise(dof_handler); 

    const std::vector<types::global_dof_index> dofs_per_component = 
      DoFTools::count_dofs_per_fe_component(dof_handler); 
    const unsigned int n_u = dofs_per_component[0], 
                       n_p = dofs_per_component[dim], 
                       n_s = dofs_per_component[dim + 1]; 

    std::cout << "Number of active cells: " << triangulation.n_active_cells() 
              << std::endl 
              << "Number of degrees of freedom: " << dof_handler.n_dofs() 
              << " (" << n_u << '+' << n_p << '+' << n_s << ')' << std::endl 
              << std::endl; 

    const unsigned int n_couplings = dof_handler.max_couplings_between_dofs(); 

    sparsity_pattern.reinit(3, 3); 
    sparsity_pattern.block(0, 0).reinit(n_u, n_u, n_couplings); 
    sparsity_pattern.block(1, 0).reinit(n_p, n_u, n_couplings); 
    sparsity_pattern.block(2, 0).reinit(n_s, n_u, n_couplings); 
    sparsity_pattern.block(0, 1).reinit(n_u, n_p, n_couplings); 
    sparsity_pattern.block(1, 1).reinit(n_p, n_p, n_couplings); 
    sparsity_pattern.block(2, 1).reinit(n_s, n_p, n_couplings); 
    sparsity_pattern.block(0, 2).reinit(n_u, n_s, n_couplings); 
    sparsity_pattern.block(1, 2).reinit(n_p, n_s, n_couplings); 
    sparsity_pattern.block(2, 2).reinit(n_s, n_s, n_couplings); 

    sparsity_pattern.collect_sizes(); 

    DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern); 
    sparsity_pattern.compress(); 

    system_matrix.reinit(sparsity_pattern); 

    solution.reinit(3); 
    solution.block(0).reinit(n_u); 
    solution.block(1).reinit(n_p); 
    solution.block(2).reinit(n_s); 
    solution.collect_sizes(); 

    old_solution.reinit(3); 
    old_solution.block(0).reinit(n_u); 
    old_solution.block(1).reinit(n_p); 
    old_solution.block(2).reinit(n_s); 
    old_solution.collect_sizes(); 

    system_rhs.reinit(3); 
    system_rhs.block(0).reinit(n_u); 
    system_rhs.block(1).reinit(n_p); 
    system_rhs.block(2).reinit(n_s); 
    system_rhs.collect_sizes(); 
  } 



  template <int dim> 
  void TwoPhaseFlowProblem<dim>::assemble_system() 
  { 
    system_matrix = 0; 
    system_rhs    = 0; 

    QGauss<dim>     quadrature_formula(degree + 2); 
    QGauss<dim - 1> face_quadrature_formula(degree + 2); 

    FEValues<dim>     fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 
    FEFaceValues<dim> fe_face_values(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_normal_vectors | 
                                       update_quadrature_points | 
                                       update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 

    const unsigned int n_q_points      = quadrature_formula.size(); 
    const unsigned int n_face_q_points = face_quadrature_formula.size(); 

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     local_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    const PressureRightHandSide<dim>  pressure_right_hand_side; 
    const PressureBoundaryValues<dim> pressure_boundary_values; 
    const RandomMedium::KInverse<dim> k_inverse; 

    std::vector<double>         pressure_rhs_values(n_q_points); 
    std::vector<double>         boundary_values(n_face_q_points); 
    std::vector<Tensor<2, dim>> k_inverse_values(n_q_points); 

    std::vector<Vector<double>>              old_solution_values(n_q_points, 
                                                                 Vector<double>(dim + 2)); 
    std::vector<std::vector<Tensor<1, dim>>> old_solution_grads( 
      n_q_points, std::vector<Tensor<1, dim>>(dim + 2)); 

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure(dim); 
    const FEValuesExtractors::Scalar saturation(dim + 1); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_values.reinit(cell); 
        local_matrix = 0; 
        local_rhs    = 0; 


        fe_values.get_function_values(old_solution, old_solution_values); 


        pressure_right_hand_side.value_list(fe_values.get_quadrature_points(), 
                                            pressure_rhs_values); 
        k_inverse.value_list(fe_values.get_quadrature_points(), 
                             k_inverse_values); 


        for (unsigned int q = 0; q < n_q_points; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              const double old_s = old_solution_values[q](dim + 1); 

              const Tensor<1, dim> phi_i_u = fe_values[velocities].value(i, q); 
              const double div_phi_i_u = fe_values[velocities].divergence(i, q); 
              const double phi_i_p     = fe_values[pressure].value(i, q); 
              const double phi_i_s     = fe_values[saturation].value(i, q); 

              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                { 
                  const Tensor<1, dim> phi_j_u = 
                    fe_values[velocities].value(j, q); 
                  const double div_phi_j_u = 
                    fe_values[velocities].divergence(j, q); 
                  const double phi_j_p = fe_values[pressure].value(j, q); 
                  const double phi_j_s = fe_values[saturation].value(j, q); 

                  local_matrix(i, j) += 
                    (phi_i_u * k_inverse_values[q] * 
                       mobility_inverse(old_s, viscosity) * phi_j_u - 
                     div_phi_i_u * phi_j_p - phi_i_p * div_phi_j_u + 
                     phi_i_s * phi_j_s) * 
                    fe_values.JxW(q); 
                } 

              local_rhs(i) += 
                (-phi_i_p * pressure_rhs_values[q]) * fe_values.JxW(q); 
            } 


        for (const auto &face : cell->face_iterators()) 
          if (face->at_boundary()) 
            { 
              fe_face_values.reinit(cell, face); 

              pressure_boundary_values.value_list( 
                fe_face_values.get_quadrature_points(), boundary_values); 

              for (unsigned int q = 0; q < n_face_q_points; ++q) 
                for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                  { 
                    const Tensor<1, dim> phi_i_u = 
                      fe_face_values[velocities].value(i, q); 

                    local_rhs(i) += 
                      -(phi_i_u * fe_face_values.normal_vector(q) * 
                        boundary_values[q] * fe_face_values.JxW(q)); 
                  } 
            } 


        cell->get_dof_indices(local_dof_indices); 
        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
            system_matrix.add(local_dof_indices[i], 
                              local_dof_indices[j], 
                              local_matrix(i, j)); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          system_rhs(local_dof_indices[i]) += local_rhs(i); 
      } 
  } 




  template <int dim> 
  void TwoPhaseFlowProblem<dim>::assemble_rhs_S() 
  { 
    QGauss<dim>       quadrature_formula(degree + 2); 
    QGauss<dim - 1>   face_quadrature_formula(degree + 2); 
    FEValues<dim>     fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 
    FEFaceValues<dim> fe_face_values(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_normal_vectors | 
                                       update_quadrature_points | 
                                       update_JxW_values); 
    FEFaceValues<dim> fe_face_values_neighbor(fe, 
                                              face_quadrature_formula, 
                                              update_values); 

    const unsigned int dofs_per_cell   = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points      = quadrature_formula.size(); 
    const unsigned int n_face_q_points = face_quadrature_formula.size(); 

    Vector<double> local_rhs(dofs_per_cell); 

    std::vector<Vector<double>> old_solution_values(n_q_points, 
                                                    Vector<double>(dim + 2)); 
    std::vector<Vector<double>> old_solution_values_face(n_face_q_points, 
                                                         Vector<double>(dim + 
                                                                        2)); 
    std::vector<Vector<double>> old_solution_values_face_neighbor( 
      n_face_q_points, Vector<double>(dim + 2)); 
    std::vector<Vector<double>> present_solution_values(n_q_points, 
                                                        Vector<double>(dim + 
                                                                       2)); 
    std::vector<Vector<double>> present_solution_values_face( 
      n_face_q_points, Vector<double>(dim + 2)); 

    std::vector<double>                  neighbor_saturation(n_face_q_points); 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    SaturationBoundaryValues<dim> saturation_boundary_values; 

    const FEValuesExtractors::Scalar saturation(dim + 1); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        local_rhs = 0; 
        fe_values.reinit(cell); 

        fe_values.get_function_values(old_solution, old_solution_values); 
        fe_values.get_function_values(solution, present_solution_values); 


        for (unsigned int q = 0; q < n_q_points; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              const double   old_s = old_solution_values[q](dim + 1); 
              Tensor<1, dim> present_u; 
              for (unsigned int d = 0; d < dim; ++d) 
                present_u[d] = present_solution_values[q](d); 

              const double         phi_i_s = fe_values[saturation].value(i, q); 
              const Tensor<1, dim> grad_phi_i_s = 
                fe_values[saturation].gradient(i, q); 

              local_rhs(i) += 
                (time.get_next_step_size() * fractional_flow(old_s, viscosity) * 
                   present_u * grad_phi_i_s + 
                 old_s * phi_i_s) * 
                fe_values.JxW(q); 
            } 



        for (const auto face_no : cell->face_indices()) 
          { 
            fe_face_values.reinit(cell, face_no); 

            fe_face_values.get_function_values(old_solution, 
                                               old_solution_values_face); 
            fe_face_values.get_function_values(solution, 
                                               present_solution_values_face); 

            if (cell->at_boundary(face_no)) 
              saturation_boundary_values.value_list( 
                fe_face_values.get_quadrature_points(), neighbor_saturation); 
            else 
              { 
                const auto         neighbor = cell->neighbor(face_no); 
                const unsigned int neighbor_face = 
                  cell->neighbor_of_neighbor(face_no); 

                fe_face_values_neighbor.reinit(neighbor, neighbor_face); 

                fe_face_values_neighbor.get_function_values( 
                  old_solution, old_solution_values_face_neighbor); 

                for (unsigned int q = 0; q < n_face_q_points; ++q) 
                  neighbor_saturation[q] = 
                    old_solution_values_face_neighbor[q](dim + 1); 
              } 

            for (unsigned int q = 0; q < n_face_q_points; ++q) 
              { 
                Tensor<1, dim> present_u_face; 
                for (unsigned int d = 0; d < dim; ++d) 
                  present_u_face[d] = present_solution_values_face[q](d); 

                const double normal_flux = 
                  present_u_face * fe_face_values.normal_vector(q); 

                const bool is_outflow_q_point = (normal_flux >= 0); 

                for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                  local_rhs(i) -= 
                    time.get_next_step_size() * normal_flux * 
                    fractional_flow((is_outflow_q_point == true ? 
                                       old_solution_values_face[q](dim + 1) : 
                                       neighbor_saturation[q]), 
                                    viscosity) * 
                    fe_face_values[saturation].value(i, q) * 
                    fe_face_values.JxW(q); 
              } 
          } 

        cell->get_dof_indices(local_dof_indices); 
        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          system_rhs(local_dof_indices[i]) += local_rhs(i); 
      } 
  } 



  template <int dim> 
  void TwoPhaseFlowProblem<dim>::solve() 
  { 
    const InverseMatrix<SparseMatrix<double>> m_inverse( 
      system_matrix.block(0, 0)); 
    Vector<double> tmp(solution.block(0).size()); 
    Vector<double> schur_rhs(solution.block(1).size()); 
    Vector<double> tmp2(solution.block(2).size()); 


    { 
      m_inverse.vmult(tmp, system_rhs.block(0)); 
      system_matrix.block(1, 0).vmult(schur_rhs, tmp); 
      schur_rhs -= system_rhs.block(1); 

      SchurComplement schur_complement(system_matrix, m_inverse); 

      ApproximateSchurComplement approximate_schur_complement(system_matrix); 

      InverseMatrix<ApproximateSchurComplement> preconditioner( 
        approximate_schur_complement); 

      SolverControl            solver_control(solution.block(1).size(), 
                                   1e-12 * schur_rhs.l2_norm()); 
      SolverCG<Vector<double>> cg(solver_control); 

      cg.solve(schur_complement, solution.block(1), schur_rhs, preconditioner); 

      std::cout << "   " << solver_control.last_step() 
                << " CG Schur complement iterations for pressure." << std::endl; 
    } 


    { 
      system_matrix.block(0, 1).vmult(tmp, solution.block(1)); 
      tmp *= -1; 
      tmp += system_rhs.block(0); 

      m_inverse.vmult(solution.block(0), tmp); 
    } 



    time.set_desired_next_step_size(std::pow(0.5, double(n_refinement_steps)) / 
                                    get_maximal_velocity()); 


    assemble_rhs_S(); 
    { 
      SolverControl            solver_control(system_matrix.block(2, 2).m(), 
                                   1e-8 * system_rhs.block(2).l2_norm()); 
      SolverCG<Vector<double>> cg(solver_control); 
      cg.solve(system_matrix.block(2, 2), 
               solution.block(2), 
               system_rhs.block(2), 
               PreconditionIdentity()); 

      project_back_saturation(); 

      std::cout << "   " << solver_control.last_step() 
                << " CG iterations for saturation." << std::endl; 
    } 

    old_solution = solution; 
  } 



  template <int dim> 
  void TwoPhaseFlowProblem<dim>::output_results() const 
  { 
    if (time.get_step_number() % 5 != 0) 
      return; 

    std::vector<std::string> solution_names; 
    switch (dim) 
      { 
        case 2: 
          solution_names = {"u", "v", "p", "S"}; 
          break; 

        case 3: 
          solution_names = {"u", "v", "w", "p", "S"}; 
          break; 

        default: 
          Assert(false, ExcNotImplemented()); 
      } 

    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, solution_names); 

    data_out.build_patches(degree + 1); 

    std::ofstream output("solution-" + 
                         Utilities::int_to_string(time.get_step_number(), 4) + 
                         ".vtk"); 
    data_out.write_vtk(output); 
  } 




  template <int dim> 
  void TwoPhaseFlowProblem<dim>::project_back_saturation() 
  { 
    for (unsigned int i = 0; i < solution.block(2).size(); ++i) 
      if (solution.block(2)(i) < 0) 
        solution.block(2)(i) = 0; 
      else if (solution.block(2)(i) > 1) 
        solution.block(2)(i) = 1; 
  } 


  template <int dim> 
  double TwoPhaseFlowProblem<dim>::get_maximal_velocity() const 
  { 
    QGauss<dim>        quadrature_formula(degree + 2); 
    const unsigned int n_q_points = quadrature_formula.size(); 

    FEValues<dim> fe_values(fe, quadrature_formula, update_values); 
    std::vector<Vector<double>> solution_values(n_q_points, 
                                                Vector<double>(dim + 2)); 
    double                      max_velocity = 0; 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_values.reinit(cell); 
        fe_values.get_function_values(solution, solution_values); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            Tensor<1, dim> velocity; 
            for (unsigned int i = 0; i < dim; ++i) 
              velocity[i] = solution_values[q](i); 

            max_velocity = std::max(max_velocity, velocity.norm()); 
          } 
      } 

    return max_velocity; 
  } 



  template <int dim> 
  void TwoPhaseFlowProblem<dim>::run() 
  { 
    make_grid_and_dofs(); 

    { 
      AffineConstraints<double> constraints; 
      constraints.close(); 

      VectorTools::project(dof_handler, 
                           constraints, 
                           QGauss<dim>(degree + 2), 
                           InitialValues<dim>(), 
                           old_solution); 
    } 

    do 
      { 
        std::cout << "Timestep " << time.get_step_number() + 1 << std::endl; 

        assemble_system(); 

        solve(); 

        output_results(); 

        time.advance_time(); 
        std::cout << "   Now at t=" << time.get_current_time() 
                  << ", dt=" << time.get_previous_step_size() << '.' 
                  << std::endl 
                  << std::endl; 
      } 
    while (time.is_at_end() == false); 
  } 
} // namespace Step21 


int main() 
{ 
  try 
    { 
      using namespace Step21; 

      TwoPhaseFlowProblem<2> two_phase_flow_problem(0); 
      two_phase_flow_problem.run(); 
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


