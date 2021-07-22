

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2016 - 2021 by the deal.II authors 
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
 * Author: Timo Heister, Clemson University, 2016 
 */ 



#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/timer.h> 

// 下面这块出场代码与 step-40 相同，可以在PETSc和Trilinos之间切换。

#include <deal.II/lac/generic_linear_algebra.h> 

/* #define FORCE_USE_OF_TRILINOS */ 



namespace LA 
{ 
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \ 
  !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS)) 
  using namespace dealii::LinearAlgebraPETSc; 
#  define USE_PETSC_LA 
#elif defined(DEAL_II_WITH_TRILINOS) 
  using namespace dealii::LinearAlgebraTrilinos; 
#else 
#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required 
#endif 
} // namespace LA 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/solver_gmres.h> 
#include <deal.II/lac/solver_minres.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 

#include <deal.II/lac/petsc_sparse_matrix.h> 
#include <deal.II/lac/petsc_vector.h> 
#include <deal.II/lac/petsc_solver.h> 
#include <deal.II/lac/petsc_precondition.h> 

#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/manifold_lib.h> 
#include <deal.II/grid/grid_tools.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

#include <deal.II/base/utilities.h> 
#include <deal.II/base/conditional_ostream.h> 
#include <deal.II/base/index_set.h> 
#include <deal.II/lac/sparsity_tools.h> 
#include <deal.II/distributed/tria.h> 
#include <deal.II/distributed/grid_refinement.h> 

#include <cmath> 
#include <fstream> 
#include <iostream> 

namespace Step55 
{ 
  using namespace dealii; 
// @sect3{Linear solvers and preconditioners}  

// 我们需要一些辅助类来表示我们在介绍中描述的求解器策略。

  namespace LinearSolvers 
  { 

// 这个类暴露了通过函数 InverseMatrix::vmult(). 应用给定矩阵的逆的动作，在内部，逆不是显式形成的。相反，一个带有CG的线性求解器被执行。这个类扩展了 step-22 中的InverseMatrix类，增加了一个指定预处理程序的选项，并允许在vmult函数中使用不同的矢量类型。

    template <class Matrix, class Preconditioner> 
    class InverseMatrix : public Subscriptor 
    { 
    public: 
      InverseMatrix(const Matrix &m, const Preconditioner &preconditioner); 

      template <typename VectorType> 
      void vmult(VectorType &dst, const VectorType &src) const; 

    private: 
      const SmartPointer<const Matrix> matrix; 
      const Preconditioner &           preconditioner; 
    }; 

    template <class Matrix, class Preconditioner> 
    InverseMatrix<Matrix, Preconditioner>::InverseMatrix( 
      const Matrix &        m, 
      const Preconditioner &preconditioner) 
      : matrix(&m) 
      , preconditioner(preconditioner) 
    {} 

    template <class Matrix, class Preconditioner> 
    template <typename VectorType> 
    void 
    InverseMatrix<Matrix, Preconditioner>::vmult(VectorType &      dst, 
                                                 const VectorType &src) const 
    { 
      SolverControl solver_control(src.size(), 1e-8 * src.l2_norm()); 
      SolverCG<LA::MPI::Vector> cg(solver_control); 
      dst = 0; 

      try 
        { 
          cg.solve(*matrix, dst, src, preconditioner); 
        } 
      catch (std::exception &e) 
        { 
          Assert(false, ExcMessage(e.what())); 
        } 
    } 

// 该类是一个简单的2x2矩阵的块状对角线预处理器的模板类。

    template <class PreconditionerA, class PreconditionerS> 
    class BlockDiagonalPreconditioner : public Subscriptor 
    { 
    public: 
      BlockDiagonalPreconditioner(const PreconditionerA &preconditioner_A, 
                                  const PreconditionerS &preconditioner_S); 

      void vmult(LA::MPI::BlockVector &      dst, 
                 const LA::MPI::BlockVector &src) const; 

    private: 
      const PreconditionerA &preconditioner_A; 
      const PreconditionerS &preconditioner_S; 
    }; 

    template <class PreconditionerA, class PreconditionerS> 
    BlockDiagonalPreconditioner<PreconditionerA, PreconditionerS>:: 
      BlockDiagonalPreconditioner(const PreconditionerA &preconditioner_A, 
                                  const PreconditionerS &preconditioner_S) 
      : preconditioner_A(preconditioner_A) 
      , preconditioner_S(preconditioner_S) 
    {} 

    template <class PreconditionerA, class PreconditionerS> 
    void BlockDiagonalPreconditioner<PreconditionerA, PreconditionerS>::vmult( 
      LA::MPI::BlockVector &      dst, 
      const LA::MPI::BlockVector &src) const 
    { 
      preconditioner_A.vmult(dst.block(0), src.block(0)); 
      preconditioner_S.vmult(dst.block(1), src.block(1)); 
    } 

 
// @sect3{Problem setup}  

// 下面的类代表测试问题的右手边和精确解。

  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    RightHandSide() 
      : Function<dim>(dim + 1) 
    {} 

    virtual void vector_value(const Point<dim> &p, 
                              Vector<double> &  value) const override; 
  }; 

  template <int dim> 
  void RightHandSide<dim>::vector_value(const Point<dim> &p, 
                                        Vector<double> &  values) const 
  { 
    const double R_x = p[0]; 
    const double R_y = p[1]; 

    const double pi  = numbers::PI; 
    const double pi2 = pi * pi; 
    values[0] = 
      -1.0L / 2.0L * (-2 * sqrt(25.0 + 4 * pi2) + 10.0) * 
        exp(R_x * (-2 * sqrt(25.0 + 4 * pi2) + 10.0)) - 
      0.4 * pi2 * exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * cos(2 * R_y * pi) + 
      0.1 * pow(-sqrt(25.0 + 4 * pi2) + 5.0, 2) * 
        exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * cos(2 * R_y * pi); 
    values[1] = 0.2 * pi * (-sqrt(25.0 + 4 * pi2) + 5.0) * 
                  exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * sin(2 * R_y * pi) - 
                0.05 * pow(-sqrt(25.0 + 4 * pi2) + 5.0, 3) * 
                  exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * sin(2 * R_y * pi) / 
                  pi; 
    values[2] = 0; 
  } 

  template <int dim> 
  class ExactSolution : public Function<dim> 
  { 
  public: 
    ExactSolution() 
      : Function<dim>(dim + 1) 
    {} 

    virtual void vector_value(const Point<dim> &p, 
                              Vector<double> &  value) const override; 
  }; 

  template <int dim> 
  void ExactSolution<dim>::vector_value(const Point<dim> &p, 
                                        Vector<double> &  values) const 
  { 
    const double R_x = p[0]; 
    const double R_y = p[1]; 

    const double pi  = numbers::PI; 
    const double pi2 = pi * pi; 
    values[0] = 
      -exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * cos(2 * R_y * pi) + 1; 
    values[1] = (1.0L / 2.0L) * (-sqrt(25.0 + 4 * pi2) + 5.0) * 
                exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * sin(2 * R_y * pi) / 
                pi; 
    values[2] = 
      -1.0L / 2.0L * exp(R_x * (-2 * sqrt(25.0 + 4 * pi2) + 10.0)) - 
      2.0 * 
        (-6538034.74494422 + 
         0.0134758939981709 * exp(4 * sqrt(25.0 + 4 * pi2))) / 
        (-80.0 * exp(3 * sqrt(25.0 + 4 * pi2)) + 
         16.0 * sqrt(25.0 + 4 * pi2) * exp(3 * sqrt(25.0 + 4 * pi2))) - 
      1634508.68623606 * exp(-3.0 * sqrt(25.0 + 4 * pi2)) / 
        (-10.0 + 2.0 * sqrt(25.0 + 4 * pi2)) + 
      (-0.00673794699908547 * exp(sqrt(25.0 + 4 * pi2)) + 
       3269017.37247211 * exp(-3 * sqrt(25.0 + 4 * pi2))) / 
        (-8 * sqrt(25.0 + 4 * pi2) + 40.0) + 
      0.00336897349954273 * exp(1.0 * sqrt(25.0 + 4 * pi2)) / 
        (-10.0 + 2.0 * sqrt(25.0 + 4 * pi2)); 
  } 

//  @sect3{The main program}  

// 主类与  step-40  非常相似，只是矩阵和向量现在是块状的，而且我们为拥有的和相关的DoF存储一个  std::vector<IndexSet>  ，而不是一个IndexSet。我们正好有两个IndexSets，一个用于所有速度未知数，一个用于所有压力未知数。

  template <int dim> 
  class StokesProblem 
  { 
  public: 
    StokesProblem(unsigned int velocity_degree); 

    void run(); 

  private: 
    void make_grid(); 
    void setup_system(); 
    void assemble_system(); 
    void solve(); 
    void refine_grid(); 
    void output_results(const unsigned int cycle) const; 

    unsigned int velocity_degree; 
    double       viscosity; 
    MPI_Comm     mpi_communicator; 

    FESystem<dim>                             fe; 
    parallel::distributed::Triangulation<dim> triangulation; 
    DoFHandler<dim>                           dof_handler; 

    std::vector<IndexSet> owned_partitioning; 
    std::vector<IndexSet> relevant_partitioning; 

    AffineConstraints<double> constraints; 

    LA::MPI::BlockSparseMatrix system_matrix; 
    LA::MPI::BlockSparseMatrix preconditioner_matrix; 
    LA::MPI::BlockVector       locally_relevant_solution; 
    LA::MPI::BlockVector       system_rhs; 

    ConditionalOStream pcout; 
    TimerOutput        computing_timer; 
  }; 

  template <int dim> 
  StokesProblem<dim>::StokesProblem(unsigned int velocity_degree) 
    : velocity_degree(velocity_degree) 
    , viscosity(0.1) 
    , mpi_communicator(MPI_COMM_WORLD) 
    , fe(FE_Q<dim>(velocity_degree), dim, FE_Q<dim>(velocity_degree - 1), 1) 
    , triangulation(mpi_communicator, 
                    typename Triangulation<dim>::MeshSmoothing( 
                      Triangulation<dim>::smoothing_on_refinement | 
                      Triangulation<dim>::smoothing_on_coarsening)) 
    , dof_handler(triangulation) 
    , pcout(std::cout, 
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)) 
    , computing_timer(mpi_communicator, 
                      pcout, 
                      TimerOutput::summary, 
                      TimerOutput::wall_times) 
  {} 

// Kovasnay流定义在域[-0.5, 1.5]^2上，我们通过将最小和最大值传递给 GridGenerator::hyper_cube. 来创建这个域。
  template <int dim> 
  void StokesProblem<dim>::make_grid() 
  { 
    GridGenerator::hyper_cube(triangulation, -0.5, 1.5); 
    triangulation.refine_global(3); 
  } 
// @sect3{System Setup}  

// 与 step-40 相比，块矩阵和向量的构造是新的，与 step-22 这样的串行代码相比也是不同的，因为我们需要提供属于我们处理器的行的集合。

  template <int dim> 
  void StokesProblem<dim>::setup_system() 
  { 
    TimerOutput::Scope t(computing_timer, "setup"); 

    dof_handler.distribute_dofs(fe); 

// 将所有的昏暗速度放入0区块，压力放入1区块，然后按区块重新排列未知数。最后计算每块有多少个未知数。

    std::vector<unsigned int> stokes_sub_blocks(dim + 1, 0); 
    stokes_sub_blocks[dim] = 1; 
    DoFRenumbering::component_wise(dof_handler, stokes_sub_blocks); 

    const std::vector<types::global_dof_index> dofs_per_block = 
      DoFTools::count_dofs_per_fe_block(dof_handler, stokes_sub_blocks); 

    const unsigned int n_u = dofs_per_block[0]; 
    const unsigned int n_p = dofs_per_block[1]; 

    pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << " (" 
          << n_u << '+' << n_p << ')' << std::endl; 

// 我们根据我们想要创建块状矩阵和向量的方式，将本地拥有的和本地相关的DoF的IndexSet分割成两个IndexSets。

    owned_partitioning.resize(2); 
    owned_partitioning[0] = dof_handler.locally_owned_dofs().get_view(0, n_u); 
    owned_partitioning[1] = 
      dof_handler.locally_owned_dofs().get_view(n_u, n_u + n_p); 

    IndexSet locally_relevant_dofs; 
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs); 
    relevant_partitioning.resize(2); 
    relevant_partitioning[0] = locally_relevant_dofs.get_view(0, n_u); 
    relevant_partitioning[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p); 

// 设置边界条件和悬挂节点的约束与  step-40  相同。尽管我们没有任何悬空节点，因为我们只进行全局细化，但把这个函数调用放进去仍然是个好主意，以备以后引入自适应细化。

    { 
      constraints.reinit(locally_relevant_dofs); 

      FEValuesExtractors::Vector velocities(0); 
      DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
      VectorTools::interpolate_boundary_values(dof_handler, 
                                               0, 
                                               ExactSolution<dim>(), 
                                               constraints, 
                                               fe.component_mask(velocities)); 
      constraints.close(); 
    } 

// 现在我们根据BlockDynamicSparsityPattern来创建系统矩阵。我们知道我们不会有不同速度分量之间的耦合（因为我们使用的是拉普拉斯而不是变形张量），也不会有压力与其测试函数之间的耦合，所以我们使用一个表来将这个耦合信息传达给  DoFTools::make_sparsity_pattern.  。
    { 
      system_matrix.clear(); 

      Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1); 
      for (unsigned int c = 0; c < dim + 1; ++c) 
        for (unsigned int d = 0; d < dim + 1; ++d) 
          if (c == dim && d == dim) 
            coupling[c][d] = DoFTools::none; 
          else if (c == dim || d == dim || c == d) 
            coupling[c][d] = DoFTools::always; 
          else 
            coupling[c][d] = DoFTools::none; 

      BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block); 

      DoFTools::make_sparsity_pattern( 
        dof_handler, coupling, dsp, constraints, false); 

      SparsityTools::distribute_sparsity_pattern( 
        dsp, 
        dof_handler.locally_owned_dofs(), 
        mpi_communicator, 
        locally_relevant_dofs); 

      system_matrix.reinit(owned_partitioning, dsp, mpi_communicator); 
    } 

// 先决条件矩阵有不同的耦合（我们只在1,1块中填入质量矩阵），否则这段代码与上面的system_matrix的构造是相同的。

    { 
      preconditioner_matrix.clear(); 

      Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1); 
      for (unsigned int c = 0; c < dim + 1; ++c) 
        for (unsigned int d = 0; d < dim + 1; ++d) 
          if (c == dim && d == dim) 
            coupling[c][d] = DoFTools::always; 
          else 
            coupling[c][d] = DoFTools::none; 

      BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block); 

      DoFTools::make_sparsity_pattern( 
        dof_handler, coupling, dsp, constraints, false); 
      SparsityTools::distribute_sparsity_pattern( 
        dsp, 
        Utilities::MPI::all_gather(mpi_communicator, 
                                   dof_handler.locally_owned_dofs()), 
        mpi_communicator, 
        locally_relevant_dofs); 
      preconditioner_matrix.reinit(owned_partitioning, 

// owned_partitioning。

                                   dsp, 
                                   mpi_communicator); 
    } 

// 最后，我们以正确的尺寸构建块状向量。带有两个 std::vector<IndexSet> 的函数调用将创建一个重影向量。

    locally_relevant_solution.reinit(owned_partitioning, 
                                     relevant_partitioning, 
                                     mpi_communicator); 
    system_rhs.reinit(owned_partitioning, mpi_communicator); 
  } 

//  @sect3{Assembly}  

// 这个函数将系统矩阵、预处理矩阵和右手边集合起来。其代码非常标准。

  template <int dim> 
  void StokesProblem<dim>::assemble_system() 
  { 
    TimerOutput::Scope t(computing_timer, "assembly"); 

    system_matrix         = 0; 
    preconditioner_matrix = 0; 
    system_rhs            = 0; 

    const QGauss<dim> quadrature_formula(velocity_degree + 1); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    FullMatrix<double> cell_matrix2(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 

    const RightHandSide<dim>    right_hand_side; 
    std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim + 1)); 

    std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell); 
    std::vector<double>         div_phi_u(dofs_per_cell); 
    std::vector<double>         phi_p(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
    const FEValuesExtractors::Vector     velocities(0); 
    const FEValuesExtractors::Scalar     pressure(dim); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          cell_matrix  = 0; 
          cell_matrix2 = 0; 
          cell_rhs     = 0; 

          fe_values.reinit(cell); 
          right_hand_side.vector_value_list(fe_values.get_quadrature_points(), 
                                            rhs_values); 
          for (unsigned int q = 0; q < n_q_points; ++q) 
            { 
              for (unsigned int k = 0; k < dofs_per_cell; ++k) 
                { 
                  grad_phi_u[k] = fe_values[velocities].gradient(k, q); 
                  div_phi_u[k]  = fe_values[velocities].divergence(k, q); 
                  phi_p[k]      = fe_values[pressure].value(k, q); 
                } 

              for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                { 
                  for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                    { 
                      cell_matrix(i, j) += 
                        (viscosity * 
                           scalar_product(grad_phi_u[i], grad_phi_u[j]) - 
                         div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) * 
                        fe_values.JxW(q); 

                      cell_matrix2(i, j) += 1.0 / viscosity * phi_p[i] * 
                                            phi_p[j] * fe_values.JxW(q); 
                    } 

                  const unsigned int component_i = 
                    fe.system_to_component_index(i).first; 
                  cell_rhs(i) += fe_values.shape_value(i, q) * 
                                 rhs_values[q](component_i) * fe_values.JxW(q); 
                } 
            } 

          cell->get_dof_indices(local_dof_indices); 
          constraints.distribute_local_to_global(cell_matrix, 
                                                 cell_rhs, 
                                                 local_dof_indices, 
                                                 system_matrix, 
                                                 system_rhs); 

          constraints.distribute_local_to_global(cell_matrix2, 
                                                 local_dof_indices, 
                                                 preconditioner_matrix); 
        } 

    system_matrix.compress(VectorOperation::add); 
    preconditioner_matrix.compress(VectorOperation::add); 
    system_rhs.compress(VectorOperation::add); 
  } 

//  @sect3{Solving}  

// 这个函数用MINRES求解线性系统，如介绍中所述，对两个对角线块使用块状对角线预处理和AMG。预处理程序对0,0块应用v循环，对1,1块应用质量矩阵的CG（Schur补充）。

  template <int dim> 
  void StokesProblem<dim>::solve() 
  { 
    TimerOutput::Scope t(computing_timer, "solve"); 

    LA::MPI::PreconditionAMG prec_A; 
    { 
      LA::MPI::PreconditionAMG::AdditionalData data; 

#ifdef USE_PETSC_LA 
      data.symmetric_operator = true; 
#endif 
      prec_A.initialize(system_matrix.block(0, 0), data); 
    } 

    LA::MPI::PreconditionAMG prec_S; 
    { 
      LA::MPI::PreconditionAMG::AdditionalData data; 

#ifdef USE_PETSC_LA 
      data.symmetric_operator = true; 
#endif 
      prec_S.initialize(preconditioner_matrix.block(1, 1), data); 
    } 

// InverseMatrix用于解决质量矩阵的问题。

    using mp_inverse_t = LinearSolvers::InverseMatrix<LA::MPI::SparseMatrix, 
                                                      LA::MPI::PreconditionAMG>; 
    const mp_inverse_t mp_inverse(preconditioner_matrix.block(1, 1), prec_S); 

// 这是在上面定义的各个块的预处理的基础上构造的块预处理。

    const LinearSolvers::BlockDiagonalPreconditioner<LA::MPI::PreconditionAMG, 
                                                     mp_inverse_t> 
      preconditioner(prec_A, mp_inverse); 

// 有了这些，我们终于可以设置线性求解器并求解该系统。

    SolverControl solver_control(system_matrix.m(), 
                                 1e-10 * system_rhs.l2_norm()); 

    SolverMinRes<LA::MPI::BlockVector> solver(solver_control); 

    LA::MPI::BlockVector distributed_solution(owned_partitioning, 
                                              mpi_communicator); 

    constraints.set_zero(distributed_solution); 

    solver.solve(system_matrix, 
                 distributed_solution, 
                 system_rhs, 
                 preconditioner); 

    pcout << "   Solved in " << solver_control.last_step() << " iterations." 
          << std::endl; 

    constraints.distribute(distributed_solution); 

// 像在  step-56  中一样，我们减去平均压力，以便与我们的参考解决方案进行误差计算，该解决方案的平均值为零。

    locally_relevant_solution = distributed_solution; 
    const double mean_pressure = 
      VectorTools::compute_mean_value(dof_handler, 
                                      QGauss<dim>(velocity_degree + 2), 
                                      locally_relevant_solution, 
                                      dim); 
    distributed_solution.block(1).add(-mean_pressure); 
    locally_relevant_solution.block(1) = distributed_solution.block(1); 
  } 

//  @sect3{The rest}  

// 其余处理网格细化、输出和主循环的代码非常标准。

  template <int dim> 
  void StokesProblem<dim>::refine_grid() 
  { 
    TimerOutput::Scope t(computing_timer, "refine"); 

    triangulation.refine_global(); 
  } 

  template <int dim> 
  void StokesProblem<dim>::output_results(const unsigned int cycle) const 
  { 
    { 
      const ComponentSelectFunction<dim> pressure_mask(dim, dim + 1); 
      const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim), 
                                                       dim + 1); 

      Vector<double> cellwise_errors(triangulation.n_active_cells()); 
      QGauss<dim>    quadrature(velocity_degree + 2); 

      VectorTools::integrate_difference(dof_handler, 
                                        locally_relevant_solution, 
                                        ExactSolution<dim>(), 
                                        cellwise_errors, 
                                        quadrature, 
                                        VectorTools::L2_norm, 
                                        &velocity_mask); 

      const double error_u_l2 = 
        VectorTools::compute_global_error(triangulation, 
                                          cellwise_errors, 
                                          VectorTools::L2_norm); 

      VectorTools::integrate_difference(dof_handler, 
                                        locally_relevant_solution, 
                                        ExactSolution<dim>(), 
                                        cellwise_errors, 
                                        quadrature, 
                                        VectorTools::L2_norm, 
                                        &pressure_mask); 

      const double error_p_l2 = 
        VectorTools::compute_global_error(triangulation, 
                                          cellwise_errors, 
                                          VectorTools::L2_norm); 

      pcout << "error: u_0: " << error_u_l2 << " p_0: " << error_p_l2 
            << std::endl; 
    } 

    std::vector<std::string> solution_names(dim, "velocity"); 
    solution_names.emplace_back("pressure"); 
    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      data_component_interpretation( 
        dim, DataComponentInterpretation::component_is_part_of_vector); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 

    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(locally_relevant_solution, 
                             solution_names, 
                             DataOut<dim>::type_dof_data, 
                             data_component_interpretation); 

    LA::MPI::BlockVector interpolated; 
    interpolated.reinit(owned_partitioning, MPI_COMM_WORLD); 
    VectorTools::interpolate(dof_handler, ExactSolution<dim>(), interpolated); 

    LA::MPI::BlockVector interpolated_relevant(owned_partitioning, 
                                               relevant_partitioning, 
                                               MPI_COMM_WORLD); 
    interpolated_relevant = interpolated; 
    { 
      std::vector<std::string> solution_names(dim, "ref_u"); 
      solution_names.emplace_back("ref_p"); 
      data_out.add_data_vector(interpolated_relevant, 
                               solution_names, 
                               DataOut<dim>::type_dof_data, 
                               data_component_interpretation); 
    } 

    Vector<float> subdomain(triangulation.n_active_cells()); 
    for (unsigned int i = 0; i < subdomain.size(); ++i) 
      subdomain(i) = triangulation.locally_owned_subdomain(); 
    data_out.add_data_vector(subdomain, "subdomain"); 

    data_out.build_patches(); 

    data_out.write_vtu_with_pvtu_record( 
      "./", "solution", cycle, mpi_communicator, 2); 
  } 

  template <int dim> 
  void StokesProblem<dim>::run() 
  { 
#ifdef USE_PETSC_LA 
    pcout << "Running using PETSc." << std::endl; 
#else 
    pcout << "Running using Trilinos." << std::endl; 
#endif 
    const unsigned int n_cycles = 5; 
    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle) 
      { 
        pcout << "Cycle " << cycle << ':' << std::endl; 

        if (cycle == 0) 
          make_grid(); 
        else 
          refine_grid(); 

        setup_system(); 

        assemble_system(); 
        solve(); 

        if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32) 
          { 
            TimerOutput::Scope t(computing_timer, "output"); 
            output_results(cycle); 
          } 

        computing_timer.print_summary(); 
        computing_timer.reset(); 

        pcout << std::endl; 
      } 
  } 
} // namespace Step55 

int main(int argc, char *argv[]) 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step55; 

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 

      StokesProblem<2> problem(2); 
      problem.run(); 
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
