

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
 * Author: Marc Fehling, Colorado State University, 2021
 *         Peter Munch, Technical University of Munich and Helmholtz-Zentrum
 *                      hereon, 2021
 *         Wolfgang Bangerth, Colorado State University, 2021
 */


// @sect3{Include files}

// 在以前的教程程序中，特别是在 step-27 和 step-40
// 中，已经使用和讨论了以下包含文件。

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_series.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/refinement.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/smoothness_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <fstream>
#include <iostream>

// 为了实现负载平衡，我们将在单元格上分配单独的权重，为此我们将使用类
// parallel::CellWeights.  。
#include <deal.II/distributed/cell_weights.h>

// 求解函数需要从直角坐标到极坐标的转换。 GeometricUtilities::Coordinates
// 命名空间提供了必要的工具。

#include <deal.II/base/function.h>
#include <deal.II/base/geometric_utilities.h>

// 以下包含的文件将启用MatrixFree功能。

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tools.h>

// 我们将使用 LinearAlgebra::distributed::Vector 进行线性代数操作。

#include <deal.II/lac/la_parallel_vector.h>

// 我们剩下的就是包含多网格求解器所需的文件。

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/multigrid.h>

namespace Step75
{
  using namespace dealii;
  // @sect3{The <code>Solution</code> class template}

  // 我们有一个分析性的方案可以使用。我们将用这个解来为问题的数值解施加边界条件。解决方案的表述需要转换为极坐标。为了从笛卡尔坐标转换到球面坐标，我们将使用
  // GeometricUtilities::Coordinates
  // 命名空间的一个辅助函数。这个转换的前两个坐标对应于x-y面的极坐标。

  template <int dim>
  class Solution : public Function<dim>
  {
  public:
    Solution()
      : Function<dim>()
    {}

    virtual double
    value(const Point<dim> &p, const unsigned int /*component*/) const override
    {
      const std::array<double, dim> p_sphere =
        GeometricUtilities::Coordinates::to_spherical(p);

      constexpr const double alpha = 2. / 3.;
      return std::pow(p_sphere[0], alpha) * std::sin(alpha * p_sphere[1]);
    }
  };

  //  @sect3{Parameters}

  // 在本教程中，我们将使用一个简化的参数集。这里也可以使用ParameterHandler类，但为了使本教程简短，我们决定使用简单的结构。所有这些参数的实际意图将在接下来的类中描述，在它们各自使用的位置。

  // 下面的参数集控制着多网格机制的粗网格求解器、平滑器和网格间传输方案。我们用默认参数来填充它。

  struct MultigridParameters
  {
    struct
    {
      std::string  type            = "cg_with_amg";
      unsigned int maxiter         = 10000;
      double       abstol          = 1e-20;
      double       reltol          = 1e-4;
      unsigned int smoother_sweeps = 1;
      unsigned int n_cycles        = 1;
      std::string  smoother_type   = "ILU";
    } coarse_solver;

    struct
    {
      std::string  type                = "chebyshev";
      double       smoothing_range     = 20;
      unsigned int degree              = 5;
      unsigned int eig_cg_n_iterations = 20;
    } smoother;

    struct
    {
      MGTransferGlobalCoarseningTools::PolynomialCoarseningSequenceType
        p_sequence = MGTransferGlobalCoarseningTools::
          PolynomialCoarseningSequenceType::decrease_by_one;
      bool perform_h_transfer = true;
    } transfer;
  };

  // 这是该问题类的一般参数结构。你会发现这个结构分为几个类别，包括一般的运行时参数、级别限制、细化和粗化分数，以及单元加权的参数。它还包含一个上述结构的实例，用于多网格参数，这些参数将被传递给多网格算法。

  struct Parameters
  {
    unsigned int n_cycles         = 8;
    double       tolerance_factor = 1e-12;

    MultigridParameters mg_data;

    unsigned int min_h_level            = 5;
    unsigned int max_h_level            = 12;
    unsigned int min_p_degree           = 2;
    unsigned int max_p_degree           = 6;
    unsigned int max_p_level_difference = 1;

    double refine_fraction    = 0.3;
    double coarsen_fraction   = 0.03;
    double p_refine_fraction  = 0.9;
    double p_coarsen_fraction = 0.9;

    double weighting_factor   = 1e6;
    double weighting_exponent = 1.;
  };

  //  @sect3{Matrix-free Laplace operator}

  // 这是一个无矩阵的拉普拉斯算子的实现，基本上将接管其他教程中的`assemble_system()`函数的部分。所有成员函数的含义将在后面的定义中解释。

  // 我们将使用FEEvaluation类来评估正交点的解向量并进行积分。与其他教程不同的是，模板参数`度数`被设置为
  // $-1$  ，`一维正交数`被设置为  $0$
  // 。在这种情况下，FEEvaluation会动态地选择正确的多项式度数和正交点的数量。在这里，我们为FEEvaluation引入一个带有正确模板参数的别名，这样我们以后就不用担心这些参数了。

  template <int dim, typename number>
  class LaplaceOperator : public Subscriptor
  {
  public:
    using VectorType = LinearAlgebra::distributed::Vector<number>;

    using FECellIntegrator = FEEvaluation<dim, -1, 0, 1, number>;

    LaplaceOperator() = default;

    LaplaceOperator(const hp::MappingCollection<dim> &mapping,
                    const DoFHandler<dim> &           dof_handler,
                    const hp::QCollection<dim> &      quad,
                    const AffineConstraints<number> & constraints,
                    VectorType &                      system_rhs);

    void
    reinit(const hp::MappingCollection<dim> &mapping,
           const DoFHandler<dim> &           dof_handler,
           const hp::QCollection<dim> &      quad,
           const AffineConstraints<number> & constraints,
           VectorType &                      system_rhs);

    types::global_dof_index
    m() const;

    number
    el(unsigned int, unsigned int) const;

    void
    initialize_dof_vector(VectorType &vec) const;

    void
    vmult(VectorType &dst, const VectorType &src) const;

    void
    Tvmult(VectorType &dst, const VectorType &src) const;

    const TrilinosWrappers::SparseMatrix &
    get_system_matrix() const;

    void
    compute_inverse_diagonal(VectorType &diagonal) const;

  private:
    void
    do_cell_integral_local(FECellIntegrator &integrator) const;

    void
    do_cell_integral_global(FECellIntegrator &integrator,
                            VectorType &      dst,
                            const VectorType &src) const;

    void
    do_cell_integral_range(
      const MatrixFree<dim, number> &              matrix_free,
      VectorType &                                 dst,
      const VectorType &                           src,
      const std::pair<unsigned int, unsigned int> &range) const;

    MatrixFree<dim, number> matrix_free;

    // 为了用AMG预处理程序解决最粗层次的方程系统，我们需要一个最粗层次的实际系统矩阵。为此，我们提供了一种机制，可以选择从无矩阵公式中计算出一个矩阵，为此我们引入了一个专门的SparseMatrix对象。在默认情况下，这个矩阵保持为空。一旦`get_system_matrix()`被调用，这个矩阵就会被填充（懒惰分配）。由于这是一个
    // "const "函数，我们需要在这里使用 "mutable
    // "关键字。我们还需要一个约束对象来构建矩阵。

    AffineConstraints<number>              constraints;
    mutable TrilinosWrappers::SparseMatrix system_matrix;
  };

  // 下面的部分包含了初始化和重新初始化该类的函数。特别是，这些函数初始化了内部的MatrixFree实例。为了简单起见，我们还计算了系统右侧的向量。

  template <int dim, typename number>
  LaplaceOperator<dim, number>::LaplaceOperator(
    const hp::MappingCollection<dim> &mapping,
    const DoFHandler<dim> &           dof_handler,
    const hp::QCollection<dim> &      quad,
    const AffineConstraints<number> & constraints,
    VectorType &                      system_rhs)
  {
    this->reinit(mapping, dof_handler, quad, constraints, system_rhs);
  }

  template <int dim, typename number>
  void
  LaplaceOperator<dim, number>::reinit(
    const hp::MappingCollection<dim> &mapping,
    const DoFHandler<dim> &           dof_handler,
    const hp::QCollection<dim> &      quad,
    const AffineConstraints<number> & constraints,
    VectorType &                      system_rhs)
  {
    // 清除内部数据结构（在操作者被重复使用的情况下）。

    this->system_matrix.clear();

    // 复制约束条件，因为以后在计算系统矩阵时可能需要它们。

    this->constraints.copy_from(constraints);

    // 设置MatrixFree。在正交点，我们只需要评估解的梯度，并用形状函数的梯度进行测试，所以我们只需要设置标志`update_gradients`。

    typename MatrixFree<dim, number>::AdditionalData data;
    data.mapping_update_flags = update_gradients;

    matrix_free.reinit(mapping, dof_handler, constraints, quad, data);

    // 计算右手边的向量。为此，我们设置了第二个MatrixFree实例，它使用一个修改过的AffineConstraints，不包含由于Dirichlet-边界条件的约束。这个修改过的算子被应用于一个只设置了迪里希特值的向量。其结果是负的右手边向量。

    {
      AffineConstraints<number> constraints_without_dbc;

      IndexSet locally_relevant_dofs;
      DoFTools::extract_locally_relevant_dofs(dof_handler,
                                              locally_relevant_dofs);
      constraints_without_dbc.reinit(locally_relevant_dofs);

      DoFTools::make_hanging_node_constraints(dof_handler,
                                              constraints_without_dbc);
      constraints_without_dbc.close();

      VectorType b, x;

      this->initialize_dof_vector(system_rhs);

      MatrixFree<dim, number> matrix_free;
      matrix_free.reinit(
        mapping, dof_handler, constraints_without_dbc, quad, data);

      matrix_free.initialize_dof_vector(b);
      matrix_free.initialize_dof_vector(x);

      constraints.distribute(x);

      matrix_free.cell_loop(&LaplaceOperator::do_cell_integral_range,
                            this,
                            b,
                            x);

      constraints.set_zero(b);

      system_rhs -= b;
    }
  }

  // 以下函数是多网格算法隐含需要的，包括平滑器。

  // 由于我们没有矩阵，所以要向DoFHandler查询自由度的数量。

  template <int dim, typename number>
  types::global_dof_index
  LaplaceOperator<dim, number>::m() const
  {
    return matrix_free.get_dof_handler().n_dofs();
  }

  // 访问矩阵中的一个特定元素。这个函数既不需要也没有实现，但是，在编译程序时需要它。

  template <int dim, typename number>
  number
  LaplaceOperator<dim, number>::el(unsigned int, unsigned int) const
  {
    Assert(false, ExcNotImplemented());
    return 0;
  }

  // 初始化给定的向量。我们只是把这个任务委托给同名的MatrixFree函数。

  template <int dim, typename number>
  void
  LaplaceOperator<dim, number>::initialize_dof_vector(VectorType &vec) const
  {
    matrix_free.initialize_dof_vector(vec);
  }

  // 在MatrixFree的帮助下，通过在所有单元中循环进行运算评估，并评估单元积分的效果（参见。`do_cell_integral_local()`和`do_cell_integral_global()`）。)

  template <int dim, typename number>
  void
  LaplaceOperator<dim, number>::vmult(VectorType &      dst,
                                      const VectorType &src) const
  {
    this->matrix_free.cell_loop(
      &LaplaceOperator::do_cell_integral_range, this, dst, src, true);
  }

  // 执行转置的运算符评估。由于我们考虑的是对称的
  // "矩阵"，这个函数可以简单地将其任务委托给vmult()。

  template <int dim, typename number>
  void
  LaplaceOperator<dim, number>::Tvmult(VectorType &      dst,
                                       const VectorType &src) const
  {
    this->vmult(dst, src);
  }

  // 由于我们没有一个系统矩阵，我们不能循环计算矩阵的对角线项。相反，我们通过对单位基向量进行一连串的运算符评估来计算对角线。为此，我们使用了MatrixFreeTools命名空间中的一个优化函数。之后再手动进行反转。

  template <int dim, typename number>
  void
  LaplaceOperator<dim, number>::compute_inverse_diagonal(
    VectorType &diagonal) const
  {
    MatrixFreeTools::compute_diagonal(matrix_free,
                                      diagonal,
                                      &LaplaceOperator::do_cell_integral_local,
                                      this);

    for (auto &i : diagonal)
      i = (std::abs(i) > 1.0e-10) ? (1.0 / i) : 1.0;
  }

  // 在无矩阵的情况下，在这个类的初始化过程中没有设置系统矩阵。因此，如果需要的话，它必须在这里被计算出来。由于矩阵在本教程中只对线性元素进行计算（在粗略的网格上），这一点是可以接受的。矩阵的条目是通过运算符的评估序列得到的。为此，使用了优化函数
  // MatrixFreeTools::compute_matrix()
  // 。矩阵只有在尚未设置的情况下才会被计算（懒惰分配）。

  template <int dim, typename number>
  const TrilinosWrappers::SparseMatrix &
  LaplaceOperator<dim, number>::get_system_matrix() const
  {
    if (system_matrix.m() == 0 && system_matrix.n() == 0)
      {
        const auto &dof_handler = this->matrix_free.get_dof_handler();

        TrilinosWrappers::SparsityPattern dsp(
          dof_handler.locally_owned_dofs(),
          dof_handler.get_triangulation().get_communicator());

        DoFTools::make_sparsity_pattern(dof_handler, dsp, this->constraints);

        dsp.compress();
        system_matrix.reinit(dsp);

        MatrixFreeTools::compute_matrix(
          matrix_free,
          constraints,
          system_matrix,
          &LaplaceOperator::do_cell_integral_local,
          this);
      }

    return this->system_matrix;
  }

  // 对一个单元格批处理进行单元格积分，不需要收集和分散数值。MatrixFreeTools函数需要这个函数，因为这些函数直接对FEEvaluation的缓冲区进行操作。

  template <int dim, typename number>
  void
  LaplaceOperator<dim, number>::do_cell_integral_local(
    FECellIntegrator &integrator) const
  {
    integrator.evaluate(EvaluationFlags::gradients);

    for (unsigned int q = 0; q < integrator.n_q_points; ++q)
      integrator.submit_gradient(integrator.get_gradient(q), q);

    integrator.integrate(EvaluationFlags::gradients);
  }

  // 与上述相同，但可以访问全局向量。

  template <int dim, typename number>
  void
  LaplaceOperator<dim, number>::do_cell_integral_global(
    FECellIntegrator &integrator,
    VectorType &      dst,
    const VectorType &src) const
  {
    integrator.gather_evaluate(src, EvaluationFlags::gradients);

    for (unsigned int q = 0; q < integrator.n_q_points; ++q)
      integrator.submit_gradient(integrator.get_gradient(q), q);

    integrator.integrate_scatter(EvaluationFlags::gradients, dst);
  }

  // 这个函数在一个单元格批次范围内的所有单元格批次上循环，并调用上述函数。

  template <int dim, typename number>
  void
  LaplaceOperator<dim, number>::do_cell_integral_range(
    const MatrixFree<dim, number> &              matrix_free,
    VectorType &                                 dst,
    const VectorType &                           src,
    const std::pair<unsigned int, unsigned int> &range) const
  {
    FECellIntegrator integrator(matrix_free, range);

    for (unsigned cell = range.first; cell < range.second; ++cell)
      {
        integrator.reinit(cell);

        do_cell_integral_global(integrator, dst, src);
      }
  }

  //  @sect3{Solver and preconditioner}
  // @sect4{Conjugate-gradient solver with multigrid preconditioner}

  // 这个函数用一连串提供的多网格对象来解决方程组。它的目的是为了尽可能的通用，因此有许多模板参数。

  template <typename VectorType,
            int dim,
            typename SystemMatrixType,
            typename LevelMatrixType,
            typename MGTransferType>
  static void
  mg_solve(SolverControl &            solver_control,
           VectorType &               dst,
           const VectorType &         src,
           const MultigridParameters &mg_data,
           const DoFHandler<dim> &    dof,
           const SystemMatrixType &   fine_matrix,
           const MGLevelObject<std::unique_ptr<LevelMatrixType>> &mg_matrices,
           const MGTransferType &                                 mg_transfer)
  {
    AssertThrow(mg_data.coarse_solver.type == "cg_with_amg",
                ExcNotImplemented());
    AssertThrow(mg_data.smoother.type == "chebyshev", ExcNotImplemented());

    const unsigned int min_level = mg_matrices.min_level();
    const unsigned int max_level = mg_matrices.max_level();

    using SmootherPreconditionerType = DiagonalMatrix<VectorType>;
    using SmootherType               = PreconditionChebyshev<LevelMatrixType,
                                               VectorType,
                                               SmootherPreconditionerType>;
    using PreconditionerType = PreconditionMG<dim, VectorType, MGTransferType>;

    // 我们在这里初始化电平运算符和切比雪夫平滑器。

    mg::Matrix<VectorType> mg_matrix(mg_matrices);

    MGLevelObject<typename SmootherType::AdditionalData> smoother_data(
      min_level, max_level);

    for (unsigned int level = min_level; level <= max_level; level++)
      {
        smoother_data[level].preconditioner =
          std::make_shared<SmootherPreconditionerType>();
        mg_matrices[level]->compute_inverse_diagonal(
          smoother_data[level].preconditioner->get_vector());
        smoother_data[level].smoothing_range = mg_data.smoother.smoothing_range;
        smoother_data[level].degree          = mg_data.smoother.degree;
        smoother_data[level].eig_cg_n_iterations =
          mg_data.smoother.eig_cg_n_iterations;
      }

    MGSmootherPrecondition<LevelMatrixType, SmootherType, VectorType>
      mg_smoother;
    mg_smoother.initialize(mg_matrices, smoother_data);

    // 接下来，我们初始化粗略网格求解器。我们使用共轭梯度法和AMG作为预处理程序。

    ReductionControl coarse_grid_solver_control(mg_data.coarse_solver.maxiter,
                                                mg_data.coarse_solver.abstol,
                                                mg_data.coarse_solver.reltol,
                                                false,
                                                false);
    SolverCG<VectorType> coarse_grid_solver(coarse_grid_solver_control);

    std::unique_ptr<MGCoarseGridBase<VectorType>> mg_coarse;

    TrilinosWrappers::PreconditionAMG                 precondition_amg;
    TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
    amg_data.smoother_sweeps = mg_data.coarse_solver.smoother_sweeps;
    amg_data.n_cycles        = mg_data.coarse_solver.n_cycles;
    amg_data.smoother_type   = mg_data.coarse_solver.smoother_type.c_str();

    precondition_amg.initialize(mg_matrices[min_level]->get_system_matrix(),
                                amg_data);

    mg_coarse =
      std::make_unique<MGCoarseGridIterativeSolver<VectorType,
                                                   SolverCG<VectorType>,
                                                   LevelMatrixType,
                                                   decltype(precondition_amg)>>(
        coarse_grid_solver, *mg_matrices[min_level], precondition_amg);

    // 最后，我们创建Multigrid对象，将其转换为预处理程序，并在共轭梯度求解器中使用它来解决线性方程组。

    Multigrid<VectorType> mg(
      mg_matrix, *mg_coarse, mg_transfer, mg_smoother, mg_smoother);

    PreconditionerType preconditioner(dof, mg, mg_transfer);

    SolverCG<VectorType>(solver_control)
      .solve(fine_matrix, dst, src, preconditioner);
  }

  //  @sect4{Hybrid polynomial/geometric-global-coarsening multigrid
  //  preconditioner}

  // 上述函数处理给定的多网格对象序列的实际解决方案。这个函数创建了实际的多重网格层次，特别是运算符，以及作为MGTransferGlobalCoarsening对象的转移运算符。

  template <typename VectorType, typename OperatorType, int dim>
  void
  solve_with_gmg(SolverControl &                  solver_control,
                 const OperatorType &             system_matrix,
                 VectorType &                     dst,
                 const VectorType &               src,
                 const MultigridParameters &      mg_data,
                 const hp::MappingCollection<dim> mapping_collection,
                 const DoFHandler<dim> &          dof_handler,
                 const hp::QCollection<dim> &     quadrature_collection)
  {
    // 为每个多网格层次创建一个DoFHandler和操作符，以及，创建转移操作符。为了能够设置运算符，我们需要一组DoFHandler，通过p或h的全局粗化来创建。

    // 如果没有要求h-transfer，我们为`emplace_back()`函数提供一个空的删除器，因为我们的DoFHandler的Triangulation是一个外部字段，其析构器在其他地方被调用。

    MGLevelObject<DoFHandler<dim>>                     dof_handlers;
    MGLevelObject<std::unique_ptr<OperatorType>>       operators;
    MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> transfers;

    std::vector<std::shared_ptr<const Triangulation<dim>>>
      coarse_grid_triangulations;
    if (mg_data.transfer.perform_h_transfer)
      coarse_grid_triangulations =
        MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
          dof_handler.get_triangulation());
    else
      coarse_grid_triangulations.emplace_back(
        const_cast<Triangulation<dim> *>(&(dof_handler.get_triangulation())),
        [](auto &) {});

    // 确定多栅格操作的总层数，并为所有层数分配足够的内存。

    const unsigned int n_h_levels = coarse_grid_triangulations.size() - 1;

    const auto get_max_active_fe_degree = [&](const auto &dof_handler) {
      unsigned int max = 0;

      for (auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
          max =
            std::max(max, dof_handler.get_fe(cell->active_fe_index()).degree);

      return Utilities::MPI::max(max, MPI_COMM_WORLD);
    };

    const unsigned int n_p_levels =
      MGTransferGlobalCoarseningTools::create_polynomial_coarsening_sequence(
        get_max_active_fe_degree(dof_handler), mg_data.transfer.p_sequence)
        .size();

    std::map<unsigned int, unsigned int> fe_index_for_degree;
    for (unsigned int i = 0; i < dof_handler.get_fe_collection().size(); ++i)
      {
        const unsigned int degree = dof_handler.get_fe(i).degree;
        Assert(fe_index_for_degree.find(degree) == fe_index_for_degree.end(),
               ExcMessage("FECollection does not contain unique degrees."));
        fe_index_for_degree[degree] = i;
      }

    unsigned int minlevel   = 0;
    unsigned int minlevel_p = n_h_levels;
    unsigned int maxlevel   = n_h_levels + n_p_levels - 1;

    dof_handlers.resize(minlevel, maxlevel);
    operators.resize(minlevel, maxlevel);
    transfers.resize(minlevel, maxlevel);

    // 从最小（最粗）到最大（最细）级别的循环，并相应地设置DoFHandler。我们从h层开始，在这里我们分布在越来越细的网格上的线性元素。

    for (unsigned int l = 0; l < n_h_levels; ++l)
      {
        dof_handlers[l].reinit(*coarse_grid_triangulations[l]);
        dof_handlers[l].distribute_dofs(dof_handler.get_fe_collection());
      }

    // 在我们达到最细的网格后，我们将调整每一层的多项式度数。我们反向迭代我们的数据结构，从包含所有活动FE指数信息的最细网格开始。然后我们逐级降低每个单元的多项式度数。

    for (unsigned int i = 0, l = maxlevel; i < n_p_levels; ++i, --l)
      {
        dof_handlers[l].reinit(dof_handler.get_triangulation());

        if (l == maxlevel) // finest level
          {
            auto &dof_handler_mg = dof_handlers[l];

            auto cell_other = dof_handler.begin_active();
            for (auto &cell : dof_handler_mg.active_cell_iterators())
              {
                if (cell->is_locally_owned())
                  cell->set_active_fe_index(cell_other->active_fe_index());
                cell_other++;
              }
          }
        else // coarse level
          {
            auto &dof_handler_fine   = dof_handlers[l + 1];
            auto &dof_handler_coarse = dof_handlers[l + 0];

            auto cell_other = dof_handler_fine.begin_active();
            for (auto &cell : dof_handler_coarse.active_cell_iterators())
              {
                if (cell->is_locally_owned())
                  {
                    const unsigned int next_degree =
                      MGTransferGlobalCoarseningTools::
                        create_next_polynomial_coarsening_degree(
                          cell_other->get_fe().degree,
                          mg_data.transfer.p_sequence);
                    Assert(fe_index_for_degree.find(next_degree) !=
                             fe_index_for_degree.end(),
                           ExcMessage("Next polynomial degree in sequence "
                                      "does not exist in FECollection."));

                    cell->set_active_fe_index(fe_index_for_degree[next_degree]);
                  }
                cell_other++;
              }
          }

        dof_handlers[l].distribute_dofs(dof_handler.get_fe_collection());
      }

    // 接下来，我们将在每个多重网格层面上创建所有额外需要的数据结构。这涉及到确定具有同质Dirichlet边界条件的约束，并像在活动层上一样建立运算器。

    MGLevelObject<AffineConstraints<typename VectorType::value_type>>
      constraints(minlevel, maxlevel);

    for (unsigned int level = minlevel; level <= maxlevel; ++level)
      {
        const auto &dof_handler = dof_handlers[level];
        auto &      constraint  = constraints[level];

        IndexSet locally_relevant_dofs;
        DoFTools::extract_locally_relevant_dofs(dof_handler,
                                                locally_relevant_dofs);
        constraint.reinit(locally_relevant_dofs);

        DoFTools::make_hanging_node_constraints(dof_handler, constraint);
        VectorTools::interpolate_boundary_values(mapping_collection,
                                                 dof_handler,
                                                 0,
                                                 Functions::ZeroFunction<dim>(),
                                                 constraint);
        constraint.close();

        VectorType dummy;

        operators[level] = std::make_unique<OperatorType>(mapping_collection,
                                                          dof_handler,
                                                          quadrature_collection,
                                                          constraint,
                                                          dummy);
      }

    //根据多网格求解器类的需要，在单个算子中设置网格间算子和收集转移算子。

    for (unsigned int level = minlevel; level < minlevel_p; ++level)
      transfers[level + 1].reinit_geometric_transfer(dof_handlers[level + 1],
                                                     dof_handlers[level],
                                                     constraints[level + 1],
                                                     constraints[level]);

    for (unsigned int level = minlevel_p; level < maxlevel; ++level)
      transfers[level + 1].reinit_polynomial_transfer(dof_handlers[level + 1],
                                                      dof_handlers[level],
                                                      constraints[level + 1],
                                                      constraints[level]);

    MGTransferGlobalCoarsening<dim, VectorType> transfer(
      transfers, [&](const auto l, auto &vec) {
        operators[l]->initialize_dof_vector(vec);
      });

    // 最后，继续用多网格法解决问题。

    mg_solve(solver_control,
             dst,
             src,
             mg_data,
             dof_handler,
             system_matrix,
             operators,
             transfer);
  }

  //  @sect3{The <code>LaplaceProblem</code> class template}

  // 现在，我们将最后声明这个程序的主类，它在随后的精炼函数空间上求解拉普拉斯方程。它的结构看起来很熟悉，因为它与
  // step-27  和  step-40  的主类类似。基本上只增加了两个。

  // -
  // 持有系统矩阵的SparseMatrix对象已经被MatrixFree公式中的LaplaceOperator类对象所取代。

  // - 加入了一个 parallel::CellWeights, 的对象，它将帮助我们实现负载平衡。

  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem(const Parameters &parameters);

    void
    run();

  private:
    void
    initialize_grid();
    void
    setup_system();
    void
    print_diagnostics();
    void
    solve_system();
    void
    compute_indicators();
    void
    adapt_resolution();
    void
    output_results(const unsigned int cycle);

    MPI_Comm mpi_communicator;

    const Parameters prm;

    parallel::distributed::Triangulation<dim> triangulation;
    DoFHandler<dim>                           dof_handler;

    hp::MappingCollection<dim> mapping_collection;
    hp::FECollection<dim>      fe_collection;
    hp::QCollection<dim>       quadrature_collection;
    hp::QCollection<dim - 1>   face_quadrature_collection;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    AffineConstraints<double> constraints;

    LaplaceOperator<dim, double>               laplace_operator;
    LinearAlgebra::distributed::Vector<double> locally_relevant_solution;
    LinearAlgebra::distributed::Vector<double> system_rhs;

    std::unique_ptr<FESeries::Legendre<dim>>    legendre;
    std::unique_ptr<parallel::CellWeights<dim>> cell_weights;

    Vector<float> estimated_error_per_cell;
    Vector<float> hp_decision_indicators;

    ConditionalOStream pcout;
    TimerOutput        computing_timer;
  };

  //  @sect3{The <code>LaplaceProblem</code> class implementation}
  // @sect4{Constructor}

  // 构造函数以一个初始化器列表开始，该列表看起来与  step-40
  // 的列表相似。我们再次准备好ConditionalOStream对象，只允许第一个进程在控制台输出任何东西，并正确初始化计算计时器。

  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem(const Parameters &parameters)
    : mpi_communicator(MPI_COMM_WORLD)
    , prm(parameters)
    , triangulation(mpi_communicator)
    , dof_handler(triangulation)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
  {
    Assert(prm.min_h_level <= prm.max_h_level,
           ExcMessage(
             "Triangulation level limits have been incorrectly set up."));
    Assert(prm.min_p_degree <= prm.max_p_degree,
           ExcMessage("FECollection degrees have been incorrectly set up."));

    // 我们需要在构造函数的实际主体中为hp-functionality准备数据结构，并在参数结构的指定范围内为每个度数创建相应的对象。由于我们只处理非扭曲的矩形单元，在这种情况下，一个线性映射对象就足够了。

    // 在参数结构中，我们为函数空间以合理的分辨率运行的层级提供范围。多网格算法需要在最粗的层次上使用线性元素。所以我们从最低的多项式度数开始，用连续的高度数填充集合，直到达到用户指定的最大值。

    mapping_collection.push_back(MappingQ1<dim>());

    for (unsigned int degree = 1; degree <= prm.max_p_degree; ++degree)
      {
        fe_collection.push_back(FE_Q<dim>(degree));
        quadrature_collection.push_back(QGauss<dim>(degree + 1));
        face_quadrature_collection.push_back(QGauss<dim - 1>(degree + 1));
      }

    // 由于我们的FECollection包含的有限元比我们想用于求解的有限元近似值要多，我们想限制活动FE指数可以操作的范围。为此，FECollection类允许注册一个层次结构，在p-精简和p-粗化的情况下，分别决定后续的和前面的有限元。
    // hp::Refinement
    // 命名空间中的所有函数都会参考这个层次结构来确定未来的FE指数。我们将注册这样一个层次结构，它只对建议范围内的多项式程度的有限元起作用
    // <code>[min_p_degree, max_p_degree]</code>  。

    const unsigned int min_fe_index = prm.min_p_degree - 1;
    fe_collection.set_hierarchy(

      /*下一个_index=  */
      [](const typename hp::FECollection<dim> &fe_collection,
         const unsigned int                    fe_index) -> unsigned int {
        return ((fe_index + 1) < fe_collection.size()) ? fe_index + 1 :
                                                         fe_index;
      },
      /*上一页_index=  */
      [min_fe_index](const typename hp::FECollection<dim> &,
                     const unsigned int fe_index) -> unsigned int {
        Assert(fe_index >= min_fe_index,
               ExcMessage("Finite element is not part of hierarchy!"));
        return (fe_index > min_fe_index) ? fe_index - 1 : fe_index;
      });

    // 我们以默认配置初始化 FESeries::Legendre 对象，以便进行平滑度估计。

    legendre = std::make_unique<FESeries::Legendre<dim>>(
      SmoothnessEstimator::Legendre::default_fe_series(fe_collection));

    // 接下来的部分会很棘手。在执行细化的过程中，有几个hp-算法需要干扰三角形对象上的实际细化过程。我们通过将几个函数连接到
    // Triangulation::Signals:
    // 信号，在实际细化过程中的不同阶段被调用，并触发所有连接的函数来做到这一点。我们需要这个功能来实现负载平衡和限制相邻单元的多项式度数。

    // 对于前者，我们希望给每个单元分配一个权重，这个权重与它未来的有限元的自由度数成正比。该库提供了一个类
    // parallel::CellWeights
    // ，允许在细化过程中的正确位置轻松地附加单个权重，即在所有细化和粗化标志被正确设置为hp-adaptation之后，以及在即将发生的负载平衡的重新划分之前。可以注册一些函数，这些函数将以
    // $a (n_\text{dofs})^b$  提供的一对参数的形式附加权重  $(a,b)$
    // 。我们在下文中注册了这样一个函数。每个单元在创建时将被赋予一个恒定的权重，这个值是1000（见
    // Triangulation::Signals::cell_weight).  ）。

    // 为了实现负载平衡，像我们使用的高效求解器应该与拥有的自由度数量成线性比例。此外，为了增加我们想要附加的权重的影响，确保单个权重将超过这个基础权重的数量级。我们相应地设置单元加权的参数。大的加权系数为
    // $10^6$ ，指数为 $1$  。

    cell_weights = std::make_unique<parallel::CellWeights<dim>>(
      dof_handler,
      parallel::CellWeights<dim>::ndofs_weighting(
        {prm.weighting_factor, prm.weighting_exponent}));

    // 在h-adaptive应用中，我们通过限制相邻单元的细化水平的差异为1来确保2:1的网格平衡。通过下面代码片段中的第二个调用，我们将确保相邻单元的p级数也是如此：未来有限元的级数不允许相差超过指定的差值。函数
    // hp::Refinement::limit_p_level_difference
    // 可以处理这个问题，但需要与并行环境中的一个非常特殊的信号相连。问题是，我们需要知道网格的实际细化情况，以便相应地设置未来的FE指数。由于我们要求p4est神谕进行细化，我们需要确保Triangulation已经先用神谕的适应标志进行了更新。
    // parallel::distributed::TemporarilyMatchRefineFlags
    // 的实例化在其生命期内正是如此。因此，我们将在限制p级差之前创建这个类的对象，并将相应的lambda函数连接到信号
    // Triangulation::Signals::post_p4est_refinement,
    // 上，该信号将在神谕被完善之后，但在三角法被完善之前被触发。此外，我们指定这个函数将被连接到信号的前面，以确保修改在连接到同一信号的任何其他函数之前进行。

    triangulation.signals.post_p4est_refinement.connect(
      [&, min_fe_index]() {
        const parallel::distributed::TemporarilyMatchRefineFlags<dim>
          refine_modifier(triangulation);
        hp::Refinement::limit_p_level_difference(dof_handler,
                                                 prm.max_p_level_difference,

                                                 //包含=  */

                                                 min_fe_index）。)
      },
      boost::signals2::at_front);
  }

  //  @sect4{LaplaceProblem::initialize_grid}

  // 对于L型域，我们可以使用 GridGenerator::hyper_L() 这个函数，如 step-50
  // 中所演示的。然而在二维的情况下，该函数只去除第一象限，而在我们的方案中我们需要去除第四象限。因此，我们将使用一个不同的函数
  // GridGenerator::subdivided_hyper_L()
  // ，它给我们更多的选择来创建网格。此外，我们在制定该函数时，也会生成一个三维网格：二维L型域基本上会在正Z方向上拉长1。

  // 我们首先假装建立一个  GridGenerator::subdivided_hyper_rectangle().
  // 我们需要提供的参数是左下角和右上角的点对象，以及基本网格在每个方向的重复次数。我们为前两个维度提供这些参数，对更高的第三维度单独处理。

  // 为了创建一个L型域，我们需要去除多余的单元。为此，我们相应地指定
  // <code>cells_to_remove</code>
  // 。我们希望从负方向的每一个单元格中移除一个单元格，但从正的x方向移除一个。

  // 最后，我们提供与所提供的最小网格细化水平相对应的初始细化数。此外，我们相应地设置初始活动FE指数。

  template <int dim>
  void
  LaplaceProblem<dim>::initialize_grid()
  {
    TimerOutput::Scope t(computing_timer, "initialize grid");

    std::vector<unsigned int> repetitions(dim);
    Point<dim>                bottom_left, top_right;
    for (unsigned int d = 0; d < dim; ++d)
      if (d < 2)
        {
          repetitions[d] = 2;
          bottom_left[d] = -1.;
          top_right[d]   = 1.;
        }
      else
        {
          repetitions[d] = 1;
          bottom_left[d] = 0.;
          top_right[d]   = 1.;
        }

    std::vector<int> cells_to_remove(dim, 1);
    cells_to_remove[0] = -1;

    GridGenerator::subdivided_hyper_L(
      triangulation, repetitions, bottom_left, top_right, cells_to_remove);

    triangulation.refine_global(prm.min_h_level);

    const unsigned int min_fe_index = prm.min_p_degree - 1;
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        cell->set_active_fe_index(min_fe_index);
  }

  //  @sect4{LaplaceProblem::setup_system}

  // 这个函数看起来和 step-40
  // 的函数完全一样，但是你会注意到没有系统矩阵以及围绕它的脚手架。相反，我们将在这里初始化
  // <code>laplace_operator</code>
  // 中的MatrixFree公式。对于边界条件，我们将使用本教程前面介绍的Solution类。

  template <int dim>
  void
  LaplaceProblem<dim>::setup_system()
  {
    TimerOutput::Scope t(computing_timer, "setup system");

    dof_handler.distribute_dofs(fe_collection);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    locally_relevant_solution.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);
    system_rhs.reinit(locally_owned_dofs, mpi_communicator);

    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(
      mapping_collection, dof_handler, 0, Solution<dim>(), constraints);
    constraints.close();

    laplace_operator.reinit(mapping_collection,
                            dof_handler,
                            quadrature_collection,
                            constraints,
                            system_rhs);
  }

  //  @sect4{LaplaceProblem::print_diagnostics}

  // 这是一个打印关于方程组及其划分的额外诊断的函数。除了通常的全局活动单元数和自由度外，我们还输出它们的局部等价物。为了规范输出，我们将用
  // Utilities::MPI::gather
  // 操作将局部数量传达给第一个进程，然后由该进程输出所有信息。本地量的输出只限于前8个进程，以避免终端的杂乱。

  // 此外，我们想打印数值离散化中的多项式度数的频率。由于这些信息只存储在本地，我们将计算本地拥有的单元上的有限元，随后通过
  // Utilities::MPI::sum. 进行交流。
  template <int dim>
  void
  LaplaceProblem<dim>::print_diagnostics()
  {
    const unsigned int first_n_processes =
      std::min<unsigned int>(8,
                             Utilities::MPI::n_mpi_processes(mpi_communicator));
    const bool output_cropped =
      first_n_processes < Utilities::MPI::n_mpi_processes(mpi_communicator);

    {
      pcout << "   Number of active cells:       "
            << triangulation.n_global_active_cells() << std::endl
            << "     by partition:              ";

      std::vector<unsigned int> n_active_cells_per_subdomain =
        Utilities::MPI::gather(mpi_communicator,
                               triangulation.n_locally_owned_active_cells());
      for (unsigned int i = 0; i < first_n_processes; ++i)
        pcout << ' ' << n_active_cells_per_subdomain[i];
      if (output_cropped)
        pcout << " ...";
      pcout << std::endl;
    }

    {
      pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl
            << "     by partition:              ";

      std::vector<types::global_dof_index> n_dofs_per_subdomain =
        Utilities::MPI::gather(mpi_communicator,
                               dof_handler.n_locally_owned_dofs());
      for (unsigned int i = 0; i < first_n_processes; ++i)
        pcout << ' ' << n_dofs_per_subdomain[i];
      if (output_cropped)
        pcout << " ...";
      pcout << std::endl;
    }

    {
      std::vector<types::global_dof_index> n_constraints_per_subdomain =
        Utilities::MPI::gather(mpi_communicator, constraints.n_constraints());

      pcout << "   Number of constraints:        "
            << std::accumulate(n_constraints_per_subdomain.begin(),
                               n_constraints_per_subdomain.end(),
                               0)
            << std::endl
            << "     by partition:              ";
      for (unsigned int i = 0; i < first_n_processes; ++i)
        pcout << ' ' << n_constraints_per_subdomain[i];
      if (output_cropped)
        pcout << " ...";
      pcout << std::endl;
    }

    {
      std::vector<unsigned int> n_fe_indices(fe_collection.size(), 0);
      for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
          n_fe_indices[cell->active_fe_index()]++;

      Utilities::MPI::sum(n_fe_indices, mpi_communicator, n_fe_indices);

      pcout << "   Frequencies of poly. degrees:";
      for (unsigned int i = 0; i < fe_collection.size(); ++i)
        if (n_fe_indices[i] > 0)
          pcout << ' ' << fe_collection[i].degree << ":" << n_fe_indices[i];
      pcout << std::endl;
    }
  }

  //  @sect4{LaplaceProblem::solve_system}

  // 围绕解决方案的脚手架与  step-40
  // 的类似。我们准备一个符合MatrixFree要求的向量，并收集本地相关的自由度，我们解决了方程系统。解决方法是通过前面介绍的函数进行的。

  template <int dim>
  void
  LaplaceProblem<dim>::solve_system()
  {
    TimerOutput::Scope t(computing_timer, "solve system");

    LinearAlgebra::distributed::Vector<double> completely_distributed_solution;
    laplace_operator.initialize_dof_vector(completely_distributed_solution);

    SolverControl solver_control(system_rhs.size(),
                                 prm.tolerance_factor * system_rhs.l2_norm());

    solve_with_gmg(solver_control,
                   laplace_operator,
                   completely_distributed_solution,
                   system_rhs,
                   prm.mg_data,
                   mapping_collection,
                   dof_handler,
                   quadrature_collection);

    pcout << "   Solved in " << solver_control.last_step() << " iterations."
          << std::endl;

    constraints.distribute(completely_distributed_solution);

    locally_relevant_solution.copy_locally_owned_data_from(
      completely_distributed_solution);
    locally_relevant_solution.update_ghost_values();
  }

  //  @sect4{LaplaceProblem::compute_indicators}

  // 这个函数只包含其他教程中典型的 <code>refine_grid</code>
  // 函数的一部分，在这个意义上是新的。在这里，我们将只计算与实际细化网格相适应的所有指标。我们这样做的目的是将所有的指标写到文件系统中，以便为以后储存。

  // 由于我们处理的是一个椭圆问题，我们将再次利用KellyErrorEstimator，但有一点不同。修改底层面积分的缩放系数，使其取决于相邻元素的实际多项式程度，这对hp-adaptive应用是有利的
  // @cite davydov2017hp
  // 。我们可以通过指定你所注意到的附加参数中的最后一个参数来做到这一点。其他的实际上只是默认的。

  // 为了hp-adaptation的目的，我们将用教程介绍中的策略来计算平滑度估计，并使用
  // SmoothnessEstimator::Legendre. 中的实现
  // 在参数结构中，我们将最小多项式度数设置为2，因为似乎平滑度估计算法在处理线性元素时有问题。

  template <int dim>
  void
  LaplaceProblem<dim>::compute_indicators()
  {
    TimerOutput::Scope t(computing_timer, "compute indicators");

    estimated_error_per_cell.grow_or_shrink(triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate(
      dof_handler,
      face_quadrature_collection,
      std::map<types::boundary_id, const Function<dim> *>(),
      locally_relevant_solution,
      estimated_error_per_cell,
      /*component_mask=  */
      ComponentMask(),
      /*coefficients=  */
      nullptr,
      /*n_threads=*/
      numbers::invalid_unsigned_int,
      /*subdomain_id=*/
      numbers::invalid_subdomain_id,
      /*material_id=*/
      numbers::invalid_material_id,
      /*策略=  */
      KellyErrorEstimator<dim>::Strategy::face_diameter_over_twice_max_degree);

    hp_decision_indicators.grow_or_shrink(triangulation.n_active_cells());
    SmoothnessEstimator::Legendre::coefficient_decay(*legendre,
                                                     dof_handler,
                                                     locally_relevant_solution,
                                                     hp_decision_indicators);
  }

  //  @sect4{LaplaceProblem::adapt_resolution}

  // 有了之前计算出的指标，我们最终将标记所有单元进行适应，同时在这个函数中执行细化。和以前的教程一样，我们将使用
  // "固定数字 "策略，但现在是针对hp-adaptation。

  template <int dim>
  void
  LaplaceProblem<dim>::adapt_resolution()
  {
    TimerOutput::Scope t(computing_timer, "adapt resolution");

    // 首先，我们将根据每个单元的误差估计值来设置细化和粗化标志。这里没有什么新东西。

    // 我们将使用在其他deal.II教程中阐述过的一般细化和粗化比例：使用固定数字策略，我们将标记所有单元中的30%进行细化，3%进行粗化，如参数结构中提供的。

    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
      triangulation,
      estimated_error_per_cell,
      prm.refine_fraction,
      prm.coarsen_fraction);

    // 接下来，我们将对hp-adaptation进行所有调整。我们想细化和粗化那些在上一步中被标记的单元，但需要决定是通过调整网格分辨率还是调整多项式程度来实现。

    // 下一个函数调用根据之前计算的平滑度指标设置未来的FE指数，作为p-adaptation指标。这些指数将只设置在那些分配了细化或粗化标志的单元上。

    // 对于p-adaptation分数，我们将采取一个有根据的猜测。由于我们只期望在我们的方案中出现一个单一的奇点，即在域的原点，而在其他任何地方都有一个平滑的解决方案，所以我们希望强烈倾向于使用p-adaptation而不是h-adaptation。这反映在我们对p-精简和p-粗化都选择了90%的分数。

    hp::Refinement::p_adaptivity_fixed_number(dof_handler,
                                              hp_decision_indicators,
                                              prm.p_refine_fraction,
                                              prm.p_coarsen_fraction);

    // 在这个阶段，我们既有未来的FE指数，也有经典的细化和粗化标志，后者将由
    // Triangulation::execute_coarsening_and_refinement()
    // 解释为h-适应性。我们希望只对细胞施加一种适应，这就是下一个函数将为我们解决的问题。简而言之，在分配有两种类型指标的单元格上，我们将倾向于p-适应的那一种，并删除h-适应的那一种。

    hp::Refinement::choose_p_over_h(dof_handler);

    // 设置完所有指标后，我们将删除那些超过参数结构中提供的水平范围的指定限制的指标。由于提供的有限元数量有限，这种限制自然会出现在p-adaptation中。此外，我们在构造函数中为p-adaptation注册了一个自定义层次结构。现在，我们需要像
    // step-31  中那样，在h-adaptive的上下文中手动完成。

    // 我们将遍历指定的最小和最大层次上的所有单元格，并删除相应的标志。作为一种选择，我们也可以通过相应地设置未来的FE指数来标记这些单元的p适应性，而不是简单地清除细化和粗化的标志。

    Assert(triangulation.n_levels() >= prm.min_h_level + 1 &&
             triangulation.n_levels() <= prm.max_h_level + 1,
           ExcInternalError());

    if (triangulation.n_levels() > prm.max_h_level)
      for (const auto &cell :
           triangulation.active_cell_iterators_on_level(prm.max_h_level))
        cell->clear_refine_flag();

    for (const auto &cell :
         triangulation.active_cell_iterators_on_level(prm.min_h_level))
      cell->clear_coarsen_flag();

    // 最后，我们就剩下执行粗化和细化了。在这里，不仅网格会被更新，而且所有以前的未来FE指数也会变得活跃。

    // 记得我们在构造函数中为三角化信号附加了函数，将在这个函数调用中被触发。所以会有更多的事情发生：加权重新分区将被执行以确保负载平衡，以及我们将限制相邻单元之间的p级差。

    triangulation.execute_coarsening_and_refinement();
  }

  //  @sect4{LaplaceProblem::output_results}

  // 在并行应用中向文件系统写入结果的工作方式与  step-40
  // 中完全相同。除了我们在整个教程中准备的数据容器外，我们还想写出网格上每个有限元的多项式程度，以及每个单元所属的子域。我们在这个函数的范围内为此准备必要的容器。

  template <int dim>
  void
  LaplaceProblem<dim>::output_results(const unsigned int cycle)
  {
    TimerOutput::Scope t(computing_timer, "output results");

    Vector<float> fe_degrees(triangulation.n_active_cells());
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        fe_degrees(cell->active_cell_index()) = cell->get_fe().degree;

    Vector<float> subdomain(triangulation.n_active_cells());
    for (auto &subd : subdomain)
      subd = triangulation.locally_owned_subdomain();

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(locally_relevant_solution, "solution");
    data_out.add_data_vector(fe_degrees, "fe_degree");
    data_out.add_data_vector(subdomain, "subdomain");
    data_out.add_data_vector(estimated_error_per_cell, "error");
    data_out.add_data_vector(hp_decision_indicators, "hp_indicator");
    data_out.build_patches(mapping_collection);

    data_out.write_vtu_with_pvtu_record(
      "./", "solution", cycle, mpi_communicator, 2, 1);
  }

  //  @sect4{LaplaceProblem::run}

  // 实际的运行函数看起来又和  step-40
  // 非常相似。唯一增加的是实际循环之前的括号内的部分。在这里，我们将预先计算Legendre变换矩阵。一般来说，每当需要某个矩阵时，这些矩阵将通过懒惰分配的方式进行实时计算。然而，出于计时的目的，我们希望在实际的时间测量开始之前，一次性地计算它们。因此，我们将把它们的计算指定为自己的范围。

  template <int dim>
  void
  LaplaceProblem<dim>::run()
  {
    pcout << "Running with Trilinos on "
          << Utilities::MPI::n_mpi_processes(mpi_communicator)
          << " MPI rank(s)..." << std::endl;

    {
      pcout << "Calculating transformation matrices..." << std::endl;
      TimerOutput::Scope t(computing_timer, "calculate transformation");
      legendre->precalculate_all_transformation_matrices();
    }

    for (unsigned int cycle = 0; cycle < prm.n_cycles; ++cycle)
      {
        pcout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          initialize_grid();
        else
          adapt_resolution();

        setup_system();

        print_diagnostics();

        solve_system();

        compute_indicators();

        if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32)
          output_results(cycle);

        computing_timer.print_summary();
        computing_timer.reset();

        pcout << std::endl;
      }
  }
} // namespace Step75

//  @sect4{main()}

// 最后一个函数是 <code>main</code>
// 函数，它将最终创建并运行一个LaplaceOperator实例。它的结构与其他大多数教程程序相似。

int
main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step75;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      Parameters        prm;
      LaplaceProblem<2> laplace_problem(prm);
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
