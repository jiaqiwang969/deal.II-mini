CCTest_file/step-66.cc

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
 * Authors: Fabian Castelli, Karlsruhe Institute of Technology (KIT) 
 */ 



// 首先我们包括本教程所需的deal.II库的典型头文件。

#include <deal.II/base/function.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/base/vectorization.h> 

#include <deal.II/dofs/dof_accessor.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/mapping_q_generic.h> 

#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/manifold_lib.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/tria_accessor.h> 
#include <deal.II/grid/tria_iterator.h> 

#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/vector.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 

// 特别是，我们需要包括无矩阵框架的头文件。

#include <deal.II/matrix_free/fe_evaluation.h> 
#include <deal.II/matrix_free/matrix_free.h> 
#include <deal.II/matrix_free/operators.h> 
#include <deal.II/matrix_free/tools.h> 

// 由于我们要使用几何多网格预处理程序，所以我们还需要多级头文件。

#include <deal.II/multigrid/mg_coarse.h> 
#include <deal.II/multigrid/mg_constrained_dofs.h> 
#include <deal.II/multigrid/mg_matrix.h> 
#include <deal.II/multigrid/mg_smoother.h> 
#include <deal.II/multigrid/mg_tools.h> 
#include <deal.II/multigrid/mg_transfer_matrix_free.h> 
#include <deal.II/multigrid/multigrid.h> 

// 最后是一些常用的C++头文件，用于输入和输出。

#include <fstream> 
#include <iostream> 

namespace Step66 
{ 
  using namespace dealii; 

//  @sect3{Matrix-free JacobianOperator}  

// 在开始时，我们定义了雅各布系数的无矩阵算子。作为指导，我们遵循教程  step-37  和  step-48  ，其中广泛记录了  MatrixFreeOperators::Base  类的精确接口。

// 由于我们希望将雅各布（Jacobian）作为系统矩阵使用，并将其传递给线性求解器以及多级预处理类，我们从 MatrixFreeOperators::Base 类派生出 <code>JacobianOperator</code> 类，这样我们就有了正确的接口。我们需要从基类中覆盖的两个函数是 MatrixFreeOperators::Base::apply_add() 和 MatrixFreeOperators::Base::compute_diagonal() 函数。为了允许用浮动精度进行预处理，我们将数字类型定义为模板参数。

// 正如在介绍中提到的，我们需要在最后一个牛顿步骤 $u_h^n$ 中评估雅各布 $F'$ ，以便计算牛顿更新 $s_h^n$ 。为了获得最后一个牛顿步骤 $u_h^n$ 的信息，我们的做法与 step-37 基本相同，在使用无矩阵算子之前，我们将一个系数函数的值存储在一个表中 <code>nonlinear_values</code> 。我们在这里实现的不是一个函数  <code>evaluate_coefficient()</code>  ，而是一个函数  <code>evaluate_newton_step()</code>  。

// 作为 <code>JacobianOperator</code> 的额外私有成员函数，我们实现了 <code>local_apply()</code> 和 <code>local_compute_diagonal()</code> 函数。第一个是矩阵-向量应用的实际工作函数，我们在 <code>apply_add()</code> 函数中将其传递给 MatrixFree::cell_loop() 。后面一个是计算对角线的工作函数，我们把它传递给 MatrixFreeTools::compute_diagonal() 函数。

// 为了提高源代码的可读性，我们进一步为FEEvaluation对象定义了一个别名。

  template <int dim, int fe_degree, typename number> 
  class JacobianOperator 
    : public MatrixFreeOperators:: 
        Base<dim, LinearAlgebra::distributed::Vector<number>> 
  { 
  public: 
    using value_type = number; 

    using FECellIntegrator = 
      FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number>; 

    JacobianOperator(); 

    virtual void clear() override; 

    void evaluate_newton_step( 
      const LinearAlgebra::distributed::Vector<number> &newton_step); 

    virtual void compute_diagonal() override; 

  private: 
    virtual void apply_add( 
      LinearAlgebra::distributed::Vector<number> &      dst, 
      const LinearAlgebra::distributed::Vector<number> &src) const override; 

    void 
    local_apply(const MatrixFree<dim, number> &                   data, 
                LinearAlgebra::distributed::Vector<number> &      dst, 
                const LinearAlgebra::distributed::Vector<number> &src, 
                const std::pair<unsigned int, unsigned int> &cell_range) const; 

    void local_compute_diagonal(FECellIntegrator &integrator) const; 

    Table<2, VectorizedArray<number>> nonlinear_values; 
  }; 

//  <code>JacobianOperator</code> 的构造函数只是调用基类 MatrixFreeOperators::Base, 的构造函数，而基类本身就是派生于Subscriptor类。

  template <int dim, int fe_degree, typename number> 
  JacobianOperator<dim, fe_degree, number>::JacobianOperator() 
    : MatrixFreeOperators::Base<dim, 
                                LinearAlgebra::distributed::Vector<number>>() 
  {} 

//  <code>clear()</code> 函数重置了保存非线性值的表格，并调用基类的 <code>clear()</code> 函数。

  template <int dim, int fe_degree, typename number> 
  void JacobianOperator<dim, fe_degree, number>::clear() 
  { 
    nonlinear_values.reinit(0, 0); 
    MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<number>>:: 
      clear(); 
  } 

//  @sect4{Evaluation of the old Newton step}  

// 下面的  <code>evaluate_newton_step()</code>  函数是基于  step-37  的  <code>evaluate_coefficient()</code>  函数。然而，它并不评估一个函数对象，而是评估一个代表有限元函数的向量，即雅各布系数所需的最后一个牛顿步骤。因此，我们设置了一个FEEvaluation对象，用 FEEvaluation::read_dof_values_plain() 和 FEEvaluation::evaluate() 函数评估正交点的有限元函数。我们将有限元函数的评估值直接存储在 <code>nonlinear_values</code> 表中。

//这样做会很好，在 <code>local_apply()</code> 函数中我们可以使用存储在表中的值来应用矩阵-向量乘积。然而，我们也可以在这个阶段优化雅各布系数的实现。我们可以直接评估非线性函数 <code>std::exp(newton_step[q])</code> 并将这些值存储在表中。这就跳过了在每次调用 <code>vmult()</code> 函数时对非线性的所有评估。

  template <int dim, int fe_degree, typename number> 
  void JacobianOperator<dim, fe_degree, number>::evaluate_newton_step( 
    const LinearAlgebra::distributed::Vector<number> &newton_step) 
  { 
    const unsigned int n_cells = this->data->n_cell_batches(); 
    FECellIntegrator   phi(*this->data); 

    nonlinear_values.reinit(n_cells, phi.n_q_points); 

    for (unsigned int cell = 0; cell < n_cells; ++cell) 
      { 
        phi.reinit(cell); 
        phi.read_dof_values_plain(newton_step); 
        phi.evaluate(EvaluationFlags::values); 

        for (unsigned int q = 0; q < phi.n_q_points; ++q) 
          { 
            nonlinear_values(cell, q) = std::exp(phi.get_value(q)); 
          } 
      } 
  } 

//  @sect4{Nonlinear matrix-free operator application}  

// 现在在  <code>local_apply()</code>  函数中，实际上实现了系统矩阵的单元格动作，我们可以使用存储在表  <code>nonlinear_values</code>  中的最后一个牛顿步骤的信息。这个函数的其余部分与  step-37  中的基本相同。我们设置 FEEvaluation 对象，收集并评估输入向量的值和梯度  <code>src</code>  ，根据雅各布的形式提交值和梯度，最后调用  FEEvaluation::integrate_scatter()  进行单元积分，将局部贡献分配到全局向量  <code> dst</code>  。

  template <int dim, int fe_degree, typename number> 
  void JacobianOperator<dim, fe_degree, number>::local_apply( 
    const MatrixFree<dim, number> &                   data, 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src, 
    const std::pair<unsigned int, unsigned int> &     cell_range) const 
  { 
    FECellIntegrator phi(data); 

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) 
      { 
        AssertDimension(nonlinear_values.size(0), 
                        phi.get_matrix_free().n_cell_batches()); 
        AssertDimension(nonlinear_values.size(1), phi.n_q_points); 

        phi.reinit(cell); 

        phi.gather_evaluate(src, 
                            EvaluationFlags::values | 
                              EvaluationFlags::gradients); 

        for (unsigned int q = 0; q < phi.n_q_points; ++q) 
          { 
            phi.submit_value(-nonlinear_values(cell, q) * phi.get_value(q), q); 
            phi.submit_gradient(phi.get_gradient(q), q); 
          } 

        phi.integrate_scatter(EvaluationFlags::values | 
                                EvaluationFlags::gradients, 
                              dst); 
      } 
  } 

// 接下来我们使用 MatrixFree::cell_loop() 对所有单元进行实际循环，计算单元对矩阵-向量积的贡献。

  template <int dim, int fe_degree, typename number> 
  void JacobianOperator<dim, fe_degree, number>::apply_add( 
    LinearAlgebra::distributed::Vector<number> &      dst, 
    const LinearAlgebra::distributed::Vector<number> &src) const 
  { 
    this->data->cell_loop(&JacobianOperator::local_apply, this, dst, src); 
  } 

//  @sect4{Diagonal of the JacobianOperator}  

// 用于计算对角线的内部工作函数  <code>local_compute_diagonal()</code>  与上述工作函数  <code>local_apply()</code>  类似。然而，作为主要区别，我们不从输入向量中读取数值，也不将任何局部结果分配给输出向量。相反，唯一的输入参数是使用的FEEvaluation对象。

  template <int dim, int fe_degree, typename number> 
  void JacobianOperator<dim, fe_degree, number>::local_compute_diagonal( 
    FECellIntegrator &phi) const 
  { 
    AssertDimension(nonlinear_values.size(0), 
                    phi.get_matrix_free().n_cell_batches()); 
    AssertDimension(nonlinear_values.size(1), phi.n_q_points); 

    const unsigned int cell = phi.get_current_cell_index(); 

    phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients); 

    for (unsigned int q = 0; q < phi.n_q_points; ++q) 
      { 
        phi.submit_value(-nonlinear_values(cell, q) * phi.get_value(q), q); 
        phi.submit_gradient(phi.get_gradient(q), q); 
      } 

    phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients); 
  } 

// 最后我们覆盖  MatrixFreeOperators::Base::compute_diagonal()  的基类的  <code>JacobianOperator</code>  的函数。虽然这个函数的名字表明只是计算对角线，但这个函数的作用更大。因为我们实际上只需要矩阵对角线元素的逆值，用于多网格预处理器的切比雪夫平滑器，我们计算对角线并存储逆值元素。因此我们首先初始化 <code>inverse_diagonal_entries</code>  。然后我们通过将工作函数 <code>local_compute_diagonal()</code> 传递给 MatrixFreeTools::compute_diagonal() 函数来计算对角线。最后，我们在对角线上循环，用手反转这些元素。注意，在这个循环过程中，我们捕捉受限的DOF，并手动将其设置为1。

  template <int dim, int fe_degree, typename number> 
  void JacobianOperator<dim, fe_degree, number>::compute_diagonal() 
  { 
    this->inverse_diagonal_entries.reset( 
      new DiagonalMatrix<LinearAlgebra::distributed::Vector<number>>()); 
    LinearAlgebra::distributed::Vector<number> &inverse_diagonal = 
      this->inverse_diagonal_entries->get_vector(); 
    this->data->initialize_dof_vector(inverse_diagonal); 

    MatrixFreeTools::compute_diagonal(*this->data, 
                                      inverse_diagonal, 
                                      &JacobianOperator::local_compute_diagonal, 
                                      this); 

    for (auto &diagonal_element : inverse_diagonal) 
      { 
        diagonal_element = (std::abs(diagonal_element) > 1.0e-10) ? 
                             (1.0 / diagonal_element) : 
                             1.0; 
      } 
  } 

//  @sect3{GelfandProblem class}  

// 在实现了无矩阵运算符之后，我们现在可以为<i>Gelfand problem</i>定义求解器类。这个类是基于之前所有教程程序的共同结构，特别是它是基于 step-15 ，解决的也是一个非线性问题。由于我们使用的是无矩阵框架，所以我们不再需要assemble_system函数，相反，在每次调用 <code>vmult()</code> 函数时都会重建矩阵的信息。然而，对于牛顿方案的应用，我们需要组装线性化问题的右手边并计算残差。因此，我们实现了一个额外的函数 <code>evaluate_residual()</code> ，后来我们在 <code>assemble_rhs()</code> and the <code>compute_residual()</code> 函数中调用了它。最后，这里典型的 <code>solve()</code> 函数实现了牛顿方法，而线性化系统的解是在 <code>compute_update()</code> 函数中计算的。由于MatrixFree框架将拉格朗日有限元方法的多项式程度作为一个模板参数来处理，我们也将其作为问题求解器类的模板参数来声明。

  template <int dim, int fe_degree> 
  class GelfandProblem 
  { 
  public: 
    GelfandProblem(); 

    void run(); 

  private: 
    void make_grid(); 

    void setup_system(); 

    void evaluate_residual( 
      LinearAlgebra::distributed::Vector<double> &      dst, 
      const LinearAlgebra::distributed::Vector<double> &src) const; 

    void local_evaluate_residual( 
      const MatrixFree<dim, double> &                   data, 
      LinearAlgebra::distributed::Vector<double> &      dst, 
      const LinearAlgebra::distributed::Vector<double> &src, 
      const std::pair<unsigned int, unsigned int> &     cell_range) const; 

    void assemble_rhs(); 

    double compute_residual(const double alpha); 

    void compute_update(); 

    void solve(); 

    double compute_solution_norm() const; 

    void output_results(const unsigned int cycle) const; 

// 对于并行计算，我们定义了一个  parallel::distributed::Triangulation.  由于计算域在二维是一个圆，在三维是一个球，我们除了为边界单元分配SphericalManifold外，还为内部单元的映射分配了一个TransfiniteInterpolationManifold对象，它负责处理内部单元的映射。在这个例子中，我们使用了一个等参数的有限元方法，因此使用了MappingQGeneric类。注意，我们也可以创建一个MappingQ类的实例，并在构造函数调用中设置 <code>use_mapping_q_on_all_cells</code> 标志为 <code>true</code>  。关于MappingQ和MappingQGeneric连接的进一步细节，你可以阅读这些类的详细描述。

    parallel::distributed::Triangulation<dim> triangulation; 
    const MappingQGeneric<dim>                mapping; 

// 像往常一样，我们接着定义拉格朗日有限元FE_Q和一个DoFHandler。

    FE_Q<dim>       fe; 
    DoFHandler<dim> dof_handler; 

// 对于线性化的离散系统，我们定义一个AffineConstraints对象和 <code>system_matrix</code>  ，在本例中它被表示为一个无矩阵算子。

    AffineConstraints<double> constraints; 
    using SystemMatrixType = JacobianOperator<dim, fe_degree, double>; 
    SystemMatrixType system_matrix; 

// 多级对象也是基于雅各布系数的无矩阵算子。由于我们需要用最后一个牛顿步骤来评估雅各布，所以我们也需要用最后一个牛顿步骤来评估预处理器的水平算子。因此，除了 <code>mg_matrices</code> 之外，我们还需要一个MGLevelObject来存储每一级的插值解向量。与 step-37 一样，我们对预处理程序使用浮点精度。此外，我们将MGTransferMatrixFree对象定义为一个类变量，因为我们只需要在三角形变化时设置一次，然后可以在每个牛顿步骤中再次使用它。

    MGConstrainedDoFs mg_constrained_dofs; 
    using LevelMatrixType = JacobianOperator<dim, fe_degree, float>; 
    MGLevelObject<LevelMatrixType>                           mg_matrices; 
    MGLevelObject<LinearAlgebra::distributed::Vector<float>> mg_solution; 
    MGTransferMatrixFree<dim, float>                         mg_transfer; 

// 当然，我们还需要持有  <code>solution</code>  ,  <code>newton_update</code> and the <code>system_rhs</code>  的向量。这样，我们就可以一直将上一个牛顿步存储在解的向量中，只需添加更新就可以得到下一个牛顿步。

    LinearAlgebra::distributed::Vector<double> solution; 
    LinearAlgebra::distributed::Vector<double> newton_update; 
    LinearAlgebra::distributed::Vector<double> system_rhs; 

// 最后我们有一个变量，用来表示线性求解器的迭代次数。

    unsigned int linear_iterations; 

// 对于与MPI并行运行的程序中的输出，我们使用ConditionalOStream类来避免不同MPI等级对同一数据的多次输出。

    ConditionalOStream pcout; 

// 最后，对于时间测量，我们使用一个TimerOutput对象，它在程序结束后将每个函数的耗时CPU和墙体时间打印在一个格式良好的表格中。

    TimerOutput computing_timer; 
  }; 

//  <code>GelfandProblem</code> 的构造函数初始化了类的变量。特别是，我们为 parallel::distributed::Triangulation, 设置了多级支持，将映射度设为有限元度，初始化ConditionalOStream，并告诉TimerOutput，我们只想在需求时看到墙体时间。

  template <int dim, int fe_degree> 
  GelfandProblem<dim, fe_degree>::GelfandProblem() 
    : triangulation(MPI_COMM_WORLD, 
                    Triangulation<dim>::limit_level_difference_at_vertices, 
                    parallel::distributed::Triangulation< 
                      dim>::construct_multigrid_hierarchy) 
    , mapping(fe_degree) 
    , fe(fe_degree) 
    , dof_handler(triangulation) 
    , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
    , computing_timer(MPI_COMM_WORLD, 
                      pcout, 
                      TimerOutput::never, 
                      TimerOutput::wall_times) 
  {} 

//  @sect4{GelfandProblem::make_grid}  

// 作为计算域，我们使用 <code>dim</code>  -维的单位球。我们按照TransfiniteInterpolationManifold类的说明，也为边界指定了一个SphericalManifold。最后，我们将初始网格细化为3

// -  <code>dim</code> 次全局。

  template <int dim, int fe_degree> 
  void GelfandProblem<dim, fe_degree>::make_grid() 
  { 
    TimerOutput::Scope t(computing_timer, "make grid"); 

    SphericalManifold<dim>                boundary_manifold; 
    TransfiniteInterpolationManifold<dim> inner_manifold; 

    GridGenerator::hyper_ball(triangulation); 

 
    triangulation.set_all_manifold_ids_on_boundary(0); 

    triangulation.set_manifold(0, boundary_manifold); 

    inner_manifold.initialize(triangulation); 
    triangulation.set_manifold(1, inner_manifold); 

    triangulation.refine_global(3 - dim); 
  } 

//  @sect4{GelfandProblem::setup_system}  

//  <code>setup_system()</code> 函数与  step-37  中的函数基本相同。唯一的区别显然是时间测量只有一个 TimerOutput::Scope ，而不是单独测量每个部分，更重要的是对前一个牛顿步骤的内插解向量的MGLevelObject的初始化。另一个重要的变化是MGTransferMatrixFree对象的设置，我们可以在每个牛顿步骤中重复使用它，因为 <code>triangulation</code> 不会被改变。

// 注意我们如何在 <code>JacobianOperator</code> 和多网格预处理程序中两次使用同一个MatrixFree对象。

  template <int dim, int fe_degree> 
  void GelfandProblem<dim, fe_degree>::setup_system() 
  { 
    TimerOutput::Scope t(computing_timer, "setup system"); 

    system_matrix.clear(); 
    mg_matrices.clear_elements(); 

    dof_handler.distribute_dofs(fe); 
    dof_handler.distribute_mg_dofs(); 

    IndexSet locally_relevant_dofs; 
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs); 

    constraints.clear(); 
    constraints.reinit(locally_relevant_dofs); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             0, 
                                             Functions::ZeroFunction<dim>(), 
                                             constraints); 
    constraints.close(); 

    { 
      typename MatrixFree<dim, double>::AdditionalData additional_data; 
      additional_data.tasks_parallel_scheme = 
        MatrixFree<dim, double>::AdditionalData::partition_color; 
      additional_data.mapping_update_flags = 
        (update_values | update_gradients | update_JxW_values | 
         update_quadrature_points); 
      auto system_mf_storage = std::make_shared<MatrixFree<dim, double>>(); 
      system_mf_storage->reinit(mapping, 
                                dof_handler, 
                                constraints, 
                                QGauss<1>(fe.degree + 1), 
                                additional_data); 

      system_matrix.initialize(system_mf_storage); 
    } 

    system_matrix.initialize_dof_vector(solution); 
    system_matrix.initialize_dof_vector(newton_update); 
    system_matrix.initialize_dof_vector(system_rhs); 

    const unsigned int nlevels = triangulation.n_global_levels(); 
    mg_matrices.resize(0, nlevels - 1); 
    mg_solution.resize(0, nlevels - 1); 

    std::set<types::boundary_id> dirichlet_boundary; 
    dirichlet_boundary.insert(0); 
    mg_constrained_dofs.initialize(dof_handler); 
    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, 
                                                       dirichlet_boundary); 

    mg_transfer.initialize_constraints(mg_constrained_dofs); 
    mg_transfer.build(dof_handler); 

    for (unsigned int level = 0; level < nlevels; ++level) 
      { 
        IndexSet relevant_dofs; 
        DoFTools::extract_locally_relevant_level_dofs(dof_handler, 
                                                      level, 
                                                      relevant_dofs); 

        AffineConstraints<double> level_constraints; 
        level_constraints.reinit(relevant_dofs); 
        level_constraints.add_lines( 
          mg_constrained_dofs.get_boundary_indices(level)); 
        level_constraints.close(); 

        typename MatrixFree<dim, float>::AdditionalData additional_data; 
        additional_data.tasks_parallel_scheme = 
          MatrixFree<dim, float>::AdditionalData::partition_color; 
        additional_data.mapping_update_flags = 
          (update_values | update_gradients | update_JxW_values | 
           update_quadrature_points); 
        additional_data.mg_level = level; 
        auto mg_mf_storage_level = std::make_shared<MatrixFree<dim, float>>(); 
        mg_mf_storage_level->reinit(mapping, 
                                    dof_handler, 
                                    level_constraints, 
                                    QGauss<1>(fe.degree + 1), 
                                    additional_data); 

        mg_matrices[level].initialize(mg_mf_storage_level, 
                                      mg_constrained_dofs, 
                                      level); 
        mg_matrices[level].initialize_dof_vector(mg_solution[level]); 
      } 
  } 

//  @sect4{GelfandProblem::evaluate_residual}  

// 接下来我们实现一个函数，该函数对给定的输入向量评估非线性离散残差（  $\texttt{dst} = F(\texttt{src})$  ）。这个函数随后被用于组装线性化系统的右手边，随后用于计算下一个牛顿步骤的残差，以检查我们是否已经达到了误差容忍度。由于这个函数不应该影响任何类别的变量，我们把它定义为一个常数函数。在内部，我们通过FEEvaluation类和类似于 MatrixFree::cell_loop(), 的 <code>apply_add()</code> function of the <code>JacobianOperator</code> 来利用快速有限元评估。

// 首先我们创建一个指向MatrixFree对象的指针，它被存储在  <code>system_matrix</code>  中。然后，我们将用于残差的单元评估的工作函数  <code>local_evaluate_residual()</code>  以及输入和输出向量传递给  MatrixFree::cell_loop().  此外，我们在循环中启用输出向量的清零，这比之前单独调用<code>dst = 0.0</code>更有效率。

// 注意，使用这种方法，我们不必关心MPI相关的数据交换，因为所有的记账工作都是由  MatrixFree::cell_loop().  完成的。
  template <int dim, int fe_degree> 
  void GelfandProblem<dim, fe_degree>::evaluate_residual( 
    LinearAlgebra::distributed::Vector<double> &      dst, 
    const LinearAlgebra::distributed::Vector<double> &src) const 
  { 
    auto matrix_free = system_matrix.get_matrix_free(); 

    matrix_free->cell_loop( 
      &GelfandProblem::local_evaluate_residual, this, dst, src, true); 
  } 

//  @sect4{GelfandProblem::local_evaluate_residual}  

// 这是用于评估残差的内部工作函数。本质上它与  <code>JacobianOperator</code>  的  <code>local_apply()</code>  函数具有相同的结构，在给定的单元格集合  <code>cell_range</code>  上对输入向量  <code>src</code>  进行残差评估。与上述 <code>local_apply()</code> 函数不同的是，我们将 FEEvaluation::gather_evaluate() 函数分成 FEEvaluation::read_dof_values_plain() 和 FEEvaluation::evaluate(), ，因为输入向量可能有受限的DOF。

  template <int dim, int fe_degree> 
  void GelfandProblem<dim, fe_degree>::local_evaluate_residual( 
    const MatrixFree<dim, double> &                   data, 
    LinearAlgebra::distributed::Vector<double> &      dst, 
    const LinearAlgebra::distributed::Vector<double> &src, 
    const std::pair<unsigned int, unsigned int> &     cell_range) const 
  { 
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, double> phi(data); 

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) 
      { 
        phi.reinit(cell); 

        phi.read_dof_values_plain(src); 
        phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients); 

        for (unsigned int q = 0; q < phi.n_q_points; ++q) 
          { 
            phi.submit_value(-std::exp(phi.get_value(q)), q); 
            phi.submit_gradient(phi.get_gradient(q), q); 
          } 

        phi.integrate_scatter(EvaluationFlags::values | 
                                EvaluationFlags::gradients, 
                              dst); 
      } 
  } 

//  @sect4{GelfandProblem::assemble_rhs}  

// 使用上述函数 <code>evaluate_residual()</code> 来评估非线性残差，组装线性化系统的右手边现在变得非常容易。我们只需调用 <code>evaluate_residual()</code> 函数并将结果乘以减一。

// 经验表明，使用FEEvaluation类要比使用FEValues和co的经典实现快得多。

  template <int dim, int fe_degree> 
  void GelfandProblem<dim, fe_degree>::assemble_rhs() 
  { 
    TimerOutput::Scope t(computing_timer, "assemble right hand side"); 

    evaluate_residual(system_rhs, solution); 

    system_rhs *= -1.0; 
  } 

//  @sect4{GelfandProblem::compute_residual}  

// 根据 step-15 ，下面的函数在 $u_h^n + \alpha s_h^n$ 函数的帮助下计算出解的非线性残差的规范。如果我们使用牛顿方法的自适应版本，牛顿步长 $\alpha$ 就变得很重要。例如，我们将计算不同步长的残差并比较残差。然而，对于我们的问题，使用 $\alpha=1$ 的完整牛顿步长是我们能做的最好的。如果我们没有好的初始值，牛顿方法的自适应版本就变得有趣了。请注意，在理论上，牛顿方法是以二次方顺序收敛的，但只有当我们有一个合适的初始值时才会收敛。对于不合适的初始值，牛顿方法甚至在二次方程下也会发散。一个常见的方法是使用阻尼版本 $\alpha<1$ ，直到牛顿步骤足够好，可以进行完整的牛顿步骤。这在  step-15  中也有讨论。

  template <int dim, int fe_degree> 
  double GelfandProblem<dim, fe_degree>::compute_residual(const double alpha) 
  { 
    TimerOutput::Scope t(computing_timer, "compute residual"); 

    LinearAlgebra::distributed::Vector<double> residual; 
    LinearAlgebra::distributed::Vector<double> evaluation_point; 

    system_matrix.initialize_dof_vector(residual); 
    system_matrix.initialize_dof_vector(evaluation_point); 

    evaluation_point = solution; 
    if (alpha > 1e-12) 
      { 
        evaluation_point.add(alpha, newton_update); 
      } 

    evaluate_residual(residual, evaluation_point); 

    return residual.l2_norm(); 
  } 

//  @sect4{GelfandProblem::compute_update}  

// 为了计算每个牛顿步骤中的牛顿更新，我们用CG算法和一个几何多网格预处理程序来解决线性系统。为此，我们首先像在  step-37  中那样，用切比雪夫平滑器设置PreconditionMG对象。

  template <int dim, int fe_degree> 
  void GelfandProblem<dim, fe_degree>::compute_update() 
  { 
    TimerOutput::Scope t(computing_timer, "compute update"); 

// 我们记得，雅各布系数取决于存储在解决方案向量中的最后一个牛顿步骤。所以我们更新牛顿步骤的鬼魂值，并将其传递给 <code>JacobianOperator</code> 来存储信息。

    solution.update_ghost_values(); 

    system_matrix.evaluate_newton_step(solution); 

// 接下来我们还要将最后一个牛顿步骤传递给多级运算符。因此，我们需要将牛顿步骤插值到三角形的所有层面。这是用 MGTransferMatrixFree::interpolate_to_mg(). 来完成的。
    mg_transfer.interpolate_to_mg(dof_handler, mg_solution, solution); 

// 现在我们可以设置预处理程序了。我们定义平滑器并将牛顿步的内插向量传递给多级运算器。

    using SmootherType = 
      PreconditionChebyshev<LevelMatrixType, 
                            LinearAlgebra::distributed::Vector<float>>; 
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
            smoother_data[level].degree              = 4; 
            smoother_data[level].eig_cg_n_iterations = 10; 
          } 
        else 
          { 
            smoother_data[0].smoothing_range = 1e-3; 
            smoother_data[0].degree          = numbers::invalid_unsigned_int; 
            smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m(); 
          } 

        mg_matrices[level].evaluate_newton_step(mg_solution[level]); 
        mg_matrices[level].compute_diagonal(); 

        smoother_data[level].preconditioner = 
          mg_matrices[level].get_matrix_diagonal_inverse(); 
      } 
    mg_smoother.initialize(mg_matrices, smoother_data); 

    MGCoarseGridApplySmoother<LinearAlgebra::distributed::Vector<float>> 
      mg_coarse; 
    mg_coarse.initialize(mg_smoother); 

    mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_matrix( 
      mg_matrices); 

    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>> 
      mg_interface_matrices; 
    mg_interface_matrices.resize(0, triangulation.n_global_levels() - 1); 
    for (unsigned int level = 0; level < triangulation.n_global_levels(); 
         ++level) 
      { 
        mg_interface_matrices[level].initialize(mg_matrices[level]); 
      } 
    mg::Matrix<LinearAlgebra::distributed::Vector<float>> mg_interface( 
      mg_interface_matrices); 

    Multigrid<LinearAlgebra::distributed::Vector<float>> mg( 
      mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother); 
    mg.set_edge_matrices(mg_interface, mg_interface); 

    PreconditionMG<dim, 
                   LinearAlgebra::distributed::Vector<float>, 
                   MGTransferMatrixFree<dim, float>> 
      preconditioner(dof_handler, mg, mg_transfer); 

// 最后我们设置了SolverControl和SolverCG来解决当前牛顿更新的线性化问题。实现SolverCG或SolverGMRES的一个重要事实是，持有线性系统解决方案的向量（这里是 <code>newton_update</code>  ）可以用来传递一个起始值。为了使迭代求解器总是以零向量开始，我们在调用 SolverCG::solve(). 之前明确地重置了 <code>newton_update</code> ，然后我们分配了存储在 <code>constraints</code> 中的Dirichlet边界条件，并为以后的输出存储了迭代的步数。

    SolverControl solver_control(100, 1.e-12); 
    SolverCG<LinearAlgebra::distributed::Vector<double>> cg(solver_control); 

    newton_update = 0.0; 

    cg.solve(system_matrix, newton_update, system_rhs, preconditioner); 

    constraints.distribute(newton_update); 

    linear_iterations = solver_control.last_step(); 

// 然后，为了记账，我们将幽灵值清零。

    solution.zero_out_ghost_values(); 
  } 

//  @sect4{GelfandProblem::solve}  

// 现在我们实现非线性问题的实际牛顿求解器。

  template <int dim, int fe_degree> 
  void GelfandProblem<dim, fe_degree>::solve() 
  { 
    TimerOutput::Scope t(computing_timer, "solve"); 

// 我们定义了牛顿步骤的最大数量和收敛标准的公差。通常情况下，如果有好的起始值，牛顿方法在三到六步内就能收敛，所以最大的十步应该是完全足够的。作为公差，我们使用 $\|F(u^n_h)\|<\text{TOL}_f = 10^{-12}$ 作为残差的规范， $\|s_h^n\| < \text{TOL}_x = 10^{-10}$ 作为牛顿更新的规范。这似乎有点过头了，但我们将看到，对于我们的例子，我们将在几步之后达到这些公差。

    const unsigned int itmax = 10; 
    const double       TOLf  = 1e-12; 
    const double       TOLx  = 1e-10; 

    Timer solver_timer; 
    solver_timer.start(); 

// 现在我们开始实际的牛顿迭代。

    for (unsigned int newton_step = 1; newton_step <= itmax; ++newton_step) 
      { 

// 我们将线性化问题的右侧集合起来，计算牛顿更新。

        assemble_rhs(); 
        compute_update(); 

// 然后，我们计算误差，即牛顿更新的规范和残差。注意，在这一点上，我们可以通过改变compute_residual函数的输入参数 $\alpha$ 来加入牛顿方法的步长控制。然而，在这里我们只是使用 $\alpha$ 等于1来进行普通的牛顿迭代。

        const double ERRx = newton_update.l2_norm(); 
        const double ERRf = compute_residual(1.0); 

// 接下来我们通过将牛顿更新添加到当前的牛顿步骤中来推进牛顿步骤。

        solution.add(1.0, newton_update); 

// 一个简短的输出将告知我们当前的牛顿步数。

        pcout << "   Nstep " << newton_step << ", errf = " << ERRf 
              << ", errx = " << ERRx << ", it = " << linear_iterations 
              << std::endl; 

// 在每个牛顿步骤之后，我们检查收敛标准。如果其中至少有一个得到满足，我们就完成了，并结束循环。如果我们在牛顿迭代的最大数量之后还没有找到一个满意的解决方案，我们就会通知用户这个缺点。

        if (ERRf < TOLf || ERRx < TOLx) 
          { 
            solver_timer.stop(); 

            pcout << "Convergence step " << newton_step << " value " << ERRf 
                  << " (used wall time: " << solver_timer.wall_time() << " s)" 
                  << std::endl; 

            break; 
          } 
        else if (newton_step == itmax) 
          { 
            solver_timer.stop(); 
            pcout << "WARNING: No convergence of Newton's method after " 
                  << newton_step << " steps." << std::endl; 

            break; 
          } 
      } 
  } 

//  @sect4{GelfandProblem::compute_solution_norm}  

// 解的H1-seminorm的计算可以用与 step-59 相同的方法进行。我们更新幽灵值并使用函数  VectorTools::integrate_difference().  最后我们收集所有MPI行列的所有计算，并返回规范。

  template <int dim, int fe_degree> 
  double GelfandProblem<dim, fe_degree>::compute_solution_norm() const 
  { 
    solution.update_ghost_values(); 

    Vector<float> norm_per_cell(triangulation.n_active_cells()); 

    VectorTools::integrate_difference(mapping, 
                                      dof_handler, 
                                      solution, 
                                      Functions::ZeroFunction<dim>(), 
                                      norm_per_cell, 
                                      QGauss<dim>(fe.degree + 2), 
                                      VectorTools::H1_seminorm); 

    solution.zero_out_ghost_values(); 

    return VectorTools::compute_global_error(triangulation, 
                                             norm_per_cell, 
                                             VectorTools::H1_seminorm); 
  } 

//  @sect4{GelfandProblem::output_results}  

// 我们通过调用  DataOut::write_vtu_with_pvtu_record()  函数，以与  step-37  中相同的方式，一次性生成 vtu 格式的图形输出文件和 pvtu 主文件。此外，与  step-40  一样，我们查询每个单元的  types::subdomain_id  并将三角形在MPI行列中的分布写进输出文件。最后，我们通过调用 DataOut::build_patches(). 生成解决方案的补丁。然而，由于我们的计算域有一个弯曲的边界，我们另外传递 <code>mapping</code> 和有限元度作为细分的数量。但这仍然不足以正确表示解决方案，例如在ParaView中，因为我们将TransfiniteInterpolationManifold附在内部单元上，这导致内部的单元是弯曲的。因此，我们将 DataOut::curved_inner_cells 选项作为第三个参数，这样，内部单元也会使用相应的流形描述来构建补丁。

// 注意，我们可以用标志 DataOutBase::VtkFlags::write_higher_order_cells. 来处理高阶元素，但是由于对ParaView以前版本的兼容性有限，而且VisIt也不支持，所以我们把这个选项留给未来的版本。

  template <int dim, int fe_degree> 
  void 
  GelfandProblem<dim, fe_degree>::output_results(const unsigned int cycle) const 
  { 
    if (triangulation.n_global_active_cells() > 1e6) 
      return; 

    solution.update_ghost_values(); 

    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "solution"); 

    Vector<float> subdomain(triangulation.n_active_cells()); 
    for (unsigned int i = 0; i < subdomain.size(); ++i) 
      { 
        subdomain(i) = triangulation.locally_owned_subdomain(); 
      } 
    data_out.add_data_vector(subdomain, "subdomain"); 

    data_out.build_patches(mapping, 
                           fe.degree, 
                           DataOut<dim>::curved_inner_cells); 

    DataOutBase::VtkFlags flags; 
    flags.compression_level = DataOutBase::VtkFlags::best_speed; 
    data_out.set_flags(flags); 
    data_out.write_vtu_with_pvtu_record( 
      "./", "solution_" + std::to_string(dim) + "d", cycle, MPI_COMM_WORLD, 3); 

    solution.zero_out_ghost_values(); 
  } 

//  @sect4{GelfandProblem::run}  

// <i>Gelfand
 problem</i>的求解器类的最后一个缺失的函数是运行函数。在开始的时候，我们打印关于系统规格和我们使用的有限元空间的信息。该问题在一个连续细化的网格上被多次求解。

  template <int dim, int fe_degree> 
  void GelfandProblem<dim, fe_degree>::run() 
  { 
    { 
      const unsigned int n_ranks = 
        Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); 
      const unsigned int n_vect_doubles = VectorizedArray<double>::size(); 
      const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles; 

      std::string DAT_header = "START DATE: " + Utilities::System::get_date() + 
                               ", TIME: " + Utilities::System::get_time(); 
      std::string MPI_header = "Running with " + std::to_string(n_ranks) + 
                               " MPI process" + (n_ranks > 1 ? "es" : ""); 
      std::string VEC_header = 
        "Vectorization over " + std::to_string(n_vect_doubles) + 
        " doubles = " + std::to_string(n_vect_bits) + " bits (" + 
        Utilities::System::get_current_vectorization_level() + 
        "), VECTORIZATION_LEVEL=" + 
        std::to_string(DEAL_II_COMPILER_VECTORIZATION_LEVEL); 
      std::string SOL_header = "Finite element space: " + fe.get_name(); 

      pcout << std::string(80, '=') << std::endl; 
      pcout << DAT_header << std::endl; 
      pcout << std::string(80, '-') << std::endl; 

      pcout << MPI_header << std::endl; 
      pcout << VEC_header << std::endl; 
      pcout << SOL_header << std::endl; 

      pcout << std::string(80, '=') << std::endl; 
    } 

    for (unsigned int cycle = 0; cycle < 9 - dim; ++cycle) 
      { 
        pcout << std::string(80, '-') << std::endl; 
        pcout << "Cycle " << cycle << std::endl; 
        pcout << std::string(80, '-') << std::endl; 

// 实际解决问题的第一项任务是生成或完善三角图。

        if (cycle == 0) 
          { 
            make_grid(); 
          } 
        else 
          { 
            triangulation.refine_global(1); 
          } 

// 现在我们建立了系统并解决这个问题。这些步骤都伴随着时间测量和文本输出。

        Timer timer; 

        pcout << "Set up system..." << std::endl; 
        setup_system(); 

        pcout << "   Triangulation: " << triangulation.n_global_active_cells() 
              << " cells" << std::endl; 
        pcout << "   DoFHandler:    " << dof_handler.n_dofs() << " DoFs" 
              << std::endl; 
        pcout << std::endl; 

        pcout << "Solve using Newton's method..." << std::endl; 
        solve(); 
        pcout << std::endl; 

        timer.stop(); 
        pcout << "Time for setup+solve (CPU/Wall) " << timer.cpu_time() << "/" 
              << timer.wall_time() << " s" << std::endl; 
        pcout << std::endl; 

// 在问题被解决后，我们计算出解决方案的法线，并生成图形输出文件。

        pcout << "Output results..." << std::endl; 
        const double norm = compute_solution_norm(); 
        output_results(cycle); 

        pcout << "  H1 seminorm: " << norm << std::endl; 
        pcout << std::endl; 

// 最后在每个周期后，我们打印计时信息。

        computing_timer.print_summary(); 
        computing_timer.reset(); 
      } 
  } 
} // namespace Step66 

//  @sect3{The <code>main</code> function}  

// 作为使用MPI并行运行的典型程序，我们设置了MPI框架，并通过限制线程数为1来禁用共享内存并行化。最后，为了运行<i>Gelfand problem</i>的求解器，我们创建一个 <code>GelfandProblem</code> 类的对象并调用运行函数。例如，我们用四阶拉格朗日有限元在二维和三维中各解决一次问题。

int main(int argc, char *argv[]) 
{ 
  try 
    { 
      using namespace Step66; 

      Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1); 

      { 
        GelfandProblem<2, 4> gelfand_problem; 
        gelfand_problem.run(); 
      } 

      { 
        GelfandProblem<3, 4> gelfand_problem; 
        gelfand_problem.run(); 
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
    } 

  return 0; 
} 


