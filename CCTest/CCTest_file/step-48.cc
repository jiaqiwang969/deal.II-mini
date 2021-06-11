

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2011 - 2021 by the deal.II authors 
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
 * Author: Katharina Kormann, Martin Kronbichler, Uppsala University, 2011-2012 
 */ 



// deal.II库中的必要文件。

#include <deal.II/base/logstream.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/conditional_ostream.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/distributed/tria.h> 

// 这包括用于有效实现无矩阵方法的数据结构。

#include <deal.II/lac/la_parallel_vector.h> 
#include <deal.II/matrix_free/matrix_free.h> 
#include <deal.II/matrix_free/fe_evaluation.h> 

#include <fstream> 
#include <iostream> 
#include <iomanip> 

namespace Step48 
{ 
  using namespace dealii; 

// 我们首先定义了两个全局变量，以便在一个地方收集所有需要改变的参数。一个是尺寸，一个是有限元度。维度在主函数中是作为实际类的模板参数使用的（就像所有其他deal.II程序一样），而有限元的度数则更为关键，因为它是作为模板参数传递给Sine-Gordon算子的实现。因此，它需要成为一个编译时常数。

  const unsigned int dimension = 2; 
  const unsigned int fe_degree = 4; 
// @sect3{SineGordonOperation}  

//  <code>SineGordonOperation</code> 类实现了每个时间步骤中需要的基于单元的操作。这个非线性操作可以在 <code>MatrixFree</code> 类的基础上直接实现，与线性操作在这个实现的有限元算子应用中的处理方式相同。我们对该类应用了两个模板参数，一个是尺寸，一个是有限元的程度。这与deal.II中的其他函数不同，其中只有维度是模板参数。这对于为 @p FEEvaluation 中的内循环提供关于循环长度等的信息是必要的，这对于效率是至关重要的。另一方面，这使得将度数作为一个运行时参数来实现更具挑战性。

  template <int dim, int fe_degree> 
  class SineGordonOperation 
  { 
  public: 
    SineGordonOperation(const MatrixFree<dim, double> &data_in, 
                        const double                   time_step); 

    void apply(LinearAlgebra::distributed::Vector<double> &dst, 
               const std::vector<LinearAlgebra::distributed::Vector<double> *> 
                 &src) const; 

  private: 
    const MatrixFree<dim, double> &            data; 
    const VectorizedArray<double>              delta_t_sqr; 
    LinearAlgebra::distributed::Vector<double> inv_mass_matrix; 

    void local_apply( 
      const MatrixFree<dim, double> &                                  data, 
      LinearAlgebra::distributed::Vector<double> &                     dst, 
      const std::vector<LinearAlgebra::distributed::Vector<double> *> &src, 
      const std::pair<unsigned int, unsigned int> &cell_range) const; 
  }; 

//  @sect4{SineGordonOperation::SineGordonOperation}  

// 这是SineGordonOperation类的构造函数。它接收一个对MatrixFree的引用，该引用持有问题信息和时间步长作为输入参数。初始化程序设置了质量矩阵。由于我们使用Gauss-Lobatto元素，质量矩阵是一个对角矩阵，可以存储为一个矢量。利用FEEvaluation提供的数据结构，质量矩阵对角线的计算很容易实现。只要在所有的单元格批次上循环，即由于SIMD矢量化的单元格集合，并通过使用 <code>integrate</code> 函数与 @p true 参数在数值的槽上对所有正交点上常一的函数进行积分。最后，我们将对角线条目进行反转，以便在每个时间步长中直接获得反质量矩阵。

  template <int dim, int fe_degree> 
  SineGordonOperation<dim, fe_degree>::SineGordonOperation( 
    const MatrixFree<dim, double> &data_in, 
    const double                   time_step) 
    : data(data_in) 
    , delta_t_sqr(make_vectorized_array(time_step * time_step)) 
  { 
    data.initialize_dof_vector(inv_mass_matrix); 

    FEEvaluation<dim, fe_degree> fe_eval(data); 
    const unsigned int           n_q_points = fe_eval.n_q_points; 

    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell) 
      { 
        fe_eval.reinit(cell); 
        for (unsigned int q = 0; q < n_q_points; ++q) 
          fe_eval.submit_value(make_vectorized_array(1.), q); 
        fe_eval.integrate(EvaluationFlags::values); 
        fe_eval.distribute_local_to_global(inv_mass_matrix); 
      } 

    inv_mass_matrix.compress(VectorOperation::add); 
    for (unsigned int k = 0; k < inv_mass_matrix.locally_owned_size(); ++k) 
      if (inv_mass_matrix.local_element(k) > 1e-15) 
        inv_mass_matrix.local_element(k) = 
          1. / inv_mass_matrix.local_element(k); 
      else 
        inv_mass_matrix.local_element(k) = 1; 
  } 

//  @sect4{SineGordonOperation::local_apply}  

// 这个算子实现了程序的核心操作，即对正弦-戈登问题的非线性算子进行单元范围的积分。其实现是基于  step-37  中的FEEvaluation类。由于Gauss-Lobatto元素的特殊结构，某些操作变得更加简单，特别是正交点上的形状函数值的评估，这只是单元自由度值的注入。MatrixFree类在初始化时检测了正交点上有限元的可能结构，然后由FEEvaluation自动用于选择最合适的数值核。

// 我们要为时间步进例程评估的非线性函数包括当前时间的函数值 @p current 以及前一个时间步进的值 @p old. 这两个值都在源向量集合 @p src, 中传递给运算器，该集合只是一个指向实际解向量的 <tt>std::vector</tt> 指针。这种将多个源向量收集到一起的结构是必要的，因为 @p MatrixFree 中的单元格循环正好需要一个源向量和一个目的向量，即使我们碰巧使用了很多向量，比如本例中的两个。请注意，单元格循环接受任何有效的输入和输出类，这不仅包括向量，还包括一般的数据类型。 然而，只有在遇到收集这些向量的 LinearAlgebra::distributed::Vector<Number> 或 <tt>std::vector</tt> 时，它才会在循环的开始和结束时调用由于MPI而交换幽灵数据的函数。在单元格的循环中，我们首先要读入与本地值相关的向量中的值。 然后，我们评估当前求解向量的值和梯度以及正交点的旧向量的值。接下来，我们在正交点的循环中结合方案中的条款。最后，我们将结果与测试函数进行积分，并将结果累积到全局解向量 @p  dst。

  template <int dim, int fe_degree> 
  void SineGordonOperation<dim, fe_degree>::local_apply( 
    const MatrixFree<dim> &                                          data, 
    LinearAlgebra::distributed::Vector<double> &                     dst, 
    const std::vector<LinearAlgebra::distributed::Vector<double> *> &src, 
    const std::pair<unsigned int, unsigned int> &cell_range) const 
  { 
    AssertDimension(src.size(), 2); 
    FEEvaluation<dim, fe_degree> current(data), old(data); 
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) 
      { 
        current.reinit(cell); 
        old.reinit(cell); 

        current.read_dof_values(*src[0]); 
        old.read_dof_values(*src[1]); 

        current.evaluate(EvaluationFlags::values | EvaluationFlags::gradients); 
        old.evaluate(EvaluationFlags::values); 

        for (unsigned int q = 0; q < current.n_q_points; ++q) 
          { 
            const VectorizedArray<double> current_value = current.get_value(q); 
            const VectorizedArray<double> old_value     = old.get_value(q); 

            current.submit_value(2. * current_value - old_value - 
                                   delta_t_sqr * std::sin(current_value), 
                                 q); 
            current.submit_gradient(-delta_t_sqr * current.get_gradient(q), q); 
          } 

        current.integrate(EvaluationFlags::values | EvaluationFlags::gradients); 
        current.distribute_local_to_global(dst); 
      } 
  } 

//  @sect4{SineGordonOperation::apply}  

// 该函数根据单元本地策略执行时间步进例程。请注意，在添加当前时间步长的积分贡献之前，我们需要将目标向量设置为零（通过 FEEvaluation::distribute_local_to_global() 调用）。在本教程中，我们通过传递给 MatrixFree::cell_loop. 的第五个`true`参数让单元格循环进行归零操作。 循环可以将归零操作安排在更接近对支持的向量项的操作，从而可能提高数据的定位性（首先被归零的向量项后来在`distribute_local_to_global()`调用中重新使用）。单元循环的结构是在单元有限元运算器类中实现的。在每个单元上，它应用定义为类 <code>local_apply()</code> 方法的例程  <code>SineGordonOperation</code>, i.e., <code>this</code>  。我们也可以提供一个具有相同签名的、不属于类的函数。最后，积分的结果要乘以质量矩阵的逆值。

  template <int dim, int fe_degree> 
  void SineGordonOperation<dim, fe_degree>::apply( 
    LinearAlgebra::distributed::Vector<double> &                     dst, 
    const std::vector<LinearAlgebra::distributed::Vector<double> *> &src) const 
  { 
    data.cell_loop( 
      &SineGordonOperation<dim, fe_degree>::local_apply, this, dst, src, true); 
    dst.scale(inv_mass_matrix); 
  } 

//  @sect3{Equation data}  

// 我们定义了一个随时间变化的函数，作为初始值使用。通过改变起始时间，可以得到不同的解决方案。这个函数取自 step-25 ，将代表一维中所有时间的分析解，但在这里只是用来设置一些感兴趣的起始解。在  step-25  中给出了可以测试该程序收敛性的更详细的选择。

  template <int dim> 
  class InitialCondition : public Function<dim> 
  { 
  public: 
    InitialCondition(const unsigned int n_components = 1, 
                     const double       time         = 0.) 
      : Function<dim>(n_components, time) 
    {} 
    virtual double value(const Point<dim> &p, 
                         const unsigned int /*component*/) const override 
    { 
      double t = this->get_time(); 

      const double m  = 0.5; 
      const double c1 = 0.; 
      const double c2 = 0.; 
      const double factor = 
        (m / std::sqrt(1. - m * m) * std::sin(std::sqrt(1. - m * m) * t + c2)); 
      double result = 1.; 
      for (unsigned int d = 0; d < dim; ++d) 
        result *= -4. * std::atan(factor / std::cosh(m * p[d] + c1)); 
      return result; 
    } 
  }; 

//  @sect3{SineGordonProblem class}  

// 这是在  step-25  中的类基础上的主类。 然而，我们用MatrixFree类代替了SparseMatrix<double>类来存储几何数据。另外，我们在这个例子中使用了一个分布式三角形。

  template <int dim> 
  class SineGordonProblem 
  { 
  public: 
    SineGordonProblem(); 
    void run(); 

  private: 
    ConditionalOStream pcout; 

    void make_grid_and_dofs(); 
    void output_results(const unsigned int timestep_number); 

#ifdef DEAL_II_WITH_P4EST 
    parallel::distributed::Triangulation<dim> triangulation; 
#else 
    Triangulation<dim> triangulation; 
#endif 
    FE_Q<dim>       fe; 
    DoFHandler<dim> dof_handler; 

    MappingQ1<dim> mapping; 

    AffineConstraints<double> constraints; 
    IndexSet                  locally_relevant_dofs; 

    MatrixFree<dim, double> matrix_free_data; 

    LinearAlgebra::distributed::Vector<double> solution, old_solution, 
      old_old_solution; 

    const unsigned int n_global_refinements; 
    double             time, time_step; 
    const double       final_time; 
    const double       cfl_number; 
    const unsigned int output_timestep_skip; 
  }; 
// @sect4{SineGordonProblem::SineGordonProblem}  

// 这是SineGordonProblem类的构造函数。时间间隔和时间步长在此定义。此外，我们使用在程序顶部定义的有限元的程度来初始化一个基于Gauss-Lobatto支持点的FE_Q有限元。这些点很方便，因为与同阶的QGauss-Lobatto正交规则相结合，它们可以得到一个对角线质量矩阵，而不会太影响精度（注意，虽然积分是不精确的），也可以参见介绍中的讨论。请注意，FE_Q默认选择Gauss-Lobatto结点，因为它们相对于等距结点有更好的条件。为了使事情更加明确，我们还是要说明节点的选择。

  template <int dim> 
  SineGordonProblem<dim>::SineGordonProblem() 
    : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
    , 
#ifdef DEAL_II_WITH_P4EST 
    triangulation(MPI_COMM_WORLD) 
    , 
#endif 
    fe(QGaussLobatto<1>(fe_degree + 1)) 
    , dof_handler(triangulation) 
    , n_global_refinements(10 - 2 * dim) 
    , time(-10) 
    , time_step(10.) 
    , final_time(10.) 
    , cfl_number(.1 / fe_degree) 
    , output_timestep_skip(200) 
  {} 
// @sect4{SineGordonProblem::make_grid_and_dofs}  

// 和 step-25 一样，这个函数在 <code>dim</code> 维度上设置了一个范围为 $[-15,15]$ 的立方体网格。我们在域的中心更多的细化网格，因为解决方案都集中在那里。我们首先细化所有中心在半径为11的单元，然后再细化一次半径为6的单元。 这种简单的临时细化可以通过在时间步进过程中使用误差估计器来适应网格，并使用 parallel::distributed::SolutionTransfer 将解决方案转移到新的网格中来完成。

  template <int dim> 
  void SineGordonProblem<dim>::make_grid_and_dofs() 
  { 
    GridGenerator::hyper_cube(triangulation, -15, 15); 
    triangulation.refine_global(n_global_refinements); 
    { 
      typename Triangulation<dim>::active_cell_iterator 
        cell     = triangulation.begin_active(), 
        end_cell = triangulation.end(); 
      for (; cell != end_cell; ++cell) 
        if (cell->is_locally_owned()) 
          if (cell->center().norm() < 11) 
            cell->set_refine_flag(); 
      triangulation.execute_coarsening_and_refinement(); 

      cell     = triangulation.begin_active(); 
      end_cell = triangulation.end(); 
      for (; cell != end_cell; ++cell) 
        if (cell->is_locally_owned()) 
          if (cell->center().norm() < 6) 
            cell->set_refine_flag(); 
      triangulation.execute_coarsening_and_refinement(); 
    } 

    pcout << "   Number of global active cells: " 
#ifdef DEAL_II_WITH_P4EST 
          << triangulation.n_global_active_cells() 
#else 
          << triangulation.n_active_cells() 
#endif 
          << std::endl; 

    dof_handler.distribute_dofs(fe); 

    pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
          << std::endl; 

// 我们生成悬挂节点约束，以确保解决方案的连续性。如同在 step-40 中，我们需要为约束矩阵配备本地相关自由度的IndexSet，以避免它在大问题中消耗过多的内存。接下来，问题的<code>MatrixFree</code>对象被设置。请注意，我们为共享内存并行化指定了一个特定的方案（因此，人们会使用多线程来实现节点内的并行化，而不是MPI；我们在这里选择了标准选项&mdash；如果我们想在程序中有一个以上的TBB线程的情况下禁用共享内存并行化，我们会选择 MatrixFree::AdditionalData::TasksParallelScheme::none).  另外请注意，我们没有使用默认的QGauss正交参数，而是提供一个QGaussLobatto正交公式来实现期望的行为。最后，三个求解向量被初始化。MatrixFree期望有一个特定的鬼魂索引布局（因为它在MPI本地数字中处理索引访问，需要在向量和MatrixFree之间匹配），所以我们只是要求它初始化向量，以确保鬼魂交换得到正确处理。

    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs); 
    constraints.clear(); 
    constraints.reinit(locally_relevant_dofs); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
    constraints.close(); 

    typename MatrixFree<dim>::AdditionalData additional_data; 
    additional_data.tasks_parallel_scheme = 
      MatrixFree<dim>::AdditionalData::TasksParallelScheme::partition_partition; 

    matrix_free_data.reinit(mapping, 
                            dof_handler, 
                            constraints, 
                            QGaussLobatto<1>(fe_degree + 1), 
                            additional_data); 

    matrix_free_data.initialize_dof_vector(solution); 
    old_solution.reinit(solution); 
    old_old_solution.reinit(solution); 
  } 

//  @sect4{SineGordonProblem::output_results}  

// 这个函数打印出解的规范，并将解的向量写到一个文件中。法线是标准的（除了我们需要累积所有处理器上的法线，用于并行网格，我们通过  VectorTools::compute_global_error()  函数来做），第二项类似于我们在  step-40  或  step-37  . 请注意，我们可以使用与计算过程中使用的相同的向量进行输出。无矩阵框架中的向量总是提供所有本地拥有的单元的全部信息（这也是本地评估中需要的），包括这些单元上的鬼向量条目。这是 VectorTools::integrate_difference() 函数以及DataOut中唯一需要的数据。这时唯一要做的就是确保在我们从矢量中读取数据之前更新其鬼魂值，并在完成后重置鬼魂值。这是一个只存在于 LinearAlgebra::distributed::Vector 类中的特性。另一方面，带有PETSc和Trilinos的分布式向量需要被复制到包括ghost值的特殊向量（见 step-40 中的相关章节 ）。如果我们还想访问幽灵单元上的所有自由度（例如，当计算使用单元边界上的解的跳跃的误差估计时），我们将需要更多的信息，并创建一个初始化了本地相关自由度的向量，就像在  step-40  中一样。还请注意，我们需要为输出分配约束条件

// --它们在计算过程中不被填充（相反，它们在无矩阵的方法中被实时插值  FEEvaluation::read_dof_values()).  
  template <int dim> 
  void 
  SineGordonProblem<dim>::output_results(const unsigned int timestep_number) 
  { 
    constraints.distribute(solution); 

    Vector<float> norm_per_cell(triangulation.n_active_cells()); 
    solution.update_ghost_values(); 
    VectorTools::integrate_difference(mapping, 
                                      dof_handler, 
                                      solution, 
                                      Functions::ZeroFunction<dim>(), 
                                      norm_per_cell, 
                                      QGauss<dim>(fe_degree + 1), 
                                      VectorTools::L2_norm); 
    const double solution_norm = 
      VectorTools::compute_global_error(triangulation, 
                                        norm_per_cell, 
                                        VectorTools::L2_norm); 

    pcout << "   Time:" << std::setw(8) << std::setprecision(3) << time 
          << ", solution norm: " << std::setprecision(5) << std::setw(7) 
          << solution_norm << std::endl; 

    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "solution"); 
    data_out.build_patches(mapping); 

    data_out.write_vtu_with_pvtu_record( 
      "./", "solution", timestep_number, MPI_COMM_WORLD, 3); 

    solution.zero_out_ghost_values(); 
  } 
// @sect4{SineGordonProblem::run}  

// 这个函数被主函数调用，并步入类的子程序中。

// 在打印了一些关于并行设置的信息后，第一个动作是设置网格和单元运算器。然后，根据构造函数中给出的CFL编号和最细的网格尺寸计算出时间步长。最细的网格尺寸计算为三角形中最后一个单元的直径，也就是网格中最细层次上的最后一个单元。这只适用于一个层次上的所有元素都具有相同尺寸的网格，否则就需要对所有单元进行循环。请注意，我们需要查询所有处理器的最细单元，因为不是所有的处理器都可能持有网格处于最细级别的区域。然后，我们重新调整一下时间步长，以准确地达到最后的时间。

  template <int dim> 
  void SineGordonProblem<dim>::run() 
  { 
    { 
      pcout << "Number of MPI ranks:            " 
            << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) << std::endl; 
      pcout << "Number of threads on each rank: " 
            << MultithreadInfo::n_threads() << std::endl; 
      const unsigned int n_vect_doubles = VectorizedArray<double>::size(); 
      const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles; 
      pcout << "Vectorization over " << n_vect_doubles 
            << " doubles = " << n_vect_bits << " bits (" 
            << Utilities::System::get_current_vectorization_level() << ")" 
            << std::endl 
            << std::endl; 
    } 
    make_grid_and_dofs(); 

    const double local_min_cell_diameter = 
      triangulation.last()->diameter() / std::sqrt(dim); 
    const double global_min_cell_diameter = 
      -Utilities::MPI::max(-local_min_cell_diameter, MPI_COMM_WORLD); 
    time_step = cfl_number * global_min_cell_diameter; 
    time_step = (final_time - time) / (int((final_time - time) / time_step)); 
    pcout << "   Time step size: " << time_step 
          << ", finest cell: " << global_min_cell_diameter << std::endl 
          << std::endl; 

// 接下来是初始值的设置。由于我们有一个两步的时间步进方法，我们还需要一个在时间步进时的解的值。为了得到准确的结果，需要根据初始时间的解的时间导数来计算，但是在这里我们忽略了这个困难，只是将其设置为该人工时间的初始值函数。

// 然后，我们继续将初始状态写入文件，并将两个初始解收集到 <tt>std::vector</tt> 的指针中，这些指针随后被 SineGordonOperation::apply() 函数消耗。接下来，根据文件顶部指定的有限元程度，建立一个 <code> SineGordonOperation class </code> 的实例。

    VectorTools::interpolate(mapping, 
                             dof_handler, 
                             InitialCondition<dim>(1, time), 
                             solution); 
    VectorTools::interpolate(mapping, 
                             dof_handler, 
                             InitialCondition<dim>(1, time - time_step), 
                             old_solution); 
    output_results(0); 

    std::vector<LinearAlgebra::distributed::Vector<double> *> 
      previous_solutions({&old_solution, &old_old_solution}); 

    SineGordonOperation<dim, fe_degree> sine_gordon_op(matrix_free_data, 
                                                       time_step); 

// 现在在时间步骤上循环。在每个迭代中，我们将解的向量移动一个，并调用`正弦戈登运算器'类的`应用'函数。然后，我们将解决方案写到一个文件中。我们对所需的计算时间和创建输出所需的时间进行计时，并在时间步长结束后报告这些数字。

// 注意这个交换是如何实现的。我们只是在两个向量上调用了交换方法，只交换了一些指针，而不需要复制数据，这在显式时间步进方法中是比较昂贵的操作。让我们来看看发生了什么。首先，我们交换 <code>old_solution</code> with <code>old_old_solution</code> ，这意味着 <code>old_old_solution</code> 得到 <code>old_solution</code> ，这就是我们所期望的。同样，在下一步中， <code>old_solution</code> gets the content from <code>solution</code> 也是如此。在这之后， <code>solution</code> 持有 <code>old_old_solution</code> ，但这将在这一步被覆盖。

    unsigned int timestep_number = 1; 

    Timer  timer; 
    double wtime       = 0; 
    double output_time = 0; 
    for (time += time_step; time <= final_time; 
         time += time_step, ++timestep_number) 
      { 
        timer.restart(); 
        old_old_solution.swap(old_solution); 
        old_solution.swap(solution); 
        sine_gordon_op.apply(solution, previous_solutions); 
        wtime += timer.wall_time(); 

        timer.restart(); 
        if (timestep_number % output_timestep_skip == 0) 
          output_results(timestep_number / output_timestep_skip); 

        output_time += timer.wall_time(); 
      } 
    timer.restart(); 
    output_results(timestep_number / output_timestep_skip + 1); 
    output_time += timer.wall_time(); 

    pcout << std::endl 
          << "   Performed " << timestep_number << " time steps." << std::endl; 

    pcout << "   Average wallclock time per time step: " 
          << wtime / timestep_number << "s" << std::endl; 

    pcout << "   Spent " << output_time << "s on output and " << wtime 
          << "s on computations." << std::endl; 
  } 
} // namespace Step48 

//  @sect3{The <code>main</code> function}  

// 与 step-40 中一样，我们在程序开始时初始化MPI。由于我们一般会将MPI并行化与线程混合在一起，所以我们也将MPI_InitFinalize中控制线程数量的第三个参数设置为无效数字，这意味着TBB库会自动选择线程的数量，通常为系统中可用的内核数量。作为一种选择，如果你想设置一个特定的线程数（例如，当需要只使用MPI时），你也可以手动设置这个数字。

int main(int argc, char **argv) 
{ 
  using namespace Step48; 
  using namespace dealii; 

  Utilities::MPI::MPI_InitFinalize mpi_initialization( 
    argc, argv, numbers::invalid_unsigned_int); 

  try 
    { 
      SineGordonProblem<dimension> sg_problem; 
      sg_problem.run(); 
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

