CCTest_file/step-27.cc

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
 * Authors: Wolfgang Bangerth, Texas A&M University, 2006, 2007; 
 *          Denis Davydov, University of Erlangen-Nuremberg, 2016; 
 *          Marc Fehling, Colorado State University, 2020. 
 */ 


// @sect3{Include files}  

// 前面几个文件已经在前面的例子中讲过了，因此不再做进一步的评论。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

// 这些是我们需要的新文件。第一个和第二个提供了FECollection和<i>hp</i>版本的FEValues类，如本程序介绍中所述。下一个文件提供了自动 $hp$ 适应的功能，为此我们将使用基于衰减系列扩展系数的估计算法，这是最后两个文件的一部分。

#include <deal.II/hp/fe_collection.h> 
#include <deal.II/hp/fe_values.h> 
#include <deal.II/hp/refinement.h> 
#include <deal.II/fe/fe_series.h> 
#include <deal.II/numerics/smoothness_estimator.h> 

// 最后一组包含文件是标准的C++头文件。

#include <fstream> 
#include <iostream> 

// 最后，这和以前的程序一样。

namespace Step27 
{ 
  using namespace dealii; 
// @sect3{The main class}  

// 这个程序的主类看起来非常像前几个教程程序中已经使用过的，例如  step-6  中的那个。主要的区别是我们将refine_grid和output_results函数合并为一个，因为我们还想输出一些用于决定如何细化网格的量（特别是估计的解决方案的平滑度）。

// 就成员变量而言，我们使用与 step-6 中相同的结构，但我们需要集合来代替单个的有限元、正交和面状正交对象。我们将在类的构造函数中填充这些集合。最后一个变量， <code>max_degree</code> ，表示所用形状函数的最大多项式程度。

  template <int dim> 
  class LaplaceProblem 
  { 
  public: 
    LaplaceProblem(); 
    ~LaplaceProblem(); 

    void run(); 

  private: 
    void setup_system(); 
    void assemble_system(); 
    void solve(); 
    void create_coarse_grid(); 
    void postprocess(const unsigned int cycle); 

    Triangulation<dim> triangulation; 

 
 
    hp::QCollection<dim>     quadrature_collection; 
    hp::QCollection<dim - 1> face_quadrature_collection; 

    AffineConstraints<double> constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> solution; 
    Vector<double> system_rhs; 

    const unsigned int max_degree; 
  }; 

//  @sect3{Equation data}  

// 接下来，让我们为这个问题定义右手边的函数。它在1d中是 $x+1$ ，在2d中是 $(x+1)(y+1)$ ，以此类推。

  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component) const override; 
  }; 

  template <int dim> 
  double RightHandSide<dim>::value(const Point<dim> &p, 
                                   const unsigned int /*component*/) const 
  { 
    double product = 1; 
    for (unsigned int d = 0; d < dim; ++d) 
      product *= (p[d] + 1); 
    return product; 
  } 

//  @sect3{Implementation of the main class}  
// @sect4{LaplaceProblem::LaplaceProblem constructor}  

// 这个类的构造函数是相当直接的。它将DoFHandler对象与三角形相关联，然后将最大多项式度数设置为7（在1d和2d中）或5（在3d及以上）。我们这样做是因为使用高阶多项式度数会变得非常昂贵，尤其是在更高的空间维度上。

// 在这之后，我们填充有限元、单元和面的四分法对象集合。我们从二次元开始，每个正交公式的选择都是为了适合 hp::FECollection 对象中的匹配有限元。

  template <int dim> 
  LaplaceProblem<dim>::LaplaceProblem() 
    : dof_handler(triangulation) 
    , max_degree(dim <= 2 ? 7 : 5) 
  { 
    for (unsigned int degree = 2; degree <= max_degree; ++degree) 
      { 
        fe_collection.push_back(FE_Q<dim>(degree)); 
        quadrature_collection.push_back(QGauss<dim>(degree + 1)); 
        face_quadrature_collection.push_back(QGauss<dim - 1>(degree + 1)); 
      } 
  } 
// @sect4{LaplaceProblem::~LaplaceProblem destructor}  

// 解构器与我们在  step-6  中已经做过的没有变化。

  template <int dim> 
  LaplaceProblem<dim>::~LaplaceProblem() 
  { 
    dof_handler.clear(); 
  } 
// @sect4{LaplaceProblem::setup_system}  

// 这个函数又是对我们在  step-6  中已经做过的事情的逐字复制。尽管函数调用的名称和参数完全相同，但内部使用的算法在某些方面是不同的，因为这里的dof_handler变量是在  $hp$  -mode。

  template <int dim> 
  void LaplaceProblem<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe_collection); 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 

    constraints.clear(); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             0, 
                                             Functions::ZeroFunction<dim>(), 
                                             constraints); 
    constraints.close(); 

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false); 
    sparsity_pattern.copy_from(dsp); 

    system_matrix.reinit(sparsity_pattern); 
  } 

//  @sect4{LaplaceProblem::assemble_system}  

// 这是一个从每个单元的局部贡献中集合全局矩阵和右侧向量的函数。它的主要工作与之前许多教程中描述的一样。重要的差异是<i>hp</i>有限元方法所需要的。特别是，我们需要使用FEValues对象的集合（通过 hp::FEValues 类实现），并且在将局部贡献复制到全局对象时，我们必须消除受限自由度。这两点在本程序的介绍中都有详细解释。

// 还有一个小问题是，由于我们在不同的单元格中使用了不同的多项式度数，持有局部贡献的矩阵和向量在所有单元格中的大小不尽相同。因此，在所有单元的循环开始时，我们每次都必须将它们的大小调整到正确的大小（由 <code>dofs_per_cell</code> 给出）。因为这些类的实现方式是减少矩阵或向量的大小不会释放当前分配的内存（除非新的大小为零），所以在循环开始时调整大小的过程只需要在最初几次迭代中重新分配内存。一旦我们在一个单元中找到了最大的有限元度，就不会再发生重新分配，因为所有后续的 <code>reinit</code> 调用只会将大小设置为适合当前分配的内存。这一点很重要，因为分配内存是很昂贵的，而且每次我们访问一个新的单元时都这样做会花费大量的计算时间。

  template <int dim> 
  void LaplaceProblem<dim>::assemble_system() 
  { 
    hp::FEValues<dim> hp_fe_values(fe_collection, 
                                   quadrature_collection, 
                                   update_values | update_gradients | 
                                     update_quadrature_points | 
                                     update_JxW_values); 

    RightHandSide<dim> rhs_function; 

    FullMatrix<double> cell_matrix; 
    Vector<double>     cell_rhs; 

    std::vector<types::global_dof_index> local_dof_indices; 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell(); 

        cell_matrix.reinit(dofs_per_cell, dofs_per_cell); 
        cell_matrix = 0; 

        cell_rhs.reinit(dofs_per_cell); 
        cell_rhs = 0; 

        hp_fe_values.reinit(cell); 

        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values(); 

        std::vector<double> rhs_values(fe_values.n_quadrature_points); 
        rhs_function.value_list(fe_values.get_quadrature_points(), rhs_values); 

        for (unsigned int q_point = 0; q_point < fe_values.n_quadrature_points; 
             ++q_point) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                cell_matrix(i, j) += 
                  (fe_values.shape_grad(i, q_point) * // grad phi_i(x_q) 
                   fe_values.shape_grad(j, q_point) * // grad phi_j(x_q) 
                   fe_values.JxW(q_point));           // dx 

              cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q) 
                              rhs_values[q_point] *               // f(x_q) 
                              fe_values.JxW(q_point));            // dx 
            } 

        local_dof_indices.resize(dofs_per_cell); 
        cell->get_dof_indices(local_dof_indices); 

        constraints.distribute_local_to_global( 
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs); 
      } 
  } 

//  @sect4{LaplaceProblem::solve}  

// 解决线性系统的函数与之前的例子完全没有变化。我们只是试图将初始残差（相当于右手边的 $l_2$ 准则）减少一定的系数。

  template <int dim> 
  void LaplaceProblem<dim>::solve() 
  { 
    SolverControl            solver_control(system_rhs.size(), 
                                 1e-12 * system_rhs.l2_norm()); 
    SolverCG<Vector<double>> cg(solver_control); 

    PreconditionSSOR<SparseMatrix<double>> preconditioner; 
    preconditioner.initialize(system_matrix, 1.2); 

    cg.solve(system_matrix, solution, system_rhs, preconditioner); 

    constraints.distribute(solution); 
  } 

//  @sect4{LaplaceProblem::postprocess}  

// 解完线性系统后，我们要对解进行后处理。在这里，我们所做的就是估计误差，估计解的局部平滑度，如介绍中所述，然后写出图形输出，最后根据之前计算的指标细化 $h$ 和 $p$ 中的网格。我们在同一个函数中完成这一切，因为我们希望估计的误差和平滑度指标不仅用于细化，而且还包括在图形输出中。

  template <int dim> 
  void LaplaceProblem<dim>::postprocess(const unsigned int cycle) 
  { 

// 让我们开始计算估计的误差和平滑度指标，这两个指标对于我们三角测量的每个活动单元来说都是一个数字。对于误差指标，我们一如既往地使用KellyErrorEstimator类。

    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      face_quadrature_collection, 
      std::map<types::boundary_id, const Function<dim> *>(), 
      solution, 
      estimated_error_per_cell); 

// 估计平滑度是用介绍中所述的衰减膨胀系数的方法进行的。我们首先需要创建一个对象，能够将每一个单元上的有限元解转化为一串傅里叶级数系数。SmoothnessEstimator命名空间为这样一个 FESeries::Fourier 对象提供了一个工厂函数，它为估计平滑度的过程进行了优化。然后在最后一个函数中实际确定每个单独单元上的傅里叶系数的衰减情况。

    Vector<float> smoothness_indicators(triangulation.n_active_cells()); 
    FESeries::Fourier<dim> fourier = 
      SmoothnessEstimator::Fourier::default_fe_series(fe_collection); 
    SmoothnessEstimator::Fourier::coefficient_decay(fourier, 
                                                    dof_handler, 
                                                    solution, 
                                                    smoothness_indicators); 

// 接下来我们要生成图形输出。除了上面得出的两个估计量之外，我们还想输出网格上每个元素所使用的有限元的多项式程度。

// 要做到这一点，我们需要在所有单元上循环，用  <code>cell-@>active_fe_index()</code>  轮询它们的活动有限元索引。然后我们使用这个操作的结果，在有限元集合中查询具有该索引的有限元，最后确定该元素的多项式程度。我们将结果放入一个矢量，每个单元有一个元素。DataOut类要求这是一个 <code>float</code> or <code>double</code> 的向量，尽管我们的值都是整数，所以我们就用这个向量。

    { 
      Vector<float> fe_degrees(triangulation.n_active_cells()); 
      for (const auto &cell : dof_handler.active_cell_iterators()) 
        fe_degrees(cell->active_cell_index()) = 
          fe_collection[cell->active_fe_index()].degree; 

// 现在有了所有的数据向量--解决方案、估计误差和平滑度指标以及有限元度--我们创建一个用于图形输出的DataOut对象并附加所有数据。

      DataOut<dim> data_out; 

      data_out.attach_dof_handler(dof_handler); 
      data_out.add_data_vector(solution, "solution"); 
      data_out.add_data_vector(estimated_error_per_cell, "error"); 
      data_out.add_data_vector(smoothness_indicators, "smoothness"); 
      data_out.add_data_vector(fe_degrees, "fe_degree"); 
      data_out.build_patches(); 

// 生成输出的最后一步是确定一个文件名，打开文件，并将数据写入其中（这里，我们使用VTK格式）。

      const std::string filename = 
        "solution-" + Utilities::int_to_string(cycle, 2) + ".vtk"; 
      std::ofstream output(filename); 
      data_out.write_vtk(output); 
    } 

// 在这之后，我们想在 $h$ 和 $p$ 两个地方实际细化网格。我们要做的是：首先，我们用估计的误差来标记那些误差最大的单元，以便进行细化。这就是我们一直以来的做法。

    { 
      GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                      estimated_error_per_cell, 
                                                      0.3, 
                                                      0.03); 

// 接下来我们要弄清楚哪些被标记为细化的单元格实际上应该增加 $p$ 而不是减少 $h$ 。我们在这里选择的策略是，我们查看那些被标记为细化的单元格的平滑度指标，并为那些平滑度大于某个相对阈值的单元格增加 $p$ 。换句话说，对于每一个(i)细化标志被设置，(ii)平滑度指标大于阈值，以及(iii)我们在有限元集合中仍有一个多项式度数高于当前度数的有限元的单元，我们将分配一个未来的FE指数，对应于一个比当前度数高一的多项式。下面的函数正是能够做到这一点。在没有更好的策略的情况下，我们将通过在标记为细化的单元上的最小和最大平滑度指标之间进行插值来设置阈值。由于角部奇点具有很强的局部性，我们将支持 $p$ 。

// - 而不是 $h$  - 精细化的数量。我们通过设置0.2的小插值系数，以低门槛实现这一点。用同样的方法，我们处理那些要被粗化的单元，当它们的平滑度指标低于在要粗化的单元上确定的相应阈值时，减少它们的多项式程度。

      hp::Refinement::p_adaptivity_from_relative_threshold( 
        dof_handler, smoothness_indicators, 0.2, 0.2); 

// 上面的函数只决定了多项式程度是否会通过未来的FE指数发生变化，但并没有操作 $h$  -细化标志。因此，对于被标记为两个细化类别的单元格，我们更倾向于 $p$  。

// 而不是 $h$  -细化。下面的函数调用确保只有 $p$ 中的一个

// - 或  $h$  - 精炼中的一种，而不是同时实施两种。

      hp::Refinement::choose_p_over_h(dof_handler); 

// 对于网格自适应细化，我们通过调用 Triangulation::prepare_coarsening_and_refinement(). 将相邻单元的细化水平差限制为1来确保2:1的网格平衡。 我们希望对相邻单元的p水平实现类似的效果：未来有限元的水平差不允许超过指定的差。通过其默认参数，调用 hp::Refinement::limit_p_level_difference() 可以确保它们的级差被限制在1以内。这不一定会减少域中的悬挂节点的数量，但可以确保高阶多项式不会被限制在面的低得多的多项式上，例如五阶多项式到二阶多项式。

      triangulation.prepare_coarsening_and_refinement(); 
      hp::Refinement::limit_p_level_difference(dof_handler); 

// 在这个过程结束后，我们再细化网格。在这个过程中，正在进行分割的单元的子单元会继承其母单元的有限元索引。此外，未来的有限元指数将变成活动的，因此新的有限元将在下一次调用 DoFHandler::distribute_dofs(). 后被分配给单元。
      triangulation.execute_coarsening_and_refinement(); 
    } 
  } 
// @sect4{LaplaceProblem::create_coarse_grid}  

// 在创建初始网格时，会用到下面这个函数。我们想要创建的网格实际上与 step-14 中的网格类似，即中间有方孔的方形域。它可以由完全相同的函数生成。然而，由于它的实现只是2d情况下的一种特殊化，我们将介绍一种不同的方法来创建这个域，它是独立于维度的。

// 我们首先创建一个有足够单元的超立方体三角形，这样它就已经包含了我们想要的域 $[-1,1]^d$ ，并细分为 $4^d$ 单元。然后，我们通过测试每个单元上顶点的坐标值来移除域中心的那些单元。最后，我们像往常一样对如此创建的网格进行全局细化。

  template <int dim> 
  void LaplaceProblem<dim>::create_coarse_grid() 
  { 
    Triangulation<dim> cube; 
    GridGenerator::subdivided_hyper_cube(cube, 4, -1., 1.); 

    std::set<typename Triangulation<dim>::active_cell_iterator> cells_to_remove; 
    for (const auto &cell : cube.active_cell_iterators()) 
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v) 
        if (cell->vertex(v).square() < .1) 
          cells_to_remove.insert(cell); 

    GridGenerator::create_triangulation_with_removed_cells(cube, 
                                                           cells_to_remove, 
                                                           triangulation); 

    triangulation.refine_global(3); 
  } 

//  @sect4{LaplaceProblem::run}  

// 这个函数实现了程序的逻辑，就像以前大多数程序中的相应函数一样，例如见  step-6  。

// 基本上，它包含了自适应循环：在第一次迭代中创建一个粗略的网格，然后建立线性系统，对其进行组合，求解，并对解进行后处理，包括网格细化。然后再重新开始。同时，也为那些盯着屏幕试图弄清楚程序是干什么的人输出一些信息。

  template <int dim> 
  void LaplaceProblem<dim>::run() 
  { 
    for (unsigned int cycle = 0; cycle < 6; ++cycle) 
      { 
        std::cout << "Cycle " << cycle << ':' << std::endl; 

        if (cycle == 0) 
          create_coarse_grid(); 

        setup_system(); 

        std::cout << "   Number of active cells      : " 
                  << triangulation.n_active_cells() << std::endl 
                  << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
                  << std::endl 
                  << "   Number of constraints       : " 
                  << constraints.n_constraints() << std::endl; 

        assemble_system(); 
        solve(); 
        postprocess(cycle); 
      } 
  } 
} // namespace Step27 
// @sect3{The main function}  

// 主函数仍然是我们之前的版本：将创建和运行一个主类的对象包装成一个 <code>try</code> 块，并捕捉任何抛出的异常，从而在出现问题时产生有意义的输出。

int main() 
{ 
  try 
    { 
      using namespace Step27; 

      LaplaceProblem<2> laplace_problem; 
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

