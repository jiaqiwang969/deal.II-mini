

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2000 - 2021 by the deal.II authors 
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
 * Author: Wolfgang Bangerth, University of Heidelberg, 2000 
 */ 


// @sect3{Include files}  

// 像往常一样，前几个include文件已经知道了，所以我们将不再评论它们。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/tensor.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

// 在这个例子中，我们需要矢量值的有限元。对这些的支持可以在下面的include文件中找到。

#include <deal.II/fe/fe_system.h> 

// 我们将用常规的Q1元素组成矢量值的有限元素，这些元素可以在这里找到，像往常一样。

#include <deal.II/fe/fe_q.h> 

// 这又是C++语言。

#include <fstream> 
#include <iostream> 

// 最后一步和以前的程序一样。特别是，就像在 step-7 中一样，我们把这个程序所特有的一切都打包到一个自己的命名空间中。

namespace Step8 
{ 
  using namespace dealii; 
// @sect3{The <code>ElasticProblem</code> class template}  

// 主类除了名称外，与 step-6 的例子相比几乎没有变化。

// 唯一的变化是为 <code>fe</code> 变量使用了一个不同的类。我们现在使用的不是FE_Q这样具体的有限元类，而是一个更通用的类，FESystem。事实上，FESystem本身并不是一个真正的有限元，因为它没有实现自己的形状函数。相反，它是一个可以用来将其他几个元素堆叠在一起形成一个矢量值的有限元的类。在我们的例子中，我们将组成 <code>FE_Q(1)</code> 对象的矢量值元素，如下所示，在这个类的构造函数中。

  template <int dim> 
  class ElasticProblem 
  { 
  public: 
    ElasticProblem(); 
    void run(); 

  private: 
    void setup_system(); 
    void assemble_system(); 
    void solve(); 
    void refine_grid(); 
    void output_results(const unsigned int cycle) const; 

    Triangulation<dim> triangulation; 
    DoFHandler<dim>    dof_handler; 

    FESystem<dim> fe; 

    AffineConstraints<double> constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> solution; 
    Vector<double> system_rhs; 
  }; 
// @sect3{Right hand side values}  

// 在进入主类的实现之前，我们声明并定义描述右手边的函数。这一次，右手边是向量值，解决方案也是如此，所以我们将更详细地描述为此所需的变化。

// 为了防止出现返回向量没有被设置成正确大小的情况，我们对这种情况进行了测试，否则将在函数的开始部分抛出一个异常。请注意，强制输出参数已经具有正确的大小是deal.II中的一个惯例，并且几乎在所有地方都强制执行。原因是，否则我们将不得不在函数开始时检查，并可能改变输出向量的大小。这很昂贵，而且几乎总是不必要的（对函数的第一次调用会将向量设置为正确的大小，随后的调用只需要做多余的检查）。此外，如果我们不能依赖向量已经具有正确大小的假设，那么检查和可能调整向量大小的操作是不能被删除的；这与Assert调用是一个契约，如果程序在优化模式下编译，Assert调用将被完全删除。

// 同样，如果由于某种意外，有人试图在只有一个空间维度的情况下编译和运行程序（在这种情况下，弹性方程没有什么意义，因为它们还原为普通的拉普拉斯方程），我们在第二个断言中终止程序。然而，该程序在三维空间中也能正常工作。

  template <int dim> 
  void right_hand_side(const std::vector<Point<dim>> &points, 
                       std::vector<Tensor<1, dim>> &  values) 
  { 
    Assert(values.size() == points.size(), 
           ExcDimensionMismatch(values.size(), points.size())); 
    Assert(dim >= 2, ExcNotImplemented()); 

// 该函数的其余部分实现了计算力值。我们将使用一个位于(0.5,0)和(-0.5,0)点周围的两个小圆圈（或球体，在3D中）的X方向的恒定（单位）力，以及位于原点周围的Y方向的力；在3D中，这些中心的Z分量也是零。

// 为此，让我们首先定义两个对象，表示这些区域的中心。请注意，在构建点对象时，所有的分量都被设置为零。

    Point<dim> point_1, point_2; 
    point_1(0) = 0.5; 
    point_2(0) = -0.5; 

    for (unsigned int point_n = 0; point_n < points.size(); ++point_n) 
      { 

// 如果 <code>points[point_n]</code> 处于围绕这些点之一的半径为0.2的圆（球）中，那么将X方向的力设置为1，否则为0。

        if (((points[point_n] - point_1).norm_square() < 0.2 * 0.2) || 
            ((points[point_n] - point_2).norm_square() < 0.2 * 0.2)) 
          values[point_n][0] = 1.0; 
        else 
          values[point_n][0] = 0.0; 

// 同样地，如果 <code>points[point_n]</code> 在原点附近，那么将y力设置为1，否则为0。

        if (points[point_n].norm_square() < 0.2 * 0.2) 
          values[point_n][1] = 1.0; 
        else 
          values[point_n][1] = 0.0; 
      } 
  } 

//  @sect3{The <code>ElasticProblem</code> class implementation}  
// @sect4{ElasticProblem::ElasticProblem constructor}  

// 下面是主类的构造函数。如前所述，我们想构造一个由多个标量有限元组成的矢量值有限元（即，我们想构造矢量值元素，使其每个矢量成分都由一个标量元素的形状函数组成）。当然，我们想堆叠在一起的标量有限元的数量等于解函数的分量数量，由于我们考虑每个空间方向上的位移，所以是 <code>dim</code> 。FESystem类可以处理这个问题：我们传递给它我们想组成系统的有限元，以及它的重复频率。

  template <int dim> 
  ElasticProblem<dim>::ElasticProblem() 
    : dof_handler(triangulation) 
    , fe(FE_Q<dim>(1), dim) 
  {} 

// 事实上，FESystem类还有几个构造函数，可以进行更复杂的操作，而不仅仅是将几个相同类型的标量有限元堆叠在一起；我们将在后面的例子中了解这些可能性。

//  @sect4{ElasticProblem::setup_system}  

// 设置方程组与 step-6 例子中使用的函数相同。DoFHandler类和这里使用的所有其他类都完全知道我们要使用的有限元是矢量值的，并且照顾到了有限元本身的矢量值。(事实上，它们不知道，但这不需要困扰你：因为它们只需要知道每个顶点、直线和单元有多少个自由度，它们不问它们代表什么，也就是说，考虑的有限元是矢量值的，还是例如在每个顶点上有几个自由度的标量Hermite元)。

  template <int dim> 
  void ElasticProblem<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 
    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 

    constraints.clear(); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             0, 
                                             Functions::ZeroFunction<dim>(dim), 
                                             constraints); 
    constraints.close(); 

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, 
                                    dsp, 
                                    constraints, 
                                    /*keep_constrained_dofs =  */ false);

    sparsity_pattern.copy_from(dsp); 

    system_matrix.reinit(sparsity_pattern); 
  } 
// @sect4{ElasticProblem::assemble_system}  

// 这个程序中最大的变化是创建矩阵和右手边，因为它们是取决于问题的。我们将一步一步地完成这个过程  step-  ，因为它比以前的例子要复杂一些。

// 然而，这个函数的前几部分和以前一样：设置一个合适的正交公式，为我们使用的（矢量值）有限元以及正交对象初始化一个FEValues对象，并声明了一些辅助数组。此外，我们还声明了永远相同的两个缩写。  <code>n_q_points</code>  和  <code>dofs_per_cell</code>  。每个单元的自由度数量，我们现在显然是从组成的有限元中询问，而不是从底层的标量Q1元中询问。在这里，它是 <code>dim</code> 乘以Q1元素的每个单元的自由度数，尽管这不是我们需要关心的明确知识。

  template <int dim> 
  void ElasticProblem<dim>::assemble_system() 
  { 
    QGauss<dim> quadrature_formula(fe.degree + 1); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

// 正如前面的例子所示，我们需要一个地方来存储单元格上所有正交点的系数值。在目前的情况下，我们有两个系数，lambda和mu。

    std::vector<double> lambda_values(n_q_points); 
    std::vector<double> mu_values(n_q_points); 

// 好吧，我们也可以省略上面的两个数组，因为我们将对lambda和mu使用常数系数，可以这样声明。它们都代表函数总是返回常量值1.0。尽管我们可以在矩阵的组合中省略各自的系数，但为了演示，我们在这里使用它们。

    Functions::ConstantFunction<dim> lambda(1.), mu(1.); 

// 和上面的两个常量函数一样，我们将在每个单元格中只调用一次函数right_hand_side，以使事情更简单。

    std::vector<Tensor<1, dim>> rhs_values(n_q_points); 

// 现在我们可以开始对所有单元格进行循环。

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_matrix = 0; 
        cell_rhs    = 0; 

        fe_values.reinit(cell); 

// 接下来我们得到正交点的系数值。同样，对于右手边也是如此。

        lambda.value_list(fe_values.get_quadrature_points(), lambda_values); 
        mu.value_list(fe_values.get_quadrature_points(), mu_values); 
        right_hand_side(fe_values.get_quadrature_points(), rhs_values); 

// 然后将局部刚度矩阵的条目和右手边的向量组合起来。这几乎是一对一地遵循本例介绍中描述的模式。 在位的几个评论之一是，我们可以计算数字  <code>comp(i)</code>  ，即使用下面的  <code>fe.system_to_component_index(i).first</code>  函数调用形状函数  <code>i</code>  的唯一非零向量成分的索引。

//（通过访问 <code>system_to_component_index</code> 函数返回值的 <code>first</code> 变量，你可能已经猜到其中还有更多的内容。事实上，该函数返回一个 <code>std::pair@<unsigned int，无符号int @></code>, ，其中第一个元素是 <code>comp(i)</code> ，第二个元素是介绍中也指出的值 <code>base(i)</code> ，即这个形状函数在这个组件中所有非零的形状函数中的索引，即介绍中的字典 <code>base(i)</code> 。不过，这不是我们通常感兴趣的数字）。)

// 有了这些知识，我们就可以把局部矩阵的贡献集合起来。

        for (const unsigned int i : fe_values.dof_indices()) 
          { 
            const unsigned int component_i = 
              fe.system_to_component_index(i).first; 

            for (const unsigned int j : fe_values.dof_indices()) 
              { 
                const unsigned int component_j = 
                  fe.system_to_component_index(j).first; 

                for (const unsigned int q_point : 
                     fe_values.quadrature_point_indices()) 
                  { 
                    cell_matrix(i, j) += 

// 第一个项是  $\lambda \partial_i u_i, \partial_j v_j) + (\mu \partial_i u_j, \partial_j v_i)$  。注意， <code>shape_grad(i,q_point)</code> 返回正交点q_point处第i个形状函数的唯一非零分量的梯度。梯度的分量 <code>comp(i)</code> 是第i个形状函数的唯一非零矢量分量相对于comp(i)th坐标的导数，由附加的括号访问。

                      (                                                  // 
                        (fe_values.shape_grad(i, q_point)[component_i] * // 
                         fe_values.shape_grad(j, q_point)[component_j] * // 
                         lambda_values[q_point])                         // 
                        +                                                // 
                        (fe_values.shape_grad(i, q_point)[component_j] * // 
                         fe_values.shape_grad(j, q_point)[component_i] * // 
                         mu_values[q_point])                             // 
                        +                                                // 

// 第二个项是  $(\mu \nabla u_i, \nabla v_j)$  。我们不需要访问梯度的具体分量，因为我们只需要计算两个梯度的标量乘积，这个问题由<tt>operator*</tt>的重载版本来负责，就像前面的例子一样。                            注意，通过使用<tt>?:</tt>操作符，我们只在<tt>component_i</tt>等于<tt>component_j</tt>时才这样做，否则会加上一个零（编译器会将其优化掉）。

                        ((component_i == component_j) ?        // 
                           (fe_values.shape_grad(i, q_point) * // 
                            fe_values.shape_grad(j, q_point) * // 
                            mu_values[q_point]) :              // 
                           0)                                  // 
                        ) *                                    // 
                      fe_values.JxW(q_point);                  // 
                  } 
              } 
          } 

// 组装右手边也和介绍中讨论的一样。

        for (const unsigned int i : fe_values.dof_indices()) 
          { 
            const unsigned int component_i = 
              fe.system_to_component_index(i).first; 

            for (const unsigned int q_point : 
                 fe_values.quadrature_point_indices()) 
              cell_rhs(i) += fe_values.shape_value(i, q_point) * 
                             rhs_values[q_point][component_i] * 
                             fe_values.JxW(q_point); 
          } 

// 从局部自由度到全局矩阵和右手向量的转移不取决于所考虑的方程，因此与之前所有的例子相同。

        cell->get_dof_indices(local_dof_indices); 
        constraints.distribute_local_to_global( 
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs); 
      } 
  } 

//  @sect4{ElasticProblem::solve}  

// 解算器并不关心方程组的来源，只要它保持正定和对称（这是使用CG解算器的要求），而这个方程组确实是这样。因此，我们不需要改变任何东西。

  template <int dim> 
  void ElasticProblem<dim>::solve() 
  { 
    SolverControl            solver_control(1000, 1e-12); 
    SolverCG<Vector<double>> cg(solver_control); 

    PreconditionSSOR<SparseMatrix<double>> preconditioner; 
    preconditioner.initialize(system_matrix, 1.2); 

    cg.solve(system_matrix, solution, system_rhs, preconditioner); 

    constraints.distribute(solution); 
  } 
// @sect4{ElasticProblem::refine_grid}  

// 对网格进行细化的函数与 step-6 的例子相同。正交公式再次适应了线性元素。请注意，误差估计器默认情况下是将从有限元解的所有分量中得到的估计值相加，也就是说，它使用所有方向的位移，权重相同。如果我们希望网格只适应x方向的位移，我们可以给函数传递一个额外的参数，告诉它这样做，而不考虑其他所有方向的位移作为误差指标。然而，对于目前的问题，似乎应该考虑所有的位移分量，而且权重相同。

  template <int dim> 
  void ElasticProblem<dim>::refine_grid() 
  { 
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

    KellyErrorEstimator<dim>::estimate(dof_handler, 
                                       QGauss<dim - 1>(fe.degree + 1), 
                                       {}, 
                                       solution, 
                                       estimated_error_per_cell); 

    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    estimated_error_per_cell, 
                                                    0.3, 
                                                    0.03); 

    triangulation.execute_coarsening_and_refinement(); 
  } 
// @sect4{ElasticProblem::output_results}  

// 输出的情况与之前的例子中已经显示过的差不多了。唯一的区别是，求解函数是矢量值的。DataOut类会自动处理这个问题，但我们必须给求解向量的每个分量一个不同的名字。

// 为了做到这一点， DataOut::add_vector() 函数想要一个字符串的向量。由于分量的数量与我们工作的维数相同，我们使用下面的 <code>switch</code> 语句。

// 我们注意到，一些图形程序对变量名称中允许的字符有限制。因此，deal.II只支持所有程序都支持的这些字符的最小子集。基本上，这些字符是字母、数字、下划线和其他一些字符，但特别是没有空格和减号/横线。否则该库将抛出一个异常，至少在调试模式下是这样。

// 在列出了1d、2d和3d的情况后，如果我们遇到一个我们没有考虑到的情况，让程序死亡是一种很好的风格。请记住，如果第一个参数中的条件没有得到满足，Assert宏会产生一个异常。当然，条件 <code>false</code> 永远不可能被满足，所以只要程序运行到默认语句，就会中止。

  template <int dim> 
  void ElasticProblem<dim>::output_results(const unsigned int cycle) const 
  { 
    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 

    std::vector<std::string> solution_names; 
    switch (dim) 
      { 
        case 1: 
          solution_names.emplace_back("displacement"); 
          break; 
        case 2: 
          solution_names.emplace_back("x_displacement"); 
          solution_names.emplace_back("y_displacement"); 
          break; 
        case 3: 
          solution_names.emplace_back("x_displacement"); 
          solution_names.emplace_back("y_displacement"); 
          solution_names.emplace_back("z_displacement"); 
          break; 
        default: 
          Assert(false, ExcNotImplemented()); 
      } 

// 在为解向量的不同组成部分设置了名称之后，我们可以将解向量添加到计划输出的数据向量列表中。请注意，下面的函数需要一个字符串向量作为第二个参数，而我们在以前所有例子中使用的函数在那里接受一个字符串。(事实上，我们之前使用的函数会将单个字符串转换成只有一个元素的向量，并将其转发给另一个函数)。

    data_out.add_data_vector(solution, solution_names); 
    data_out.build_patches(); 

    std::ofstream output("solution-" + std::to_string(cycle) + ".vtk"); 
    data_out.write_vtk(output); 
  } 

//  @sect4{ElasticProblem::run}  

//  <code>run</code> 函数所做的事情与 step-6 中的相同，比如说。这一次，我们使用平方[-1,1]^d作为域，在开始第一次迭代之前，我们在全局上对其进行了四次细化。

// 细化的原因有点意外：我们使用QGauss正交公式，在每个方向上有两个点用于整合右手边；这意味着每个单元上有四个正交点（在二维）。如果我们只对初始网格进行一次全局细化，那么在域上每个方向上就只有四个正交点。然而，右侧函数被选择为相当局部的，在这种情况下，纯属偶然，恰好所有的正交点都位于右侧函数为零的点上（用数学术语来说，正交点恰好在右侧函数的<i>support</i>之外的点上）。这样一来，用正交计算的右手向量将只包含零（尽管如果我们完全用积分计算右手向量的话，它当然会是非零的），方程组的解就是零向量，也就是一个处处为零的有限元函数。从某种意义上说，我们不应该对这种情况的发生感到惊讶，因为我们选择了一个完全不适合手头问题的初始网格。

// 不幸的是，如果离散解是常数，那么KellyErrorEstimator类计算的误差指标对每个单元来说也是零，对 Triangulation::refine_and_coarsen_fixed_number() 的调用将不会标记任何单元进行细化（如果每个单元的指示误差是零，为什么要这样做？因此，下一次迭代中的网格也将只由四个单元组成，同样的问题再次发生。

// 结论是：虽然我们当然不会把初始网格选择得非常适合问题的精确解决，但我们至少必须选择它，使它有机会捕捉到解决方案的重要特征。在这种情况下，它需要能够看到右手边的情况。因此，我们进行了四次全局细化。(任何更大的全局细化步骤当然也可以。)

  template <int dim> 
  void ElasticProblem<dim>::run() 
  { 
    for (unsigned int cycle = 0; cycle < 8; ++cycle) 
      { 
        std::cout << "Cycle " << cycle << ':' << std::endl; 

        if (cycle == 0) 
          { 
            GridGenerator::hyper_cube(triangulation, -1, 1); 
            triangulation.refine_global(4); 
          } 
        else 
          refine_grid(); 

        std::cout << "   Number of active cells:       " 
                  << triangulation.n_active_cells() << std::endl; 

        setup_system(); 

        std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
                  << std::endl; 

        assemble_system(); 
        solve(); 
        output_results(cycle); 
      } 
  } 
} // namespace Step8 
// @sect3{The <code>main</code> function}  

// 在上面最后一行关闭了 <code>Step8</code> 命名空间后，下面是程序的主要功能，又和 step-6 中一模一样（当然，除了改变了类名）。

int main() 
{ 
  try 
    { 
      Step8::ElasticProblem<2> elastic_problem_2d; 
      elastic_problem_2d.run(); 
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

