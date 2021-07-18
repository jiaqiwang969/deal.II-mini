

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2000 - 2020 by the deal.II authors 
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

// 前面几个文件已经在前面的例子中讲过了，因此不再做进一步的评论。

#include <deal.II/base/quadrature_lib.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_values.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 

#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/vector.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 

#include <fstream> 

// 从下面的include文件中我们将导入H1-conforming有限元形状函数的声明。这个有限元系列被称为  <code>FE_Q</code>  ，在之前的所有例子中已经被用来定义通常的双线性或三线性元素，但我们现在将用它来定义双二次元元素。

#include <deal.II/fe/fe_q.h> 

// 我们不会像前面的例子那样从文件中读取网格，而是使用库中的一个函数生成网格。然而，我们将希望在每一步中写出局部细化的网格（只是网格，而不是解决方案），所以我们需要以下的include文件，而不是 <code>grid_in.h</code>  。

#include <deal.II/grid/grid_out.h> 

// 当使用局部细化网格时，我们会得到所谓的<code>悬空节点</code>。然而，标准的有限元方法假定离散的解空间是连续的，所以我们需要确保悬挂节点上的自由度符合一些约束条件，这样全局解是连续的。我们也要在这个对象中存储边界条件。下面的文件包含一个用来处理这些约束条件的类。

#include <deal.II/lac/affine_constraints.h> 

// 为了在本地细化我们的网格，我们需要一个来自库的函数，根据我们计算的误差指标来决定哪些单元需要细化或粗化。这个函数被定义在这里。

#include <deal.II/grid/grid_refinement.h> 

// 最后，我们需要一个简单的方法来实际计算基于某种误差估计的细化指标。虽然一般来说，适应性是非常具体的问题，但以下文件中的误差指标通常会对一大类问题产生相当好的适应网格。

#include <deal.II/numerics/error_estimator.h> 

// 最后，这和以前的程序一样。

using namespace dealii; 
// @sect3{The <code>Step6</code> class template}  

// 主类又是几乎没有变化的。然而，我们增加了两项内容：我们增加了 <code>refine_grid</code> 函数，该函数用于自适应地细化网格（而不是之前例子中的全局细化），还有一个变量，它将保存约束条件。

template <int dim> 
class Step6 
{ 
public: 
  Step6(); 

  void run(); 

private: 
  void setup_system(); 
  void assemble_system(); 
  void solve(); 
  void refine_grid(); 
  void output_results(const unsigned int cycle) const; 

  Triangulation<dim> triangulation; 

  FE_Q<dim>       fe; 
  DoFHandler<dim> dof_handler; 

// 这是主类中的新变量。我们需要一个对象，它持有一个约束条件的列表，以保持悬挂节点和边界条件。

  AffineConstraints<double> constraints; 

  SparseMatrix<double> system_matrix; 
  SparsityPattern      sparsity_pattern; 

  Vector<double> solution; 
  Vector<double> system_rhs; 
}; 
// @sect3{Nonconstant coefficients}  

//非恒定系数的实现是逐字复制自  step-5  。

template <int dim> 
double coefficient(const Point<dim> &p) 
{ 
  if (p.square() < 0.5 * 0.5) 
    return 20; 
  else 
    return 1; 
} 

//  @sect3{The <code>Step6</code> class implementation}  
// @sect4{Step6::Step6}  

// 这个类的构造函数与之前的基本相同，但这一次我们要使用二次元。为此，我们只需用所需的多项式度数（这里是 <code>2</code> ）替换构造函数参数（在之前的所有例子中是 <code>1</code> ）。

template <int dim> 
Step6<dim>::Step6() 
  : fe(2) 
  , dof_handler(triangulation) 
{} 

//  @sect4{Step6::setup_system}  

// 下一个函数设置了所有描述线性有限元问题的变量，如DoFHandler、矩阵和向量。与我们在 step-5 中所做的不同的是，我们现在还必须处理悬挂节点约束。这些约束几乎完全由库来处理，也就是说，你只需要知道它们的存在以及如何获得它们，但你不需要知道它们是如何形成的，也不需要知道对它们到底做了什么。

// 在函数的开头，你会发现所有与 step-5 中相同的东西：设置自由度（这次我们有二次元，但从用户代码的角度看与线性--或任何其他程度的情况没有区别），生成稀疏模式，并初始化解和右手向量。请注意，现在每行的稀疏模式将有更多的条目，因为现在每个单元有9个自由度（而不是只有4个），它们可以相互耦合。

template <int dim> 
void Step6<dim>::setup_system() 
{ 
  dof_handler.distribute_dofs(fe); 

  solution.reinit(dof_handler.n_dofs()); 
  system_rhs.reinit(dof_handler.n_dofs()); 

// 我们现在可以用悬挂节点的约束来填充AffineConstraints对象。由于我们将在一个循环中调用这个函数，所以我们首先清除上一个系统中的当前约束集，然后计算新的约束。

  constraints.clear(); 
  DoFTools::make_hanging_node_constraints(dof_handler, constraints); 

// 现在我们准备用指标0（整个边界）来插值边界值，并将得到的约束存储在我们的 <code>constraints</code> 对象中。请注意，我们并不像在前面的步骤中那样，在装配后应用边界条件：相反，我们将所有的约束条件放在AffineConstraints对象中的我们的函数空间。我们可以按任何顺序向AffineConstraints对象添加约束：如果两个约束发生冲突，那么约束矩阵要么中止，要么通过Assert宏抛出一个异常。

  VectorTools::interpolate_boundary_values(dof_handler, 
                                           0, 
                                           Functions::ZeroFunction<dim>(), 
                                           constraints); 

// 在所有约束条件被添加之后，需要对它们进行排序和重新排列，以便更有效地执行一些操作。这种后处理是用 <code>close()</code> 函数完成的，之后就不能再添加任何约束了。

  constraints.close(); 

// 现在我们首先建立我们的压缩稀疏模式，就像我们在前面的例子中做的那样。然而，我们并没有立即将其复制到最终的稀疏度模式中。 请注意，我们调用了make_sparsity_pattern的一个变体，它把AffineConstraints对象作为第三个参数。我们通过将参数 <code>keep_constrained_dofs</code> 设置为false（换句话说，我们永远不会写入矩阵中对应于受限自由度的条目），让该例程知道我们永远不会写入 <code>constraints</code> 所给的位置。如果我们在装配后对约束进行压缩，我们就必须通过 <code>true</code> 来代替，因为这样我们就会先写进这些位置，然后在压缩过程中再将它们设置为零。

  DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
  DoFTools::make_sparsity_pattern(dof_handler, 
                                  dsp, 
                                  constraints, 
                                  /*keep_constrained_dofs =  */ false);

// 现在，矩阵的所有非零条目都是已知的（即那些来自定期组装矩阵的条目和那些通过消除约束引入的条目）。我们可以将我们的中间对象复制到稀疏模式中。

  sparsity_pattern.copy_from(dsp); 

// 我们现在可以，最后，初始化稀疏矩阵。

  system_matrix.reinit(sparsity_pattern); 
} 
// @sect4{Step6::assemble_system}  

// 接下来，我们要对矩阵进行组装。然而，为了将每个单元上的本地矩阵和向量复制到全局系统中，我们不再使用手写的循环。相反，我们使用 AffineConstraints::distribute_local_to_global() ，在内部执行这个循环，同时对对应于受限自由度的行和列进行高斯消除。

// 构成局部贡献的其余代码保持不变。然而，值得注意的是，在引擎盖下，有几件事与以前不同。首先，变量 <code>dofs_per_cell</code> 和返回值 <code>quadrature_formula.size()</code> 现在各为9，以前是4。引入这样的变量作为缩写是一个很好的策略，可以使代码在不同的元素下工作，而不需要改变太多的代码。其次， <code>fe_values</code> 对象当然也需要做其他事情，因为现在的形状函数是二次的，而不是线性的，在每个坐标变量中。不过，这也是完全由库来处理的事情。

template <int dim> 
void Step6<dim>::assemble_system() 
{ 
  const QGauss<dim> quadrature_formula(fe.degree + 1); 

  FEValues<dim> fe_values(fe, 
                          quadrature_formula, 
                          update_values | update_gradients | 
                            update_quadrature_points | update_JxW_values); 

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
  Vector<double>     cell_rhs(dofs_per_cell); 

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

  for (const auto &cell : dof_handler.active_cell_iterators()) 
    { 
      cell_matrix = 0; 
      cell_rhs    = 0; 

      fe_values.reinit(cell); 

      for (const unsigned int q_index : fe_values.quadrature_point_indices()) 
        { 
          const double current_coefficient = 
            coefficient(fe_values.quadrature_point(q_index)); 
          for (const unsigned int i : fe_values.dof_indices()) 
            { 
              for (const unsigned int j : fe_values.dof_indices()) 
                cell_matrix(i, j) += 
                  (current_coefficient *              // a(x_q) 
                   fe_values.shape_grad(i, q_index) * // grad phi_i(x_q) 
                   fe_values.shape_grad(j, q_index) * // grad phi_j(x_q) 
                   fe_values.JxW(q_index));           // dx 

              cell_rhs(i) += (1.0 *                               // f(x) 
                              fe_values.shape_value(i, q_index) * // phi_i(x_q) 
                              fe_values.JxW(q_index));            // dx 
            } 
        } 

// 最后，将 @p cell_matrix 和 @p cell_rhs 中的贡献转移到全局对象中。

      cell->get_dof_indices(local_dof_indices); 
      constraints.distribute_local_to_global( 
        cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs); 
    } 

// 现在我们已经完成了线性系统的组装。约束矩阵处理了应用边界条件的问题，也消除了悬挂的节点约束。受约束的节点仍然在线性系统中（在矩阵的对角线上有一个非零条目，选择的方式是使矩阵具有良好的条件，并且这一行的所有其他条目都被设置为零），但是计算出来的值是无效的（也就是说， <code>system_rhs</code> 中的相应条目目前是没有意义的）。我们在 <code>solve</code> 函数的最后为这些节点计算出正确的值。

} 
// @sect4{Step6::solve}  

// 我们继续逐步改进。解决线性系统的函数再次使用了SSOR预处理程序，除了我们必须加入悬空节点约束外，其他的都没有改变。如上所述，通过对矩阵的行和列进行特殊处理，从AffineConstraints对象中删除了对应于悬挂节点约束和边界值的自由度。这样一来，这些自由度的值在求解线性系统后就有了错误的、但定义明确的值。然后我们要做的就是利用约束条件给它们分配它们应该有的值。这个过程被称为 <code>distributing</code> 约束，从无约束的节点的值中计算出约束节点的值，只需要一个额外的函数调用，你可以在这个函数的末尾找到。

template <int dim> 
void Step6<dim>::solve() 
{ 
  SolverControl            solver_control(1000, 1e-12); 
  SolverCG<Vector<double>> solver(solver_control); 

  PreconditionSSOR<SparseMatrix<double>> preconditioner; 
  preconditioner.initialize(system_matrix, 1.2); 

  solver.solve(system_matrix, solution, system_rhs, preconditioner); 

  constraints.distribute(solution); 
} 
// @sect4{Step6::refine_grid}  

// 我们使用一个复杂的误差估计方案来细化网格，而不是全局细化。我们将使用KellyErrorEstimator类，该类实现了拉普拉斯方程的误差估计器；原则上它可以处理可变系数，但我们不会使用这些高级功能，而是使用其最简单的形式，因为我们对定量结果不感兴趣，只对生成局部细化网格的快速方法感兴趣。

// 尽管Kelly等人得出的误差估计器最初是为拉普拉斯方程开发的，但我们发现它也很适合于为一类广泛的问题快速生成局部细化网格。这个误差估计器使用了解梯度在单元面上的跳跃（这是一个测量二阶导数的方法），并将其按单元的大小进行缩放。因此，它是对每个单元的解的局部平滑性的测量，因此可以理解，它对双曲运输问题或波浪方程也能产生合理的网格，尽管这些网格与专门针对该问题的方法相比肯定是次优的。因此，这个误差估计器可以理解为测试自适应程序的一种快速方法。

// 估算器的工作方式是将描述自由度的 <code>DoFHandler</code> 对象和每个自由度的数值向量作为输入，为三角剖分的每个活动单元计算一个指标值（即每个活动单元一个数值）。为此，它需要两个额外的信息：一个面部正交公式，即 <code>dim-1</code> 维物体上的正交公式。我们再次使用3点高斯法则，这个选择与本程序中的双二次方有限元形状函数是一致和合适的。当然，什么是合适的正交规则取决于对误差估计器评估解场的方式的了解。如上所述，梯度的跳跃在每个面上都是集成的，对于本例中使用的二次元元素来说，这将是每个面上的二次元函数。然而，事实上，它是梯度跳动的平方，正如该类文件中所解释的那样，这是一个二次函数，对于它来说，3点高斯公式就足够了，因为它可以精确地整合5阶以下的多项式。)

// 其次，该函数需要一个边界指示器的列表，用于那些我们施加了 $\partial_n u(\mathbf x) = h(\mathbf x)$ 类诺伊曼值的边界，以及每个此类边界的函数 $h(\mathbf x)$ 。这些信息由一个从边界指标到描述诺伊曼边界值的函数对象的映射来表示。在本例程序中，我们不使用诺伊曼边界值，所以这个映射是空的，实际上是在函数调用期望得到相应函数参数的地方使用映射的默认构造器构造的。

// 输出是一个所有活动单元的值的向量。虽然非常精确地计算一个解的自由度的<b>value</b>可能是有意义的，但通常没有必要特别精确地计算一个单元上的解对应的<b>error indicator</b>。因此，我们通常使用一个浮点数的向量而不是一个双数的向量来表示误差指标。

template <int dim> 
void Step6<dim>::refine_grid() 
{ 
  Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

  KellyErrorEstimator<dim>::estimate(dof_handler, 
                                     QGauss<dim - 1>(fe.degree + 1), 
                                     {}, 
                                     solution, 
                                     estimated_error_per_cell); 

// 上述函数为 <code>estimated_error_per_cell</code> 数组中的每个单元格返回一个错误指标值。现在的细化工作如下：细化那些误差值最高的30%的单元，粗化那些误差值最低的3%的单元。

// 人们可以很容易地验证，如果第二个数字为零，这大约会导致在两个空间维度上的每一步的细胞翻倍，因为对于每一个30%的细胞，四个新的将被替换，而其余70%的细胞保持不动。在实践中，通常会产生一些更多的单元，因为不允许一个单元被精炼两次而相邻的单元没有被精炼；在这种情况下，相邻的单元也会被精炼。

// 在许多应用中，被粗化的单元格数量将被设置为大于3%的数值。一个非零的值是很有用的，特别是当初始（粗）网格由于某种原因已经相当精细时。在这种情况下，可能有必要在某些区域进行细化，而在另一些区域进行粗化是有用的。在我们这里，初始网格是非常粗的，所以粗化只需要在一些可能发生过度细化的区域。因此，一个小的、非零的值在这里是合适的。

// 下面的函数现在接受这些细化指标，并使用上述方法对三角形的一些单元进行细化或粗化标记。它来自一个实现了几种不同算法的类，可以根据单元的误差指标来细化三角形。

  GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                  estimated_error_per_cell, 
                                                  0.3, 
                                                  0.03); 

// 在前一个函数退出后，一些单元被标记为细化，另一些单元被标记为粗化。然而，细化或粗化本身并没有被执行，因为有些情况下，进一步修改这些标志是有用的。在这里，我们不想做任何这样的事情，所以我们可以告诉三角计算执行单元格被标记的动作。

  triangulation.execute_coarsening_and_refinement(); 
} 
// @sect4{Step6::output_results}  

// 在每个网格的计算结束后，在我们继续下一个网格细化周期之前，我们要输出这个周期的结果。

// 我们已经在 step-1 中看到了如何实现对网格本身的输出。在这里，我们改变一些东西。  <ol>  
// <li>  我们使用两种不同的格式。gnuplot和VTU。 </li>  
// <li>  我们在输出文件名中嵌入了周期号。 </li>  
// <li>  对于gnuplot输出，我们设置了一个 GridOutFlags::Gnuplot 对象，以提供一些额外的可视化参数，使边缘看起来是弯曲的。这在  step-10  中有进一步的详细解释。 </li>  
// </ol>  
template <int dim> 
void Step6<dim>::output_results(const unsigned int cycle) const 
{ 
  { 
    GridOut               grid_out; 
    std::ofstream         output("grid-" + std::to_string(cycle) + ".gnuplot"); 
    GridOutFlags::Gnuplot gnuplot_flags(false, 5); 
    grid_out.set_flags(gnuplot_flags); 
    MappingQGeneric<dim> mapping(3); 
    grid_out.write_gnuplot(triangulation, output, &mapping); 
  } 

  { 
    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "solution"); 
    data_out.build_patches(); 

    std::ofstream output("solution-" + std::to_string(cycle) + ".vtu"); 
    data_out.write_vtu(output); 
  } 
} 
// @sect4{Step6::run}  

//  <code>main()</code> 之前的最后一个函数又是该类的主要驱动，  <code>run()</code>  。它与  step-5  的函数类似，只是我们在程序中再次生成一个文件，而不是从磁盘中读取，我们自适应地而不是全局地细化网格，并且我们在本函数中输出最终网格上的解决方案。

// 该函数主循环的第一个块是处理网格生成。如果这是该程序的第一个循环，我们现在不是像上一个例子那样从磁盘上的文件中读取网格，而是再次使用库函数来创建它。域还是一个圆，中心在原点，半径为1（这是函数的两个隐藏参数，有默认值）。

// 你会注意到粗略的网格比我们在前面的例子中从文件中读出的网格质量要差：单元格的形成不太平均。然而，使用库函数，这个程序在任何空间维度上都可以工作，而以前不是这样的。

// 如果我们发现这不是第一个周期，我们要细化网格。与上一个例子程序中采用的全局细化不同，我们现在使用上述的自适应程序。

// 循环的其余部分看起来和以前一样。

template <int dim> 
void Step6<dim>::run() 
{ 
  for (unsigned int cycle = 0; cycle < 8; ++cycle) 
    { 
      std::cout << "Cycle " << cycle << ':' << std::endl; 

      if (cycle == 0) 
        { 
          GridGenerator::hyper_ball(triangulation); 
          triangulation.refine_global(1); 
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
// @sect3{The <code>main</code> function}  

// 主函数的功能与之前的例子相比没有改变，但我们采取了额外的谨慎措施。有时，会出现一些问题（比如写输出文件时磁盘空间不足，试图分配向量或矩阵时内存不足，或者由于某种原因我们无法从文件中读取或写入文件），在这些情况下，库会抛出异常。由于这些是运行时的问题，而不是可以一劳永逸的编程错误，这种异常在优化模式下不会被关闭，与我们用来测试编程错误的 <code>Assert</code> 宏相反。如果没有被捕获，这些异常会传播到 <code>main</code> 函数的调用树上，如果它们在那里也没有被捕获，程序就会被中止。在很多情况下，比如内存或磁盘空间不足，我们什么也做不了，但我们至少可以打印一些文字，试图解释程序失败的原因。下面显示了一种方法。以这种方式编写任何较大的程序当然是有用的，你可以通过或多或少地复制这个函数来做到这一点，但 <code>try</code> 块除外，它实际上编码了本应用程序所特有的功能。

int main() 
{ 

// 这个函数布局的总体思路如下：让我们试着像以前那样运行程序......

  try 
    { 
      Step6<2> laplace_problem_2d; 
      laplace_problem_2d.run(); 
    } 

// ......如果这应该是失败的，尽量收集尽可能多的信息。具体来说，如果被抛出的异常是一个从C++标准类派生出来的对象  <code>exception</code>, then we can use the <code>what</code>  成员函数，以获得一个描述异常被抛出原因的字符串。

// deal.II的异常类都是从标准类派生出来的，特别是 <code>exc.what()</code> 函数将返回与使用 <code>Assert</code> 宏抛出的异常所产生的字符串大致相同。在前面的例子中，你已经看到了这种异常的输出，然后你知道它包含了异常发生的文件和行号，以及其他一些信息。这也是下面的语句会打印的内容。

// 除此以外，除了用错误代码退出程序（这就是 <code>return 1;</code> 的作用），我们能做的并不多。

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

// 如果在某处抛出的异常不是从标准 <code>exception</code> 类派生出来的对象，那么我们根本无法做任何事情。那么我们就简单地打印一个错误信息并退出。

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

// 如果我们走到这一步，就没有任何异常传播到主函数上（可能有异常，但它们在程序或库的某个地方被捕获）。因此，程序按预期执行，我们可以无误返回。

  return 0; 
} 


