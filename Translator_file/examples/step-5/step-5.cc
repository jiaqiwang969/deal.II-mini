

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 1999 - 2021 by the deal.II authors 
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
 * Author: Wolfgang Bangerth, University of Heidelberg, 1999 
 */ 


// @sect3{Include files}  

// 同样，前几个include文件已经知道了，所以我们不会对它们进行评论。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 

// 这个是新的。我们想从磁盘上读取一个三角图，做这个的类在下面的文件中声明。

#include <deal.II/grid/grid_in.h> 

// 我们将使用一个圆形域，而描述其边界的对象来自这个文件。

#include <deal.II/grid/manifold_lib.h> 

// 这是C++ ...

#include <fstream> 
#include <iostream> 

// 最后，这在以前的教程程序中已经讨论过了。

using namespace dealii; 
// @sect3{The <code>Step5</code> class template}  

// 主类大部分和前面的例子一样。最明显的变化是删除了 <code>make_grid</code> 函数，因为现在创建网格是在 <code>run</code> 函数中完成的，其余功能都在 <code>setup_system</code> 中。除此以外，一切都和以前一样。

template <int dim> 
class Step5 
{ 
public: 
  Step5(); 
  void run(); 

private: 
  void setup_system(); 
  void assemble_system(); 
  void solve(); 
  void output_results(const unsigned int cycle) const; 

  Triangulation<dim> triangulation; 
  FE_Q<dim>          fe; 
  DoFHandler<dim>    dof_handler; 

  SparsityPattern      sparsity_pattern; 
  SparseMatrix<double> system_matrix; 

  Vector<double> solution; 
  Vector<double> system_rhs; 
}; 
// @sect3{Working with nonconstant coefficients}  

// 在  step-4  中，我们展示了如何使用非恒定边界值和右手边。 在这个例子中，我们想在椭圆算子中使用一个可变系数来代替。由于我们有一个只取决于空间中的点的函数，我们可以做得更简单一些，使用一个普通的函数而不是继承自Function。

// 这是对单点的系数函数的实现。如果与原点的距离小于0.5，我们让它返回20，否则返回1。

template <int dim> 
double coefficient(const Point<dim> &p) 
{ 
  if (p.square() < 0.5 * 0.5) 
    return 20; 
  else 
    return 1; 
} 
// @sect3{The <code>Step5</code> class implementation}  
// @sect4{Step5::Step5}  

// 这个函数和以前一样。

template <int dim> 
Step5<dim>::Step5() 
  : fe(1) 
  , dof_handler(triangulation) 
{} 

//  @sect4{Step5::setup_system}  

// 这是前面例子中的函数 <code>make_grid</code> ，减去了生成网格的部分。其他一切都没有变化。

template <int dim> 
void Step5<dim>::setup_system() 
{ 
  dof_handler.distribute_dofs(fe); 

  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
            << std::endl; 

  DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
  DoFTools::make_sparsity_pattern(dof_handler, dsp); 
  sparsity_pattern.copy_from(dsp); 

  system_matrix.reinit(sparsity_pattern); 

  solution.reinit(dof_handler.n_dofs()); 
  system_rhs.reinit(dof_handler.n_dofs()); 
} 

//  @sect4{Step5::assemble_system}  

// 和前面的例子一样，这个函数在功能上没有太大变化，但仍有一些优化，我们将展示这些优化。对此，需要注意的是，如果使用高效的求解器（如预设条件的CG方法），组装矩阵和右手边会花费相当的时间，你应该考虑在某些地方使用一到两个优化。

// 该函数的前几部分与之前完全没有变化。

template <int dim> 
void Step5<dim>::assemble_system() 
{ 
  QGauss<dim> quadrature_formula(fe.degree + 1); 

  FEValues<dim> fe_values(fe, 
                          quadrature_formula, 
                          update_values | update_gradients | 
                            update_quadrature_points | update_JxW_values); 

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
  Vector<double>     cell_rhs(dofs_per_cell); 

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

// 接下来是对所有单元的典型循环，以计算局部贡献，然后将它们转移到全局矩阵和向量中。与 step-4 相比，这部分的唯一变化是我们将使用上面定义的 <code>coefficient()</code> 函数来计算每个正交点的系数值。

  for (const auto &cell : dof_handler.active_cell_iterators()) 
    { 
      cell_matrix = 0.; 
      cell_rhs    = 0.; 

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

              cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q) 
                              1.0 *                               // f(x_q) 
                              fe_values.JxW(q_index));            // dx 
            } 
        } 

      cell->get_dof_indices(local_dof_indices); 
      for (const unsigned int i : fe_values.dof_indices()) 
        { 
          for (const unsigned int j : fe_values.dof_indices()) 
            system_matrix.add(local_dof_indices[i], 
                              local_dof_indices[j], 
                              cell_matrix(i, j)); 

          system_rhs(local_dof_indices[i]) += cell_rhs(i); 
        } 
    } 

// 有了这样构建的矩阵，我们再次使用零边界值。

  std::map<types::global_dof_index, double> boundary_values; 
  VectorTools::interpolate_boundary_values(dof_handler, 
                                           0, 
                                           Functions::ZeroFunction<dim>(), 
                                           boundary_values); 
  MatrixTools::apply_boundary_values(boundary_values, 
                                     system_matrix, 
                                     solution, 
                                     system_rhs); 
} 
// @sect4{Step5::solve}  

// 求解过程看起来又和前面的例子差不多。然而，我们现在将使用一个预设条件的共轭梯度算法。做出这种改变并不难。事实上，我们唯一需要改变的是，我们需要一个作为预处理程序的对象。我们将使用SSOR（对称连续过度放松），放松系数为1.2。为此， <code>SparseMatrix</code> 类有一个函数可以做一个SSOR步骤，我们需要把这个函数的地址和它应该作用的矩阵（也就是要反转的矩阵）以及松弛因子打包成一个对象。 <code>PreconditionSSOR</code> 类为我们做了这个。(  <code>PreconditionSSOR</code>  类需要一个模板参数，表示它应该工作的矩阵类型。默认值是 <code>SparseMatrix@<double@></code> ，这正是我们在这里需要的，所以我们只需坚持使用默认值，不在角括号中指定任何东西。)

// 请注意，在目前的情况下，SSOR的表现并不比其他大多数预处理程序好多少（尽管比没有预处理好）。在下一个教程程序  step-6  的结果部分，将对不同的预处理进行简要比较。

// 有了这个，函数的其余部分就很简单了：我们现在使用我们声明的预处理程序，而不是之前创建的 <code>PreconditionIdentity</code> 对象，CG求解器将为我们完成其余的工作。

template <int dim> 
void Step5<dim>::solve() 
{ 
  SolverControl            solver_control(1000, 1e-12); 
  SolverCG<Vector<double>> solver(solver_control); 

  PreconditionSSOR<SparseMatrix<double>> preconditioner; 
  preconditioner.initialize(system_matrix, 1.2); 

  solver.solve(system_matrix, solution, system_rhs, preconditioner); 

  std::cout << "   " << solver_control.last_step() 
            << " CG iterations needed to obtain convergence." << std::endl; 
} 
// @sect4{Step5::output_results and setting output flags}  

// 将输出写入文件的方法与上一个教程中的基本相同。唯一不同的是，我们现在需要为每个细化周期构建一个不同的文件名。

// 这个函数以VTU格式写入输出，这是VTK格式的一个变种，因为它压缩了数据，所以需要更少的磁盘空间。当然，如果你希望使用一个不理解VTK或VTU的可视化程序，DataOut类还支持许多其他格式。

template <int dim> 
void Step5<dim>::output_results(const unsigned int cycle) const 
{ 
  DataOut<dim> data_out; 

  data_out.attach_dof_handler(dof_handler); 
  data_out.add_data_vector(solution, "solution"); 

  data_out.build_patches(); 

  std::ofstream output("solution-" + std::to_string(cycle) + ".vtu"); 
  data_out.write_vtu(output); 
} 

//  @sect4{Step5::run}  

// 在这个程序中，倒数第二件事是对 <code>run()</code> 函数的定义。与之前的程序不同，我们将在一个网格序列上进行计算，在每次迭代后都会进行全局细化。因此，该函数由6个周期的循环组成。在每个循环中，我们首先打印循环编号，然后决定如何处理网格。如果这不是第一个周期，我们就简单地对现有的网格进行一次全局精炼。然而，在运行这些循环之前，我们必须先生成一个网格。

// 在前面的例子中，我们已经使用了 <code>GridGenerator</code> 类中的一些函数。在这里，我们想从一个存储单元的文件中读取网格，这个文件可能来自其他人，也可能是一个网格生成工具的产物。

// 为了从文件中读取网格，我们生成一个数据类型为GridIn的对象，并将三角剖分与之相关联（也就是说，当我们要求它读取文件时，我们告诉它要填充我们的三角剖分对象）。然后我们打开相应的文件，用文件中的数据初始化三角剖分。

template <int dim> 
void Step5<dim>::run() 
{ 
  GridIn<dim> grid_in; 
  grid_in.attach_triangulation(triangulation); 
  std::ifstream input_file("circle-grid.inp"); 

// 我们现在想读取该文件。但是，输入文件只针对二维三角图，而这个函数是一个任意维度的模板。由于这只是一个演示程序，我们不会为不同的维度使用不同的输入文件，而是在不在二维的情况下迅速杀死整个程序。当然，由于下面的主函数假定我们是在二维空间工作，我们可以跳过这个检查，在这个版本的程序中，不会有任何不良影响。

// 事实证明，90%以上的编程错误都是无效的函数参数，如无效的数组大小等，所以我们在整个deal.II中大量使用断言来捕捉此类错误。为此， <code>Assert</code> 宏是一个很好的选择，因为它确保作为第一个参数的条件是有效的，如果不是，就抛出一个异常（它的第二个参数），通常会终止程序，并给出错误发生的位置和原因的信息。关于 @p Assert 宏的具体作用，可以在 @ref Exceptions "异常文档模块 "中找到更详细的讨论）。这通常会大大减少发现编程错误的时间，我们发现断言是快速编程的宝贵手段。

// 另一方面，如果你想做大的计算，所有这些检查（目前库中有超过10000个）不应该使程序太慢。为此， <code>Assert</code> 宏只在调试模式下使用，如果在优化模式下则扩展为零。因此，当你在小问题上测试你的程序并进行调试时，断言会告诉你问题出在哪里。一旦你的程序稳定了，你可以关闭调试，程序将在没有断言的情况下以最大速度运行你的实际计算。更准确地说：通过在优化模式下编译你的程序，关闭库中的所有检查（这些检查可以防止你用错误的参数调用函数，从数组中走出来，等等），通常可以使程序的运行速度提高四倍左右。即使优化后的程序性能更高，我们仍然建议在调试模式下开发，因为它允许库自动发现许多常见的编程错误。对于那些想尝试的人来说。从调试模式切换到优化模式的方法是用<code>make release</code>命令重新编译你的程序。现在 <code>make</code> 程序的输出应该向你表明，该程序现在是以优化模式编译的，以后也会被链接到已经为优化模式编译的库。为了切换回调试模式，只需用  <code>make debug</code>  命令重新编译。

  Assert(dim == 2, ExcInternalError()); 

// ExcInternalError是一个全局定义的异常，每当出现严重的错误时就会抛出。通常，人们希望使用更具体的异常，特别是在这种情况下，如果 <code>dim</code> 不等于2，人们当然会尝试做其他事情，例如使用库函数创建一个网格。终止程序通常不是一个好主意，断言实际上只应该用于不应该发生的特殊情况，但由于程序员、用户或其他人的愚蠢而可能发生。上面的情况并不是对Assert的巧妙使用，但是再次强调：这是一个教程，也许值得展示一下什么是不应该做的，毕竟。

// 所以，如果我们通过了断言，我们就知道dim==2，现在我们就可以真正地读取网格。它的格式是UCD（非结构化单元数据）（尽管惯例是使用UCD文件的后缀 <code>inp</code> ）。

  grid_in.read_ucd(input_file); 

// 如果你想使用其他输入格式，你必须使用其他 <code>grid_in.read_xxx</code> 函数之一。(参见  <code>GridIn</code>  类的文档，以了解目前支持哪些输入格式)。

// 文件中的网格描述了一个圆。因此，我们必须使用一个流形对象，告诉三角计算在细化网格时将边界上的新点放在哪里。与 step-1 不同的是，由于GridIn不知道域的边界是圆形的（与 GridGenerator::hyper_shell) 不同的是，我们必须在创建三角网格后明确地将流形附加到边界上，以便在细化网格时获得正确的结果。

  const SphericalManifold<dim> boundary; 
  triangulation.set_all_manifold_ids_on_boundary(0); 
  triangulation.set_manifold(0, boundary); 

  for (unsigned int cycle = 0; cycle < 6; ++cycle) 
    { 
      std::cout << "Cycle " << cycle << ':' << std::endl; 

      if (cycle != 0) 
        triangulation.refine_global(1); 

// 现在我们有了一个确定的网格，我们写一些输出，做所有我们在前面的例子中已经看到的事情。

      std::cout << "   Number of active cells: "  // 
                << triangulation.n_active_cells() // 
                << std::endl                      // 
                << "   Total number of cells: "   // 
                << triangulation.n_cells()        // 
                << std::endl; 

      setup_system(); 
      assemble_system(); 
      solve(); 
      output_results(cycle); 
    } 
} 
// @sect3{The <code>main</code> function}  

// 主函数看起来和前面的例子中的函数差不多，所以我们就不进一步评论了。

int main() 
{ 
  Step5<2> laplace_problem_2d; 
  laplace_problem_2d.run(); 
  return 0; 
} 

