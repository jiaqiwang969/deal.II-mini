CCTest_file/step-4.cc

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

// 前面几个（很多）include文件已经在前面的例子中使用过了，所以我们在这里不再解释它们的含义。

#include <deal.II/grid/tria.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 

#include <deal.II/numerics/data_out.h> 
#include <fstream> 
#include <iostream> 

// 这是新的，但是：在前面的例子中，我们从线性求解器得到了一些不需要的输出。如果我们想抑制它，我们必须包括这个文件，并在程序的某个地方添加一行字（见下面的main()函数）。

#include <deal.II/base/logstream.h> 

// 最后一步，和以前的程序一样，是将所有deal.II的类和函数名导入全局命名空间中。

using namespace dealii; 
// @sect3{The <code>Step4</code> class template}  

// 这又是前面例子中的 <code>Step4</code> 类。唯一不同的是，我们现在把它声明为一个带有模板参数的类，而模板参数当然是我们要解决拉普拉斯方程的空间维度。当然，几个成员变量也取决于这个维度，特别是Triangulation类，它必须分别表示四边形或六面体。除此以外，一切都和以前一样。

template <int dim> 
class Step4 
{ 
public: 
  Step4(); 
  void run(); 

private: 
  void make_grid(); 
  void setup_system(); 
  void assemble_system(); 
  void solve(); 
  void output_results() const; 

  Triangulation<dim> triangulation; 
  FE_Q<dim>          fe; 
  DoFHandler<dim>    dof_handler; 

  SparsityPattern      sparsity_pattern; 
  SparseMatrix<double> system_matrix; 

  Vector<double> solution; 
  Vector<double> system_rhs; 
}; 
// @sect3{Right hand side and boundary values}  

// 在下文中，我们又声明了两个类，表示右手边和非均质的Dirichlet边界值。两者都是一个二维空间变量的函数，所以我们也将它们声明为模板。

// 这些类中的每一个都是从一个共同的、抽象的基类Function派生出来的，它声明了所有函数都必须遵循的共同接口。特别是，具体的类必须重载 <code>value</code> 函数，该函数接收二维空间中的一个点作为参数，并将该点的值作为 <code>double</code> 变量返回。

//  <code>value</code> 函数需要第二个参数，我们在这里将其命名为 <code>component</code>  : 这只适用于矢量值函数，你可能想访问点 <code>p</code> 处的矢量的某个分量。然而，我们的函数是标量的，所以我们不需要担心这个参数，在函数的实现中也不会使用它。在库的头文件中，Function基类对 <code>value</code> 函数的声明中，分量的默认值为0，所以我们在访问右侧的 <code>value</code> 函数时，只需要一个参数，即我们要评估函数的点。然后，对于标量函数，可以简单地省略分量的值。

// 函数对象在库中很多地方都有使用（例如，在 step-3 中我们使用了一个 Functions::ZeroFunction 实例作为 VectorTools::interpolate_boundary_values) 的参数，这是我们定义一个继承自Function的新类的第一个教程。由于我们只调用 Function::value(), ，我们可以只用一个普通的函数（这就是 step-5 中的做法），但由于这是一个教程，为了举例说明，我们继承了Function。

template <int dim> 
class RightHandSide : public Function<dim> 
{ 
public: 
  virtual double value(const Point<dim> & p, 
                       const unsigned int component = 0) const override; 
}; 

template <int dim> 
class BoundaryValues : public Function<dim> 
{ 
public: 
  virtual double value(const Point<dim> & p, 
                       const unsigned int component = 0) const override; 
}; 

// 如果你不熟悉上述函数声明中的关键字 "virtual "和 "override "是什么意思，你可能会想看看你最喜欢的C++书籍或在线教程，如http:www.cplusplus.com/doc/tutorial/polymorphism/ 。从本质上讲，这里发生的事情是Function<dim>是一个 "抽象 "基类，它声明了某种 "接口"--一组可以在这类对象上调用的函数。但它实际上并没有*实现*这些函数：它只是说 "Function对象是这样的"，但它实际上是什么样的函数，则留给实现了`value()`函数的派生类。

// 从另一个类中派生出一个类，通常称为 "is-a "关系函数。在这里，`RightHandSide`类 "是一个 "函数类，因为它实现了Function基类所描述的接口。("value() "函数的实际实现在下面的代码块中)。那么`virtual`关键字意味着 "是的，这里的函数可以被派生类覆盖"，而`override`关键字意味着 "是的，这实际上是一个我们知道已经被声明为基类一部分的函数"。覆盖 "关键字不是严格必要的，但它是防止打字错误的一个保险。如果我们把函数的名字或一个参数的类型弄错了，编译器会警告我们说："你说这个函数覆盖了基类中的一个函数，但实际上我不知道有任何这样的函数有这个名字和这些参数。"

// 但回到这里的具体案例。在本教程中，我们选择2D中的函数 $4(x^4+y^4)$ ，或者3D中的 $4(x^4+y^4+z^4)$ 作为右手边。我们可以用空间维度上的if语句来写这个区别，但这里有一个简单的方法，通过使用一个短循环，也允许我们在一维（或四维，如果你想这样做）中使用相同的函数。 幸运的是，编译器在编译时就知道循环的大小（记住，在你定义模板时，编译器不知道 <code>dim</code> 的值，但当它后来遇到语句或声明 <code>RightHandSide@<2@></code> 时，它将采取模板，用2替换所有出现的dim，并编译出结果函数）。 换句话说，在编译这个函数的时候，主体将被执行的次数是已知的，编译器可以将循环所需的开销降到最低；结果将和我们马上使用上面的公式一样快。

// 最后要注意的是， <code>Point@<dim@></code> 表示二维空间中的一个点，它的各个组成部分（即 $x$ 、 $y$ 、...坐标）可以像C和C++中一样用（）运算符访问（事实上，[]运算符也同样有效），索引从0开始。

template <int dim> 
double RightHandSide<dim>::value(const Point<dim> &p, 
                                 const unsigned int /*component*/) const 
{ 
  double return_value = 0.0; 
  for (unsigned int i = 0; i < dim; ++i) 
    return_value += 4.0 * std::pow(p(i), 4.0); 

  return return_value; 
} 

// 作为边界值，我们选择二维的 $x^2+y^2$ ，三维的 $x^2+y^2+z^2$ 。这恰好等于从原点到我们想评估函数的点的矢量的平方，而不考虑维度。所以这就是我们的返回值。

template <int dim> 
double BoundaryValues<dim>::value(const Point<dim> &p, 
                                  const unsigned int /*component*/) const 
{ 
  return p.square(); 
} 

//  @sect3{Implementation of the <code>Step4</code> class}  

// 接下来是利用上述函数的类模板的实现。和以前一样，我们将把所有东西写成模板，这些模板有一个形式参数 <code>dim</code> ，在我们定义模板函数时，我们假设这个参数是未知的。只有在以后，编译器才会发现 <code>Step4@<2@></code> (in the <code>main</code> 函数的声明，实际上），并在编译整个类时将 <code>dim</code> 替换成2，这个过程被称为 "模板的实例化"。这样做的时候，它也会用 <code>RightHandSide@<dim@></code> 的实例替换 <code>RightHandSide@<2@></code> ，并从类模板中实例化后一个类。

// 事实上，编译器也会在 <code>main()</code> 中找到一个 <code>Step4@<3@></code> 声明。这将导致它再次回到一般的 <code>Step4@<dim@></code> 模板，替换所有出现的 <code>dim</code> ，这次是3，并第二次编译这个类。注意这两个实例  <code>Step4@<2@></code>  和  <code>Step4@<3@></code>  是完全独立的类；它们唯一的共同特征是它们都是从同一个通用模板中实例化出来的，但是它们不能相互转换，例如，它们没有共享代码（两个实例都是完全独立编译的）。

//  @sect4{Step4::Step4}  

// 在这个介绍之后，这里是  <code>Step4</code>  类的构造函数。它指定了所需的有限元素的多项式程度，并将DoFHandler与三角形关联起来，就像在前面的例子程序中一样，  step-3  。

template <int dim> 
Step4<dim>::Step4() 
  : fe(1) 
  , dof_handler(triangulation) 
{} 
// @sect4{Step4::make_grid}  

// 网格的创建在本质上是与维度有关的东西。然而，只要领域在二维或三维中足够相似，库就可以为你抽象。在我们的例子中，我们想再次在二维的正方形 $[-1,1]\times [-1,1]$ 上求解，或者在三维的立方体 $[-1,1] \times [-1,1] \times [-1,1]$ 上求解；两者都可以被称为 GridGenerator::hyper_cube(), ，因此我们可以在任何维度上使用同一个函数。当然，在二维和三维中创建超立方体的函数有很大的不同，但这是你不需要关心的事情。让库来处理这些困难的事情。

template <int dim> 
void Step4<dim>::make_grid() 
{ 
  GridGenerator::hyper_cube(triangulation, -1, 1); 
  triangulation.refine_global(4); 

  std::cout << "   Number of active cells: " << triangulation.n_active_cells() 
            << std::endl 
            << "   Total number of cells: " << triangulation.n_cells() 
            << std::endl; 
} 
// @sect4{Step4::setup_system}  

// 这个函数看起来和前面的例子完全一样，尽管它执行的动作在细节上有很大的不同，如果 <code>dim</code> 刚好是3。从用户的角度来看，唯一显著的区别是所产生的单元格数量，在三个空间维度中比两个空间维度中要高得多

template <int dim> 
void Step4<dim>::setup_system() 
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
// @sect4{Step4::assemble_system}  

// 与前面的例子不同，我们现在想使用一个非恒定的右侧函数和非零边界值。这两个任务都是很容易实现的，只需在矩阵和右手边的组合中增加几行代码即可。

// 更有趣的是，我们将矩阵和右手边的向量维度独立组装起来的方式：与二维的情况根本没有区别。由于这个函数中使用的重要对象（正交公式、FEValues）也通过模板参数的方式依赖于维度，它们可以为这个函数所编译的维度正确设置一切。通过使用模板参数声明所有可能依赖于维度的类，库可以为你完成几乎所有的工作，你不需要关心大多数事情。

template <int dim> 
void Step4<dim>::assemble_system() 
{ 
  QGauss<dim> quadrature_formula(fe.degree + 1); 

// 我们希望有一个非恒定的右手，所以我们使用上面声明的类的一个对象来生成必要的数据。由于这个右侧对象只在本函数中局部使用，所以我们在这里把它声明为一个局部变量。

  RightHandSide<dim> right_hand_side; 

// 与之前的例子相比，为了评估非恒定右手函数，我们现在还需要我们目前所在单元上的正交点（之前，我们只需要FEValues对象中的形状函数的值和梯度，以及正交权重， FEValues::JxW() ）。我们可以通过给FEValues对象添加#update_quadrature_points标志来让它为我们做事。

  FEValues<dim> fe_values(fe, 
                          quadrature_formula, 
                          update_values | update_gradients | 
                            update_quadrature_points | update_JxW_values); 

// 然后我们再次定义与前面程序中相同的缩写。这个变量的值当然取决于我们现在使用的维度，但是FiniteElement类为你做了所有必要的工作，你不需要关心与维度有关的部分。

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
  Vector<double>     cell_rhs(dofs_per_cell); 

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

// 接下来，我们又要在所有的单元格上进行循环，并汇集局部贡献。 请注意，一个单元在两个空间维度上是一个四边形，但在三维上是一个六面体。事实上， <code>active_cell_iterator</code> 的数据类型是不同的，这取决于我们所处的维度，但对外界来说，它们看起来是一样的，你可能永远不会看到区别。在任何情况下，真正的类型是通过使用`auto`来隐藏的。

  for (const auto &cell : dof_handler.active_cell_iterators()) 
    { 
      fe_values.reinit(cell); 
      cell_matrix = 0; 
      cell_rhs    = 0; 

// 现在我们要把本地矩阵和右手边组合起来。这个过程和前面的例子完全一样，但是现在我们重新调整循环的顺序（我们可以安全地这样做，因为它们是相互独立的），并尽可能地合并本地矩阵和本地向量的循环，使事情变得更快。

// 组装右手边与我们在 step-3 中的做法有唯一的区别：我们没有使用值为1的常数右手边，而是使用代表右手边的对象并在正交点对其进行评估。

      for (const unsigned int q_index : fe_values.quadrature_point_indices()) 
        for (const unsigned int i : fe_values.dof_indices()) 
          { 
            for (const unsigned int j : fe_values.dof_indices()) 
              cell_matrix(i, j) += 
                (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q) 
                 fe_values.shape_grad(j, q_index) * // grad phi_j(x_q) 
                 fe_values.JxW(q_index));           // dx 

            const auto &x_q = fe_values.quadrature_point(q_index); 
            cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q) 
                            right_hand_side.value(x_q) *        // f(x_q) 
                            fe_values.JxW(q_index));            // dx 
          } 

// 作为对这些循环的最后说明：当我们将局部贡献集合到 <code>cell_matrix(i,j)</code> 时，我们必须将形状函数 $i$ 和 $j$ 在点号q_index的梯度相乘并与标量权重JxW相乘。这就是实际发生的情况。  <code>fe_values.shape_grad(i,q_index)</code> 返回一个 <code>dim</code> 维向量，由 <code>Tensor@<1,dim@></code> 对象表示，将其与 <code>fe_values.shape_grad(j,q_index)</code> 的结果相乘的运算器*确保两个向量的 <code>dim</code> 分量被适当收缩，结果是一个标量浮点数，然后与权重相乘。在内部，这个操作符*确保对向量的所有 <code>dim</code> 分量都能正确发生，无论 <code>dim</code> 是2、3还是其他空间维度；从用户的角度来看，这并不值得费心，然而，如果想独立编写代码维度，事情就会简单很多。

// 随着本地系统的组装，转移到全局矩阵和右手边的工作与之前完全一样，但在这里我们再次合并了一些循环以提高效率。

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

// 作为这个函数的最后一步，我们希望在这个例子中拥有非均质的边界值，与之前的例子不同。这是一个简单的任务，我们只需要用一个描述我们想使用的边界值的类的对象（即上面声明的 <code>BoundaryValues</code> 类）来替换那里使用的 Functions::ZeroFunction 。

// 函数 VectorTools::interpolate_boundary_values() 只对标有边界指标0的面起作用（因为我们在下面的第二个参数中说该函数应该对其起作用）。如果有的面的边界指标不是0，那么函数interpolate_boundary_values将对这些面不起作用。对于拉普拉斯方程来说，什么都不做相当于假设在边界的这些部分，零诺伊曼边界条件成立。

  std::map<types::global_dof_index, double> boundary_values; 
  VectorTools::interpolate_boundary_values(dof_handler, 
                                           0, 
                                           BoundaryValues<dim>(), 
                                           boundary_values); 
  MatrixTools::apply_boundary_values(boundary_values, 
                                     system_matrix, 
                                     solution, 
                                     system_rhs); 
} 
// @sect4{Step4::solve}  

// 解决线性方程组是在大多数程序中看起来几乎相同的事情。特别是，它与维度无关，所以这个函数是从前面的例子中逐字复制的。

template <int dim> 
void Step4<dim>::solve() 
{ 
  SolverControl            solver_control(1000, 1e-12); 
  SolverCG<Vector<double>> solver(solver_control); 
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity()); 

// 不过我们做了一个补充：由于我们抑制了线性求解器的输出，我们必须手工打印迭代次数。

  std::cout << "   " << solver_control.last_step() 
            << " CG iterations needed to obtain convergence." << std::endl; 
} 
// @sect4{Step4::output_results}  

// 这个函数也做了  step-3  中各自的工作。这里也没有改变维度的独立性。

// 由于程序将同时运行拉普拉斯求解器的2D和3D版本，我们使用文件名中的维度为每次运行生成不同的文件名（在一个更好的程序中，我们将检查 <code>dim</code> 是否可以有2或3以外的其他值，但为了简洁起见，我们在这里忽略了这一点）。

template <int dim> 
void Step4<dim>::output_results() const 
{ 
  DataOut<dim> data_out; 

  data_out.attach_dof_handler(dof_handler); 
  data_out.add_data_vector(solution, "solution"); 

  data_out.build_patches(); 

  std::ofstream output(dim == 2 ? "solution-2d.vtk" : "solution-3d.vtk"); 
  data_out.write_vtk(output); 
} 

//  @sect4{Step4::run}  

// 这是一个对所有事情都有最高级别控制的函数。除了一行额外的输出外，它与前面的例子相同。

template <int dim> 
void Step4<dim>::run() 
{ 
  std::cout << "Solving problem in " << dim << " space dimensions." 
            << std::endl; 

  make_grid(); 
  setup_system(); 
  assemble_system(); 
  solve(); 
  output_results(); 
} 
// @sect3{The <code>main</code> function}  

// 这是主函数。它看起来也大多像 step-3 中的内容，但如果你看下面的代码，注意我们是如何首先创建一个 <code>Step4@<2@></code> 类型的变量（迫使编译器用 <code>dim</code> replaced by <code>2</code> 编译类模板）并运行一个2d模拟，然后我们用3d做整个事情。

// 在实践中，这可能不是你经常做的事情（你可能要么想解决一个2D的问题，要么想解决一个3D的问题，但不会同时解决这两个问题）。然而，它展示了一种机制，我们可以在一个地方简单地改变我们想要的维度，从而迫使编译器为我们要求的维度重新编译独立的类模板。这里的重点在于，我们只需要改变一个地方。这使得在计算速度较快的2D环境下调试程序变得非常简单，然后将一个地方切换到3，在3D环境下运行计算量大得多的程序，进行 "真实 "的计算。

// 这两个区块中的每一个都用大括号括起来，以确保 <code>laplace_problem_2d</code> 这个变量在我们继续为3D情况分配内存之前就已经超出了范围（并释放了它所持有的内存）。如果没有额外的大括号， <code>laplace_problem_2d</code> 变量只会在函数结束时被销毁，也就是在运行完3d问题后被销毁，而且会在3d运行时不必要地占用内存，而实际使用它。

int main() 
{ 
  { 
    Step4<2> laplace_problem_2d; 
    laplace_problem_2d.run(); 
  } 

  { 
    Step4<3> laplace_problem_3d; 
    laplace_problem_3d.run(); 
  } 

  return 0; 
} 


