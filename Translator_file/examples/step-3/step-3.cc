CCTest_file/step-3.cc

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
 * Authors: Wolfgang Bangerth, 1999, 
 *          Guido Kanschat, 2011 
 */ 


// @sect3{Many new include files}  

// 这些包含文件已经为你所知。它们声明了处理三角形和自由度枚举的类。

#include <deal.II/grid/tria.h> 
#include <deal.II/dofs/dof_handler.h> 

// 在这个文件中声明了创建网格的函数。

#include <deal.II/grid/grid_generator.h> 

// 这个文件包含了对拉格朗日插值有限元的描述。

#include <deal.II/fe/fe_q.h> 

// 而这个文件是创建稀疏矩阵的稀疏模式所需要的，如前面的例子中所示。

#include <deal.II/dofs/dof_tools.h> 

// 接下来的两个文件是在每个单元上使用正交法组装矩阵所需要的。下面将对其中声明的类进行解释。

#include <deal.II/fe/fe_values.h> 
#include <deal.II/base/quadrature_lib.h> 

// 以下是我们在处理边界值时需要的三个包含文件。

#include <deal.II/base/function.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 

// 我们现在几乎到了终点。第二组到最后一组include文件是用于线性代数的，我们用它来解决拉普拉斯方程的有限元离散化所产生的方程组。我们将使用向量和全矩阵在每个单元中组装方程组，并将结果转移到稀疏矩阵中。然后我们将使用共轭梯度求解器来解决这个问题，为此我们需要一个预处理程序（在这个程序中，我们使用身份预处理程序，它没有任何作用，但我们还是需要包括这个文件）。

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 

// 最后，这是为了输出到文件和控制台。

#include <deal.II/numerics/data_out.h> 
#include <fstream> 
#include <iostream> 

// ...这是为了将deal.II命名空间导入到全局范围。

using namespace dealii; 
// @sect3{The <code>Step3</code> class}  

// 在这个程序中，我们没有采用以前例子中的程序化编程，而是将所有东西都封装到一个类中。这个类由一些函数组成，这些函数分别执行有限元程序的某些方面，一个`main`函数控制先做什么和后做什么，还有一个成员变量列表。

// 该类的公共部分相当简短：它有一个构造函数和一个从外部调用的函数`run`，其作用类似于`main`函数：它协调该类的哪些操作应以何种顺序运行。该类中的其他东西，即所有真正做事情的函数，都在该类的私有部分。

class Step3 
{ 
public: 
  Step3(); 

  void run(); 

// 然后，还有一些成员函数，它们主要是做它们名字所暗示的事情，在介绍中已经讨论过了。由于它们不需要从外部调用，所以它们是本类的私有函数。

private: 
  void make_grid(); 
  void setup_system(); 
  void assemble_system(); 
  void solve(); 
  void output_results() const; 

// 最后我们还有一些成员变量。有一些变量描述了三角形和自由度的全局编号（我们将在这个类的构造函数中指定有限元的确切多项式程度）...

  Triangulation<2> triangulation; 
  FE_Q<2>          fe; 
  DoFHandler<2>    dof_handler; 

// ...拉普拉斯方程离散化产生的系统矩阵的稀疏模式和数值的变量...

  SparsityPattern      sparsity_pattern; 
  SparseMatrix<double> system_matrix; 

// .......以及用于保存右手边和解决方案向量的变量。

  Vector<double> solution; 
  Vector<double> system_rhs; 
}; 
// @sect4{Step3::Step3}  

// 这里是构造函数。它除了首先指定我们需要双线性元素（由有限元对象的参数表示，它表示多项式的程度），并将dof_handler变量与我们使用的三角形相关联之外，没有做更多的工作。(注意，目前三角结构并没有设置网格，但是DoFHandler并不关心：它只想知道它将与哪个三角结构相关联，只有当你使用distribution_dofs()函数试图在网格上分布自由度时，它才开始关心实际的网格。) Step3类的所有其他成员变量都有一个默认的构造函数，它可以完成我们想要的一切。

Step3::Step3() 
  : fe(1) 
  , dof_handler(triangulation) 
{} 
// @sect4{Step3::make_grid}  

// 现在，我们要做的第一件事是生成我们想在其上进行计算的三角形，并对每个顶点进行自由度编号。我们之前在 step-1 和 step-2 中分别看到过这两个步骤。

// 这个函数做的是第一部分，创建网格。 我们创建网格并对所有单元格进行五次细化。由于初始网格（也就是正方形 $[-1,1] \times [-1,1]$ ）只由一个单元组成，所以最终的网格有32乘以32个单元，总共是1024个。

// 不确定1024是否是正确的数字？我们可以通过使用三角形上的 <code>n_active_cells()</code> 函数输出单元格的数量来检查。

void Step3::make_grid() 
{ 
  GridGenerator::hyper_cube(triangulation, -1, 1); 
  triangulation.refine_global(5); 

  std::cout << "Number of active cells: " << triangulation.n_active_cells() 
            << std::endl; 
} 
// @note  我们调用 Triangulation::n_active_cells() 函数，而不是 Triangulation::n_cells(). 这里，<i>active</i>指的是没有进一步提炼的单元。我们强调 "活跃 "这个形容词，因为还有更多的单元，即最细的单元的父单元，它们的父单元等等，直到构成初始网格的一个单元为止。当然，在下一个更粗的层次上，单元格的数量是最细层次上的单元格的四分之一，即256，然后是64、16、4和1。如果你在上面的代码中调用 <code>triangulation.n_cells()</code> ，你会因此得到一个1365的值。另一方面，单元格的数量（相对于活动单元格的数量）通常没有什么意义，所以没有很好的理由去打印它。

//  @sect4{Step3::setup_system}  

// 接下来我们列举所有的自由度，并建立矩阵和向量对象来保存系统数据。枚举是通过使用 DoFHandler::distribute_dofs(), 来完成的，我们在 step-2 的例子中已经看到了。由于我们使用了FE_Q类，并且在构造函数中设置了多项式的度数为1，即双线性元素，这就将一个自由度与每个顶点联系起来。当我们在生成输出时，让我们也看看有多少自由度被生成。

void Step3::setup_system() 
{ 
  dof_handler.distribute_dofs(fe); 
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs() 
            << std::endl; 

// 每个顶点应该有一个DoF。因为我们有一个32乘以32的网格，所以DoFs的数量应该是33乘以33，即1089。

// 正如我们在前面的例子中所看到的，我们通过首先创建一个临时结构，标记那些可能为非零的条目，然后将数据复制到SparsityPattern对象中，然后可以被系统矩阵使用，来设置一个稀疏模式。

  DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
  DoFTools::make_sparsity_pattern(dof_handler, dsp); 
  sparsity_pattern.copy_from(dsp); 

// 注意，SparsityPattern对象并不保存矩阵的值，它只保存条目所在的位置。条目本身存储在SparseMatrix类型的对象中，我们的变量system_matrix就是其中之一。

// 稀疏模式和矩阵之间的区别是为了让几个矩阵使用相同的稀疏模式。这在这里似乎并不重要，但是当你考虑到矩阵的大小，以及建立稀疏模式可能需要一些时间时，如果你必须在程序中存储几个矩阵，这在大规模问题中就变得很重要了。

  system_matrix.reinit(sparsity_pattern); 

// 在这个函数中要做的最后一件事是将右侧向量和解向量的大小设置为正确的值。

  solution.reinit(dof_handler.n_dofs()); 
  system_rhs.reinit(dof_handler.n_dofs()); 
} 
// @sect4{Step3::assemble_system}  

// 下一步是计算形成线性系统的矩阵和右手边的条目，我们从中计算出解决方案。这是每一个有限元程序的核心功能，我们在介绍中已经讨论了主要步骤。

// 组装矩阵和向量的一般方法是在所有单元上循环，并在每个单元上通过正交计算该单元对全局矩阵和右侧的贡献。现在要认识到的一点是，我们需要实心单元上正交点位置的形状函数值。然而，有限元形状函数和正交点都只定义在参考单元上。因此，它们对我们帮助不大，事实上，我们几乎不会直接从这些对象中查询有关有限元形状函数或正交点的信息。

// 相反，我们需要的是一种将这些数据从参考单元映射到实际单元的方法。能够做到这一点的类都是由Mapping类派生出来的，尽管人们常常不必直接与它们打交道：库中的许多函数都可以将映射对象作为参数，但当它被省略时，它们只是简单地诉诸于标准的双线性Q1映射。我们将走这条路，暂时不打扰它（我们将在 step-10 、 step-11 和 step-12 中再讨论这个问题）。

// 所以我们现在有三个类的集合来处理：有限元、正交、和映射对象。这就太多了，所以有一种类型的类可以协调这三者之间的信息交流：FEValues类。如果给这三个对象各一个实例（或两个，以及一个隐式线性映射），它就能为你提供实心单元上正交点的形状函数值和梯度的信息。

// 利用所有这些，我们将把这个问题的线性系统组装在以下函数中。

void Step3::assemble_system() 
{ 

// 好的，我们开始吧：我们需要一个正交公式来计算每个单元格的积分。让我们采用一个高斯公式，每个方向有两个正交点，即总共有四个点，因为我们是在二维。这个正交公式可以准确地积分三度以下的多项式（在一维）。很容易检查出，这对目前的问题来说是足够的。

  QGauss<2> quadrature_formula(fe.degree + 1); 

// 然后我们初始化我们在上面简单谈及的对象。它需要被告知我们要使用哪个有限元，以及正交点和它们的权重（由一个正交对象共同描述）。如前所述，我们使用隐含的Q1映射，而不是自己明确指定一个。最后，我们必须告诉它我们希望它在每个单元上计算什么：我们需要正交点的形状函数值（对于右手 $(\varphi_i,f)$ ），它们的梯度（对于矩阵条目 $(\nabla \varphi_i, \nabla \varphi_j)$ ），以及正交点的权重和从参考单元到实际单元的雅各布变换的行列式。

// 我们实际需要的信息列表是作为FEValues构造函数的第三个参数的标志集合给出的。由于这些值必须重新计算，或者说更新，每次我们进入一个新的单元时，所有这些标志都以前缀 <code>update_</code> 开始，然后指出我们想要更新的实际内容。如果我们想要计算形状函数的值，那么给出的标志是#update_values；对于梯度，它是#update_gradients。雅各布的行列式和正交权重总是一起使用的，所以只计算乘积（雅各布乘以权重，或者简称 <code>JxW</code> ）；由于我们需要它们，我们必须同时列出#update_JxW_values。

  FEValues<2> fe_values(fe, 
                        quadrature_formula, 
                        update_values | update_gradients | update_JxW_values); 

// 这种方法的优点是，我们可以指定每个单元上究竟需要什么样的信息。很容易理解的是，这种方法可以大大加快有限元计算的速度，相比之下，所有的东西，包括二阶导数、单元的法向量等都在每个单元上计算，不管是否需要它们。

//  @note  <code>update_values | update_gradients | update_JxW_values</code>的语法对于那些不习惯用C语言编程多年的位操作的人来说不是很明显。首先， <code>operator|</code> 是<i>bitwise or operator</i>，也就是说，它接受两个整数参数，这些参数被解释为比特模式，并返回一个整数，其中每个比特都被设置，因为在两个参数中至少有一个的对应比特被设置。例如，考虑操作 <code>9|10</code>. In binary, <code>9=0b1001</code> （其中前缀 <code>0b</code> 表示该数字将被解释为二进制数字）和 <code>10=0b1010</code>  。通过每个比特，看它是否在其中一个参数中被设置，我们得出 <code>0b1001|0b1010=0b1011</code> ，或者用十进制符号表示， <code>9|10=11</code>  。你需要知道的第二个信息是，各种 <code>update_*</code> 标志都是有<i>exactly one bit set</i>的整数。例如，假设  <code>update_values=0b00001=1</code>  ,  <code>update_gradients=0b00010=2</code>  ,  <code>update_JxW_values=0b10000=16</code>  。那么<code>update_values | update_gradients | update_JxW_values = 0b10011 = 19</code>。换句话说，我们得到一个数字，即<i>encodes a binary mask representing all of the operations you want to happen</i>，其中每个操作正好对应于整数中的一个位，如果等于1，意味着每个单元格上应该更新一个特定的片断，如果是0，意味着我们不需要计算它。换句话说，即使 <code>operator|</code> 是<i>bitwise OR operation</i>，它真正代表的是<i>I want this AND that AND the other</i>。这样的二进制掩码在C语言编程中很常见，但在C++这样的高级语言中也许不是这样，但对当前的目的有很好的作用。

// 为了在下文中进一步使用，我们为一个将被频繁使用的值定义了一个快捷方式。也就是每个单元的自由度数的缩写（因为我们是在二维，自由度只与顶点相关，所以这个数字是4，但是我们更希望在写这个变量的定义时，不妨碍我们以后选择不同的有限元，每个单元有不同的自由度数，或者在不同的空间维度工作）。

// 一般来说，使用符号名称而不是硬编码这些数字是个好主意，即使你知道它们，因为例如，你可能想在某个时候改变有限元。改变元素就必须在不同的函数中进行，而且很容易忘记在程序的另一部分做相应的改变。最好不要依赖自己的计算，而是向正确的对象索取信息。在这里，我们要求有限元告诉我们每个单元的自由度数，无论我们在程序中的其他地方选择什么样的空间尺寸或多项式程度，我们都会得到正确的数字。

// 这里定义的快捷方式主要是为了讨论基本概念，而不是因为它节省了大量的输入，然后会使下面的循环更容易阅读。在大型程序中，你会在很多地方看到这样的快捷方式，`dofs_per_cell`就是一个或多或少是这类对象的传统名称。

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 

// 现在，我们说我们想逐个单元地组装全局矩阵和向量。我们可以将结果直接写入全局矩阵，但是这样做的效率并不高，因为对稀疏矩阵元素的访问是很慢的。相反，我们首先在一个小矩阵中计算每个单元的贡献，并在这个单元的计算结束后将其转移到全局矩阵中。我们对右手边的向量也是这样做的。所以我们首先分配这些对象（这些是局部对象，所有的自由度都与所有其他的自由度耦合，我们应该使用一个完整的矩阵对象，而不是一个用于局部操作的稀疏矩阵；以后所有的东西都将转移到全局的稀疏矩阵中）。

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
  Vector<double>     cell_rhs(dofs_per_cell); 

// 在集合每个单元的贡献时，我们用自由度的局部编号（即从零到dofs_per_cell-1的编号）来做。然而，当我们将结果转移到全局矩阵时，我们必须知道自由度的全局编号。当我们查询它们时，我们需要为这些数字建立一个从头开始的（临时）数组（关于这里使用的类型， types::global_dof_index, ，见介绍末尾的讨论）。

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

// 现在是所有单元格的循环。我们之前已经看到这对一个三角形是如何工作的。DoFHandler的单元格迭代器与Triangulation的迭代器完全类似，但有关于你所使用的有限元的自由度的额外信息。在自由度处理程序的活动单元上进行循环操作的方法与三角法相同。

// 注意，这次我们将单元的类型声明为`const auto &`，而不是`auto`。在第1步中，我们通过用细化指标标记来修改三角形的单元。在这里，我们只检查单元格而不修改它们，所以把`cell`声明为`const`是很好的做法，以便执行这个不变性。

  for (const auto &cell : dof_handler.active_cell_iterators()) 
    { 

// 我们现在坐在一个单元上，我们希望计算形状函数的值和梯度，以及参考单元和真实单元之间映射的雅各布矩阵的行列式，在正交点上。由于所有这些值都取决于单元格的几何形状，我们必须让FEValues对象在每个单元格上重新计算它们。

      fe_values.reinit(cell); 

// 接下来，在我们填充之前，将本地单元对全局矩阵和全局右手边的贡献重置为零。

      cell_matrix = 0; 
      cell_rhs    = 0; 

// 现在是时候开始对单元进行积分了，我们通过对所有的正交点进行循环来完成，我们将用q_index来编号。

      for (const unsigned int q_index : fe_values.quadrature_point_indices()) 
        { 

// 首先组装矩阵。对于拉普拉斯问题，每个单元格上的矩阵是形状函数i和j的梯度的积分。由于我们不进行积分，而是使用正交，所以这是在所有正交点的积分之和乘以正交点的雅各布矩阵的行列式乘以这个正交点的权重。你可以通过使用 <code>fe_values.shape_grad(i,q_index)</code> 得到形状函数 $i$ 在数字q_index的正交点上的梯度；这个梯度是一个二维向量（事实上它是张量 @<1,dim@>, 类型，这里dim=2），两个这样的向量的乘积是标量乘积，即两个shape_grad函数调用的积是点乘。这又要乘以雅各布行列式和正交点权重（通过调用 FEValues::JxW() 得到）。最后，对所有形状函数 $i$ 和 $j$ 重复上述操作。

          for (const unsigned int i : fe_values.dof_indices()) 
            for (const unsigned int j : fe_values.dof_indices()) 
              cell_matrix(i, j) += 
                (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q) 
                 fe_values.shape_grad(j, q_index) * // grad phi_j(x_q) 
                 fe_values.JxW(q_index));           // dx 

// 然后我们对右手边做同样的事情。在这里，积分是对形状函数i乘以右手边的函数，我们选择的是常值为1的函数（更有趣的例子将在下面的程序中考虑）。

          for (const unsigned int i : fe_values.dof_indices()) 
            cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q) 
                            1. *                                // f(x_q) 
                            fe_values.JxW(q_index));            // dx 
        } 

// 现在我们有了这个单元的贡献，我们必须把它转移到全局矩阵和右手边。为此，我们首先要找出这个单元上的自由度有哪些全局数字。让我们简单地询问该单元的信息。

      cell->get_dof_indices(local_dof_indices); 

// 然后再次循环所有形状函数i和j，并将局部元素转移到全局矩阵中。全局数字可以用local_dof_indices[i]获得。

      for (const unsigned int i : fe_values.dof_indices()) 
        for (const unsigned int j : fe_values.dof_indices()) 
          system_matrix.add(local_dof_indices[i], 
                            local_dof_indices[j], 
                            cell_matrix(i, j)); 

// 再来，我们对右边的向量做同样的事情。

      for (const unsigned int i : fe_values.dof_indices()) 
        system_rhs(local_dof_indices[i]) += cell_rhs(i); 
    } 

// 现在，几乎所有的东西都为离散系统的求解做好了准备。然而，我们还没有照顾到边界值（事实上，没有迪里切特边界值的拉普拉斯方程甚至不是唯一可解的，因为你可以在离散解中加入一个任意的常数）。因此，我们必须对这种情况做一些处理。

// 为此，我们首先获得边界上的自由度列表以及形状函数在那里的值。为了简单起见，我们只对边界值函数进行插值，而不是将其投影到边界上。库中有一个函数正是这样做的。  VectorTools::interpolate_boundary_values(). 它的参数是（省略存在默认值而我们不关心的参数）：DoFHandler对象，用于获取边界上自由度的全局数字；边界上边界值应被内插的部分；边界值函数本身；以及输出对象。

// 边界分量的含义如下：在很多情况下，你可能只想在边界的一部分施加某些边界值。例如，在流体力学中，你可能有流入和流出的边界，或者在身体变形计算中，身体的夹紧和自由部分。那么你就想用指标来表示边界的这些不同部分，并告诉interpolate_boundary_values函数只计算边界的某一部分（例如夹住的部分，或流入的边界）的边界值。默认情况下，所有的边界都有一个0的边界指标，除非另有规定。如果边界的部分有不同的边界条件，你必须用不同的边界指示器为这些部分编号。然后，下面的函数调用将只确定那些边界指标实际上是作为第二个参数指定的0的边界部分的边界值。

// 描述边界值的函数是一个Function类型的对象或一个派生类的对象。其中一个派生类是 Functions::ZeroFunction, ，它描述了一个到处都是零的函数（并不意外）。我们就地创建这样一个对象，并将其传递给 VectorTools::interpolate_boundary_values() 函数。

// 最后，输出对象是一对全局自由度数（即边界上的自由度数）和它们的边界值（这里所有条目都是零）的列表。这种自由度数到边界值的映射是由 <code>std::map</code> 类完成的。

  std::map<types::global_dof_index, double> boundary_values; 
  VectorTools::interpolate_boundary_values(dof_handler, 
                                           0, 
                                           Functions::ZeroFunction<2>(), 
                                           boundary_values); 

// 现在我们得到了边界DoF的列表和它们各自的边界值，让我们用它们来相应地修改方程组。这可以通过以下函数调用来实现。

  MatrixTools::apply_boundary_values(boundary_values, 
                                     system_matrix, 
                                     solution, 
                                     system_rhs); 
} 
// @sect4{Step3::solve}  

// 下面的函数简单地求解了离散化的方程。由于该系统对于高斯消除或LU分解等直接求解器来说是一个相当大的系统，我们使用共轭梯度算法。你应该记住，这里的变量数量（只有1089个）对于有限元计算来说是一个非常小的数字，而100.000是一个比较常见的数字。 对于这个数量的变量，直接方法已经不能使用了，你不得不使用CG这样的方法。

void Step3::solve() 
{ 

// 首先，我们需要有一个对象，知道如何告诉CG算法何时停止。这是通过使用SolverControl对象来实现的，作为停止标准，我们说：在最多1000次迭代后停止（这远远超过了1089个变量的需要；见结果部分以了解真正使用了多少次），如果残差的规范低于 $10^{-12}$ 就停止。在实践中，后一个标准将是停止迭代的一个标准。

  SolverControl solver_control(1000, 1e-12); 

// 然后，我们需要解算器本身。SolverCG类的模板参数是向量的类型，留下空的角括号将表明我们采取的是默认参数（即 <code>Vector@<double@></code>  ）。然而，我们明确地提到了模板参数。

  SolverCG<Vector<double>> solver(solver_control); 

// 现在求解方程组。CG求解器的第四个参数是一个预处理程序。我们觉得还没有准备好深入研究这个问题，所以我们告诉它使用身份运算作为预处理。

  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity()); 

// 现在求解器已经完成了它的工作，求解变量包含了求解函数的结点值。

} 
// @sect4{Step3::output_results}  

// 典型的有限元程序的最后一部分是输出结果，也许会做一些后处理（例如计算边界处的最大应力值，或者计算整个流出物的平均通量，等等）。我们这里没有这样的后处理，但是我们想把解决方案写到一个文件里。

void Step3::output_results() const 
{ 

// 为了将输出写入文件，我们需要一个知道输出格式等的对象。这就是DataOut类，我们需要一个该类型的对象。

  DataOut<2> data_out; 

// 现在我们必须告诉它从哪里获取它要写的值。我们告诉它使用哪个DoFHandler对象，以及求解向量（以及求解变量在输出文件中的名称）。如果我们有不止一个我们想在输出中查看的向量（例如右手边，每个单元格的错误，等等），我们也要把它们加进去。

  data_out.attach_dof_handler(dof_handler); 
  data_out.add_data_vector(solution, "solution"); 

// 在DataOut对象知道它要处理哪些数据后，我们必须告诉它把它们处理成后端可以处理的数据。原因是我们将前端（知道如何处理DoFHandler对象和数据向量）与后端（知道许多不同的输出格式）分开，使用一种中间数据格式将数据从前端传输到后端。数据通过以下函数转换为这种中间格式。

  data_out.build_patches(); 

// 现在我们已经为实际输出做好了一切准备。只要打开一个文件，用VTK格式把数据写进去就可以了（在我们这里使用的DataOut类中还有很多其他函数，可以把数据写成postscript、AVS、GMV、Gnuplot或其他一些文件格式）。

  std::ofstream output("solution.vtk"); 
  data_out.write_vtk(output); 
} 
// @sect4{Step3::run}  

// 最后，这个类的最后一个函数是主函数，调用 <code>Step3</code> 类的所有其他函数。这样做的顺序类似于大多数有限元程序的工作顺序。由于这些名字大多是不言自明的，所以没有什么可评论的。

void Step3::run() 
{ 
  make_grid(); 
  setup_system(); 
  assemble_system(); 
  solve(); 
  output_results(); 
} 
// @sect3{The <code>main</code> function}  

// 这是程序的主函数。由于主函数的概念大多是C++编程之前的面向对象时代的遗留物，所以它通常不做更多的事情，只是创建一个顶层类的对象并调用其原理函数。

// 最后，函数的第一行是用来启用deal.II可以生成的一些诊断程序的输出。  @p deallog 变量（代表deal-log，而不是de-allog）代表一个流，库的某些部分将输出写入其中。例如，迭代求解器将产生诊断程序（起始残差、求解器步骤数、最终残差），在运行这个教程程序时可以看到。

//  @p deallog 的输出可以写到控制台，也可以写到文件，或者两者都写。两者在默认情况下都是禁用的，因为多年来我们已经知道，一个程序只应该在用户明确要求的时候才产生输出。但这是可以改变的，为了解释如何做到这一点，我们需要解释 @p deallog 是如何工作的。当库的个别部分想要记录输出时，它们会打开一个 "上下文 "或 "部分"，这个输出将被放入其中。在想要写输出的部分结束时，人们再次退出这个部分。由于一个函数可以在这个输出部分打开的范围内调用另一个函数，所以输出实际上可以分层嵌套到这些部分。LogStream类（ @p deallog 是一个变量）将这些部分中的每一个称为 "前缀"，因为所有的输出都以这个前缀打印在行的左端，前缀由冒号分隔。总是有一个默认的前缀叫做 "DEAL"（暗示了deal.II的历史，它是以前一个叫做 "DEAL "的库的继承者，LogStream类是被带入deal.II的少数代码之一）。

// 默认情况下， @p logstream 只输出前缀为零的行--也就是说，所有的输出都是禁用的，因为默认的 "DEAL "前缀总是存在的。但人们可以为应该输出的行设置不同的最大前缀数，以达到更大的效果，事实上在这里我们通过调用 LogStream::depth_console(). 将其设置为两个。这意味着对于所有的屏幕输出，在默认的 "DEAL "之外再推一个前缀的上下文被允许将其输出打印到屏幕上（"控制台"），而所有进一步嵌套的部分将有三个或更多的前缀被激活，会写到 @p deallog, ，但 @p deallog 并不转发这个输出到屏幕。因此，运行这个例子（或者看 "结果 "部分），你会看到解算器的统计数据前缀为 "DEAL:CG"，这是两个前缀。这对于当前程序的上下文来说已经足够了，但是你将在以后看到一些例子（例如，在 step-22 中），其中求解器嵌套得更深，你可能通过设置更高的深度来获得有用的信息。

int main() 
{ 
  deallog.depth_console(2); 

  Step3 laplace_problem; 
  laplace_problem.run(); 

  return 0; 
} 


