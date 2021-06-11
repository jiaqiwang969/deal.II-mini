

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



// 前面几个包括的内容和前面的程序一样，所以不需要额外的注释。

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

// 然而，下一个文件是新的。我们需要这个包含文件来将自由度（DoF）与顶点、直线和单元联系起来。

#include <deal.II/dofs/dof_handler.h>

// 以下文件包含了对双线性有限元的描述，包括它在三角形的每个顶点上有一个自由度，但在面和单元内部没有自由度。

// (事实上，该文件包含了对拉格朗日元素的一般描述，即还有二次、三次等版本，而且不仅是2d，还有1d和3d。)

#include <deal.II/fe/fe_q.h>

// 在下面的文件中，可以找到几个操作自由度的工具。

#include <deal.II/dofs/dof_tools.h>

// 我们将使用一个稀疏矩阵来可视化自由度在网格上的分布所产生的非零条目模式。这个类可以在这里找到。

#include <deal.II/lac/sparse_matrix.h>

// 我们还需要使用一个中间的稀疏模式结构，可以在这个文件中找到。

#include <deal.II/lac/dynamic_sparsity_pattern.h>

// 我们希望使用一种特殊的算法来重新计算自由度。它被声明在这里。

#include <deal.II/dofs/dof_renumbering.h>

// 而这又是C++输出所需要的。

#include <fstream>

// 最后，和 step-1 一样，我们将deal.II命名空间导入到全局范围。

using namespace dealii;
// @sect3{Mesh generation}

// 这就是前面 step-1
// 例子程序中产生圆形网格的函数，细化步骤较少。唯一不同的是，它通过其参数返回它所产生的网格。

void
make_grid(Triangulation<2> &triangulation)
{
  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 5);

  for (unsigned int step = 0; step < 3; ++step)
    {
      for (auto &cell : triangulation.active_cell_iterators())
        for (const auto v : cell->vertex_indices())
          {
            const double distance_from_center =
              center.distance(cell->vertex(v));

            if (std::fabs(distance_from_center - inner_radius) <=
                1e-6 * inner_radius)
              {
                cell->set_refine_flag();
                break;
              }
          }

      triangulation.execute_coarsening_and_refinement();
    }
}
// @sect3{Creation of a DoFHandler}

// 到目前为止，我们只有一个网格，即一些几何信息（顶点的位置）和一些拓扑信息（顶点如何与线相连，线与单元格相连，以及哪些单元格与哪些其他单元格相邻）。要使用数值算法，还需要一些逻辑信息：我们希望将自由度数字与每个顶点（或线，或单元，如果我们使用高阶元素的话）联系起来，以便以后生成描述三角形上有限元场的矩阵和矢量。

// 这个函数显示了如何做到这一点。要考虑的对象是 <code>DoFHandler</code> 类模板。
// 然而，在这之前，我们首先需要一些东西来描述这些对象中的每一个要与多少个自由度相关联。由于这是有限元空间定义的一个方面，有限元基类存储了这个信息。在目前的情况下，我们因此创建了一个描述拉格朗日元素的派生类
// <code>FE_Q</code>
// 的对象。它的构造函数需要一个参数，说明元素的多项式程度，这里是1（表示一个双线性元素）；这就对应于每个顶点的一个自由度，而线和四边形内部没有自由度。如果给构造函数的值是3，我们就会得到一个双立方体元素，每个顶点有一个自由度，每条线有两个自由度，单元内有四个自由度。一般来说，
// <code>FE_Q</code>
// 表示具有完整多项式（即张量积多项式）的连续元素家族，直到指定的顺序。

// 我们首先需要创建一个这个类的对象，然后把它传递给 <code>DoFHandler</code>
// 对象，为自由度分配存储空间（用deal.II的行话说：我们<i>distribute degrees of
// freedom</i>）。

void
distribute_dofs(DoFHandler<2> &dof_handler)
{
  const FE_Q<2> finite_element(1);
  dof_handler.distribute_dofs(finite_element);

  // 现在我们已经将自由度与每个顶点的全局数字联系起来，我们想知道如何将其可视化？
  // 没有简单的方法可以直接将与每个顶点相关的自由度数字可视化。然而，这样的信息几乎不会真正重要，因为编号本身或多或少是任意的。还有更重要的因素，我们将在下文中展示其中一个。

  // 与三角形的每个顶点相关的是一个形状函数。假设我们想解决类似拉普拉斯方程的问题，那么不同的矩阵条目将是每对这样的形状函数的梯度的积分。显然，由于形状函数只在与它们相关的顶点相邻的单元格上是非零的，所以只有当与该列和行%号相关的形状函数的支持相交时，矩阵条目才是非零的。这只是相邻形状函数的情况，因此也只是相邻顶点的情况。现在，由于顶点被上述函数
  // (DoFHandler::distribute_dofs),
  // 或多或少地随机编号，矩阵中非零项的模式将有些参差不齐，我们现在就来看看它。

  // 首先，我们要创建一个结构，用来存储非零元素的位置。然后，这个结构可以被一个或多个稀疏矩阵对象使用，这些对象在这个稀疏模式所存储的位置上存储条目的值。存储这些位置的类是SparsityPattern类。然而，事实证明，当我们试图立即填充这个类时，它有一些缺点：它的数据结构的设置方式是，我们需要对我们可能希望在每一行的最大条目数有一个估计。在两个空间维度上，通过 DoFHandler::max_couplings_between_dofs() 函数可以得到合理的估计值，但是在三个维度上，该函数几乎总是严重高估真实的数字，导致大量的内存浪费，有时对于所使用的机器来说太多，即使未使用的内存可以在计算稀疏模式后立即释放。为了避免这种情况，我们使用了一个中间对象DynamicSparsityPattern，该对象使用了一个不同的%内部数据结构，我们可以随后将其复制到SparsityPattern对象中，而不需要太多的开销。关于这些数据结构的一些更多信息可以在 @ref Sparsity 模块中找到）。为了初始化这个中间数据结构，我们必须给它提供矩阵的大小，在我们的例子中，矩阵是正方形的，行和列的数量与网格上的自由度相同。

  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());

  // 然后我们在这个对象中填入非零元素的位置，考虑到目前自由度的编号。

  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);

  // 现在我们已经准备好创建实际的稀疏模式了，以后我们可以用在我们的矩阵上。它将包含已经在DynamicSparsityPattern中集合的数据。

  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);

  // 有了这个，我们现在可以把结果写到一个文件里。

  std::ofstream out("sparsity_pattern1.svg");
  sparsity_pattern.print_svg(out);

  // 结果被存储在一个 <code>.svg</code>
  // 文件中，矩阵中的每个非零条目都对应于图像中的一个红色方块。输出结果将显示如下。

  // 如果你看一下，你会注意到稀疏性模式是对称的。这不应该是一个惊喜，因为我们没有给
  // <code>DoFTools::make_sparsity_pattern</code>
  // 任何信息，表明我们的双线性形式可能以非对称的方式耦合形状函数。你还会注意到它有几个明显的区域，这源于编号从最粗的单元开始，然后到较细的单元；由于它们都是围绕原点对称分布的，这在稀疏模式中再次显示出来。
}
// @sect3{Renumbering of DoFs}

// 在上面产生的稀疏模式中，非零条目在对角线上延伸得很远。对于某些算法来说，例如不完全LU分解或Gauss-Seidel预处理，这是不利的，我们将展示一个简单的方法来改善这种情况。

// 请记住，为了使矩阵中的一个条目 $(i,j)$
// 不为零，形状函数i和j的支持需要相交（否则在积分中，积分将到处为零，因为在某个点上，一个或另一个形状函数为零）。然而，形状函数的支撑点只有在彼此相邻的情况下才会相交，所以为了使非零条目聚集在对角线周围（其中
// $i$ 等于 $j$ ），我们希望相邻的形状函数的索引（DoF编号）相差不大。

// 这可以通过一个简单的前行算法来实现，即从一个给定的顶点开始，给它的索引为0。然后，依次对其邻居进行编号，使其指数接近于原始指数。然后，他们的邻居，如果还没有被编号，也被编号，以此类推。

// 有一种算法沿着这些思路增加了一点复杂性，那就是Cuthill和McKee的算法。我们将在下面的函数中使用它来对自由度进行重新编号，从而使产生的稀疏模式在对角线周围更加本地化。该函数唯一有趣的部分是对
// <code>DoFRenumbering::Cuthill_McKee</code>
// 的第一次调用，其余部分基本上与以前一样。

void
renumber_dofs(DoFHandler<2> &dof_handler)
{
  DoFRenumbering::Cuthill_McKee(dof_handler);

  DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(),
                                                  dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);

  SparsityPattern sparsity_pattern;
  sparsity_pattern.copy_from(dynamic_sparsity_pattern);

  std::ofstream out("sparsity_pattern2.svg");
  sparsity_pattern.print_svg(out);
}

// 再次，输出如下。请注意，非零项在对角线附近的聚类情况要比以前好得多。这种效果对于较大的矩阵来说更加明显（目前的矩阵有1260行和列，但是大的矩阵往往有几十万行）。

// 值得注意的是， <code>DoFRenumbering</code>
// 类也提供了一些其他的算法来重新编号自由度。例如，如果所有的耦合都在矩阵的下三角或上三角部分，那当然是最理想的，因为这样的话，解决线性系统就只需要向前或向后替换。当然，这对于对称稀疏模式来说是无法实现的，但在一些涉及传输方程的特殊情况下，通过列举从流入边界沿流线到流出边界的自由度，这是可能的。毫不奇怪，
// <code>DoFRenumbering</code> 也有这方面的算法。

//  @sect3{The main function}

// 最后，这是主程序。它所做的唯一一件事就是分配和创建三角形，然后创建一个
// <code>DoFHandler</code> 对象并将其与三角形相关联，最后对其调用上述两个函数。

int
main()
{
  Triangulation<2> triangulation;
  make_grid(triangulation);

  DoFHandler<2> dof_handler(triangulation);

  distribute_dofs(dof_handler);
  renumber_dofs(dof_handler);
}
