//include/deal.II-translator/dofs/dof_renumbering_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_dof_renumbering_h
#define dealii_dof_renumbering_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/hp/dof_handler.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 实现一些对三角形自由度的重新编号算法。这个命名空间中的函数为DoFHandler对象的每个自由度计算新的索引，然后调用
 * DoFHandler::renumber_dofs().  。   <h3>Cuthill-McKee like algorithms</h3>
 * 在这个类中，Cuthill-McKee算法被实现。它从一个自由度开始，在其他自由度中搜索那些与我们开始的自由度相耦合的自由度，并以一定的方式对这些自由度进行编号。然后它找到第二层的自由度，即那些与上一层的自由度（也就是与初始自由度耦合的自由度）耦合的自由度，并对这些自由度进行编号。以此类推。关于该算法的细节，特别是每个层次内的编号，请参见H.R.Schwarz："Methode
 * der finiten
 * Elemente"。反向的Cuthill-McKee算法也做同样的工作，但对所有元素的编号顺序是相反的。
 * 这些算法有一个主要的缺点：它们需要一个好的起点，即自由度指数，将得到一个新的指数为零。因此，重新编号函数允许调用者指定这样一个初始自由度，例如通过利用域的实际拓扑结构知识。也可以给出几个起始指数，可以用来模拟简单的上游编号（通过给流入的DoF作为起始值）或使预处理更快（通过让Dirichlet边界指数作为起始点）。
 * 如果没有给出起始指数，会自动选择一个，即协调数最小的一个（协调数是这个道口与其他道口的耦合数）。这个道夫通常位于域的边界上。然而，在使用本库中使用的分层网格时，这一点存在很大的模糊性，因为在大多数情况下，计算域不是通过倾斜和变形元素以及在顶点插入可变数量的元素来逼近的，而是通过分层细化。因此，有大量的协调数相等的道夫。因此，重新编号的算法将不会得到最佳的结果。
 * 在施瓦茨（H.R.Schwarz: Methode der finiten
 * Elemente）的书中，建议测试许多起点，如果可能的话，所有的协调数都是最小的，还有那些略高的协调数。然而，这似乎只适用于20世纪80年代初的小型工程问题（第二版于1984年出版）中的最多几十个或几百个元素的网格，但肯定不适用于本库中的那些具有几万至几十万个元素的网格。
 *
 *  <h4>Implementation of renumbering schemes</h4>
 * 重新编号的算法需要相当多的内存，因为他们必须为每个道夫存储它与哪些其他道夫的耦合。这是用一个SparsityPattern对象完成的，用来存储矩阵的稀疏模式。用户在分配道夫和重新编号之间做任何事情都是没有用的，也就是说，对
 * DoFHandler::distribute_dofs 和 DoFHandler::renumber_dofs
 * 的调用应该紧接着进行。如果你试图在两者之间创建一个稀疏模式或其他任何东西，这些之后将是无效的。
 * 重新编号可能会照顾到仅由消除约束引起的dof-to-of耦合。除了上面提到的内存消耗外，这也需要相当多的计算时间，但它可以在调用
 * @p renumber_dofs
 * 函数时被关闭。这样就会得到较差的结果，因为图中的结（代表道夫）不会被发现是邻居，即使它们在浓缩后会是邻居。
 * 重新编号的算法是在纯代数的基础上工作的，这是因为算法所依据的图论基础和稀疏模式所代表的二进制矩阵（其条目为二进制值的矩阵）之间存在同构性。特别是，这些算法并不试图利用拓扑学知识（如角检测）来寻找适当的起点。然而，这样一来，他们在任意的空间维度上工作。
 * 如果你想给起点，你可以给一个dof指数的列表，这将形成重新编号的第一步。列表中的道夫将从零开始连续编号，也就是说，这个列表不是根据节点的协调数来重新编号的。不在允许范围内的索引被删除。如果不允许有索引，算法将搜索自己的起点。
 *
 *  <h4>Results of renumbering</h4>
 * 上面提到的重新编号方案并没有导致最佳结果。然而，毕竟没有任何算法能在合理的时间内完成这个任务。在有些情况下，缺乏最优性甚至会导致比原来粗略的逐级编号方案更差的结果；其中一个例子是一个由四个单元组成的网格，其中总是精炼那些与中心相邻的单元（你可以称这种网格为
 * "放大 "网格）。在这样一个例子中，带宽增加了大约50%。
 * 在其他大多数情况下，带宽会明显减少。网格的结构越少，减少的效果就越好。有一个网格的单元是按照随机驱动的算法进行细化的，带宽减少了6倍。
 * 使用约束信息通常会使带宽减少10%或20%，但对于一些非常非结构化的网格，也可能导致带宽增加。你必须权衡在你的情况下减少的带宽与使用约束信息所花费的时间，后者通常比
 * "纯 "重新编号的算法长几倍。
 * 在几乎所有情况下，重新编号方案都会找到一个角来开始。由于大多数网格中有不止一个角，而且即使是内部自由度也可能是一个更好的起点，如果你有一个简单的方案来推导出一个合适的点，那么由用户给出起点可能是一个可行的方法（例如，如果你想得到左上角的顶点，可以依次取最粗层次的单元的第三个孩子，取其第三个顶点和其道夫指数）。如果你事先不知道你的网格会是什么样子（例如，当使用自适应算法时），寻找一个最佳的起始点可能会很困难，然而，在许多情况下，这并不值得努力。
 *
 *  <h3>Component-wise and block-wise numberings</h3>
 * 对于使用FESystem类的几个基本元素组成的有限元，或者对于本身提供了几个组件的元素，按组件对DoF指数进行排序可能是有意义的。这样就能显示出块状矩阵结构，因为否则自由度是按单元编号的，而没有考虑到它们可能属于不同的构件。例如，我们可能希望对斯托克斯离散化的自由度进行排序，这样我们首先得到所有的速度，然后是所有的压力，这样得到的矩阵就会自然地分解为
 * $2\times 2$ 系统。
 * 这种编号可以通过调用本类的component_wise()函数获得。由于它没有触及每个分量中的指数顺序，因此可能值得首先使用Cuthill-
 * McKee或类似的算法进行重新编号，然后再按分量重新编号。这将带出矩阵结构，并在每个块内有一个良好的编号。
 * component_wise()函数不仅允许基于向量成分的枚举，还允许使用各种 DoFRenumbering::component_wise() 函数的默认参数将向量成分组合成 "块"
 * （见 @ref GlossComponent 与 @ref GlossBlock
 * 的区别说明）。通过这个参数指定的块可以，但不一定要等于有限元报告的块。例如，一个典型的斯托克斯元素将是
 *
 * @code
 * FESystem<dim> stokes_fe (FE_Q<dim>(2), dim,   // dim velocities
 *                          FE_Q<dim>(1), 1);    // one pressure
 * @endcode
 * 这个元素有 <code>dim+1</code>
 * 个矢量分量和同样多的块。然而，人们可能希望将速度视为一个逻辑块，这样所有的速度自由度都以同样的方式列举，与它们是
 * $x$  - 还是 $y$  - 速度无关。例如，在 step-20 和 step-22
 * 以及其他几个教程程序中就这样做了。
 * 另一方面，如果你真的想使用有限元本身报告的块状结构（如果你的有限元有多个矢量分量，例如FE_RaviartThomas或FE_Nedelec元，这种情况经常发生），那么你可以使用
 * DoFRenumbering::block_wise 而不是 DoFRenumbering::component_wise
 * 函数。
 *
 *  <h3>Cell-wise numbering</h3>
 * 给定一个有序的单元格向量，函数cell_wise()对自由度进行排序，使该向量中较早的单元格上的自由度出现在较晚的单元格上的自由度之前。
 * 对于非连续Galerkin方法（FE_DGP,
 * FE_DGQ），这一规则产生了一个明确的排序。对于连续方法，我们使用额外的规则，即每个自由度根据它所属的有序向量中的第一个单元来排序。
 * 这个方案的应用是downstream()和clock_wise_dg()。第一个是根据下游方向对单元进行排序，然后应用cell_wise()。
 *
 *
 * @note
 * 对于DG元素，每个单元的内部编号保持不受影响。这对于连续元素来说已经不能保证了，因为与前面的单元共享的自由度会被另一个单元所占。
 *
 *  <h3>Random renumbering</h3>
 * random()函数对自由度进行随机编号。这个函数可能很少使用，除了检查求解器（迭代或直接）对自由度编号的依赖性。
 *
 *  <h3>A comparison of reordering strategies</h3>
 * 作为比较的基准，让我们考虑一下，当使用通常用于斯托克斯方程离散化的
 * $Q_2^d\times Q_1$ 元素组合时，在 step-22
 * 中得到的网格上，经过三维自适应网格细化后，各种算法产生的不同稀疏性模式是什么。空间维度和耦合有限元导致了一个相当密集的系统矩阵，平均每行大约有180个非零条目。在应用下面所示的每一种重排策略后，自由度也用
 * DoFRenumbering::component_wise
 * 排序为速度组和压力组；这产生了下面的 $2\times 2$
 * 块结构，大的速度-速度块在左上方，小的压力-压力块在右下方，而耦合块在右上方和左下方。
 * 重排策略的目的是为了改进预处理程序。在 step-22
 * 中，我们使用SparseILU对左上方的速度-速度块进行预处理。预处理程序的质量可以通过解决这个块的线性系统所需的CG迭代次数来衡量。对于下面的一些重排策略，我们记录了自适应细化周期3的这个数字，有93176个自由度；因为我们用舒尔补数中的同一个矩阵来解决几个线性系统，所以报告了平均迭代次数。数字越小，预处理程序就越好，因此自由度的重新编号就越适合这项任务。我们还说明了程序的运行时间，部分由所需的迭代次数决定，在我们的一台机器上的前4个周期。请注意，报告的时间对应于整个程序的运行时间，而不仅仅是受影响的求解器；如果一个程序用一个特定的排序比用另一个排序运行快一倍，那么这意味着实际的求解器实际上要快几倍。
 * <table> <tr> <td>
 @image html "reorder_sparsity_step_31_original.png"
 * </td> <td>
 @image html "reorder_sparsity_step_31_random.png"
 * </td> <td>
 @image html "reorder_sparsity_step_31_deal_cmk.png"
 * </td> </tr> <tr> <td> Enumeration as produced by deal.II's
 * DoFHandler::distribute_dofs function and no further reordering apart from
 * the component-wise one.
 *
 * With this renumbering, we needed an average of 92.2 iterations for the
 * testcase outlined above, and a runtime of 7min53s. </td> <td> Random
 * enumeration as produced by applying DoFRenumbering::random after calling
 * DoFHandler::distribute_dofs. This enumeration produces nonzero entries in
 * matrices pretty much everywhere, appearing here as an entirely unstructured
 * matrix.
 *
 * With this renumbering, we needed an average of 71 iterations for the
 * testcase outlined above, and a runtime of 10min55s. The longer runtime
 * despite less iterations compared to the default ordering may be due to the
 * fact that computing and applying the ILU requires us to jump back and forth
 * all through memory due to the lack of localization of matrix entries around
 * the diagonal; this then leads to many cache misses and consequently bad
 * timings. </td> <td> Cuthill-McKee enumeration as produced by calling the
 * deal.II implementation of the algorithm provided by
 * DoFRenumbering::Cuthill_McKee after DoFHandler::distribute_dofs.
 *
 * With this renumbering, we needed an average of 57.3 iterations for the
 * testcase outlined above, and a runtime of 6min10s. </td> </td> </tr>
 *
 * <tr> <td>
 @image html "reorder_sparsity_step_31_boost_cmk.png"
 * </td> <td>
 @image html "reorder_sparsity_step_31_boost_king.png"
 * </td> <td>
 @image html "reorder_sparsity_step_31_boost_md.png"
 * </td> </tr> <tr> <td> Cuthill- McKee enumeration as produced by calling the
 * BOOST implementation of the algorithm provided by
 * DoFRenumbering::boost::Cuthill_McKee after DoFHandler::distribute_dofs.
 *
 * With this renumbering, we needed an average of 51.7 iterations for the
 * testcase outlined above, and a runtime of 5min52s. </td> <td> King
 * enumeration as produced by calling the BOOST implementation of the
 * algorithm provided by DoFRenumbering::boost::king_ordering after
 * DoFHandler::distribute_dofs. The sparsity pattern appears denser than with
 * BOOST's Cuthill-McKee algorithm; however, this is only an illusion: the
 * number of nonzero entries is the same, they are simply not as well
 * clustered.
 *
 * With this renumbering, we needed an average of 51.0 iterations for the
 * testcase outlined above, and a runtime of 5min03s. Although the number of
 * iterations is only slightly less than with BOOST's Cuthill-McKee
 * implementation, runtime is significantly less. This, again, may be due to
 * cache effects. As a consequence, this is the algorithm best suited to the
 * testcase, and is in fact used in step-22. </td> <td> Minimum degree
 * enumeration as produced by calling the BOOST implementation of the
 * algorithm provided by DoFRenumbering::boost::minimum_degree after
 * DoFHandler::distribute_dofs. The minimum degree algorithm does not attempt
 * to minimize the bandwidth of a matrix but to minimize the amount of fill-in
 * a LU decomposition would produce, i.e. the number of places in the matrix
 * that would be occupied by elements of an LU decomposition that are not
 * already occupied by elements of the original matrix. The resulting sparsity
 * pattern obviously has an entirely different structure than the ones
 * produced by algorithms trying to minimize the bandwidth.
 *
 * With this renumbering, we needed an average of 58.9 iterations for the
 * testcase outlined above, and a runtime of 6min11s. </td> </tr>
 *
 * <tr> <td>
 @image html "reorder_sparsity_step_31_downstream.png"
 * </td> <td> </td> <td> </td> </tr> <tr> <td> Downstream enumeration using
 * DoFRenumbering::downstream using a direction that points diagonally through
 * the domain.
 *
 * With this renumbering, we needed an average of 90.5 iterations for the
 * testcase outlined above, and a runtime of 7min05s. </td> <td> </td> <td>
 * </td> </tr> </table>   <h3>Multigrid DoF numbering</h3>
 * 上面列出的大多数算法也适用于多网格自由度的编号。请参考实际的函数声明以获得更多这方面的信息。
 *
 *
 * @ingroup dofs
 *
 *
 */
namespace DoFRenumbering
{
  /**
   * 基于方向的单元格迭代器的比较器：如果第二个单元格的中心在第一个单元格的中心的下游，相对于构造函数所给的方向，它返回
   * @p true 。
   *
   */
  template <class Iterator, int dim>
  struct CompareDownstream
  {
    /**
     * 构造函数。
     *
     */
    CompareDownstream(const Tensor<1, dim> &dir)
      : dir(dir)
    {}
    /**
     * 如果c1小于c2，返回true。
     *
     */
    bool
    operator()(const Iterator &c1, const Iterator &c2) const
    {
      const Tensor<1, dim> diff = c2->center() - c1->center();
      return (diff * dir > 0);
    }

  private:
    /**
     * 流动方向。
     *
     */
    const Tensor<1, dim> dir;
  };


  /**
   * 基于点的下游方向的比较器：如果第二个点相对于第一个点的下游方向是给构造函数的，则返回
   * @p true
   * 。如果这两个点相对于下游方向是相同的，那么具有较低DoF数的点被认为是较小的。
   *
   */
  template <int dim>
  struct ComparePointwiseDownstream
  {
    /**
     * 构造函数。
     *
     */
    ComparePointwiseDownstream(const Tensor<1, dim> &dir)
      : dir(dir)
    {}
    /**
     * 如果c1小于c2，返回true。
     *
     */
    bool
    operator()(const std::pair<Point<dim>, types::global_dof_index> &c1,
               const std::pair<Point<dim>, types::global_dof_index> &c2) const
    {
      const Tensor<1, dim> diff = c2.first - c1.first;
      return (diff * dir > 0 || (diff * dir == 0 && c1.second < c2.second));
    }

  private:
    /**
     * 流动方向。
     *
     */
    const Tensor<1, dim> dir;
  };



  /**
   * 一个命名空间，用于实现一些基于Jeremy Siek等人在Boost
   * Graph Library（BGL）中实现的重新编号算法。
   * 虽然计算速度通常稍慢，但使用BOOST的算法往往会导致矩阵的带宽更小，因此基于这种编号的稀疏ILU更有效率。
   * 关于这些算法与DoFRenumbering中定义的算法的比较，请参见DoFRenumbering命名空间的文档中的比较部分。
   *
   */
  namespace boost
  {
    /**
     * 根据Cuthill-McKee方法对自由度重新编号，最终使用反向编号方案。
     * 关于不同方法的细节，请参见父类的一般文档。
     * 作为这种算法结果的一个例子，看看DoFRenumbering命名空间的文档中各种算法的比较。
     *
     */
    template <int dim, int spacedim>
    void
    Cuthill_McKee(DoFHandler<dim, spacedim> &dof_handler,
                  const bool                 reversed_numbering = false,
                  const bool                 use_constraints    = false);

    /**
     * 计算Cuthill_McKee()函数所需的重新编号向量。
     * 不对DoFHandler
     * dofs进行重新编号，但返回重新编号的向量。
     *
     */
    template <int dim, int spacedim>
    void
    compute_Cuthill_McKee(std::vector<types::global_dof_index> &new_dof_indices,
                          const DoFHandler<dim, spacedim> &,
                          const bool reversed_numbering = false,
                          const bool use_constraints    = false);

    /**
     * 根据King算法的BOOST实现对自由度进行重新编号。这通常会导致比Cuthill-McKee算法稍大（百分之几）的带宽，但是稀疏ILU通常是稍好（也是百分之几）的前置条件。
     * 作为这种算法结果的一个例子，看看DoFRenumbering命名空间的文档中各种算法的比较。
     * 这种算法在  step-22  中使用。
     *
     */
    template <int dim, int spacedim>
    void
    king_ordering(DoFHandler<dim, spacedim> &dof_handler,
                  const bool                 reversed_numbering = false,
                  const bool                 use_constraints    = false);

    /**
     * 计算King算法的重新编号，但实际上不对DoF处理参数中的自由度进行重新编号。
     *
     */
    template <int dim, int spacedim>
    void
    compute_king_ordering(std::vector<types::global_dof_index> &new_dof_indices,
                          const DoFHandler<dim, spacedim> &,
                          const bool reversed_numbering = false,
                          const bool use_constraints    = false);

    /**
     * 根据最小度算法的BOOST实现对自由度进行重新编号。与Cuthill-McKee算法不同，该算法并不试图最小化矩阵的带宽，而是试图在进行LU分解时最小化填充量。由于这一特性，它有时可能会产生更好的ILU。
     * 作为这种算法结果的一个例子，看看DoFRenumbering命名空间的文档中各种算法的比较。
     *
     */
    template <int dim, int spacedim>
    void
    minimum_degree(DoFHandler<dim, spacedim> &dof_handler,
                   const bool                 reversed_numbering = false,
                   const bool                 use_constraints    = false);

    /**
     * 计算最小度数算法的重新编号，但实际上不对DoF处理参数中的自由度进行重新编号。
     *
     */
    template <int dim, int spacedim>
    void
    compute_minimum_degree(
      std::vector<types::global_dof_index> &new_dof_indices,
      const DoFHandler<dim, spacedim> &,
      const bool reversed_numbering = false,
      const bool use_constraints    = false);
  } // namespace boost

  /**
   * 根据Cuthill-McKee方法对自由度重新编号，可能使用反向编号方案。    关于不同方法的细节，请参见这个类的一般文档。    作为这种算法结果的一个例子，看看DoFRenumbering命名空间的文档中各种算法的比较。      @param  dof_handler 要工作的DoFHandler对象。    @param  reversed_numbering 是否使用原始的Cuthill-McKee算法，或者颠倒排序。    @param  use_constraints 在确定自由度的重新排序时是否使用悬挂节点约束。    @param  starting_indices 一组自由度，构成重新编号的第一层自由度。如果这个集合是空的，那么就会在那些具有最小数量的与之耦合的其他自由度中自动选择一个起始条目。    <h4> Operation in parallel </h4> 如果给定的DoFHandler使用分布式三角形（即，如果dof_handler.local_owned()不是完整的索引集），重新编号将在每个处理器的自由度上单独执行，处理器之间没有任何通信。换句话说，由此产生的重新编号是分别对<i>each diagonal block of the matrix corresponding to one processor</i>的带宽进行最小化的尝试，而没有对全局矩阵的带宽进行最小化的尝试。此外，重新编号完全重复使用每个处理器之前使用的同一组DoF指数。换句话说，如果一个处理器上先前的DoF编号使用了DoF指数的连续范围，那么在重新编号后，该处理器上的DoF也将如此，它们将占据同样的范围。如果一个处理器上的DoF以前的编号是由一些索引范围或单个索引组成的，情况也是如此：重新编号后，该处理器上的本地拥有的DoF将使用完全相同的索引，只是顺序不同。    此外，如果DoFHandler是建立在平行三角形上的，那么在每个处理器上，重新编号的起始指数需要是 @ref GlossLocallyActiveDof "本地活动自由度 "
   * 的一个（可能是空）子集。一般来说，这些起始指数在每个处理器上都是不同的（当然，除非你传递一个默认的空列表），每个处理器将使用它们作为该处理器上本地重新编号的起始指数。
   * 起始指数必须是本地活动的自由度，但该函数将只对本地拥有的自由度子集进行重新编号。该函数接受来自最大的本地活动自由度集合的起始索引，因为用这个函数进行的典型的重新编号操作是从位于边界上的索引开始的
   *
   * - 在当前函数的情况下，这将是处理器子域之间的边界。由于位于子域界面上的自由度可能被拥有相邻子域的两个处理器中的任何一个所拥有，所以要确定本地拥有的起始索引并不总是容易。另一方面，子域界面上的所有自由度都是本地活动的，因此该函数接受它们作为起始索引，尽管它只能在给定的处理器上对它们进行重新编号，如果它们也是本地拥有的。
   *
   */
  template <int dim, int spacedim>
  void
  Cuthill_McKee(DoFHandler<dim, spacedim> &dof_handler,
                const bool                 reversed_numbering = false,
                const bool                 use_constraints    = false,
                const std::vector<types::global_dof_index> &starting_indices =
                  std::vector<types::global_dof_index>());

  /**
   * 计算Cuthill_McKee()函数所需的重新编号向量。
   * 这个函数不对DoFHandler
   * DoF进行重新编号，只返回重新编号的向量。
   * 如果一个有效的级别作为参数被传递，那么这个网格级别的重新编号向量将被返回。
   * 其他参数的解释见Cuthill_McKee()函数。
   *
   */
  template <int dim, int spacedim>
  void
  compute_Cuthill_McKee(
    std::vector<types::global_dof_index> &new_dof_indices,
    const DoFHandler<dim, spacedim> &,
    const bool                                  reversed_numbering = false,
    const bool                                  use_constraints    = false,
    const std::vector<types::global_dof_index> &starting_indices =
      std::vector<types::global_dof_index>(),
    const unsigned int level = numbers::invalid_unsigned_int);

  /**
   * 根据Cuthill-McKee方法对自由度进行重新编号，最终使用反向编号方案，在这种情况下为自由度的多网格编号。
   * 你可以给出一个三角形水平，这个函数将被应用于此。
   * 因为在逐级编号的情况下，没有悬空的节点，所以不能使用约束，所以前一个函数的相应参数被省略了。
   * 有关不同方法的细节，请参见该类的一般文档。
   *
   */
  template <int dim, int spacedim>
  void
  Cuthill_McKee(DoFHandler<dim, spacedim> &dof_handler,
                const unsigned int         level,
                const bool                 reversed_numbering = false,
                const std::vector<types::global_dof_index> &starting_indices =
                  std::vector<types::global_dof_index>());

  /**
   * @name 分量上的编号  
     * @{ 
   *
   */

  /**
   * 按向量分量对自由度进行排序。每个分量内的编号不被触及，所以一个自由度的索引
   * $i$  ，属于某个分量，而另一个自由度的索引  $j$
   * 属于同一个分量，将被分配新的索引  $n(i)$  和  $n(j)$
   * ，如果  $i<j$  ，则  $n(i)>n(j)$  。
   * 你可以指定组件的排序方式与你使用的FESystem对象所建议的方式不同。为此，设置向量
   * @p target_component ，使索引 @p i
   * 处的条目表示FESystem中具有组件 @p i
   * 的目标组件的编号。
   * 对同一目标组件进行多次命名是可能的，其结果是将几个组件阻塞成一个。这将在
   * step-22
   * 中讨论。如果你省略了这个参数，就会使用有限元给出的相同顺序。
   * 如果这里考虑的全局有限元所来自的基础有限元之一是一个非原始的有限元，即它的形状函数有一个以上的非零分量，那么就不可能将这些自由度与单个矢量分量联系起来。在这种情况下，它们被与它们所属的第一个矢量分量相关联。
   * 对于只有一个分量的有限元，或单一的非原始基元，这个函数是身份操作。
   *
   */
  template <int dim, int spacedim>
  void
  component_wise(DoFHandler<dim, spacedim> &      dof_handler,
                 const std::vector<unsigned int> &target_component =
                   std::vector<unsigned int>());


  /**
   * 按分量对自由度进行排序。它的作用与上述函数相同，只是它对多级离散化中的单级进行了操作。DoFHandler的非多网格部分没有被触及。
   *
   */
  template <int dim, int spacedim>
  void
  component_wise(DoFHandler<dim, spacedim> &      dof_handler,
                 const unsigned int               level,
                 const std::vector<unsigned int> &target_component =
                   std::vector<unsigned int>());

  /**
   * 计算component_wise()函数所需的重新编号矢量。
   * 不对DoFHandler
   * dofs进行重新编号，但返回重新编号的向量。
   *
   */
  template <int dim, int spacedim, typename CellIterator>
  types::global_dof_index
  compute_component_wise(std::vector<types::global_dof_index> &new_dof_indices,
                         const CellIterator &                  start,
                         const typename identity<CellIterator>::type &end,
                         const std::vector<unsigned int> &target_component,
                         const bool                       is_level_operation);

  /**
   * @}
   *
   */

  /**
   * @name  块状的编号  @{
   *
   */

  /**
   * 按矢量块对自由度进行排序。每个块内的编号不被触及，所以一个自由度的索引
   * $i$ ，属于某个块，而另一个自由度的索引 $j$
   * 属于同一个块，将被分配新的索引 $n(i)$ 和 $n(j)$ ，如果
   * $i<j$ ，则 $n(i)>n(j)$ ，如果 $i>j$  。
   * @note  只有当附加到DoFHandler参数的 hp::FECollection 中的每个元素都有完全相同的块数时，该函数才会成功（更多信息见 @ref GlossBlock "术语表"
   * ）。请注意，这并不总是给定的：虽然 hp::FECollection
   * 类确保其所有元素具有相同数量的向量分量，但它们不需要有相同数量的块。同时，这里的这个函数需要跨元素匹配单个块，因此要求元素具有相同的块数，并且一个元素中的后续块与另一个元素中的意义相同。
   *
   */
  template <int dim, int spacedim>
  void
  block_wise(DoFHandler<dim, spacedim> &dof_handler);

  /**
   * 按矢量块对自由度进行排序。它的作用与上述函数相同，只是它对多级离散化中的一个单级进行排序。DoFHandler的非多网格部分不被触及。
   *
   */
  template <int dim, int spacedim>
  void
  block_wise(DoFHandler<dim, spacedim> &dof_handler, const unsigned int level);

  /**
   * 计算block_wise()函数所需的重新编号矢量。  不在DoFHandler
   * dofs上执行重新编号，但返回重新编号向量。
   *
   */
  template <int dim, int spacedim, class ITERATOR, class ENDITERATOR>
  types::global_dof_index
  compute_block_wise(std::vector<types::global_dof_index> &new_dof_indices,
                     const ITERATOR &                      start,
                     const ENDITERATOR &                   end,
                     bool                                  is_level_operation);

  /**
   * @}
   *
   */

  /**
   * @name  各种单元格的编号 
     * @{ 
   *
   */

  /**
   * 通过以 @ref GlossZOrder "Z顺序 "
   * 遍历单元格，逐个单元格重新编号。
   * 使用这个函数有两个原因。
   *
   *
   *
   *
   *
   *
   * - 它产生一个可预测的自由度排序，与你究竟是如何得到一个网格的无关。    特别是，一般来说，一个网格的单元的顺序取决于在网格所经历的细化周期中，单元被标记为细化和粗化的顺序。另一方面，单元格的Z轴顺序与网格的历史无关，因此产生一个可预测的DoF编号。
   *
   *
   *
   *
   *
   *
   * - 对于基于 parallel::distributed::Triangulation, 的网格，每个MPI进程的 @ref GlossLocallyOwnedCell "本地拥有的单元 "按Z顺序是连续的。这意味着通过按Z顺序访问单元来对自由度进行编号可以得到 @ref GlossLocallyOwnedDof "本地拥有的DoF指数"，它由每个进程的连续范围组成。这对于这种三角形上的DoF的默认排序也是如此，但是默认排序产生的枚举也取决于有多少个处理器参与到网格中，而这个函数产生的枚举是在一个特定的单元上的自由度的指数，无论网格被分割成多少个进程都是一样的。    对于基于 parallel::shared::Triangulation, 的网格，情况更为复杂。在这里，本地拥有的单元的集合是由分区算法决定的（通过传递一个 parallel::shared::Triangulation::Settings 类型的对象给三角化的构造器来选择），一般来说，这些分区算法可能会根据与Z顺序无关的决定将单元分配到 @ref GlossSubdomainId 的 "子域"。(尽管有可能以一种方式选择这些标志，从而使分区使用Z顺序)。因此，一个子域的单元在Z顺序中是不连续的，如果根据单元的Z顺序对自由度进行重新编号，通常会导致每个处理器上的自由度指数不形成一个连续的范围。  这通常是不方便的（例如，因为PETSc不能存储本地拥有的指数集不连续的向量和矩阵），因此这个函数对 parallel::shared::Triangulation 对象使用以下算法。
   *
   *
   *
   *
   *
   *
   * - 它决定了每个处理器拥有多少个自由度。    这是在重新编号下的一个不变因素，因此我们可以使用每个处理器在当前函数开始时拥有多少个自由度。让我们把这个数字称为 $n_P$ ，用于处理器 $P$  。
   *
   *
   *
   *
   *
   *
   * - 它为每个处理器确定一个新的DoF指数 $[b_P,e_P)$ 的连续范围，以便 $e_P-b_P=n_P$ 、 $b_0=0$ 和 $b_P=e_{P-1}$  。
   *
   *
   *
   *
   *
   *
   * - 它按Z顺序遍历<i>locally owned cells</i>，并对这些单元上本地拥有的自由度重新编号，使新的数字符合区间 $[b_P,e_P)$ 。  换句话说，每个处理器上的<i>locally owned degrees of freedom</i>是根据它们所在的本地拥有的单元的Z顺序进行排序的，但这一属性可能在全局上、在各单元之间不成立。这是因为分区算法可能已经决定，例如，处理器0拥有的单元格比分配给处理器1的一个单元格的Z顺序来<i>later</i>。  另一方面，上述算法在这个单元上分配的自由度<i>earlier</i>指数比处理器1拥有的所有指数都要高。
   * @note
   * 这个函数产生的排序与之前自由度的编号无关。换句话说，任何可能由先前调用重新编号函数产生的信息都被忽略。
   *
   */
  template <int dim, int spacedim>
  void
  hierarchical(DoFHandler<dim, spacedim> &dof_handler);

  /**
   * 按单元重新编号自由度。该函数接收一个单元格迭代器的向量（需要列出<i>all</i>DoF处理对象的本地拥有的活动单元格），并将根据自由度所在的单元格在给定的单元格列表中的位置给予自由度新的指数。存在于两个或多个单元间界面的自由度将在首先遇到时被编号。
   * 在同一单元上首先遇到的自由度在重新编号步骤前保留其原始顺序。
   * @param[in,out]  dof_handler
   * 其自由度将被重新编号的DoFHandler。    @param[in]  cell_order
   * 包含定义自由度应被重新编号的单元格的顺序的向量。
   * @pre 对于串行三角形 @p cell_order 必须有
   * <code>dof_handler.get_triangulation().n_active_cells()</code>
   * 的大小，而对于并行三角形，其大小应该是
   * parallel::TriangulationBase::n_locally_owned_active_cells().
   * 该三角形的每个活动单元迭代器需要正好出现在 @p
   * cell_order 中一次。
   *
   */
  template <int dim, int spacedim>
  void
  cell_wise(
    DoFHandler<dim, spacedim> &dof_handler,
    const std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator>
      &cell_order);

  /**
   * 按单元计算自由度的重新编号。该函数接收一个单元格迭代器的向量（需要列出<i>all</i>DoF处理对象的本地拥有的活动单元格），并将根据自由度所在的单元格在给定的单元格列表中的位置，给自由度以新的索引。存在于两个或多个单元间界面的自由度将在首先遇到时被编号。
   * 在同一单元上首先遇到的自由度在重新编号步骤之前保留其原始排序。
   * @param[out]  重新编号 一个长度为
   * <code>dof_handler.n_locally_owned_dofs()</code>
   * 的向量，包含每个自由度（以其当前编号）的未来DoF索引。因此，这个向量呈现了当前DoF指数的（非常特殊的）<i>permutation</i>。
   * @param[out]  inverse_renumbering
   * 前一个参数中返回的排列组合的反向。在
   * parallel::TriangulationBase
   * 的情况下，逆向是在本地拥有的DoF内。    @param[in]
   * dof_handler 其自由度要重新编号的DoFHandler。    @param[in]
   * cell_order
   * 包含定义自由度应被重新编号的单元的顺序的向量。
   * @pre  对于串行三角计算 @p cell_order 必须有大小
   * <code>dof_handler.get_triangulation().n_active_cells()</code>
   * ，而对于并行三角计算，其大小应该是
   * parallel::TriangulationBase::n_locally_owned_active_cells().
   * 该三角计算的每个活动单元迭代器需要在 @p
   * cell_order中精确出现一次。  @post  对于零和
   * <code>dof_handler.n_locally_owned_dofs()</code> 之间的每个 @p i
   * ，<code>renumbering[inverse_renumbering[i]] ==
   * dof_handler.local_owned_dofs().nth_index_in_set(i)</code>条件将成立。
   *
   */
  template <int dim, int spacedim>
  void
  compute_cell_wise(
    std::vector<types::global_dof_index> &renumbering,
    std::vector<types::global_dof_index> &inverse_renumbering,
    const DoFHandler<dim, spacedim> &     dof_handler,
    const std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator>
      &cell_order);

  /**
   * 像其他cell_wise()函数一样，但对于多级自由度枚举中的一级。
   *
   */
  template <int dim, int spacedim>
  void
  cell_wise(
    DoFHandler<dim, spacedim> &dof_handler,
    const unsigned int         level,
    const std::vector<typename DoFHandler<dim, spacedim>::level_cell_iterator>
      &cell_order);

  /**
   * 像其他compute_cell_wise()函数一样，但对多级自由度枚举中的一级而言。
   *
   */
  template <int dim, int spacedim>
  void
  compute_cell_wise(
    std::vector<types::global_dof_index> &renumbering,
    std::vector<types::global_dof_index> &inverse_renumbering,
    const DoFHandler<dim, spacedim> &     dof_handler,
    const unsigned int                    level,
    const std::vector<typename DoFHandler<dim, spacedim>::level_cell_iterator>
      &cell_order);

  /**
   * @}
   *
   */

  /**
   * @name 方向性编号  
     * @{ 
   *
   */

  /**
   * 相对于恒定流向的下游编号。如果附加参数 @p
   * dof_wise_renumbering 被设置为 @p false,
   * ，则按单元进行编号，否则根据支持点的位置进行编号。
   * 单元的排序是这样的：相对于常数向量 @p direction
   * ，编号较高的单元中心比编号较低的单元中心更靠下游。即使这在相当普遍的网格中产生了一个相对于边缘通量的下游编号，这也不一定能保证所有网格的通量。
   * 如果 @p dof_wise_renumbering 参数被设置为 @p false,
   * ，该函数会产生一个网格单元的下游排序，并调用cell_wise()。
   * 因此，在这种情况下，输出结果只对非连续加尔金有限元有意义（在这种情况下，所有的自由度都必须与单元的内部相关联）。
   * 如果 @p dof_wise_renumbering 被设置为 @p true,
   * ，自由度将根据各个自由度的支持点位置重新编号（显然，有限元需要定义支持点才能发挥作用）。下游位置相同的点的编号（例如那些与流动方向平行的点，或者一个FES系统内的几个自由度）将不受影响。
   *
   */
  template <int dim, int spacedim>
  void
  downstream(DoFHandler<dim, spacedim> &dof_handler,
             const Tensor<1, spacedim> &direction,
             const bool                 dof_wise_renumbering = false);


  /**
   * 相对于多网格层次结构中某一层的恒定流向，进行单元的下游编号。参见同名的其他函数。
   *
   */
  template <int dim, int spacedim>
  void
  downstream(DoFHandler<dim, spacedim> &dof_handler,
             const unsigned int         level,
             const Tensor<1, spacedim> &direction,
             const bool                 dof_wise_renumbering = false);

  /**
   * 计算downstream()函数所需的重新编号索引集。
   * 不对DoFHandler
   * dofs进行重新编号，但返回重新编号的向量。
   *
   */
  template <int dim, int spacedim>
  void
  compute_downstream(std::vector<types::global_dof_index> &new_dof_indices,
                     std::vector<types::global_dof_index> &reverse,
                     const DoFHandler<dim, spacedim> &     dof_handler,
                     const Tensor<1, spacedim> &           direction,
                     const bool dof_wise_renumbering);

  /**
   * 计算下游()函数所需的重新编号索引集。  不对DoFHandler
   * dofs进行重新编号，但返回重新编号的向量。
   *
   */
  template <int dim, int spacedim>
  void
  compute_downstream(std::vector<types::global_dof_index> &new_dof_indices,
                     std::vector<types::global_dof_index> &reverse,
                     const DoFHandler<dim, spacedim> &     dof_handler,
                     const unsigned int                    level,
                     const Tensor<1, spacedim> &           direction,
                     const bool dof_wise_renumbering);

  /**
   * 单元的顺时针编号。    这个函数产生一个相对于中心 @p
   * center
   * 的网格单元的（反）顺时针排序，并调用cell_wise()。
   * 因此，它只适用于不连续的Galerkin有限元，即所有的自由度都必须与单元的内部相关。
   *
   */
  template <int dim, int spacedim>
  void
  clockwise_dg(DoFHandler<dim, spacedim> &dof_handler,
               const Point<spacedim> &    center,
               const bool                 counter = false);

  /**
   * 在多网格层次结构的一个层面上进行单元的顺时针编号。参见同名的其他函数。
   *
   */
  template <int dim, int spacedim>
  void
  clockwise_dg(DoFHandler<dim, spacedim> &dof_handler,
               const unsigned int         level,
               const Point<spacedim> &    center,
               const bool                 counter = false);

  /**
   * 计算clockwise_dg()函数所需的重新编号向量。
   * 不对DoFHandler
   * dofs进行重新编号，但返回重新编号的向量。
   *
   */
  template <int dim, int spacedim>
  void
  compute_clockwise_dg(std::vector<types::global_dof_index> &new_dof_indices,
                       const DoFHandler<dim, spacedim> &     dof_handler,
                       const Point<spacedim> &               center,
                       const bool                            counter);

  /**
   * @}
   *
   */

  /**
   * @name  选择性的和随机的编号  
     * @{ 
   *
   */

  /**
   * 将那些在 @p selected_dofs数组中被标记为 @p true
   * 的自由度排序到自由度编号的后面。这种排序是稳定的，也就是说，被标记的自由度中的相对顺序被保留，未被标记的自由度中的相对顺序也被保留。
   * @pre  @p selected_dofs 数组的元素数量必须与 @p
   * Dof_handler的自由度一样多。
   *
   */
  template <int dim, int spacedim>
  void
  sort_selected_dofs_back(DoFHandler<dim, spacedim> &dof_handler,
                          const std::vector<bool> &  selected_dofs);

  /**
   * 将那些在 @p selected_dofs数组中被标记为 @p true
   * 的自由度排序到DoF数字的后面。
   * 这种排序是稳定的，也就是说，被标记的自由度内的相对顺序被保留，未被标记的自由度内的相对顺序也被保留。
   * @pre  @p selected_dofs 数组的元素数量必须与 @p
   * Dof_handler在给定层次上的自由度相同。
   *
   */
  template <int dim, int spacedim>
  void
  sort_selected_dofs_back(DoFHandler<dim, spacedim> &dof_handler,
                          const std::vector<bool> &  selected_dofs,
                          const unsigned int         level);

  /**
   * 计算sort_selected_dofs_back()函数所需的重新编号向量。不对DoFHandler
   * dofs进行重新编号，但返回重新编号的向量。      @pre  @p
   * selected_dofs 数组的元素数必须与 @p
   * dof_handler的自由度数一样多。
   *
   */
  template <int dim, int spacedim>
  void
  compute_sort_selected_dofs_back(
    std::vector<types::global_dof_index> &new_dof_indices,
    const DoFHandler<dim, spacedim> &     dof_handler,
    const std::vector<bool> &             selected_dofs);

  /**
   * 这个函数计算sort_selected_dofs_back()函数所需要的每一层的重新编号向量。不对DoFHandler
   * dofs进行重新编号，只计算重新编号并返回重新编号向量。
   * @pre  @p selected_dofs 数组的元素数量必须与 @p
   * Dof_handler在给定层次上的自由度相同。
   *
   */
  template <int dim, int spacedim>
  void
  compute_sort_selected_dofs_back(
    std::vector<types::global_dof_index> &new_dof_indices,
    const DoFHandler<dim, spacedim> &     dof_handler,
    const std::vector<bool> &             selected_dofs,
    const unsigned int                    level);

  /**
   * 以随机的方式对自由度重新编号。这个函数的结果是可重复的，即同一程序的两次运行将产生相同的结果。这是通过每次输入这个函数时用固定的种子创建一个新的随机数发生器来实现的。特别是，该函数因此不依赖于外部随机数生成器，因为它在该函数之前被调用的次数很重要（或者，就这一点而言，与该函数同时运行的其他线程是否也在绘制随机数）。
   *
   */
  template <int dim, int spacedim>
  void
  random(DoFHandler<dim, spacedim> &dof_handler);

  /**
   * 以随机的方式对自由度重新编号。它的作用与上述函数相同，只是它是针对多级离散化中的一个单级进行的。DoFHandler的非多网格部分不被触及。
   *
   */
  template <int dim, int spacedim>
  void
  random(DoFHandler<dim, spacedim> &dof_handler, const unsigned int level);

  /**
   * 计算随机()函数所需的重新编号向量。关于计算的随机重新编号的更多信息，见那里。
   * 这个函数并不对DoFHandler
   * dofs进行重新编号，而是返回重新编号的向量。
   *
   */
  template <int dim, int spacedim>
  void
  compute_random(std::vector<types::global_dof_index> &new_dof_indices,
                 const DoFHandler<dim, spacedim> &     dof_handler);

  /**
   * 计算随机()函数所需的重新编号向量。与上述函数相同，但用于多级离散化的单级。
   *
   */
  template <int dim, int spacedim>
  void
  compute_random(std::vector<types::global_dof_index> &new_dof_indices,
                 const DoFHandler<dim, spacedim> &     dof_handler,
                 const unsigned int                    level);

  /**
   * @}
   *
   */

  /**
   * @name  基于单元格属性的编号  
     * @{ 
   *
   */

  /**
   * 对自由度进行重新编号，使其与它们所在的单元格的子域ID相关联，即首先对属于子域0的单元格的所有自由度进行编号，然后是子域1的所有自由度，等等。这在使用分区器分配子域ID后用标准三角法进行并行计算时非常有用（见
   * GridTools::partition_triangulation 函数）。当使用
   * parallel::shared::Triangulation 或 parallel::distributed::Triangulation,
   * 时，调用这个函数是不必要的，因为自由度已经根据MPI进程ID被列举出来。因此，如果底层三角结构是这种类型，那么将产生一个错误。
   * 请注意，与面、边和顶点相关的自由度可能与多个子域相关，如果它们位于分区边界。因此，它们必须与哪个子域相关联是不确定的。对于这一点，我们使用从
   * DoFTools::get_subdomain_association 函数中得到的东西。
   * 该算法是稳定的，即如果两个道夫i,j有<tt>i<j</tt>并且属于同一个子域，那么它们在重新排序后也将是这个顺序。
   *
   */
  template <int dim, int spacedim>
  void
  subdomain_wise(DoFHandler<dim, spacedim> &dof_handler);

  /**
   * 计算subdomain_wise()函数所需的重新编号的向量。  不对 @p
   * DoFHandler Dofs进行重新编号，但返回重新编号的向量。
   *
   */
  template <int dim, int spacedim>
  void
  compute_subdomain_wise(std::vector<types::global_dof_index> &new_dof_indices,
                         const DoFHandler<dim, spacedim> &     dof_handler);

  /**
   * @}
   *
   */



  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(ExcDoFHandlerNotInitialized,
                   "The DoFHandler on which this function should work has not "
                   "been initialized, i.e., it doesn't appear that DoF indices "
                   "have been distributed on it.");

  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcInvalidComponentOrder);

  /**
   * 该函数仅对不连续的Galerkin有限元素实现。
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcNotDGFEM);
} // namespace DoFRenumbering


DEAL_II_NAMESPACE_CLOSE

#endif


