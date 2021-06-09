//include/deal.II-translator/lac/sparsity_tools_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2020 by the deal.II authors
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

#ifndef dealii_sparsity_tools_h
#define dealii_sparsity_tools_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <memory>
#include <vector>

#ifdef DEAL_II_WITH_MPI
#  include <deal.II/base/index_set.h>

#  include <mpi.h>
#endif

DEAL_II_NAMESPACE_OPEN


/*!   @addtogroup Sparsity  
     * @{  

* 
*
*/

/**
 * 一个命名空间，用于处理可以在稀疏性模式上做的事情的函数，例如根据连通性对行和列（或者自由度，如果你想的话）重新编号，或者对自由度进行分区。
 *
 *
 */
namespace SparsityTools
{
  /**
   * 带有分区器选项的枚举器
   *
   */
  enum class Partitioner
  {
    /**
     * 使用METIS分区器。
     *
     */
    metis = 0,
    /**
     * 使用ZOLTAN分区器。
     *
     */
    zoltan
  };


  /**
   * 使用一个分区算法来生成这个稀疏模式所代表的自由度的分区。实际上，我们把这种稀疏模式看作是各种自由度之间的连接图，其中稀疏模式中的每个非零条目对应于连接图中两个节点之间的一条边。然后，我们的目标是将这个图分解成节点组，使节点组之间的边界切割出最小数量的边。这种分割是由METIS或ZOLTAN完成的，这取决于第四个参数中选择的分割器。默认是METIS。请注意，METIS和ZOLTAN只能划分对称的稀疏模式，当然，稀疏模式必须是方形的。我们不检查稀疏性模式的对称性，因为这是一个昂贵的操作，而是让这个函数的调用者负责。
   * 调用这个函数后，输出数组的每个节点（即矩阵的行或列）的值将在0和
   * @p n_partitions-1 之间。
   * 如果deal.II没有与ZOLTAN或METIS软件包一起安装，当选择相应的分区方法时，该函数将产生一个错误，除非
   * @p n_partitions
   * 是一个。也就是说，你可以写一个程序，使其在单处理器单分区的情况下运行，而不安装软件包，只有在需要多分区时才需要安装它们。
   * 注意，稀疏模式本身不会因为调用这个函数而改变。然而，你可能会使用调用这个函数产生的信息来重新编号自由度，之后你当然要重新生成稀疏模式。
   * 这个函数很少被单独调用，因为在有限元方法中，你想对网格进行分割，而不是对矩阵进行分割。这可以通过调用
   * @p GridTools::partition_triangulation. 来实现。
   *
   */
  void
  partition(const SparsityPattern &    sparsity_pattern,
            const unsigned int         n_partitions,
            std::vector<unsigned int> &partition_indices,
            const Partitioner          partitioner = Partitioner::metis);


  /**
   * 这个函数执行的操作与上面的函数相同，只是它考虑到了一组
   * @p cell_weights, ，这些 @p cell_weights,
   * 允许分区器在考虑到每个单元上所花费的计算力的同时平衡图形。
   * @note  如果 @p cell_weights
   * 向量为空，则不考虑加权。如果不是，那么这个向量的大小必须等于三角形中有效单元的数量。
   *
   */
  void
  partition(const SparsityPattern &          sparsity_pattern,
            const std::vector<unsigned int> &cell_weights,
            const unsigned int               n_partitions,
            std::vector<unsigned int> &      partition_indices,
            const Partitioner                partitioner = Partitioner::metis);

  /**
   使用ZOLTAN提供的着色算法为节点着色，这些节点的连接是用SparsityPattern对象表示的。实际上，我们把这个稀疏模式看作是一个节点之间的连接图，其中 @p sparsity_pattern 中的每个非零条目对应于两个节点之间的一条边。我们的目标是给每个节点分配一个颜色指数，使得没有两个直接连接的节点具有相同的颜色。  所分配的颜色在 @p color_indices 中列出，从一开始就有索引，该函数返回所用颜色的总数。ZOLTAN着色算法是由每个处理器串行运行的。因此，所有处理器都有所有节点的着色信息。GraphColoring命名空间中也有这个函数的封装函数。    请注意，这个函数要求SparsityPattern是对称的（因此是方形的），即无向图表示。我们不检查稀疏性模式的对称性，因为这是一个昂贵的操作，而是让这个函数的调用者负责。      @image html color_undirected_graph.png  上图说明了该函数的用法，显示了五个从0到4编号的节点，显示的连接是双向的。  着色后，很明显，没有两个直接连接的节点被赋予相同的颜色。    如果deal.II没有和包ZOLTAN一起安装，这个函数将产生一个错误。
   * @note 当前函数是 GraphColoring::make_graph_coloring()
   * 的替代品，该函数是为矩阵的共享内存并行装配中产生的图形着色而定制的。
   *
   */
  unsigned int
  color_sparsity_pattern(const SparsityPattern &    sparsity_pattern,
                         std::vector<unsigned int> &color_indices);

  /**
   * 对于一个给定的稀疏模式，根据Cuthill-McKee的算法计算行/列索引的重新列举。
   * 该算法是一种图的重新编号算法，我们试图根据图中所有节点与其他节点的连接性（即连接节点的边）找到一个新的编号。这种连接性在这里是由稀疏模式表示的。在库中的许多情况下，节点代表自由度，边是矩阵中的非零条目，即通过双线性形式的作用耦合的自由度对。
   * 该算法从一个节点开始，在其他节点中寻找那些与我们开始的那个节点相耦合的节点，并以某种方式为这些节点编号。然后它找到第二层的节点，即那些与上一层的节点（即与初始节点耦合的节点）耦合的节点，并对这些节点编号。以此类推。关于该算法的细节，特别是每一层内的编号，我们请读者参考施瓦茨的书（H.
   * R. Schwarz: Methode der finiten Elemente）。
   * 这些算法有一个主要的缺点：它们需要一个好的起始节点，即在输出阵列中数字为零的节点。因此，形成初始层节点的起始节点可以由用户给出，例如，通过利用领域的实际拓扑知识。
   * 也可以给出几个起始指数，可以用来模拟简单的上游编号（通过给流入节点作为起始值）或使预处理更快（通过让Dirichlet边界指数作为起始点）。
   * 如果没有给出起始指数，就会自动选择一个，即协调数最小的一个节点（协调数是这个节点与其他节点的耦合数）。这个节点通常位于域的边界上。然而，在使用本库中使用的分层网格时，这一点存在很大的模糊性，因为在大多数情况下，计算域不是通过倾斜和变形的元素以及在顶点插入可变数量的元素来逼近的，而是通过分层细化来实现的。因此有大量的节点具有相等的协调数。因此，重新编号的算法将不会得到最佳结果。
   * 如果图有两个或更多的不相连的组件，如果没有给出起始指数，算法将对每个组件连续编号。然而，这需要为每个组件确定一个起始索引；因此，如果给出了起始索引，算法将产生一个异常，认为后者表明函数的调用者希望覆盖算法中选择起始索引的部分。
   *
   */
  void
  reorder_Cuthill_McKee(
    const DynamicSparsityPattern &                        sparsity,
    std::vector<DynamicSparsityPattern::size_type> &      new_indices,
    const std::vector<DynamicSparsityPattern::size_type> &starting_indices =
      std::vector<DynamicSparsityPattern::size_type>());

  /**
   * 对于一个给定的稀疏模式，以分层的方式计算行/列索引的重新列举，类似于
   * DoFRenumbering::hierarchical
   * 对分层细化网格上的自由度的计算。
   * 该算法首先选择一个邻居数量最少的节点，并将该节点和它的直接邻居放入一个块中。接下来，它选择已经选择的节点的一个邻居，将该节点和它的直接邻居中不属于前一个块的节点，加入到下一个块中。在这次扫描之后，相邻的节点被分组在一起。
   * 为了确保在更大的范围内进行类似的分组，这种分组在如此形成的组上被递归调用。当不可能进一步分组时，递归停止。最终，通过这种方法得到的排序以类似于Z的方式通过稀疏模式中的指数。
   * 如果该图有两个或更多的不相连的组件，该算法将对每个组件连续编号，从节点数最少的组件开始。
   *
   */
  void
  reorder_hierarchical(
    const DynamicSparsityPattern &                  sparsity,
    std::vector<DynamicSparsityPattern::size_type> &new_indices);

#ifdef DEAL_II_WITH_MPI
  /**
   * 通过MPI以动态稀疏模式进行行的通信。      @param[in,out]
   * dsp
   * 一个在本地建立的动态稀疏模式，我们需要与其他处理器交换条目，以确保每个处理器知道它所存储的、最终可能被写入的矩阵的所有行的元素。这种稀疏模式将因这个函数而改变。属于不同处理器的行中的所有条目都会被发送到他们那里并添加到那里。
   * @param  locally_owned_rows
   * 一个描述由调用MPI进程拥有的行的索引集。该索引集在通信器中的处理器之间应是一对一的。
   * @param  mpi_comm
   * 参与此操作的处理器之间共享的MPI通信器。      @param
   * locally_relevant_rows
   * 存储在本地MPI进程上的元素范围。这应该是在DynamicSparsityPattern的构造函数中使用的，也应该是本地相关集。只有包含在这个集合中的行才会在dsp中被检查以进行传输。这个函数需要和
   * PETScWrappers::MPI::SparseMatrix
   * 一起使用，才能在并行计算中正确工作。
   *
   */
  void
  distribute_sparsity_pattern(DynamicSparsityPattern &dsp,
                              const IndexSet &        locally_owned_rows,
                              const MPI_Comm &        mpi_comm,
                              const IndexSet &        locally_relevant_rows);

  /**
   * 通过MPI以动态稀疏模式交流行，与上面的类似，但使用一个向量`rows_per_cpu`，包含每个CPU的行数来确定所有权。这通常是由
   * Utilities::MPI::all_gather(MPI_Comm,   DoFHandler::locally_owned_dofs())
   * 返回的值。
   *
   * - 考虑到这个函数的输入构造涉及到全对全的通信，对于超过一千的进程，它通常比上面的函数慢（对于小规模的进程也足够快）。
   *
   */
  void
  distribute_sparsity_pattern(
    DynamicSparsityPattern &                              dsp,
    const std::vector<DynamicSparsityPattern::size_type> &rows_per_cpu,
    const MPI_Comm &                                      mpi_comm,
    const IndexSet &                                      myrange);

  /**
   * 类似于上面的函数，但对BlockDynamicSparsityPattern而言，则是。
   * @param[in,out]  dsp 要修改的本地构建的稀疏性模式。
   * @param  locally_owned_rows
   * 描述由调用MPI进程拥有的行的索引集。该索引集在通信器中的处理器之间应是一对一的。
   * @param  mpi_comm 要使用的MPI通信器。      @param
   * locally_relevant_rows 通常是本地相关的DoF。
   *
   */
  void
  distribute_sparsity_pattern(BlockDynamicSparsityPattern &dsp,
                              const IndexSet &             locally_owned_rows,
                              const MPI_Comm &             mpi_comm,
                              const IndexSet &locally_relevant_rows);

  /**
   * @deprecated
   * 使用distribut_sparsity_pattern()，只为当前的MPI进程设置一个索引。
   *
   */
  void
  distribute_sparsity_pattern(BlockDynamicSparsityPattern &dsp,
                              const std::vector<IndexSet> &owned_set_per_cpu,
                              const MPI_Comm &             mpi_comm,
                              const IndexSet &             myrange);

  /**
   * 通过MPI以动态稀疏模式收集行。  该函数类似于
   * SparsityTools::distribute(),
   * ，但是该函数将从其他MPI进程中收集稀疏度，并将其添加到本地DynamicSparsityPattern中，而不是在这个MPI进程中分配存储在非所有权行中的稀疏度。
   * @param  dsp
   * 一个已经在本地建立的动态稀疏模式，我们需要根据存储在其他MPI进程上的行的稀疏度来扩展它。
   * @param  locally_owned_rows
   * 一个描述由调用MPI进程拥有的行的索引集。该索引集在通信器中的处理器之间应是一对一的。
   * @param  mpi_comm
   * 参与此操作的处理器之间共享的MPI通信器。      @param
   * local_relevant_rows
   * 这个MPI进程需要收集的行的范围。只有不包括在本地拥有的行中的部分将被使用。
   *
   */
  void
  gather_sparsity_pattern(DynamicSparsityPattern &dsp,
                          const IndexSet &        locally_owned_rows,
                          const MPI_Comm &        mpi_comm,
                          const IndexSet &        locally_relevant_rows);

  /**
   * @deprecated
   * 使用gather_sparsity_pattern()方法，只为当前处理器设置索引。
   *
   */
  DEAL_II_DEPRECATED void
  gather_sparsity_pattern(DynamicSparsityPattern &     dsp,
                          const std::vector<IndexSet> &owned_rows_per_processor,
                          const MPI_Comm &             mpi_comm,
                          const IndexSet &             ghost_range);

#endif


  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcMETISNotInstalled,
                   "The function you called requires METIS, but you did not "
                   "configure deal.II with METIS.");

  /**
   * 异常情况
   *
   */
  DeclException1(ExcInvalidNumberOfPartitions,
                 int,
                 << "The number of partitions you gave is " << arg1
                 << ", but must be greater than zero.");

  /**
   * 异常情况
   *
   */
  DeclException1(ExcMETISError,
                 int,
                 << "    An error with error number " << arg1
                 << " occurred while calling a METIS function");

  /**
   * 异常情况
   *
   */
  DeclException2(ExcInvalidArraySize,
                 int,
                 int,
                 << "The array has size " << arg1 << " but should have size "
                 << arg2);
  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(
    ExcZOLTANNotInstalled,
    "The function you called requires ZOLTAN, but you did not "
    "configure deal.II with ZOLTAN or zoltan_cpp.h is not available.");
} // namespace SparsityTools

/**
 * @}
 *
 *
 */

DEAL_II_NAMESPACE_CLOSE

#endif


