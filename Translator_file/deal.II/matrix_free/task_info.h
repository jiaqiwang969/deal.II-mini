//include/deal.II-translator/matrix_free/task_info_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2021 by the deal.II authors
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


#ifndef dealii_matrix_free_task_info_h
#define dealii_matrix_free_task_info_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>


DEAL_II_NAMESPACE_OPEN



namespace internal
{
  /**
   * 工作者对象的接口，该对象在无矩阵循环中运行我们要执行的各种操作。
   *
   */
  struct MFWorkerInterface
  {
  public:
    virtual ~MFWorkerInterface() = default;

    /// Starts the communication for the update ghost values operation
    virtual void
    vector_update_ghosts_start() = 0;

    /// Finishes the communication for the update ghost values operation
    virtual void
    vector_update_ghosts_finish() = 0;

    /// Starts the communication for the vector compress operation
    virtual void
    vector_compress_start() = 0;

    /// Finishes the communication for the vector compress operation
    virtual void
    vector_compress_finish() = 0;

    /// Zeros part of the vector according to a given range as stored in
    /// DoFInfo
    virtual void
    zero_dst_vector_range(const unsigned int range_index) = 0;

    virtual void
    cell_loop_pre_range(const unsigned int range_index) = 0;

    virtual void
    cell_loop_post_range(const unsigned int range_index) = 0;

    /// Runs the cell work specified by MatrixFree::loop or
    /// MatrixFree::cell_loop
    virtual void
    cell(const std::pair<unsigned int, unsigned int> &cell_range) = 0;

    /// Runs the cell work specified by MatrixFree::loop or
    /// MatrixFree::cell_loop
    virtual void
    cell(const unsigned int range_index) = 0;

    /// Runs the body of the work on interior faces specified by
    /// MatrixFree::loop
    virtual void
    face(const unsigned int range_index) = 0;

    /// Runs the body of the work on boundary faces specified by
    /// MatrixFree::loop
    virtual void
    boundary(const unsigned int range_index) = 0;
  };



  namespace MatrixFreeFunctions
  {
    /**
     * 一个收集与线程并行化有关的所有信息的结构。工作被细分为可以独立完成的任务。
     *
     */
    struct TaskInfo
    {
      // enum for choice of how to build the task graph. Odd add versions with
      // preblocking and even versions with postblocking. partition_partition
      // and partition_color are deprecated but kept for backward
      // compatibility.
      enum TasksParallelScheme
      {
        none,
        partition_partition,
        partition_color,
        color
      };

      /**
       * 构造函数。
       *
       */
      TaskInfo();

      /**
       * 清除所有的数据字段并将其重置为零。
       *
       */
      void
      clear();

      /**
       * 运行无矩阵循环。
       *
       */
      void
      loop(MFWorkerInterface &worker) const;

      /**
       * 使只能在通信重叠中处理的单元格的数量被矢量化长度所除。
       *
       */
      void
      make_boundary_cells_divisible(std::vector<unsigned int> &boundary_cells);

      /**
       * 根据输入参数控制的选项，设置运行单元格循环的块。
       * @param  cells_with_comm
       * 一个需要在执行计算前交换数据的单元的列表。这些将在分区中被赋予一定的id，以确保与通信重叠的单元格循环有准备好的幽灵数据。
       * @param  dofs_per_cell
       * 给出一个单元上自由度数量的预期值，用于确定交错单元和面积分的块大小。
       * @param  categories_are_hp
       * 定义`cell_vectorization_categories`是源于具有可变多项式程度的hp自适应计算还是用户定义的变体。
       * @param  cell_vectorization_categories
       * 这组类别定义了在向量数组的通道内应该被分组的单元。这可以是hp-元素中的多项式程度，也可以是用户提供的分组。
       * @param  cell_vectorization_categories_strict
       * 定义前面的变量所定义的类别是否应该被严格分开，或者是否允许将较低的类别插入到下一个较高的类别。
       * @param  parent_relation
       * 这个数据字段用于指定哪些单元格有相同的父单元格。具有相同祖先的单元格被分组到同一批次（es），并在单元格之间进行矢量处理。
       * @param  重新编号
       * 当离开这个函数时，向量包含一个新的单元格编号，与存储在这个类中的分组一致。
       * @param  incompletely_filled_vectorization
       * 鉴于本类的矢量布局，一些单元格批次可能在矢量数组（SIMD通道）中有未被使用的组件，并且不承载有效数据。这个数组根据这个函数返回的重新编号，指出发生这种情况的单元批处理。
       *
       */
      void
      create_blocks_serial(
        const std::vector<unsigned int> &cells_with_comm,
        const unsigned int               dofs_per_cell,
        const bool                       categories_are_hp,
        const std::vector<unsigned int> &cell_vectorization_categories,
        const bool                       cell_vectorization_categories_strict,
        const std::vector<unsigned int> &parent_relation,
        std::vector<unsigned int> &      renumbering,
        std::vector<unsigned char> &     incompletely_filled_vectorization);

      /**
       * 任务并行分块设置的分块创建的第一步。
       * @param  boundary_cells
       * 一个需要在执行计算前交换数据的单元格的列表。这些单元在分区中会被赋予一定的ID。
       * @param  重新编号
       * 当离开这个函数时，向量包含一个新的单元格编号，与存储在这个类中的分组一致（在实际创建任务之前）。
       * @param  incompletely_filled_vectorization
       * 考虑到这个类的矢量布局，一些单元格批次可能在矢量数组（SIMD通道）中有未被使用的组件，并不承载有效数据。这个数组根据这个函数返回的重新编号，指出发生这种情况的单元批处理。
       *
       */
      void
      initial_setup_blocks_tasks(
        const std::vector<unsigned int> &boundary_cells,
        std::vector<unsigned int> &      renumbering,
        std::vector<unsigned char> &     incompletely_filled_vectorization);

      /**
       * 如果用户决定不通过 MatrixFree::AdditionalData.
       * 强制确定一个块的大小，这个辅助函数确定一个块的大小，这是根据系统的硬件线程数和我们应该工作的宏单元数计算出来的。
       *
       */
      void
      guess_block_size(const unsigned int dofs_per_cell);

      /**
       * 这个方法通过所有已经填入 @p
       * dof_indices的单元，找出哪些单元可以独立工作，哪些单元是相邻的，在并行使用时需要在不同时间进行。
       * 该策略是基于两层的方法。外层被细分为类似于Cuthill-McKee中邻居类型的分区，内层则通过颜色进行细分（对于同一颜色内的块，可以独立工作）。一个任务由一个单元块来代表。单元块在细分为分区和颜色之前就已经形成。
       * @param  连通性（in/out）
       * 判断单元格`i`和`j`是否冲突，用位置（i,j）的条目表示。
       * @param  重新编号（in/out）
       * 在输出时，该变量的元素j给出了单元格的原始编号，由于线程图的排序而被重新排序到j的位置。
       * @param  irregular_cells (in/out)
       * 通知当前函数，对于给定的单元格批次索引，VectorizedArray中的一些SIMD通道是否不会被填充。
       * @param  hp_bool 定义我们是否在hp模式下。
       *
       */
      void
      make_thread_graph_partition_color(
        DynamicSparsityPattern &    connectivity,
        std::vector<unsigned int> & renumbering,
        std::vector<unsigned char> &irregular_cells,
        const bool                  hp_bool);

      /**
       * 这个函数会浏览所有已经填入 @p
       * dof_indices的单元格，找出哪些单元格可以独立工作，哪些单元格是相邻的，在并行使用时需要在不同时间进行。
       * 该策略是基于两层的方法。外层被细分为类似于Cuthill-McKee中邻居类型的分区，内层又被细分为类似Cuthill-McKee的分区（级别相差2以上的分区可以独立工作）。一个任务由一个单元块来代表。细胞块是在细分为两级分区后形成的。
       * @param  cell_active_fe_index
       * 与所有单元格索引列表中的单个索引相对应的活动FE索引，以便能够不将具有不同索引的单元格放置在具有矢量化的同一单元格批次中。
       * @param  connectivity (in/out)
       * 确定单元格`i`和`j`是否冲突，由位置（i,j）的条目表示。
       * @param  重新编号（in/out）
       * 在输出时，该变量的元素j给出了单元格的原始编号，由于线程图的排序而被重新排序到j的位置。
       * @param  irregular_cells (in/out)
       * 通知当前函数，对于给定的单元格批次索引，VectorizedArray中的一些SIMD通道是否不会被填充。
       * @param  hp_bool 定义我们是否在hp模式下。
       *
       */
      void
      make_thread_graph_partition_partition(
        const std::vector<unsigned int> &cell_active_fe_index,
        DynamicSparsityPattern &         connectivity,
        std::vector<unsigned int> &      renumbering,
        std::vector<unsigned char> &     irregular_cells,
        const bool                       hp_bool);

      /**
       * 要么调用make_thread_graph_partition_color()，要么调用外部可访问的make_thread_graph_partition_partition()，这取决于数据结构中的设置。
       * @param  cell_active_fe_index
       * 与所有单元格索引列表中的单个索引相对应的活动FE索引，以便能够不将具有不同索引的单元格放入同一个具有矢量的单元格批次中。
       * @param  connectivity (in/out)
       * 确定单元格`i`和`j`是否冲突，用位置（i,j）的条目表示。
       * @param  重新编号（in/out）
       * 在输出时，该变量的元素j给出了单元格的原始编号，由于线程图的排序而被重新排序到j的位置。
       * @param  irregular_cells (in/out)
       * 通知当前函数，对于给定的单元格批次索引，VectorizedArray中的一些SIMD通道是否不会被填充。
       * @param  hp_bool 定义我们是否在hp模式下。
       *
       */
      void
      make_thread_graph(const std::vector<unsigned int> &cell_active_fe_index,
                        DynamicSparsityPattern &         connectivity,
                        std::vector<unsigned int> &      renumbering,
                        std::vector<unsigned char> &     irregular_cells,
                        const bool                       hp_bool);

      /**
       * 这个函数从单个细胞之间的连接性计算出细胞块之间的连接性。
       *
       */
      void
      make_connectivity_cells_to_blocks(
        const std::vector<unsigned char> &irregular_cells,
        const DynamicSparsityPattern &    connectivity_cells,
        DynamicSparsityPattern &          connectivity_blocks) const;

      /**
       * 在每个分区内的第二层上创建着色的%函数。
       *
       */
      void
      make_coloring_within_partitions_pre_blocked(
        const DynamicSparsityPattern &   connectivity,
        const unsigned int               partition,
        const std::vector<unsigned int> &cell_partition,
        const std::vector<unsigned int> &partition_list,
        const std::vector<unsigned int> &partition_size,
        std::vector<unsigned int> &      partition_color_list);

      /**
       * 在每个分区的第二层上创建分区的函数。
       *
       */
      void
      make_partitioning_within_partitions_post_blocked(
        const DynamicSparsityPattern &   connectivity,
        const std::vector<unsigned int> &cell_active_fe_index,
        const unsigned int               partition,
        const unsigned int               cluster_size,
        const bool                       hp_bool,
        const std::vector<unsigned int> &cell_partition,
        const std::vector<unsigned int> &partition_list,
        const std::vector<unsigned int> &partition_size,
        std::vector<unsigned int> &      partition_partition_list,
        std::vector<unsigned char> &     irregular_cells);

      /**
       * 这个函数根据提供的连接图创建分区。
       * @param  connectivity （细胞块）之间的连通性  @param
       * cluster_size
       * 每个分区的细胞数应该是cluster_size的倍数（用于以后的阻塞）
       * @param  cell_partition
       * 保存每个（细胞块）的分区，该块属于哪个分区
       * @param  ] partition_list partition_list[j]
       * 给出了由于分区而应该重新编号为j的块的旧编号
       * @param  partition_size
       * 指向每个分区的开始的矢量（在输出时）  @param
       * partition 创建的分区数
       *
       */
      void
      make_partitioning(const DynamicSparsityPattern &connectivity,
                        const unsigned int            cluster_size,
                        std::vector<unsigned int> &   cell_partition,
                        std::vector<unsigned int> &   partition_list,
                        std::vector<unsigned int> &   partition_size,
                        unsigned int &                partition) const;

      /**
       * 为make_thread_graph中设置的任务图更新任务信息的字段。
       *
       */
      void
      update_task_info(const unsigned int partition);

      /**
       * 从连接结构中创建一个任务图。
       *
       */
      void
      create_flow_graph();

      /**
       * 返回该类的内存消耗。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 打印MPI进程的最小、平均和最大的内存消耗。
       *
       */
      template <typename StreamType>
      void
      print_memory_statistics(StreamType &out, std::size_t data_length) const;

      /**
       * 网格中的物理单元的数量，而不是矢量化后的单元批次
       *
       */
      unsigned int n_active_cells;

      /**
       * 网格中的物理幽灵单元的数量，这些单元要进行特殊处理，不应包括在循环中。
       *
       */
      unsigned int n_ghost_cells;

      /**
       * SIMD阵列中用于矢量化的通道数量
       *
       */
      unsigned int vectorization_length;

      /**
       * 用于多线程的块大小信息
       *
       */
      unsigned int block_size;

      /**
       * 用于多线程的块的数量
       *
       */
      unsigned int n_blocks;

      /**
       * 多线程应用的并行方案
       *
       */
      TasksParallelScheme scheme;

      /**
       * 块是由矢量的概念组织的，这个数据字段 @p
       * partition_row_index 在所有数据的线性存储中存储了一个
       * "矢量 "到下一个的距离，以两级划分。
       *
       */
      std::vector<unsigned int> partition_row_index;

      /**
       * 这是所有分区的线性存储，在MatrixFree中所有单元的整数列表中建立一个形式为cell_partition_data[idx]到cell_partition_data[idx+1]的索引范围，通过
       * @p partition_row_index. 细分为大块。
       *
       */
      std::vector<unsigned int> cell_partition_data;

      /**
       * 像cell_partition_data一样，但对每个活动的fe索引有预先计算的子范围。分区的起点和终点由cell_partition_data_hp_ptr给出。
       *
       */
      std::vector<unsigned int> cell_partition_data_hp;

      /**
       * cell_partition_data_hp内的指针，表示一个分区的开始和结束。
       *
       */
      std::vector<unsigned int> cell_partition_data_hp_ptr;

      /**
       * 这是对所有内面的分区的线性存储，在MatrixFree中所有内面的整数列表中建立一个形式为face_partition_data[idx]到face_partition_data[idx+1]的索引范围，按
       * @p partition_row_index细分为大块。
       *
       */
      std::vector<unsigned int> face_partition_data;

      /**
       * 像face_partition_data一样，但对每个活动fe索引对都有预先计算的子范围。分区的起点和终点由face_partition_data_hp_ptr给出。
       *
       */
      std::vector<unsigned int> face_partition_data_hp;

      /**
       * face_partition_data_hp内的指针，表示一个分区的开始和结束。
       *
       */
      std::vector<unsigned int> face_partition_data_hp_ptr;

      /**
       * 这是对边界面所有分区的线性存储，在MatrixFree中所有边界面的整数列表中建立一个形式为boundary_partition_data[idx]到boundary_partition_data[idx+1]的索引范围，按
       * @p partition_row_index细分为大块。
       *
       */
      std::vector<unsigned int> boundary_partition_data;

      /**
       * 像boundary_partition_data一样，但对每个活动的fe索引有预先计算的子范围。分区的起点和终点是由boundary_partition_data_hp_ptr给出的。
       *
       */
      std::vector<unsigned int> boundary_partition_data_hp;

      /**
       * boundary_partition_data_hp内的指针，表示一个分区的开始和结束。
       *
       */
      std::vector<unsigned int> boundary_partition_data_hp_ptr;

      /**
       * 这是对边界上所有内部面的分区的线性存储，以其他处理器不在本地使用的方式，在MatrixFree中所有此类面的整数列表中建立一个形式为ghost_face_partition_data[idx]到ghost_face_partition_data[idx+1]的索引范围，按
       * @p  partition_row_index细分为大块。
       *
       */
      std::vector<unsigned int> ghost_face_partition_data;

      /**
       * 这是一个线性存储多网格级别的所有面的分区，这些面有一个更粗的邻居，只包括在某些剩余计算中，但不包括在平滑中，在MatrixFree中所有这些面的整数列表中建立一个形式为refinement_edge_face_partition_data[idx]到refinement_edge_face_partition_data[idx+1]的索引范围，按
       * @p  partition_row_index细分为几块。
       *
       */
      std::vector<unsigned int> refinement_edge_face_partition_data;

      /**
       * 将交给动态任务调度器的线程信息（从哪个块开始
       * "偶数 "分区）。
       *
       */
      std::vector<unsigned int> partition_evens;

      /**
       * 将线程信息（从哪个区块启动 "奇数
       * "分区）交给动态任务调度器
       *
       */
      std::vector<unsigned int> partition_odds;

      /**
       * 关于移交给动态任务调度器的分区的依赖性的线程信息
       *
       */
      std::vector<unsigned int> partition_n_blocked_workers;

      /**
       * 关于移交给动态任务调度器的分区的依赖关系的线程信息
       *
       */
      std::vector<unsigned int> partition_n_workers;

      /**
       * 在字段 @p 中积累的偶数分区的数量partitions_even
       *
       */
      unsigned int evens;

      /**
       * 场上累计的奇数分区数量  @p  partitions_odd
       *
       */
      unsigned int odds;

      /**
       * 在该领域积累的被封锁工人的数量  @p
       * partition_n_blocked_workers
       *
       */
      unsigned int n_blocked_workers;

      /**
       * 在该领域积累的工人数量  @p partition_n_workers  。
       *
       */
      unsigned int n_workers;

      /**
       * 存储一个特定的任务是否处于MPI边界并需要数据交换
       *
       */
      std::vector<unsigned char> task_at_mpi_boundary;

      /**
       * MPI通信器
       *
       */
      MPI_Comm communicator;

      /**
       * 共享内存的MPI通信器
       *
       */
      MPI_Comm communicator_sm;

      /**
       * MPI进程的等级
       *
       */
      unsigned int my_pid;

      /**
       * 当前通信器的MPI等级的数量
       *
       */
      unsigned int n_procs;
    };

  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif


