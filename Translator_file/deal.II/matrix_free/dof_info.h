//include/deal.II-translator/matrix_free/dof_info_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2020 by the deal.II authors
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


#ifndef dealii_matrix_free_dof_info_h
#define dealii_matrix_free_dof_info_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/matrix_free/face_info.h>
#include <deal.II/matrix_free/mapping_info.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/task_info.h>
#include <deal.II/matrix_free/vector_data_exchange.h>

#include <array>
#include <memory>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * 一个结构，它接收描述约束的条目，并将其放入一个排序的列表中，其中重复的内容被过滤掉。
     *
     */
    template <typename Number>
    struct ConstraintValues
    {
      ConstraintValues();

      /**
       * 这个函数在所有值的集合中插入一些受约束的条目。它将道夫的（重新排序的）编号（根据与函数相匹配的排序）存储在new_indices中，并返回存储位置的双数组，以便以后访问。
       *
       */
      template <typename number2>
      unsigned short
      insert_entries(
        const std::vector<std::pair<types::global_dof_index, number2>>
          &entries);

      std::vector<std::pair<types::global_dof_index, double>>
                                           constraint_entries;
      std::vector<types::global_dof_index> constraint_indices;

      std::pair<std::vector<Number>, types::global_dof_index> next_constraint;
      std::map<std::vector<Number>,
               types::global_dof_index,
               FPArrayComparator<Number>>
        constraints;
    };

    /**
     * 存储所有单元格自由度的索引的类。本质上，这是一个DoFHandler风格的智能数字缓存，同时也将约束的描述直接嵌入到单元格层面，而不需要引用外部AffineConstraints对象。
     * 这个类只存储索引关系。悬挂节点约束的权重被存储在一个不同的字段中。这是因为不同的字段允许在不同的DoFHandlers上对矢量值问题进行相同的压缩权重数据。在那里，指数可能在不同的分量上受到不同的约束（例如，Dirichlet条件只在选定的分量上），而来自悬挂节点的权重是相同的，只需要存储一次。这种组合将在MatrixFree类中处理。
     * @ingroup matrixfree
     *
     */
    struct DoFInfo
    {
      /**
       * 这个值用来定义向量中的子范围，我们可以在
       * MatrixFree::loop()
       * 的调用中归零。我们的目标是每次只清除向量的一部分，以保持缓存中被清零的值，为应用这个方法的情况节省一个全局向量访问，而不是`vector
       * = 0. ;`。            我们将颗粒度设置为64
       *
       * - 这是一个足够大的数字，以尽量减少工作中的循环剥离开销（并与高达16的矢量化长度兼容），并且足够小，不会浪费在单个块的大小上。
       *
       */
      static const unsigned int chunk_size_zero_vector = 64;

      /**
       * 默认的空构造函数。
       *
       */
      DoFInfo();

      /**
       * 复制构造函数。
       *
       */
      DoFInfo(const DoFInfo &) = default;

      /**
       * 移动构造函数。
       *
       */
      DoFInfo(DoFInfo &&) noexcept = default;

      /**
       * 解除构造器。
       *
       */
      ~DoFInfo() = default;

      /**
       * 复制赋值操作符。
       *
       */
      DoFInfo &
      operator=(const DoFInfo &) = default;

      /**
       * 移动赋值运算符。
       *
       */
      DoFInfo &
      operator=(DoFInfo &&) noexcept = default;

      /**
       * 清除该类中的所有数据字段。
       *
       */
      void
      clear();

      /**
       * 返回一个给定的有限元度的FE索引。如果不是在hp-模式下，这个函数总是返回索引0。如果在hp-模式下没有找到索引，则返回
       * numbers::invalid_unsigned_int.  。
       *
       */
      unsigned int
      fe_index_from_degree(const unsigned int first_selected_component,
                           const unsigned int fe_degree) const;

      /**
       用存储在单元块 @p cell. 上的本地拥有的自由度填充向量 @p locall_indices  如果 @p with_constraints 是`true`，那么返回的向量将包含解决约束条件所需的索引。            下面的图片说明了这个函数对单元块0和1的输出，域的底部是零Dirichlet边界条件。请注意，由于约束条件的存在，在 "with_constraints = true "的情况下，该函数返回的DoFs不是单元块上每个单元DoFs的简单联合  @p cell.   @image html dofinfo_get_dof_indices.png  。
       * @note  返回的索引可能包含重复的。使用 `std::sort()`
       * 和 `std::unique()` 可以获得唯一的集合。
       *
       */
      void
      get_dof_indices_on_cell_batch(std::vector<unsigned int> &locall_indices,
                                    const unsigned int         cell,
                                    const bool with_constraints = true) const;

      /**
       * 这个内部方法获取单元格上的局部索引，并将其填入该类。它解决了约束条件并分配了结果。鬼魂索引，即位于另一个处理器上的索引，通过这个函数得到一个临时数字，以后在所有鬼魂索引被调用
       * @p assign_ghosts. 收集后，将被分配到正确的索引。
       *
       */
      template <typename number>
      void
      read_dof_indices(
        const std::vector<types::global_dof_index> &local_indices,
        const std::vector<unsigned int> &           lexicographic_inv,
        const dealii::AffineConstraints<number> &   constraints,
        const unsigned int                          cell_number,
        ConstraintValues<double> &                  constraint_values,
        bool &                                      cell_at_boundary);

      /**
       * 这种方法从 @p read_dof_indices
       * 函数采用的临时编号中为鬼魂指数分配了正确的指数。这些数字相对于MPI进程来说是本地化的，而鬼魂则从本地拥有的范围的末端开始。这样，我们就可以直接访问所有的向量条目。
       *
       */
      void
      assign_ghosts(const std::vector<unsigned int> &boundary_cells,
                    const MPI_Comm &                 communicator_sm,
                    const bool use_vector_data_exchanger_full);

      /**
       * 这个方法根据给定的单元格的重新编号，重新排列单元格的走法。它还将
       * @p vectorization_length
       * 单元格放在一起，并将它们解释为只有一个单元格，这是矢量化的需要。
       *
       */
      void
      reorder_cells(const TaskInfo &                  task_info,
                    const std::vector<unsigned int> & renumbering,
                    const std::vector<unsigned int> & constraint_pool_row_index,
                    const std::vector<unsigned char> &irregular_cells);

      /**
       * 查找可能的单元格索引压缩，我们可以应用它来提高效率。在reorder_cells的最后运行。
       *
       */
      void
      compute_cell_index_compression(
        const std::vector<unsigned char> &irregular_cells);

      /**
       * 查找可能的压缩面的指数，我们可以应用这些指数来提高效率。在reorder_cells的结尾处运行。
       *
       */
      template <int length>
      void
      compute_face_index_compression(
        const std::vector<FaceToCellTopology<length>> &faces);

      /**
       * 这个函数计算当前存储的指数在各个单元之间的连接性，并将结构填充为稀疏模式。
       *
       */
      void
      make_connectivity_graph(const TaskInfo &                 task_info,
                              const std::vector<unsigned int> &renumbering,
                              DynamicSparsityPattern &connectivity) const;

      /**
       * 在启用面积分的情况下，找出某些对未知数的循环是否只访问我们在主分区器中保存的所有鬼魂道夫的子集。
       *
       */
      void
      compute_tight_partitioners(
        const Table<2, ShapeInfo<double>> &       shape_info,
        const unsigned int                        n_owned_cells,
        const unsigned int                        n_lanes,
        const std::vector<FaceToCellTopology<1>> &inner_faces,
        const std::vector<FaceToCellTopology<1>> &ghosted_faces,
        const bool                                fill_cell_centric,
        const MPI_Comm &                          communicator_sm,
        const bool use_vector_data_exchanger_full);

      /**
       * 给出 @p cell_indices_contiguous_sm
       * 包含宏面（内/外）和宏面的单元格的局部索引，计算dof_indices_contiguous_sm。
       *
       */
      void
      compute_shared_memory_contiguous_indices(
        std::array<std::vector<std::pair<unsigned int, unsigned int>>, 3>
          &cell_indices_contiguous_sm);

      /**
       * 计算自由度的重新编号，以改善该类的数据访问模式，可以被IndexStorageVariants枚举中的类别所利用。例如，对于典型的DG元素，可以通过将成批的单元格的自由度交错排列来改善索引排序，这就避免了
       * IndexStorageVariants::contiguous. 中明确的数据转置。
       * 目前，这些更高级的功能没有实现，所以这个函数的价值有限。
       *
       */
      void
      compute_dof_renumbering(
        std::vector<types::global_dof_index> &renumbering);

      /**
       * 填充数组，该数组定义了如何在单元格循环内将结果向量中的选定范围归零，填充两个成员变量
       * @p  vector_zero_range_list_index和 @p vector_zero_range_list.
       * 这个模式的意图是将向量条目在时间上与第一次访问相近的地方归零，从而将向量条目保留在缓存中。
       *
       */
      template <int length>
      void
      compute_vector_zero_access_pattern(
        const TaskInfo &                               task_info,
        const std::vector<FaceToCellTopology<length>> &faces);

      /**
       * 返回该类的内存消耗，以字节为单位。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 在给定的输出流中打印出该类不同结构的内存消耗的详细摘要。
       *
       */
      template <typename StreamType>
      void
      print_memory_consumption(StreamType &    out,
                               const TaskInfo &size_info) const;

      /**
       * 在给定的输出流中打印该类中的索引表示。
       *
       */
      template <typename Number>
      void
      print(const std::vector<Number> &      constraint_pool_data,
            const std::vector<unsigned int> &constraint_pool_row_index,
            std::ostream &                   out) const;

      /**
       * 用于索引的各种存储变体的枚举。这种存储格式用于在底层数据结构允许的情况下实现更有效的索引方案，并告知
       * FEEvaluationBase::read_write_operation()
       * 中的访问函数要从哪个数组中获取数据。更有效存储的一个例子是枚举值
       * IndexStorageVariants::contiguous,
       * ，这意味着只需读取每个单元的第一个索引就可以获得单元的所有自由度的索引，而所有后续的索引仅仅是第一个索引的一个偏移。
       *
       */
      enum class IndexStorageVariants : unsigned char
      {
        /**
         * 这个值表示没有找到索引压缩，唯一有效的存储是访问单元上存在的所有索引，可能包括约束。对于这种索引类型的单元格/面，FEEvaluationBase中的数据访问被引向数组
         * @p
         * dof_indices，索引为`row_starts[cell_index*n_vectorization*n_components].first`。
         *
         */
        full,
        /**
         * 这个值表示指数是交错的，以便用矢量的聚集和分散操作来访问。这种存储变体在单元格没有约束条件，且该批单元格中的指数没有指向矢量化数组不同槽位中的同一个全局指数的情况下是可能的（为了支持散点操作）。对于这种索引类型的单元格/面，FEEvaluationBase中的数据访问是指向索引为
         * "row_starts[cell_index*n_vectorization*n_components].first
         * "的阵列`dof_indices_interleaved`。
         *
         */
        interleaved,
        /**
         * 这个值表示一个单元格内的索引都是连续的，人们可以通过读取该单元格批中每个单元格的这个单一值来获得该单元格的索引。对于这种索引类型的单元格/面，FEEvaluationBase中的数据访问是指向索引为
         * "cell_index*n_vectorization*n_components
         * "的数组`dof_indices_contiguous`。
         *
         */
        contiguous,
        /**
         * 这个值表示与一个单元的索引是连续的，并交错进行矢量化，即一个单元上的第一个DoF索引到矢量化批次中的四个或八个单元先到，比第二个DoF索引，以此类推。此外，单元格之间的交错意味着只有用于矢量化的批次可以被有效地访问，而对于只获得部分条目的访问则是一种分层访问。
         * 两个额外的类别`interleaved_contiguous_strided`和`interleaved_contiguous_mixed_strides`是这种存储类型的结果。前者适用于相邻的两个面中至少有一个面会因交错存储而断开。那么我们就必须按照下一个类别的描述进行串联访问。最后一个类别`interleaved_contiguous_mixed_strides`出现在ghost层，见下面对该类别的详细描述。
         * 同样，一旦我们在单元格之间交错索引，这也是一般情况下无法避免的事情。
         * 对于这种索引类型的单元格/面，FEEvaluationBase中的数据访问是指向索引为`cell_index*n_vectorization*n_components`的`dof_indices_contiguous`数组。
         *
         */
        interleaved_contiguous,
        /**
         * 与Interleaved_contiguous存储类似，但适用于Interleaved指数只在自由度内连续，而不是在矢量数组的组件上连续的情况。
         * 这种情况通常发生在有DG的面，其中的单元有`交错_连续`存储，但面的编号与单元的编号不一样。对于这种索引类型的单元格/面，FEEvaluationBase中的数据访问被引导到索引为`cell_index*n_vectorization*n_components'的`dof_indices_contiguous'阵列。
         *
         */
        interleaved_contiguous_strided,
        /**
         * 类似于interleaved_contiguous_separate存储，但针对的是交错索引不相隔`n_vectorization`的情况。这种情况通常发生在DG的幽灵层中，远程所有者应用了交错存储，而当前的处理器只看到其中的一些单元。对于这种索引类型的单元/面，FEEvaluationBase中的数据访问被引向索引为`cell_index*n_vectorization*n_components'的`dof_indices_contiguous'数组，包括实际跨度信息的`dof_indices_interleave_strides`数组。
         *
         */
        interleaved_contiguous_mixed_strides
      };

      /**
       * 用于区分单元格和面的矢量化类型的数据阵列的枚举。
       *
       */
      enum DoFAccessIndex : unsigned char
      {
        /**
         * 指定为内部的面的数据索引
         *
         */
        dof_access_face_interior = 0,
        /**
         * 指定为外部的面的数据索引
         *
         */
        dof_access_face_exterior = 1,
        /**
         * 单元的数据索引
         *
         */
        dof_access_cell = 2
      };

      /**
       * 存储底层DoFHandler的尺寸。由于索引不是模板化的，这是一个变量，可以在这个类中需要的（很少）情况下访问尺寸。
       *
       */
      unsigned int dimension;

      /**
       * 出于效率的考虑，总是将具有类似属性的单元格的固定数量放在一起。这个变量控制着被分在一起的单元格的数量。相对于其他类在数字类型上的模板化，这个类作为一个纯粹的索引容器是没有模板化的，所以我们需要保留否则包含在
       * VectorizedArrayType::size(). 的信息。
       *
       */
      unsigned int vectorization_length;

      /**
       * 存储所有单元和面批的索引存储变体。
       * 这里给出的三个数组根据CellOrFaceAccess解决作为内部装饰的面（0）、作为外部装饰的面（1）和单元（2）的类型。
       *
       */
      std::array<std::vector<IndexStorageVariants>, 3> index_storage_variants;

      /**
       * 在 @p  dof_indices和 @p constraint_indicator
       * 字段中存储压缩行存储的行起始指数。这两个字段总是被一起访问，所以只为它们保留一个变量会更简单。这也避免了保持两个行开始向量的同步。
       *
       */
      std::vector<std::pair<unsigned int, unsigned int>> row_starts;

      /**
       * 存储每个单元的自由度指数。这些指数是在MPI本地索引空间中计算的，即每个处理器将本地拥有的指数存储为<tt>0</tt>和<tt>n_locally_owned_dofs-1</tt>之间的数字，以及<tt>n_locally_owned_dofs</tt>到<tt>n_locally_owned_dofs+n_ghost_dofs</tt>范围内的幽灵指数。这个MPI本地索引空间和全局自由度编号之间的转换被存储在
       * @p vector_partitioner 数据结构中。
       * 这个数组还包括来自约束的间接贡献，由 @p
       * constraint_indicator
       * 字段描述。由于行的长度可变，这将是一个矢量的矢量。
       * 然而，我们使用一个连续的内存区域，并将行开始存储在变量
       * @p row_starts. 中
       *
       */
      std::vector<unsigned int> dof_indices;

      /**
       * 这个变量以单元上自由度的局部编号来描述约束的位置。第一个数字存储了从一个约束自由度到下一个约束自由度的距离。当我们从矢量中读出或写入单元的局部自由度时，这可以确定受约束自由度的位置。第二个数字存储的是约束权重的索引，存储在变量constraint_pool_data中。
       *
       */
      std::vector<std::pair<unsigned short, unsigned short>>
        constraint_indicator;

      /**
       * 为 `IndexStorageVariants::interleaved`.
       * 重新排序的索引存储。
       *
       */
      std::vector<unsigned int> dof_indices_interleaved;

      /**
       * 压缩的索引存储，比通过 @p
       * dof_indices使用的索引存储更快，根据IndexStorageVariants中的描述。
       * 这里给出的三个数组根据CellOrFaceAccess处理内部装饰的面（0）、外部装饰的面（1）和单元格（2）的类型。
       *
       */
      std::array<std::vector<unsigned int>, 3> dof_indices_contiguous;

      /**
       * 与上述相同，但用于共享内存的使用。对的第一个值是识别拥有的进程，第二个是该进程的本地拥有数据中的索引。
       * @note
       * 这个数据结构只有在index_storage_variants[2]中的所有条目都是
       * IndexStorageVariants::contiguous. 时才会设置。
       *
       */
      std::array<std::vector<std::pair<unsigned int, unsigned int>>, 3>
        dof_indices_contiguous_sm;

      /**
       * 压缩索引存储，比通过根据IndexStorageVariants中的描述使用的
       * @p  dof_indices更快地访问。
       * 这里给出的三个数组解决了装饰为减号的面（0）、装饰为加号的面（1）和细胞（2）的类型。
       *
       */
      std::array<std::vector<unsigned int>, 3> dof_indices_interleave_strides;

      /**
       * 缓存矢量化时填充的索引数。这个信息可以隐含地从row_starts数据字段中推导出来，但是这个字段允许更快的访问。
       * 这里给出的三个数组根据CellOrFaceAccess来处理内部装饰的面（0）、外部装饰的面（1）和单元格（2）的类型。
       *
       */
      std::array<std::vector<unsigned char>, 3> n_vectorization_lanes_filled;

      /**
       * 这存储了可用于设置向量的平行分区。分区器包括对向量中局部范围的描述，也包括鬼魂的样子。这使得基于DoFInfo字段的向量的初始化成为可能。
       *
       */
      std::shared_ptr<const Utilities::MPI::Partitioner> vector_partitioner;

      /**
       * 与vector_partitioner兼容的向量交换器。
       *
       */
      std::shared_ptr<
        const internal::MatrixFreeFunctions::VectorDataExchange::Base>
        vector_exchanger;

      /**
       * 与分区器兼容的向量交换器，该分区器选择了存储在
       * @p vector_partitioner.
       * 中的完整向量分区器的鬼魂索引子集。这些分区器用于专门的循环，只导入鬼魂区域的一部分，以减少通信量。有五种变体的分区器被初始化。
       *
       *
       *
       *
       *
       * - 一个是只查询单元格的值。
       *
       * - 一个额外描述在相关面上评估函数值的指数。
       *
       *
       *
       *
       *
       *
       * - 一个描述用于评估函数值和邻近本地所有单元的相关面上的梯度的指数。
       *
       *
       *
       * - 一个额外描述在所有面上评估函数值的索引，以及
       *
       *
       *
       *
       *
       *
       * - 一个描述用于评估函数值和梯度的指数，在与本地所有单元相邻的所有面上。
       *
       */
      std::array<
        std::shared_ptr<
          const internal::MatrixFreeFunctions::VectorDataExchange::Base>,
        5>
        vector_exchanger_face_variants;

      /**
       * 这存储了所有本地拥有的自由度的（排序）列表，这些自由度被约束。
       *
       */
      std::vector<unsigned int> constrained_dofs;

      /**
       * 在 @p
       * plain_dof_indices字段中存储压缩行存储的行开始指数。
       *
       */
      std::vector<unsigned int> row_starts_plain_indices;

      /**
       * 存储每个单元的自由度指数。这个数组不包括约束条件的间接贡献，这些贡献包括在
       * @p dof_indices. 中。
       * 由于行的长度可变，这将是一个矢量的矢量。然而，我们使用一个连续的内存区域，并将行开始存储在变量
       * @p  row_starts_plain_indices中。
       *
       */
      std::vector<unsigned int> plain_dof_indices;

      /**
       * 以所有DoFInfo对象的基数元素的数量来存储偏移量。
       *
       */
      unsigned int global_base_element_offset;

      /**
       * 存储从DoFHandler中读出指数的基本元素的数量。
       *
       */
      unsigned int n_base_elements;

      /**
       * 储存有限元中每个基元的分量数，索引是从那里读来的。
       *
       */
      std::vector<unsigned int> n_components;

      /**
       * 这个向量的第1个条目存储了给定基元的分量号。
       *
       */
      std::vector<unsigned int> start_components;

      /**
       * 对于FES系统中的一个给定的分量，这个变量告诉人们该索引属于哪个基元。
       *
       */
      std::vector<unsigned int> component_to_base_index;

      /**
       * 对于一个矢量值元素，这给出了从给定组件开始的自由度数量中的常数偏移，因为自由度是按自由度编号的。这个数据结构没有考虑到可能的约束，因此，更短或更长的列表。这一信息直接编码在row_starts变量中。
       * 外围向量在hp情况下经过各种FE指数，与 @p dofs_per_cell
       * 变量类似。
       *
       */
      std::vector<std::vector<unsigned int>> component_dof_indices_offset;

      /**
       * 存储每个单元的自由度数。
       *
       */
      std::vector<unsigned int> dofs_per_cell;

      /**
       * 存储每个面的自由度数。
       *
       */
      std::vector<unsigned int> dofs_per_face;

      /**
       * 告知平原指数是否被缓存。
       *
       */
      bool store_plain_indices;

      /**
       * 存储hp情况下的活动有限元的索引。
       *
       */
      std::vector<unsigned int> cell_active_fe_index;

      /**
       * 存储hp-情况下不同有限元的最大度数。
       *
       */
      unsigned int max_fe_index;

      /**
       * 对hp-adaptive
       * case中的每个槽，内向量存储相应的元素度。这被FEEvaluationBase的构造器用来识别hp-案例中的正确数据槽。
       *
       */
      std::vector<std::vector<unsigned int>> fe_index_conversion;

      /**
       * 在设置过程中临时存储鬼魂的数量。在调用 @p
       * assign_ghosts.
       * 时被清除。然后，所有的信息都由分区器收集。
       *
       */
      std::vector<types::global_dof_index> ghost_dofs;

      /**
       * 为TaskInfo中的每个分区存储一个整数，表明如果用户在
       * MatrixFree::loop.
       * 中用相应的参数要求，是否清除结果向量中的某些部分。
       *
       */
      std::vector<unsigned int> vector_zero_range_list_index;

      /**
       * 存储要清除的向量中的实际范围。
       *
       */
      std::vector<std::pair<unsigned int, unsigned int>> vector_zero_range_list;

      /**
       * 为TaskInfo中的每个分区存储一个整数，表明何时安排操作，这些操作将在对向量项的任何访问之前完成。
       *
       */
      std::vector<unsigned int> cell_loop_pre_list_index;

      /**
       * 在对向量项进行任何访问之前，存储操作的实际范围。
       *
       */
      std::vector<std::pair<unsigned int, unsigned int>> cell_loop_pre_list;

      /**
       * 为TaskInfo中的每个分区存储一个整数，表明何时安排将在所有访问向量项之后进行的操作。
       *
       */
      std::vector<unsigned int> cell_loop_post_list_index;

      /**
       * 存储在所有访问向量项之后的操作的实际范围。
       *
       */
      std::vector<std::pair<unsigned int, unsigned int>> cell_loop_post_list;
    };


     /*-------------------------- Inline functions ---------------------------*/ 

#ifndef DOXYGEN


    inline unsigned int
    DoFInfo::fe_index_from_degree(const unsigned int first_selected_component,
                                  const unsigned int fe_degree) const
    {
      const unsigned int n_indices = fe_index_conversion.size();
      if (n_indices <= 1)
        return 0;
      for (unsigned int i = 0; i < n_indices; ++i)
        if (fe_index_conversion[i][first_selected_component] == fe_degree)
          return i;
      return numbers::invalid_unsigned_int;
    }

#endif // ifndef DOXYGEN

  } // end of namespace MatrixFreeFunctions
} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif


