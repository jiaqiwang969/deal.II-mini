//include/deal.II-translator/matrix_free/matrix_free_0.txt
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


#ifndef dealii_matrix_free_h
#define dealii_matrix_free_h

#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/thread_local_storage.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_operation.h>

#include <deal.II/matrix_free/dof_info.h>
#include <deal.II/matrix_free/mapping_info.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/task_info.h>
#include <deal.II/matrix_free/type_traits.h>

#include <cstdlib>
#include <limits>
#include <list>
#include <memory>


DEAL_II_NAMESPACE_OPEN



/**
 * 这个类收集了所有为矩阵自由实现而存储的数据。这个存储方案是针对用相同数据进行的几个循环而定制的，也就是说，通常在同一个网格上做许多矩阵-向量乘积或残差计算。该类在
 * step-37  和  step-48  中使用。
 * 该类不实现任何涉及有限元基函数的操作，即关于对单元的操作。对于这些操作，FEEvaluation类被设计为使用该类中收集的数据。
 * 存储的数据可以细分为三个主要部分。
 *
 *
 *
 * - DoFInfo：它存储了局部自由度与全局自由度的关系。它包括对约束条件的描述，这些约束条件被评估为穿过一个单元上的所有局部自由度。
 *
 *
 *
 * - MappingInfo：它存储了从实数到单位单元的转换，这些转换对于建立有限元函数的导数和找到物理空间中正交权重的位置是必要的。
 *
 *
 * - ShapeInfo：它包含有限元的形状函数，在单元格上进行评估。
 * 除了初始化程序，这个类只实现了一个单一的操作，即所有单元的循环（cell_loop()）。这个循环的安排方式是，共享自由度的单元不会同时工作，这意味着可以并行地写入向量（或矩阵），而不必明确地同步访问这些向量和矩阵。这个类没有实现任何形状值，它所做的只是缓存相应的数据。要实现有限元操作，请使用FEEvaluation类（或一些相关的类）。
 * 这个类以不同的顺序遍历单元，与deal.II中通常的Triangulation类不同，以便在共享内存和矢量化的并行化方面具有灵活性。
 * 矢量化是通过将几个拓扑单元合并成一个所谓的宏观单元来实现的。这使得几个单元的所有相关操作都可以用一条CPU指令来实现，这也是这个框架的主要特点之一。
 * 关于这个类的使用细节，见FEEvaluation的描述或 @ref matrixfree "无矩阵模块"
 * 。
 *
 *
 * @ingroup matrixfree
 *
 */

template <int dim,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
class MatrixFree : public Subscriptor
{
  static_assert(
    std::is_same<Number, typename VectorizedArrayType::value_type>::value,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  /**
   * 由模板参数指定的底层数字类型的别名。
   *
   */
  using value_type            = Number;
  using vectorized_value_type = VectorizedArrayType;

  /**
   * 由模板参数`dim`设定的尺寸。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 收集用于MatrixFree类初始化的选项。第一个参数指定要使用的MPI通信器，第二个参数是共享内存中的并行化选项（基于任务的并行化，可以选择无并行化和避免访问相同向量项的单元同时被访问的三种方案），第三个参数是任务并行调度的块大小，第四个参数是这个类应该存储的更新标志。
   * 第五个参数指定了三角形中的级别，索引将从该级别开始使用。如果级别设置为
   * numbers::invalid_unsigned_int,
   * ，将遍历活动单元，否则将遍历指定级别的单元。这个选项在给出DoFHandler的情况下没有影响。
   * 参数 @p initialize_plain_indices
   * 表示DoFInfo类是否也应该允许访问向量而不解决约束。
   * 两个参数`initialize_indices`和`initialize_mapping`允许用户停用一些初始化过程。例如，如果只需要找到避免同时触及相同向量/矩阵指数的调度，就不需要初始化映射。同样，如果映射从一个迭代到下一个迭代发生了变化，但拓扑结构没有变化（比如使用MappingQEulerian的变形网格时），只需要初始化映射就足够了。
   * 两个参数`cell_vectorization_categories`和`cell_vectorization_categories_strict`控制在几个单元上形成向量化的批次。它在使用hp-adaptivity时被隐含使用，但在其他情况下也很有用，比如在本地时间步进中，人们想控制哪些元素一起形成一批单元。数组`cell_vectorization_categories`在`mg_level`设置为
   * `numbers::invalid_unsigned_int`
   * 的情况下，通过cell->active_cell_index()访问活动单元，对于水平单元则通过cell->index()访问。默认情况下，`cell_vectorization_category`中的不同类别可以混合使用，如果算法内部需要，允许算法将较低的类别数字与下一个较高的类别合并，以尽可能避免部分填充的SIMD通道。这样可以更好地利用矢量，但可能需要特殊处理，特别是对面积分。如果设置为
   * @p true,
   * ，算法反而会将不同的类别分开，而不是将它们混合在一个矢量化数组中。
   *
   */
  struct AdditionalData
  {
    /**
     * 收集任务并行的选项。参见成员变量
     * MatrixFree::AdditionalData::tasks_parallel_scheme
     * 的文档以获得详尽的描述。
     *
     */
    enum TasksParallelScheme
    {
      /**
       * 以串行方式执行应用。
       *
       */
      none = internal::MatrixFreeFunctions::TaskInfo::none,
      /**
       * 将单元格分成两层，之后形成块状。
       *
       */
      partition_partition =
        internal::MatrixFreeFunctions::TaskInfo::partition_partition,
      /**
       * 在全局层面进行分区，并对分区内的单元格进行着色。
       *
       */
      partition_color =
        internal::MatrixFreeFunctions::TaskInfo::partition_color,
      /**
       * 使用传统的着色算法：这就像
       * TasksParallelScheme::partition_color, ，但只使用一个分区。
       *
       */
      color = internal::MatrixFreeFunctions::TaskInfo::color
    };

    /**
     * 附加数据的构造函数。
     *
     */
    AdditionalData(
      const TasksParallelScheme tasks_parallel_scheme = partition_partition,
      const unsigned int        tasks_block_size      = 0,
      const UpdateFlags         mapping_update_flags  = update_gradients |
                                               update_JxW_values,
      const UpdateFlags  mapping_update_flags_boundary_faces = update_default,
      const UpdateFlags  mapping_update_flags_inner_faces    = update_default,
      const UpdateFlags  mapping_update_flags_faces_by_cells = update_default,
      const unsigned int mg_level            = numbers::invalid_unsigned_int,
      const bool         store_plain_indices = true,
      const bool         initialize_indices  = true,
      const bool         initialize_mapping  = true,
      const bool         overlap_communication_computation    = true,
      const bool         hold_all_faces_to_owned_cells        = false,
      const bool         cell_vectorization_categories_strict = false)
      : tasks_parallel_scheme(tasks_parallel_scheme)
      , tasks_block_size(tasks_block_size)
      , mapping_update_flags(mapping_update_flags)
      , mapping_update_flags_boundary_faces(mapping_update_flags_boundary_faces)
      , mapping_update_flags_inner_faces(mapping_update_flags_inner_faces)
      , mapping_update_flags_faces_by_cells(mapping_update_flags_faces_by_cells)
      , mg_level(mg_level)
      , store_plain_indices(store_plain_indices)
      , initialize_indices(initialize_indices)
      , initialize_mapping(initialize_mapping)
      , overlap_communication_computation(overlap_communication_computation)
      , hold_all_faces_to_owned_cells(hold_all_faces_to_owned_cells)
      , cell_vectorization_categories_strict(
          cell_vectorization_categories_strict)
      , communicator_sm(MPI_COMM_SELF)
    {}

    /**
     * 复制构造器。
     *
     */
    AdditionalData(const AdditionalData &other)
      : tasks_parallel_scheme(other.tasks_parallel_scheme)
      , tasks_block_size(other.tasks_block_size)
      , mapping_update_flags(other.mapping_update_flags)
      , mapping_update_flags_boundary_faces(
          other.mapping_update_flags_boundary_faces)
      , mapping_update_flags_inner_faces(other.mapping_update_flags_inner_faces)
      , mapping_update_flags_faces_by_cells(
          other.mapping_update_flags_faces_by_cells)
      , mg_level(other.mg_level)
      , store_plain_indices(other.store_plain_indices)
      , initialize_indices(other.initialize_indices)
      , initialize_mapping(other.initialize_mapping)
      , overlap_communication_computation(
          other.overlap_communication_computation)
      , hold_all_faces_to_owned_cells(other.hold_all_faces_to_owned_cells)
      , cell_vectorization_category(other.cell_vectorization_category)
      , cell_vectorization_categories_strict(
          other.cell_vectorization_categories_strict)
      , communicator_sm(other.communicator_sm)
    {}

    /**
     * 拷贝赋值。
     *
     */
    AdditionalData &
    operator=(const AdditionalData &other)
    {
      tasks_parallel_scheme = other.tasks_parallel_scheme;
      tasks_block_size      = other.tasks_block_size;
      mapping_update_flags  = other.mapping_update_flags;
      mapping_update_flags_boundary_faces =
        other.mapping_update_flags_boundary_faces;
      mapping_update_flags_inner_faces = other.mapping_update_flags_inner_faces;
      mapping_update_flags_faces_by_cells =
        other.mapping_update_flags_faces_by_cells;
      mg_level            = other.mg_level;
      store_plain_indices = other.store_plain_indices;
      initialize_indices  = other.initialize_indices;
      initialize_mapping  = other.initialize_mapping;
      overlap_communication_computation =
        other.overlap_communication_computation;
      hold_all_faces_to_owned_cells = other.hold_all_faces_to_owned_cells;
      cell_vectorization_category   = other.cell_vectorization_category;
      cell_vectorization_categories_strict =
        other.cell_vectorization_categories_strict;
      communicator_sm = other.communicator_sm;

      return *this;
    }

    /**
     * 设置任务并行的方案。有四个选项可用。
     * 如果设置为 @p none,
     * ，运算器的应用是以串行方式进行的，没有共享内存并行。如果这个类和MPI一起使用，并且MPI也用于节点内的并行，那么这个标志应该设置为
     * @p none. 默认值是 @p partition_partition,
     * ，即我们实际上用下面描述的第一个选项使用多线程。
     * 第一个选项 @p partition_partition
     * 是在洋葱皮一样的分区中对单元格进行两级分割，并在分割后形成tasks_block_size的块。分区找到独立的单元组，使其能够在不同时访问相同的向量项的情况下并行工作。
     * 第二种选择 @p partition_color
     * 是在全局层面上使用分区，并在分区内为单元格着色（在一个颜色内的所有块是独立的）。在这里，将单元格细分为若干块是在分区之前完成的，如果自由度没有被重新编号，可能会得到更差的分区，但缓存性能更好。
     * 第三种方案 @p color
     * 是在全局层面使用传统的着色算法。这种方案是第二种方案的一个特例，即只存在一个分区。请注意，对于有悬空节点的问题，有相当多的颜色（3D中为50种或更多），这可能会降低并行性能（不良的缓存行为，许多同步点）。
     * @note
     * 对于进行内面积分的情况，线程支持目前是实验性的，如果可能的话，建议使用MPI并行。虽然该方案已经验证了在通常的DG元素的情况下可以使用`partition_partition`选项，但对于更一般的元素系统，如连续和不连续元素的组合，对所有项添加面积分的情况，没有进行全面的测试。
     *
     */
    TasksParallelScheme tasks_parallel_scheme;

    /**
     * 设置应该形成一个分区的所谓宏单元的数量。如果给定的尺寸为零，该类会尝试根据
     * MultithreadInfo::n_threads()
     * 和存在的单元数为块找到一个好的尺寸。否则，将使用给定的数字。如果给定的数字大于总单元格数的三分之一，这意味着没有并行性。请注意，在使用矢量化的情况下，一个宏单元由一个以上的物理单元组成。
     *
     */
    unsigned int tasks_block_size;

    /**
     * 这个标志决定了单元格上的映射数据被缓存起来。该类可以缓存梯度计算（反雅各布）、雅各布行列式（JxW）、正交点以及Hessians（雅各布导数）的数据。默认情况下，只有梯度和雅各布行列式乘以正交权重JxW的数据被缓存。如果需要正交点或二次导数，它们必须由这个字段指定（即使这里没有设置这个选项，二次导数仍然可以在笛卡尔单元上评估，因为那里的雅各布系数完全描述了映射）。
     *
     */
    UpdateFlags mapping_update_flags;

    /**
     * 这个标志决定了边界面上的映射数据要被缓存。注意MatrixFree对面的积分使用单独的循环布局，以便在悬空节点（需要在两边有不同的子面设置）或者在一个VectorizedArray的批次中的一些单元与边界相邻而其他单元不相邻的情况下也能有效地进行矢量。
     * 如果设置为不同于update_general的值（默认），面的信息将被显式建立。目前，MatrixFree支持缓存以下面的数据：逆雅各布，雅各布行列式（JxW），正交点，Hessians数据（雅各布的导数），以及法向量。
     * @note 为了能够在 MatrixFree::loop()`, 中执行 "face_operation
     * "或 "boundary_operation"，这个字段或 @p
     * mapping_update_flags_inner_faces 必须被设置为与
     * UpdateFlags::update_default. 不同的值。
     *
     */
    UpdateFlags mapping_update_flags_boundary_faces;

    /**
     * 这个标志决定了要缓存的内部面的映射数据。请注意，MatrixFree对面的积分使用单独的循环布局，以便在悬空节点（需要在两边有不同的子面设置）或者在一个VectorizedArray的批中的一些单元与边界相邻而其他单元不相邻的情况下也能有效地进行矢量化。
     * 如果设置为不同于update_general的值（默认），面的信息将被显式建立。目前，MatrixFree支持缓存以下面的数据：逆雅各布，雅各布行列式（JxW），正交点，Hessians数据（雅各布的导数），以及法向量。
     * @note 为了能够在 MatrixFree::loop()`, 中执行 "face_operation
     * "或 "boundary_operation"，这个字段或 @p
     * mapping_update_flags_boundary_faces 必须被设置为与
     * UpdateFlags::update_default. 不同的值。
     *
     */
    UpdateFlags mapping_update_flags_inner_faces;

    /**
     * 这个标志决定了不同布局的面的映射数据与矢量化的关系。当`mapping_update_flags_inner_faces`和`mapping_update_flags_boundary_faces`触发以面为中心的方式建立数据并进行适当的矢量化时，当前的数据域将面的信息附加到单元格及其矢量化的方式。这只在特殊情况下需要，例如在块状Jacobi方法中，单元格的全部算子包括其面都被评估。该数据由
     * <code>FEFaceEvaluation::reinit(cell_batch_index,
     * face_number)</code>访问。然而，目前用这种方法不能计算到邻居的耦合项，因为邻居不是由VectorizedArray数据布局的数组-结构-数组类型的数据结构来布置的。
     * 注意，你应该只在真正需要的情况下计算这个数据域，因为它使面的映射数据所需的内存增加了一倍以上。
     * 如果设置为不同于update_general的值（默认值），就会显式地建立面的信息。目前，MatrixFree支持缓存以下面的数据：反Jacobian，Jacobian行列式（JxW），正交点，Hessians（Jacobian的导数）数据，以及法向量。
     *
     */
    UpdateFlags mapping_update_flags_faces_by_cells;

    /**
     * 这个选项可以用来定义我们是否在网格的某一层工作，而不是活动单元。如果设置为invalid_unsigned_int
     * (这是默认值)，就会对活动单元进行处理，否则就按这个参数给定的级别处理。请注意，如果你指定在一个层次上工作，其道夫必须通过使用
     * <code>dof_handler.distribute_mg_dofs(fe);</code> 来分配。
     *
     */
    unsigned int mg_level;

    /**
     * 控制是否启用从向量读取而不解决约束，即只读取向量的局部值。默认情况下，这个选项是启用的。如果你想使用
     * FEEvaluationBase::read_dof_values_plain,
     * ，这个标志需要被设置。
     *
     */
    bool store_plain_indices;

    /**
     * 选项用于控制是否应该读取存储在DoFHandler中的索引，以及是否应该在MatrixFree的初始化方法中设置任务并行的模式。默认值是true。如果映射需要重新计算（例如，当使用通过MappingEulerian描述的变形网格时），但单元格的拓扑结构保持不变，可以禁用。
     *
     */
    bool initialize_indices;

    /**
     * 控制是否应该在MatrixFree的初始化方法中计算映射信息的选项。默认值为true。当只需要设置一些索引时可以禁用（例如，当只需要计算一组独立的单元格时）。
     *
     */
    bool initialize_mapping;

    /**
     * 如果传递给循环的向量支持非阻塞数据交换，可以选择控制循环是否应该尽可能地重叠通信和计算。在大多数情况下，如果要发送的数据量超过几千字节，重叠会更快。如果发送的数据较少，在大多数集群上，通信是有延迟的（按照2016年的标准，在好的集群上，点到点的延迟大约是1微秒）。根据MPI实现和结构的不同，不重叠并等待数据的到来可能会更快。默认情况下是true，即通信和计算是重叠的。
     *
     */
    bool overlap_communication_computation;

    /**
     * 默认情况下，面的部分将只保存那些要被本地处理的面（和面后面的鬼元素）。如果MatrixFree需要访问本地拥有的单元上的所有邻域，这个选项可以在面域的末端添加相应的面。
     *
     */
    bool hold_all_faces_to_owned_cells;

    /**
     * 这个数据结构允许在建立矢量化信息时将一部分单元分配给不同的类别。在使用hp-adaptivity时，它被隐含地使用，但在其他情况下也是有用的，比如在本地时间步进中，人们想控制哪些元素共同构成一批单元。
     * 在 @p mg_level 设置为 numbers::invalid_unsigned_int
     * 的活动单元上工作时，通过cell->active_cell_index()给出的数字访问这个数组，对于水平单元则通过cell->index()访问。
     * @note
     * 这个字段在构建AdditionalData时是空的。用户有责任在填充数据时将此字段调整为`triangulation.n_active_cells()`或`triangulation.n_cells(level)`。
     *
     */
    std::vector<unsigned int> cell_vectorization_category;

    /**
     * 默认情况下， @p cell_vectorization_category
     * 中的不同类别可以混合，如果在算法内部有必要，允许算法将较低的类别与次高的类别合并。这样可以更好地利用矢量，但可能需要特殊处理，特别是对于面积分。如果设置为
     * @p true,
     * ，算法会将不同的类别分开，而不是将它们混合在一个矢量化的数组中。
     *
     */
    bool cell_vectorization_categories_strict;

    /**
     * 共享内存的MPI通信器。默认值。MPI_COMM_SELF。
     *
     */
    MPI_Comm communicator_sm;
  };

  /**
   * @name  1：构造和初始化
   *
   */
  //@{
  /**
   * 默认的空构造函数。什么都不做。
   *
   */
  MatrixFree();

  /**
   * 复制构造函数，调用copy_from
   *
   */
  MatrixFree(const MatrixFree<dim, Number, VectorizedArrayType> &other);

  /**
   * 解构器。
   *
   */
  ~MatrixFree() override = default;

  /**
   * 提取在单元格上进行循环所需的信息。DoFHandler和AffineConstraints对象描述了自由度的布局，DoFHandler和映射描述了从单元到实数单元的转换，而DoFHandler底层的有限元与正交公式一起描述了局部操作。请注意，DoFHandler底层的有限元必须是标量的，或者包含同一元素的几个副本。不允许将几个不同的元素混入一个FES系统。在这种情况下，使用带有几个DoFHandler参数的初始化函数。
   *
   */
  template <typename QuadratureType, typename number2, typename MappingType>
  void
  reinit(const MappingType &               mapping,
         const DoFHandler<dim> &           dof_handler,
         const AffineConstraints<number2> &constraint,
         const QuadratureType &            quad,
         const AdditionalData &            additional_data = AdditionalData());

  /**
   * 初始化数据结构。和上面一样，但使用 $Q_1$ 的映射。
   * @deprecated  使用重载来代替Mapping对象。
   *
   */
  template <typename QuadratureType, typename number2>
  DEAL_II_DEPRECATED void
  reinit(const DoFHandler<dim> &           dof_handler,
         const AffineConstraints<number2> &constraint,
         const QuadratureType &            quad,
         const AdditionalData &            additional_data = AdditionalData());

  /**
   * 提取在单元格上进行循环所需的信息。DoFHandler和AffineConstraints对象描述了自由度的布局，DoFHandler和映射描述了从单元到实数单元的转换，DoFHandler的基础有限元与正交公式一起描述了局部操作。与其他初始化函数处理的标量情况不同，这个函数允许有两个或多个不同的有限元的问题。每个元素的DoFHandlers必须作为指针传递给初始化函数。另外，一个由多个构件组成的系统也可以由一个带有FESystem元素的单一DoFHandler来表示。这种情况的前提条件是FESystem的每个基本元素必须与本类兼容，如FE_Q或FE_DGQ类。
   * 这个函数也允许使用多个正交公式，例如当描述中包含不同程度的元素的独立积分时。然而，不同正交公式的数量可以独立于DoFHandlers的数量进行设置，当几个元素总是用同一个正交公式进行积分时。
   *
   */
  template <typename QuadratureType, typename number2, typename MappingType>
  void
  reinit(const MappingType &                                    mapping,
         const std::vector<const DoFHandler<dim> *> &           dof_handler,
         const std::vector<const AffineConstraints<number2> *> &constraint,
         const std::vector<QuadratureType> &                    quad,
         const AdditionalData &additional_data = AdditionalData());

  /**
   * 初始化数据结构。与上述相同，但使用DoFHandlerType。
   * @deprecated 使用获取DoFHandler对象的重载来代替。
   *
   */
  template <typename QuadratureType,
            typename number2,
            typename DoFHandlerType,
            typename MappingType>
  DEAL_II_DEPRECATED void
  reinit(const MappingType &                                    mapping,
         const std::vector<const DoFHandlerType *> &            dof_handler,
         const std::vector<const AffineConstraints<number2> *> &constraint,
         const std::vector<QuadratureType> &                    quad,
         const AdditionalData &additional_data = AdditionalData());

  /**
   * 初始化数据结构。和上面一样，但使用 $Q_1$ 的映射。
   * @deprecated  使用获取Mapping对象的重载来代替。
   *
   */
  template <typename QuadratureType, typename number2>
  DEAL_II_DEPRECATED void
  reinit(const std::vector<const DoFHandler<dim> *> &           dof_handler,
         const std::vector<const AffineConstraints<number2> *> &constraint,
         const std::vector<QuadratureType> &                    quad,
         const AdditionalData &additional_data = AdditionalData());

  /**
   * 初始化数据结构。和上面一样，但使用DoFHandlerType。
   * @deprecated  使用获取DoFHandler对象的重载来代替。
   *
   */
  template <typename QuadratureType, typename number2, typename DoFHandlerType>
  DEAL_II_DEPRECATED void
  reinit(const std::vector<const DoFHandlerType *> &            dof_handler,
         const std::vector<const AffineConstraints<number2> *> &constraint,
         const std::vector<QuadratureType> &                    quad,
         const AdditionalData &additional_data = AdditionalData());

  /**
   * 初始化数据结构。和以前一样，但现在本地拥有的自由度范围的索引集描述取自DoFHandler。此外，只使用一个正交公式，当一个矢量值问题中的几个分量基于同一个正交公式被整合在一起时，这可能是必要的。
   *
   */
  template <typename QuadratureType, typename number2, typename MappingType>
  void
  reinit(const MappingType &                                    mapping,
         const std::vector<const DoFHandler<dim> *> &           dof_handler,
         const std::vector<const AffineConstraints<number2> *> &constraint,
         const QuadratureType &                                 quad,
         const AdditionalData &additional_data = AdditionalData());

  /**
   * 初始化数据结构。与上述相同，但使用DoFHandlerType。
   * @deprecated 使用获取DoFHandler对象的重载来代替。
   *
   */
  template <typename QuadratureType,
            typename number2,
            typename DoFHandlerType,
            typename MappingType>
  DEAL_II_DEPRECATED void
  reinit(const MappingType &                                    mapping,
         const std::vector<const DoFHandlerType *> &            dof_handler,
         const std::vector<const AffineConstraints<number2> *> &constraint,
         const QuadratureType &                                 quad,
         const AdditionalData &additional_data = AdditionalData());

  /**
   * 初始化数据结构。和上面一样，但使用 $Q_1$ 的映射。
   * @deprecated  使用获取映射对象的重载来代替。
   *
   */
  template <typename QuadratureType, typename number2>
  DEAL_II_DEPRECATED void
  reinit(const std::vector<const DoFHandler<dim> *> &           dof_handler,
         const std::vector<const AffineConstraints<number2> *> &constraint,
         const QuadratureType &                                 quad,
         const AdditionalData &additional_data = AdditionalData());

  /**
   * 初始化数据结构。和上面一样，但使用DoFHandlerType。
   * @deprecated  使用获取DoFHandler对象的重载来代替。
   *
   */
  template <typename QuadratureType, typename number2, typename DoFHandlerType>
  DEAL_II_DEPRECATED void
  reinit(const std::vector<const DoFHandlerType *> &            dof_handler,
         const std::vector<const AffineConstraints<number2> *> &constraint,
         const QuadratureType &                                 quad,
         const AdditionalData &additional_data = AdditionalData());

  /**
   * 复制函数。创建一个所有数据结构的深度拷贝。通常情况下，为不同的操作保留一次数据就足够了，所以不应该经常需要这个函数。
   *
   */
  void
  copy_from(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free_base);

  /**
   * 当底层几何体发生变化（例如，通过像MappingFEField这样可以通过空间配置的变化而变形的映射）而网格和未知数的拓扑结构保持不变时，重新刷新存储在MappingInfo字段中的几何体数据。与reinit()相比，这个操作只需要重新生成几何体阵列，因此可以大大降低成本（取决于评估几何体的成本）。
   *
   */
  void
  update_mapping(const Mapping<dim> &mapping);

  /**
   * 与上述相同，但有 hp::MappingCollection. 。
   *
   */
  void
  update_mapping(const std::shared_ptr<hp::MappingCollection<dim>> &mapping);

  /**
   * 清除所有的数据字段，使类进入类似于调用了默认构造函数后的状态。
   *
   */
  void
  clear();

  //@}

  /**
   * 这个类定义了循环（）中的面积分的数据访问类型，它被传递给并行向量的`update_ghost_values`和`compress`函数，目的是能够减少必须交换的数据量。数据交换是一个真正的瓶颈，特别是对于高等级的DG方法，因此更严格的交换方式显然是有益的。请注意，这种选择只适用于分配给访问
   * `FaceToCellTopology::exterior_cells`
   * 的单元格外部的FEFaceEvaluation对象；所有<i>interior</i>对象在任何情况下都是可用的。
   *
   */
  enum class DataAccessOnFaces
  {
    /**
     * 循环不涉及任何FEFaceEvaluation对邻居的访问，就像只有边界积分（但没有内部面积分）或在类似
     * MatrixFree::cell_loop() 的设置中做质量矩阵时的情况。
     *
     */
    none,

    /**
     * 循环只涉及到FEFaceEvaluation通过函数值访问邻居，如
     * FEFaceEvaluation::gather_evaluate() 的参数 EvaluationFlags::values,
     * ，但没有访问形状函数导数（通常需要访问更多数据）。对于只有部分形状函数在面上有支持的FiniteElement类型，如FE_DGQ元素，其节点在元素表面的拉格朗日多项式，数据交换从`(k+1)^dim`减少到`(k+1)^(dim-1)`。
     *
     */
    values,

    /**
     * 与上述相同。如果FEFaceEvaluation通过提供单元批号和一个面的编号而被重新初始化，则必须从外部面访问数据时使用。这种配置在以单元为中心的循环中很有用。
     * @pre   AdditionalData::hold_all_faces_to_owned_cells  必须启用。
     *
     */
    values_all_faces,

    /**
     * 循环确实涉及到FEFaceEvaluation通过函数值和梯度访问到邻居，但没有二阶导数，比如
     * FEFaceEvaluation::gather_evaluate() 设置了 EvaluationFlags::values
     * 和 EvaluationFlags::gradients
     * 。对于只有部分形状函数在一个面上有非零值和一阶导数的FiniteElement类型，例如FE_DGQHermite元素，数据交换会减少，例如从`(k+1)^dim`到`2(k+1)^(dim-1)`。请注意，对于不具备这种特殊属性的基数，无论如何都要发送完整的邻接数据。
     *
     */
    gradients,

    /**
     * 与上述相同。如果FEFaceEvaluation通过提供单元格批号和一个面的编号而被重新初始化的话，要用于必须从外部面访问数据。这种配置在以单元为中心的循环中很有用。
     * @pre   AdditionalData::hold_all_faces_to_owned_cells  必须启用。
     *
     */
    gradients_all_faces,

    /**
     * 用户不想做限制的一般设置。这通常比其他选项更昂贵，但也是最保守的一种，因为要在本地计算的面后面的元素的全部数据将被交换。
     *
     */
    unspecified
  };

  /**
   * @name  2：无矩阵循环
   *
   */
  //@{
  /**
   * 这种方法在所有单元格上运行循环（并行），并在源向量和目的向量上执行MPI数据交换。
   * @param  cell_operation  `std::function` ，签名为<tt>cell_operation
   * (const MatrixFree<dim,Number> &, OutVector &, InVector &,
   * std::pair<unsigned  int,unsigned int>
   * &)</tt>，第一个参数传递调用类的数据，最后一个参数定义应该被处理的单元的范围（通常应该处理多个单元以减少开销）。
   * 如果一个对象有一个具有正确参数集的`operator()`，人们可以在这个地方传递一个指针，因为这样的指针可以被转换为函数对象。
   * @param  dst 保存结果的目标向量。如果向量是
   * LinearAlgebra::distributed::Vector 类型（或其复合对象，如
   * LinearAlgebra::distributed::BlockVector),
   * ，循环在内部调用结束时调用
   * LinearAlgebra::distributed::Vector::compress()
   * 。对于其他向量，包括平行的Trilinos或PETSc向量，不发出这样的调用。请注意，Trilinos/Epetra或PETSc向量目前不能并行工作，因为本类使用MPI本地索引寻址，而不是这些外部库所暗示的全局寻址。
   * @param  src 输入向量。如果向量是
   * LinearAlgebra::distributed::Vector 类型（或其复合对象，如
   * LinearAlgebra::distributed::BlockVector),
   * ，则循环在内部调用开始时调用
   * LinearAlgebra::distributed::Vector::update_ghost_values()
   * ，以确保所有必要的数据在本地可用。然而，请注意，在循环结束时，向量会被重置为原始状态，即如果向量在进入循环时没有被重影，那么在完成循环时也不会被重影。
   * @param  zero_dst_vector
   * 如果这个标志被设置为`true`，向量`dst`将在循环内被设置为零。在你对矩阵对象进行典型的`vmult()`操作时使用这种情况，因为它通常会比在循环之前单独调用`dst
   * =
   * 0;`更快。这是因为向量项只在向量的子范围内被设置为零，确保向量项尽可能地留在缓存中。
   *
   */
  template <typename OutVector, typename InVector>
  void
  cell_loop(const std::function<void(
              const MatrixFree<dim, Number, VectorizedArrayType> &,
              OutVector &,
              const InVector &,
              const std::pair<unsigned int, unsigned int> &)> &cell_operation,
            OutVector &                                        dst,
            const InVector &                                   src,
            const bool zero_dst_vector = false) const;

  /**
   * 这是第二个变种，在所有单元格上运行循环，现在提供一个指向`CLASS`类成员函数的函数指针。这种方法避免了定义lambda函数或调用
   * std::bind
   * 将类绑定到给定的函数中，以防本地函数需要访问类中的数据（即，它是一个非静态成员函数）。
   * @param  cell_operation
   * 指向`CLASS`的成员函数，其签名为<tt>cell_operation (const
   * MatrixFree<dim,Number> &, OutVector &, InVector &,  std::pair<unsigned
   * int,unsigned int>
   * &)</tt>，其中第一个参数传递调用类的数据，最后一个参数定义应该被处理的单元范围（通常应该处理一个以上的单元以减少开销）。
   * @param  owning_class
   * 提供`cell_operation`调用的对象。为了与该接口兼容，该类必须允许调用`owning_class->cell_operation(..)`。
   * @param  dst 保存结果的目标向量。如果向量是
   * LinearAlgebra::distributed::Vector 类型（或其复合对象，如
   * LinearAlgebra::distributed::BlockVector),
   * ，循环在内部调用结束时调用
   * LinearAlgebra::distributed::Vector::compress()
   * 。对于其他向量，包括平行的Trilinos或PETSc向量，不发出这样的调用。请注意，Trilinos/Epetra或PETSc向量目前不能并行工作，因为本类使用MPI本地索引寻址，而不是那些外部库所暗示的全局寻址。
   * @param  src 输入向量。如果向量是
   * LinearAlgebra::distributed::Vector 类型（或其复合对象，如
   * LinearAlgebra::distributed::BlockVector),
   * ），循环在内部调用开始时调用
   * LinearAlgebra::distributed::Vector::update_ghost_values()
   * ，以确保所有必要的数据在本地可用。然而，请注意，在循环结束时，向量会被重置为原始状态，即如果向量在进入循环时没有被重影，那么在完成循环时也不会被重影。
   * @param  zero_dst_vector
   * 如果这个标志被设置为`true`，向量`dst`将在循环中被设置为零。在你对矩阵对象进行典型的`vmult()`操作时使用这种情况，因为它通常会比在循环之前单独调用`dst
   * =
   * 0;`更快。这是因为向量项只在向量的子范围内被设置为零，确保向量项尽可能地留在缓存中。
   *
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  cell_loop(void (CLASS::*cell_operation)(
              const MatrixFree &,
              OutVector &,
              const InVector &,
              const std::pair<unsigned int, unsigned int> &) const,
            const CLASS *   owning_class,
            OutVector &     dst,
            const InVector &src,
            const bool      zero_dst_vector = false) const;

  /**
   * 与上述相同，但对于非const的类成员函数。
   *
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  cell_loop(void (CLASS::*cell_operation)(
              const MatrixFree &,
              OutVector &,
              const InVector &,
              const std::pair<unsigned int, unsigned int> &),
            CLASS *         owning_class,
            OutVector &     dst,
            const InVector &src,
            const bool      zero_dst_vector = false) const;

  /**
   * 这个函数类似于cell_loop，用一个 std::function
   * 对象来指定要对单元格进行的操作，但增加了两个额外的向量，在计算单元格积分前后执行一些额外的工作。
   * 这两个额外的向量在自由度范围内工作，以MPI本地索引中选定的DoFHandler`dof_handler_index_pre_post`的自由度编号表示。向量的参数代表自由度范围，粒度为
   * internal::MatrixFreeFunctions::DoFInfo::chunk_size_zero_vector
   * 个条目（除了最后一个块被设置为本地拥有的条目数），形式为`[first,
   * last)`。这些向量的想法是使向量上的操作更接近于它们在无矩阵循环中的访问点，目的是通过时间上的定位来增加缓存的点击。这个循环保证了`operation_before_loop`在cell_operation（包括MPI数据交换）中第一次接触到所有相关的未知数之前，就已经击中了这些未知数，允许执行一些`src`向量依赖的向量更新。循环后的操作
   * "是类似的
   *
   *
   *
   *
   * - 一旦该范围内的所有自由度被`cell_operation`最后一次触及（包括MPI数据交换），它就开始在该范围的自由度上执行，允许例如计算一些取决于当前`dst`中的单元循环结果的向量操作，或者想修改`src`。缓存的效率取决于自由度的编号，因为范围的粒度不同。      @param  cell_operation 指向`CLASS`的成员函数，签名为<tt>cell_operation (const MatrixFree<dim,Number> &, OutVector &, InVector &,  std::pair<unsigned  int,unsigned int> &)</tt>，第一个参数传递调用类的数据，最后一个参数定义应该被处理的单元范围（通常应该处理多个单元，以减少开销）。      @param  owning_class 提供`cell_operation`调用的对象。为了与该接口兼容，该类必须允许调用`owning_class->cell_operation(..)`。      @param  dst 保存结果的目标向量。如果向量是 LinearAlgebra::distributed::Vector 类型（或其复合对象，如 LinearAlgebra::distributed::BlockVector), ，循环在内部调用结束时调用 LinearAlgebra::distributed::Vector::compress() 。对于其他向量，包括平行的Trilinos或PETSc向量，不发出这样的调用。请注意，Trilinos/Epetra或PETSc向量目前不能并行工作，因为本类使用MPI本地索引寻址，而不是这些外部库所暗示的全局寻址。      @param  src 输入向量。如果向量是 LinearAlgebra::distributed::Vector 类型（或其复合对象，如 LinearAlgebra::distributed::BlockVector), ），循环在内部调用开始时调用 LinearAlgebra::distributed::Vector::update_ghost_values() ，以确保所有必要的数据在本地可用。然而，请注意，在循环结束时，向量被重置为其原始状态，即如果向量在进入循环时没有被重影，那么在完成循环时也不会被重影。      @param  operation_before_loop 这个函数可以用来对`src'和`dst'向量（或其他向量）的条目进行操作，在对单元的操作第一次接触到特定的DoF之前，根据上面文本中的一般描述。这个函数被传递给选定的`dof_handler_index_pre_post`（用MPI本地编号）上的本地拥有的自由度范围。      @param  operation_after_loop 这个函数可以用来对`src'和`dst'向量（或其他向量）的条目进行操作，在对单元的操作最后触及一个特定的DoF之后，根据上面文字的一般描述。这个函数被传递给选定的`dof_handler_index_pre_post`（以MPI本地编号）上的本地拥有的自由度范围。      @param  dof_handler_index_pre_post 由于MatrixFree可以用DoFHandler对象的矢量初始化，一般来说，每个对象都会有矢量大小，因此返回给`operation_before_loop`和`operation_after_loop`的范围也不同。使用这个变量来指定索引范围应该与哪一个DoFHandler对象相关。默认为`dof_handler_index` 0。
   * @note
   * `operation_before_loop`和`operation_after_loop`的近距离定位目前只在仅MPI的情况下实现。在启用线程的情况下，由于复杂的依赖关系，完整的`operation_before_loop`被安排在并行循环之前，而`operation_after_loop`被严格安排在之后。
   *
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  cell_loop(void (CLASS::*cell_operation)(
              const MatrixFree &,
              OutVector &,
              const InVector &,
              const std::pair<unsigned int, unsigned int> &) const,
            const CLASS *   owning_class,
            OutVector &     dst,
            const InVector &src,
            const std::function<void(const unsigned int, const unsigned int)>
              &operation_before_loop,
            const std::function<void(const unsigned int, const unsigned int)>
              &                operation_after_loop,
            const unsigned int dof_handler_index_pre_post = 0) const;

  /**
   * 与上述相同，但对于非const的类成员函数。
   *
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  cell_loop(void (CLASS::*cell_operation)(
              const MatrixFree &,
              OutVector &,
              const InVector &,
              const std::pair<unsigned int, unsigned int> &),
            CLASS *         owning_class,
            OutVector &     dst,
            const InVector &src,
            const std::function<void(const unsigned int, const unsigned int)>
              &operation_before_loop,
            const std::function<void(const unsigned int, const unsigned int)>
              &                operation_after_loop,
            const unsigned int dof_handler_index_pre_post = 0) const;

  /**
   * 同上，但取一个 `std::function`
   * 作为`cell_operation`，而不是类成员函数。
   *
   */
  template <typename OutVector, typename InVector>
  void
  cell_loop(const std::function<void(
              const MatrixFree<dim, Number, VectorizedArrayType> &,
              OutVector &,
              const InVector &,
              const std::pair<unsigned int, unsigned int> &)> &cell_operation,
            OutVector &                                        dst,
            const InVector &                                   src,
            const std::function<void(const unsigned int, const unsigned int)>
              &operation_before_loop,
            const std::function<void(const unsigned int, const unsigned int)>
              &                operation_after_loop,
            const unsigned int dof_handler_index_pre_post = 0) const;

  /**
   * 这个方法在所有单元格上运行一个循环（并行），并在源向量和目的向量上执行MPI数据交换。与其他只在单元格上运行一个函数的变体不同，这个方法还分别以内部面和边界面的函数为参数。
   * @param  cell_operation  `std::function`  的签名为<tt>cell_operation
   * (const MatrixFree<dim,Number> &, OutVector &, InVector &,
   * std::pair<unsigned  int,unsigned int>
   * &)</tt>，第一个参数传递调用类的数据，最后一个参数定义应该被处理的单元的范围（通常应该处理一个以上的单元以减少开销）。如果有一个
   * <code>operator()</code>
   * 具有正确的参数集，人们可以在这个地方传递一个对象的指针，因为这样的指针可以被转换为函数对象。
   * @param  face_operation  `std::function`  的签名是<tt>face_operation
   * (const MatrixFree<dim,Number> &, OutVector &, InVector &,
   * std::pair<unsigned  int,unsigned int>
   * &)</tt>，类似于`cell_operation`，但现在是与内部面的工作有关的部分。请注意，MatrixFree框架将周期性面视为内部面，因此，在调用face_operation时应用周期性约束后，它们将被分配到正确的邻居。
   * @param  boundary_operation  `std::function`
   * 签名为<tt>boundary_operation (const MatrixFree<dim,Number> &,
   * OutVector &, InVector &,  std::pair<unsigned  int,unsigned int>
   * &)</tt>，与`cell_operation`和`face_operation`类似，但现在是与边界面的工作相关的部分。边界面是由它们的`boundary_id'分开的，可以用
   * MatrixFree::get_boundary_id().
   * 来查询这个id。注意，内部和面都使用相同的编号，内部的面被分配的编号比边界面低。
   * @param  dst 保存结果的目的向量。如果向量是
   * LinearAlgebra::distributed::Vector 类型（或其复合对象，如
   * LinearAlgebra::distributed::BlockVector),
   * ，循环在内部调用结束时调用
   * LinearAlgebra::distributed::Vector::compress() 。      @param  src
   * 输入向量。如果向量是 LinearAlgebra::distributed::Vector
   * 类型（或其复合对象，如
   * LinearAlgebra::distributed::BlockVector),
   * ，则循环在内部调用开始时调用
   * LinearAlgebra::distributed::Vector::update_ghost_values()
   * 以确保所有必要的数据在本地可用。然而，请注意，在循环结束时，向量会被重置为其原始状态，即如果向量在进入循环时没有被重影，那么在完成循环时也不会被重影。
   * @param  zero_dst_vector
   * 如果这个标志被设置为`true`，向量`dst`将在循环中被设置为零。在你对矩阵对象进行典型的`vmult()`操作时使用这种情况，因为它通常会比在循环之前单独调用`dst
   * =
   * 0;`更快。这是因为向量项只在向量的子范围内被设置为0，确保向量项尽可能地留在缓存中。
   * @param  dst_vector_face_access
   * 设置对向量`dst`的访问类型，该访问将发生在 @p
   * face_operation
   * 函数内部。正如在DataAccessOnFaces结构的描述中所解释的，这个选择的目的是减少必须通过MPI网络（如果在节点的共享内存区域内，则通过`memcpy`）交换的数据量以获得性能。请注意，没有办法与FEFaceEvaluation类沟通这一设置，因此除了在`face_operation`函数中实现的内容外，这一选择必须在这个地方进行。因此，也没有办法检查传递给这个调用的设置是否与后来`FEFaceEvaluation`所做的一致，确保数据的正确性是用户的责任。
   * @param  src_vector_face_access
   * 设置对向量`src`的访问类型，将在 @p face_operation
   * 函数体内发生，与`dst_vector_face_access`相类似。
   *
   */
  template <typename OutVector, typename InVector>
  void
  loop(const std::function<
         void(const MatrixFree<dim, Number, VectorizedArrayType> &,
              OutVector &,
              const InVector &,
              const std::pair<unsigned int, unsigned int> &)> &cell_operation,
       const std::function<
         void(const MatrixFree<dim, Number, VectorizedArrayType> &,
              OutVector &,
              const InVector &,
              const std::pair<unsigned int, unsigned int> &)> &face_operation,
       const std::function<void(
         const MatrixFree<dim, Number, VectorizedArrayType> &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &)> &boundary_operation,
       OutVector &                                        dst,
       const InVector &                                   src,
       const bool              zero_dst_vector = false,
       const DataAccessOnFaces dst_vector_face_access =
         DataAccessOnFaces::unspecified,
       const DataAccessOnFaces src_vector_face_access =
         DataAccessOnFaces::unspecified) const;

  /**
   * 这是第二个变体，在所有的单元格、内部面和边界面上运行循环，现在提供了三个指向
   * @p CLASS 类成员函数的函数指针，其签名为<code>operation
   * (const MatrixFree<dim,Number> &, OutVector &, InVector &,
   * std::pair<unsigned  int,unsigned int>&)
   * const</code>。如果本地函数需要访问类中的数据（即，它是一个非静态成员函数），该方法就不需要定义lambda函数或调用
   * std::bind 将类绑定到给定函数中。      @param  cell_operation
   * 指向`CLASS`的成员函数，其签名为<tt>cell_operation (const
   * MatrixFree<dim,Number> &, OutVector &, InVector &,  std::pair<unsigned
   * int,unsigned int>
   * &)</tt>，其中第一个参数传递调用类的数据，最后一个参数定义应该被处理的单元范围（通常应该处理一个以上的单元以减少开销）。注意，这个循环通常会将`cell_range'分割成更小的部分，并交替工作于`cell_operation'、`face_operation'和`boundary_operation'，以增加缓存中向量项的潜在重用。
   * @param  face_operation
   * 指向`CLASS`的成员函数，其签名为<tt>face_operation (const
   * MatrixFree<dim,Number> &, OutVector &, InVector &,  std::pair<unsigned
   * int,unsigned int>
   * &)</tt>，与`cell_operation`相类似，但现在是与内部面的工作相关的部分。请注意，MatrixFree框架将周期性面视为内部面，因此，在调用face_operation时应用周期性约束后，它们将被分配给正确的邻居。
   * @param  boundary_operation
   * 指向`CLASS`的成员函数，其签名为<tt>boundary_operation (const
   * MatrixFree<dim,Number> &, OutVector &, InVector &,  std::pair<unsigned
   * int,unsigned int>
   * &)</tt>，与`cell_operation`和`face_operation`类似，但现在是与边界面的工作相关的部分。边界面由它们的`boundary_id'分开，可以使用
   * MatrixFree::get_boundary_id().
   * 查询该id。注意，内部和面都使用相同的编号，内部的面被分配的编号比边界面低。
   * @param  owning_class
   * 提供`cell_operation`调用的对象。为了与该接口兼容，该类必须允许调用`owning_class->cell_operation(...)`,
   * `owning_class->face_operation(...)`,
   * 和`owning_class->boundary_operation(...)`。      @param  dst
   * 保存结果的目标向量。如果向量是
   * LinearAlgebra::distributed::Vector 类型（或其复合对象，如
   * LinearAlgebra::distributed::BlockVector),
   * ，循环在内部调用结束时调用
   * LinearAlgebra::distributed::Vector::compress() 。      @param  src
   * 输入向量。如果向量是 LinearAlgebra::distributed::Vector
   * 类型（或其复合对象，如
   * LinearAlgebra::distributed::BlockVector),
   * ，则循环在内部调用开始时调用
   * LinearAlgebra::distributed::Vector::update_ghost_values()
   * 以确保所有必要的数据在本地可用。然而，请注意，在循环结束时，向量会被重置为其原始状态，即如果向量在进入循环时没有被重影，那么在完成循环时也不会被重影。
   * @param  zero_dst_vector
   * 如果这个标志被设置为`true`，向量`dst`将在循环中被设置为零。在你对矩阵对象进行典型的`vmult()`操作时使用这种情况，因为它通常会比在循环之前单独调用`dst
   * =
   * 0;`更快。这是因为向量项只在向量的子范围内被设置为0，确保向量项尽可能地留在缓存中。
   * @param  dst_vector_face_access
   * 设置对向量`dst`的访问类型，将在 @p face_operation
   * 函数内部发生。正如在DataAccessOnFaces结构的描述中所解释的，这个选择的目的是减少必须通过MPI网络（如果在节点的共享内存区域内，则通过`memcpy`）交换的数据量以获得性能。请注意，没有办法与FEFaceEvaluation类沟通这一设置，因此除了在`face_operation`函数中实现的内容外，这一选择必须在这个地方进行。因此，也没有办法检查传递给这个调用的设置是否与后来`FEFaceEvaluation`所做的一致，确保数据的正确性是用户的责任。
   * @param  src_vector_face_access
   * 设置对向量`src`的访问类型，将在 @p face_operation
   * 函数体内发生，与`dst_vector_face_access`相类似。
   *
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  loop(
    void (CLASS::*cell_operation)(const MatrixFree &,
                                  OutVector &,
                                  const InVector &,
                                  const std::pair<unsigned int, unsigned int> &)
      const,
    void (CLASS::*face_operation)(const MatrixFree &,
                                  OutVector &,
                                  const InVector &,
                                  const std::pair<unsigned int, unsigned int> &)
      const,
    void (CLASS::*boundary_operation)(
      const MatrixFree &,
      OutVector &,
      const InVector &,
      const std::pair<unsigned int, unsigned int> &) const,
    const CLASS *           owning_class,
    OutVector &             dst,
    const InVector &        src,
    const bool              zero_dst_vector = false,
    const DataAccessOnFaces dst_vector_face_access =
      DataAccessOnFaces::unspecified,
    const DataAccessOnFaces src_vector_face_access =
      DataAccessOnFaces::unspecified) const;

  /**
   * 和上面一样，但对于非const的类成员函数。
   *
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  loop(void (CLASS::*cell_operation)(
         const MatrixFree &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &),
       void (CLASS::*face_operation)(
         const MatrixFree &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &),
       void (CLASS::*boundary_operation)(
         const MatrixFree &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &),
       CLASS *                 owning_class,
       OutVector &             dst,
       const InVector &        src,
       const bool              zero_dst_vector = false,
       const DataAccessOnFaces dst_vector_face_access =
         DataAccessOnFaces::unspecified,
       const DataAccessOnFaces src_vector_face_access =
         DataAccessOnFaces::unspecified) const;

  /**
   * 这个方法与cell_loop()的做法类似，在所有单元格上运行循环（并行）。然而，这个函数的目的是用于面和边界积分也应该被评估的情况。与loop()相反，用户只需提供一个单一的函数，该函数应包含一个单元（或向量化时的一批单元）的单元积分和所有面的面和边界积分。这在文献中被称为
   * "以元素为中心的循环 "或 "以单元为中心的循环"。
   * 为了能够评估所有的面积分（用来自相邻单元的值或梯度），相邻单元的所有幽灵值都要更新。使用
   * FEFaceEvalution::reinit(cell,
   * face_no)来访问一个单元的任意面和各自的邻居的数量。
   * @param  cell_operation
   * 指向`CLASS`的成员函数，其签名为<tt>cell_operation (const
   * MatrixFree<dim,Number> &, OutVector &, InVector &,  std::pair<unsigned
   * int,unsigned int>
   * &)</tt>，第一个参数传递调用类的数据，最后一个参数定义应该被处理的单元范围（通常从循环中传递多个单元以减少开销）。
   * @param  owning_class
   * 提供`cell_operation`调用的对象。为了与该接口兼容，该类必须允许调用`owning_class->cell_operation(..)`。
   * @param  dst 保存结果的目标向量。如果向量是
   * LinearAlgebra::distributed::Vector 类型（或其复合对象，如
   * LinearAlgebra::distributed::BlockVector),
   * ，循环在内部调用结束时调用
   * LinearAlgebra::distributed::Vector::compress() 。      @param  src
   * 输入向量。如果向量是 LinearAlgebra::distributed::Vector
   * 类型（或其复合对象，如
   * LinearAlgebra::distributed::BlockVector),
   * ，循环在内部调用开始时调用
   * LinearAlgebra::distributed::Vector::update_ghost_values()
   * ，以确保所有必要的数据在本地可用。然而，请注意，在循环结束时，向量会被重置为原始状态，即如果向量在进入循环时没有被重影，那么在完成循环时也不会被重影。
   * @param  zero_dst_vector
   * 如果这个标志被设置为`true`，向量`dst`将在循环中被设置为零。在你对矩阵对象进行典型的`vmult()`操作时使用这种情况，因为它通常会比在循环之前单独调用`dst
   * =
   * 0;`更快。这是因为向量项只在向量的子范围内被设置为0，确保向量项尽可能地留在缓存中。
   * @param  src_vector_face_access
   * 设置对向量`src`的访问类型，这将在面积分过程中发生在
   * @p cell_operation 函数的内部。
   * 正如在DataAccessOnFaces结构的描述中所解释的，这个选择的目的是减少必须通过MPI网络（如果在节点的共享内存区域内，则通过`memcpy`）交换的数据量以获得性能。请注意，没有办法与FEFaceEvaluation类沟通这一设置，因此除了在`face_operation`函数中实现的内容外，这一选择必须在这个地方进行。因此，也没有办法检查传递给这个调用的设置是否与后来`FEFaceEvaluation`所做的一致，确保数据的正确性是用户的责任。
   *
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  loop_cell_centric(void (CLASS::*cell_operation)(
                      const MatrixFree &,
                      OutVector &,
                      const InVector &,
                      const std::pair<unsigned int, unsigned int> &) const,
                    const CLASS *           owning_class,
                    OutVector &             dst,
                    const InVector &        src,
                    const bool              zero_dst_vector = false,
                    const DataAccessOnFaces src_vector_face_access =
                      DataAccessOnFaces::unspecified) const;

  /**
   * 同上，但对于类的成员函数，它是非const的。
   *
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  loop_cell_centric(void (CLASS::*cell_operation)(
                      const MatrixFree &,
                      OutVector &,
                      const InVector &,
                      const std::pair<unsigned int, unsigned int> &),
                    CLASS *                 owning_class,
                    OutVector &             dst,
                    const InVector &        src,
                    const bool              zero_dst_vector = false,
                    const DataAccessOnFaces src_vector_face_access =
                      DataAccessOnFaces::unspecified) const;

  /**
   * 与上述相同，但有 std::function. 。
   *
   */
  template <typename OutVector, typename InVector>
  void
  loop_cell_centric(
    const std::function<void(const MatrixFree &,
                             OutVector &,
                             const InVector &,
                             const std::pair<unsigned int, unsigned int> &)>
      &                     cell_operation,
    OutVector &             dst,
    const InVector &        src,
    const bool              zero_dst_vector = false,
    const DataAccessOnFaces src_vector_face_access =
      DataAccessOnFaces::unspecified) const;

  /**
   * 在hp-adaptive情况下，在单元格循环中计算的单元格子范围可能包含不同程度的元素。使用这个函数来计算一个单独的有限元度的子范围是什么。有限元度与函数调用中给出的矢量分量相关。
   *
   */
  std::pair<unsigned int, unsigned int>
  create_cell_subrange_hp(const std::pair<unsigned int, unsigned int> &range,
                          const unsigned int fe_degree,
                          const unsigned int dof_handler_index = 0) const;

  /**
   * 在hp-adaptive情况下，在单元循环中计算的单元子范围可能包含不同程度的元素。使用这个函数来计算给定索引的子范围的hp-finite元素，而不是其他函数中的有限元素程度。
   *
   */
  std::pair<unsigned int, unsigned int>
  create_cell_subrange_hp_by_index(
    const std::pair<unsigned int, unsigned int> &range,
    const unsigned int                           fe_index,
    const unsigned int                           dof_handler_index = 0) const;

  /**
   * 在hp自适应情况下，返回active_fe_indices的数量。
   *
   */
  unsigned int
  n_active_fe_indices() const;

  /**
   * 在hp-adaptive情况下，返回单元格范围的active_fe_index。
   *
   */
  unsigned int
  get_cell_active_fe_index(
    const std::pair<unsigned int, unsigned int> range) const;

  /**
   * 在hp-adaptive的情况下，返回一个面域的active_fe_index。
   *
   */
  unsigned int
  get_face_active_fe_index(const std::pair<unsigned int, unsigned int> range,
                           const bool is_interior_face = true) const;

  //@}

  /**
   * @name  3: 向量的初始化
   *
   */
  //@{
  /**
   *
   */
  template <typename VectorType>
  void
  initialize_dof_vector(VectorType &       vec,
                        const unsigned int dof_handler_index = 0) const;

  /**
   *
   */
  template <typename Number2>
  void
  initialize_dof_vector(LinearAlgebra::distributed::Vector<Number2> &vec,
                        const unsigned int dof_handler_index = 0) const;

  /**
   * 返回代表本地拥有的数据的分区器，以及单元格循环需要访问的幽灵索引。分区器是由各自字段给出的本地拥有的道夫和幽灵道夫构建的。如果你想知道这些对象的具体信息，你可以用各自的访问函数来查询它们。如果你只是想初始化一个（平行）向量，你通常应该更喜欢这种数据结构，因为数据交换信息可以从一个向量重复使用到另一个向量。
   *
   */
  const std::shared_ptr<const Utilities::MPI::Partitioner> &
  get_vector_partitioner(const unsigned int dof_handler_index = 0) const;

  /**
   * 返回由处理器拥有的单元格集合。
   *
   */
  const IndexSet &
  get_locally_owned_set(const unsigned int dof_handler_index = 0) const;

  /**
   * 返回需要但不为处理器所拥有的幽灵单元的集合。
   *
   */
  const IndexSet &
  get_ghost_set(const unsigned int dof_handler_index = 0) const;

  /**
   * 返回一个所有被约束的自由度的列表。该列表是在本地拥有的向量范围的MPI本地索引空间中返回的，而不是跨越所有MPI处理器的全局MPI索引空间。要获得全局索引空间的数字，请在向量的一个条目上调用<tt>get_vector_partitioner()->local_to_global</tt>。此外，它只返回本地拥有的自由度的指数，而不是鬼魂的指数。
   *
   */
  const std::vector<unsigned int> &
  get_constrained_dofs(const unsigned int dof_handler_index = 0) const;

  /**
   * 根据给定的数据布局，计算自由度的重新编号，使其更符合MatrixFree中的数据布局。注意，这个函数并不重新排列存储在这个类中的信息，而是创建一个重新编号以消耗
   * DoFHandler::renumber_dofs.
   * 为了产生任何效果，必须使用重新编号的DoFHandler和AffineConstraints再次设置MatrixFree对象。注意，如果DoFHandler调用
   * DoFHandler::renumber_dofs, ，MatrixFree中的所有信息都会失效。
   *
   */
  void
  renumber_dofs(std::vector<types::global_dof_index> &renumbering,
                const unsigned int                    dof_handler_index = 0);

  //@}

  /**
   * @name  4：一般信息
   *
   */
  //@{
  /**
   * 返回一个给定的FiniteElement  @p fe  是否被这个类所支持。
   *
   */
  template <int spacedim>
  static bool
  is_supported(const FiniteElement<dim, spacedim> &fe);

  /**
   * 返回初始化时指定的不同DoFHandlers的数量。
   *
   */
  unsigned int
  n_components() const;

  /**
   * 对于由 @p
   * dof_handler_index指定的DoFHandler的基础有限元，返回基础元的数量。
   *
   */
  unsigned int
  n_base_elements(const unsigned int dof_handler_index) const;

  /**
   * 返回这个结构所基于的单元数。如果你使用的是一个通常的DoFHandler，它对应于（本地拥有的）活动单元的数量。请注意，这个类中的大多数数据结构并不直接作用于这个数字，而是作用于n_cell_batches()，它给出了用矢量化将几个单元拼凑在一起时看到的单元的数量。
   *
   */
  unsigned int
  n_physical_cells() const;

  /**
   * @deprecated  用n_cell_batches()代替。
   *
   */
  DEAL_II_DEPRECATED unsigned int
  n_macro_cells() const;

  /**
   * 返回该结构所处理的单元格批次的数量。批次是通过在一般的几个单元上应用矢量化而形成的。
   * @p cell_loop
   * 中的细胞范围从零到n_cell_batches()（独占），所以如果你想为所有要处理的细胞存储数据数组，这是一个合适的大小。这个数字大约是
   * `n_physical_cells()/VectorizedArray::%size()`
   * （取决于有多少细胞批没有被完全填满）。
   *
   */
  unsigned int
  n_cell_batches() const;

  /**
   * 返回该结构为面层集成而保留的额外单元格批次的数量。请注意，并不是所有在三角形中被重影的单元格都被保留在这个数据结构中，而是只有那些对评估两边的面积分有必要的单元格。
   *
   */
  unsigned int
  n_ghost_cell_batches() const;

  /**
   * 返回这个结构所处理的内部面批的数量。
   * 这些批次是通过在一般的几个面上应用矢量化而形成的。
   * @p loop
   * 中的面的范围从零到n_inner_face_batches()（独占），所以如果你想为所有要处理的内部面存储数据的数组，这就是合适的大小。
   *
   */
  unsigned int
  n_inner_face_batches() const;

  /**
   * 返回这个结构所处理的边界面批的数量。
   * 这些批次是通过在一般的几个面上应用矢量化而形成的。
   * @p loop
   * 中的面的范围从n_inner_face_batches()到n_inner_face_batches()+n_boundary_face_batches()（独占），所以如果你需要存储所有边界面而不是内部面的数据的数组，这个数字给出适当的大小。
   *
   */
  unsigned int
  n_boundary_face_batches() const;

  /**
   * 返回未在本地处理但属于本地拥有的面的数量。
   *
   */
  unsigned int
  n_ghost_inner_face_batches() const;

  /**
   * 为了对边界的不同部分应用不同的运算符，这个方法可以用来查询面孔自己在VectorizedArray中按车道排序的边界ID。只对表示边界面的索引有效。
   *
   */
  types::boundary_id
  get_boundary_id(const unsigned int macro_face) const;

  /**
   * 返回一个单元格内的面的边界ID，使用单元格在VectorizedArray中按车道的排序。
   *
   */
  std::array<types::boundary_id, VectorizedArrayType::size()>
  get_faces_by_cells_boundary_id(const unsigned int cell_batch_index,
                                 const unsigned int face_number) const;

  /**
   * 返回DoFHandler，其索引与reinit()函数中各自的 `std::vector`
   * 参数一样。
   *
   */
  const DoFHandler<dim> &
  get_dof_handler(const unsigned int dof_handler_index = 0) const;

  /**
   * 返回DoFHandler，其索引与reinit()函数中相应的 `std::vector`
   * 参数一样。注意，如果你想用不同于默认的模板参数来调用这个函数，你需要在函数调用前使用`template`，也就是说，你会有类似`matrix_free.template
   * get_dof_handler<hp::DoFHandler<dim>>()`.   @deprecated
   * 使用这个函数的非模板化等价物。
   *
   */
  template <typename DoFHandlerType>
  DEAL_II_DEPRECATED const DoFHandlerType &
                           get_dof_handler(const unsigned int dof_handler_index = 0) const;

  /**
   * 返回deal.II中的单元格迭代器讲到一个给定的单元格批次（在一个VectorizedArray中填充几个车道）和在这个结构的重新编号中跨单元格的矢量化中的车道索引。
   * 请注意，deal.II中的单元格迭代器与本类的单元格循环的处理方式不同。这是因为几个单元被一起处理（跨单元的矢量化），而且在访问远程数据和与计算重叠的通信时，在不同的MPI处理器上有邻居的单元需要在某个时间被访问。
   *
   */
  typename DoFHandler<dim>::cell_iterator
  get_cell_iterator(const unsigned int cell_batch_index,
                    const unsigned int lane_index,
                    const unsigned int dof_handler_index = 0) const;

  /**
   * 这将返回由get_cell_iterator()对相同参数`cell_batch_index`和`lane_index`所返回的单元的级别和索引。
   *
   */
  std::pair<int, int>
  get_cell_level_and_index(const unsigned int cell_batch_index,
                           const unsigned int lane_index) const;

  /**
   * 在deal.II中返回单元格的迭代器，该迭代器与一个面的内部/外部单元格在一对面批和巷道索引中。这一对中的第二个元素是面的编号，这样就可以访问面的迭代器。
   * `pair.first()->face(pair.second());`注意deal.II中的面迭代器通过单元格的方式与本类的面/边界循环的方式不同。这是因为几个面是一起工作的（矢量化），而且在访问远程数据和与计算重叠的通信时，在不同的MPI处理器上有相邻单元的面需要在某个时间被访问。
   *
   */
  std::pair<typename DoFHandler<dim>::cell_iterator, unsigned int>
  get_face_iterator(const unsigned int face_batch_index,
                    const unsigned int lane_index,
                    const bool         interior     = true,
                    const unsigned int fe_component = 0) const;

  /**
   * @copydoc   MatrixFree::get_cell_iterator()   @deprecated
   * 使用get_cell_iterator()代替。
   *
   */
  DEAL_II_DEPRECATED typename DoFHandler<dim>::active_cell_iterator
  get_hp_cell_iterator(const unsigned int cell_batch_index,
                       const unsigned int lane_index,
                       const unsigned int dof_handler_index = 0) const;

  /**
   * 由于该类使用的是矢量数据类型，数据域中通常有多个值，因此可能会出现矢量类型的某些分量与网格中的实际单元不对应的情况。当只使用这个类时，通常不需要理会这个事实，因为这些值是用零填充的。然而，当这个类与访问单元格的deal.II混合时，需要注意。如果不是所有给定的`cell_batch_index`的`n_lanes`单元都对应于网格中的实际单元，有些只是为了填充的原因而出现，则该函数返回
   * @p true
   * 。要知道有多少单元被实际使用，可以使用函数n_active_entries_per_cell_batch()。
   *
   */
  bool
  at_irregular_cell(const unsigned int cell_batch_index) const;

  /**
   * @deprecated  使用n_active_entries_per_cell_batch()代替。
   *
   */
  DEAL_II_DEPRECATED unsigned int
  n_components_filled(const unsigned int cell_batch_number) const;

  /**
   * 这个查询返回在一个单元格批中的
   * `VectorizedArrayType::size()`
   * 个单元格中，有多少单元格是网格中的实际单元格，而不是因为填充的原因而出现的。对于大多数给定的n_cell_batches()中的单元格批次，这个数字等于
   * `VectorizedArrayType::size()`,
   * ，但在网格中可能有一个或几个单元格批次（数字不相加），其中只有一个批次中的一些单元格被使用，由函数at_irregular_cell()表示。
   *
   */
  unsigned int
  n_active_entries_per_cell_batch(const unsigned int cell_batch_index) const;

  /**
   * 使用这个函数找出在矢量化数据类型的长度上有多少个面对应于网格中的真实面（包括内部和边界面，因为这些面使用相同的索引，但范围不同）。对于大多数在n_inner_faces_batches()和n_boundary_face_batches()中给定的索引，这只是
   * @p vectorization_length
   * 个，但可能有一个或几个网格（数字不相加），其中有更少的这种道的填充。
   *
   */
  unsigned int
  n_active_entries_per_face_batch(const unsigned int face_batch_index) const;

  /**
   * 返回给定hp-index的每个单元的自由度数量。
   *
   */
  unsigned int
  get_dofs_per_cell(const unsigned int dof_handler_index  = 0,
                    const unsigned int hp_active_fe_index = 0) const;

  /**
   * 返回给定hp-index的每个单元的正交点的数量。
   *
   */
  unsigned int
  get_n_q_points(const unsigned int quad_index         = 0,
                 const unsigned int hp_active_fe_index = 0) const;

  /**
   * 返回给定hp-index的单元格每个面上的自由度数量。
   *
   */
  unsigned int
  get_dofs_per_face(const unsigned int dof_handler_index  = 0,
                    const unsigned int hp_active_fe_index = 0) const;

  /**
   * 返回给定hp-index的单元格每个面上的正交点的数量。
   *
   */
  unsigned int
  get_n_q_points_face(const unsigned int quad_index         = 0,
                      const unsigned int hp_active_fe_index = 0) const;

  /**
   * 返回给定hp-index的正交规则。
   *
   */
  const Quadrature<dim> &
  get_quadrature(const unsigned int quad_index         = 0,
                 const unsigned int hp_active_fe_index = 0) const;

  /**
   * 返回给定hp-index的正交规则。
   *
   */
  const Quadrature<dim - 1> &
  get_face_quadrature(const unsigned int quad_index         = 0,
                      const unsigned int hp_active_fe_index = 0) const;

  /**
   * 返回当前批次的电池被分配到的类别。对于非hp-DoFHandler类型，类别在字段
   * AdditionalData::cell_vectorization_category
   * 中的给定值之间运行，在hp-adaptive情况下返回活动的FE指数。
   *
   */
  unsigned int
  get_cell_category(const unsigned int cell_batch_index) const;

  /**
   * 返回当前一批面的两边的单元格上的类别。
   *
   */
  std::pair<unsigned int, unsigned int>
  get_face_category(const unsigned int macro_face) const;

  /**
   * 查询是否已经设置了索引。
   *
   */
  bool
  indices_initialized() const;

  /**
   * 查询单元格的几何相关信息是否已被设置。
   *
   */
  bool
  mapping_initialized() const;

  /**
   * 返回要处理的网格的级别。如果在活动单元上工作，返回
   * numbers::invalid_unsigned_int 。
   *
   */
  unsigned int
  get_mg_level() const;

  /**
   * 返回该类的内存消耗量的近似值，单位为字节。
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
  print_memory_consumption(StreamType &out) const;

  /**
   * 在给定的输出流中打印这个类的摘要。它集中在索引上，并不打印所有存储的数据。
   *
   */
  void
  print(std::ostream &out) const;

  //@}

  /**
   * @name  5: 访问内部数据结构 注意：专家模式，接口在不同版本之间不稳定。
   *
   */
  //@{
  /**
   * 返回任务图的信息。
   *
   */
  const internal::MatrixFreeFunctions::TaskInfo &
  get_task_info() const;

  /*返回单元格上与几何有关的信息。 
*
*/
  const internal::MatrixFreeFunctions::
    MappingInfo<dim, Number, VectorizedArrayType> &
    get_mapping_info() const;

  /**
   * 返回关于索引自由度的信息。
   *
   */
  const internal::MatrixFreeFunctions::DoFInfo &
  get_dof_info(const unsigned int dof_handler_index_component = 0) const;

  /**
   * 返回约束池中的权重数量。
   *
   */
  unsigned int
  n_constraint_pool_entries() const;

  /**
   * 返回一个指向约束池数据中第一个数字的指针，索引为
   * @p pool_index （与 @p constraint_pool_end()). 一起使用
   *
   */
  const Number *
  constraint_pool_begin(const unsigned int pool_index) const;

  /**
   * 返回一个指向约束池数据中最后一个数字的指针，索引为
   * @p pool_index （与 @p  constraint_pool_begin()一起使用）。
   *
   */
  const Number *
  constraint_pool_end(const unsigned int pool_index) const;

  /**
   * 返回给定hp-index的单元格信息。
   *
   */
  const internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &
  get_shape_info(const unsigned int dof_handler_index_component = 0,
                 const unsigned int quad_index                  = 0,
                 const unsigned int fe_base_element             = 0,
                 const unsigned int hp_active_fe_index          = 0,
                 const unsigned int hp_active_quad_index        = 0) const;

  /**
   * 返回一个面的连接信息。
   *
   */
  const internal::MatrixFreeFunctions::FaceToCellTopology<
    VectorizedArrayType::size()> &
  get_face_info(const unsigned int face_batch_index) const;


  /**
   * 返回将宏观单元格号、单元格内的面的索引和向量的单元格批内的索引这三者转化为面数组内的索引的表格。
   *
   */
  const Table<3, unsigned int> &
  get_cell_and_face_to_plain_faces() const;

  /**
   * 获取一个内部使用的抓取数据对象。请确保事后将你从这个对象获得的指针传递给release_scratch_data()函数，以释放它。这个接口被FEEvaluation对象用来存储其数据结构。    内部数据结构的组织是一个线程本地存储的向量列表。多个线程将各自获得一个单独的存储域和单独的向量，确保线程安全。获取和释放对象的机制类似于WorkStream的本地贡献机制，见 @ref workstream_paper "WorkStream论文"
   * 。
   *
   */
  AlignedVector<VectorizedArrayType> *
  acquire_scratch_data() const;

  /**
   * 使得抓取板的对象再次可用。
   *
   */
  void
  release_scratch_data(const AlignedVector<VectorizedArrayType> *memory) const;

  /**
   * 获取一个用于内部使用的划痕数据对象。确保事后将你从这个对象获得的指针传递给release_scratch_data_non_threadsafe()函数，以释放它。注意，与acquisition_scratch_data()相反，这个方法一次只能由一个线程调用，但与acquisition_scratch_data()相反，释放scratch数据的线程也有可能与获取它的线程不同。
   *
   */
  AlignedVector<Number> *
  acquire_scratch_data_non_threadsafe() const;

  /**
   * 使得从头数据的对象再次可用。
   *
   */
  void
  release_scratch_data_non_threadsafe(
    const AlignedVector<Number> *memory) const;

  //@}

private:
  /**
   * 这是实际的reinit函数，为DoFHandler的情况设置索引。
   *
   */
  template <typename number2, int q_dim>
  void
  internal_reinit(
    const std::shared_ptr<hp::MappingCollection<dim>> &    mapping,
    const std::vector<const DoFHandler<dim, dim> *> &      dof_handlers,
    const std::vector<const AffineConstraints<number2> *> &constraint,
    const std::vector<IndexSet> &                          locally_owned_set,
    const std::vector<hp::QCollection<q_dim>> &            quad,
    const AdditionalData &                                 additional_data);

  /**
   * 初始化DoFInfo中的字段和约束池，约束池中保存了所有不同的权重（不是DoFInfo的一部分，因为几个DoFInfo类可以有相同的权重，因此只需要存储一次）。
   *
   */
  template <typename number2>
  void
  initialize_indices(
    const std::vector<const AffineConstraints<number2> *> &constraint,
    const std::vector<IndexSet> &                          locally_owned_set,
    const AdditionalData &                                 additional_data);

  /**
   * 基于DoFHandler<dim>参数初始化DoFHandlers。
   *
   */
  void
  initialize_dof_handlers(
    const std::vector<const DoFHandler<dim, dim> *> &dof_handlers,
    const AdditionalData &                           additional_data);

  /**
   * 指向当前问题基础的DoFHandlers的指针。
   *
   */
  std::vector<SmartPointer<const DoFHandler<dim>>> dof_handlers;

  /**
   * 包含各个单元格上的自由度和约束信息。
   *
   */
  std::vector<internal::MatrixFreeFunctions::DoFInfo> dof_info;

  /**
   * 包含存储在DoFInfo中的约束条件的权重。填充到一个单独的字段中，因为几个向量组件可能共享类似的权重，这样可以减少内存消耗。此外，它省去了DoFInfo上的模板参数，使之成为一个只有索引的普通字段。
   *
   */
  std::vector<Number> constraint_pool_data;

  /**
   * 包含一个指向约束池数据中第i个索引开始的指标。
   *
   */
  std::vector<unsigned int> constraint_pool_row_index;

  /**
   * 保存单元格从参考单元格到实数单元格的转换信息，这是评估积分所需要的。
   *
   */
  internal::MatrixFreeFunctions::MappingInfo<dim, Number, VectorizedArrayType>
    mapping_info;

  /**
   * 包含单元格的形状值信息。
   *
   */
  Table<4, internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType>>
    shape_info;

  /**
   * 描述了单元格是如何被穿过的。有了单元格级别（该字段的第一个索引）和级别内的索引，就可以重建一个deal.II单元格迭代器，并使用deal.II提供的所有传统的单元格迭代器的功能。
   *
   */
  std::vector<std::pair<unsigned int, unsigned int>> cell_level_index;


  /**
   * 对于非连续Galerkin，cell_level_index包括不在本地处理器上的单元，但需要用来计算单元积分。在cell_level_index_end_local中，我们存储本地单元的数量。
   *
   */
  unsigned int cell_level_index_end_local;

  /**
   * 存储要处理的单元和面的基本布局，包括共享内存并行化的任务布局以及与MPI的通信和计算之间可能的重叠。
   *
   */
  internal::MatrixFreeFunctions::TaskInfo task_info;

  /**
   * 持有面孔信息的向量。只在build_face_info=true时初始化。
   *
   */
  internal::MatrixFreeFunctions::FaceInfo<VectorizedArrayType::size()>
    face_info;

  /**
   * 存储索引是否已被初始化。
   *
   */
  bool indices_are_initialized;

  /**
   * 存储索引是否已被初始化。
   *
   */
  bool mapping_is_initialized;

  /**
   * 用于评估的刮板内存。我们允许超过一个评估对象附加到这个字段（这个，外
   * std::vector),
   * 所以我们需要跟踪某个数据字段是否已经被使用（第一部分对）并保持一个对象的列表。
   *
   */
  mutable Threads::ThreadLocalStorage<
    std::list<std::pair<bool, AlignedVector<VectorizedArrayType>>>>
    scratch_pad;

  /**
   * 用于评估和其他情况下的Scratchpad内存，非线程安全的变体。
   *
   */
  mutable std::list<std::pair<bool, AlignedVector<Number>>>
    scratch_pad_non_threadsafe;

  /**
   * 存储了要处理的网格的级别。
   *
   */
  unsigned int mg_level;
};



 /*----------------------- Inline functions ----------------------------------*/ 

#ifndef DOXYGEN



template <int dim, typename Number, typename VectorizedArrayType>
template <typename VectorType>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::initialize_dof_vector(
  VectorType &       vec,
  const unsigned int comp) const
{
  AssertIndexRange(comp, n_components());
  vec.reinit(dof_info[comp].vector_partitioner->size());
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename Number2>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::initialize_dof_vector(
  LinearAlgebra::distributed::Vector<Number2> &vec,
  const unsigned int                           comp) const
{
  AssertIndexRange(comp, n_components());
  vec.reinit(dof_info[comp].vector_partitioner, task_info.communicator_sm);
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const std::shared_ptr<const Utilities::MPI::Partitioner> &
MatrixFree<dim, Number, VectorizedArrayType>::get_vector_partitioner(
  const unsigned int comp) const
{
  AssertIndexRange(comp, n_components());
  return dof_info[comp].vector_partitioner;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const std::vector<unsigned int> &
MatrixFree<dim, Number, VectorizedArrayType>::get_constrained_dofs(
  const unsigned int comp) const
{
  AssertIndexRange(comp, n_components());
  return dof_info[comp].constrained_dofs;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_components() const
{
  AssertDimension(dof_handlers.size(), dof_info.size());
  return dof_handlers.size();
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_base_elements(
  const unsigned int dof_no) const
{
  AssertDimension(dof_handlers.size(), dof_info.size());
  AssertIndexRange(dof_no, dof_handlers.size());
  return dof_handlers[dof_no]->get_fe().n_base_elements();
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const internal::MatrixFreeFunctions::TaskInfo &
MatrixFree<dim, Number, VectorizedArrayType>::get_task_info() const
{
  return task_info;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_macro_cells() const
{
  return *(task_info.cell_partition_data.end() - 2);
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_physical_cells() const
{
  return task_info.n_active_cells;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_cell_batches() const
{
  return *(task_info.cell_partition_data.end() - 2);
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_ghost_cell_batches() const
{
  return *(task_info.cell_partition_data.end() - 1) -
         *(task_info.cell_partition_data.end() - 2);
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_inner_face_batches() const
{
  if (task_info.face_partition_data.size() == 0)
    return 0;
  return task_info.face_partition_data.back();
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_boundary_face_batches() const
{
  if (task_info.face_partition_data.size() == 0)
    return 0;
  return task_info.boundary_partition_data.back() -
         task_info.face_partition_data.back();
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_ghost_inner_face_batches() const
{
  if (task_info.face_partition_data.size() == 0)
    return 0;
  return face_info.faces.size() - task_info.boundary_partition_data.back();
}



template <int dim, typename Number, typename VectorizedArrayType>
inline types::boundary_id
MatrixFree<dim, Number, VectorizedArrayType>::get_boundary_id(
  const unsigned int macro_face) const
{
  Assert(macro_face >= task_info.boundary_partition_data[0] &&
           macro_face < task_info.boundary_partition_data.back(),
         ExcIndexRange(macro_face,
                       task_info.boundary_partition_data[0],
                       task_info.boundary_partition_data.back()));
  return types::boundary_id(face_info.faces[macro_face].exterior_face_no);
}



template <int dim, typename Number, typename VectorizedArrayType>
inline std::array<types::boundary_id, VectorizedArrayType::size()>
MatrixFree<dim, Number, VectorizedArrayType>::get_faces_by_cells_boundary_id(
  const unsigned int cell_batch_index,
  const unsigned int face_number) const
{
  AssertIndexRange(cell_batch_index, n_cell_batches());
  AssertIndexRange(face_number, GeometryInfo<dim>::faces_per_cell);
  Assert(face_info.cell_and_face_boundary_id.size(0) >= n_cell_batches(),
         ExcNotInitialized());
  std::array<types::boundary_id, VectorizedArrayType::size()> result;
  result.fill(numbers::invalid_boundary_id);
  for (unsigned int v = 0;
       v < n_active_entries_per_cell_batch(cell_batch_index);
       ++v)
    result[v] =
      face_info.cell_and_face_boundary_id(cell_batch_index, face_number, v);
  return result;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const internal::MatrixFreeFunctions::
  MappingInfo<dim, Number, VectorizedArrayType> &
  MatrixFree<dim, Number, VectorizedArrayType>::get_mapping_info() const
{
  return mapping_info;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const internal::MatrixFreeFunctions::DoFInfo &
MatrixFree<dim, Number, VectorizedArrayType>::get_dof_info(
  const unsigned int dof_index) const
{
  AssertIndexRange(dof_index, n_components());
  return dof_info[dof_index];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_constraint_pool_entries() const
{
  return constraint_pool_row_index.size() - 1;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const Number *
MatrixFree<dim, Number, VectorizedArrayType>::constraint_pool_begin(
  const unsigned int row) const
{
  AssertIndexRange(row, constraint_pool_row_index.size() - 1);
  return constraint_pool_data.empty() ?
           nullptr :
           constraint_pool_data.data() + constraint_pool_row_index[row];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const Number *
MatrixFree<dim, Number, VectorizedArrayType>::constraint_pool_end(
  const unsigned int row) const
{
  AssertIndexRange(row, constraint_pool_row_index.size() - 1);
  return constraint_pool_data.empty() ?
           nullptr :
           constraint_pool_data.data() + constraint_pool_row_index[row + 1];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline std::pair<unsigned int, unsigned int>
MatrixFree<dim, Number, VectorizedArrayType>::create_cell_subrange_hp(
  const std::pair<unsigned int, unsigned int> &range,
  const unsigned int                           degree,
  const unsigned int                           dof_handler_component) const
{
  if (dof_info[dof_handler_component].cell_active_fe_index.empty())
    {
      AssertDimension(
        dof_info[dof_handler_component].fe_index_conversion.size(), 1);
      AssertDimension(
        dof_info[dof_handler_component].fe_index_conversion[0].size(), 1);
      if (dof_info[dof_handler_component].fe_index_conversion[0][0] == degree)
        return range;
      else
        return {range.second, range.second};
    }

  const unsigned int fe_index =
    dof_info[dof_handler_component].fe_index_from_degree(0, degree);
  if (fe_index >= dof_info[dof_handler_component].max_fe_index)
    return {range.second, range.second};
  else
    return create_cell_subrange_hp_by_index(range,
                                            fe_index,
                                            dof_handler_component);
}



template <int dim, typename Number, typename VectorizedArrayType>
inline bool
MatrixFree<dim, Number, VectorizedArrayType>::at_irregular_cell(
  const unsigned int cell_batch_index) const
{
  AssertIndexRange(cell_batch_index, task_info.cell_partition_data.back());
  return VectorizedArrayType::size() > 1 &&
         cell_level_index[(cell_batch_index + 1) * VectorizedArrayType::size() -
                          1] == cell_level_index[(cell_batch_index + 1) *
                                                   VectorizedArrayType::size() -
                                                 2];
}



template <int dim, typename Number, typename VectorizedArrayType>
unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_active_fe_indices() const
{
  return shape_info.size(2);
}


template <int dim, typename Number, typename VectorizedArrayType>
unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_cell_active_fe_index(
  const std::pair<unsigned int, unsigned int> range) const
{
  const auto &fe_indices = dof_info[0].cell_active_fe_index;

  if (fe_indices.empty() == true)
    return 0;

  const auto index = fe_indices[range.first];

  for (unsigned int i = range.first; i < range.second; ++i)
    AssertDimension(index, fe_indices[i]);

  return index;
}



template <int dim, typename Number, typename VectorizedArrayType>
unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_face_active_fe_index(
  const std::pair<unsigned int, unsigned int> range,
  const bool                                  is_interior_face) const
{
  const auto &fe_indices = dof_info[0].cell_active_fe_index;

  if (fe_indices.empty() == true)
    return 0;

  if (is_interior_face)
    {
      const unsigned int index =
        fe_indices[face_info.faces[range.first].cells_interior[0] /
                   VectorizedArrayType::size()];

      for (unsigned int i = range.first; i < range.second; ++i)
        AssertDimension(index,
                        fe_indices[face_info.faces[i].cells_interior[0] /
                                   VectorizedArrayType::size()]);

      return index;
    }
  else
    {
      const unsigned int index =
        fe_indices[face_info.faces[range.first].cells_exterior[0] /
                   VectorizedArrayType::size()];

      for (unsigned int i = range.first; i < range.second; ++i)
        AssertDimension(index,
                        fe_indices[face_info.faces[i].cells_exterior[0] /
                                   VectorizedArrayType::size()]);

      return index;
    }
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_components_filled(
  const unsigned int cell_batch_index) const
{
  return n_active_entries_per_cell_batch(cell_batch_index);
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_active_entries_per_cell_batch(
  const unsigned int cell_batch_index) const
{
  AssertIndexRange(cell_batch_index, task_info.cell_partition_data.back());
  unsigned int n_lanes = VectorizedArrayType::size();
  while (n_lanes > 1 &&
         cell_level_index[cell_batch_index * VectorizedArrayType::size() +
                          n_lanes - 1] ==
           cell_level_index[cell_batch_index * VectorizedArrayType::size() +
                            n_lanes - 2])
    --n_lanes;
  AssertIndexRange(n_lanes - 1, VectorizedArrayType::size());
  return n_lanes;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_active_entries_per_face_batch(
  const unsigned int face_batch_index) const
{
  AssertIndexRange(face_batch_index, face_info.faces.size());
  unsigned int n_lanes = VectorizedArrayType::size();
  while (n_lanes > 1 &&
         face_info.faces[face_batch_index].cells_interior[n_lanes - 1] ==
           numbers::invalid_unsigned_int)
    --n_lanes;
  AssertIndexRange(n_lanes - 1, VectorizedArrayType::size());
  return n_lanes;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_dofs_per_cell(
  const unsigned int dof_handler_index,
  const unsigned int active_fe_index) const
{
  return dof_info[dof_handler_index].dofs_per_cell[active_fe_index];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_n_q_points(
  const unsigned int quad_index,
  const unsigned int active_fe_index) const
{
  AssertIndexRange(quad_index, mapping_info.cell_data.size());
  return mapping_info.cell_data[quad_index]
    .descriptor[active_fe_index]
    .n_q_points;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_dofs_per_face(
  const unsigned int dof_handler_index,
  const unsigned int active_fe_index) const
{
  return dof_info[dof_handler_index].dofs_per_face[active_fe_index];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_n_q_points_face(
  const unsigned int quad_index,
  const unsigned int active_fe_index) const
{
  AssertIndexRange(quad_index, mapping_info.face_data.size());
  return mapping_info.face_data[quad_index]
    .descriptor[active_fe_index]
    .n_q_points;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const IndexSet &
MatrixFree<dim, Number, VectorizedArrayType>::get_locally_owned_set(
  const unsigned int dof_handler_index) const
{
  return dof_info[dof_handler_index].vector_partitioner->locally_owned_range();
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const IndexSet &
MatrixFree<dim, Number, VectorizedArrayType>::get_ghost_set(
  const unsigned int dof_handler_index) const
{
  return dof_info[dof_handler_index].vector_partitioner->ghost_indices();
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const internal::MatrixFreeFunctions::ShapeInfo<VectorizedArrayType> &
MatrixFree<dim, Number, VectorizedArrayType>::get_shape_info(
  const unsigned int dof_handler_index,
  const unsigned int index_quad,
  const unsigned int index_fe,
  const unsigned int active_fe_index,
  const unsigned int active_quad_index) const
{
  AssertIndexRange(dof_handler_index, dof_info.size());
  const unsigned int ind =
    dof_info[dof_handler_index].global_base_element_offset + index_fe;
  AssertIndexRange(ind, shape_info.size(0));
  AssertIndexRange(index_quad, shape_info.size(1));
  AssertIndexRange(active_fe_index, shape_info.size(2));
  AssertIndexRange(active_quad_index, shape_info.size(3));
  return shape_info(ind, index_quad, active_fe_index, active_quad_index);
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const internal::MatrixFreeFunctions::FaceToCellTopology<
  VectorizedArrayType::size()> &
MatrixFree<dim, Number, VectorizedArrayType>::get_face_info(
  const unsigned int macro_face) const
{
  AssertIndexRange(macro_face, face_info.faces.size());
  return face_info.faces[macro_face];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const Table<3, unsigned int> &
MatrixFree<dim, Number, VectorizedArrayType>::get_cell_and_face_to_plain_faces()
  const
{
  return face_info.cell_and_face_to_plain_faces;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const Quadrature<dim> &
MatrixFree<dim, Number, VectorizedArrayType>::get_quadrature(
  const unsigned int quad_index,
  const unsigned int active_fe_index) const
{
  AssertIndexRange(quad_index, mapping_info.cell_data.size());
  return mapping_info.cell_data[quad_index]
    .descriptor[active_fe_index]
    .quadrature;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const Quadrature<dim - 1> &
MatrixFree<dim, Number, VectorizedArrayType>::get_face_quadrature(
  const unsigned int quad_index,
  const unsigned int active_fe_index) const
{
  AssertIndexRange(quad_index, mapping_info.face_data.size());
  return mapping_info.face_data[quad_index]
    .descriptor[active_fe_index]
    .quadrature;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_cell_category(
  const unsigned int cell_batch_index) const
{
  AssertIndexRange(0, dof_info.size());
  AssertIndexRange(cell_batch_index, dof_info[0].cell_active_fe_index.size());
  if (dof_info[0].cell_active_fe_index.empty())
    return 0;
  else
    return dof_info[0].cell_active_fe_index[cell_batch_index];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline std::pair<unsigned int, unsigned int>
MatrixFree<dim, Number, VectorizedArrayType>::get_face_category(
  const unsigned int macro_face) const
{
  AssertIndexRange(macro_face, face_info.faces.size());
  if (dof_info[0].cell_active_fe_index.empty())
    return std::make_pair(0U, 0U);

  std::pair<unsigned int, unsigned int> result;
  for (unsigned int v = 0; v < VectorizedArrayType::size() &&
                           face_info.faces[macro_face].cells_interior[v] !=
                             numbers::invalid_unsigned_int;
       ++v)
    result.first = std::max(
      result.first,
      dof_info[0]
        .cell_active_fe_index[face_info.faces[macro_face].cells_interior[v]]);
  if (face_info.faces[macro_face].cells_exterior[0] !=
      numbers::invalid_unsigned_int)
    for (unsigned int v = 0; v < VectorizedArrayType::size() &&
                             face_info.faces[macro_face].cells_exterior[v] !=
                               numbers::invalid_unsigned_int;
         ++v)
      result.second = std::max(
        result.first,
        dof_info[0]
          .cell_active_fe_index[face_info.faces[macro_face].cells_exterior[v]]);
  else
    result.second = numbers::invalid_unsigned_int;
  return result;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline bool
MatrixFree<dim, Number, VectorizedArrayType>::indices_initialized() const
{
  return indices_are_initialized;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline bool
MatrixFree<dim, Number, VectorizedArrayType>::mapping_initialized() const
{
  return mapping_is_initialized;
}


template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_mg_level() const
{
  return mg_level;
}



template <int dim, typename Number, typename VectorizedArrayType>
AlignedVector<VectorizedArrayType> *
MatrixFree<dim, Number, VectorizedArrayType>::acquire_scratch_data() const
{
  using list_type =
    std::list<std::pair<bool, AlignedVector<VectorizedArrayType>>>;
  list_type &data = scratch_pad.get();
  for (typename list_type::iterator it = data.begin(); it != data.end(); ++it)
    if (it->first == false)
      {
        it->first = true;
        return &it->second;
      }
  data.emplace_front(true, AlignedVector<VectorizedArrayType>());
  return &data.front().second;
}



template <int dim, typename Number, typename VectorizedArrayType>
void
MatrixFree<dim, Number, VectorizedArrayType>::release_scratch_data(
  const AlignedVector<VectorizedArrayType> *scratch) const
{
  using list_type =
    std::list<std::pair<bool, AlignedVector<VectorizedArrayType>>>;
  list_type &data = scratch_pad.get();
  for (typename list_type::iterator it = data.begin(); it != data.end(); ++it)
    if (&it->second == scratch)
      {
        Assert(it->first == true, ExcInternalError());
        it->first = false;
        return;
      }
  AssertThrow(false, ExcMessage("Tried to release invalid scratch pad"));
}



template <int dim, typename Number, typename VectorizedArrayType>
AlignedVector<Number> *
MatrixFree<dim, Number, VectorizedArrayType>::
  acquire_scratch_data_non_threadsafe() const
{
  for (typename std::list<std::pair<bool, AlignedVector<Number>>>::iterator it =
         scratch_pad_non_threadsafe.begin();
       it != scratch_pad_non_threadsafe.end();
       ++it)
    if (it->first == false)
      {
        it->first = true;
        return &it->second;
      }
  scratch_pad_non_threadsafe.push_front(
    std::make_pair(true, AlignedVector<Number>()));
  return &scratch_pad_non_threadsafe.front().second;
}



template <int dim, typename Number, typename VectorizedArrayType>
void
MatrixFree<dim, Number, VectorizedArrayType>::
  release_scratch_data_non_threadsafe(
    const AlignedVector<Number> *scratch) const
{
  for (typename std::list<std::pair<bool, AlignedVector<Number>>>::iterator it =
         scratch_pad_non_threadsafe.begin();
       it != scratch_pad_non_threadsafe.end();
       ++it)
    if (&it->second == scratch)
      {
        Assert(it->first == true, ExcInternalError());
        it->first = false;
        return;
      }
  AssertThrow(false, ExcMessage("Tried to release invalid scratch pad"));
}



// ------------------------------ reinit functions ---------------------------

namespace internal
{
  namespace MatrixFreeImplementation
  {
    template <int dim, int spacedim>
    inline std::vector<IndexSet>
    extract_locally_owned_index_sets(
      const std::vector<const ::dealii::DoFHandler<dim, spacedim> *> &dofh,
      const unsigned int                                              level)
    {
      std::vector<IndexSet> locally_owned_set;
      locally_owned_set.reserve(dofh.size());
      for (unsigned int j = 0; j < dofh.size(); j++)
        if (level == numbers::invalid_unsigned_int)
          locally_owned_set.push_back(dofh[j]->locally_owned_dofs());
        else
          locally_owned_set.push_back(dofh[j]->locally_owned_mg_dofs(level));
      return locally_owned_set;
    }
  } // namespace MatrixFreeImplementation
} // namespace internal



template <int dim, typename Number, typename VectorizedArrayType>
template <typename QuadratureType, typename number2>
void
MatrixFree<dim, Number, VectorizedArrayType>::reinit(
  const DoFHandler<dim> &           dof_handler,
  const AffineConstraints<number2> &constraints_in,
  const QuadratureType &            quad,
  const typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    &additional_data)
{
  std::vector<const DoFHandler<dim, dim> *>       dof_handlers;
  std::vector<const AffineConstraints<number2> *> constraints;
  std::vector<QuadratureType>                     quads;

  dof_handlers.push_back(&dof_handler);
  constraints.push_back(&constraints_in);
  quads.push_back(quad);

  std::vector<IndexSet> locally_owned_sets =
    internal::MatrixFreeImplementation::extract_locally_owned_index_sets(
      dof_handlers, additional_data.mg_level);

  std::vector<hp::QCollection<dim>> quad_hp;
  quad_hp.emplace_back(quad);

  internal_reinit(std::make_shared<hp::MappingCollection<dim>>(
                    StaticMappingQ1<dim>::mapping),
                  dof_handlers,
                  constraints,
                  locally_owned_sets,
                  quad_hp,
                  additional_data);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename QuadratureType, typename number2, typename MappingType>
void
MatrixFree<dim, Number, VectorizedArrayType>::reinit(
  const MappingType &               mapping,
  const DoFHandler<dim> &           dof_handler,
  const AffineConstraints<number2> &constraints_in,
  const QuadratureType &            quad,
  const typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    &additional_data)
{
  std::vector<const DoFHandler<dim, dim> *>       dof_handlers;
  std::vector<const AffineConstraints<number2> *> constraints;

  dof_handlers.push_back(&dof_handler);
  constraints.push_back(&constraints_in);

  std::vector<IndexSet> locally_owned_sets =
    internal::MatrixFreeImplementation::extract_locally_owned_index_sets(
      dof_handlers, additional_data.mg_level);

  std::vector<hp::QCollection<dim>> quad_hp;
  quad_hp.emplace_back(quad);

  internal_reinit(std::make_shared<hp::MappingCollection<dim>>(mapping),
                  dof_handlers,
                  constraints,
                  locally_owned_sets,
                  quad_hp,
                  additional_data);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename QuadratureType, typename number2>
void
MatrixFree<dim, Number, VectorizedArrayType>::reinit(
  const std::vector<const DoFHandler<dim> *> &           dof_handler,
  const std::vector<const AffineConstraints<number2> *> &constraint,
  const std::vector<QuadratureType> &                    quad,
  const typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    &additional_data)
{
  std::vector<IndexSet> locally_owned_set =
    internal::MatrixFreeImplementation::extract_locally_owned_index_sets(
      dof_handler, additional_data.mg_level);
  std::vector<hp::QCollection<dim>> quad_hp;
  for (unsigned int q = 0; q < quad.size(); ++q)
    quad_hp.emplace_back(quad[q]);

  internal_reinit(std::make_shared<hp::MappingCollection<dim>>(
                    StaticMappingQ1<dim>::mapping),
                  dof_handler,
                  constraint,
                  locally_owned_set,
                  quad_hp,
                  additional_data);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename QuadratureType,
          typename number2,
          typename DoFHandlerType,
          typename MappingType>
void
MatrixFree<dim, Number, VectorizedArrayType>::reinit(
  const MappingType &                                    mapping,
  const std::vector<const DoFHandlerType *> &            dof_handler,
  const std::vector<const AffineConstraints<number2> *> &constraint,
  const std::vector<QuadratureType> &                    quad,
  const AdditionalData &                                 additional_data)
{
  static_assert(dim == DoFHandlerType::dimension,
                "Dimension dim not equal to DoFHandlerType::dimension.");

  std::vector<const DoFHandler<dim> *> dof_handlers;

  for (const auto dh : dof_handler)
    dof_handlers.push_back(dh);

  this->reinit(mapping, dof_handlers, constraint, quad, additional_data);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename QuadratureType, typename number2>
void
MatrixFree<dim, Number, VectorizedArrayType>::reinit(
  const std::vector<const DoFHandler<dim> *> &           dof_handler,
  const std::vector<const AffineConstraints<number2> *> &constraint,
  const QuadratureType &                                 quad,
  const typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    &additional_data)
{
  std::vector<IndexSet> locally_owned_set =
    internal::MatrixFreeImplementation::extract_locally_owned_index_sets(
      dof_handler, additional_data.mg_level);
  std::vector<hp::QCollection<dim>> quad_hp;
  quad_hp.emplace_back(quad);

  internal_reinit(std::make_shared<hp::MappingCollection<dim>>(
                    StaticMappingQ1<dim>::mapping),
                  dof_handler,
                  constraint,
                  locally_owned_set,
                  quad_hp,
                  additional_data);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename QuadratureType, typename number2, typename DoFHandlerType>
void
MatrixFree<dim, Number, VectorizedArrayType>::reinit(
  const std::vector<const DoFHandlerType *> &            dof_handler,
  const std::vector<const AffineConstraints<number2> *> &constraint,
  const std::vector<QuadratureType> &                    quad,
  const AdditionalData &                                 additional_data)
{
  static_assert(dim == DoFHandlerType::dimension,
                "Dimension dim not equal to DoFHandlerType::dimension.");

  std::vector<const DoFHandler<dim> *> dof_handlers;

  for (const auto dh : dof_handler)
    dof_handlers.push_back(dof_handler);

  this->reinit(dof_handlers, constraint, quad, additional_data);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename QuadratureType, typename number2, typename MappingType>
void
MatrixFree<dim, Number, VectorizedArrayType>::reinit(
  const MappingType &                                    mapping,
  const std::vector<const DoFHandler<dim> *> &           dof_handler,
  const std::vector<const AffineConstraints<number2> *> &constraint,
  const QuadratureType &                                 quad,
  const typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    &additional_data)
{
  std::vector<IndexSet> locally_owned_set =
    internal::MatrixFreeImplementation::extract_locally_owned_index_sets(
      dof_handler, additional_data.mg_level);
  std::vector<hp::QCollection<dim>> quad_hp;
  quad_hp.emplace_back(quad);

  internal_reinit(std::make_shared<hp::MappingCollection<dim>>(mapping),
                  dof_handler,
                  constraint,
                  locally_owned_set,
                  quad_hp,
                  additional_data);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename QuadratureType, typename number2, typename MappingType>
void
MatrixFree<dim, Number, VectorizedArrayType>::reinit(
  const MappingType &                                    mapping,
  const std::vector<const DoFHandler<dim> *> &           dof_handler,
  const std::vector<const AffineConstraints<number2> *> &constraint,
  const std::vector<QuadratureType> &                    quad,
  const typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    &additional_data)
{
  std::vector<IndexSet> locally_owned_set =
    internal::MatrixFreeImplementation::extract_locally_owned_index_sets(
      dof_handler, additional_data.mg_level);
  std::vector<hp::QCollection<dim>> quad_hp;
  for (unsigned int q = 0; q < quad.size(); ++q)
    quad_hp.emplace_back(quad[q]);

  internal_reinit(std::make_shared<hp::MappingCollection<dim>>(mapping),
                  dof_handler,
                  constraint,
                  locally_owned_set,
                  quad_hp,
                  additional_data);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename QuadratureType,
          typename number2,
          typename DoFHandlerType,
          typename MappingType>
void
MatrixFree<dim, Number, VectorizedArrayType>::reinit(
  const MappingType &                                    mapping,
  const std::vector<const DoFHandlerType *> &            dof_handler,
  const std::vector<const AffineConstraints<number2> *> &constraint,
  const QuadratureType &                                 quad,
  const AdditionalData &                                 additional_data)
{
  static_assert(dim == DoFHandlerType::dimension,
                "Dimension dim not equal to DoFHandlerType::dimension.");

  std::vector<const DoFHandler<dim> *> dof_handlers;

  for (const auto dh : dof_handler)
    dof_handlers.push_back(dof_handler);

  this->reinit(mapping, dof_handlers, constraint, quad, additional_data);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename QuadratureType, typename number2, typename DoFHandlerType>
void
MatrixFree<dim, Number, VectorizedArrayType>::reinit(
  const std::vector<const DoFHandlerType *> &            dof_handler,
  const std::vector<const AffineConstraints<number2> *> &constraint,
  const QuadratureType &                                 quad,
  const AdditionalData &                                 additional_data)
{
  static_assert(dim == DoFHandlerType::dimension,
                "Dimension dim not equal to DoFHandlerType::dimension.");

  std::vector<const DoFHandler<dim> *> dof_handlers;

  for (const auto dh : dof_handler)
    dof_handlers.push_back(dof_handler);

  this->reinit(dof_handlers, constraint, quad, additional_data);
}



// ------------------------------ implementation of loops --------------------

// internal helper functions that define how to call MPI data exchange
// functions: for generic vectors, do nothing at all. For distributed vectors,
// call update_ghost_values_start function and so on. If we have collections
// of vectors, just do the individual functions of the components. In order to
// keep ghost values consistent (whether we are in read or write mode), we
// also reset the values at the end. the whole situation is a bit complicated
// by the fact that we need to treat block vectors differently, which use some
// additional helper functions to select the blocks and template magic.
namespace internal
{
  /**
   * 用于交换向量间数据的内部类。
   *
   */
  template <int dim, typename Number, typename VectorizedArrayType>
  struct VectorDataExchange
  {
    // A shift for the MPI messages to reduce the risk for accidental
    // interaction with other open communications that a user program might
    // set up (parallel vectors support unfinished communication). We let
    // the other vectors use the first 20 assigned numbers and start the
    // matrix-free communication.
    static constexpr unsigned int channel_shift = 20;



    /**
     * 构造函数。接受MF数据，DG中面部访问的标志和组件的数量。
     *
     */
    VectorDataExchange(
      const dealii::MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
      const typename dealii::MatrixFree<dim, Number, VectorizedArrayType>::
        DataAccessOnFaces vector_face_access,
      const unsigned int  n_components)
      : matrix_free(matrix_free)
      , vector_face_access(
          matrix_free.get_task_info().face_partition_data.empty() ?
            dealii::MatrixFree<dim, Number, VectorizedArrayType>::
              DataAccessOnFaces::unspecified :
            vector_face_access)
      , ghosts_were_set(false)
#  ifdef DEAL_II_WITH_MPI
      , tmp_data(n_components)
      , requests(n_components)
#  endif
    {
      (void)n_components;
      if (this->vector_face_access !=
          dealii::MatrixFree<dim, Number, VectorizedArrayType>::
            DataAccessOnFaces::unspecified)
        for (unsigned int c = 0; c < matrix_free.n_components(); ++c)
          AssertDimension(
            matrix_free.get_dof_info(c).vector_exchanger_face_variants.size(),
            5);
    }



    /**
     * 解构器。
     *
     */
    ~VectorDataExchange() // NOLINT
    {
#  ifdef DEAL_II_WITH_MPI
      for (unsigned int i = 0; i < tmp_data.size(); ++i)
        if (tmp_data[i] != nullptr)
          matrix_free.release_scratch_data_non_threadsafe(tmp_data[i]);
#  endif
    }



    /**
     * 遍历MF对象中的所有组件，选择其分区器与该组件中的分区器兼容的组件。
     *
     */
    template <typename VectorType>
    unsigned int
    find_vector_in_mf(const VectorType &vec,
                      const bool        check_global_compatibility = true) const
    {
      // case 1: vector was set up with MatrixFree::initialize_dof_vector()
      for (unsigned int c = 0; c < matrix_free.n_components(); ++c)
        if (vec.get_partitioner().get() ==
            matrix_free.get_dof_info(c).vector_partitioner.get())
          return c;

      // case 2: user provided own partitioner (compatibility mode)
      for (unsigned int c = 0; c < matrix_free.n_components(); ++c)
        if (check_global_compatibility ?
              vec.get_partitioner()->is_globally_compatible(
                *matrix_free.get_dof_info(c).vector_partitioner) :
              vec.get_partitioner()->is_compatible(
                *matrix_free.get_dof_info(c).vector_partitioner))
          return c;

      Assert(false,
             ExcNotImplemented("Could not find partitioner that fits vector"));

      return numbers::invalid_unsigned_int;
    }



    /**
     * 为给定的 @p mf_component
     * 获取分区器，同时考虑到构造函数中设置的vector_face_access。
     *
     */
    const internal::MatrixFreeFunctions::VectorDataExchange::Base &
    get_partitioner(const unsigned int mf_component) const
    {
      AssertDimension(matrix_free.get_dof_info(mf_component)
                        .vector_exchanger_face_variants.size(),
                      5);
      if (vector_face_access ==
          dealii::MatrixFree<dim, Number, VectorizedArrayType>::
            DataAccessOnFaces::none)
        return *matrix_free.get_dof_info(mf_component)
                  .vector_exchanger_face_variants[0];
      else if (vector_face_access ==
               dealii::MatrixFree<dim, Number, VectorizedArrayType>::
                 DataAccessOnFaces::values)
        return *matrix_free.get_dof_info(mf_component)
                  .vector_exchanger_face_variants[1];
      else if (vector_face_access ==
               dealii::MatrixFree<dim, Number, VectorizedArrayType>::
                 DataAccessOnFaces::gradients)
        return *matrix_free.get_dof_info(mf_component)
                  .vector_exchanger_face_variants[2];
      else if (vector_face_access ==
               dealii::MatrixFree<dim, Number, VectorizedArrayType>::
                 DataAccessOnFaces::values_all_faces)
        return *matrix_free.get_dof_info(mf_component)
                  .vector_exchanger_face_variants[3];
      else if (vector_face_access ==
               dealii::MatrixFree<dim, Number, VectorizedArrayType>::
                 DataAccessOnFaces::gradients_all_faces)
        return *matrix_free.get_dof_info(mf_component)
                  .vector_exchanger_face_variants[4];
      else
        return *matrix_free.get_dof_info(mf_component).vector_exchanger.get();
    }



    /**
     * 开始为序列向量更新_ghost_value
     *
     */
    template <typename VectorType,
              typename std::enable_if<is_serial_or_dummy<VectorType>::value,
                                      VectorType>::type * = nullptr>
    void
    update_ghost_values_start(const unsigned int  /*component_in_block_vector*/ ,
                              const VectorType &  /*vec*/ )
    {}


    /**
     * 对于不支持分割成_start()和finish()阶段的向量，开始更新_ghost_value。
     *
     */
    template <typename VectorType,
              typename std::enable_if<
                !has_update_ghost_values_start<VectorType>::value &&
                  !is_serial_or_dummy<VectorType>::value,
                VectorType>::type * = nullptr>
    void
    update_ghost_values_start(const unsigned int component_in_block_vector,
                              const VectorType & vec)
    {
      (void)component_in_block_vector;
      bool ghosts_set = vec.has_ghost_elements();
      if (ghosts_set)
        ghosts_were_set = true;

      vec.update_ghost_values();
    }



    /**
     * 对于那些支持分成_start()和finish()阶段，但不支持在DoF子集上交换的向量，开始更新_ghost_value。
     *
     */
    template <typename VectorType,
              typename std::enable_if<
                has_update_ghost_values_start<VectorType>::value &&
                  !has_exchange_on_subset<VectorType>::value,
                VectorType>::type * = nullptr>
    void
    update_ghost_values_start(const unsigned int component_in_block_vector,
                              const VectorType & vec)
    {
      (void)component_in_block_vector;
      bool ghosts_set = vec.has_ghost_elements();
      if (ghosts_set)
        ghosts_were_set = true;

      vec.update_ghost_values_start(component_in_block_vector + channel_shift);
    }



    /**
     * 最后，对于那些支持分成_start()和finish()阶段，并且支持在子集DoF上交换的向量，即
     * LinearAlgebra::distributed::Vector ，开始更新_ghost_value。
     *
     */
    template <typename VectorType,
              typename std::enable_if<
                has_update_ghost_values_start<VectorType>::value &&
                  has_exchange_on_subset<VectorType>::value,
                VectorType>::type * = nullptr>
    void
    update_ghost_values_start(const unsigned int component_in_block_vector,
                              const VectorType & vec)
    {
      static_assert(
        std::is_same<Number, typename VectorType::value_type>::value,
        "Type mismatch between VectorType and VectorDataExchange");
      (void)component_in_block_vector;
      bool ghosts_set = vec.has_ghost_elements();
      if (ghosts_set)
        ghosts_were_set = true;

      if (vec.size() != 0)
        {
#  ifdef DEAL_II_WITH_MPI
          const unsigned int mf_component = find_vector_in_mf(vec);

          const auto &part = get_partitioner(mf_component);

          if (part.n_ghost_indices() == 0 && part.n_import_indices() == 0 &&
              part.n_import_sm_procs() == 0)
            return;

          tmp_data[component_in_block_vector] =
            matrix_free.acquire_scratch_data_non_threadsafe();
          tmp_data[component_in_block_vector]->resize_fast(
            part.n_import_indices());
          AssertDimension(requests.size(), tmp_data.size());

          part.export_to_ghosted_array_start(
            component_in_block_vector * 2 + channel_shift,
            ArrayView<const Number>(vec.begin(), part.locally_owned_size()),
            vec.shared_vector_data(),
            ArrayView<Number>(const_cast<Number *>(vec.begin()) +
                                part.locally_owned_size(),
                              matrix_free.get_dof_info(mf_component)
                                .vector_partitioner->n_ghost_indices()),
            ArrayView<Number>(tmp_data[component_in_block_vector]->begin(),
                              part.n_import_indices()),
            this->requests[component_in_block_vector]);
#  endif
        }
    }



    /**
     * 对于不支持分成_start()和finish()阶段的向量和串行向量，完成update_ghost_value。
     *
     */
    template <
      typename VectorType,
      typename std::enable_if<!has_update_ghost_values_start<VectorType>::value,
                              VectorType>::type * = nullptr>
    void
    update_ghost_values_finish(const unsigned int  /*component_in_block_vector*/ ,
                               const VectorType &  /*vec*/ )
    {}



    /**
     * 完成对支持分割成_start()和finish()阶段的向量的
     * update_ghost_value，但不支持在DoF的子集上进行交换。
     *
     */
    template <typename VectorType,
              typename std::enable_if<
                has_update_ghost_values_start<VectorType>::value &&
                  !has_exchange_on_subset<VectorType>::value,
                VectorType>::type * = nullptr>
    void
    update_ghost_values_finish(const unsigned int component_in_block_vector,
                               const VectorType & vec)
    {
      (void)component_in_block_vector;
      vec.update_ghost_values_finish();
    }



    /**
     * 完成update_ghost_value，用于_支持分成_start()和finish()阶段的向量，也支持在子集DoF上的交换，即
     * LinearAlgebra::distributed::Vector
     *
     */
    template <typename VectorType,
              typename std::enable_if<
                has_update_ghost_values_start<VectorType>::value &&
                  has_exchange_on_subset<VectorType>::value,
                VectorType>::type * = nullptr>
    void
    update_ghost_values_finish(const unsigned int component_in_block_vector,
                               const VectorType & vec)
    {
      static_assert(
        std::is_same<Number, typename VectorType::value_type>::value,
        "Type mismatch between VectorType and VectorDataExchange");
      (void)component_in_block_vector;

      if (vec.size() != 0)
        {
#  ifdef DEAL_II_WITH_MPI
          AssertIndexRange(component_in_block_vector, tmp_data.size());
          AssertDimension(requests.size(), tmp_data.size());

          const unsigned int mf_component = find_vector_in_mf(vec);

          const auto &part = get_partitioner(mf_component);

          if (part.n_ghost_indices() != 0 || part.n_import_indices() != 0 ||
              part.n_import_sm_procs() != 0)
            {
              part.export_to_ghosted_array_finish(
                ArrayView<const Number>(vec.begin(), part.locally_owned_size()),
                vec.shared_vector_data(),
                ArrayView<Number>(const_cast<Number *>(vec.begin()) +
                                    part.locally_owned_size(),
                                  matrix_free.get_dof_info(mf_component)
                                    .vector_partitioner->n_ghost_indices()),
                this->requests[component_in_block_vector]);

              matrix_free.release_scratch_data_non_threadsafe(
                tmp_data[component_in_block_vector]);
              tmp_data[component_in_block_vector] = nullptr;
            }
#  endif
        }
      // let vector know that ghosts are being updated and we can read from
      // them
      vec.set_ghost_state(true);
    }



    /**
     * 对串行向量进行启动压缩
     *
     */
    template <typename VectorType,
              typename std::enable_if<is_serial_or_dummy<VectorType>::value,
                                      VectorType>::type * = nullptr>
    void
    compress_start(const unsigned int  /*component_in_block_vector*/ ,
                   VectorType &  /*vec*/ )
    {}



    /**
     * 对于不支持分割成_start()和finish()阶段的向量，开始压缩
     *
     */
    template <typename VectorType,
              typename std::enable_if<!has_compress_start<VectorType>::value &&
                                        !is_serial_or_dummy<VectorType>::value,
                                      VectorType>::type * = nullptr>
    void
    compress_start(const unsigned int component_in_block_vector,
                   VectorType &       vec)
    {
      (void)component_in_block_vector;
      Assert(vec.has_ghost_elements() == false, ExcNotImplemented());
      vec.compress(dealii::VectorOperation::add);
    }



    /**
     * 对于支持分割成_start()和finish()阶段的向量开始压缩，但不支持在DoF的一个子集上交换。
     *
     */
    template <
      typename VectorType,
      typename std::enable_if<has_compress_start<VectorType>::value &&
                                !has_exchange_on_subset<VectorType>::value,
                              VectorType>::type * = nullptr>
    void
    compress_start(const unsigned int component_in_block_vector,
                   VectorType &       vec)
    {
      (void)component_in_block_vector;
      Assert(vec.has_ghost_elements() == false, ExcNotImplemented());
      vec.compress_start(component_in_block_vector + channel_shift);
    }



    /**
     * 开始压缩那些_支持分割成_start()和finish()阶段的向量，并且也支持在子集DoF上的交换，即
     * LinearAlgebra::distributed::Vector  。
     *
     */
    template <
      typename VectorType,
      typename std::enable_if<has_compress_start<VectorType>::value &&
                                has_exchange_on_subset<VectorType>::value,
                              VectorType>::type * = nullptr>
    void
    compress_start(const unsigned int component_in_block_vector,
                   VectorType &       vec)
    {
      static_assert(
        std::is_same<Number, typename VectorType::value_type>::value,
        "Type mismatch between VectorType and VectorDataExchange");
      (void)component_in_block_vector;
      Assert(vec.has_ghost_elements() == false, ExcNotImplemented());

      if (vec.size() != 0)
        {
#  ifdef DEAL_II_WITH_MPI
          const unsigned int mf_component = find_vector_in_mf(vec);

          const auto &part = get_partitioner(mf_component);

          if (part.n_ghost_indices() == 0 && part.n_import_indices() == 0 &&
              part.n_import_sm_procs() == 0)
            return;

          tmp_data[component_in_block_vector] =
            matrix_free.acquire_scratch_data_non_threadsafe();
          tmp_data[component_in_block_vector]->resize_fast(
            part.n_import_indices());
          AssertDimension(requests.size(), tmp_data.size());

          part.import_from_ghosted_array_start(
            dealii::VectorOperation::add,
            component_in_block_vector * 2 + channel_shift,
            ArrayView<Number>(vec.begin(), part.locally_owned_size()),
            vec.shared_vector_data(),
            ArrayView<Number>(vec.begin() + part.locally_owned_size(),
                              matrix_free.get_dof_info(mf_component)
                                .vector_partitioner->n_ghost_indices()),
            ArrayView<Number>(tmp_data[component_in_block_vector]->begin(),
                              part.n_import_indices()),
            this->requests[component_in_block_vector]);
#  endif
        }
    }



    /**
     * 对不支持分成_start()和finish()阶段的向量和串行向量进行完成压缩
     *
     */
    template <typename VectorType,
              typename std::enable_if<!has_compress_start<VectorType>::value,
                                      VectorType>::type * = nullptr>
    void
    compress_finish(const unsigned int  /*component_in_block_vector*/ ,
                    VectorType &  /*vec*/ )
    {}



    /**
     * 对于支持分割成_start()和finish()阶段的向量完成压缩，但不支持在DoF的子集上交换。
     *
     */
    template <
      typename VectorType,
      typename std::enable_if<has_compress_start<VectorType>::value &&
                                !has_exchange_on_subset<VectorType>::value,
                              VectorType>::type * = nullptr>
    void
    compress_finish(const unsigned int component_in_block_vector,
                    VectorType &       vec)
    {
      (void)component_in_block_vector;
      vec.compress_finish(dealii::VectorOperation::add);
    }



    /**
     * 开始压缩那些_支持分割成_start()和finish()阶段的向量，并且也支持在子集DoF上的交换，即
     * LinearAlgebra::distributed::Vector  。
     *
     */
    template <
      typename VectorType,
      typename std::enable_if<has_compress_start<VectorType>::value &&
                                has_exchange_on_subset<VectorType>::value,
                              VectorType>::type * = nullptr>
    void
    compress_finish(const unsigned int component_in_block_vector,
                    VectorType &       vec)
    {
      static_assert(
        std::is_same<Number, typename VectorType::value_type>::value,
        "Type mismatch between VectorType and VectorDataExchange");
      (void)component_in_block_vector;
      if (vec.size() != 0)
        {
#  ifdef DEAL_II_WITH_MPI
          AssertIndexRange(component_in_block_vector, tmp_data.size());
          AssertDimension(requests.size(), tmp_data.size());

          const unsigned int mf_component = find_vector_in_mf(vec);

          const auto &part = get_partitioner(mf_component);

          if (part.n_ghost_indices() != 0 || part.n_import_indices() != 0 ||
              part.n_import_sm_procs() != 0)
            {
              part.import_from_ghosted_array_finish(
                VectorOperation::add,
                ArrayView<Number>(vec.begin(), part.locally_owned_size()),
                vec.shared_vector_data(),
                ArrayView<Number>(vec.begin() + part.locally_owned_size(),
                                  matrix_free.get_dof_info(mf_component)
                                    .vector_partitioner->n_ghost_indices()),
                ArrayView<const Number>(
                  tmp_data[component_in_block_vector]->begin(),
                  part.n_import_indices()),
                this->requests[component_in_block_vector]);

              matrix_free.release_scratch_data_non_threadsafe(
                tmp_data[component_in_block_vector]);
              tmp_data[component_in_block_vector] = nullptr;
            }

          if (Utilities::MPI::job_supports_mpi())
            {
              const int ierr =
                MPI_Barrier(matrix_free.get_task_info().communicator_sm);
              AssertThrowMPI(ierr);
            }
#  endif
        }
    }



    /**
     * 重置串行向量的所有幽灵值
     *
     */
    template <typename VectorType,
              typename std::enable_if<is_serial_or_dummy<VectorType>::value,
                                      VectorType>::type * = nullptr>
    void
    reset_ghost_values(const VectorType &  /*vec*/ ) const
    {}



    /**
     * 重置不支持在子集DoF上交换的向量的所有鬼值
     *
     */
    template <
      typename VectorType,
      typename std::enable_if<!has_exchange_on_subset<VectorType>::value &&
                                !is_serial_or_dummy<VectorType>::value,
                              VectorType>::type * = nullptr>
    void
    reset_ghost_values(const VectorType &vec) const
    {
      if (ghosts_were_set == true)
        return;

      vec.zero_out_ghost_values();
    }



    /**
     * 重置_支持在子集DoF上交换的向量的所有鬼魂值，即
     * LinearAlgebra::distributed::Vector 。
     *
     */
    template <typename VectorType,
              typename std::enable_if<has_exchange_on_subset<VectorType>::value,
                                      VectorType>::type * = nullptr>
    void
    reset_ghost_values(const VectorType &vec) const
    {
      static_assert(
        std::is_same<Number, typename VectorType::value_type>::value,
        "Type mismatch between VectorType and VectorDataExchange");
      if (ghosts_were_set == true)
        return;

      if (vec.size() != 0)
        {
#  ifdef DEAL_II_WITH_MPI
          AssertDimension(requests.size(), tmp_data.size());

          const unsigned int mf_component = find_vector_in_mf(vec);

          const auto &part = get_partitioner(mf_component);

          if (part.n_ghost_indices() > 0)
            {
              part.reset_ghost_values(ArrayView<Number>(
                const_cast<LinearAlgebra::distributed::Vector<Number> &>(vec)
                    .begin() +
                  part.locally_owned_size(),
                matrix_free.get_dof_info(mf_component)
                  .vector_partitioner->n_ghost_indices()));
            }

#  endif
        }
      // let vector know that it's not ghosted anymore
      vec.set_ghost_state(false);
    }



    /**
     * 对于_do_ support exchange on a subset of DoFs <==> begin() + ind ==
     * local_element(ind), i.e.  LinearAlgebra::distributed::Vector
     * 的向量区域清零。
     *
     */
    template <typename VectorType,
              typename std::enable_if<has_exchange_on_subset<VectorType>::value,
                                      VectorType>::type * = nullptr>
    void
    zero_vector_region(const unsigned int range_index, VectorType &vec) const
    {
      static_assert(
        std::is_same<Number, typename VectorType::value_type>::value,
        "Type mismatch between VectorType and VectorDataExchange");
      if (range_index == numbers::invalid_unsigned_int)
        vec = Number();
      else
        {
          const unsigned int mf_component = find_vector_in_mf(vec, false);
          const internal::MatrixFreeFunctions::DoFInfo &dof_info =
            matrix_free.get_dof_info(mf_component);
          Assert(dof_info.vector_zero_range_list_index.empty() == false,
                 ExcNotInitialized());

          Assert(vec.partitioners_are_compatible(*dof_info.vector_partitioner),
                 ExcInternalError());
          AssertIndexRange(range_index,
                           dof_info.vector_zero_range_list_index.size() - 1);
          for (unsigned int id =
                 dof_info.vector_zero_range_list_index[range_index];
               id != dof_info.vector_zero_range_list_index[range_index + 1];
               ++id)
            std::memset(vec.begin() + dof_info.vector_zero_range_list[id].first,
                        0,
                        (dof_info.vector_zero_range_list[id].second -
                         dof_info.vector_zero_range_list[id].first) *
                          sizeof(Number));
        }
    }



    /**
     * 对于不支持在DoFs子集上进行交换的向量，将向量区域清零
     * <==> begin() + ind == local_element(ind) 但仍然是向量类型。
     *
     */
    template <
      typename VectorType,
      typename std::enable_if<!has_exchange_on_subset<VectorType>::value,
                              VectorType>::type * = nullptr,
      typename VectorType::value_type *           = nullptr>
    void
    zero_vector_region(const unsigned int range_index, VectorType &vec) const
    {
      if (range_index == numbers::invalid_unsigned_int || range_index == 0)
        vec = typename VectorType::value_type();
    }



    /**
     * 对于非向量类型，即没有 VectorType::value_type
     * 的类，将向量区域清零。
     *
     */
    void
    zero_vector_region(const unsigned int, ...) const
    {
      Assert(false,
             ExcNotImplemented("Zeroing is only implemented for vector types "
                               "which provide VectorType::value_type"));
    }



    const dealii::MatrixFree<dim, Number, VectorizedArrayType> &matrix_free;
    const typename dealii::MatrixFree<dim, Number, VectorizedArrayType>::
      DataAccessOnFaces vector_face_access;
    bool                ghosts_were_set;
#  ifdef DEAL_II_WITH_MPI
    std::vector<AlignedVector<Number> *>  tmp_data;
    std::vector<std::vector<MPI_Request>> requests;
#  endif
  }; // VectorDataExchange

  template <typename VectorStruct>
  unsigned int
  n_components(const VectorStruct &vec);

  template <typename VectorStruct>
  unsigned int
  n_components_block(const VectorStruct &vec,
                     std::integral_constant<bool, true>)
  {
    unsigned int components = 0;
    for (unsigned int bl = 0; bl < vec.n_blocks(); ++bl)
      components += n_components(vec.block(bl));
    return components;
  }

  template <typename VectorStruct>
  unsigned int
  n_components_block(const VectorStruct &, std::integral_constant<bool, false>)
  {
    return 1;
  }

  template <typename VectorStruct>
  unsigned int
  n_components(const VectorStruct &vec)
  {
    return n_components_block(
      vec, std::integral_constant<bool, IsBlockVector<VectorStruct>::value>());
  }

  template <typename VectorStruct>
  inline unsigned int
  n_components(const std::vector<VectorStruct> &vec)
  {
    unsigned int components = 0;
    for (unsigned int comp = 0; comp < vec.size(); comp++)
      components += n_components_block(
        vec[comp],
        std::integral_constant<bool, IsBlockVector<VectorStruct>::value>());
    return components;
  }

  template <typename VectorStruct>
  inline unsigned int
  n_components(const std::vector<VectorStruct *> &vec)
  {
    unsigned int components = 0;
    for (unsigned int comp = 0; comp < vec.size(); comp++)
      components += n_components_block(
        *vec[comp],
        std::integral_constant<bool, IsBlockVector<VectorStruct>::value>());
    return components;
  }



  // A helper function to identify block vectors with many components where we
  // should not try to overlap computations and communication because there
  // would be too many outstanding communication requests.

  // default value for vectors that do not have communication_block_size
  template <
    typename VectorStruct,
    typename std::enable_if<!has_communication_block_size<VectorStruct>::value,
                            VectorStruct>::type * = nullptr>
  constexpr unsigned int
  get_communication_block_size(const VectorStruct &)
  {
    return numbers::invalid_unsigned_int;
  }



  template <
    typename VectorStruct,
    typename std::enable_if<has_communication_block_size<VectorStruct>::value,
                            VectorStruct>::type * = nullptr>
  constexpr unsigned int
  get_communication_block_size(const VectorStruct &)
  {
    return VectorStruct::communication_block_size;
  }



  // --------------------------------------------------------------------------
  // below we have wrappers to distinguish between block and non-block vectors.
  // --------------------------------------------------------------------------

  //
  // update_ghost_values_start
  //

  // update_ghost_values for block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            typename std::enable_if<IsBlockVector<VectorStruct>::value,
                                    VectorStruct>::type * = nullptr>
  void
  update_ghost_values_start(
    const VectorStruct &                                  vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    if (get_communication_block_size(vec) < vec.n_blocks())
      {
        // don't forget to set ghosts_were_set, that otherwise happens
        // inside VectorDataExchange::update_ghost_values_start()
        exchanger.ghosts_were_set = vec.has_ghost_elements();
        vec.update_ghost_values();
      }
    else
      {
        for (unsigned int i = 0; i < vec.n_blocks(); ++i)
          update_ghost_values_start(vec.block(i), exchanger, channel + i);
      }
  }



  // update_ghost_values for non-block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            typename std::enable_if<!IsBlockVector<VectorStruct>::value,
                                    VectorStruct>::type * = nullptr>
  void
  update_ghost_values_start(
    const VectorStruct &                                  vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    exchanger.update_ghost_values_start(channel, vec);
  }



  // update_ghost_values_start() for vector of vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  update_ghost_values_start(
    const std::vector<VectorStruct> &                     vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); comp++)
      {
        update_ghost_values_start(vec[comp], exchanger, component_index);
        component_index += n_components(vec[comp]);
      }
  }



  // update_ghost_values_start() for vector of pointers to vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  update_ghost_values_start(
    const std::vector<VectorStruct *> &                   vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); comp++)
      {
        update_ghost_values_start(*vec[comp], exchanger, component_index);
        component_index += n_components(*vec[comp]);
      }
  }



  //
  // update_ghost_values_finish
  //

  // for block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            typename std::enable_if<IsBlockVector<VectorStruct>::value,
                                    VectorStruct>::type * = nullptr>
  void
  update_ghost_values_finish(
    const VectorStruct &                                  vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    if (get_communication_block_size(vec) < vec.n_blocks())
      {
        // do nothing, everything has already been completed in the _start()
        // call
      }
    else
      for (unsigned int i = 0; i < vec.n_blocks(); ++i)
        update_ghost_values_finish(vec.block(i), exchanger, channel + i);
  }



  // for non-block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            typename std::enable_if<!IsBlockVector<VectorStruct>::value,
                                    VectorStruct>::type * = nullptr>
  void
  update_ghost_values_finish(
    const VectorStruct &                                  vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    exchanger.update_ghost_values_finish(channel, vec);
  }



  // for vector of vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  update_ghost_values_finish(
    const std::vector<VectorStruct> &                     vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); comp++)
      {
        update_ghost_values_finish(vec[comp], exchanger, component_index);
        component_index += n_components(vec[comp]);
      }
  }



  // for vector of pointers to vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  update_ghost_values_finish(
    const std::vector<VectorStruct *> &                   vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); comp++)
      {
        update_ghost_values_finish(*vec[comp], exchanger, component_index);
        component_index += n_components(*vec[comp]);
      }
  }



  //
  // compress_start
  //

  // for block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            typename std::enable_if<IsBlockVector<VectorStruct>::value,
                                    VectorStruct>::type * = nullptr>
  inline void
  compress_start(
    VectorStruct &                                        vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    if (get_communication_block_size(vec) < vec.n_blocks())
      vec.compress(dealii::VectorOperation::add);
    else
      for (unsigned int i = 0; i < vec.n_blocks(); ++i)
        compress_start(vec.block(i), exchanger, channel + i);
  }



  // for non-block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            typename std::enable_if<!IsBlockVector<VectorStruct>::value,
                                    VectorStruct>::type * = nullptr>
  inline void
  compress_start(
    VectorStruct &                                        vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    exchanger.compress_start(channel, vec);
  }



  // for std::vector of vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  compress_start(
    std::vector<VectorStruct> &                           vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); comp++)
      {
        compress_start(vec[comp], exchanger, component_index);
        component_index += n_components(vec[comp]);
      }
  }



  // for std::vector of pointer to vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  compress_start(
    std::vector<VectorStruct *> &                         vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); comp++)
      {
        compress_start(*vec[comp], exchanger, component_index);
        component_index += n_components(*vec[comp]);
      }
  }



  //
  // compress_finish
  //

  // for block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            typename std::enable_if<IsBlockVector<VectorStruct>::value,
                                    VectorStruct>::type * = nullptr>
  inline void
  compress_finish(
    VectorStruct &                                        vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    if (get_communication_block_size(vec) < vec.n_blocks())
      {
        // do nothing, everything has already been completed in the _start()
        // call
      }
    else
      for (unsigned int i = 0; i < vec.n_blocks(); ++i)
        compress_finish(vec.block(i), exchanger, channel + i);
  }



  // for non-block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            typename std::enable_if<!IsBlockVector<VectorStruct>::value,
                                    VectorStruct>::type * = nullptr>
  inline void
  compress_finish(
    VectorStruct &                                        vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    exchanger.compress_finish(channel, vec);
  }



  // for std::vector of vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  compress_finish(
    std::vector<VectorStruct> &                           vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); comp++)
      {
        compress_finish(vec[comp], exchanger, component_index);
        component_index += n_components(vec[comp]);
      }
  }



  // for std::vector of pointer to vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  compress_finish(
    std::vector<VectorStruct *> &                         vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); comp++)
      {
        compress_finish(*vec[comp], exchanger, component_index);
        component_index += n_components(*vec[comp]);
      }
  }



  //
  // reset_ghost_values:
  //
  // if the input vector did not have ghosts imported, clear them here again
  // in order to avoid subsequent operations e.g. in linear solvers to work
  // with ghosts all the time
  //

  // for block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            typename std::enable_if<IsBlockVector<VectorStruct>::value,
                                    VectorStruct>::type * = nullptr>
  inline void
  reset_ghost_values(
    const VectorStruct &                                  vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    // return immediately if there is nothing to do.
    if (exchanger.ghosts_were_set == true)
      return;

    for (unsigned int i = 0; i < vec.n_blocks(); ++i)
      reset_ghost_values(vec.block(i), exchanger);
  }



  // for non-block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            typename std::enable_if<!IsBlockVector<VectorStruct>::value,
                                    VectorStruct>::type * = nullptr>
  inline void
  reset_ghost_values(
    const VectorStruct &                                  vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    exchanger.reset_ghost_values(vec);
  }



  // for std::vector of vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  reset_ghost_values(
    const std::vector<VectorStruct> &                     vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    // return immediately if there is nothing to do.
    if (exchanger.ghosts_were_set == true)
      return;

    for (unsigned int comp = 0; comp < vec.size(); comp++)
      reset_ghost_values(vec[comp], exchanger);
  }



  // for std::vector of pointer to vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  reset_ghost_values(
    const std::vector<VectorStruct *> &                   vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    // return immediately if there is nothing to do.
    if (exchanger.ghosts_were_set == true)
      return;

    for (unsigned int comp = 0; comp < vec.size(); comp++)
      reset_ghost_values(*vec[comp], exchanger);
  }



  //
  // zero_vector_region
  //

  // for block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            typename std::enable_if<IsBlockVector<VectorStruct>::value,
                                    VectorStruct>::type * = nullptr>
  inline void
  zero_vector_region(
    const unsigned int                                    range_index,
    VectorStruct &                                        vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    for (unsigned int i = 0; i < vec.n_blocks(); ++i)
      exchanger.zero_vector_region(range_index, vec.block(i));
  }



  // for non-block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            typename std::enable_if<!IsBlockVector<VectorStruct>::value,
                                    VectorStruct>::type * = nullptr>
  inline void
  zero_vector_region(
    const unsigned int                                    range_index,
    VectorStruct &                                        vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    exchanger.zero_vector_region(range_index, vec);
  }



  // for std::vector of vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  zero_vector_region(
    const unsigned int                                    range_index,
    std::vector<VectorStruct> &                           vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    for (unsigned int comp = 0; comp < vec.size(); comp++)
      zero_vector_region(range_index, vec[comp], exchanger);
  }



  // for std::vector of pointers to vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  zero_vector_region(
    const unsigned int                                    range_index,
    std::vector<VectorStruct *> &                         vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    for (unsigned int comp = 0; comp < vec.size(); comp++)
      zero_vector_region(range_index, *vec[comp], exchanger);
  }



  namespace MatrixFreeFunctions
  {
    // struct to select between a const interface and a non-const interface
    // for MFWorker
    template <typename, typename, typename, typename, bool>
    struct InterfaceSelector
    {};

    // Version for constant functions
    template <typename MF,
              typename InVector,
              typename OutVector,
              typename Container>
    struct InterfaceSelector<MF, InVector, OutVector, Container, true>
    {
      using function_type = void (Container::*)(
        const MF &,
        OutVector &,
        const InVector &,
        const std::pair<unsigned int, unsigned int> &) const;
    };

    // Version for non-constant functions
    template <typename MF,
              typename InVector,
              typename OutVector,
              typename Container>
    struct InterfaceSelector<MF, InVector, OutVector, Container, false>
    {
      using function_type =
        void (Container::*)(const MF &,
                            OutVector &,
                            const InVector &,
                            const std::pair<unsigned int, unsigned int> &);
    };
  } // namespace MatrixFreeFunctions



  // A implementation class for the worker object that runs the various
  // operations we want to perform during the matrix-free loop
  template <typename MF,
            typename InVector,
            typename OutVector,
            typename Container,
            bool is_constant>
  class MFWorker : public MFWorkerInterface
  {
  public:
    // An alias to make the arguments further down more readable
    using function_type = typename MatrixFreeFunctions::
      InterfaceSelector<MF, InVector, OutVector, Container, is_constant>::
        function_type;

    // constructor, binds all the arguments to this class
    MFWorker(const MF &                           matrix_free,
             const InVector &                     src,
             OutVector &                          dst,
             const bool                           zero_dst_vector_setting,
             const Container &                    container,
             function_type                        cell_function,
             function_type                        face_function,
             function_type                        boundary_function,
             const typename MF::DataAccessOnFaces src_vector_face_access =
               MF::DataAccessOnFaces::none,
             const typename MF::DataAccessOnFaces dst_vector_face_access =
               MF::DataAccessOnFaces::none,
             const std::function<void(const unsigned int, const unsigned int)>
               &operation_before_loop = {},
             const std::function<void(const unsigned int, const unsigned int)>
               &                operation_after_loop       = {},
             const unsigned int dof_handler_index_pre_post = 0)
      : matrix_free(matrix_free)
      , container(const_cast<Container &>(container))
      , cell_function(cell_function)
      , face_function(face_function)
      , boundary_function(boundary_function)
      , src(src)
      , dst(dst)
      , src_data_exchanger(matrix_free,
                           src_vector_face_access,
                           n_components(src))
      , dst_data_exchanger(matrix_free,
                           dst_vector_face_access,
                           n_components(dst))
      , src_and_dst_are_same(PointerComparison::equal(&src, &dst))
      , zero_dst_vector_setting(zero_dst_vector_setting &&
                                !src_and_dst_are_same)
      , operation_before_loop(operation_before_loop)
      , operation_after_loop(operation_after_loop)
      , dof_handler_index_pre_post(dof_handler_index_pre_post)
    {}

    // Runs the cell work. If no function is given, nothing is done
    virtual void
    cell(const std::pair<unsigned int, unsigned int> &cell_range) override
    {
      if (cell_function != nullptr && cell_range.second > cell_range.first)
        for (unsigned int i = 0; i < matrix_free.n_active_fe_indices(); ++i)
          {
            const auto cell_subrange =
              matrix_free.create_cell_subrange_hp_by_index(cell_range, i);

            if (cell_subrange.second <= cell_subrange.first)
              continue;

            (container.*
             cell_function)(matrix_free, this->dst, this->src, cell_subrange);
          }
    }

    virtual void
    cell(const unsigned int range_index) override
    {
      process_range(cell_function,
                    matrix_free.get_task_info().cell_partition_data_hp_ptr,
                    matrix_free.get_task_info().cell_partition_data_hp,
                    range_index);
    }

    virtual void
    face(const unsigned int range_index) override
    {
      process_range(face_function,
                    matrix_free.get_task_info().face_partition_data_hp_ptr,
                    matrix_free.get_task_info().face_partition_data_hp,
                    range_index);
    }

    virtual void
    boundary(const unsigned int range_index) override
    {
      process_range(boundary_function,
                    matrix_free.get_task_info().boundary_partition_data_hp_ptr,
                    matrix_free.get_task_info().boundary_partition_data_hp,
                    range_index);
    }

  private:
    void
    process_range(const function_type &            fu,
                  const std::vector<unsigned int> &ptr,
                  const std::vector<unsigned int> &data,
                  const unsigned int               range_index)
    {
      if (fu == nullptr)
        return;

      for (unsigned int i = ptr[range_index]; i < ptr[range_index + 1]; ++i)
        (container.*fu)(matrix_free,
                        this->dst,
                        this->src,
                        std::make_pair(data[2 * i], data[2 * i + 1]));
    }

  public:
    // Starts the communication for the update ghost values operation. We
    // cannot call this update if ghost and destination are the same because
    // that would introduce spurious entries in the destination (there is also
    // the problem that reading from a vector that we also write to is usually
    // not intended in case there is overlap, but this is up to the
    // application code to decide and we cannot catch this case here).
    virtual void
    vector_update_ghosts_start() override
    {
      if (!src_and_dst_are_same)
        internal::update_ghost_values_start(src, src_data_exchanger);
    }

    // Finishes the communication for the update ghost values operation
    virtual void
    vector_update_ghosts_finish() override
    {
      if (!src_and_dst_are_same)
        internal::update_ghost_values_finish(src, src_data_exchanger);
    }

    // Starts the communication for the vector compress operation
    virtual void
    vector_compress_start() override
    {
      internal::compress_start(dst, dst_data_exchanger);
    }

    // Finishes the communication for the vector compress operation
    virtual void
    vector_compress_finish() override
    {
      internal::compress_finish(dst, dst_data_exchanger);
      if (!src_and_dst_are_same)
        internal::reset_ghost_values(src, src_data_exchanger);
    }

    // Zeros the given input vector
    virtual void
    zero_dst_vector_range(const unsigned int range_index) override
    {
      if (zero_dst_vector_setting)
        internal::zero_vector_region(range_index, dst, dst_data_exchanger);
    }

    virtual void
    cell_loop_pre_range(const unsigned int range_index) override
    {
      if (operation_before_loop)
        {
          const internal::MatrixFreeFunctions::DoFInfo &dof_info =
            matrix_free.get_dof_info(dof_handler_index_pre_post);
          if (range_index == numbers::invalid_unsigned_int)
            {
              // Case with threaded loop -> currently no overlap implemented
              dealii::parallel::apply_to_subranges(
                0U,
                dof_info.vector_partitioner->locally_owned_size(),
                operation_before_loop,
                internal::VectorImplementation::minimum_parallel_grain_size);
            }
          else
            {
              AssertIndexRange(range_index,
                               dof_info.cell_loop_pre_list_index.size() - 1);
              for (unsigned int id =
                     dof_info.cell_loop_pre_list_index[range_index];
                   id != dof_info.cell_loop_pre_list_index[range_index + 1];
                   ++id)
                operation_before_loop(dof_info.cell_loop_pre_list[id].first,
                                      dof_info.cell_loop_pre_list[id].second);
            }
        }
    }

    virtual void
    cell_loop_post_range(const unsigned int range_index) override
    {
      if (operation_after_loop)
        {
          const internal::MatrixFreeFunctions::DoFInfo &dof_info =
            matrix_free.get_dof_info(dof_handler_index_pre_post);
          if (range_index == numbers::invalid_unsigned_int)
            {
              // Case with threaded loop -> currently no overlap implemented
              dealii::parallel::apply_to_subranges(
                0U,
                dof_info.vector_partitioner->locally_owned_size(),
                operation_after_loop,
                internal::VectorImplementation::minimum_parallel_grain_size);
            }
          else
            {
              AssertIndexRange(range_index,
                               dof_info.cell_loop_post_list_index.size() - 1);
              for (unsigned int id =
                     dof_info.cell_loop_post_list_index[range_index];
                   id != dof_info.cell_loop_post_list_index[range_index + 1];
                   ++id)
                operation_after_loop(dof_info.cell_loop_post_list[id].first,
                                     dof_info.cell_loop_post_list[id].second);
            }
        }
    }

  private:
    const MF &    matrix_free;
    Container &   container;
    function_type cell_function;
    function_type face_function;
    function_type boundary_function;

    const InVector &src;
    OutVector &     dst;
    VectorDataExchange<MF::dimension,
                       typename MF::value_type,
                       typename MF::vectorized_value_type>
      src_data_exchanger;
    VectorDataExchange<MF::dimension,
                       typename MF::value_type,
                       typename MF::vectorized_value_type>
               dst_data_exchanger;
    const bool src_and_dst_are_same;
    const bool zero_dst_vector_setting;
    const std::function<void(const unsigned int, const unsigned int)>
      operation_before_loop;
    const std::function<void(const unsigned int, const unsigned int)>
                       operation_after_loop;
    const unsigned int dof_handler_index_pre_post;
  };



  /**
   * 一个内部类，将三个函数指针转换为上述带有虚拟函数的方案。
   *
   */
  template <class MF, typename InVector, typename OutVector>
  struct MFClassWrapper
  {
    using function_type =
      std::function<void(const MF &,
                         OutVector &,
                         const InVector &,
                         const std::pair<unsigned int, unsigned int> &)>;

    MFClassWrapper(const function_type cell,
                   const function_type face,
                   const function_type boundary)
      : cell(cell)
      , face(face)
      , boundary(boundary)
    {}

    void
    cell_integrator(const MF &                                   mf,
                    OutVector &                                  dst,
                    const InVector &                             src,
                    const std::pair<unsigned int, unsigned int> &range) const
    {
      if (cell)
        cell(mf, dst, src, range);
    }

    void
    face_integrator(const MF &                                   mf,
                    OutVector &                                  dst,
                    const InVector &                             src,
                    const std::pair<unsigned int, unsigned int> &range) const
    {
      if (face)
        face(mf, dst, src, range);
    }

    void
    boundary_integrator(
      const MF &                                   mf,
      OutVector &                                  dst,
      const InVector &                             src,
      const std::pair<unsigned int, unsigned int> &range) const
    {
      if (boundary)
        boundary(mf, dst, src, range);
    }

    const function_type cell;
    const function_type face;
    const function_type boundary;
  };

} // end of namespace internal



template <int dim, typename Number, typename VectorizedArrayType>
template <typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::cell_loop(
  const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                           OutVector &,
                           const InVector &,
                           const std::pair<unsigned int, unsigned int> &)>
    &             cell_operation,
  OutVector &     dst,
  const InVector &src,
  const bool      zero_dst_vector) const
{
  using Wrapper =
    internal::MFClassWrapper<MatrixFree<dim, Number, VectorizedArrayType>,
                             InVector,
                             OutVector>;
  Wrapper wrap(cell_operation, nullptr, nullptr);
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     Wrapper,
                     true>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           wrap,
           &Wrapper::cell_integrator,
           &Wrapper::face_integrator,
           &Wrapper::boundary_integrator);

  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::cell_loop(
  const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                           OutVector &,
                           const InVector &,
                           const std::pair<unsigned int, unsigned int> &)>
    &             cell_operation,
  OutVector &     dst,
  const InVector &src,
  const std::function<void(const unsigned int, const unsigned int)>
    &operation_before_loop,
  const std::function<void(const unsigned int, const unsigned int)>
    &                operation_after_loop,
  const unsigned int dof_handler_index_pre_post) const
{
  using Wrapper =
    internal::MFClassWrapper<MatrixFree<dim, Number, VectorizedArrayType>,
                             InVector,
                             OutVector>;
  Wrapper wrap(cell_operation, nullptr, nullptr);
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     Wrapper,
                     true>
    worker(*this,
           src,
           dst,
           false,
           wrap,
           &Wrapper::cell_integrator,
           &Wrapper::face_integrator,
           &Wrapper::boundary_integrator,
           DataAccessOnFaces::none,
           DataAccessOnFaces::none,
           operation_before_loop,
           operation_after_loop,
           dof_handler_index_pre_post);

  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::loop(
  const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                           OutVector &,
                           const InVector &,
                           const std::pair<unsigned int, unsigned int> &)>
    &cell_operation,
  const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                           OutVector &,
                           const InVector &,
                           const std::pair<unsigned int, unsigned int> &)>
    &face_operation,
  const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                           OutVector &,
                           const InVector &,
                           const std::pair<unsigned int, unsigned int> &)>
    &                     boundary_operation,
  OutVector &             dst,
  const InVector &        src,
  const bool              zero_dst_vector,
  const DataAccessOnFaces dst_vector_face_access,
  const DataAccessOnFaces src_vector_face_access) const
{
  using Wrapper =
    internal::MFClassWrapper<MatrixFree<dim, Number, VectorizedArrayType>,
                             InVector,
                             OutVector>;
  Wrapper wrap(cell_operation, face_operation, boundary_operation);
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     Wrapper,
                     true>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           wrap,
           &Wrapper::cell_integrator,
           &Wrapper::face_integrator,
           &Wrapper::boundary_integrator,
           src_vector_face_access,
           dst_vector_face_access);

  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::cell_loop(
  void (CLASS::*function_pointer)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &) const,
  const CLASS *   owning_class,
  OutVector &     dst,
  const InVector &src,
  const bool      zero_dst_vector) const
{
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     true>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           *owning_class,
           function_pointer,
           nullptr,
           nullptr);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::cell_loop(
  void (CLASS::*function_pointer)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &) const,
  const CLASS *   owning_class,
  OutVector &     dst,
  const InVector &src,
  const std::function<void(const unsigned int, const unsigned int)>
    &operation_before_loop,
  const std::function<void(const unsigned int, const unsigned int)>
    &                operation_after_loop,
  const unsigned int dof_handler_index_pre_post) const
{
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     true>
    worker(*this,
           src,
           dst,
           false,
           *owning_class,
           function_pointer,
           nullptr,
           nullptr,
           DataAccessOnFaces::none,
           DataAccessOnFaces::none,
           operation_before_loop,
           operation_after_loop,
           dof_handler_index_pre_post);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::loop(
  void (CLASS::*cell_operation)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &) const,
  void (CLASS::*face_operation)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &) const,
  void (CLASS::*boundary_operation)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &) const,
  const CLASS *           owning_class,
  OutVector &             dst,
  const InVector &        src,
  const bool              zero_dst_vector,
  const DataAccessOnFaces dst_vector_face_access,
  const DataAccessOnFaces src_vector_face_access) const
{
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     true>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           *owning_class,
           cell_operation,
           face_operation,
           boundary_operation,
           src_vector_face_access,
           dst_vector_face_access);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::cell_loop(
  void (CLASS::*function_pointer)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &),
  CLASS *         owning_class,
  OutVector &     dst,
  const InVector &src,
  const bool      zero_dst_vector) const
{
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     false>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           *owning_class,
           function_pointer,
           nullptr,
           nullptr);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::cell_loop(
  void (CLASS::*function_pointer)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &),
  CLASS *         owning_class,
  OutVector &     dst,
  const InVector &src,
  const std::function<void(const unsigned int, const unsigned int)>
    &operation_before_loop,
  const std::function<void(const unsigned int, const unsigned int)>
    &                operation_after_loop,
  const unsigned int dof_handler_index_pre_post) const
{
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     false>
    worker(*this,
           src,
           dst,
           false,
           *owning_class,
           function_pointer,
           nullptr,
           nullptr,
           DataAccessOnFaces::none,
           DataAccessOnFaces::none,
           operation_before_loop,
           operation_after_loop,
           dof_handler_index_pre_post);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::loop(
  void (CLASS::*cell_operation)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &),
  void (CLASS::*face_operation)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &),
  void (CLASS::*boundary_operation)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &),
  CLASS *                 owning_class,
  OutVector &             dst,
  const InVector &        src,
  const bool              zero_dst_vector,
  const DataAccessOnFaces dst_vector_face_access,
  const DataAccessOnFaces src_vector_face_access) const
{
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     false>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           *owning_class,
           cell_operation,
           face_operation,
           boundary_operation,
           src_vector_face_access,
           dst_vector_face_access);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::loop_cell_centric(
  void (CLASS::*function_pointer)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &) const,
  const CLASS *           owning_class,
  OutVector &             dst,
  const InVector &        src,
  const bool              zero_dst_vector,
  const DataAccessOnFaces src_vector_face_access) const
{
  auto src_vector_face_access_temp = src_vector_face_access;
  if (DataAccessOnFaces::gradients == src_vector_face_access_temp)
    src_vector_face_access_temp = DataAccessOnFaces::gradients_all_faces;
  else if (DataAccessOnFaces::values == src_vector_face_access_temp)
    src_vector_face_access_temp = DataAccessOnFaces::values_all_faces;

  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     true>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           *owning_class,
           function_pointer,
           nullptr,
           nullptr,
           src_vector_face_access_temp,
           DataAccessOnFaces::none);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::loop_cell_centric(
  void (CLASS::*function_pointer)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &),
  CLASS *                 owning_class,
  OutVector &             dst,
  const InVector &        src,
  const bool              zero_dst_vector,
  const DataAccessOnFaces src_vector_face_access) const
{
  auto src_vector_face_access_temp = src_vector_face_access;
  if (DataAccessOnFaces::gradients == src_vector_face_access_temp)
    src_vector_face_access_temp = DataAccessOnFaces::gradients_all_faces;
  else if (DataAccessOnFaces::values == src_vector_face_access_temp)
    src_vector_face_access_temp = DataAccessOnFaces::values_all_faces;

  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     false>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           *owning_class,
           function_pointer,
           nullptr,
           nullptr,
           src_vector_face_access_temp,
           DataAccessOnFaces::none);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::loop_cell_centric(
  const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                           OutVector &,
                           const InVector &,
                           const std::pair<unsigned int, unsigned int> &)>
    &                     cell_operation,
  OutVector &             dst,
  const InVector &        src,
  const bool              zero_dst_vector,
  const DataAccessOnFaces src_vector_face_access) const
{
  auto src_vector_face_access_temp = src_vector_face_access;
  if (DataAccessOnFaces::gradients == src_vector_face_access_temp)
    src_vector_face_access_temp = DataAccessOnFaces::gradients_all_faces;
  else if (DataAccessOnFaces::values == src_vector_face_access_temp)
    src_vector_face_access_temp = DataAccessOnFaces::values_all_faces;

  using Wrapper =
    internal::MFClassWrapper<MatrixFree<dim, Number, VectorizedArrayType>,
                             InVector,
                             OutVector>;
  Wrapper wrap(cell_operation, nullptr, nullptr);

  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     Wrapper,
                     true>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           wrap,
           &Wrapper::cell_integrator,
           &Wrapper::face_integrator,
           &Wrapper::boundary_integrator,
           src_vector_face_access_temp,
           DataAccessOnFaces::none);
  task_info.loop(worker);
}


#endif // ifndef DOXYGEN



DEAL_II_NAMESPACE_CLOSE

#endif


