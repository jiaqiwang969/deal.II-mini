//include/deal.II-translator/distributed/tria_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2021 by the deal.II authors
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

#ifndef dealii_distributed_tria_h
#define dealii_distributed_tria_h


#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/distributed/p4est_wrappers.h>
#include <deal.II/distributed/tria_base.h>

#include <deal.II/grid/tria.h>

#include <boost/range/iterator_range.hpp>

#include <functional>
#include <list>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef DEAL_II_WITH_MPI
#  include <mpi.h>
#endif

#ifdef DEAL_II_WITH_P4EST
#  include <p4est.h>
#  include <p4est_connectivity.h>
#  include <p4est_ghost.h>
#  include <p8est.h>
#  include <p8est_connectivity.h>
#  include <p8est_ghost.h>
#endif


DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_P4EST

// Forward declarations
#  ifndef DOXYGEN

namespace FETools
{
  namespace internal
  {
    template <int, int, class>
    class ExtrapolateImplementation;
  }
} // namespace FETools

// forward declaration of the data type for periodic face pairs
namespace GridTools
{
  template <typename CellIterator>
  struct PeriodicFacePair;
}

namespace parallel
{
  namespace distributed
  {
    template <int, int>
    class TemporarilyMatchRefineFlags;
  }
} // namespace parallel
#  endif

namespace parallel
{
  namespace distributed
  {
    /**
     * 这个类的作用与 dealii::Triangulation 类相似，但它在使用MPI时将网格分布在多个不同的处理器上。该类的接口与 dealii::Triangulation 类相比并没有增加很多，但在其内部有一些困难的算法，确保我们总是有一个负载平衡的、完全分布的网格。这个类的使用在  step-40  ,  step-32  ,  @ref distributed  文档模块，以及  @ref distributed_paper  中都有解释。    更多信息见那里。这个类满足了  @ref ConceptMeshType  "MeshType概念"
     * 。
     * @note
     * 这个类不支持各向异性的细化，因为它所依赖的p4est库不支持这个。试图对单元进行各向异性的细化将导致错误。
     * @note
     * 目前不支持分布1d三角形。            <h3> Interaction with
     * boundary description </h3>
     * 细化和粗化分布式三角形是一个复杂的过程，因为单元可能必须从一个处理器迁移到另一个。在一个单一的处理器上，将我们想要存储在这里的全局网格的一部分具体化，可能需要对本地存储的单元集进行多次的细化和粗化，直到我们最终从上一个三角形得到下一个。这个过程在
     * @ref distributed_paper  中有更详细的描述。
     * 不幸的是，在这个过程中，一些信息可能会丢失，这些信息是由用户代码设置的，并从母单元继承到子单元，但如果单元从一个处理器迁移到另一个处理器，这些信息就不会随单元一起移动。
     * 一个例子是边界指示器。例如，假设你从一个单一的单元开始，在全局范围内被精炼一次，产生四个子单元。如果你有四个处理器，每个处理器拥有一个单元。假设处理器1将其拥有的单元格的外部边界指标设置为42。由于处理器0并不拥有这个单元，所以它并没有设置这个单元的幽灵单元副本的边界指标。现在，假设我们做了几个网格细化循环，最后这个处理器突然发现自己是这个单元的所有者。如果边界指标42意味着我们需要沿着这个边界整合诺伊曼边界条件，那么处理器0会忘记这样做，因为它从未将这个单元的边界指标设置为42。
     * 避免这种困境的方法是确保每次精炼平行三角剖面时，设置边界指示器或材料ID等事情立即完成。这对顺序三角剖分来说是没有必要的，因为在那里，这些标志是由母单元继承到子单元的，即使它被细化了，子单元后来又被粗化了，这些标志仍然留在单元中，但这对分布式三角剖分来说是不成立的。更加困难的是，在细化一个平行分布式三角形的过程中，三角形可能会多次调用
     * dealii::Triangulation::execute_coarsening_and_refinement
     * ，这个函数需要知道边界。换句话说，<i>not</i>只是在新创建的面设置边界指标就足够了，只有<i>after</i>调用
     * <tt>distributed::parallel::TriangulationBase::execute_coarsening_and_refinement</tt>:
     * 时，实际上必须在该函数仍在运行时发生。
     * 做到这一点的方法是编写一个设置边界指标的函数，该函数将被
     * dealii::Triangulation 类所调用。
     * 三角形并没有为被调用的函数提供一个指向自身的指针，也没有任何其他信息，所以诀窍是将这些信息引入函数。C++为此提供了一个很好的机制，最好用一个例子来解释。
     * @code
     * #include <functional>
     *
     * template <int dim>
     * void set_boundary_ids (
     * parallel::distributed::Triangulation<dim> &triangulation)
     * {
     * ... set boundary indicators on the triangulation object ...
     * }
     *
     * template <int dim>
     * void
     * MyClass<dim>::create_coarse_mesh (
     * parallel::distributed::Triangulation<dim> &coarse_grid) const
     * {
     * ... create the coarse mesh ...
     *
     * coarse_grid.signals.post_refinement.connect(
     *   [&coarse_grid](){
     *     set_boundary_ids<dim>(coarse_grid);
     *   });
     * }
     * @endcode
     * 作为参数传递给 <code>connect</code>
     * 的对象是一个可以像函数一样被调用的对象，没有参数。它是通过包装一个函数来实现的，该函数实际上需要一个参数，但这个参数在创建lambda函数时被存储为对粗略网格三角的引用。在每个细化步骤之后，三角化将调用如此创建的对象，而该对象又将调用
     * <code>set_boundary_ids<dim></code>
     * ，并将粗略网格的引用作为参数。
     * 这种方法可以被推广。在上面的例子中，我们使用了一个将被调用的全局函数。然而，有时需要这个函数实际上是生成网格的类的一个成员函数，例如，因为它需要访问运行时参数。这可以通过以下方式实现：假设
     * <code>set_boundary_ids()</code> 函数已经被声明为
     * <code>MyClass</code>
     * 类的（非静态的，但可能是私有的）成员函数，那么下面的方法就可以实现。
     * @code
     * #include <functional>
     *
     * template <int dim>
     * void
     * MyClass<dim>::set_boundary_ids (
     * parallel::distributed::Triangulation<dim> &triangulation) const
     * {
     * ... set boundary indicators on the triangulation object ...
     * }
     *
     * template <int dim>
     * void
     * MyClass<dim>::create_coarse_mesh (
     * parallel::distributed::Triangulation<dim> &coarse_grid) const
     * {
     * ... create the coarse mesh ...
     *
     * coarse_grid.signals.post_refinement.connect(
     *   [this, &coarse_grid]()
     *   {
     *     this->set_boundary_ids(coarse_grid);
     *   });
     * }
     * @endcode
     * 上面的lambda函数又是一个可以像全局函数一样被调用的对象，没有参数，这个对象反过来调用当前对象的成员函数
     * <code>set_boundary_ids</code>
     * ，有一个对三角形的引用来工作。注意，由于
     * <code>create_coarse_mesh</code> 函数被声明为 <code>const</code>
     * ，所以 <code>set_boundary_ids</code> 函数也必须被声明为
     * <code>const</code> 。        <b>Note:</b>由于与
     * parallel::distributed::Triangulation
     * 的实现方式有关的原因，每次实际细化三角形时，附加到三角形细化后信号的函数会被调用一次以上，有时会被调用几次。
     * @ingroup distributed
     *
     */
    template <int dim, int spacedim = dim>
    class Triangulation
      : public dealii::parallel::DistributedTriangulationBase<dim, spacedim>
    {
    public:
      /**
       * 一个别名，用于识别单元格迭代器。迭代器的概念在 @ref Iterators "迭代器文档模块 "
       * 中有详细的讨论。
       * 当前的别名用于标识三角形中的单元格。你可以在基类自己的别名中找到它所指的确切类型，但它应该是TriaIterator<CellAccessor<dim,spacedim>>。TriaIterator类的工作原理就像一个指针，当你解除引用时，会产生一个类型为CellAccessor的对象。CellAccessor是一个标识三角形中单元格特定属性的类，但它派生（并因此继承）于TriaAccessor，后者描述了你可以对三角形中更多的一般对象（线、面以及单元格）提出的要求。
       * @ingroup Iterators
       *
       */
      using cell_iterator =
        typename dealii::Triangulation<dim, spacedim>::cell_iterator;

      /**
       * @ingroup Iterators
       *
       */
      using active_cell_iterator =
        typename dealii::Triangulation<dim, spacedim>::active_cell_iterator;

      using CellStatus =
        typename dealii::Triangulation<dim, spacedim>::CellStatus;

      /**
       * 分布式三角形的配置标志将在构造函数中设置。设置可以使用位法OR进行组合。
       *
       */
      enum Settings
      {
        /**
         * 默认设置，其他选项被禁用。
         *
         */
        default_setting = 0x0,
        /**
         * 如果设置了，每次在p4est中发生重新分区时，deal.II网格将从粗略的网格中重构出来。这可能有点昂贵，但是保证了相同的内存布局，从而保证了deal.II网格的单元排序。由于装配是在deal.II单元排序中完成的，这个标志对于在快照/恢复后获得可重复的行为是必需的。
         *
         */
        mesh_reconstruction_after_repartitioning = 0x1,
        /**
         * 这个标志需要被设置以使用几何多网格功能。这个选项需要额外的计算和通信。
         *
         */
        construct_multigrid_hierarchy = 0x2,
        /**
         * 设置这个标志将禁止在细化周期后自动重新划分单元。它可以通过调用repartition()手动执行。
         *
         */
        no_automatic_repartitioning = 0x4
      };



      /**
       * 构造函数。              @param  mpi_communicator
       * 用于三角测量的MPI通信器。              @param  smooth_grid
       * 对网格进行平滑处理的程度和种类。参见
       * dealii::Triangulation
       * 类中关于可应用的平滑操作种类的描述。
       * @param  设置 参见设置枚举器的描述。
       * 为smooth_grid提供 <code>construct_multigrid_hierarchy</code>
       * 执行 <code>Triangulation::limit_level_difference_at_vertices</code>
       * 。
       * @note  这个类目前不支持基类所提供的
       * <code>check_for_distorted_cells</code> 参数。
       * @note
       * 虽然可以将基类中列出的所有网格平滑标志传递给该类型的对象，但如果这些平滑选项需要了解不属于该处理器本地的单元上的细化/粗化标志，则并不总是可以兑现所有这些平滑选项。因此，对于其中的一些标志，平行三角形的最终单元数可能取决于它被分割成的处理器的数量。另一方面，如果不传递平滑标志，如果你总是标记相同的网格单元，你将总是得到完全相同的精炼网格，与三角形划分的处理器的数量无关。
       *
       */
      explicit Triangulation(
        const MPI_Comm &mpi_communicator,
        const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
                       smooth_grid = (dealii::Triangulation<dim, spacedim>::none),
        const Settings settings    = default_setting);

      /**
       * 解构器。
       *
       */
      virtual ~Triangulation() override;

      /**
       * 通过删除所有的数据，将这个三角网格重设为原始状态。
       * 请注意，只有当这个对象不再存在任何订阅时，才允许这个操作，例如使用它的DoFHandler对象。
       *
       */
      virtual void
      clear() override;

      /**
       * 如果支持多级层次结构并已构建，则返回。
       *
       */
      bool
      is_multilevel_hierarchy_constructed() const override;

      /**
       * 在森林间传输数据。            除了实际的 @p
       * parallel_forest,
       * 已经被细化和重新分区外，这个函数还需要关于其先前状态的信息，即每个处理器的p4est的sc_array中本地拥有的间隔。这些信息需要从旧的p4est对象中记忆复制出来，并且必须通过参数
       * @p previous_global_first_quadrant. 来提供，数据必须事先用
       * DistributedTriangulationBase::DataTransfer::pack_data(). 打包。
       *
       */
      void
      execute_transfer(
        const typename dealii::internal::p4est::types<dim>::forest
          *parallel_forest,
        const typename dealii::internal::p4est::types<dim>::gloidx
          *previous_global_first_quadrant);

      /**
       * 实现与基类中相同的功能。
       * @note 该函数可用于将一个序列三角图复制到
       * parallel::distributed::Triangulation
       * ，但仅当该序列三角图从未被细化。
       *
       */
      virtual void
      copy_triangulation(
        const dealii::Triangulation<dim, spacedim> &other_tria) override;

      /**
       * 按照基类中的记录创建一个三角剖面。
       * 这个函数还设置了各种必要的数据结构，以便将一个网格分布在多个处理器上。一旦网格被细化，这将是必要的，尽管我们将始终在所有处理器上保留这个函数生成的整个粗略网格。
       *
       */
      virtual void
      create_triangulation(const std::vector<Point<spacedim>> &vertices,
                           const std::vector<CellData<dim>> &  cells,
                           const SubCellData &subcelldata) override;

      /**
       * @copydoc   Triangulation::create_triangulation() .
       * @note  还没有实现。
       *
       */
      virtual void
      create_triangulation(
        const TriangulationDescription::Description<dim, spacedim>
          &construction_data) override;

      /**
       * 根据设置的细化和粗化标志，对网格进行粗化和细化。
       * 由于当前处理器只控制它所拥有的单元（即那些<code>cell-
       * @>subdomain_id()  == this-  @>locally_owned_subdomain()</code>),
       * 细化和粗化标志只对那些本地拥有的单元进行控制。其他单元也可以设置标志（事实上，如果你调用
       * dealii::Triangulation::prepare_coarsening_and_refinement())
       * ，也可以经常设置标志，但基本上会被忽略：细化全局网格的决定只受本地拥有的单元上设置的标志影响。
       * @note
       * 这个函数默认以这样的方式划分网格，即所有处理器上的单元数大致相等。如果你想为分区设置权重，例如，因为某些单元的计算成本比其他单元高，你可以使用信号cell_weight，如
       * dealii::Triangulation
       * 类中所记录的。这个函数将检查一个函数是否与信号相连，如果是则使用它。
       * 如果你喜欢自己只在用户定义的时间间隔内重新划分网格，你可以在构造函数中传递
       * parallel::distributed::Triangulation::no_automatic_repartitioning
       * 标志来创建你的三角形对象，这样可以确保调用当前函数只细化和粗化三角形，而不划分。然后你可以手动调用repartition()函数。
       * cell_weights信号的用法在这两种情况下是相同的，如果一个函数连接到该信号，它将被用来平衡计算的权重，否则就平衡单元格的数量。
       *
       */
      virtual void
      execute_coarsening_and_refinement() override;

      /**
       * 覆盖基类中prepare_coarsening_and_refinement的实现。如果启用了周期性边界，并且周期性边界上的顶点的水平差不能超过2:1，那么这就是必要的。
       *
       */
      virtual bool
      prepare_coarsening_and_refinement() override;

      /**
       * 在处理器之间手动重新划分活动单元。通常这种重新划分会在调用execute_coarsening_and_refinement()（或refine_global()）时自动发生，除非在构造函数中设置
       * @p no_automatic_repartitioning
       * 。设置该标志，然后调用repartition()会得到相同的结果。
       * 如果你想传输数据（使用SolutionTransfer或手动使用register_data_attach()和notify_ready_to_unpack()），你需要设置两次：一次是在调用execute_coarsening_and_refinement()时，它将处理粗化和细化，但显然不会在处理器之间传输任何数据；另一次是在调用repartition()时。
       * 在这里，不会进行粗化和细化，但信息将被打包并运送到不同的处理器。换句话说，在处理数据移动（SolutionTransfer等）方面，你可能希望以与execute_coarsening_and_refinement()相同的方式处理对repartition()的调用。
       * @note  如果没有函数连接到 dealii::Triangulation
       * 类中描述的cell_weight信号，这个函数将平衡每个处理器上的单元格数量。如果连接了一个或多个函数，它将计算权重的总和并平衡各处理器的权重。对权重的唯一要求是每个单元的权重都是正的，而且所有处理器上的所有权重之和可以用64位整数来表示。
       * 除此之外，你可以选择如何解释这些权重。
       * 一个常见的方法是认为权重与在一个单元上进行计算的成本成正比，例如，将装配和求解的时间相加。在实践中，确定这个成本当然不是小事，因为我们不是在孤立的单元上求解，而是在整个网格上。在这种情况下，例如，我们可以选择与每个单元的未知数相等的权重（在hp-finite
       * element方法的背景下），或者使用一种启发式方法来估计每个单元的成本，这取决于，例如，我们是否必须在某些单元而不是其他单元上运行一些昂贵的算法（例如，在装配期间仅在实际处于边界的单元上形成边界积分，或者仅在某些单元而不是其他单元上计算昂贵的非线性项，例如。在
       * step-42 中的弹塑性问题）。)
       *
       */
      void
      repartition();


      /**
       *
       */
      virtual bool
      has_hanging_nodes() const override;

      /**
       * 返回本地的内存消耗，以字节为单位。
       *
       */
      virtual std::size_t
      memory_consumption() const override;

      /**
       * 返回仅包含在p4est数据结构中的本地内存消耗。这已经包含在memory_consumption()中了，但为了调试的目的而单独提供。
       *
       */
      virtual std::size_t
      memory_consumption_p4est() const;

      /**
       * 一个集体的操作，产生一连串的输出文件，这些文件的基名是给定的，包含VTK格式的网格。
       * 更多的是，这个函数对于调试deal.II和p4est之间的接口非常有用。
       *
       */
      void
      write_mesh_vtk(const std::string &file_basename) const;

      /**
       * 生成三角形的检查和。
       * 这是一个集体操作，主要用于调试目的。
       *
       */
      unsigned int
      get_checksum() const;

      /**
       * 将粗略网格的细化信息保存到给定的文件中。这个文件需要在共享的网络文件系统中，从计算的所有节点都可以到达。参见SolutionTransfer类，了解如何将解决方案向量存储到这个文件中。额外的基于单元的数据可以用
       * DistributedTriangulationBase::DataTransfer::register_data_attach().
       * 来保存。
       *
       */
      virtual void
      save(const std::string &filename) const override;

      /**
       * 将用save()保存的细化信息加载回来。在调用此函数之前，网格必须包含与save()中使用的相同的粗略网格。
       * 你不需要用与保存时相同数量的MPI进程来加载。相反，如果加载的网格与保存时使用的MPI进程数量不同，则会对网格进行适当的重新划分。用
       * DistributedTriangulationBase::DataTransfer::register_data_attach()
       * 保存的基于单元的数据，可以在调用 load() 后用
       * DistributedTriangulationBase::DataTransfer::notify_ready_to_unpack()
       * 读入。            如果你使用的p4est版本>0.3.4.2， @p
       * autopartition
       * 标志会告诉p4est忽略三角图保存时的分区，并在加载时使其统一。如果
       * @p autopartition
       * 被设置为false，那么只有在需要时（即遇到不同数量的MPI进程时）才会对三角形进行重新分区。
       *
       */
      virtual void
      load(const std::string &filename,
           const bool         autopartition = true) override;

      /**
       * 从一个给定的平行森林加载细化信息。这个森林可能是通过调用
       * parallel::distributed::Triangulation::get_p4est().
       * 的函数获得的。
       *
       */
      void
      load(const typename dealii::internal::p4est::types<dim>::forest *forest);

      /**
       * 返回粗略单元被移交给p4est的顺序的排列向量。例如，这个向量中
       * $i$ 的元素的值是与p4est管理的 $i$
       * 棵树相对应的deal.II粗糙单元的索引（从begin(0)开始计算）。
       *
       */
      const std::vector<types::global_dof_index> &
      get_p4est_tree_to_coarse_cell_permutation() const;

      /**
       * 返回从粗略的交易单元到p4est树的映射的包络向量。这是get_p4est_tree_to_coarse_cell_permutation的逆过程。
       *
       */
      const std::vector<types::global_dof_index> &
      get_coarse_cell_to_p4est_tree_permutation() const;

      /**
       * 这将返回一个指向内部存储的p4est对象的指针（类型为p4est_t或p8est_t，取决于
       * @p dim).   @warning
       * 如果你修改p4est对象，内部数据结构可能会变得不一致。
       *
       */
      const typename dealii::internal::p4est::types<dim>::forest *
      get_p4est() const;

      /**
       * 除了基类Triangulation中的动作之外，这个函数还将p4est森林中的面连接起来，以获得周期性的边界条件。因此，每对面将最多相差一个细化级别，并且在这些面之间会有鬼魂邻居。
       * 该向量可以由函数  GridTools::collect_periodic_faces.
       * 填充。关于周期性边界条件的更多信息，请参见
       * GridTools::collect_periodic_faces,
       * DoFTools::make_periodicity_constraints  和  step-45  。
       * @note
       * 在使用这个函数之前，必须先初始化三角结构，并且不能进行细化。可以多次调用这个函数，但不推荐。该函数每次调用都会破坏并重建p4est森林。
       *
       */
      virtual void
      add_periodicity(
        const std::vector<dealii::GridTools::PeriodicFacePair<cell_iterator>> &)
        override;


    private:
      /**
       * 存储设置。
       *
       */
      Settings settings;

      /**
       * 一个标志，表示三角图是否有实际内容。
       *
       */
      bool triangulation_has_content;

      /**
       * 一个数据结构，用于保存树之间的连接性。由于每棵树都扎根于一个粗略的网格单元，这个数据结构持有粗略网格单元之间的连接。
       *
       */
      typename dealii::internal::p4est::types<dim>::connectivity *connectivity;

      /**
       * 一个数据结构，用于保存全局三角剖分的局部部分。
       *
       */
      typename dealii::internal::p4est::types<dim>::forest *parallel_forest;

      /**
       * 一个数据结构，用于保存三角形的幽灵单元的一些信息。
       *
       */
      typename dealii::internal::p4est::types<dim>::ghost *parallel_ghost;

      /**
       * 遍历所有p4est树，在私有成员向量local_cell_relations中记录本地拥有的p4est象限和活动的deal.II单元之间的关系。
       * 该向量包含每个本地拥有的p4est象限的活动单元的迭代器，以及描述其关系的CellStatus标志。
       * 存储的向量将按照平行森林的相应本地sc_array中的象限的出现情况排序。因此，这个向量的大小将等于parallel_forest对象中本地拥有的象限的数量。
       * 例如，这些关系将在网格细化过程中建立：在适应parallel_forest之后，但在将这些变化应用于这个三角结构之前，我们将记录单元在细化过程中的变化。有了这些信息，我们就可以相应地准备所有的缓冲区进行数据传输。
       *
       */
      virtual void
      update_cell_relations() override;

      /**
       * 两个数组，存储哪个p4est树对应于哪个粗网格单元，反之亦然。我们需要这些数组，因为p4est在建立森林时，会采用粗格单元的原始顺序，然后在每棵树内应用莫顿排序。但是，如果粗网格单元的顺序不好，这可能意味着存储在本地机器上的森林的各个部分可能被分割到几何上不相近的粗网格单元中。因此，我们根据
       * SparsityTools::reorder_hierarchical()
       * 应用分层预排序，以确保p4est存储的森林部分位于几何上接近的粗略网格单元上。
       *
       */
      std::vector<types::global_dof_index>
        coarse_cell_to_p4est_tree_permutation;
      std::vector<types::global_dof_index>
        p4est_tree_to_coarse_cell_permutation;

      /**
       * 返回一个指向属于给定的dealii_coarse_cell_index()的p4est树的指针
       *
       */
      typename dealii::internal::p4est::types<dim>::tree *
      init_tree(const int dealii_coarse_cell_index) const;

      /**
       * 计算两种数据存储方案之间的互换的函数。
       *
       */
      void
      setup_coarse_cell_to_p4est_tree_permutation();

      /**
       * 把我们新创建的三角形的内容附在上面并复制到p4est数据结构中。
       * 这个函数存在2d和3d的变体。
       *
       */
      void
      copy_new_triangulation_to_p4est(std::integral_constant<int, 2>);
      void
      copy_new_triangulation_to_p4est(std::integral_constant<int, 3>);

      /**
       * 将p4est中的精炼森林的局部部分复制到所附的三角结构中。
       *
       */
      void
      copy_local_forest_to_triangulation();

      /**
       * 内部函数通知所有注册的槽，在重新分区发生之前提供它们的权重。从execute_coarsening_and_refinement()和repartition()调用。
       * @return
       * 一个无符号整数的向量，代表细化/粗化/重新分区周期后每个单元的重量或计算负荷。注意，条目数不需要等于n_active_cells()或n_locally_owned_active_cells()，因为三角形还没有更新。权重是按照p4est在迭代时遇到的顺序排序的。
       *
       */
      std::vector<unsigned int>
      get_cell_weights() const;

      /**
       * 该方法返回一个长度为tria.n_vertices()的位向量，表示一个层面上的局部活动顶点，即在几何多网格中使用的局部拥有的层面单元所接触的顶点（可能包括由于周期性边界条件的顶点）被标记为true。
       * 由 DoFHandler::Policy::ParallelDistributed. 使用。
       *
       */
      std::vector<bool>
      mark_locally_active_vertices_on_level(const int level) const;

      virtual unsigned int
      coarse_cell_id_to_coarse_cell_index(
        const types::coarse_cell_id coarse_cell_id) const override;

      virtual types::coarse_cell_id
      coarse_cell_index_to_coarse_cell_id(
        const unsigned int coarse_cell_index) const override;

      template <int, int, class>
      friend class dealii::FETools::internal::ExtrapolateImplementation;

      template <int, int>
      friend class TemporarilyMatchRefineFlags;
    };


    /**
     * 一般模板的特化，用于1d情况。目前还不支持分布1d三角形。因此，这个类所做的只是抛出一个异常。
     *
     */
    template <int spacedim>
    class Triangulation<1, spacedim>
      : public dealii::parallel::DistributedTriangulationBase<1, spacedim>
    {
    public:
      /**
       * 哑巴设置
       *
       */
      enum Settings
      {
        default_setting                          = 0x0,
        mesh_reconstruction_after_repartitioning = 0x1,
        construct_multigrid_hierarchy            = 0x2
      };

      /**
       * 构造函数。参数表示将用于三角测量的MPI通信器。
       *
       */
      Triangulation(
        const MPI_Comm &mpi_communicator,
        const typename dealii::Triangulation<1, spacedim>::MeshSmoothing
                       smooth_grid = (dealii::Triangulation<1, spacedim>::none),
        const Settings settings    = default_setting);

      /**
       * 解构器。
       *
       */
      virtual ~Triangulation() override;

      /**
       * 返回一个关于粗略单元被移交给p4est的顺序的包络向量。例如，这个向量中的第一个元素i表示分层排序中的第一个单元是由begin(0)开始的第i个交易单元。
       *
       */
      const std::vector<types::global_dof_index> &
      get_p4est_tree_to_coarse_cell_permutation() const;

      /**
       * 这个函数没有实现，但对编译器来说需要存在。
       *
       */
      virtual void
      load(const std::string &filename,
           const bool         autopartition = true) override;

      /**
       * 这个函数没有实现，但对于编译器来说需要存在。
       *
       */
      virtual void
      save(const std::string &filename) const override;

      /**
       * 这个函数没有实现，但对于编译器来说需要存在。
       *
       */
      virtual bool
      is_multilevel_hierarchy_constructed() const override;

      /**
       * 这个函数没有实现，但对于编译器来说需要存在。
       *
       */
      virtual void
      update_cell_relations() override;

      /**
       * 虚数组。这个类是不能用的，但是编译器还是想在几个地方看到这些变量。
       *
       */
      std::vector<types::global_dof_index>
        coarse_cell_to_p4est_tree_permutation;
      std::vector<types::global_dof_index>
        p4est_tree_to_coarse_cell_permutation;

      /**
       * 这个方法只在dim =
       * 2或3的情况下实现，需要一个存根，因为它在dof_handler_policy.cc中使用。
       *
       */
      virtual std::map<unsigned int, std::set<dealii::types::subdomain_id>>
      compute_level_vertices_with_ghost_neighbors(
        const unsigned int level) const;

      /**
       * 和上面一样，这个方法只对dim =
       * 2或3实现，需要一个存根，因为它在dof_handler_policy.cc中使用。
       *
       */
      virtual std::vector<bool>
      mark_locally_active_vertices_on_level(const unsigned int level) const;

      virtual unsigned int
      coarse_cell_id_to_coarse_cell_index(
        const types::coarse_cell_id coarse_cell_id) const override;

      virtual types::coarse_cell_id
      coarse_cell_index_to_coarse_cell_id(
        const unsigned int coarse_cell_index) const override;

      template <int, int>
      friend class TemporarilyMatchRefineFlags;
    };
  } // namespace distributed
} // namespace parallel


#else // DEAL_II_WITH_P4EST

namespace parallel
{
  namespace distributed
  {
    /**
     * 如果我们没有用p4est库实际配置deal.II，编译器就会为并行分布式三角计算选择假类。这个类的存在使得我们可以在整个库中引用
     * parallel::distributed::Triangulation 对象，即使它被禁用。
     * 由于这个类的构造函数被删除，实际上不能创建这样的对象，因为考虑到p4est不可用，这将是毫无意义的。
     *
     */
    template <int dim, int spacedim = dim>
    class Triangulation
      : public dealii::parallel::DistributedTriangulationBase<dim, spacedim>
    {
    public:
      /**
       * 虚拟设置，允许定义被删除的构造函数。
       *
       */
      enum Settings
      {
        default_setting                          = 0x0,
        mesh_reconstruction_after_repartitioning = 0x1,
        construct_multigrid_hierarchy            = 0x2,
        no_automatic_repartitioning              = 0x4
      };

      /**
       * 构造函数。删除是为了确保这种类型的对象不能被构造（也可参见类的文档）。
       *
       */
      explicit Triangulation(
        const MPI_Comm &  /*mpi_communicator*/ ,
        const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
         /*smooth_grid*/ 
        = (dealii::Triangulation<dim, spacedim>::none),
        const Settings  /*settings*/  = default_setting) = delete;

      /**
       * 哑巴替换，以便在编译该类时有更好的错误信息。
       *
       */
      virtual bool
      is_multilevel_hierarchy_constructed() const override
      {
        return false;
      }

      /**
       * 哑巴替换，以便在编译这个类时有更好的错误信息。
       *
       */
      virtual void
      save(const std::string &  /*filename*/ ) const override
      {}

      /**
       * 虚伪的替换，以便在编译这个类时允许更好的错误信息。
       *
       */
      virtual void
      load(const std::string &  /*filename*/ ,
           const bool  /*autopartition*/  = true) override
      {}

      /**
       * 虚伪的替换，以便在编译这个类时允许更好的错误信息。
       *
       */
      virtual void
      update_cell_relations() override
      {}
    };
  } // namespace distributed
} // namespace parallel


#endif


namespace parallel
{
  namespace distributed
  {
    /**
     * 这个类暂时修改了所有活动单元的细化和粗化标志，以匹配p4est神谕。
     * 该修改只发生在 parallel::distributed::Triangulation
     * 对象上，并且在该类的实例化过程中持续存在。
     * TemporarilyMatchRefineFlags类应该只与
     * Triangulation::Signals::post_p4est_refinement
     * 信号结合使用。在这个阶段，p4est神谕已经被细化了，但是三角测量仍然没有改变。在修改之后，所有的refine和coarsen标志都描述了traingulation实际上将如何被细化。
     * 这个类的使用在  step-75  中得到了证明。
     *
     */
    template <int dim, int spacedim = dim>
    class TemporarilyMatchRefineFlags : public Subscriptor
    {
    public:
      /**
       * 构造函数。
       * 存储所有活动单元的细化和粗化标志，如果提供的三角图是类型
       * parallel::distributed::Triangulation.
       * 调整它们以与p4est神谕一致。
       *
       */
      TemporarilyMatchRefineFlags(dealii::Triangulation<dim, spacedim> &tria);

      /**
       * 解构器。            将 parallel::distributed::Triangulation
       * 上所有活动单元的细化和粗化标志返回到它们之前的状态。
       *
       */
      ~TemporarilyMatchRefineFlags();

    private:
      /**
       * 修改后的 parallel::distributed::Triangulation. 。
       *
       */
      const SmartPointer<
        dealii::parallel::distributed::Triangulation<dim, spacedim>>
        distributed_tria;

      /**
       * 一个向量，用于临时存储在
       * parallel::distributed::Triangulation.
       * 上被修改过的细化标志。
       *
       */
      std::vector<bool> saved_refine_flags;

      /**
       * 一个向量，在 parallel::distributed::Triangulation.
       * 上的粗化标志被修改之前，临时存储这些标志。
       *
       */
      std::vector<bool> saved_coarsen_flags;
    };
  } // namespace distributed
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif


