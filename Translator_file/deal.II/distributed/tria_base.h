//include/deal.II-translator/distributed/tria_base_0.txt
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

#ifndef dealii_distributed_tria_base_h
#define dealii_distributed_tria_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/partitioner.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/grid/tria.h>

#include <functional>
#include <list>
#include <set>
#include <utility>
#include <vector>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  /**
   * 这个类描述了所有并行工作的三角化类的接口，即
   * parallel::distributed::Triangulation,
   * parallel::fullydistributed::Triangulation, 和
   * parallel::shared::Triangulation.
   * 因此，它是一个可以用来测试三角化对象的引用指针是否指的是顺序三角化，或者三角化是否真的是并行的类。换句话说，人们可以写一个这样的函数。
   * @code
   * template <int dim, int spacedim>
   * bool is_parallel (const dealii::Triangulation<dim,spacedim> &tria)
   * {
   *   if (dynamic_cast<const parallel::TriangulationBase<dim,spacedim>*>
   *                   (&tria)
   *       != nullptr)
   *     return true;
   *   else
   *     return false;
   * }
   * @endcode
   *
   */
  template <int dim, int spacedim = dim>
  class TriangulationBase : public dealii::Triangulation<dim, spacedim>
  {
  public:
    /**
     * 构造函数。
     *
     */
    TriangulationBase(
      const MPI_Comm &mpi_communicator,
      const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
                 smooth_grid = (dealii::Triangulation<dim, spacedim>::none),
      const bool check_for_distorted_cells = false);

    /**
     * 解构器。
     *
     */
    virtual ~TriangulationBase() override;

    /**
     * 返回该三角化所使用的MPI通信器。
     *
     */
    virtual MPI_Comm
    get_communicator() const override;

    /**
     * 如果支持多级层次结构并且已经构建完成，则返回。
     *
     */
    virtual bool
    is_multilevel_hierarchy_constructed() const = 0;

    /**
     * 实现与基类中相同的函数。
     * @note
     * 这个函数复制了源三角结构的单元，但没有复制通讯器。换句话说，产生的三角形将在它所构建的通信器上操作。
     *
     */
    virtual void
    copy_triangulation(
      const dealii::Triangulation<dim, spacedim> &old_tria) override;

    /**
     * 返回三角形中本地拥有的活动单元的数量，即这些单元的subdomain_id等于local_owned_subdomain()。请注意，在目前的处理器上存储的三角形中可能有更多的活动单元，例如幽灵单元，或者离本地拥有的单元块更远的单元，但这是为了确保存储这个处理器的活动单元集的三角形仍然保持相邻单元的2:1的大小比例而需要的。
     * 由于上面的备注，这个函数的结果总是小于或等于::三角形基类中同名函数的结果，其中包括活动的幽灵和人工单元（也见
     * @ref GlossArtificialCell 和 @ref GlossGhostCell  ）。
     *
     */
    unsigned int
    n_locally_owned_active_cells() const;

    /**
     * 返回所有处理器中每个处理器所拥有的活动单元的数量之和。这等于三角形中活动单元的总数量。
     *
     */
    virtual types::global_cell_index
    n_global_active_cells() const override;

    /**
     * 返回以字节为单位的本地内存消耗。
     *
     */
    virtual std::size_t
    memory_consumption() const override;


    /**
     * 返回全局最大水平。如果当前处理器只将单元存储在域中不是很细的部分，但如果其他处理器将单元存储在域中更深的细化部分，这个数字可能大于
     * dealii::Triangulation::n_levels()
     * （这个类的基类中的一个函数）的返回值。
     *
     */
    virtual unsigned int
    n_global_levels() const override;

    /**
     * 返回那些由当前处理器拥有的单元格的子域ID。三角形中所有没有这个子域ID的单元格要么被其他处理器拥有，要么有只存在于其他处理器的子域。
     *
     */
    types::subdomain_id
    locally_owned_subdomain() const override;

    /**
     * 返回一组处理器的MPI行列，这些处理器至少有一个与本地处理器的单元格相邻的幽灵单元。换句话说，这是所有幽灵单元的subdomain_id()的集合。
     * 返回的集合是对称的，即如果 @p i 包含在处理器 @p j,
     * 的列表中，那么 @p j 也将包含在处理器 @p i.
     * 的列表中。
     *
     */
    const std::set<types::subdomain_id> &
    ghost_owners() const;

    /**
     * 返回一组MPI等级的处理器，这些处理器至少有一个与我们在几何多网格中使用的单元相邻的水平重影单元。换句话说，这是所有层面鬼单元的level_subdomain_id()的集合。
     * 返回的集合是对称的，即如果 @p i 包含在处理器 @p j,
     * 的列表中，那么 @p j 也将包含在处理器 @p i.
     * 的列表中。
     * @note
     * 只有在多网格所有权被分配的情况下（通过在构造时设置construct_multigrid_hierarchy标志），才能确定级别ghost所有者，否则返回的集合将为空。
     *
     */
    const std::set<types::subdomain_id> &
    level_ghost_owners() const;

    /**
     * 返回三角形活动层上的单元格的全局索引的分区器。
     *
     */
    const std::weak_ptr<const Utilities::MPI::Partitioner>
    global_active_cell_index_partitioner() const;

    /**
     * 返回三角形给定 @p
     * 层上的单元格的全局索引的分区器。
     *
     */
    const std::weak_ptr<const Utilities::MPI::Partitioner>
    global_level_cell_index_partitioner(const unsigned int level) const;

    /**
     * 返回一个地图，对于每个顶点，列出其子域与该顶点相邻的所有处理器。
     * @deprecated  用 GridTools::compute_vertices_with_ghost_neighbors()
     * 代替
     * parallel::TriangulationBase::compute_vertices_with_ghost_neighbors().
     * 。
     *
     */
    DEAL_II_DEPRECATED virtual std::map<unsigned int,
                                        std::set<dealii::types::subdomain_id>>
    compute_vertices_with_ghost_neighbors() const;

    /**
     * @copydoc   dealii::Triangulation::get_boundary_ids()  *
     * dealii::Triangulation::get_boundary_ids() 。
     * @note
     * 这个函数涉及到一个全局通信，收集所有进程的当前ID。
     *
     */
    virtual std::vector<types::boundary_id>
    get_boundary_ids() const override;

    /**
     * @copydoc   dealii::Triangulation::get_manifold_ids() .
     * @note
     * 这个函数涉及一个全局通信，收集所有进程的所有当前ID。
     *
     */
    virtual std::vector<types::manifold_id>
    get_manifold_ids() const override;

    /**
     * 当顶点在本地被移动时，例如使用如下代码
     * @code
     * cell->vertex(0) = new_location;
     * @endcode
     * 那么这个函数可以用来更新MPI进程之间顶点的位置。
     * 所有已经被移动的顶点和可能在进程的幽灵层中的顶点都必须在
     * @p vertex_locally_moved
     * 参数中报告。这确保了必须在进程之间发送的那部分信息被实际发送。此外，很重要的一点是，在进程之间的边界上的顶点正好在一个进程中被报告（例如，具有最高id的那个）。
     * 否则，如果多个进程以不同的方式移动一个顶点，我们可以期待不理想的结果。一个典型的策略是让处理器
     * $i$
     * 移动那些与单元相邻的顶点，这些单元的所有者包括处理器
     * $i$ ，但没有其他处理器 $j$ 与 $j<i$
     * ；换句话说，对于子域边界的顶点，子域id最低的处理器
     * "拥有 "一个顶点。
     * @note
     * 只有移动位于本地拥有的单元上的顶点或位于幽灵层上的单元才有意义。这是因为你可以确定这些顶点确实存在于所有处理器汇总的最精细的网格上，而位于人工单元但至少不在幽灵层的顶点可能存在于全局最精细的网格上，也可能不存在。因此，
     * @p vertex_locally_moved
     * 参数可能不包含至少不在幽灵单元上的顶点。
     * @note
     * 这个函数移动顶点的方式是，在每个处理器上，每个本地拥有的和鬼魂单元的顶点与其他处理器上这些单元的相应位置是一致的。另一方面，人工单元的位置一般会是错误的，因为人工单元在其他处理器上可能存在也可能不存在，因此不可能以任何方式确定它们的位置。这通常不是一个问题，因为人们从不在人工单元上做任何事情。但是，如果在以后的步骤中对移动顶点的网格进行细化，可能会导致问题。
     * 如果你想这样做，正确的方法是保存应用于每个顶点的偏移量，调用这个函数，并在细化或粗化网格之前应用相反的偏移量，并再次调用这个函数。
     * @param  vertex_locally_moved
     * 表示哪些顶点被移动的位图。这个数组的大小必须等于
     * Triangulation::n_vertices()  并且必须是
     * GridTools::get_locally_owned_vertices().   @see
     * 标记的顶点的一个子集。这个函数用于，例如，在
     * GridTools::distort_random().  中。
     *
     */
    void
    communicate_locally_moved_vertices(
      const std::vector<bool> &vertex_locally_moved);

  protected:
    /**
     * 将用于三角测量的MPI通信器。我们为这个类创建一个唯一的通信器，它是传递给构造函数的通信器的复制品。
     *
     */
    const MPI_Comm mpi_communicator;

    /**
     * 将用于当前处理器的子域ID。这就是MPI等级。
     *
     */
    types::subdomain_id my_subdomain;

    /**
     * 子域的总数（或MPI通信器的大小）。
     *
     */
    types::subdomain_id n_subdomains;

    /**
     * 一个包含分布式三角测量信息的结构。
     *
     */
    struct NumberCache
    {
      /**
       * 这个MPI等级的本地拥有的活动单元的数量。
       *
       */
      unsigned int n_locally_owned_active_cells;
      /**
       * 活动单元的总数（ @p  n_locally_owned_active_cells之和）。
       *
       */
      types::global_cell_index n_global_active_cells;
      /**
       * 全局级别数，计算为所有MPI级别的最大级别数，所以<tt>n_levels()<=n_global_levels
       * = max(n_levels() on proc i)</tt>。
       *
       */
      unsigned int n_global_levels;
      /**
       * 一个包含该处理器上幽灵单元所有者的子域_id（MPI等级）的集合。
       *
       */
      std::set<types::subdomain_id> ghost_owners;
      /**
       * 一个包含此处理器上的级别幽灵单元所有者的MPI等级的集合（对于所有级别）。
       *
       */
      std::set<types::subdomain_id> level_ghost_owners;

      /**
       * 全局活动单元索引的分区器。
       *
       */
      std::shared_ptr<const Utilities::MPI::Partitioner>
        active_cell_index_partitioner;

      /**
       * 每个级别的全局级别单元索引的分区器。
       *
       */
      std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
        level_cell_index_partitioners;

      NumberCache();
    };

    NumberCache number_cache;

    /**
     * 在网格创建或细化后更新number_cache变量。
     *
     */
    virtual void
    update_number_cache();

    /**
     * @copydoc   dealii::Triangulation::update_reference_cells() .
     *
     */
    void
    update_reference_cells() override;

    /**
     * 重置全局活动单元索引和全局水平单元索引。
     *
     */
    void
    reset_global_cell_indices();
  };



  /**
   * @code
   * template <int dim, int spacedim>
   * bool
   * is_parallel_distributed(const dealii::Triangulation<dim,spacedim> &tria)
   * {
   *   if(dynamic_cast<const
   *                   parallel::DistributedTriangulationBase<dim,spacedim>*>
   *                  (&tria)
   *      != nullptr)
   *     return true;
   *   else
   *     return false;
   * }
   * @endcode
   *
   *
   */
  template <int dim, int spacedim = dim>
  class DistributedTriangulationBase
    : public dealii::parallel::TriangulationBase<dim, spacedim>
  {
  public:
    /**
     * 构造器。
     *
     */
    DistributedTriangulationBase(
      const MPI_Comm &mpi_communicator,
      const typename dealii::Triangulation<dim, spacedim>::MeshSmoothing
                 smooth_grid = (dealii::Triangulation<dim, spacedim>::none),
      const bool check_for_distorted_cells = false);

    /**
     * 通过删除所有的数据，将这个三角剖面重设为处女状态。
     * 请注意，只有当这个对象的订阅不再存在时，才允许这个操作，比如使用它的DoFHandler对象。
     *
     */
    virtual void
    clear() override;

    using cell_iterator =
      typename dealii::Triangulation<dim, spacedim>::cell_iterator;

    using CellStatus =
      typename dealii::Triangulation<dim, spacedim>::CellStatus;

    /**
     * 将三角图保存到给定的文件中。这个文件需要在共享的网络文件系统上，从计算中的所有节点都可以到达。参见SolutionTransfer类，了解如何将解决方案向量存储到这个文件中。其他基于单元的数据可以用register_data_attach()来保存。
     *
     */
    virtual void
    save(const std::string &filename) const = 0;

    /**
     * 将用save()保存的三角图装回去。用register_data_attach()保存的基于单元的数据可以在调用load()后用notify_ready_to_unpack()读入。
     *
     */
    virtual void
    load(const std::string &filename, const bool autopartition = true) = 0;

    /**
     * 注册一个函数，可以用来将固定大小的数据附加到单元格上。这在两个方面是有用的。(i)
     * 在细化和粗化三角形时（ @a  例如在
     * parallel::distributed::Triangulation::execute_coarsening_and_refinement()),
     * 中，需要能够为每个单元存储一个或多个数据向量，以描述该单元上的求解值，这样当网格被重新划分时，这些数据可以被转移到该单元（或其父/子）的新拥有者处理器；(ii)
     * 当将计算序列化到一个文件时，有必要将数据附加到单元，以便能够被保存（
     * @a  ]，例如在 parallel::distributed::Triangulation::save())
     * 中，与单元的其他信息一起，如果有必要，以后可以从磁盘上重新加载，在处理器中对单元进行不同的细分。
     * 这个功能的工作方式是，它允许任何数量的利益方注册他们的意图，将数据附加到单元中。做到这一点的类的一个例子是
     * parallel::distributed::SolutionTransfer ，每个
     * parallel::distributed::SolutionTransfer
     * 对象在当前Triangulation对象上工作，然后需要注册其意图。
     * 每一方都注册了一个回调函数（这里的第一个参数，
     * @p pack_callback)
     * ，每当三角化的execute_coarsening_and_refinement()或save()函数被调用时，都会被调用。
     * 然后，当前函数返回一个整数句柄，对应于这里提供的回调将附加的数据集的数量。
     * 虽然这个数字可以被赋予精确的含义，但这并不重要：除了把它返回给notify_ready_to_unpack()函数外，你实际上永远不必对这个数字做任何事情。
     * 换句话说，每个感兴趣的人（即当前函数的调用者）需要在提供给notify_ready_to_unpack()的回调中存储他们各自返回的句柄，以便以后在解包数据时使用。
     * 每当 @p pack_callback
     * 再被execute_coarsening_and_refinement()或load()在给定的单元上调用时，它都会收到一些参数。特别是，传递给回调的第一个参数表示它应该附加数据的单元。这始终是一个活动单元。
     * 第二个参数，即CellStatus，提供给回调函数的参数将告诉你给定的单元是否会被粗化、细化，或保持原样。这个状态可能与设置在该单元上的细化或粗化标志不同，以适应诸如
     * "每条边一个悬空节点
     * "的规则）。这些标志需要在它们所属的p4est象限的背景下阅读，因为它们的关系被收集在local_cell_relations。
     * 具体来说，这个参数的值意味着以下几点。
     *
     *
     *
     *
     *
     *
     * - `cell_persist`: 单元不会被细化/粗化，但可能被转移到不同的处理器。如果是这种情况，回调将希望把这个单元的数据打包成一个数组，并存储在所提供的地址上，以便以后在这个单元可能出现的地方进行解包。
     *
     *
     *
     *
     *
     *
     * - `cell_refine`: 这个单元格将被细化为4个或8个单元格（分别在2d和3d中）。然而，由于这些子单元还不存在，所以在调用回调的时候，你无法访问它们。因此，在local_cell_relations中，相应的子单元的p4est象限被链接到将被精炼的deal.II单元。具体来说，只有第一个子单元被标记为 "CELL_REFINE"，而其他单元将被标记为 "CELL_INVALID"，这表明这些单元在打包或拆包过程中会被默认忽略。这确保了数据只被传输到父单元或从父单元传输一次。如果回调是以`CELL_REFINE`调用的，回调将希望把这个单元格上的数据打包成一个数组，并存储在提供的地址上，以便以后解包的方式，这样就可以把数据传输到该单元格的子单元，然后就可以使用了。换句话说，如果回调想要打包的数据对应于一个有限元字段，那么从父单元到（新）子单元的延长必须在解包时发生。
     *
     *
     *
     *
     *
     *
     * - `cell_coarsen`: 这个单元格的子代将被粗化为给定的单元格。这些子单元仍然存在，所以如果这是给回调的第二个参数的值，回调将想把数据从子单元转移到当前的父单元，并将其打包，以便以后可以在一个不再有任何子单元（也可能位于不同的处理器上）上再次解包。换句话说，如果回调想要打包的数据对应于一个有限元字段，那么它就需要在这一点上做从子单元到父单元的限制。
     *
     *
     *
     *
     *
     *
     * - `cell_invalid`: 参见 `CELL_REFINE`.
     * @note
     * 如果这个函数用于使用save()和load()的数据序列化，那么调用回调函数的单元状态参数将总是`CELL_PERSIST`。
     * 回调函数预计将返回一个格式为 `std::vector<char>`,
     * 的内存块，代表某个单元上的打包数据。
     * 第二个参数 @p returns_variable_size_data
     * 表示回调函数返回的内存区域大小是否因单元而异（<tt>=true</tt>）或在整个域中每个单元保持不变（<tt>=false</tt>）。
     * @note
     * 这个函数的目的是为了注册附加数据的意图，以便随后调用execute_coarsening_and_refinement()和notify_ready_to_unpack()、save()、load()。因此，一旦这些回调被调用，notify_ready_to_unpack()、save()和load()都会忘记已注册的回调，如果你想让它们在再次调用这些函数时被激活，你必须用三角法重新注册它们。
     *
     */
    unsigned int
    register_data_attach(
      const std::function<std::vector<char>(const cell_iterator &,
                                            const CellStatus)> &pack_callback,
      const bool returns_variable_size_data);

    /**
     * 这个函数与register_data_attach()相反。它被称为 <i>after</i>
     * execute_coarsening_and_refinement() 或 save()/load()
     * 函数，当之前将数据附加到三角结构上的类和函数准备好接收这些数据时，这些数据将被传送到其他处理器，跨越网格细化，或将数据序列化到文件。这个过程的重要部分是三角形不能在execute_coarsening_and_refinement()或load()结束后，通过先前附加的回调函数（如register_data_attach()函数）立即完成这个过程，因为最终想要回数据的类可能需要在重新创建网格的时间点和实际接收数据的时间点之间做一些设置。
     * 一个例子是 parallel::distributed::SolutionTransfer
     * 类，它只能在网格在当前处理器上完全可用后才能接收数据，而且只能在DoFHandler被重新初始化和分配自由度后才能接收。换句话说，在可以接收附加在单元上的数据的类准备好这样做之前，通常需要在用户空间进行大量的设置。
     * 当他们准备好时，他们会使用当前函数告诉三角测量对象现在是他们准备好的时候，调用当前函数。
     * 然后为每个新的本地拥有的单元调用所提供的回调函数。回调的第一个参数是一个指定单元的迭代器；第二个参数表示有关单元的状态；第三个参数通过两个迭代器定位一个内存区域，其中包含之前从提供给
     * register_data_attach() 的回调中保存的数据。
     * CellStatus将指示该单元是否被细化、粗化或持久化而未被改变。然后，回调的
     * @p cell_iterator
     * 参数将是一个活跃的、本地拥有的单元（如果该单元没有被精炼），或者是直接的父单元，如果它在execute_coarsening_and_refinement()期间被精炼。
     * 因此，与register_data_attach()期间相反，如果状态为`CELL_REFINE`，你现在可以访问子单元，但对于状态为`CELL_COARSEN`的回调，就不能再访问了。
     * 这个函数的第一个参数 "handle "对应于
     * register_data_attach()
     * 的返回值。(这个句柄的数字值应该代表什么的确切含义并不重要，你也不应该试图用它来做任何事情，除了在调用register_data_attach()和相应的调用notify_ready_to_unpack()之间传输信息。)
     *
     */
    void
    notify_ready_to_unpack(
      const unsigned int handle,
      const std::function<
        void(const cell_iterator &,
             const CellStatus,
             const boost::iterator_range<std::vector<char>::const_iterator> &)>
        &unpack_callback);

  protected:
    /**
     * 将额外的细胞附加数据保存到给定的文件中。第一个参数用于确定将缓冲区写入何处的偏移量。
     * 由  @ref save  调用。
     *
     */
    void
    save_attached_data(const unsigned int global_first_cell,
                       const unsigned int global_num_cells,
                       const std::string &filename) const;

    /**
     * 从给定的文件中加载额外的细胞附加数据，如果有任何保存的话。
     * 第一个参数用于确定从哪里读取缓冲区的偏移量。
     * 由  @ref load  调用。
     *
     */
    void
    load_attached_data(const unsigned int global_first_cell,
                       const unsigned int global_num_cells,
                       const unsigned int local_num_cells,
                       const std::string &filename,
                       const unsigned int n_attached_deserialize_fixed,
                       const unsigned int n_attached_deserialize_variable);

    /**
     * 一个记录当前活动单元的CellStatus的函数，这些单元为本地所有。这个信息对于在适应或序列化过程中在网格之间传输数据是必须的，例如，使用
     * parallel::distributed::SolutionTransfer.
     * 关系将被存储在私有成员local_cell_relations。关于CellStatus的广泛描述，请参见成员函数register_data_attach()的文档。
     *
     */
    virtual void
    update_cell_relations() = 0;

    /**
     * 用于将CellStatus分配给deal.II单元格迭代器的辅助数据结构。关于前者的广泛描述，请参见成员函数register_data_attach()的文档。
     *
     */
    using cell_relation_t = typename std::pair<cell_iterator, CellStatus>;

    /**
     * 对的向量，每个都包含一个deal.II单元格迭代器和其各自的CellStatus。要更新其内容，请使用
     * update_cell_relations() 成员函数。
     *
     */
    std::vector<cell_relation_t> local_cell_relations;

    /**
     * 一个结构，用于存储已经或将要通过register_data_attach()函数附加到单元格上的数据信息，以后通过notify_ready_to_unpack()检索。
     *
     */
    struct CellAttachedData
    {
      /**
       * 通过register_data_attach()函数附加到三角结构上的函数数量，例如SolutionTransfer。
       *
       */
      unsigned int n_attached_data_sets;

      /**
       * 从load()调用后需要解开数据的函数的数量
       *
       */
      unsigned int n_attached_deserialize;

      using pack_callback_t = std::function<std::vector<char>(
        typename dealii::Triangulation<dim, spacedim>::cell_iterator,
        typename dealii::Triangulation<dim, spacedim>::CellStatus)>;

      /**
       * 这些回调函数将按照它们在register_data_attach()函数中注册的顺序来存储。
       *
       */
      std::vector<pack_callback_t> pack_callbacks_fixed;
      std::vector<pack_callback_t> pack_callbacks_variable;
    };

    CellAttachedData cell_attached_data;

    /**
     * 这个在 parallel::DistributedTriangulationBase
     * 的私有范围内的类专门用于跨重新分区的网格和向/从文件系统的数据传输。
     * 它被设计用来存储所有用于传输的数据缓冲区。
     *
     */
    class DataTransfer
    {
    public:
      DataTransfer(const MPI_Comm &mpi_communicator);

      /**
       * 通过调用 @p cell_relations.
       * 中每个单元的打包回调函数来准备数据传输。 @p
       * pack_callbacks_fixed
       * 中所有注册的回调函数将写入固定大小的缓冲区，而
       * @p pack_callbacks_variable
       * 的每个条目将把其数据写入可变大小的缓冲区。
       *
       */
      void
      pack_data(const std::vector<cell_relation_t> &cell_relations,
                const std::vector<typename CellAttachedData::pack_callback_t>
                  &pack_callbacks_fixed,
                const std::vector<typename CellAttachedData::pack_callback_t>
                  &pack_callbacks_variable);



      /**
       * 解除 @p cell_relations.
       * 的每个条目上的CellStatus信息，数据必须事先用execute_transfer()传输或通过load()从文件系统读取。
       *
       */
      void
      unpack_cell_status(std::vector<cell_relation_t> &cell_relations) const;

      /**
       * 用提供的 @p unpack_callback 函数在 @p cell_relations
       * 中注册的每个单元上解压先前传输的数据。
       * 参数 @p handle 对应于允许 @p unpack_callback
       * 函数从存储器中读取的位置。它的值需要与之前注册的相应的pack_callback函数一致。
       * 数据必须事先用execute_transfer()传输，或通过load()从文件系统中读取。
       *
       */
      void
      unpack_data(
        const std::vector<cell_relation_t> &cell_relations,
        const unsigned int                  handle,
        const std::function<void(
          const typename dealii::Triangulation<dim, spacedim>::cell_iterator &,
          const typename dealii::Triangulation<dim, spacedim>::CellStatus &,
          const boost::iterator_range<std::vector<char>::const_iterator> &)>
          &unpack_callback) const;

      /**
       * 向文件系统传输数据。
       * 数据将被写入一个单独的文件，该文件的名称由 @p
       * filename
       * 和一个附加的标识符<tt>_fixed.data</tt>，用于固定大小的数据，<tt>_variable.data</tt>用于可变大小的数据。
       * 所有处理器通过MPIIO同时向这些文件写入。
       * 每个处理器要写入的位置将由提供的输入参数决定。
       * 数据必须事先用pack_data()打包。
       *
       */
      void
      save(const unsigned int global_first_cell,
           const unsigned int global_num_cells,
           const std::string &filename) const;

      /**
       * 从文件系统传输数据。
       * 数据将从单独的文件中读取，该文件的名称由 @p
       * filename
       * 和一个附加的标识符<tt>_fixed.data</tt>组成，用于固定大小的数据，<tt>_variable.data</tt>用于可变大小数据。
       * 需要 @p n_attached_deserialize_fixed 和 @p
       * n_attached_deserialize_variable
       * 参数来收集每个回调的内存偏移量。
       * 所有处理器通过MPIIO同时从这些文件中读取。
       * 每个处理器要读取的位置将由提供的输入参数决定。
       * 在加载之后，需要调用unpack_data()来最终将数据分配到相关的三角地带。
       *
       */
      void
      load(const unsigned int global_first_cell,
           const unsigned int global_num_cells,
           const unsigned int local_num_cells,
           const std::string &filename,
           const unsigned int n_attached_deserialize_fixed,
           const unsigned int n_attached_deserialize_variable);

      /**
       * 清除所有容器和相关数据，并将成员值重置为默认状态。
       * 完全释放了内存。
       *
       */
      void
      clear();

      /**
       * 表示可变大小的数据是否被打包的标志。
       *
       */
      bool variable_size_data_stored;

      /**
       * 那些调用register_data_attach()的函数想要附加到每个单元的累计大小，以字节为单位。这个数字只与固定大小的缓冲区有关，其中附加到每个单元的数据具有完全相同的大小。
       * 这个容器的最后一个条目对应于固定大小的缓冲区中每个单元打包的数据大小（可以调用<tt>sizes_fixed_cumulative.back()</tt>来访问）。
       *
       */
      std::vector<unsigned int> sizes_fixed_cumulative;

      /**
       * 为p4est的固定尺寸传输函数设计的连续缓冲区。
       *
       */
      std::vector<char> src_data_fixed;
      std::vector<char> dest_data_fixed;

      /**
       * 为p4est的可变大小传输函数设计的连续缓冲区。
       *
       */
      std::vector<int>  src_sizes_variable;
      std::vector<int>  dest_sizes_variable;
      std::vector<char> src_data_variable;
      std::vector<char> dest_data_variable;

    private:
      MPI_Comm mpi_communicator;
    };

    DataTransfer data_transfer;
  };

} // namespace parallel

DEAL_II_NAMESPACE_CLOSE

#endif


