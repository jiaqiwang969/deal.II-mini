//include/deal.II-translator/particles/particle_handler_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2021 by the deal.II authors
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

#ifndef dealii_particles_particle_handler_h
#define dealii_particles_particle_handler_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/bounding_box.h>
#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/particles/partitioner.h>
#include <deal.II/particles/property_pool.h>

#include <boost/range/iterator_range.hpp>
#include <boost/serialization/map.hpp>

DEAL_II_NAMESPACE_OPEN

namespace Particles
{
  /**
   * 这个类管理粒子的存储和处理。它提供了有效存储粒子所需的数据结构、迭代粒子和查找粒子的访问函数，以及在并行域中分配粒子的算法。请注意，该类的设计方式与三角化类相似。特别是，我们把本地进程域中的粒子称为本地粒子，而把属于邻居进程并生活在本地拥有的域周围的幽灵单元中的粒子称为
   * "幽灵粒子"。    这个类在  step-70  中使用。
   * @ingroup Particle
   *
   */
  template <int dim, int spacedim = dim>
  class ParticleHandler : public Subscriptor
  {
  public:
    /**
     * 一个可用于遍历域中所有粒子的类型。
     *
     */
    using particle_iterator = ParticleIterator<dim, spacedim>;

    /**
     * 一个代表粒子范围的类型。
     *
     */
    using particle_iterator_range = boost::iterator_range<particle_iterator>;

    /**
     * 默认的构造函数。
     *
     */
    ParticleHandler();

    /**
     * 构造函数，用一个给定的三角形和映射来初始化粒子处理程序。由于粒子是相对于它们周围的单元而被存储的，所以这个信息对于正确组织粒子集合是必要的。
     * 这个构造函数等同于调用默认构造函数和初始化函数。
     *
     */
    ParticleHandler(const Triangulation<dim, spacedim> &tria,
                    const Mapping<dim, spacedim> &      mapping,
                    const unsigned int                  n_properties = 0);

    /**
     * 解构函数。
     *
     */
    virtual ~ParticleHandler() override = default;

    /**
     * 初始化粒子处理程序。这个函数并不清除内部数据结构，它只是设置三角形和要使用的映射。
     *
     */
    void
    initialize(const Triangulation<dim, spacedim> &tria,
               const Mapping<dim, spacedim> &      mapping,
               const unsigned int                  n_properties = 0);

    /**
     * 将粒子处理程序 @p particle_handler
     * 的状态复制到当前对象中。这将复制所有的粒子和属性，并使这个对象成为一个与
     * @p particle_handler.
     * 相同的副本，这个对象中现有的粒子被删除。请注意，这不会复制与
     * @p particle_handler,
     * 的信号相连的函数，也不会将当前对象的成员函数与三角信号相连，如果有必要，必须由调用者完成，也就是说，如果
     * @p particle_handler 有相连的函数。
     * 这个函数很昂贵，因为它必须复制 @p particle_handler,
     * 中的所有数据并插入这个对象中，这可能是一个相当大的数据量。然而，它可以用来保存某个时间点的粒子集合的状态，并在以后的某些条件下重置这个状态，例如，如果一个时间步骤必须被撤销和重复，那么它就很有用。
     *
     */
    void
    copy_from(const ParticleHandler<dim, spacedim> &particle_handler);

    /**
     * 清除所有与粒子有关的数据。
     *
     */
    void
    clear();

    /**
     * 只清除粒子数据，但保留关于粒子数量的缓存信息。这在进程之间重组粒子数据时很有用。
     *
     */
    void
    clear_particles();

    /**
     * 更新所有内部缓存的数字。请注意，所有修改内部数据结构并作用于多个粒子的函数将自动调用这个函数（例如insert_particles），而作用于单个粒子的函数将不调用这个函数（例如insert_particle）。这样做是因为与单次操作相比，更新是很昂贵的。
     *
     */
    void
    update_cached_numbers();

    /**
     * 返回一个到第一个粒子的迭代器。
     *
     */
    particle_iterator
    begin() const;

    /**
     * 返回一个到第一个粒子的迭代器。
     *
     */
    particle_iterator
    begin();

    /**
     * 返回一个超过粒子末端的迭代器。
     *
     */
    particle_iterator
    end() const;

    /**
     * 返回一个超过粒子末端的迭代器。
     *
     */
    particle_iterator
    end();

    /**
     * 返回一个到第一个幽灵粒子的迭代器。
     *
     */
    particle_iterator
    begin_ghost() const;

    /**
     * 返回一个到第一个鬼魂粒子的迭代器。
     *
     */
    particle_iterator
    begin_ghost();

    /**
     * 返回一个超过鬼魂粒子末端的迭代器。
     *
     */
    particle_iterator
    end_ghost() const;

    /**
     * 返回一个超过鬼魂粒子末端的迭代器。
     *
     */
    particle_iterator
    end_ghost();

    /**
     * 返回住在给定单元上的粒子的数量。
     * @note 虽然这个函数在 step-19
     * 中使用，但如果粒子的数量很大，它并不是一个有效的函数。这是因为要找到位于一个单元中的粒子需要花费
     * ${\cal O}(\log N)$ ，其中 $N$
     * 是整体粒子的数量。由于你可能会对每个单元都这样做，并假设粒子的数量和单元的数量大致成正比，你最终会得到一个
     * ${\cal O}(N \log N)$
     * 的算法。一个更好的方法是利用这样一个事实：在内部，粒子是按照它们所在的活动单元的顺序排列的。换句话说，如果你遍历所有的粒子，当你在活动单元上行走时，你会以同样的顺序遇到它们。你可以利用这一点，保持对第一个单元的第一个粒子的迭代器，当你移动到下一个单元时，你也会递增粒子迭代器，直到你找到位于下一个单元的粒子。计算这个过程花了多少步，就可以得到你要找的数字，在所有单元中累积起来的代价是
     * ${\cal O}(\log N)$ 。      例如，这就是 step-70
     * 中使用的方法。这种方法在 step-19 的
     * "扩展部分的可能性 "中也有详细说明。
     *
     */
    types::particle_index
    n_particles_in_cell(
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell)
      const;

    /**
     * 返回一对粒子迭代器，标记特定单元中粒子的开始和结束。最后一个迭代器是已不在单元格中的第一个粒子。
     * 返回范围内的元素数量等于n_particles_in_cell()函数返回的内容。
     * @note  虽然这个函数在 step-19
     * 中被使用，但如果粒子的数量很大，它并不是一个有效的函数。这是因为要找到位于一个单元中的粒子需要花费
     * ${\cal O}(\log N)$ ，其中 $N$
     * 是整体粒子的数量。由于你可能会对每个单元都这样做，并假设粒子的数量和单元的数量大致成正比，你最终会得到一个
     * ${\cal O}(N \log N)$
     * 算法。一个更好的方法是利用这样一个事实：在内部，粒子是按照它们所在的活动单元的顺序排列的。换句话说，如果你遍历所有的粒子，当你走过活动单元时，你会以同样的顺序遇到它们。你可以利用这一点，在第一个单元的第一个粒子上保留一个迭代器，当你移动到下一个单元时，你也会递增粒子迭代器，直到你找到位于下一个单元的粒子。例如，这是
     * step-70
     * 中使用的方法，当累积到所有单元时，其总成本为
     * ${\cal O}(\log N)$ 。该方法在 step-19 的 "扩展可能性部分
     * "中也有详细说明。
     *
     */
    particle_iterator_range
    particles_in_cell(
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell);

    /**
     * 返回一对粒子迭代器，标记特定单元中粒子的开始和结束。最后一个迭代器是已不在单元格中的第一个粒子。
     * 返回范围内的元素数量等于n_particles_in_cell()函数返回的内容。
     * @note  虽然这个函数在 step-19
     * 中被使用，但如果粒子的数量很大，它并不是一个有效的函数。这是因为要找到位于一个单元中的粒子需要花费
     * ${\cal O}(\log N)$ ，其中 $N$
     * 是整体粒子的数量。由于你可能会对每个单元都这样做，并假设粒子的数量和单元的数量大致成正比，你最终会得到一个
     * ${\cal O}(N \log N)$
     * 的算法。一个更好的方法是利用这样一个事实：在内部，粒子是按照它们所在的活动单元的顺序排列的。换句话说，如果你遍历所有的粒子，当你在活动单元上行走时，你会以同样的顺序遇到它们。你可以利用这一点，在第一个单元的第一个粒子上保留一个迭代器，当你移动到下一个单元时，你也会递增粒子迭代器，直到你找到位于下一个单元的粒子。例如，这是在
     * step-70 中使用的方法，在所有单元中累积时，总成本为
     * ${\cal O}(\log N)$ 。该方法在 step-19 的 "扩展可能性部分
     * "中也有详细说明。
     *
     */
    particle_iterator_range
    particles_in_cell(
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell)
      const;

    /**
     * 移除迭代器所指向的一个粒子。
     *
     */
    void
    remove_particle(const particle_iterator &particle);

    /**
     * 在粒子集合中插入一个粒子。返回一个迭代器到粒子的新位置。这个函数涉及到一个粒子及其属性的副本。请注意，对于
     * $N$ 粒子来说，这个函数具有 $O(N \log N)$ 的复杂性。
     *
     */
    particle_iterator
    insert_particle(
      const Particle<dim, spacedim> &particle,
      const typename Triangulation<dim, spacedim>::active_cell_iterator &cell);

    /**
     * 在粒子集合中插入一定数量的粒子。
     * 这个函数涉及到粒子和它们的属性的一个副本。
     * 注意，这个函数的复杂度是O(n_existing_particles +
     * n_particles)。
     *
     */
    void
    insert_particles(
      const std::multimap<
        typename Triangulation<dim, spacedim>::active_cell_iterator,
        Particle<dim, spacedim>> &particles);

    /**
     * 创建并插入若干粒子到粒子集合中。
     * 这个函数接收一个位置列表，并在这些位置上创建一组粒子，然后将其添加到本地粒子集合中。请注意，这个函数目前使用的是
     * GridTools::compute_point_locations(),
     * ，它假定所有的位置都在三角形的本地部分内。如果其中一个不在本地域内，这个函数将抛出一个异常。
     *
     */
    void
    insert_particles(const std::vector<Point<spacedim>> &positions);

    /**
     * 创建并在粒子集合中插入若干粒子。
     * 这个函数接收一个位置列表，并在这些位置上创建一组粒子，然后将其分配并添加到处理器的本地粒子集合中。请注意，这个函数使用了
     * GridTools::distributed_compute_point_locations().
     * 因此，它可能需要处理器之间的密集通信。这个函数使用的是
     * step-70  。
     * 这个函数找出哪个mpi进程拥有不属于三角形的本地拥有部分的点，它把在这个进程上传递给这个函数的点发送给该进程，并从收到输入的人那里接收属于三角形的本地拥有单元的点。
     * 为了跟踪哪个 mpi
     * 进程收到了哪些点，函数会返回一个从 mpi 进程到
     * IndexSet 的映射。这个 IndexSet 包含在调用 mpi
     * 进程上传递给这个函数的点的本地索引，这些点属于这个
     * mpi 进程所拥有的三角形的一部分。        如果 @p ids
     * 的向量是空的，那么这些粒子的id是由get_next_free_particle_index()开始自动计算的。
     * 例如，如果get_next_free_particle_index()方法返回n0，在两个MPI进程分别添加n1和n2粒子的情况下调用这个函数，将导致进程0添加的n1粒子的id等于`[n0,n0+n1)`，而进程1添加的n2粒子的id为`[n0+n1,
     * n0+n1+n2)`。          @param[in]  positions
     * 一个不需要在本地处理器上的点的向量，但必须在与这个ParticleHandler对象相关的三角图中。
     * @param[in]  global_bounding_boxes 一个边界盒的向量。
     * 界限盒`global_bubes[rk]`描述网格的哪一部分是由等级`rk`的mpi进程局部拥有的。局部描述可以从
     * GridTools::compute_mesh_predicate_bounding_box(),
     * 得到，全局描述可以通过将局部的传递给
     * Utilities::MPI::all_gather().  @param[in]
     * 属性（可选）一个与每个局部点相关的属性向量。矢量的大小应该是零（没有属性会被转移，也不会附加到生成的粒子上），或者它应该是一个大小为`n_properties_per_particle()`的`positions.size()`矢量的矢量。请注意，这个函数调用将把属性从本地的mpi进程转移到将拥有每个粒子的最终mpi进程，因此它可能是通信密集型的。
     * @param[in]  ids （可选）与每个粒子相关联的id的向量。
     * 如果该向量为空，那么id将被分配为从第一个可用索引开始的连续范围，如上文所述。如果该向量不是空的，那么它的大小必须与
     * @p positions 向量的大小一致。          @return
     * 一个从所有者到IndexSet的映射，它包含在调用mpi进程上传递给这个函数的点的局部索引，并且属于这个mpi进程所拥有的三角测量部分。
     *
     */
    std::map<unsigned int, IndexSet>
    insert_global_particles(
      const std::vector<Point<spacedim>> &positions,
      const std::vector<std::vector<BoundingBox<spacedim>>>
        &                                       global_bounding_boxes,
      const std::vector<std::vector<double>> &  properties = {},
      const std::vector<types::particle_index> &ids        = {});

    /**
     * 在粒子集合中插入一定数量的粒子。这个函数接收一个我们不知道相关单元迭代器的粒子列表，并将其分配到一个处理器的正确的本地粒子集合中，方法是解压位置，通过调用
     * GridTools::distributed_compute_point_locations(),
     * 计算出粒子的发送位置，并将粒子发送到相应的进程。
     * 为了跟踪哪个mpi进程收到了什么粒子，函数会返回一个从mpi进程到IndexSet的映射。这个IndexSet包含了在调用mpi进程上传递给这个函数的粒子的本地索引，并且属于这个mpi进程所拥有的三角形的一部分。
     * @param[in]  particles
     * 一个不需要在本地处理器上的粒子的向量。
     * @param[in]  global_bounding_boxes 一个边界盒的向量。
     * 界限盒`global_bubes[rk]`描述网格的哪一部分在本地被等级为`rk`的mpi进程拥有。本地描述可以从
     * GridTools::compute_mesh_predicate_bounding_box(),
     * 获得，全局描述可以通过将本地描述传递给
     * Utilities::MPI::all_gather().   @return
     * 一个从所有者到IndexSet的映射，它包含在调用mpi进程上传递给这个函数的点的本地索引，并且属于这个mpi进程拥有的三角化部分。
     *
     */
    std::map<unsigned int, IndexSet>
    insert_global_particles(
      const std::vector<Particle<dim, spacedim>> &particles,
      const std::vector<std::vector<BoundingBox<spacedim>>>
        &global_bounding_boxes);

    /**
     * 通过使用向量 @p input_vector.  @tparam
     * 中包含的值来设置粒子的位置 VectorType
     * 该库支持的任何并行分布式向量。        矢量 @p
     * input_vector
     * 应该可以读取通过用local_owned_particle_ids()提取本地相关id创建的索引，并将其与代表范围`[0,
     * spacedim)`的索引集进行张量积，即。
     * @code
     * IndexSet ids = particle_handler.locally_owned_particle_ids().
     * tensor_product(complete_index_set(spacedim));
     * @endcode
     * 具有全局索引`id`的粒子的位置从spacedim的连续条目中读取，从`input_vector[id*spacedim]`开始。
     * 注意，不需要 @p input_vectorowns*
     * 这些索引，但是它必须有对这些索引的读取权限（也就是说，它可以是一个有鬼魂条目的分布式向量）。
     * 如果参数 @p displace_particles
     * 被设置为false，那么新的位置就会从 @p input_vector,
     * 中包含的数值中提取，以取代之前存储的粒子位置。
     * 默认情况下，粒子的位移量为 @p input_vector,
     * 中包含的量，也就是说，向量的内容被认为是偏移量*，被添加到之前的位置。
     * 在设置新的位置之后，这个函数在内部调用sort_particles_into_subdomains_and_cells()方法。你应该确保你满足该函数的要求。
     * @param[in]  input_vector
     * 一个平行分布的矢量，包含应用于每个粒子的位移，或者它们新的绝对位置。
     * @param[in]  displace_particles 控制 @p input_vector
     * 应该被解释为位移向量，还是绝对位置的向量。
     *
     */
    template <class VectorType>
    typename std::enable_if<
      std::is_convertible<VectorType *, Function<spacedim> *>::value ==
      false>::type
    set_particle_positions(const VectorType &input_vector,
                           const bool        displace_particles = true);

    /**
     * 在粒子处理程序中使用一个点的矢量来设置粒子的位置。由矢量定义的新的点集必须足够接近原始点集，以确保sort_particles_into_subdomains_and_cells()函数能够找到粒子所属的新单元。
     * 点的编号与ParticleHandler本地遍历的方式相同。使用这个方法的典型方法，是首先调用get_particle_positions()函数，然后修改产生的向量。
     * @param  [in] new_positions
     * 尺寸为particle_handler.n_locally_owned_particles()  @param  [in]
     * displace_particles
     * 当为真时，该函数将点的矢量值添加到粒子的当前位置，从而将它们置换成函数所给的数量。当为假时，粒子的位置会被矢量中的值所取代。
     *
     */
    void
    set_particle_positions(const std::vector<Point<spacedim>> &new_positions,
                           const bool displace_particles = true);


    /**
     * 使用带有spacedim组件的函数，在粒子处理程序中设置粒子的位置。由该函数定义的新的点集必须足够接近原始点集，以确保sort_particles_into_subdomains_and_cells算法能够找到粒子所属的新单元。
     * 该函数在粒子的当前位置被评估。          @param  [in]
     * function
     * 一个具有n_components==spacedim的函数，描述粒子的位移或新位置，作为粒子当前位置的函数。
     * @param  [in] displace_particles
     * 当为真时，该函数将函数的结果添加到粒子的当前位置，从而使其位移的量由函数给出。当为假时，粒子的位置会被函数的值所取代。
     *
     */
    void
    set_particle_positions(const Function<spacedim> &function,
                           const bool                displace_particles = true);

    /**
     * 读取粒子的位置，并将其存储到分布式矢量 @p
     * output_vector. 中，默认情况下， @p output_vector
     * 被此操作覆盖，但你可以通过设置 @p add_to_output_vector
     * 为`true`来增加其条目。          @tparam  VectorType
     * 该库支持的任何一个并行分布式向量。
     * 这是set_particle_positions()函数的反向操作。
     * 全局索引为`id`的粒子的位置被写入从`output_vector[id*spacedim]`开始的spacedim连续条目中。
     * 注意，如果你使用分布式矢量类型， @p output_vector
     * 没有必要拥有将被写入的索引所对应的条目。然而你应该记住，这需要一个全局通信，将上面的条目分配给它们各自的所有者。
     * @param[in,  out] output_vector
     * 一个平行分布的向量，包含粒子的位置，或者用粒子的位置更新。
     * @param[in]  add_to_output_vector 控制函数是否应该设置 @p
     * output_vector 的条目，或者是否应该向其添加。
     *
     */
    template <class VectorType>
    void
    get_particle_positions(VectorType &output_vector,
                           const bool  add_to_output_vector = false);

    /**
     * 将粒子处理程序中的粒子位置收集到一个点的矢量中。这些点的顺序与在所有（本地）粒子上进行迭代，并查询它们的位置所得到的相同。
     * @param  [in,out] positions
     * 以`particle_handler.n_locally_owned_articles`大小预先分配的向量，其点将成为本地拥有的粒子的位置
     * @param  [in] add_to_output_vector
     * 当真时，粒子的点的值被添加到位置向量。当为假时，位置向量中的点的值会被粒子的位置所取代。
     *
     */
    void
    get_particle_positions(std::vector<Point<spacedim>> &positions,
                           const bool add_to_output_vector = false);

    /**
     * 这个函数允许注册三个额外的函数，这些函数在每次粒子被转移到另一个进程时（即在分拣到单元格时，在幽灵粒子转移时，或在所有粒子的序列化时）被调用。
     * @param  size_callback
     * 一个在序列化粒子数据时被调用的函数。该函数没有得到任何参数，并且预计将返回每个粒子被序列化的额外数据的大小。请注意，目前这意味着每个粒子的数据大小都必须是相同的。
     * @param  store_callback
     * 在序列化粒子数据时，每个粒子会被调用一次的函数。该函数的参数是一个识别当前粒子的粒子迭代器和一个指向大小为_callback()的数据块的无效指针，该函数可以在其中存储额外的数据。该函数预计将返回一个指向其数据块之后的位置的无效指针。
     * @param  load_callback
     * 在反序列化粒子数据时，每个粒子都会被调用一次的函数。该函数的参数是一个识别当前粒子的粒子迭代器和一个指向大小为_callback()的数据块的无效指针，其中的额外数据被store_callback函数存储。该函数预计将返回一个指向其数据块之后的位置的无效指针。
     *
     */
    void
    register_additional_store_load_functions(
      const std::function<std::size_t()> &size_callback,
      const std::function<void *(const particle_iterator &, void *)>
        &store_callback,
      const std::function<const void *(const particle_iterator &, const void *)>
        &load_callback);

    /**
     * 返回上次调用update_cached_numbers()函数时，由这个类管理的粒子总数。
     * 如果粒子被添加或移除，实际的粒子数可能会在那之后发生变化。
     * @return 模拟中的粒子总数。
     *
     */
    types::particle_index
    n_global_particles() const;

    /**
     * 返回上次调用update_cached_numbers()函数时每个单元的最大粒子数。
     * @return  仿真中一个单元中的最大粒子数。
     *
     */
    types::particle_index
    n_global_max_particles_per_cell() const;

    /**
     * 返回三角形局部的粒子数。
     *
     */
    types::particle_index
    n_locally_owned_particles() const;

    /**
     * 返回上次调用update_cached_numbers()函数时，全局粒子集合中的下一个空闲粒子索引。
     *
     */
    types::particle_index
    get_next_free_particle_index() const;

    /**
     * 提取一个全局尺寸等于get_next_free_particle_index()的IndexSet，包含本地拥有的粒子指数。
     * 这个函数可以用来构造分布式向量和矩阵，以使用线性代数操作来操作粒子。
     * 请注意，用户有责任保证粒子索引是唯一的，而且没有进行检查以验证这一点，也没有检查每个mpi进程上所有IndexSet对象的联盟是否完整。
     * @return
     * 一个大小为get_next_free_particle_index()的IndexSet，包含n_locally_owned_particle()指数。
     * @deprecated  使用local_owned_particle_ids()代替。
     *
     */
    DEAL_II_DEPRECATED IndexSet
                       locally_relevant_ids() const;

    /**
     * 提取一个全局尺寸等于get_next_free_particle_index()的IndexSet，包含本地拥有的粒子指数。
     * 这个函数可以用来构造分布式向量和矩阵，以使用线性代数操作来操作粒子。
     * 请注意，用户有责任保证粒子索引是唯一的，而且没有进行检查以验证这一点，也没有检查每个mpi进程上所有IndexSet对象的联盟是否完整。
     * @return  一个大小为 get_next_free_particle_index() 的
     * IndexSet，包含 n_locally_owned_particle() 指数。
     *
     */
    IndexSet
    locally_owned_particle_ids() const;

    /**
     * 返回每个粒子拥有的属性数量。
     *
     */
    unsigned int
    n_properties_per_particle() const;

    /**
     * 返回一个对拥有所有粒子属性的属性池的引用，并对它们进行物理组织。
     *
     */
    PropertyPool<dim, spacedim> &
    get_property_pool() const;

    /**
     * 为所有本地拥有的粒子找到并更新包含每个粒子的单元。如果粒子移出了本地子域，它们将被发送到它们的新进程并插入那里。
     * 在这个函数调用之后，每个粒子要么在其当前进程和当前单元中，要么被删除（如果它不能找到其新进程或单元）。
     * 用户可以在信号上附加一个函数
     * Particles::ParticleHandler::Signals::particle_lost().
     * 每当一个粒子被删除时，信号就会被触发，连接的函数会被调用，并传递给有关粒子的迭代器，以及它最后已知的细胞关联。
     *
     */
    void
    sort_particles_into_subdomains_and_cells();

    /**
     * 将所有生活在细胞中的粒子换成其他进程的幽灵细胞。清除并重新填充ghost_neighbors成员变量。
     *
     */
    void
    exchange_ghost_particles(const bool enable_ghost_cache = false);

    /**
     * 将所有生活在作为幽灵细胞的细胞中的粒子更新给其他进程。在这里，更新意味着更新幽灵粒子的位置和属性，假设幽灵粒子没有改变细胞。因此，这将不会更新粒子的参考位置。
     *
     */
    void
    update_ghost_particles();

    /**
     * 回调函数，应该在每次细化前和写检查点时被调用。这个函数用于将store_particles()注册到三角结构中。这个函数在
     * step-70  中使用。
     *
     */
    void
    register_store_callback_function();

    /**
     * 回调函数，应该在每次细化后和从检查点恢复后被调用。
     * 这个函数用于将load_particles()注册到三角结构中。这个函数在
     * step-70  中使用。
     *
     */
    void
    register_load_callback_function(const bool serialization);

    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将该类的内容序列化。
     *
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int version);

    /**
     * 一个拥有 boost::signal
     * 对象的结构，用于粒子处理程序可以对自己进行的一些操作。如何在应用中使用信号，在三角形类中的
     * "当三角形发生变化时获得通知
     * "一节有更多的信息和例子解释。简而言之，这些信号允许粒子处理程序通知应用程序关于粒子处理程序内部的某些事件，例如，当一个粒子丢失时。
     * 关于信号的文档，请参见http://www.boost.org/doc/libs/release/libs/signals2
     * 。
     *
     */
    struct Signals
    {
      /**
       * 每当 ParticleHandler::sort_particles_into_subdomains_and_cells()
       * 函数遇到一个不能与细胞关联的粒子时，就会触发这个信号。这可能发生在粒子离开三角形的领域，或者在平行三角形中离开本地已知的领域（包括
       * parallel::distributed::triangulation). 的幽灵单元）。
       * 连接的函数接收到有关粒子的迭代器，以及其最后已知的单元关联。
       * 这个信号用于  step-19  。
       *
       */
      boost::signals2::signal<void(
        const typename Particles::ParticleIterator<dim, spacedim> &particle,
        const typename Triangulation<dim, spacedim>::active_cell_iterator
          &cell)>
        particle_lost;
    };

    /**
     * 粒子处理程序可以通知调用应用程序的事件的信号。
     *
     */
    mutable Signals signals;

  private:
    /**
     * 要工作的三角区的地址。
     *
     */
    SmartPointer<const Triangulation<dim, spacedim>,
                 ParticleHandler<dim, spacedim>>
      triangulation;

    /**
     * 要处理的映射的地址。
     *
     */
    SmartPointer<const Mapping<dim, spacedim>, ParticleHandler<dim, spacedim>>
      mapping;

    /**
     * 这个对象拥有并组织了所有粒子属性的内存。由于粒子引用了属性池，后者必须在*粒子被销毁之后被销毁。
     * 这可以通过确保`property_pool`成员变量在`particles`和`ghost_particles`成员的声明之前实现。
     *
     */
    std::unique_ptr<PropertyPool<dim, spacedim>> property_pool;

    /**
     * 当前生活在本地域中的粒子集合，按照它们所在的单元的级别/索引组织。
     *
     */
    std::multimap<internal::LevelInd, Particle<dim, spacedim>> particles;

    /**
     * 目前生活在本地域的幽灵单元中的粒子集合，按它们所在的单元的级别/索引组织。这些粒子相当于分布式向量中的幽灵条目。
     *
     */
    std::multimap<internal::LevelInd, Particle<dim, spacedim>> ghost_particles;

    /**
     * 这个变量存储了全局存储的粒子数量。它是由update_cached_numbers()计算的。
     *
     */
    types::particle_index global_number_of_particles;

    /**
     * 全局域中每个单元的最大粒子数。这个变量对于在解决方案的重新分区和序列化过程中存储和加载粒子数据非常重要。请注意，该变量只有在需要时才会被更新，例如，在粒子移动之后、网格细化之前/之后、创建检查点之前以及从检查点恢复之后。
     *
     */
    unsigned int global_max_particles_per_cell;

    /**
     * 这个变量存储了全局可用的下一个空闲粒子索引，以备生成新粒子。
     *
     */
    types::particle_index next_free_particle_index;

    /**
     * 一个可以通过调用register_additional_store_load_functions注册的函数。它在序列化粒子数据时被调用。该函数没有得到任何参数，并被期望返回每个粒子被序列化的额外数据的大小。请注意，目前这意味着每个粒子的数据大小必须是相同的，但并不是每个序列化过程都必须是相同的（例如，粒子运动过程中的序列化可能包括临时数据，而运动结束后的序列化则不需要传输这些数据）。
     *
     */
    std::function<std::size_t()> size_callback;

    /**
     * 一个可以通过调用register_additional_store_load_functions注册的函数。在序列化粒子数据时，每个粒子会被调用一次。该函数的参数是一个识别当前粒子的粒子迭代器和一个指向大小为_callback()的数据块的无效指针，该函数可以在其中存储额外的数据。该函数预计将返回一个指向其数据块之后的位置的无效指针。
     *
     */
    std::function<void *(const particle_iterator &, void *)> store_callback;

    /**
     * 在对粒子数据进行反序列化时，每个粒子都会被调用一次的函数。该函数的参数是一个识别当前粒子的粒子迭代器和一个指向大小为_callback()的数据块的无效指针，该函数可以从中加载额外的数据。这个数据块在序列化过程中由store_callback函数填充。这个函数被期望返回一个指向其数据块之后的位置的无效指针。
     *
     */
    std::function<const void *(const particle_iterator &, const void *)>
      load_callback;

    /**
     * 这个变量由register_store_callback_function()函数设置，并由register_load_callback_function()函数使用，以检查粒子数据在相应的三角测量对象中注册的位置。
     *
     */
    unsigned int handle;

    /**
     * GridTools::Cache
     * 用于存储顶点到单元格集合和顶点到单元格中心向量的信息，以防止每次我们sort_into_subdomain_and_cells()时重新计算它们。
     * 当三角结构发生变化时，这个缓存会自动更新。这个缓存被存储在一个唯一的指针中，因为粒子处理程序有一个构造函数，可以在没有三角结构的情况下构造它。缓存没有这样的构造函数。
     *
     */
    std::unique_ptr<GridTools::Cache<dim, spacedim>> triangulation_cache;

#ifdef DEAL_II_WITH_MPI
    /**
     * 将已经越过子域边界的粒子转移到其他处理器。
     * 所有收到的粒子和它们的新单元将被附加到 @p
     * received_particles 向量中。          @param  [in] particles_to_send
     * 所有应该被发送的粒子和它们的新子域_id都在这个图中。
     * @param  [in,out] received_particles
     * 存储所有收到的粒子的向量。请注意，不要求也不检查该列表是否为空，收到的粒子只是被附加到矢量的末端。
     * @param  [in] new_cells_for_particles
     * 可选的单元格迭代器向量，其结构与 @p particles_to_send.
     * 相同，如果给出这个参数，它应该包含每个要发送的粒子的单元格迭代器，该粒子属于其中。如果粒子迭代器的单元信息已经过期（例如，在粒子移动之后），这个参数是必要的。
     * @param  [in] enable_cache
     * 可选的bool，通过建立一个GhostParticlePartitioner类型的缓存，存储更新幽灵粒子的必要信息，使幽灵粒子不需要从头开始重建。
     * 一旦建立了这个缓存，就可以通过调用send_recv_particles_properties_and_location（）来更新幽灵粒子。
     *
     */
    void
    send_recv_particles(
      const std::map<types::subdomain_id, std::vector<particle_iterator>>
        &particles_to_send,
      std::multimap<internal::LevelInd, Particle<dim, spacedim>>
        &received_particles,
      const std::map<
        types::subdomain_id,
        std::vector<
          typename Triangulation<dim, spacedim>::active_cell_iterator>>
        &new_cells_for_particles = std::map<
          types::subdomain_id,
          std::vector<
            typename Triangulation<dim, spacedim>::active_cell_iterator>>(),
      const bool enable_cache = false);

    /**
     * 假设粒子没有改变单元格，则传输粒子的位置和属性。这个例程使用GhostParticlePartitioner作为一个缓存结构来更新粒子。
     * 它内在地假定粒子不能改变单元。
     * 所有更新的粒子将被追加到 @p received_particles 容器中。
     * @param  [in] particles_to_send
     * 所有应该发送信息的粒子和它们的新子域_id都在这个地图中。
     * @param  [in,out] received_particles
     * 一个包含所有收到的粒子的地图。注意，不要求也不检查容器是否为空，收到的粒子只是被插入到地图中。
     *
     */
    void
    send_recv_particles_properties_and_location(
      const std::map<types::subdomain_id, std::vector<particle_iterator>>
        &particles_to_send,
      std::multimap<internal::LevelInd, Particle<dim, spacedim>>
        &received_particles);


#endif

    /**
     * 缓存结构，用于存储跨处理器交换粒子信息（位置和属性）所需的元素，以便更新幽灵粒子。这个结构只用于更新鬼魂粒子。
     *
     */
    internal::GhostParticlePartitioner<dim, spacedim> ghost_particles_cache;

    /**
     * 在细化步骤之前，由Triangulation的监听函数对每个单元进行调用。所有的粒子都必须与它们的单元相连，才能被送到新的进程中去。
     *
     */
    std::vector<char>
    store_particles(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      const typename Triangulation<dim, spacedim>::CellStatus     status) const;

    /**
     * 在细化步骤之后，由听众函数调用。粒子的局部地图必须从三角形的user_pointer中读取。
     *
     */
    void
    load_particles(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell,
      const typename Triangulation<dim, spacedim>::CellStatus     status,
      const boost::iterator_range<std::vector<char>::const_iterator>
        &data_range);
  };



   /* ---------------------- inline and template functions ------------------
   */ 

  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::begin() const
  {
    return (const_cast<ParticleHandler<dim, spacedim> *>(this))->begin();
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::begin()
  {
    return particle_iterator(particles, particles.begin());
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::end() const
  {
    return (const_cast<ParticleHandler<dim, spacedim> *>(this))->end();
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::end()
  {
    return particle_iterator(particles, particles.end());
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::begin_ghost() const
  {
    return (const_cast<ParticleHandler<dim, spacedim> *>(this))->begin_ghost();
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::begin_ghost()
  {
    return particle_iterator(ghost_particles, ghost_particles.begin());
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::end_ghost() const
  {
    return (const_cast<ParticleHandler<dim, spacedim> *>(this))->end_ghost();
  }



  template <int dim, int spacedim>
  inline typename ParticleHandler<dim, spacedim>::particle_iterator
  ParticleHandler<dim, spacedim>::end_ghost()
  {
    return particle_iterator(ghost_particles, ghost_particles.end());
  }



  template <int dim, int spacedim>
  template <class Archive>
  inline void
  ParticleHandler<dim, spacedim>::serialize(Archive &ar, const unsigned int)
  {
    // Note that we do not serialize the particle data itself. Instead we
    // use the serialization functionality of the triangulation class, because
    // this guarantees that data is immediately shipped to new processes if
    // the domain is distributed differently after resuming from a checkpoint.
    ar //&particles
      &global_number_of_particles &global_max_particles_per_cell
        &                          next_free_particle_index;
  }



  template <int dim, int spacedim>
  template <class VectorType>
  inline typename std::enable_if<
    std::is_convertible<VectorType *, Function<spacedim> *>::value ==
    false>::type
  ParticleHandler<dim, spacedim>::set_particle_positions(
    const VectorType &input_vector,
    const bool        displace_particles)
  {
    AssertDimension(input_vector.size(),
                    get_next_free_particle_index() * spacedim);
    for (auto &p : *this)
      {
        auto       new_point(displace_particles ? p.get_location() :
                                            Point<spacedim>());
        const auto id = p.get_id();
        for (unsigned int i = 0; i < spacedim; ++i)
          new_point[i] += input_vector[id * spacedim + i];
        p.set_location(new_point);
      }
    sort_particles_into_subdomains_and_cells();
  }



  template <int dim, int spacedim>
  template <class VectorType>
  inline void
  ParticleHandler<dim, spacedim>::get_particle_positions(
    VectorType &output_vector,
    const bool  add_to_output_vector)
  {
    AssertDimension(output_vector.size(),
                    get_next_free_particle_index() * spacedim);
    for (const auto &p : *this)
      {
        auto       point = p.get_location();
        const auto id    = p.get_id();
        if (add_to_output_vector)
          for (unsigned int i = 0; i < spacedim; ++i)
            output_vector[id * spacedim + i] += point[i];
        else
          for (unsigned int i = 0; i < spacedim; ++i)
            output_vector[id * spacedim + i] = point[i];
      }
    if (add_to_output_vector)
      output_vector.compress(VectorOperation::add);
    else
      output_vector.compress(VectorOperation::insert);
  }



  template <int dim, int spacedim>
  inline IndexSet
  ParticleHandler<dim, spacedim>::locally_relevant_ids() const
  {
    return this->locally_owned_particle_ids();
  }

} // namespace Particles

DEAL_II_NAMESPACE_CLOSE

#endif


