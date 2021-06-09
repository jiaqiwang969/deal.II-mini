//include/deal.II-translator/grid/grid_tools_cache_0.txt
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

#ifndef dealii_grid_grid_tools_cache_h
#define dealii_grid_grid_tools_cache_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_tools_cache_update_flags.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/rtree.h>

#include <boost/signals2.hpp>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

namespace GridTools
{
  /**
   * 一个用于缓存三角测量的计算密集型信息的类。
   * 这个类将一个信号附加到构造时传递的Triangulation上，以跟踪细化的变化，并允许用户查询一些使用GridTools命名空间中的函数构造的数据结构，这些数据结构只计算一次，然后在这个类中缓存，以便在Triangulation没有变化时快速访问。
   * 请注意，这个类只注意到底层三角剖分因
   * Triangulation::Signals::any_change() 信号被触发而发生变化。
   * 如果三角剖分因其他原因而改变，例如，因为你将其与MappingQEulerian对象一起使用，该对象通过其自身的变换看到顶点，或者因为你手动改变了一些顶点位置，那么这个类中的一些结构就变得过时了，你将不得不通过手动调用mark_for_update()方法将它们标记为过时。
   *
   */
  template <int dim, int spacedim = dim>
  class Cache : public Subscriptor
  {
  public:
    /**
     * 构造函数。
     * 如果你提供了可选的`mapping`参数，那么只要需要映射，就会使用这个参数。
     * @param  tria 用于存储信息的三角形  @param  mapping
     * 计算缓存对象时使用的映射。
     *
     */
    Cache(const Triangulation<dim, spacedim> &tria,
          const Mapping<dim, spacedim> &      mapping =
            (ReferenceCells::get_hypercube<dim>()
               .template get_default_linear_mapping<dim, spacedim>()));

    /**
     * 销毁器。
     *
     */
    ~Cache() override;

    /**
     * 确保在随后调用这个类中定义的`get_*'函数时，标记为更新的对象被重新计算。
     * 注意，当你调用这个函数时，没有任何工作被执行。实际的数据结构会在你下次调用相应的`get_*`方法时进行计算。
     * @param 标志 更新的标志是什么？
     *
     */
    void
    mark_for_update(const CacheUpdateFlags &flags = update_all);


    /**
     * 返回由 GridTools::vertex_to_cell_map().
     * 计算的缓存的顶点到单元格图。
     *
     */
    const std::vector<
      std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>> &
    get_vertex_to_cell_map() const;

    /**
     * 返回由 GridTools::vertex_to_cell_centers_directions().
     * 计算的缓存的顶点到单元格中心的方向。
     *
     */
    const std::vector<std::vector<Tensor<1, spacedim>>> &
    get_vertex_to_cell_centers_directions() const;

    /**
     * 返回由 GridTools::extract_used_vertices().
     * 计算的已使用顶点的缓存图。
     *
     */
    const std::map<unsigned int, Point<spacedim>> &
    get_used_vertices() const;

    /**
     * 为顶点返回缓存的RTree对象，使用三角形的已用顶点构建。
     *
     */
    const RTree<std::pair<Point<spacedim>, unsigned int>> &
    get_used_vertices_rtree() const;

    /**
     * 返回缓存的单元格边界框的RTree对象，使用存储的三角形的活动单元格迭代器构建。对于
     * parallel::distributed::Triangulation
     * 对象，该函数也将返回鬼魂和人工单元的边界框。
     *
     */
    const RTree<
      std::pair<BoundingBox<spacedim>,
                typename Triangulation<dim, spacedim>::active_cell_iterator>> &
    get_cell_bounding_boxes_rtree() const;

    /**
     * 返回包含本地拥有的活动单元的边界盒的缓存RTree对象，该对象使用存储的三角形的活动单元迭代器构建。
     * 与前一个函数不同的是，这个函数只使用本地拥有的单元格来构建RTree，即不包括幽灵或人工单元。这两个函数在串行计算中返回相同的结果。
     *
     */
    const RTree<
      std::pair<BoundingBox<spacedim>,
                typename Triangulation<dim, spacedim>::active_cell_iterator>> &
    get_locally_owned_cell_bounding_boxes_rtree() const;


    /**
     * 返回包含每个顶点所连接的子域ID的整数集的向量。这个特征在粒子处理程序中被广泛使用，以检测必须在哪些处理器上建立幽灵粒子。
     *
     */
    const std::vector<std::set<unsigned int>> &
    get_vertex_to_neighbor_subdomain() const;

    /**
     * 返回一个对存储的三角图的引用。
     *
     */
    const Triangulation<dim, spacedim> &
    get_triangulation() const;

    /**
     * 返回一个对存储的映射的引用。
     *
     */
    const Mapping<dim, spacedim> &
    get_mapping() const;


    /**
     * 这个函数返回一个对象，该对象允许识别并行计算中的哪些进程可能拥有围绕给定点的单元。这个对象的元素是
     *
     * - 一个Rtree
     *
     * 是成对的包围盒，表示覆盖平行三角形的全部或部分本地部分的区域，以及一个无符号的int，代表拥有在给定包围盒内的单元的进程或子域。
     * 给定一个 parallel::TriangulationBase,
     * 上的点，这棵树允许识别一个或几个候选进程，对于这些进程，该点位于本地拥有的单元上。
     * 构建或更新rtree需要调用
     * GridTools::build_global_description_tree(), ，它使用
     * Utilities::MPI::all_gather(),
     * 的集体操作在所有进程之间交换边界盒。
     * 因此，这个函数必须由所有进程同时调用。
     * 本地边界盒是通过从get_locally_owned_cell_bounding_boxes_rtree()返回的rtree对象中提取指定的
     * @p level 来构造的。注意， @p level
     * 在这里不是指三角形的级别，而是指RTree对象的级别（例如，见https://en.wikipedia.org/wiki/R-tree）。
     *
     */
    const RTree<std::pair<BoundingBox<spacedim>, unsigned int>> &
    get_covering_rtree(const unsigned int level = 0) const;

  private:
    /**
     * 追踪下一步需要更新的内容。
     *
     */
    mutable CacheUpdateFlags update_flags;

    /**
     * 一个指向三角结构的指针。
     *
     */
    SmartPointer<const Triangulation<dim, spacedim>, Cache<dim, spacedim>> tria;

    /**
     * 在三角图上计算时使用的映射。
     *
     */
    SmartPointer<const Mapping<dim, spacedim>, Cache<dim, spacedim>> mapping;


    /**
     * 存储顶点到单元格的映射信息，由
     * GridTools::vertex_to_cell_map() 生成。
     *
     */
    mutable std::vector<
      std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>>
      vertex_to_cells;

    /**
     * 存储顶点到单元格的中心方向，由
     * GridTools::vertex_to_cell_centers_directions(). 生成。
     *
     */
    mutable std::vector<std::vector<Tensor<1, spacedim>>>
      vertex_to_cell_centers;

    /**
     * 一个覆盖整个网格的rtree对象的集合。
     * 地图的每个条目是由函数extract_rtree_level()应用于函数get_locally_owned_cell_bounding_boxes_rtree()返回的rtree构建的，输入的是指定级别。
     *
     */
    mutable std::map<unsigned int,
                     RTree<std::pair<BoundingBox<spacedim>, unsigned int>>>
      covering_rtree;

    /**
     * 存储由 GridTools::extract_used_vertices().
     * 生成的三角形的使用顶点。
     *
     */
    mutable std::map<unsigned int, Point<spacedim>> used_vertices;

    /**
     * 存储一个RTree对象，其中包含三角形的使用顶点。
     *
     */
    mutable RTree<std::pair<Point<spacedim>, unsigned int>> used_vertices_rtree;

    /**
     * 存储一个RTree对象，包含三角形的单元格的边界框。
     *
     */
    mutable RTree<
      std::pair<BoundingBox<spacedim>,
                typename Triangulation<dim, spacedim>::active_cell_iterator>>
      cell_bounding_boxes_rtree;

    /**
     * 存储一个RTree对象，包含三角形中本地拥有的单元格的边界框。
     *
     */
    mutable RTree<
      std::pair<BoundingBox<spacedim>,
                typename Triangulation<dim, spacedim>::active_cell_iterator>>
      locally_owned_cell_bounding_boxes_rtree;


    /**
     * 存储一个 std::vector 的 std::set
     * 的整数，包含一个顶点连接到的所有子域的id。
     *
     */
    mutable std::vector<std::set<unsigned int>> vertex_to_neighbor_subdomain;

    /**
     * 存储三角化信号的状态。
     *
     */
    boost::signals2::connection tria_signal;
  };



  // Inline functions
  template <int dim, int spacedim>
  inline const Triangulation<dim, spacedim> &
  Cache<dim, spacedim>::get_triangulation() const
  {
    return *tria;
  }



  template <int dim, int spacedim>
  inline const Mapping<dim, spacedim> &
  Cache<dim, spacedim>::get_mapping() const
  {
    return *mapping;
  }
} // namespace GridTools



DEAL_II_NAMESPACE_CLOSE

#endif


