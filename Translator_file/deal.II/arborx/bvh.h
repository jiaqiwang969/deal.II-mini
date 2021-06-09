//include/deal.II-translator/arborx/bvh_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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

#ifndef dealii_arborx_bvh_h
#define dealii_arborx_bvh_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_ARBORX
#  include <deal.II/arborx/access_traits.h>

#  include <ArborX_LinearBVH.hpp>
#  include <Kokkos_Core.hpp>

DEAL_II_NAMESPACE_OPEN

/**
 * 这个命名空间包含ArborX库的包装器。
 *
 *
 */
namespace ArborXWrappers
{
  /**
   * 该类实现了对 ArborX::BVH. BVH的包装器，BVH代表Bounding
   * Volume Hierarchy。
   * 来自[维基百科](https://en.wikipedia.org/wiki/Bounding_Volume_Hierarchy)。
   * <blockquote>
   * 边界体积层次结构（BVH）是一组几何对象上的树状结构。所有的几何对象都被包裹在构成树的叶子节点的边界体中。然后，这些节点被分组为小集，并被包围在更大的边界体积中。这些节点反过来也被分组并以递归的方式包围在其他更大的边界体中，最终形成一个树状结构，在树的顶端有一个单一的边界体。边界体积层次结构被用来有效地支持对几何对象集的若干操作，例如在碰撞检测和光线追踪中。
   * </blockquote>
   * 因为ArborX使用了Kokkos，所以在使用这个类之前需要初始化和最终确定Kokkos。
   *
   */
  class BVH
  {
  public:
    /**
     * 构造函数。使用BoundingBox @p bounding_boxes
     * 的一个向量作为基元。
     *
     */
    template <int dim, typename Number>
    BVH(const std::vector<BoundingBox<dim, Number>> &bounding_boxes);

    /**
     * 构造函数。使用一个 @p points 的向量作为基元。
     *
     */
    template <int dim, typename Number>
    BVH(const std::vector<Point<dim, Number>> &points);

    /**
     * 返回那些满足 @p queries. 的BoundingBox对象的索引 因为 @p
     * queries
     * 可以包含多个查询，所以函数返回一对索引和偏移量。
     * 下面是一段做了以下工作的示例代码。让我们假设我们有一组对象的界线盒
     *
     * --比如说，一个三角形中每个单元的边界框，或者一个平行三角形中每个部分的边界框。我们将把这些存储在下面的`bvh_bounding_boxes`数组中。
     * 然后让我们假设我们有一组其他的边界框，比方说在我们领域中移动的小物体。我们将把这些边界框放到`bb_intersect`数组中。然后我们要回答的问题是：每个bb_intersect边界盒与BVH的哪个边界盒相交？换句话说，粒子在哪个单元或分区？
     * 这个问题可以通过下面的代码来回答。
     * @code
     * const std::vector<BoundingBox<dim>> query_bounding_boxes = ...
     * ArborXWrappers::BoundingBoxIntersectPredicate
     * bb_intersect(query_bounding_boxes);
     *
     * const std::vector<BoundingBox<dim>> bvh_bounding_boxes = ...
     * ArborxWrappers::BVH bvh(bvh_bounding_boxes);
     *
     * auto [indices, offset] = bvh.query(bb_intersect);
     * @endcode
     * `bvh_bounding_boxes`中与`query_bounding_boxes`的第`j`个BoundingBox相交的元素由以下方法给出。
     * @code
     * std::vector<int> bvh_bounding_box_indices;
     * for (int i = offset[j]; i < offset[j+1]; ++i)
     * bvh_bounding_box_indices.push_back(indices[i]);
     * @endcode
     * 在许多其他应用中，我们不仅对寻找另一个边界盒所在的边界盒感兴趣，而且对个别点所在的边界盒感兴趣
     *
     * 比方说，如果我们没有物体，而是有点状的粒子在移动。在这种情况下，我们需要回答一个关于点的查询，这可以按以下方式进行。
     * @code
     * const std::vector<Point<dim>> query_points = ...
     * ArborXWrappers::PointIntersectPredicate pt_intersect(query_points);
     *
     * const std::vector<BoundingBox<dim>> bvh_bounding_boxes = ...
     * ArborxWrappers::BVH bvh(bvh_bounding_boxes);
     *
     * auto [indices, offset] = bvh.query(pt_intersect);
     * @endcode
     * 作为最后一个例子，我们要展示如何找到一组给定的点中最近的五个点。这可以按以下方式进行。
     * @code
     * const std::vector<Point<dim>> query_points = ...
     * ArborXWrappers::PointNearestPredicate pt_nearest(query_points, 5);
     *
     * const std::vector<Point<dim>> bvh_points = ...
     * ArborxWrappers::BVH bvh(bvh_points);
     *
     * auto [indices, offset] = bvh.query(pt_nearest);
     * @endcode
     *
     */
    template <typename QueryType>
    std::pair<std::vector<int>, std::vector<int>>
    query(const QueryType &queries);

  private:
    /**
     * 底层的ArborX对象。
     *
     */
    ArborX::BVH<Kokkos::HostSpace> bvh;
  };



  template <int dim, typename Number>
  BVH::BVH(const std::vector<BoundingBox<dim, Number>> &bounding_boxes)
    : bvh(Kokkos::DefaultHostExecutionSpace{}, bounding_boxes)
  {}



  template <int dim, typename Number>
  BVH::BVH(const std::vector<Point<dim, Number>> &points)
    : bvh(Kokkos::DefaultHostExecutionSpace{}, points)
  {}



  template <typename QueryType>
  std::pair<std::vector<int>, std::vector<int>>
  BVH::query(const QueryType &queries)
  {
    Kokkos::View<int *, Kokkos::HostSpace> indices("indices", 0);

    Kokkos::View<int *, Kokkos::HostSpace> offset("offset", 0);
    ArborX::query(
      bvh, Kokkos::DefaultHostExecutionSpace{}, queries, indices, offset);
    std::vector<int> indices_vector;
    indices_vector.insert(indices_vector.begin(),
                          indices.data(),
                          indices.data() + indices.extent(0));
    std::vector<int> offset_vector;
    offset_vector.insert(offset_vector.begin(),
                         offset.data(),
                         offset.data() + offset.extent(0));

    return {indices_vector, offset_vector};
  }
} // namespace ArborXWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
#endif


