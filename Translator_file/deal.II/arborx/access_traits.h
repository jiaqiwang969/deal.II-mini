//include/deal.II-translator/arborx/access_traits_0.txt
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

#ifndef dealii_arborx_access_traits_h
#define dealii_arborx_access_traits_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_ARBORX
#  include <deal.II/base/bounding_box.h>

#  include <ArborX.hpp>


DEAL_II_NAMESPACE_OPEN

namespace ArborXWrappers
{
  /**
   * 基于点的谓词的基类，为派生类提供基本功能，不应该单独使用。
   *
   */
  class PointPredicate
  {
  protected:
    /**
     * 构造函数。  @p points 是谓词所使用的点的列表。
     *
     */
    template <int dim, typename Number>
    PointPredicate(const std::vector<dealii::Point<dim, Number>> &points);

    /**
     * 存储在结构中的点的数量。
     *
     */
    std::size_t
    size() const;

    /**
     * 返回存储在对象中的第`i`个点。
     *
     */
    const dealii::Point<3, float> &
    get(unsigned int i) const;

  private:
    std::vector<dealii::Point<3, float>> points;
  };



  /**
   * 这个类定义了一个谓词，被 ArborXWrappers::BVH
   * 用来确定对于给定的点，用于建立 ArborXWrappers::BVH
   * 的边界盒与它们相交。
   * @note 该类不应该被用于多态的环境中。
   *
   */
  class PointIntersectPredicate : private PointPredicate
  {
  public:
    /**
     * 构造函数。  @p points
     * 是一个点的列表，我们有兴趣知道它们是否与
     * ArborXWrappers::BVH 边界盒相交。
     *
     */
    template <int dim, typename Number>
    PointIntersectPredicate(
      const std::vector<dealii::Point<dim, Number>> &points);

    // We need these since we inherit privately to avoid polymorphic use.
    using PointPredicate::get;
    using PointPredicate::size;
  };



  /**
   * 这个类定义了一个由 ArborXWrappers::BVH
   * 使用的谓词，以确定对于给定的点，哪些是用于建立
   * ArborXWrappers::BVH. 的最近的界线盒/点。
   * @note  该类不应该在多态环境中使用。
   *
   */
  class PointNearestPredicate : private PointPredicate
  {
  public:
    /**
     * 构造器。  @p points 是我们对 @p n_nearest_neighbors 中的
     * ArborXWrappers::BVH 边界盒/点感兴趣的点的列表。
     *
     */
    template <int dim, typename Number>
    PointNearestPredicate(const std::vector<dealii::Point<dim, Number>> &points,
                          const unsigned int n_nearest_neighbors);

    /**
     * 返回我们正在寻找的最近的邻居的数量。
     *
     */
    unsigned int
    get_n_nearest_neighbors() const;

    // We need these since we inherit privately to avoid polymorphic use.
    using PointPredicate::get;
    using PointPredicate::size;

  private:
    unsigned int n_nearest_neighbors;
  };



  /**
   * BoundingBox谓词的基类，为派生类提供基本功能，不应该单独使用。
   *
   */
  class BoundingBoxPredicate
  {
  protected:
    /**
     * 构造函数。  @p bounding_boxes
     * 是谓词所使用的边界盒的列表。
     *
     */
    template <int dim, typename Number>
    BoundingBoxPredicate(
      const std::vector<dealii::BoundingBox<dim, Number>> &bounding_boxes);

    /**
     * 存储在结构中的界线盒的数量。
     *
     */
    std::size_t
    size() const;

    /**
     * 返回对象中存储的第`i`个BoundingBox。
     *
     */
    const dealii::BoundingBox<3, float> &
    get(unsigned int i) const;

  private:
    std::vector<dealii::BoundingBox<3, float>> bounding_boxes;
  };



  /**
   * 这个类被 ArborXWrappers::BVH
   * 用来确定对于给定的边界盒，哪些用于建立
   * ArborXWrappers::BVH 的边界盒与它们相交。
   * @note 该类不应该在多态环境中使用。
   *
   */
  class BoundingBoxIntersectPredicate : private BoundingBoxPredicate
  {
  public:
    /**
     * 构造函数。  @p bounding_boxes
     * 是一个边界盒的列表，我们有兴趣知道它们是否与
     * ArborXWrappers::BVH  边界盒相交。
     *
     */
    template <int dim, typename Number>
    BoundingBoxIntersectPredicate(
      const std::vector<dealii::BoundingBox<dim, Number>> &bounding_boxes);

    // We need these since we inherit privately to avoid polymorphic use.
    using BoundingBoxPredicate::get;
    using BoundingBoxPredicate::size;
  };


  /**
   * 这个类被 ArborXWrappers::BVH
   * 用来确定对于给定的界线盒，哪些是用于构建
   * ArborXWrappers::BVH. 的界线盒/点中最近的。
   * @note  该类不应该在多态环境中使用。
   *
   */
  class BoundingBoxNearestPredicate : private BoundingBoxPredicate
  {
  public:
    /**
     * 构造函数。  @p bounding_boxes
     * 是一个界线盒的列表，对于这个列表，我们有兴趣知道用于建立
     * ArborXWrappers::BVH. 的最近的界线盒。
     *
     */
    template <int dim, typename Number>
    BoundingBoxNearestPredicate(
      const std::vector<dealii::BoundingBox<dim, Number>> &bounding_boxes,
      const unsigned int                                   n_nearest_neighbors);

    /**
     * 返回我们要找的最近的邻居的数量。
     *
     */
    unsigned int
    get_n_nearest_neighbors() const;

    // We need these since we inherit privately to avoid polymorphic use.
    using BoundingBoxPredicate::get;
    using BoundingBoxPredicate::size;

  private:
    unsigned int n_nearest_neighbors;
  };
} // namespace ArborXWrappers

DEAL_II_NAMESPACE_CLOSE

/**
 * 这个命名空间包含ArborX使用的AccessTraits的实现。
 *
 *
 */
namespace ArborX
{
  /**
   * 这个结构允许ArborX使用 std::vector<dealii::Point> 作为原语。
   *
   */
  template <int dim, typename Number>
  struct AccessTraits<std::vector<dealii::Point<dim, Number>>, PrimitivesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * 返回向量 @p v. 的大小。
     *
     */
    static std::size_t
    size(const std::vector<dealii::Point<dim, Number>> &v);

    /**
     * 从 dealii::Point  `v[i]`中返回一个 ArborX::Point 。
     *
     */
    static Point
    get(const std::vector<dealii::Point<dim, Number>> &v, std::size_t i);
  };



  /**
   * 这个结构允许ArborX使用 std::vector<dealii::BoundingBox>
   * 作为原语。
   *
   */
  template <int dim, typename Number>
  struct AccessTraits<std::vector<dealii::BoundingBox<dim, Number>>,
                      PrimitivesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * 返回向量的大小  @p v.  。
     *
     */
    static std::size_t
    size(const std::vector<dealii::BoundingBox<dim, Number>> &v);

    /**
     * 从 dealii::BoundingBox 中返回一个 ArborX::Box  `v[i]`。
     *
     */
    static Box
    get(const std::vector<dealii::BoundingBox<dim, Number>> &v, std::size_t i);
  };



  /**
   * 这个结构允许ArborX使用PointIntersectPredicate作为谓词。
   *
   */
  template <>
  struct AccessTraits<dealii::ArborXWrappers::PointIntersectPredicate,
                      PredicatesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * 存储在 @p pt_intersect. 中的点的数量。
     *
     */
    static std::size_t
    size(const dealii::ArborXWrappers::PointIntersectPredicate &pt_intersect);

    /**
     * 返回一个由存储在 @p pt_intersect. 中的第`i`个 dealii::Point
     * 构建的 Arbox::intersects(ArborX::Point) 对象。
     *
     */
    static auto
    get(const dealii::ArborXWrappers::PointIntersectPredicate &pt_intersect,
        std::size_t                                            i);
  };


  /**
   * 这个结构允许ArborX使用PointNearestPredicate作为谓词。
   *
   */
  template <>
  struct AccessTraits<dealii::ArborXWrappers::PointNearestPredicate,
                      PredicatesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * 存储在 @p pt_nearest. 中的点的数量。
     *
     */
    static std::size_t
    size(const dealii::ArborXWrappers::PointNearestPredicate &pt_nearest);

    /**
     * 返回一个由存储在 @p pt_nearest. 中的第`i`个 dealii::Point
     * 构建的 Arbox::nearest(ArborX::Point,
     * PointNearestPredicate::get_n_nearest_neighbors) 对象。
     *
     */
    static auto
    get(const dealii::ArborXWrappers::PointNearestPredicate &pt_nearest,
        std::size_t                                          i);
  };


  /**
   * 这个结构允许ArborX使用BoundingBoxIntersectPredicate作为谓词。
   *
   */
  template <>
  struct AccessTraits<dealii::ArborXWrappers::BoundingBoxIntersectPredicate,
                      PredicatesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * 存储在 @p bb_intersect. 中的BoundingBox的数量。
     *
     */
    static std::size_t
    size(const dealii::ArborXWrappers::BoundingBoxIntersectPredicate
           &bb_intersect);

    /**
     * 返回一个由存储在 @p bb_intersect. 中的第`i`个
     * dealii::BoundingBox 构建的 Arbox::intersects(ArborX::Box) 对象。
     *
     */
    static auto
    get(
      const dealii::ArborXWrappers::BoundingBoxIntersectPredicate &bb_intersect,
      std::size_t                                                  i);
  };


  /**
   * 这个结构允许ArborX使用BoundingBoxNearstPredicate作为谓词。
   *
   */
  template <>
  struct AccessTraits<dealii::ArborXWrappers::BoundingBoxNearestPredicate,
                      PredicatesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * 存储在 @p bb_nearest. 的BoundingBox的数量。
     *
     */
    static std::size_t
    size(const dealii::ArborXWrappers::BoundingBoxNearestPredicate &bb_nearest);

    /**
     * 返回一个由存储在 @p bb_nearest. 中的第`i`个
     * dealii::BoundingBox 构建的 Arbox::nearest(ArborX::Box,
     * BoundingBoxtNearestPredicate::get_n_nearest_neighbors) 对象。
     *
     */
    static auto
    get(const dealii::ArborXWrappers::BoundingBoxNearestPredicate &bb_nearest,
        std::size_t                                                i);
  };

  // ------------------------------- Inline ----------------------------------//

  // The implementation of AccessTraits<..., PredicatesTag> needs to be in the
  // header file otherwise the return type of auto get() cannot be determined.
  // We use auto because ArborX does not expose the type of intersects

  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::PointIntersectPredicate, PredicatesTag>::
    size(const dealii::ArborXWrappers::PointIntersectPredicate &pt_intersect)
  {
    return pt_intersect.size();
  }



  inline auto
  AccessTraits<dealii::ArborXWrappers::PointIntersectPredicate, PredicatesTag>::
    get(const dealii::ArborXWrappers::PointIntersectPredicate &pt_intersect,
        std::size_t                                            i)
  {
    const auto dealii_point = pt_intersect.get(i);
    return intersects(Point{dealii_point[0], dealii_point[1], dealii_point[2]});
  }



  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::PointNearestPredicate, PredicatesTag>::
    size(const dealii::ArborXWrappers::PointNearestPredicate &pt_nearest)
  {
    return pt_nearest.size();
  }



  inline auto
  AccessTraits<dealii::ArborXWrappers::PointNearestPredicate, PredicatesTag>::
    get(const dealii::ArborXWrappers::PointNearestPredicate &pt_nearest,
        std::size_t                                          i)
  {
    const auto dealii_point = pt_nearest.get(i);
    return nearest(Point{dealii_point[0], dealii_point[1], dealii_point[2]},
                   pt_nearest.get_n_nearest_neighbors());
  }



  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::BoundingBoxIntersectPredicate,
               PredicatesTag>::
    size(
      const dealii::ArborXWrappers::BoundingBoxIntersectPredicate &bb_intersect)
  {
    return bb_intersect.size();
  }



  inline auto
  AccessTraits<dealii::ArborXWrappers::BoundingBoxIntersectPredicate,
               PredicatesTag>::
    get(
      const dealii::ArborXWrappers::BoundingBoxIntersectPredicate &bb_intersect,
      std::size_t                                                  i)
  {
    const auto boundary_points = bb_intersect.get(i).get_boundary_points();
    const dealii::Point<3, float> min_corner = boundary_points.first;
    const dealii::Point<3, float> max_corner = boundary_points.second;

    return intersects(Box{{min_corner[0], min_corner[1], min_corner[2]},
                          {max_corner[0], max_corner[1], max_corner[2]}});
  }



  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::BoundingBoxNearestPredicate,
               PredicatesTag>::
    size(const dealii::ArborXWrappers::BoundingBoxNearestPredicate &bb_nearest)
  {
    return bb_nearest.size();
  }



  inline auto
  AccessTraits<dealii::ArborXWrappers::BoundingBoxNearestPredicate,
               PredicatesTag>::
    get(const dealii::ArborXWrappers::BoundingBoxNearestPredicate &bb_nearest,
        std::size_t                                                i)
  {
    const auto boundary_points = bb_nearest.get(i).get_boundary_points();
    const dealii::Point<3, float> min_corner = boundary_points.first;
    const dealii::Point<3, float> max_corner = boundary_points.second;

    return nearest(Box{{min_corner[0], min_corner[1], min_corner[2]},
                       {max_corner[0], max_corner[1], max_corner[2]}},
                   bb_nearest.get_n_nearest_neighbors());
  }
} // namespace ArborX

#endif

#endif


