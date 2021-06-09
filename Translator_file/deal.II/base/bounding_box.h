//include/deal.II-translator/base/bounding_box_0.txt
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

#ifndef dealii_base_bounding_box_h
#define dealii_base_bounding_box_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 枚举器NeighborType描述了两个界线盒之间的相邻关系。
 *
 *
 */
enum class NeighborType
{
  /**
   * 不相邻：相交处为空。
   *
   */
  not_neighbors = 0,

  /**
   * 简单相邻：盒子与一个维度最多为`spacedim的交点相交。
   *
   * - 2`. 例如，在2d中，这意味着两个盒子在每个盒子的一个角上相接触。
   *
   */
  simple_neighbors = 1,

  /**
   * 附着的邻居：与`维度>spacedim的相交的邻居
   *
   * - 2`. 例如，在2d中，这意味着两个盒子沿着一条边接触。
   *
   */
  attached_neighbors = 2,

  /**
   * 可合并的邻居：可以用一个BoundingBox来表示的邻居，比如说
   * @code
   * .--V--W    .-----V
   * |  |  | =  |     |
   * V--W--.    V-----.
   * @endcode
   * 或者一个在另一个里面
   *
   */
  mergeable_neighbors = 3
};

/**
 * 一个表示任意尺寸<tt>spacedim</tt>且边与坐标轴平行的盒子的类，也就是一个区域
 * @f[
 * [x_0^L, x_0^U] \times ... \times [x_{spacedim-1}^L, x_{spacedim-1}^U],
 * @f]
 * 其中 $(x_0^L , ..., x_{spacedim-1}^L) and $  (x_0^U , ...,
 * x_{spacedim-1}^U)表示两个顶点（左下和右上），用于表示盒子。
 * 从几何学上看，一个边界盒就是这样。
 *
 *
 *
 * - 1D：一个线段（由其顶点按适当顺序表示
 *
 *
 *
 * - 二维：一个长方形（由左下角和右上角的顶点V代表）。
 *
 * @code
 * .--------V
 * |        |
 * V--------.
 * @endcode
 *
 *
 *
 *
 * - 三维：长方体（在这种情况下，两个顶点V遵循惯例，不属于同一个面）。
 *
 * @code
 * .------V
 * /      /|
 * .------. |
 * |      | /
 * |      |/
 * V------.
 * @endcode
 *
 * 例如，边界盒在平行分布的网格中很有用，可以对网格的每一部分的所有者进行一般描述。
 * 将BoundingBox<spacedim>的横截面与给定的方向正交，可以得到一个低一维的盒子：BoundingBox<spacedim
 *
 * - 1>. 在三维中，BoundingBox<3>横截面的两个坐标可以用两种不同的方式排序。也就是说，如果我们将横截面与y方向正交，我们可以将3D坐标排序为 $(x,z)$ 或 $(z,x)$  的2D坐标。本类使用第二种约定，对应于坐标的循环排序 $x \rightarrow y \rightarrow z \rightarrow x \rightarrow ... $  准确地说，如果我们取一个横截面。
 * 正交于 * 横截面坐标排序为 |
 * |:-------------:|:------------------------------------:| | x | (y, z) | y |
 * (z, x) | z | (x, y) 这是根据函数
 * <code>coordinate_to_one_dim_higher</code> 所设定的惯例。
 *
 *
 */
template <int spacedim, typename Number = double>
class BoundingBox
{
public:
  /**
   * 标准构造函数。创建一个对应于空盒子的对象，即一个两点都是原点的退化盒子。
   *
   */
  BoundingBox() = default;

  /**
   * 非空盒子的标准构造函数：它使用一对描述盒子的点：一个是底角，一个是顶角。
   *
   */
  BoundingBox(const std::pair<Point<spacedim, Number>, Point<spacedim, Number>>
                &boundary_points);

  /**
   * 构建包围给定容器中所有点的包围盒。
   * 该构造函数支持任何为Point<spacedim,
   * Number>元素提供begin()和end()遍历器的容器。
   *
   */
  template <class Container>
  BoundingBox(const Container &points);

  /**
   * 返回一个对边界点的引用
   *
   */
  std::pair<Point<spacedim, Number>, Point<spacedim, Number>> &
  get_boundary_points();

  /**
   * 返回一个对边界点的常量引用
   *
   */
  const std::pair<Point<spacedim, Number>, Point<spacedim, Number>> &
  get_boundary_points() const;

  /**
   * 测试是否相等。
   *
   */
  bool
  operator==(const BoundingBox<spacedim, Number> &box) const;

  /**
   * 测试不等式。
   *
   */
  bool
  operator!=(const BoundingBox<spacedim, Number> &box) const;

  /**
   * 检查当前对象和 @p other_bbox
   * 是否为邻居，也就是说，如果盒子的尺寸为spacedim，检查它们的交集是否为非空。
   * 返回一个NeighborType类型的枚举器。
   *
   */
  NeighborType
  get_neighbor_type(const BoundingBox<spacedim, Number> &other_bbox) const;

  /**
   * 扩大当前对象，使其包含  @p other_bbox  。
   * 如果当前对象已经包含 @p other_bbox
   * ，那么它不会被这个函数所改变。
   *
   */
  void
  merge_with(const BoundingBox<spacedim, Number> &other_bbox);

  /**
   * 如果点在包围盒内，返回真，否则返回假。参数 @p
   * tolerance
   * 是一个系数，相对于包围盒的尺寸来说，包围盒被放大了，以便以一种数字上稳健的方式确定该点是否在里面。
   *
   */
  bool
  point_inside(
    const Point<spacedim, Number> &p,
    const double tolerance = std::numeric_limits<Number>::epsilon()) const;

  /**
   * 按给定的数量增加（或减少）边界框的大小。
   * 调用此方法后，边界框的左下角的每个坐标将减少 @p
   * amount, ，边界框的右上角的每个坐标将增加 @p amount.
   * 如果你用一个负数调用此方法，并且原始边界框的一个轴小于
   * amount/2，此方法将触发一个断言。
   *
   */
  void
  extend(const Number &amount);

  /**
   * 计算BoundingBox的体积（即二维度量）。
   *
   */
  double
  volume() const;

  /**
   * 返回盒子中心的点。
   *
   */
  Point<spacedim, Number>
  center() const;

  /**
   * 返回盒子的边长，单位是 @p direction. 。
   *
   */
  Number
  side_length(const unsigned int direction) const;

  /**
   * 返回盒子的下限，单位为 @p direction. 。
   *
   */
  Number
  lower_bound(const unsigned int direction) const;

  /**
   * 返回 @p direction. 中的盒子的上界
   *
   */
  Number
  upper_bound(const unsigned int direction) const;

  /**
   * 返回 @p direction,
   * 中的盒子的边界，作为一个一维的盒子。
   *
   */
  BoundingBox<1, Number>
  bounds(const unsigned int direction) const;

  /**
   * 返回盒子的第索引顶点。顶点的含义与单元格的含义相同，因此，
   * @p index   $\in [0, 2^{\text{dim}}
   *
   * - 1]$  .
   *
   */
  Point<spacedim, Number>
  vertex(const unsigned int index) const;

  /**
   * 返回盒子的第1个子节点。子项的含义与单元格的含义相同。
   *
   */
  BoundingBox<spacedim, Number>
  child(const unsigned int index) const;

  /**
   * 返回正交于 @p direction.
   * 的盒子的横截面，这是一个低一维的盒子。
   * @note 在一维中调用此方法将导致一个异常，因为
   * <code>BoundingBox&lt;0&gt;</code> 没有实现。
   *
   */
  BoundingBox<spacedim - 1, Number>
  cross_section(const unsigned int direction) const;

  /**
   * 应用仿射变换，将此BoundingBox转换为一个单元BoundingBox对象。
   * 如果 $B$ 是这个边界盒，而 $\hat{B}$
   * 是单位边界盒，计算满足 $G(B) = \hat{B}$
   * 的仿生映射并将其应用于 @p point. 。
   *
   */
  Point<spacedim, Number>
  real_to_unit(const Point<spacedim, Number> &point) const;

  /**
   * 应用将单位BoundingBox对象转换为该对象的仿射变换。
   * 如果 $B$ 是这个边界盒，而 $\hat{B}$
   * 是单位边界盒，计算满足 $F(\hat{B}) = B$
   * 的仿生映射并将其应用于 @p point. 。
   *
   */
  Point<spacedim, Number>
  unit_to_real(const Point<spacedim, Number> &point) const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写入或读出到一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

private:
  std::pair<Point<spacedim, Number>, Point<spacedim, Number>> boundary_points;
};

/**
 * 这个类的存在是为了实现与尺寸无关的编程，但在其构造函数中无条件地抛出一个异常。
 *
 *
 */
template <typename Number>
class BoundingBox<0, Number>
{
public:
  /**
   * 默认构造函数。抛出一个异常。
   *
   */
  BoundingBox();

  /**
   * 等价的两点式构造函数。抛出一个异常。
   *
   */
  BoundingBox(const std::pair<Point<0, Number>, Point<0, Number>> &);

  /**
   * 等效的容器构造函数。抛出一个异常。
   *
   */
  template <class Container>
  BoundingBox(const Container &);
};


/**
 * 返回单位盒子  $[0,1]^\text{dim}$  。
 * @relates  BoundingBox
 *
 *
 */
template <int dim, typename Number = double>
BoundingBox<dim, Number>
create_unit_bounding_box();


namespace internal
{
  /**
   * 这个函数定义了一个惯例，当dim+1维中的一个坐标被锁定为一个给定值时，dim维中的坐标应该如何转换为dim+1维中的坐标。
   * 这个约定是这样的。从锁定的坐标开始，我们连续存储低维的坐标，并在越过维度时进行环绕。这种关系是，在二维中，|锁定在二维中|一维坐标|二维坐标||:------------:|:-------------:|:-------------:||x0|(a)|(x0,
   * a)||x1|(a)|(a ,
   * x1)|，在三维中，|锁定在三维|二维坐标|三维坐标||:-------------|。
   * --------------:|:--------------:| | x0 | (a, b) | (x0, a, b) | | x1 | (a,
   * b) | ( b, x1, a) | | x2 | (a, b) | ( a, b, x2) |
   * 给定一个锁定坐标，这个函数将dim维度的坐标索引映射到dim+1维度的坐标索引。
   * @param  locked_coordinate应该在[0, dim+1]范围内。    @param
   * coordinate_in_dim应该在[0, dim]范围内。    @return  在[0,
   * dim+1]范围内的一个坐标索引  @relates  BoundingBox
   *
   */
  template <int dim>
  inline int
  coordinate_to_one_dim_higher(const int locked_coordinate,
                               const int coordinate_in_dim)
  {
    AssertIndexRange(locked_coordinate, dim + 1);
    AssertIndexRange(coordinate_in_dim, dim);
    return (locked_coordinate + coordinate_in_dim + 1) % (dim + 1);
  }

} // namespace internal

 /*------------------------ Inline functions: BoundingBox --------------------*/ 

#ifndef DOXYGEN


template <int spacedim, typename Number>
inline BoundingBox<spacedim, Number>::BoundingBox(
  const std::pair<Point<spacedim, Number>, Point<spacedim, Number>>
    &boundary_points)
{
  // We check the Bounding Box is not degenerate
  for (unsigned int i = 0; i < spacedim; ++i)
    Assert(boundary_points.first[i] <= boundary_points.second[i],
           ExcMessage("Bounding Box can't be created: the points' "
                      "order should be bottom left, top right!"));

  this->boundary_points = boundary_points;
}



template <int spacedim, typename Number>
template <class Container>
inline BoundingBox<spacedim, Number>::BoundingBox(const Container &points)
{
  // Use the default constructor in case points is empty instead of setting
  // things to +oo and -oo
  if (points.size() > 0)
    {
      auto &min = boundary_points.first;
      auto &max = boundary_points.second;
      std::fill(min.begin_raw(),
                min.end_raw(),
                std::numeric_limits<Number>::infinity());
      std::fill(max.begin_raw(),
                max.end_raw(),
                -std::numeric_limits<Number>::infinity());

      for (const Point<spacedim, Number> &point : points)
        for (unsigned int d = 0; d < spacedim; ++d)
          {
            min[d] = std::min(min[d], point[d]);
            max[d] = std::max(max[d], point[d]);
          }
    }
}



template <int spacedim, typename Number>
inline std::pair<Point<spacedim, Number>, Point<spacedim, Number>> &
BoundingBox<spacedim, Number>::get_boundary_points()
{
  return this->boundary_points;
}



template <int spacedim, typename Number>
inline const std::pair<Point<spacedim, Number>, Point<spacedim, Number>> &
BoundingBox<spacedim, Number>::get_boundary_points() const
{
  return this->boundary_points;
}



template <int spacedim, typename Number>
inline bool
BoundingBox<spacedim, Number>::
operator==(const BoundingBox<spacedim, Number> &box) const
{
  return boundary_points == box.boundary_points;
}



template <int spacedim, typename Number>
inline bool
BoundingBox<spacedim, Number>::
operator!=(const BoundingBox<spacedim, Number> &box) const
{
  return boundary_points != box.boundary_points;
}



template <int spacedim, typename Number>
inline void
BoundingBox<spacedim, Number>::extend(const Number &amount)
{
  for (unsigned int d = 0; d < spacedim; ++d)
    {
      boundary_points.first[d] -= amount;
      boundary_points.second[d] += amount;
      Assert(boundary_points.first[d] <= boundary_points.second[d],
             ExcMessage("Bounding Box can't be shrunk this much: the points' "
                        "order should remain bottom left, top right."));
    }
}


template <int spacedim, typename Number>
template <class Archive>
void
BoundingBox<spacedim, Number>::serialize(Archive &ar,
                                         const unsigned int  /*version*/ )
{
  ar &boundary_points;
}



template <typename Number>
inline BoundingBox<0, Number>::BoundingBox()
{
  AssertThrow(false, ExcImpossibleInDim(0));
}



template <typename Number>
inline BoundingBox<0, Number>::BoundingBox(
  const std::pair<Point<0, Number>, Point<0, Number>> &)
{
  AssertThrow(false, ExcImpossibleInDim(0));
}



template <typename Number>
template <class Container>
inline BoundingBox<0, Number>::BoundingBox(const Container &)
{
  AssertThrow(false, ExcImpossibleInDim(0));
}



#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif


