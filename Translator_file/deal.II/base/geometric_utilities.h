//include/deal.II-translator/base/geometric_utilities_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

#ifndef dealii_geometric_utilities_h
#define dealii_geometric_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>

#include <array>


DEAL_II_NAMESPACE_OPEN


/**
 * 一个几何实用函数的命名空间，这些函数并不特别针对有限元计算或数值程序，但在编写应用程序时，在各种情况下都需要。
 *
 *
 * @ingroup utilities
 *
 *
 */
namespace GeometricUtilities
{
  /**
   * 一个用于坐标变换的命名空间。
   *
   */
  namespace Coordinates
  {
    /**
     * 返回笛卡尔点的球面坐标  @p point.
     * 返回的数组充满了半径、方位角  $\in [0,2 \pi)$
     * 和极地/倾角  $ \in [0,\pi]$  （在2D中省略）。
     * 在三维中，转换的方式是
     * @f{align*}{
     * r &= \sqrt{x^2+y^2+z^2} \\
     * \theta &= {\rm atan}(y/x) \\
     * \phi &= {\rm acos} (z/r)
     * @f}
     * 该函数的使用在 step-75 中演示。
     *
     */
    template <int dim>
    std::array<double, dim>
    to_spherical(const Point<dim> &point);

    /**
     * 返回由 @p scoord
     * 定义的球面点的直角坐标，该球面点的半径为 $r \in
     * [0,\infty)$ ，方位角为 $\theta \in [0,2 \pi)$ ，极角/倾角为
     * $\phi \in [0,\pi]$ （2D中省略）。
     * 在三维中，转换是由
     * @f{align*}{
     * x &= r\, \cos(\theta) \, \sin(\phi) \\
     * y &= r\, \sin(\theta) \, \sin(\phi) \\
     * z &= r\, \cos(\phi)
     * @f}
     *
     */
    template <std::size_t dim>
    Point<dim>
    from_spherical(const std::array<double, dim> &scoord);

  } // namespace Coordinates
} // namespace GeometricUtilities

DEAL_II_NAMESPACE_CLOSE

#endif


