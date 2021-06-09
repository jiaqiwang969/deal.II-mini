//include/deal.II-translator/base/function_restriction_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2021 by the deal.II authors
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

#ifndef dealii_function_restriction_h
#define dealii_function_restriction_h

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/tensor.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  /**
   * 这个类在`dim +
   * 1`维度上接收一个函数，并通过将其中一个坐标限制在给定的值上，在低一个维度上创建一个新函数。在数学上，这相当于取一个函数
   * $f = f(x, y, z)$ ，一个固定值， $Z$
   * ，并定义一个新的函数（限制） $g = g(x, y) = f(x, y, Z)$ 。
   * 使用这个类，这可以转化为
   * @code
   * Function<3> &            function             = ...
   * double                   z                    = ...
   * unsigned int             restricted_direction = 2;
   * CoordinateRestriction<2> restriction(function, restricted_direction, z);
   * @endcode
   * 限制上的`dim`维坐标从限制的（`dim +
   * 1`）坐标开始排序。特别是，这意味着如果 $y$
   * -坐标在三维中被锁定为 $Y$
   * ，那么在限制上的坐标被排序为 $(z, x)$  。    $g = g(z, x)
   * = f(x, Y, z)$  .   这与 BoundingBox::cross_section.
   * 中的惯例相同。
   *
   */
  template <int dim>
  class CoordinateRestriction : public Function<dim>
  {
  public:
    /**
     * 构造函数，接受(`dim +
     * 1`)坐标方向和传入函数应被限制的值。
     * 一个指向传入函数的指针在内部存储，所以该函数的寿命必须比创建的限制值长。
     *
     */
    CoordinateRestriction(const Function<dim + 1> &function,
                          const unsigned int       direction,
                          const double             coordinate_value);

    double
    value(const Point<dim> &point, const unsigned int component) const override;

    Tensor<1, dim>
    gradient(const Point<dim> & point,
             const unsigned int component) const override;

    SymmetricTensor<2, dim>
    hessian(const Point<dim> & point,
            const unsigned int component) const override;

  private:
    // The higher-dimensional function that has been restricted.
    const SmartPointer<const Function<dim + 1>> function;

    // The (`dim + 1`)-coordinate direction that has been restricted.
    const unsigned int restricted_direction;

    // Value of the restricted coordinate.
    const double coordinate_value;
  };



  /**
   * 该类通过将坐标值的`dim`限制在一个给定的点上，从`dim+1`维函数中创建一个一维函数。
   * 在数学上，这相当于从一个函数 $f = f(x, y, z)$ 和一个点
   * $(Y, Z)$ ，定义了一个新的函数 $g = g(x) = f(x, Y, Z)$ 。
   * 使用这个类，这可以转化为
   * @code
   * Function<3> &       function = ...
   * Point<2>            point(y, z);
   * unsigned int        open_direction = 0;
   * PointRestriction<2> restriction(function, open_direction, point);
   * @endcode
   * 点的坐标将在高维函数坐标中展开，从开放方向开始（并环绕）。特别是，如果我们限制到一个点
   * $(Z, X)$
   * ，并选择保持y方向的开放性，那么创建的限制就是函数
   * $g(y) = f(X, y, Z)$  。  这与 BoundingBox::cross_section.
   * 中的惯例一致
   *
   */
  template <int dim>
  class PointRestriction : public Function<1>
  {
  public:
    /**
     * 构造函数，接受传入函数应该被限制的点，以及哪个（`dim
     * + 1`）维坐标方向应该被保持 "开放"。
     * 一个指向传入函数的指针在内部存储，所以该函数的寿命必须比创建的限制更长。
     *
     */
    PointRestriction(const Function<dim + 1> &function,
                     const unsigned int       open_direction,
                     const Point<dim> &       point);

    double
    value(const Point<1> &point, const unsigned int component) const override;

    Tensor<1, 1>
    gradient(const Point<1> &   point,
             const unsigned int component) const override;

    SymmetricTensor<2, 1>
    hessian(const Point<1> &point, const unsigned int component) const override;

  private:
    // The higher-dimensional function that has been restricted.
    const SmartPointer<const Function<dim + 1>> function;

    // The (`dim + 1`)-coordinate direction that is kept "open"
    const unsigned int open_direction;

    // The point that we have restricted the above function to.
    const Point<dim> point;
  };

} // namespace Functions


namespace internal
{
  /**
   * 创建一个(`dim +
   * 1`)维的点，方法是复制传入的`dim`维的点的坐标，并将
   * "缺失的"(`dim + 1`)维分量设置为传入坐标值。
   * 例如，给定输入 $\{(x, y), 2, z \}$ ，该函数创建点 $(x, y,
   * z)$  。
   * `dim`维的点的坐标按照函数coordinate_to_one_dim_higher给出的约定顺序写到(`dim
   * +
   * 1`)维的点的坐标。因此，低维点上的坐标顺序没有被保留。
   * $\{(z, x), 1, y \}$  创建点  $(x, y, z)$  .
   *
   */
  template <int dim>
  Point<dim + 1>
  create_higher_dim_point(const Point<dim> & point,
                          const unsigned int component_in_dim_plus_1,
                          const double       coordinate_value);
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif  /* dealii_function_restriction_h */ 


