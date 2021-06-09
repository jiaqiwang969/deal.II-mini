//include/deal.II-translator/base/geometry_info_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#ifndef dealii_geometry_info_h
#define dealii_geometry_info_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/ndarray.h>
#include <deal.II/base/point.h>
#include <deal.II/base/std_cxx20/iota_view.h>

#include <array>
#include <cstdint>



DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
namespace internal
{
  namespace GeometryInfoHelper
  {
    // A struct that holds the values for all the arrays we want to initialize
    // in GeometryInfo
    template <int dim>
    struct Initializers;

    template <>
    struct Initializers<1>
    {
      static constexpr std::array<unsigned int, 2>
      ucd_to_deal()
      {
        return {{0, 1}};
      }

      static constexpr std::array<unsigned int, 2>
      unit_normal_direction()
      {
        return {{0, 0}};
      }

      static constexpr std::array<int, 2>
      unit_normal_orientation()
      {
        return {{-1, 1}};
      }

      static constexpr std::array<Tensor<1, 1>, 2>
      unit_normal_vector()
      {
        return {{Tensor<1, 1>{{-1}}, Tensor<1, 1>{{1}}}};
      }

      static constexpr dealii::ndarray<Tensor<1, 1>, 2, 0>
      unit_tangential_vectors()
      {
        return {{{{}}, {{}}}};
      }

      static constexpr std::array<unsigned int, 2>
      opposite_face()
      {
        return {{1, 0}};
      }

      static constexpr std::array<unsigned int, 2>
      dx_to_deal()
      {
        return {{0, 1}};
      }

      static constexpr dealii::ndarray<unsigned int, 2, 1>
      vertex_to_face()
      {
        return {{{{0}}, {{1}}}};
      }
    };

    template <>
    struct Initializers<2>
    {
      static constexpr std::array<unsigned int, 4>
      ucd_to_deal()
      {
        return {{0, 1, 3, 2}};
      }

      static constexpr std::array<unsigned int, 4>
      unit_normal_direction()
      {
        return {{0, 0, 1, 1}};
      }

      static constexpr std::array<int, 4>
      unit_normal_orientation()
      {
        return {{-1, 1, -1, 1}};
      }

      static constexpr std::array<Tensor<1, 2>, 4>
      unit_normal_vector()
      {
        return {{Tensor<1, 2>{{-1., 0.}},
                 Tensor<1, 2>{{1., 0.}},
                 Tensor<1, 2>{{0., -1.}},
                 Tensor<1, 2>{{0., 1.}}}};
      }

      static constexpr dealii::ndarray<Tensor<1, 2>, 4, 1>
      unit_tangential_vectors()
      {
        return {{{{Tensor<1, 2>{{0, -1}}}},
                 {{Tensor<1, 2>{{0, 1}}}},
                 {{Tensor<1, 2>{{1, 0}}}},
                 {{Tensor<1, 2>{{-1, 0}}}}}};
      }

      static constexpr std::array<unsigned int, 4>
      opposite_face()
      {
        return {{1, 0, 3, 2}};
      }

      static constexpr std::array<unsigned int, 4>
      dx_to_deal()
      {
        return {{0, 2, 1, 3}};
      }

      static constexpr dealii::ndarray<unsigned int, 4, 2>
      vertex_to_face()
      {
        return {{{{0, 2}}, {{1, 2}}, {{0, 3}}, {{1, 3}}}};
      }
    };

    template <>
    struct Initializers<3>
    {
      static constexpr std::array<unsigned int, 8>
      ucd_to_deal()
      {
        return {{0, 1, 5, 4, 2, 3, 7, 6}};
      }

      static constexpr std::array<unsigned int, 6>
      unit_normal_direction()
      {
        return {{0, 0, 1, 1, 2, 2}};
      }

      static constexpr std::array<int, 6>
      unit_normal_orientation()
      {
        return {{-1, 1, -1, 1, -1, 1}};
      }

      static constexpr std::array<Tensor<1, 3>, 6>
      unit_normal_vector()
      {
        return {{Tensor<1, 3>{{-1, 0, 0}},
                 Tensor<1, 3>{{1, 0, 0}},
                 Tensor<1, 3>{{0, -1, 0}},
                 Tensor<1, 3>{{0, 1, 0}},
                 Tensor<1, 3>{{0, 0, -1}},
                 Tensor<1, 3>{{0, 0, 1}}}};
      }

      static constexpr dealii::ndarray<Tensor<1, 3>, 6, 2>
      unit_tangential_vectors()
      {
        return {{{{Tensor<1, 3>{{0, -1, 0}}, Tensor<1, 3>{{0, 0, 1}}}},
                 {{Tensor<1, 3>{{0, 1, 0}}, Tensor<1, 3>{{0, 0, 1}}}},
                 {{Tensor<1, 3>{{0, 0, -1}}, Tensor<1, 3>{{1, 0, 0}}}},
                 {{Tensor<1, 3>{{0, 0, 1}}, Tensor<1, 3>{{1, 0, 0}}}},
                 {{Tensor<1, 3>{{-1, 0, 0}}, Tensor<1, 3>{{0, 1, 0}}}},
                 {{Tensor<1, 3>{{1, 0, 0}}, Tensor<1, 3>{{0, 1, 0}}}}}};
      }

      static constexpr std::array<unsigned int, 6>
      opposite_face()
      {
        return {{1, 0, 3, 2, 5, 4}};
      }

      static constexpr std::array<unsigned int, 8>
      dx_to_deal()
      {
        return {{0, 4, 2, 6, 1, 5, 3, 7}};
      }

      static constexpr dealii::ndarray<unsigned int, 8, 3>
      vertex_to_face()
      {
        return {{{{0, 2, 4}},
                 {{1, 2, 4}},
                 {{0, 3, 4}},
                 {{1, 3, 4}},
                 {{0, 2, 5}},
                 {{1, 2, 5}},
                 {{0, 3, 5}},
                 {{1, 3, 5}}}};
      }
    };

    template <>
    struct Initializers<4>
    {
      static constexpr std::array<unsigned int, 16>
      ucd_to_deal()
      {
        return {{numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int}};
      }

      static constexpr std::array<unsigned int, 8>
      unit_normal_direction()
      {
        return {{0, 0, 1, 1, 2, 2, 3, 3}};
      }

      static constexpr std::array<int, 8>
      unit_normal_orientation()
      {
        return {{-1, 1, -1, 1, -1, 1, -1, 1}};
      }

      static constexpr std::array<Tensor<1, 4>, 8>
      unit_normal_vector()
      {
        return {{Tensor<1, 4>{{-1, 0, 0, 0}},
                 Tensor<1, 4>{{1, 0, 0, 0}},
                 Tensor<1, 4>{{0, -1, 0, 0}},
                 Tensor<1, 4>{{0, 1, 0, 0}},
                 Tensor<1, 4>{{0, 0, -1, 0}},
                 Tensor<1, 4>{{0, 0, 1, 0}},
                 Tensor<1, 4>{{0, 0, 0, -1}},
                 Tensor<1, 4>{{0, 0, 0, 1}}}};
      }

      static constexpr dealii::ndarray<Tensor<1, 4>, 8, 3>
      unit_tangential_vectors()
      {
        return {{{{Tensor<1, 4>{{0, -1, 0, 0}},
                   Tensor<1, 4>{{0, 0, 1, 0}},
                   Tensor<1, 4>{{0, 0, 0, 1}}}},
                 {{Tensor<1, 4>{{0, 1, 0, 0}},
                   Tensor<1, 4>{{0, 0, 1, 0}},
                   Tensor<1, 4>{{0, 0, 0, 1}}}},
                 {{Tensor<1, 4>{{0, 0, -1, 0}},
                   Tensor<1, 4>{{0, 0, 0, 1}},
                   Tensor<1, 4>{{1, 0, 0, 0}}}},
                 {{Tensor<1, 4>{{0, 0, 1, 0}},
                   Tensor<1, 4>{{0, 0, 0, 1}},
                   Tensor<1, 4>{{1, 0, 0, 0}}}},
                 {{Tensor<1, 4>{{0, 0, 0, -1}},
                   Tensor<1, 4>{{1, 0, 0, 0}},
                   Tensor<1, 4>{{0, 1, 0, 0}}}},
                 {{Tensor<1, 4>{{0, 0, 0, 1}},
                   Tensor<1, 4>{{1, 0, 0, 0}},
                   Tensor<1, 4>{{0, 1, 0, 0}}}},
                 {{Tensor<1, 4>{{-1, 0, 0, 0}},
                   Tensor<1, 4>{{0, 1, 0, 0}},
                   Tensor<1, 4>{{0, 0, 1, 0}}}},
                 {{Tensor<1, 4>{{1, 0, 0, 0}},
                   Tensor<1, 4>{{0, 1, 0, 0}},
                   Tensor<1, 4>{{0, 0, 1, 0}}}}}};
      }

      static constexpr std::array<unsigned int, 8>
      opposite_face()
      {
        return {{1, 0, 3, 2, 5, 4, 7, 6}};
      }

      static constexpr std::array<unsigned int, 16>
      dx_to_deal()
      {
        return {{numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int,
                 numbers::invalid_unsigned_int}};
      }

      static constexpr dealii::ndarray<unsigned int, 16, 4>
      vertex_to_face()
      {
        return {{{{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}},
                 {{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}},
                 {{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}},
                 {{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}},
                 {{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}},
                 {{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}},
                 {{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}},
                 {{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}},
                 {{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}},
                 {{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}},
                 {{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}},
                 {{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}},
                 {{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}},
                 {{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}},
                 {{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}},
                 {{numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int,
                   numbers::invalid_unsigned_int}}}};
      }
    };
  } // namespace GeometryInfoHelper
} // namespace internal
#endif // DOXYGEN


/**
 * 一个可以表示三角形所组成的各种对象的类：顶点、线、四边形和六边形。
 * 这个类是相当原始的：它只存储一个整数，代表所代表对象的维度。换句话说，这个类主要是作为一种传递对象的方式，它的数据类型解释了它的作用（不像只是传递一个整数），并为这些对象提供符号名称，如
 * GeometryPrimitive::vertex 而不是一个整数0。
 * 由于能够用所代表的对象的积分维度来识别这些对象，这个类提供了与无符号整数的转换操作符。
 *
 *
 */
class GeometryPrimitive
{
public:
  /**
   * 一个枚举，为可由该类表示的对象提供符号名称。这些符号名称的数值等于代表对象的几何维度，以使整数变量的转换更加简单。
   *
   */
  enum Object
  {
    /**
     * 一个顶点。
     *
     */
    vertex = 0,
    /**
     * 一条线。
     *
     */
    line = 1,
    /**
     * 一个四边形。
     *
     */
    quad = 2,
    /**
     * 一个六面体。
     *
     */
    hex = 3
  };

  /**
   * 构造函数。用给定的参数初始化对象，代表一个顶点、线等。
   *
   */
  GeometryPrimitive(const Object object);

  /**
   * 构造函数。用一个整数初始化对象，该整数应代表有关几何对象的尺寸。这通常是一个介于零（顶点）和三（六面体）之间的数字。
   *
   */
  GeometryPrimitive(const unsigned int object_dimension);

  /**
   * 返回当前表示的对象的积分维度，即0代表顶点，1代表线，等等。
   *
   */
  operator unsigned int() const;

private:
  /**
   * 当前代表的对象。
   *
   */
  Object object;
};



/**
 * 一个提供当前空间维度下各向同性和各向异性细化标志的可能选择的类。
 * 除了在一些奇怪的模板结构中，这个一般的模板是不用的。然而，实际使用的是专业化的
 * <code>RefinementPossibilities@<1@></code> ,
 * <code>RefinementPossibilities@<2@></code>  , 和
 * <code>RefinementPossibilities@<3@></code>  。
 *
 *
 * @ingroup aniso
 *
 *
 */
template <int dim>
struct RefinementPossibilities
{
  /**
   * 在当前维度的细化情况下，可能的数值。
   * 注意数值的构造：最低位描述的是X轴的切面，第二至最低位对应的是Y轴的切面，第三至最低位对应的是Z轴的切面。因此，以下关系成立（除其他外）。
   * @code
   * cut_xy  == cut_x  | cut_y
   * cut_xyz == cut_xy | cut_xz
   * cut_x   == cut_xy & cut_xz
   * @endcode
   * 当然，只提供那些在特定空间维度上合理的切割。
   * 此外，标签 <code>isotropic_refinement</code>
   * 表示在这个类的模板参数所选择的空间维度上的各向同性细化。
   * 如果你选择了各向异性的细化，例如通过传递
   * CellIterator::set_refine_flag() 的参数 RefinementCase::cut_x,
   * RefinementCase::cut_y,  RefinementCase::cut_z,
   * 中的一个标志或这些标志的组合，那么请记住，X、Y或Z方向的细化是针对单元的
   * <em> 局部 </em>
   * 坐标系进行的。换句话说，这些标志决定了单元格的哪些边和面将被切割成新的边和面。另一方面，这个过程与细胞在
   * <em> 全局 </em>
   * 坐标系中的方向无关，你不应该假定细胞的局部坐标系在它所在空间的全局坐标系中的任何特定方向。
   *
   */
  enum Possibilities
  {
    /**
     * 不要进行细化。
     *
     */
    no_refinement = 0,

    /**
     * 执行各向同性的细化。这意味着在所有坐标方向进行细化。对于当前的一般模板类
     *
     * - 因为有针对1d、2d和3d情况的特殊性，所以从未使用过
     *
     * - ，我们简单地将这个数字设置为一个有所有比特设置的值。RefinementPossibilities<1>、RefinementPossibilities<2>和RefinementPossibilities<3>中的特殊化将相应的`enum`元素设置为更合理的值。
     *
     */
    isotropic_refinement = 0xFF
  };
};



/**
 * 一个提供当前空间维度下各向同性和各向异性细化标志的可能选择的类。
 * 这个特化用于 <code>dim=1</code> ，它提供X方向的细化。
 *
 *
 * @ingroup aniso
 *
 *
 */
template <>
struct RefinementPossibilities<1>
{
  /**
   * 在当前维度的细化情况下，可能的数值。
   * 注意数值的构造：最低位描述的是x轴的切面，第二至最低位对应的是y轴的切面，第三至最低位对应的是z轴的切面。因此，以下关系成立（除其他外）。
   * @code
   * cut_xy  == cut_x  | cut_y
   * cut_xyz == cut_xy | cut_xz
   * cut_x   == cut_xy & cut_xz
   * @endcode
   * 当然，只提供那些在特定空间维度上合理的切割。
   * 此外，标签 <code>isotropic_refinement</code>
   * 表示在这个类的模板参数所选择的空间维度上的各向同性细化。
   * 如果你选择各向异性的细化，例如通过传递
   * CellIterator::set_refine_flag() 的一个参数作为标志
   * RefinementCase::cut_x,   RefinementCase::cut_y,   RefinementCase::cut_z,
   * 或这些标志的组合，那么请记住，X、Y或Z方向的细化是在单元的
   * <em> 局部 </em>
   * 坐标系中进行。换句话说，这些标志决定了单元格的哪些边和面将被切割成新的边和面。另一方面，这个过程与细胞在
   * <em> 全局 </em>
   * 坐标系中的方向无关，你不应该假定细胞的局部坐标系在它所处空间的全局坐标系中的任何特定方向。
   *
   */
  enum Possibilities
  {
    /**
     * 不要细化。
     *
     */
    no_refinement = 0,
    /**
     * 在X方向上进行切割。
     *
     */
    cut_x = 1,
    /**
     * 进行各向同性的精细化处理。
     *
     */
    isotropic_refinement = cut_x
  };
};



/**
 * 一个提供在当前空间维度上各向同性和各向异性细化标志的可能选择的类。
 * 这个特殊化用于 <code>dim=2</code>
 * ，它提供了在X方向和Y方向分别进行细化，以及同时在两个方向进行各向同性的细化。
 *
 *
 * @ingroup aniso
 *
 *
 */
template <>
struct RefinementPossibilities<2>
{
  /**
   * 在当前维度的细化情况下，可能的数值。
   * 注意数值的构造：最低位描述的是x轴的切面，第二至最低位对应的是y轴的切面，第三至最低位对应的是z轴的切面。因此，以下关系成立（除其他外）。
   * @code
   * cut_xy  == cut_x  | cut_y
   * cut_xyz == cut_xy | cut_xz
   * cut_x   == cut_xy & cut_xz
   * @endcode
   * 当然，只提供那些在特定空间维度上合理的切割。
   * 此外，标签 <code>isotropic_refinement</code>
   * 表示在这个类的模板参数所选择的空间维度上的各向同性细化。
   * 如果你选择了各向异性的细化，例如通过传递
   * CellIterator::set_refine_flag() 的参数 RefinementCase::cut_x,
   * RefinementCase::cut_y,  RefinementCase::cut_z,
   * 中的一个标志或这些标志的组合，那么请记住，X、Y或Z方向的细化是针对单元的
   * <em> 局部 </em>
   * 坐标系进行的。换句话说，这些标志决定了单元格的哪些边和面将被切割成新的边和面。另一方面，这个过程与细胞在
   * <em> 全局 </em>
   * 坐标系中的方向无关，你不应该假定细胞的局部坐标系在它所处空间的全局坐标系中的任何特定方向。
   *
   */
  enum Possibilities
  {
    /**
     * 不要细化。
     *
     */
    no_refinement = 0,
    /**
     * 在X方向上进行切割。
     *
     */
    cut_x = 1,
    /**
     * 在Y方向上进行切割。
     *
     */
    cut_y = 2,
    /**
     * 在x方向和y方向上进行切割。
     *
     */
    cut_xy = cut_x | cut_y,

    /**
     * 执行各向同性的细化。
     *
     */
    isotropic_refinement = cut_xy
  };
};



/**
 * 一个提供在当前空间维度上各向同性和各向异性细化标志的可能选择的类。
 * 这个特殊化用于 <code>dim=3</code>
 * ，它提供了在x、y和z方向上分别进行细化，以及这些细化的组合和同时在所有方向上进行各向同性的细化。
 *
 *
 * @ingroup aniso
 *
 *
 */
template <>
struct RefinementPossibilities<3>
{
  /**
   * 在当前维度的细化情况下，可能的数值。
   * 注意数值的构造：最低位描述的是x轴的切面，第二至最低位对应的是y轴的切面，第三至最低位对应的是z轴的切面。因此，以下关系成立（除其他外）。
   * @code
   * cut_xy  == cut_x  | cut_y
   * cut_xyz == cut_xy | cut_xz
   * cut_x   == cut_xy & cut_xz
   * @endcode
   * 当然，只提供那些在特定空间维度上合理的切割。
   * 此外，标签 <code>isotropic_refinement</code>
   * 表示在这个类的模板参数所选择的空间维度上的各向同性细化。
   * 如果你选择了各向异性的细化，例如通过传递
   * CellIterator::set_refine_flag() 的参数 RefinementCase::cut_x,
   * RefinementCase::cut_y,  RefinementCase::cut_z,
   * 中的一个标志或这些标志的组合，那么请记住，X、Y或Z方向的细化是在单元的
   * <em> 局部 </em>
   * 坐标系下进行的。换句话说，这些标志决定了单元格的哪些边和面将被切割成新的边和面。另一方面，这个过程与细胞在
   * <em> 全局 </em>
   * 坐标系中的方向无关，你不应该假定细胞的局部坐标系在它所处空间的全局坐标系中的任何特定方向。
   *
   */
  enum Possibilities
  {
    /**
     * 不要细化。
     *
     */
    no_refinement = 0,
    /**
     * 在X方向上进行切割。
     *
     */
    cut_x = 1,
    /**
     * 在Y方向上进行切割。
     *
     */
    cut_y = 2,
    /**
     * 在x和y方向上进行切割。
     *
     */
    cut_xy = cut_x | cut_y,
    /**
     * 在Z方向上进行切割。
     *
     */
    cut_z = 4,
    /**
     * 在x方向和y方向上进行切割。
     *
     */
    cut_xz = cut_x | cut_z,
    /**
     * 在X和Y方向上进行切割。
     *
     */
    cut_yz = cut_y | cut_z,
    /**
     * 在X、Y和Z方向上进行切割。
     *
     */
    cut_xyz = cut_x | cut_y | cut_z,

    /**
     * 进行各向同性的精细化处理。
     *
     */
    isotropic_refinement = cut_xyz
  };
};



/**
 * 一个存储具有 <code>dim</code>
 * 维度的对象的可能的各向异性和各向同性的细化情况的类（例如，对于一条直线
 * <code>dim=1</code> 在任何空间维度，对于一个四边形
 * <code>dim=2</code>
 * ，等等）。这个类的可能值是在基类中声明的枚举中列出的；更多信息请看那里。
 *
 *
 * @ingroup aniso
 *
 *
 */
template <int dim>
class RefinementCase : public RefinementPossibilities<dim>
{
public:
  /**
   * 默认构造函数。用no_refinement初始化细化情况。
   *
   */
  RefinementCase();

  /**
   * 构造函数。从基类中指定的可能的细化列表中获取并存储一个表示特定细化的值。
   *
   */
  RefinementCase(
    const typename RefinementPossibilities<dim>::Possibilities refinement_case);

  /**
   * 构造函数。取并存储一个表示特定细化的值，作为一个位域。为了避免隐式转换到积分值或从积分值转换，这个构造函数被标记为显式。
   *
   */
  explicit RefinementCase(const std::uint8_t refinement_case);

  /**
   * 返回这个类所存储的数字值。虽然这个操作符的存在看起来很危险，但在人们希望有类似<tt>switch
   * (refinement_flag)... case  RefinementCase<dim>::cut_x:  ...
   * 的代码的情况下，它是很有用的。</tt>，它可以写成
   * <code>switch (static_cast@<std::uint8_t@>(refinement_flag)</code>  。
   * 另一个应用是使用当前类型的对象作为数组的索引；然而，这种用法已被废弃，因为它假定了从RefinementPossibilities基类中定义的符号标志到实际数值（数组索引）的某种映射。
   *
   */
  operator std::uint8_t() const;

  /**
   * 返回当前对象所代表的细化标志和作为参数给出的细化标志的联合。
   *
   */
  RefinementCase
  operator|(const RefinementCase &r) const;

  /**
   * 返回当前对象所代表的细化标志和作为参数给出的细化标志的交集。
   *
   */
  RefinementCase operator&(const RefinementCase &r) const;

  /**
   * 返回当前对象所代表的细化标志的否定值。例如，在2d中，如果当前对象持有标志
   * <code>cut_x</code>, then the returned value will be <code>cut_y</code>
   * ；如果当前值是 <code>isotropic_refinement</code>
   * ，那么结果将是 <code>no_refinement</code> ；等等。
   *
   */
  RefinementCase
  operator~() const;


  /**
   * 返回与沿参数给出的轴线切割单元格相对应的标志。例如，如果
   * <code>i=0</code> ，则返回值为
   * <tt>RefinementPossibilities<dim>::cut_x</tt>. 。
   *
   */
  static RefinementCase
  cut_axis(const unsigned int i);

  /**
   * 返回该类型的对象所占用的内存量。
   *
   */
  static std::size_t
  memory_consumption();

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

  /**
   * 异常情况。
   *
   */
  DeclException1(
    ExcInvalidRefinementCase,
    int,
    << "The refinement flags given (" << arg1
    << ") contain set bits that do not "
    << "make sense for the space dimension of the object to which they are applied.");

private:
  /**
   * 将细化情况存储为一个比特字段，在任何给定的维度上都有必要的比特。
   *
   */
  std::uint8_t value : (dim > 0 ? dim : 1);
};


namespace internal
{
  /**
   * 一个提供一个面的所有可能情况的类（在当前空间维度
   * @p dim) 可能被细分为子面。对于 <code>dim=1</code> and
   * <code>dim=2</code> ，它们对应于
   * <code>RefinementPossibilities@<dim-1@></code>
   * 中给出的情况。然而， <code>SubfacePossibilities@<3@></code>
   * 包括 <code>RefinementPossibilities@<2@></code>
   * 的细化情况，但另外还有一些子面的可能性，一个面可能被细分为，这是在两个相邻单元中的一个上重复进行的各向异性细化步骤。
   * 除了一些奇怪的模板结构，这个一般的模板是不用的。
   * 然而，实际情况是，对 <code>SubfacePossibilities@<1@></code> 、
   * <code>SubfacePossibilities@<2@></code> 和
   * <code>SubfacePossibilities@<3@></code> 进行了专门化。
   * @ingroup aniso
   *
   */
  template <int dim>
  struct SubfacePossibilities
  {
    /**
     * 面被细分为子面的可能情况。
     *
     */
    enum Possibilities
    {
      /**
       * 不要细化。
       *
       */
      case_none = 0,

      /**
       * 各向同性地进行细化。
       *
       */
      case_isotropic = static_cast<std::uint8_t>(-1)
    };
  };


  /**
   * 一个提供所有可能情况的面（在当前空间维度 @p dim)
   * 中可能被细分为子面的类。    对于 <code>dim=0</code>
   * 我们只提供一个假的实现。
   * @ingroup aniso
   *
   */
  template <>
  struct SubfacePossibilities<0>
  {
    /**
     * 面被细分为子面的可能情况。        假的实现。
     *
     */
    enum Possibilities
    {
      /**
       * 不要细化。
       *
       */
      case_none = 0,

      /**
       * 各向同性地进行细化。
       *
       */
      case_isotropic = case_none
    };
  };



  /**
   * 一个提供所有可能情况的面（在当前空间维度 @p dim)
   * 可能被细分为子面的类。    对于 <code>dim=1</code>
   * ，没有面。因此，没有子面的可能性。
   * @ingroup aniso
   *
   */
  template <>
  struct SubfacePossibilities<1>
  {
    /**
     * 面被细分为子面的可能情况。
     * 在1d中没有面，因此没有子面的可能性。
     *
     */
    enum Possibilities
    {
      /**
       * 不要细化。
       *
       */
      case_none = 0,

      /**
       * 各向同性地进行细化。
       *
       */
      case_isotropic = case_none
    };
  };



  /**
   * 一个提供一个面（在当前空间维度 @p dim)
   * 中）可能被细分为子面的所有可能情况的类。
   * 这个特殊化用于 <code>dim=2</code>
   * ，它提供了以下可能性：一个面（线）被细化（
   * <code>case_x</code>) or not refined (<code>case_no</code> ）。
   * @ingroup aniso
   *
   */
  template <>
  struct SubfacePossibilities<2>
  {
    /**
     * 面被细分为子面的可能情况。
     * 在2D中，有以下可能性：一个面（线）被细化（
     * <code>case_x</code>) or not refined (<code>case_no</code>  ）。
     *
     */
    enum Possibilities
    {
      /**
       * 不进行细化。
       *
       */
      case_none = 0,
      /**
       * 在X方向上切割。
       *
       */
      case_x = 1,
      /**
       * 各向同性地细化。
       *
       */
      case_isotropic = case_x
    };
  };



  /**
   * 一个提供一个面（在当前空间维度 @p dim)
   * 中）可能被细分为子面的所有可能情况的类。
   * 这个特殊化用于dim=3，它提供了以下可能性：一个面（四边形）在x方向或y方向（在面-内坐标系中）分别被细化，（
   * <code>case_x</code>  或（  <code>case_y</code>), and in both directions
   * (<code>case_x</code>  相当于（  <code>case_isotropic</code>
   * ）。此外，它还提供了一个面通过在两个相邻单元中的一个单元上执行的重复各向异性细化步骤所能得到的可能性。
   * 例如，有可能一个面（四边形）用 <code>cut_x</code>
   * 进行细化，之后左边的孩子又用 <code>cut_y</code>
   * 进行细化，这样就有了三个活跃的子面。然而，请注意，只允许在细化的情况下，两个六边形之间的面的每条线不超过一个悬挂节点。此外，不允许两个相邻的六面体被细化，使得其中一个六面体用
   * <code>cut_x</code> 细化共同的面，而另一个六面体用
   * <code>cut_y</code> 细化该面。事实上，
   * Triangulation::prepare_coarsening_and_refinement
   * 照顾到了这种情况，并确保精炼单元的每个面都完全包含在相邻单元的一个面中。
   * 下图解释了SubfacePossibilities并给出了相应的子面编号。
   * @code
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * -------*
   * |       |
   * |   0   |    case_none
   * |       |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * -------*
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*
   * |   |   |
   * | 0 | 1 |    case_x
   * |   |   |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*
   * | 1 |   |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---* 2 |    case_x1y
   * | 0 |   |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*
   * |   | 2 |
   * | 0---*    case_x2y
   * |   | 1 |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*
   * | 1 | 3 |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*    case_x1y2y   (successive refinement: first cut_x, then cut_y for both children)
   * | 0 | 2 |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * -------*
   * |   1   |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * -------*    case_y
   * |   0   |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * -------*
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * -------*
   * |   2   |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*    case_y1x
   * | 0 | 1 |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*
   * | 1 | 2 |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*    case_y2x
   * |   0   |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * -------*
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*
   * | 2 | 3 |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*    case_y1x2x   (successive refinement: first cut_y, then cut_x for both children)
   * | 0 | 1 |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---+---*
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*
   * | 2 | 3 |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*    case_xy      (one isotropic refinement step)
   * | 0 | 1 |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * ---*---*
   *
   * @endcode
   *
   * @ingroup aniso
   *
   */
  template <>
  struct SubfacePossibilities<3>
  {
    /**
     * 面被细分为子面的可能情况。
     * 关于子面可能性的更多细节，请参见SubfacePossibilities<3>的文档。
     *
     */
    enum Possibilities
    {
      case_none  = 0,
      case_x     = 1,
      case_x1y   = 2,
      case_x2y   = 3,
      case_x1y2y = 4,
      case_y     = 5,
      case_y1x   = 6,
      case_y2x   = 7,
      case_y1x2x = 8,
      case_xy    = 9,

      case_isotropic = case_xy
    };
  };



  /**
   * 一个提供一个面（在当前空间维度 @p dim)
   * ）可能被细分为子面的所有可能情况的类。
   * @ingroup aniso
   *
   */
  template <int dim>
  class SubfaceCase : public SubfacePossibilities<dim>
  {
  public:
    /**
     * 构造函数。在基类指定的可能情况列表中，取并存储一个表示特定子面可能性的值。
     *
     */
    SubfaceCase(const typename SubfacePossibilities<dim>::Possibilities
                  subface_possibility);

    /**
     * 返回该类所存储的数值。虽然这个操作符的存在看起来很危险，但在人们希望有<code>switch
     * (subface_case)... case  SubfaceCase::case_x:  ...
     * 这样的代码的情况下，它是很有用的。</code>，可以写成<code>switch
     * (static_cast@<std::uint8_t@>(subface_case)</code>.
     * 另一个应用是使用当前类型的对象作为数组的索引；然而，这种用法已被废弃，因为它假定了从SubfacePossibilities基类中定义的符号标志到实际数值（数组索引）的某种映射。
     *
     */
    operator std::uint8_t() const;

    /**
     * 返回该类型的对象所占用的内存量。
     *
     */
    static constexpr std::size_t
    memory_consumption();

    /**
     * 异常情况。
     *
     */
    DeclException1(
      ExcInvalidSubfaceCase,
      int,
      << "The subface case given (" << arg1 << ") does not make sense "
      << "for the space dimension of the object to which they are applied.");

  private:
    /**
     * 将细化情况存储为一个比特字段，在任何给定的维度上都有必要的比特。
     *
     */
    std::uint8_t value : (dim == 3 ? 4 : 1);
  };

} // namespace internal



template <int dim>
struct GeometryInfo;



/**
 * 这个类提供了一个零维单元的描述。它已经被ReferenceCell类所取代。
 *
 * - 更多信息见那里。
 * 零维单元的拓扑描述，即点。这个类可能看起来不是很有用，但是如果在某个维度上，我们想查询比现在低一个维度的对象的信息，例如关于面的信息，往往是有用的。
 * 这个类包含了 @p dim-dimensional
 * 网格单元的顶点和面的信息作为静态成员。该接口对所有维度都是一样的。如果一个值在低维单元中没有用处，它将被（正确地）设置为零，例如1d中的#max_children_per_face。
 * 这个信息应该总是取代顶点、邻居等的硬编码数字，因为它可以独立使用维度。
 *
 *
 * @ingroup grid geomprimitives aniso
 *
 *
 */
template <>
struct GeometryInfo<0>
{
  /**
   * 一个单元格的最大子女数，即一个各向同性细化的单元格的子女数。
   * 如果一个单元格是各向异性精炼的，实际的子代数可能少于这里给出的值。
   *
   */
  static constexpr unsigned int max_children_per_cell = 1;

  /**
   * 一个单元有多少个面。
   *
   */
  static constexpr unsigned int faces_per_cell = 0;

  /**
   * 返回一个对象，它可以被认为是一个包含从零到`面孔_per_cell`所有索引的数组。这允许使用基于范围的for循环来编写以下类型的代码。
   * @code
   * for (auto &cell : triangulation.active_cell_iterators())
   *   for (auto face_index : GeometryInfo<dim>::face_indices())
   *     if (cell->face(face_index)->at_boundary())
   *       ... do something ...
   * @endcode
   * 这里，我们在所有单元格的所有面进行循环，`face_index`使用所有有效的索引。
   * 当然，由于这个类是针对`dim==0`的情况，返回的对象实际上是一个空数组。
   *
   */
  static std::array<unsigned int, 0>
  face_indices();

  /**
   * 一个精炼面的最大子女数，即一个各向同性的精炼面的子女数。
   * 如果一个单元被各向异性地细化，实际的子代数可能少于这里给出的值。
   *
   */
  static constexpr unsigned int max_children_per_face = 0;

  /**
   * 返回用<tt>ref_case</tt>精炼的单元格（或面）的子数。因为我们在这里关注的是点，所以子节点的数量等于1。
   *
   */
  static unsigned int
  n_children(const RefinementCase<0> &refinement_case);

  /**
   * 一个单元格的顶点数。
   *
   */
  static constexpr unsigned int vertices_per_cell = 1;

  /**
   * 返回一个对象，它可以被认为是一个包含从零到`顶点_per_cell'所有索引的数组。这允许使用基于范围的for循环来编写以下类型的代码。
   * @code
   * for (auto &cell : triangulation.active_cell_iterators())
   *   for (auto vertex_index : GeometryInfo<dim>::vertex_indices())
   *     if (cell->vertex(vertex_index) satisfies some condition)
   *       ... do something ...
   * @endcode
   * 这里，我们在所有单元格的所有顶点上循环，`vertex_index`使用所有有效的索引。
   * 当然，由于这个类是针对`dim==0`的情况，返回的对象是一个只有一个条目的数组：0。这是因为维度为0的对象实际上只是一个单点，对应于一个顶点本身。
   *
   */
  static std::array<unsigned int, vertices_per_cell>
  vertex_indices();

  /**
   * 将面的顶点编号映射到单元的顶点编号，即给出面的<tt>顶点</tt>第1个顶点的单元顶点编号，例如 <tt>GeometryInfo<2>::face_to_cell_vertices(3,0)=2</tt>, 见本类文档2d部分N4点下的图片。    通过<tt>face_orientation</tt>, <tt>face_flip</tt> 和 <tt>face_rotation</tt> 参数，这个函数可以处理以标准和非标准方向的面。<tt>face_orientation</tt>默认为<tt>true</tt>，<tt>face_flip</tt>和<tt>face_rotation</tt>默认为<tt>false</tt>（标准方向）。在2d中，只有<tt>face_flip</tt>被考虑。更多信息请参见 @ref GlossFaceOrientation 的 "词汇表 "
   * 文章。
   * 由于单元格的子节点是根据单元格的顶点排序的，这个调用被传递给child_cell_on_face()函数。
   * 因此，这个函数只是child_cell_on_face()的一个包装器，给它起了一个暗示性的名字。
   * 当然，由于这个类是针对`dim==0'的情况，这个函数没有被实现。
   *
   */
  static unsigned int
  face_to_cell_vertices(const unsigned int face,
                        const unsigned int vertex,
                        const bool         face_orientation = true,
                        const bool         face_flip        = false,
                        const bool         face_rotation    = false);

  /**
   * 将面的行数映射到单元格行数，即给出面的<tt>行</tt>第1行的单元格行数，例如：
   * <tt>GeometryInfo<3>::face_to_cell_lines(5,0)=4</tt>.
   * 通过<tt>面的方向</tt>，<tt>面的翻转</tt>和<tt>面的旋转</tt>参数，这个函数处理以标准和非标准方向为方向的面。<tt>face_orientation</tt>默认为<tt>true</tt>，<tt>face_flip</tt>和<tt>face_rotation</tt>默认为<tt>false</tt>（标准方向），在2d中没有影响。
   * 当然，由于这个类是针对`dim==0'的情况，这个函数没有实现。
   *
   */
  static unsigned int
  face_to_cell_lines(const unsigned int face,
                     const unsigned int line,
                     const bool         face_orientation = true,
                     const bool         face_flip        = false,
                     const bool         face_rotation    = false);

  /**
   * 每个面所拥有的顶点数量。由于这在一维中没有用，我们提供了一个无用的数字（希望编译器在看到类似<tt>for
   * (i=0; i<vertices_per_face;
   * ++i)</tt>的结构时可以发出警告，至少在 @p i
   * 是一个<tt>无符号int</tt>的情况下。
   *
   */
  static constexpr unsigned int vertices_per_face = 0;

  /**
   * 每个面有多少条线。
   *
   */
  static constexpr unsigned int lines_per_face = 0;

  /**
   * 每个面的四边形的数量。
   *
   */
  static constexpr unsigned int quads_per_face = 0;

  /**
   * 一个单元格的线数。
   *
   */
  static constexpr unsigned int lines_per_cell = 0;

  /**
   * 一个单元格的四边形的数量。
   *
   */
  static constexpr unsigned int quads_per_cell = 0;

  /**
   * 一个单元格的六面体的数量。
   *
   */
  static constexpr unsigned int hexes_per_cell = 0;

  /**
   * 为UCD输出重新排列顶点。
   * 对于以UCD格式写入的单元格，该字段的每个条目包含<tt>deal.II</tt>中与该位置的UCD编号相对应的一个顶点的编号。
   * 典型的例子：写一个单元并安排顶点，这样UCD就能理解它们。
   * @code
   * for (i=0; i< n_vertices; ++i)
   * out << cell->vertex(ucd_to_deal[i]);
   * @endcode
   * 由于deal.II版本<=5.1中的顶点编号恰好与UCD的编号一致，这个字段也可以像old_to_lexicographic映射一样使用。
   *
   */
  static const std::array<unsigned int, vertices_per_cell> ucd_to_deal;

  /**
   * 为OpenDX输出重新排列顶点。
   * 对于一个被写成OpenDX格式的单元格，这个字段的每个条目都包含<tt>deal.II</tt>中的一个顶点的编号，它与这个位置的DX编号相对应。
   * 典型的例子：写一个单元并安排顶点，这样OpenDX就能理解它们。
   * @code
   * for (i=0; i< n_vertices; ++i)
   * out << cell->vertex(dx_to_deal[i]);
   * @endcode
   *
   *
   */
  static const std::array<unsigned int, vertices_per_cell> dx_to_deal;
};



/**
 * 该类为构成单元的所有拓扑结构提供独立的维度信息，或 @ref GlossReferenceCell "参考单元"
 * 。这个类已经被ReferenceCell类所取代。
 *
 * - 更多信息见那里。
 *
 * 它是库中的一个中心点，关于参考单元的顶点、线或面的编号的信息被收集起来。因此，这个类的信息被广泛用于Triangulation对象的几何描述，以及代码的其他各个部分。特别是，它也是以独立于维度的方式编写代码的重点；例如，在2D中不写顶点0<=v<4的循环，而是写成
 * 0<=v<GeometryInfo<dim>::vertices_per_cell,
 * ，从而使代码在3D中也能工作，而无需改变。
 * 该类中最常用的部分是它的静态成员，如vertices_per_cell、faces_per_cell等。然而，该类也提供了关于更抽象的问题的信息，如面的方向等。下面的文档给出了许多这些概念的文字描述。
 *
 *  <h3>Implementation conventions for two spatial dimensions</h3>
 * 从5.2版本开始，deal.II基于一个编号方案，尽可能使用lexicographic排序（x跑得最快），因此试图采用一种'经典'排序。
 * 2d中顶点和面（线）的排序定义为
 *
 *
 * - 顶点的编号是按词法排序的
 *
 * - 面（2d中的线）：首先是在x方向和y方向上有法线的两个面。对于每两个面：首先是法线为负坐标方向的面，然后是法线为正坐标方向的面，也就是说，面的排序是根据它们的法线指向
 *
 * - , x,
 *
 * - , y方向的法线排序。
 *
 *
 *
 * - 线条的方向由点0对点1的方向表示，并且总是在其中一个坐标方向上
 *
 *
 *
 * - 3d中的面线是有序的，这样诱导的2D局部坐标系（x,y）意味着（右手规则）面线法线方向的法线，见N2/。
 * 由此产生的2d中顶点和面（线）的编号，以及线的方向，如下所示。
 *
 * @verbatim
 *     3
 *  2-->--3
 *  |     |
 * 0^     ^1
 *  |     |
 *  0-->--1
 *      2
 * @endverbatim
 *
 * 请注意，在构建网格时，线的方向必须是正确的；但是，在细化时，它将自动被保留。
 * 此外，我们定义子线的方向与父线相同，即<tt>line->child(0)->vertex(0)==line->vertex(0)</tt>，<tt>line->child(1)->vertex(1)==line->vertex(1)</tt>。这也意味着，第一条子线（<tt>line->child(0)</tt>）是位于旧线的顶点(0)的那条。
 * 同样地，我们定义，一个四边形的四个子线与旧四边形的顶点相邻，并具有相同的编号。
 * 请注意，关于这些约定的几个信息可以在运行或编译时从本类的成员函数和变量中提取出来。
 *
 *  <h4>Coordinate systems</h4>
 * 当单元格中的点需要明确的坐标时（例如用于正交公式或试用函数的定义点），我们为单元格定义以下坐标系统。
 *
 * @verbatim
 * y^   2-----3
 * |   |     |
 * |   |     |
 * |   |     |
 * |   0-----1
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ------------>x
 * @endverbatim
 *
 * 这里，顶点0是坐标系的原点，顶点1的坐标是<tt>(1,0)</tt>，顶点2在<tt>(0,1)</tt>，顶点3在<tt>(1,1)</tt>。
 * GeometryInfo<dim>::unit_cell_vertex()
 * 函数可以用来在运行时查询这些信息。
 *
 *  <h3>Implementation conventions for three spatial dimensions</h3>
 * 按照惯例，我们将对三空间维度的六面体的顶点、线和面使用以下编号惯例。在给出这些约定之前，我们声明以下草图是绘制六面体的三维图片的标准方式。
 *
 * @verbatim
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 *                    /|       |       /       /|
 *                   / |       |      /       / |
 * z                 /  |       |     /       /  |
 * ^                  |       |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*   |
 * |   ^y           |
 *
 *
 *
 *
 *
 * -------*    |       |
 * |  /             |  /       /     |       |  /
 * | /              | /       /      |       | /
 * |/               |/       /       |       |/
 *
 *
 *
 *
 *
 * ------>x
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 * @endverbatim
 * 图片的左边部分显示的是立方体的左面、底部和背面，而右边的应是顶部、右侧和正面。你可以通过将这两部分移到一起恢复整个立方体。
 * 请再次注意，以下几个约定的信息可以在运行或编译时从本类的成员函数和变量中提取。
 * <h4>Vertices</h4>
 * 3d中顶点的排序是由与2D情况相同的规则定义的。特别是，以下情况仍然是真实的。
 *
 *
 *
 * - 顶点的编号是按词法排序的。
 * 因此，顶点的编号方式如下
 *
 * @verbatim
 *     6-------7        6-------7
 *    /|       |       /       /|
 *   / |       |      /       / |
 *  /  |       |     /       /  |
 * 4   |       |    4-------5   |
 * |   2-------3    |       |   3
 * |  /       /     |       |  /
 * | /       /      |       | /
 * |/       /       |       |/
 * 0-------1        0-------1
 * @endverbatim
 *
 * 我们注意到，首先，底面的顶点（z=0）的编号方式与四边形的顶点完全相同。然后将底面的顶点(z=1)移到顶面，进行类似的编号。同样，
 * GeometryInfo<dim>::unit_cell_vertex()
 * 函数可以用来在运行时查询这些信息。
 *
 *  <h4>Lines</h4> 这里，与顶点的情况相同。
 *
 *
 *
 * - 3d中的线序。  <ul>   <li>  首先是面(z=0)的线，按2D线排序， <li>  然后是面(z=1)的线，按2D线排序， <li>  最后是Z方向的线，按lexicographic排序 </ul>  。
 * @verbatim
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---7---*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---7---*
 *    /|       |       /       /|
 *   4 |       11     4       5 11
 *  /  10      |     /       /  |
 *   |       |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---6---*   |
 * |
 *
 *
 *
 *
 *
 * ---3---*    |       |
 * |  /       /     |       9  /
 * 8 0       1      8       | 1
 * |/       /       |       |/
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---2---*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---2---*
 * @endverbatim
 * 如同在2d中，线是以坐标方向指向的。
 *
 * @verbatim
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * --->---*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * --->---*
 *    /|       |       /       /|
 *   ^ |       ^      ^       ^ ^
 *  /  ^       |     /       /  |
 *   |       |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * --->---*   |
 * |
 *
 *
 *
 *
 *
 * --->---*    |       |
 * |  /       /     |       ^  /
 * ^ ^       ^      ^       | ^
 * |/       /       |       |/
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * --->---*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * --->---*
 * @endverbatim
 *
 * 边缘（就像顶点和面）是以其自身的权利存储的实体，而不是在每次需要时从单元格中构造出来，这一事实意味着相邻的单元格实际上有指向边缘的指针，因此它们之间是共享的。这意味着平行边的集合具有平行方向的约定不仅仅是一个局部条件。在单元格列表被传递给Triangulation类的对象以创建三角形之前，你必须确保单元格的方向是兼容的，这样边的方向才是全局的，符合上述惯例。然而，GridReordering类可以为你做到这一点，它可以对输入单元的任意列表中的单元和边进行重新定向，这些单元不需要已经被排序。
 * <h4>Faces</h4>
 * 3D中面的编号是由一个类似2D的规则来定义的。
 *
 *
 *
 * - 面（3d中的四边形）：首先是在x-方向有法线的两个面，然后是y-和z-方向。对于每两个面：首先是法线在负坐标方向的面，然后是法线在正方向的面，也就是说，面的顺序是根据它们的法线指向
 *
 * - , x,
 *
 * - , y,
 *
 *
 *
 * - ，Z方向。
 * 因此，面的编号顺序为：左、右、前、后、底、顶面。
 * @verbatim
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 *    /|       |       /       /|
 *   / |   3   |      /   5   / |
 *  /  |       |     /       /  |
 *   |       |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*   |
 * | 0-------*    |       | 1
 * |  /       /     |       |  /
 * | /   4   /      |   2   | /
 * |/       /       |       |/
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 * @endverbatim
 *
 * 面的 <em> 标准 </em> 方向是这样的，即诱导的2d局部坐标系（x,y）意味着（右手规则）面法线方向的法线，见N2a）。 在下文中，我们展示了局部坐标系和面线的编号。  <ul>   <li>  面的0和1。
 * @verbatim
 *        Face 0           Face 1
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 *     /|       |       /       /|
 *    3 1       |      /       3 1
 *  y/  |       |     /      y/  |
 *    |x      |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*   |x
 *  |
 *
 *
 *
 *
 *
 * -------*    |       |
 *  0  /       /     |       0  /
 *  | 2       /      |       | 2
 *  |/       /       |       |/
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 * @endverbatim
 * <li>  脸部2和3。
 * @verbatim
 *      x Face 3           Face 2
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---1---*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 *     /|       |       /       /|
 *    / |       3      /       / |
 *   /  2       |    x/       /  |
 *    |       |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---1---*   |
 *  |
 *
 *
 *
 *
 *
 * ---0---*y   |       |
 *  |  /       /     |       3  /
 *  | /       /      2       | /
 *  |/       /       |       |/
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---0---*y
 * @endverbatim
 *
 * <li>  脸部4和5。
 * @verbatim
 *        Face 4         y Face 5
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---3---*
 *     /|       |       /       /|
 *    / |       |      0       1 |
 *   /  |       |     /       /  |
 *    |y      |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---2---* x |
 *  |
 *
 *
 *
 *
 *
 * ---3---*    |       |
 *  |  /       /     |       |  /
 *  | 0       1      |       | /
 *  |/       /       |       |/
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---2---* x
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 * @endverbatim
 * </ul>
 * 面的行号（0,1,2,3）对应于以下细胞行号。  <ul>   <li>  面0：第8、10、0、4行；  <li>  面1：第9、11、1、5行；  <li>  面2：第2、6、8、9行；  <li>  ] 面部3：第3、7、10、11行；  <li>  面部4：第0、1、2、3行；  <li>  面部5：第4、5、6、7行；  </ul>  你可以用 GeometryInfo<3>::face_to_cell_lines() 函数得到这些数字。
 * 脸部法线可以通过应用右手规则（x,y
 *
 * -> 法线）。)
 * 我们注意到，在2D的标准方向上，面0和面2的法线指向单元格内，面1和面3的法线指向外。在3D中，面0、2和4的法线指向单元格内，而面1、3和5的法线指向外面。这些信息同样可以从
 * GeometryInfo<dim>::unit_normal_orientation. 中查询到。
 * 然而，事实证明，大量的三维网格不能满足这个约定。这是由于一个单元的面的约定已经暗示了相邻单元的东西，因为它们共享一个共同的面，对第一个单元的固定也固定了两个单元的相对面的法向量。很容易构建单元格循环的案例，对于这些案例，我们无法为所有面找到与该约定一致的方向。
 * 由于这个原因，上述惯例只是我们所说的 <em> 标准方向 </em>  ......II实际上允许3d中的面要么有标准方向，要么有其相反的方向，在这种情况下，构成单元格的线会有还原的顺序，上述线的等价关系就不再成立了。你可以通过调用<tt>cell->face_orientation(face_no)</tt>来询问一个单元格的某个面是否有标准方向：如果结果是 @p true, ，那么这个面就有标准方向，否则它的法向量就指向另一个方向。在应用程序中，你需要这些信息的地方其实并不多，但库中有几个地方用到了这个信息。注意，在2d中，结果总是 @p true. 关于这个主题的更多信息可以在这个 @ref GlossFaceOrientation "词汇表 "
 * 文章中找到。
 * 为了允许3D中的各种网格，包括 <em>  Moebius </em>  -loops，一个面甚至可能从一个单元看是旋转的，而从共享该特定面的邻近单元看则是按照标准。为了解决这个问题，有两个标志<tt>face_flip</tt>和<tt>face_rotation</tt>，分别代表180度和90度的旋转。设置这两个标志相当于270度的旋转（都是逆时针）。你可以像询问<tt>面的方向</tt>一样询问单元格的这些标志。为了启用旋转的面，甚至线条也可以偏离它们在3D中的标准方向。这个信息可以作为<tt>line_orientation</tt>标志提供给三维中的单元和面。同样，这应该是库内部的东西，应用程序可能永远都不需要去管它。更多的信息请参见  @ref GlossFaceOrientation  "这个词汇表条目"
 * 。
 *
 *  <h4>Children</h4>
 * 一个各向同性的细化单元的八个子单元是根据它们相邻的顶点来编号的。
 *
 * @verbatim
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 *    /| 6  |  7 |       / 6  /  7 /|
 *  6|    |    |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*7|
 *  /|----*----*     / 4  /  5 /|
 * |/|    |    |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----* |/|
 * |4* | 2  |  3 |    | 4  |  5 |5*3|
 * |/|2*----*----*    |    |    |/|
 * |/ 2  /  3 /
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----* |/
 * |0*----*----*      |    |    |1*
 * |/0   /  1 /       | 0  |  1 |/
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 * @endverbatim
 * 考虑到面的方向，下列子节点与各自的面相邻。  <ul>   <li>  面0：子0、2、4、6；  <li>  面1：子1、3、5、7；  <li>  面2：子0、4、1、5；  <li>  面3：子2、6、3、7；  <li>  面4：子0、1、2、3；  <li>  面5：子4、5、6、7。  </ul>  你可以用 GeometryInfo<3>::child_cell_on_face() 函数得到这些数字。由于每个孩子都与具有相同数字的顶点相邻，这些数字也可以通过 GeometryInfo<3>::face_to_cell_vertices() 函数得到。
 * 请注意，上述列表只适用于标准方向的面。如果一个面不在标准方向上，那么在位置1和2（从0到3算起）的孩子会被交换。事实上，这就是GeometryInfo<3>的child_cell_on_face和face_to_cell_vertices函数在调用<tt>face_orientation=false</tt>参数时的作用。
 * 哪个子单元在哪个面的哪个位置的信息，在计算有悬挂节点的面之间的跳跃项时，最常使用的是FESubfaceValues类型的对象。坐在一个单元格上，你会看一个面，然后算出邻居的哪个子单元是坐在现在和邻居单元之间的一个给定的子表面上。为了避免在这种情况下每次都要查询两个单元格面的标准方向，你应该使用类似<tt>cell->neighbor_child_on_subface(face_no,subface_no)</tt>的函数调用，它在2D（面的方向不重要）和3D（需要使用面的方向作为
 * <tt>GeometryInfo<3>::child_cell_on_face</tt>).
 * 的附加参数）中都能返回正确的结果。
 * 对于各向异性的细化，子单元不能根据相邻的顶点进行编号，因此使用了以下惯例。
 *
 * @verbatim
 *          RefinementCase<3>::cut_x
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 *    /|    |    |       /    /    /|
 *   / |    |    |      / 0  /  1 / |
 *  /  | 0  |  1 |     /    /    /  |
 *   |    |    |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*   |
 * | 0 |    |    |    |    |    | 1 |
 * |
 *
 *
 *
 *
 *
 * ----*----*    |    |    |
 * |  /    /    /     | 0  | 1  |  /
 * | / 0  /  1 /      |    |    | /
 * |/    /    /       |    |    |/
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 * @endverbatim
 *
 *
 *
 * @verbatim
 *          RefinementCase<3>::cut_y
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------*
 *    /|         |       /    1    /|
 *   |         |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------* |
 *  /| |    1    |     /    0    /| |
 * |1|         |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------* |1|
 * | | |         |    |         | | |
 * |0|---------*    |         |0|
 * | |/    1    /     |    0    | |/
 * |---------*      |         |
 * |/    0    /       |         |/
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------*
 * @endverbatim
 *
 * @verbatim
 *          RefinementCase<3>::cut_z
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------*
 *    /|    1    |       /         /|
 *   / |         |      /    1    / |
 *  /
 *
 * ---------*     /         /
 * 1/|         |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------* 1/|
 * | / |    0    |    |    1    | / |
 * |/
 *
 * ---------*    |         |/
 * 0/         /
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------* 0/
 * | /    0    /      |         | /
 * |/         /       |    0    |/
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------*
 * @endverbatim
 * @verbatim
 *          RefinementCase<3>::cut_xy
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 *    /|    |    |       / 2  /  3 /|
 *   |    |    |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----* |
 *  /| | 2  |  3 |     / 0  /  1 /| |
 * |2|    |    |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----* |3|
 * | | |    |    |    |    |    | | |
 * |0|----*----*    |    |    |1|
 * | |/ 2  /  3 /     | 0  |  1 | |/
 * |----*----*      |    |    |
 * |/ 0  /  1 /       |    |    |/
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 * @endverbatim
 *
 * @verbatim
 *          RefinementCase<3>::cut_xz
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 *    /| 1  |  3 |       /    /    /|
 *   / |    |    |      / 1  /  3 / |
 *  /
 *
 * ----*----*     /    /    /
 * 1/|    |    |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----* 3/|
 * | / | 0  |  2 |    | 1  |  3 | / |
 * |/
 *
 * ----*----*    |    |    |/
 * 0/    /    /
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----* 2/
 * | / 0  /  2 /      |    |    | /
 * |/    /    /       | 0  |  2 |/
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ----*----*
 * @endverbatim
 * @verbatim
 *          RefinementCase<3>::cut_yz
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------*
 *    /|    3    |       /    3    /|
 *   |         |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------* |
 *  /|3*---------*     /    2    /|3*
 * |/|         |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------* |/|
 * |2* |    1    |    |    2    |2* |
 * |/|1*---------*    |         |/|1*
 * |/    1    /
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------* |/
 * |0*---------*      |         |0*
 * |/    0    /       |    0    |/
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------*
 * @endverbatim
 *
 * 这一信息也可以通过 <tt>GeometryInfo<3>::child_cell_on_face</tt>
 * 函数获得。 <h4>Coordinate systems</h4>
 * 我们为单元格顶点的显性坐标定义了以下坐标系。
 *
 * @verbatim
 *                     6-------7        6-------7
 *                    /|       |       /       /|
 *                   / |       |      /       / |
 * z                 /  |       |     /       /  |
 * ^                4   |       |    4-------5   |
 * |   ^y           |   2-------3    |       |   3
 * |  /             |  /       /     |       |  /
 * | /              | /       /      |       | /
 * |/               |/       /       |       |/
 *
 *
 *
 *
 *
 * ------>x        0-------1        0-------1
 * @endverbatim
 *
 * 根据上面的惯例，顶点的坐标如下（排位法，x跑得最快）。  <ul>   <li>  顶点0：<tt>(0,0,0)</tt>;  <li>  顶点1：<tt>(1,0,0)</tt>;  <li>  顶点2：<tt>（0,1,0）</tt>;  <li>  顶点3：<tt>（1,1,0）</tt>;  <li>  顶点4。<tt>(0,0,1)</tt>;  <li>  顶点5：<tt>(1,0,1)</tt>;  <li>  顶点6。<tt>(0,1,1)</tt>;  <li>  顶点7：<tt>(1,1,1)</tt>。  </ul>
 *
 *
 *
 * @note
 * 这个模板的实例是为维度1,2,3,4提供的，还有一个针对dim=0的特殊化（见手册中的
 * @ref Instantiations 部分）。
 *
 *
 * @ingroup grid geomprimitives aniso
 *
 */
template <int dim>
struct GeometryInfo
{
  /**
   * 精炼单元的最大子女数，即各向同性精炼单元的子女数。
   * 如果一个单元被各向异性地精炼，实际的子女数可能少于这里给出的值。
   *
   */
  static constexpr unsigned int max_children_per_cell = 1 << dim;

  /**
   * 一个单元的面的数量。
   *
   */
  static constexpr unsigned int faces_per_cell = 2 * dim;

  /**
   * 返回一个对象，它可以被认为是一个包含从零到`面孔_per_cell`所有索引的数组。这允许使用基于范围的for循环来编写以下类型的代码。
   * @code
   * for (auto &cell : triangulation.active_cell_iterators())
   *   for (auto face_index : GeometryInfo<dim>::face_indices())
   *     if (cell->face(face_index)->at_boundary())
   *       ... do something ...
   * @endcode
   * 这里，我们正在循环所有单元格的所有面，`face_index`接收所有有效的面的索引（1d中的0和1，2D中的0到3，以及3D中的0到5）。
   * @see  CPP11
   *
   */
  static std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  face_indices();

  /**
   * 一个精炼面的最大子女数，即一个各向同性的精炼面的子女数。
   * 如果一个单元被各向异性地细化，实际的子女数可能少于这里给出的值。
   *
   */
  static constexpr unsigned int max_children_per_face =
    GeometryInfo<dim - 1>::max_children_per_cell;

  /**
   * 一个单元格的顶点数。
   *
   */
  static constexpr unsigned int vertices_per_cell = 1 << dim;

  /**
   * 返回一个对象，它可以被认为是一个包含从零到`顶点_per_cell'所有索引的数组。这样就可以使用基于范围的for循环来编写以下类型的代码。
   * @code
   * for (auto &cell : triangulation.active_cell_iterators())
   *   for (auto vertex_index : GeometryInfo<dim>::vertex_indices())
   *     if (cell->vertex(vertex_index) satisfies some condition)
   *       ... do something ...
   * @endcode
   * 这里，我们在所有单元格的所有顶点上循环，`vertex_index`使用所有有效的索引。
   * @see  CPP11
   *
   */
  static std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  vertex_indices();

  /**
   * 每个面上的顶点数量。
   *
   */
  static constexpr unsigned int vertices_per_face =
    GeometryInfo<dim - 1>::vertices_per_cell;

  /**
   * 每个面的线数。
   *
   */
  static constexpr unsigned int lines_per_face =
    GeometryInfo<dim - 1>::lines_per_cell;

  /**
   * 每个面上的四边形数目。
   *
   */
  static constexpr unsigned int quads_per_face =
    GeometryInfo<dim - 1>::quads_per_cell;

  /**
   * 一个单元格的行数。
   * 计算公式利用了这样一个事实：当从一个维度到下一个维度时，低维度的对象被复制一次（因此是旧的线数的两倍），然后在旧对象的每个顶点和副本中的相应顶点之间插入一条新线。
   *
   */
  static constexpr unsigned int lines_per_cell =
    (2 * GeometryInfo<dim - 1>::lines_per_cell +
     GeometryInfo<dim - 1>::vertices_per_cell);

  /**
   * 一个单元格的四边形的数量。
   * 这个数字与前一个数字一样是递归计算的，不同的是，新的四边形是由连接原始线和它的副本产生的。
   *
   */
  static constexpr unsigned int quads_per_cell =
    (2 * GeometryInfo<dim - 1>::quads_per_cell +
     GeometryInfo<dim - 1>::lines_per_cell);

  /**
   * 一个单元格的六面体的数量。
   *
   */
  static constexpr unsigned int hexes_per_cell =
    (2 * GeometryInfo<dim - 1>::hexes_per_cell +
     GeometryInfo<dim - 1>::quads_per_cell);

  /**
   * 为UCD输出重新排列顶点。
   * 对于正在以UCD格式写入的单元格，该字段的每个条目包含<tt>deal.II</tt>中与该位置的UCD编号相对应的一个顶点的编号。
   * 典型的例子：写一个单元并安排顶点，这样UCD就能理解它们。
   * @code
   * for (i=0; i< n_vertices; ++i)
   * out << cell->vertex(ucd_to_deal[i]);
   * @endcode
   * 由于deal.II版本<=5.1中的顶点编号恰好与UCD的编号一致，这个字段也可以像old_to_lexicographic映射一样使用。
   *
   */
  static constexpr std::array<unsigned int, vertices_per_cell> ucd_to_deal =
    internal::GeometryInfoHelper::Initializers<dim>::ucd_to_deal();

  /**
   * 为OpenDX输出重新排列顶点。
   * 对于一个被写成OpenDX格式的单元格，这个字段的每个条目都包含<tt>deal.II</tt>中的一个顶点的编号，它与这个位置的DX编号相对应。
   * 典型的例子：写一个单元并安排顶点，这样OpenDX就能理解它们。
   * @code
   * for (i=0; i< n_vertices; ++i)
   * out << cell->vertex(dx_to_deal[i]);
   * @endcode
   *
   *
   */
  static constexpr std::array<unsigned int, vertices_per_cell> dx_to_deal =
    internal::GeometryInfoHelper::Initializers<dim>::dx_to_deal();

  /**
   * 这个字段为每个顶点存储它所属的面。在任何给定的维度中，面的数量都等于维度。这个二维数组的第一个索引是所有顶点的，第二个索引是顶点所属的
   * @p dim 面。
   * 每个顶点的面的顺序是这样的：第一个列出的面在<i>x</i>方向上包围参考单元，第二个在<i>y</i>方向上，以此类推。
   *
   */
  static constexpr ndarray<unsigned int, vertices_per_cell, dim>
    vertex_to_face =
      internal::GeometryInfoHelper::Initializers<dim>::vertex_to_face();

  /**
   * 返回用<tt>ref_case</tt>提炼的单元格（或面）的子节点数量。
   *
   */
  static unsigned int
  n_children(const RefinementCase<dim> &refinement_case);

  /**
   * 返回根据 internal::SubfaceCase  @p face_ref_case.
   * 精炼的一个面的子面的数量。
   *
   */
  static unsigned int
  n_subfaces(const internal::SubfaceCase<dim> &subface_case);

  /**
   * 给定参考元素上的一个具有
   * <code>internal::SubfaceCase@<dim@></code>  @p face_refinement_case
   * 的面，该函数返回 @p subface_no
   * 第1个子面的面积与该面的面积（=1）之间的比率。
   * 例如，对于 internal::SubfaceCase::cut_xy
   * ，每个子面的比例是1/4。
   *
   */
  static double
  subface_ratio(const internal::SubfaceCase<dim> &subface_case,
                const unsigned int                subface_no);

  /**
   * 给定一个用 <code>RefinementCase</code>   @p
   * cell_refinement_case精炼的单元格，返回 @p  face_no th面的
   * <code>SubfaceCase</code> 。
   *
   */
  static RefinementCase<dim - 1>
  face_refinement_case(const RefinementCase<dim> &cell_refinement_case,
                       const unsigned int         face_no,
                       const bool                 face_orientation = true,
                       const bool                 face_flip        = false,
                       const bool                 face_rotation    = false);

  /**
   * 给出第 @p face_no 个面的SubfaceCase @p face_refinement_case
   * ，返回最小的单元格的RefinementCase，它对应于该面的细化。
   *
   */
  static RefinementCase<dim>
  min_cell_refinement_case_for_face_refinement(
    const RefinementCase<dim - 1> &face_refinement_case,
    const unsigned int             face_no,
    const bool                     face_orientation = true,
    const bool                     face_flip        = false,
    const bool                     face_rotation    = false);

  /**
   * 给定一个用RefinementCase @p cell_refinement_case
   * 精炼的单元格，返回第 @p line_no 个面的RefinementCase。
   *
   */
  static RefinementCase<1>
  line_refinement_case(const RefinementCase<dim> &cell_refinement_case,
                       const unsigned int         line_no);

  /**
   * 返回单元格的最小/最小的RefinementCase，确保细化线 @p
   * line_no. 。
   *
   */
  static RefinementCase<dim>
  min_cell_refinement_case_for_line_refinement(const unsigned int line_no);

  /**
   * 这个字段存储了哪些子单元与母单元的某个面相邻。
   * 例如，在2D中，一个单元格的布局如下。
   * @verbatim
   * .      3
   * .   2-->--3
   * .   |     |
   * . 0 ^     ^ 1
   * .   |     |
   * .   0-->--1
   * .      2
   * @endverbatim
   * 顶点和面用它们的数字表示，面也用它们的方向表示。
   * 现在，在细化时，布局是这样的。
   * @verbatim
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * --*--*
   * | 2|3 |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * --*--*
   * | 0|1 |
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * --*--*
   * @endverbatim
   * 因此，面0上的子单元是（按面的方向排序）0和2，面3上是2和3，等等。    对于三个空间维度，子单元的确切顺序在这个类的一般文件中规定了。    通过<tt>face_orientation</tt>, <tt>face_flip</tt> 和 <tt>face_rotation</tt> 参数，这个函数可以处理以标准和非标准方向的面。<tt>face_orientation</tt>默认为<tt>true</tt>，<tt>face_flip</tt>和<tt>face_rotation</tt>默认为<tt>false</tt>（标准方向），在2d中没有影响。  面部方向的概念在这个 @ref GlossFaceOrientation "词汇表 "
   * 条目中得到了解释。
   * 在各向异性细化单元和面的情况下，面的 @p
   * RefinementCase，<tt>face_ref_case</tt>，可能对哪个子面在哪个给定子面后面有影响，因此这是一个附加参数，默认为面的各向同性细化。
   *
   */
  static unsigned int
  child_cell_on_face(const RefinementCase<dim> &    ref_case,
                     const unsigned int             face,
                     const unsigned int             subface,
                     const bool                     face_orientation = true,
                     const bool                     face_flip        = false,
                     const bool                     face_rotation    = false,
                     const RefinementCase<dim - 1> &face_refinement_case =
                       RefinementCase<dim - 1>::isotropic_refinement);

  /**
   * 将线的顶点编号映射到单元格顶点编号，即给出<tt>线</tt>的第1个顶点的单元格顶点编号，例如：
   * <tt>GeometryInfo<2>::line_to_cell_vertices(3,0)=2</tt>.
   * 线的顺序，以及它们的方向（这又决定了哪一个是线上的第一个顶点，哪一个是第二个顶点）是deal.II中的典范，如该类的一般文档中所述。
   * 对于<tt>dim=2</tt>，这个调用被简单地传递给face_to_cell_vertices()函数。
   *
   */
  static unsigned int
  line_to_cell_vertices(const unsigned int line, const unsigned int vertex);

  /**
   * 将面的顶点编号映射到单元格的顶点编号，即给出面的<tt>顶点</tt>第1个顶点的单元格顶点编号，例如， <tt>GeometryInfo<2>::face_to_cell_vertices(3,0)=2</tt>, 见本类文档2d部分N4点下的图片。    通过<tt>face_orientation</tt>, <tt>face_flip</tt>和<tt>face_rotation</tt>参数，这个函数可以处理以标准和非标准方向的面。<tt>face_orientation</tt>默认为<tt>true</tt>，<tt>face_flip</tt>和<tt>face_rotation</tt>默认为<tt>false</tt>（标准方向）。在2d中，只有<tt>face_flip</tt>被考虑。更多信息请参见这篇 @ref GlossFaceOrientation "词汇表 "
   * 文章。
   * 由于单元格的子节点是根据单元格的顶点排序的，这个调用被传递给child_cell_on_face()函数。
   * 因此，这个函数只是child_cell_on_face()的一个包装器，给它起了一个暗示性的名字。
   *
   */
  static unsigned int
  face_to_cell_vertices(const unsigned int face,
                        const unsigned int vertex,
                        const bool         face_orientation = true,
                        const bool         face_flip        = false,
                        const bool         face_rotation    = false);

  /**
   * 将面的行数映射到单元格行数，即给出面的<tt>行</tt>第1行的单元格行数，例如：
   * <tt>GeometryInfo<3>::face_to_cell_lines(5,0)=4</tt>.
   * 通过<tt>face_orientation</tt>,
   * <tt>face_flip</tt>和<tt>face_rotation</tt>参数，这个函数处理以标准和非标准方向为导向的面。<tt>face_orientation</tt>默认为<tt>true</tt>，<tt>face_flip</tt>和<tt>face_rotation</tt>默认为<tt>false</tt>（标准方向），在2d中没有影响。
   *
   */
  static unsigned int
  face_to_cell_lines(const unsigned int face,
                     const unsigned int line,
                     const bool         face_orientation = true,
                     const bool         face_flip        = false,
                     const bool         face_rotation    = false);

  /**
   * 将标准方向的面的顶点索引 @p vertex 映射到具有任意 @p
   * face_orientation,  @p face_flip 和 @p
   * 面旋转的面。这三个标志的值分别默认为<tt>true</tt>,
   * <tt>false</tt>和<tt>false</tt>。这个组合描述了一个标准方向的面。
   * 这个函数只在三维中实现。
   *
   */
  static unsigned int
  standard_to_real_face_vertex(const unsigned int vertex,
                               const bool         face_orientation = true,
                               const bool         face_flip        = false,
                               const bool         face_rotation    = false);

  /**
   * 将具有任意 @p  face_orientation、 @p face_flip 和 @p face_rotation
   * 的面的顶点索引 @p vertex
   * 映射到一个标准方向的面。这三个标志的值分别默认为<tt>true</tt>,
   * <tt>false</tt>和<tt>false</tt>。这个组合描述了一个标准方向的面。
   * 这个函数只在三维中实现。
   *
   */
  static unsigned int
  real_to_standard_face_vertex(const unsigned int vertex,
                               const bool         face_orientation = true,
                               const bool         face_flip        = false,
                               const bool         face_rotation    = false);

  /**
   * 将一个标准方向的面的线条索引 @p line
   * 映射到一个具有任意 @p face_orientation,  @p face_flip 和 @p
   * face_rotation的面。这三个标志的值默认分别为<tt>true</tt>,
   * <tt>false</tt>和<tt>false</tt>。这个组合描述了一个标准方向的面。
   * 这个函数只在三维中实现。
   *
   */
  static unsigned int
  standard_to_real_face_line(const unsigned int line,
                             const bool         face_orientation = true,
                             const bool         face_flip        = false,
                             const bool         face_rotation    = false);

  /**
   * 将标准方向的线的顶点索引 @p vertex 映射到一个任意 @p
   * line_orientation.
   * 的面的顶点索引，这个标志的值默认为<tt>true</tt>。
   *
   */
  static unsigned int
  standard_to_real_line_vertex(const unsigned int vertex,
                               const bool         line_orientation = true);

  /**
   * 将四边形中的顶点索引分解为一对线的索引和该线中的顶点索引。
   * @note 哪条线被选中并不重要（也不是故意暴露）。
   *
   */
  static std::array<unsigned int, 2>
  standard_quad_vertex_to_line_vertex_index(const unsigned int vertex);

  /**
   * 将一个六边形中的顶点索引分解为一对四边形索引和这个四边形中的顶点索引。
   * @note 选择哪个象限并不重要（也不是故意暴露的）。
   *
   */
  static std::array<unsigned int, 2>
  standard_hex_vertex_to_quad_vertex_index(const unsigned int vertex);

  /**
   * 将六边形中的线索引分解成一对四边形索引和该四边形中的线索引。
   * @note 选择哪个象限并不重要（也不是故意暴露的）。
   *
   */
  static std::array<unsigned int, 2>
  standard_hex_line_to_quad_line_index(const unsigned int line);

  /**
   * 将具有任意 @p face_orientation,  @p face_flip 和 @p face_rotation
   * 的面的线指数 @p line
   * 映射到一个标准方向的面。这三个标志的值默认分别为<tt>true</tt>,
   * <tt>false</tt>和<tt>false</tt>。这个组合描述了一个标准方向的面。
   * 这个函数只在三维中实现。
   *
   */
  static unsigned int
  real_to_standard_face_line(const unsigned int line,
                             const bool         face_orientation = true,
                             const bool         face_flip        = false,
                             const bool         face_rotation    = false);

  /**
   * 返回单元格上 @p ith
   * 顶点的位置。顶点的顺序是deal.II中的典范顺序，如该类的一般文档中所述。
   *
   */
  static Point<dim>
  unit_cell_vertex(const unsigned int vertex);

  /**
   * 给出一个单位坐标的点 @p p
   * ，返回它所在的子单元的编号。如果该点位于两个子单元的界面上，则返回它们中的任何一个指数。结果总是小于
   * GeometryInfo<dimension>::max_children_per_cell.
   * ，子单元的顺序在该类的一般文档中描述。
   *
   */
  static unsigned int
  child_cell_from_point(const Point<dim> &p);

  /**
   * 给出单元格上的坐标 @p p
   * ，返回该点在给定子单元格的坐标系中的坐标值。
   * 原始坐标和返回的坐标实际上都不需要在单元格内，我们只是简单地进行了一个缩放和移位操作，移位的程度取决于子单元的数量。
   *
   */
  static Point<dim>
  cell_to_child_coordinates(const Point<dim> &        p,
                            const unsigned int        child_index,
                            const RefinementCase<dim> refine_case =
                              RefinementCase<dim>::isotropic_refinement);

  /**
   * 与上面的函数相反：在子单元的坐标系中取一个点，并将其转换为母单元的坐标系。
   *
   */
  static Point<dim>
  child_to_cell_coordinates(const Point<dim> &        p,
                            const unsigned int        child_index,
                            const RefinementCase<dim> refine_case =
                              RefinementCase<dim>::isotropic_refinement);

  /**
   * 如果给定的点在当前空间维度的单元格内，则返回真。
   *
   */
  static bool
  is_inside_unit_cell(const Point<dim> &p);

  /**
   * 如果给定的点在当前空间维度的单元格内，则返回真。这个函数接受一个额外的参数，指定点的位置实际上可能在真正的单元格之外的程度。这是很有用的，因为在实践中，我们往往不能准确地计算出参考坐标中的点的坐标，而只能是数字上的舍入。
   * 容差参数可以小于零，表示该点应该安全地在单元格内。
   *
   */
  static bool
  is_inside_unit_cell(const Point<dim> &p, const double eps);

  /**
   * 将一个给定的点投射到单元格上，也就是说，[0...1]以外的每个坐标被修改为位于该区间内。
   *
   */
  template <typename Number = double>
  static Point<dim, Number>
  project_to_unit_cell(const Point<dim, Number> &p);

  /**
   * 返回单元格外的给定点 @p p
   * 到最近的单元格边界之间的向量的无穷大规范。对于单元格内的点，这被定义为零。
   *
   */
  static double
  distance_to_unit_cell(const Point<dim> &p);

  /**
   * 计算 $i$ -第 $d$
   * -线性（即（双，三）线性）形状函数在位置 $\xi$
   * 的值。
   *
   */
  static double
  d_linear_shape_function(const Point<dim> &xi, const unsigned int i);

  /**
   * 计算 $i$ -第 $d$
   * -线性（即（双，三）线性）形状函数在位置 $\xi$
   * 的梯度。
   *
   */
  static Tensor<1, dim>
  d_linear_shape_function_gradient(const Point<dim> &xi, const unsigned int i);

  /**
   * 对于从参考单元、面或边到由给定顶点指定的物体的（双、三）线性映射，计算转换后的单位向量顶点的交替形式。对于维数为
   * @p dim, 的物体，有 @p dim 个矢量，每个矢量有 @p spacedim
   * 个分量，交替形式是等级为spacedim-dim的张量，对应于 @p
   * dim
   * 个单位矢量的楔形积，它对应于从参考元素到顶点描述的元素映射的体积和法向量。
   * 例如，如果dim==spacedim==2，那么交替形式是一个标量（因为spacedim-dim=0），其值等于
   * $\mathbf v_1\wedge \mathbf v_2=\mathbf v_1^\perp \cdot\mathbf v_2$
   * ，其中 $\mathbf v_1^\perp$ 是一个从 $\mathbf v_1$
   * 向右旋转90度的矢量。如果dim==spacedim==3，那么结果又是一个标量，其值为
   * $\mathbf v_1\wedge \mathbf v_2 \wedge \mathbf v_3 = (\mathbf v_1\times
   * \mathbf v_2)\cdot \mathbf v_3$ ，其中 $\mathbf v_1, \mathbf v_2,
   * \mathbf v_3$
   * 是单位dim维单元的一个顶点的单位向量在转化为spacedim维空间的dim维单元时的图像。在这两种情况下，即对于dim==2或3，其结果恰好等于实空间中从参考单元到单元的映射的Jacobian的行列式。请注意，这是实际的行列式，而不是它的绝对值，因为它经常用于将积分从一个坐标系转换到另一个坐标系。特别是，如果顶点指定的对象是一个平行四边形（即参考单元的线性变换），那么计算值在所有顶点都是相同的，并且等于该单元的（有符号的）面积；同样，对于平行四边形，它是该单元的体积。
   * 同样，如果我们有dim==spacedim-1（例如，我们有一个三维空间的四边形，或者一条二维空间的线），那么交替乘积表示每个顶点的物体的法向量（即一个秩-1张量，因为spacedim-dim=1），其中法向量的大小表示从参考物体到顶点所给物体的转换的面积元素。特别是，如果从参考物体到这里考虑的物体的映射是线性的（不是双线或三线），那么返回的向量都是%平行的，垂直于顶点描述的映射物体，并且其大小等于映射物体的面积/体积。如果dim=1，spacedim=2，那么返回值为
   * $\mathbf v_1^\perp$  ，其中 $\mathbf v_1$
   * 是映射到顶点给出的2D中的直线的唯一单位向量的图像；如果dim=2，spacedim=3，那么返回值为
   * $\mathbf v_1 \wedge \mathbf v_2=\mathbf v_1 \times \mathbf v_2$ ，其中
   * $\mathbf
   * v_1,\mathbf v_2$ 是与映射到三维空间的四边形相切的两个三维向量。    这个函数用于确定一个单元的扭曲程度（见术语表中 @ref GlossDistorted "扭曲的单元 "
   * 条目）。
   *
   */
  template <int spacedim>
  static void
  alternating_form_at_vertices
#ifndef DEAL_II_CXX14_CONSTEXPR_BUG
    (const Point<spacedim> (&vertices)[vertices_per_cell],
     Tensor<spacedim - dim, spacedim> (&forms)[vertices_per_cell]);
#else
    (const Point<spacedim> *vertices, Tensor<spacedim - dim, spacedim> *forms);
#endif

  /**
   * 对于参考单元的每个面，这个字段存储其法向量指向的坐标方向。在<tt>dim</tt>维度中，这些是<tt>2*dim</tt>的第一个条目<tt>{0,0,1,1,2,2,3,3}</tt>。
   * 请注意，这只是坐标数。法向量的实际方向是通过将这个方向的单位向量与#unit_normal_orientation相乘而得到的。
   *
   */
  static constexpr std::array<unsigned int, faces_per_cell>
    unit_normal_direction =
      internal::GeometryInfoHelper::Initializers<dim>::unit_normal_direction();

  /**
   * 参考单元的一个面的单位法向量的方向。在<tt>dim</tt>维度中，这些是<tt>2*dim</tt><tt>{-1,1,-1,1,-1,-1,1}</tt>的第一个条目。    每个值都是<tt>1</tt>或<tt>-1</tt>，分别对应于指向正坐标或负坐标方向的法向量。    注意，这只是面的 <em> 标准方向 </em> 。至少在3D中，三角形中单元格的实际面也可以有相反的方向，这取决于人们可以从它所属的单元格中查询到的一个标志。欲了解更多信息，请参见 @ref GlossFaceOrientation "词汇表 "
   * 中关于面的方向的条目。
   *
   */
  static constexpr std::array<int, faces_per_cell> unit_normal_orientation =
    internal::GeometryInfoHelper::Initializers<dim>::unit_normal_orientation();

  /**
   * 参考单元的一个面的单位法向量（Point<dim>）。    注意，这只是面的 <em> 标准方向 </em> 。至少在3D中，三角形中单元格的实际面也可以有相反的方向，这取决于一个可以从它所属的单元格中查询的标志。欲了解更多信息，请参见 @ref GlossFaceOrientation "词汇表 "
   * 中关于面的方向的条目。
   *
   */
  static constexpr std::array<Tensor<1, dim>, faces_per_cell>
    unit_normal_vector =
      internal::GeometryInfoHelper::Initializers<dim>::unit_normal_vector();

  /**
   * 参考单元面的单位切向矢量（Point<dim>的`dim-1`元素数组），排列在右手坐标系中，使两个矢量之间的交积返回单位法向量。    注意，这只是面的 <em> 标准方向 </em> 。至少在3D中，三角形中单元格的实际面也可以有相反的方向，这取决于一个可以从它所属的单元格中查询的标志。欲了解更多信息，请参见 @ref GlossFaceOrientation "词汇表 "
   * 中关于面的方向的条目。
   *
   */
  static constexpr ndarray<Tensor<1, dim>, faces_per_cell, dim - 1>
    unit_tangential_vectors = internal::GeometryInfoHelper::Initializers<
      dim>::unit_tangential_vectors();

  /**
   * 数字列表，表示哪个面与给定的面相对。它的条目是<tt>{
   * 1, 0, 3, 2, 5, 4, 7, 6}</tt>的第一个<tt>2*dim</tt>条目。
   *
   */
  static constexpr std::array<unsigned int, faces_per_cell> opposite_face =
    internal::GeometryInfoHelper::Initializers<dim>::opposite_face();


  /**
   * 异常情况
   *
   */
  DeclException1(ExcInvalidCoordinate,
                 double,
                 << "The coordinates must satisfy 0 <= x_i <= 1, "
                 << "but here we have x_i=" << arg1);

  /**
   * 异常情况
   *
   */
  DeclException3(ExcInvalidSubface,
                 int,
                 int,
                 int,
                 << "RefinementCase<dim> " << arg1 << ": face " << arg2
                 << " has no subface " << arg3);
};



#ifndef DOXYGEN


 /* -------------- declaration of explicit specializations ------------- */ 


template <>
Tensor<1, 1>
GeometryInfo<1>::d_linear_shape_function_gradient(const Point<1> &   xi,
                                                  const unsigned int i);
template <>
Tensor<1, 2>
GeometryInfo<2>::d_linear_shape_function_gradient(const Point<2> &   xi,
                                                  const unsigned int i);
template <>
Tensor<1, 3>
GeometryInfo<3>::d_linear_shape_function_gradient(const Point<3> &   xi,
                                                  const unsigned int i);



 /* -------------- inline functions ------------- */ 


inline GeometryPrimitive::GeometryPrimitive(const Object object)
  : object(object)
{}



inline GeometryPrimitive::GeometryPrimitive(const unsigned int object_dimension)
  : object(static_cast<Object>(object_dimension))
{}


inline GeometryPrimitive::operator unsigned int() const
{
  return static_cast<unsigned int>(object);
}



namespace internal
{
  template <int dim>
  inline SubfaceCase<dim>::SubfaceCase(
    const typename SubfacePossibilities<dim>::Possibilities subface_possibility)
    : value(subface_possibility)
  {}


  template <int dim>
  inline SubfaceCase<dim>::operator std::uint8_t() const
  {
    return value;
  }


} // namespace internal


template <int dim>
inline RefinementCase<dim>
RefinementCase<dim>::cut_axis(const unsigned int)
{
  Assert(false, ExcInternalError());
  return static_cast<std::uint8_t>(-1);
}


template <>
inline RefinementCase<1>
RefinementCase<1>::cut_axis(const unsigned int i)
{
  AssertIndexRange(i, 1);

  const RefinementCase options[1] = {RefinementPossibilities<1>::cut_x};
  return options[i];
}



template <>
inline RefinementCase<2>
RefinementCase<2>::cut_axis(const unsigned int i)
{
  AssertIndexRange(i, 2);

  const RefinementCase options[2] = {RefinementPossibilities<2>::cut_x,
                                     RefinementPossibilities<2>::cut_y};
  return options[i];
}



template <>
inline RefinementCase<3>
RefinementCase<3>::cut_axis(const unsigned int i)
{
  AssertIndexRange(i, 3);

  const RefinementCase options[3] = {RefinementPossibilities<3>::cut_x,
                                     RefinementPossibilities<3>::cut_y,
                                     RefinementPossibilities<3>::cut_z};
  return options[i];
}



template <int dim>
inline RefinementCase<dim>::RefinementCase()
  : value(RefinementPossibilities<dim>::no_refinement)
{}



template <int dim>
inline RefinementCase<dim>::RefinementCase(
  const typename RefinementPossibilities<dim>::Possibilities refinement_case)
  : value(refinement_case)
{
  // check that only those bits of
  // the given argument are set that
  // make sense for a given space
  // dimension
  Assert((refinement_case &
          RefinementPossibilities<dim>::isotropic_refinement) ==
           refinement_case,
         ExcInvalidRefinementCase(refinement_case));
}



template <int dim>
inline RefinementCase<dim>::RefinementCase(const std::uint8_t refinement_case)
  : value(refinement_case)
{
  // check that only those bits of
  // the given argument are set that
  // make sense for a given space
  // dimension
  Assert((refinement_case &
          RefinementPossibilities<dim>::isotropic_refinement) ==
           refinement_case,
         ExcInvalidRefinementCase(refinement_case));
}



template <int dim>
inline RefinementCase<dim>::operator std::uint8_t() const
{
  return value;
}



template <int dim>
inline RefinementCase<dim>
RefinementCase<dim>::operator|(const RefinementCase<dim> &r) const
{
  return RefinementCase<dim>(static_cast<std::uint8_t>(value | r.value));
}



template <int dim>
inline RefinementCase<dim> RefinementCase<dim>::
                           operator&(const RefinementCase<dim> &r) const
{
  return RefinementCase<dim>(static_cast<std::uint8_t>(value & r.value));
}



template <int dim>
inline RefinementCase<dim>
RefinementCase<dim>::operator~() const
{
  return RefinementCase<dim>(static_cast<std::uint8_t>(
    (~value) & RefinementPossibilities<dim>::isotropic_refinement));
}



template <int dim>
inline std::size_t
RefinementCase<dim>::memory_consumption()
{
  return sizeof(RefinementCase<dim>);
}



template <int dim>
template <class Archive>
inline void
RefinementCase<dim>::serialize(Archive &ar, const unsigned int)
{
  // serialization can't deal with bitfields, so copy from/to a full sized
  // std::uint8_t
  std::uint8_t uchar_value = value;
  ar &         uchar_value;
  value = uchar_value;
}



template <>
inline Point<1>
GeometryInfo<1>::unit_cell_vertex(const unsigned int vertex)
{
  AssertIndexRange(vertex, vertices_per_cell);

  return Point<1>(static_cast<double>(vertex));
}



template <>
inline Point<2>
GeometryInfo<2>::unit_cell_vertex(const unsigned int vertex)
{
  AssertIndexRange(vertex, vertices_per_cell);

  return {static_cast<double>(vertex % 2), static_cast<double>(vertex / 2)};
}



template <>
inline Point<3>
GeometryInfo<3>::unit_cell_vertex(const unsigned int vertex)
{
  AssertIndexRange(vertex, vertices_per_cell);

  return {static_cast<double>(vertex % 2),
          static_cast<double>(vertex / 2 % 2),
          static_cast<double>(vertex / 4)};
}



inline std::array<unsigned int, 0>
GeometryInfo<0>::face_indices()
{
  return {{}};
}



inline std::array<unsigned int, 1>
GeometryInfo<0>::vertex_indices()
{
  return {{0}};
}



template <int dim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
GeometryInfo<dim>::face_indices()
{
  return {0U, faces_per_cell};
}



template <int dim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
GeometryInfo<dim>::vertex_indices()
{
  return {0U, vertices_per_cell};
}



template <int dim>
inline Point<dim>
GeometryInfo<dim>::unit_cell_vertex(const unsigned int)
{
  Assert(false, ExcNotImplemented());

  return {};
}



template <>
inline unsigned int
GeometryInfo<1>::child_cell_from_point(const Point<1> &p)
{
  Assert((p[0] >= 0) && (p[0] <= 1), ExcInvalidCoordinate(p[0]));

  return (p[0] <= 0.5 ? 0 : 1);
}



template <>
inline unsigned int
GeometryInfo<2>::child_cell_from_point(const Point<2> &p)
{
  Assert((p[0] >= 0) && (p[0] <= 1), ExcInvalidCoordinate(p[0]));
  Assert((p[1] >= 0) && (p[1] <= 1), ExcInvalidCoordinate(p[1]));

  return (p[0] <= 0.5 ? (p[1] <= 0.5 ? 0 : 2) : (p[1] <= 0.5 ? 1 : 3));
}



template <>
inline unsigned int
GeometryInfo<3>::child_cell_from_point(const Point<3> &p)
{
  Assert((p[0] >= 0) && (p[0] <= 1), ExcInvalidCoordinate(p[0]));
  Assert((p[1] >= 0) && (p[1] <= 1), ExcInvalidCoordinate(p[1]));
  Assert((p[2] >= 0) && (p[2] <= 1), ExcInvalidCoordinate(p[2]));

  return (p[0] <= 0.5 ?
            (p[1] <= 0.5 ? (p[2] <= 0.5 ? 0 : 4) : (p[2] <= 0.5 ? 2 : 6)) :
            (p[1] <= 0.5 ? (p[2] <= 0.5 ? 1 : 5) : (p[2] <= 0.5 ? 3 : 7)));
}


template <int dim>
inline unsigned int
GeometryInfo<dim>::child_cell_from_point(const Point<dim> &)
{
  Assert(false, ExcNotImplemented());

  return 0;
}



template <>
inline Point<1>
GeometryInfo<1>::cell_to_child_coordinates(const Point<1> &        p,
                                           const unsigned int      child_index,
                                           const RefinementCase<1> refine_case)

{
  AssertIndexRange(child_index, 2);
  Assert(refine_case == RefinementCase<1>::cut_x, ExcInternalError());
  (void)refine_case; // removes -Wunused-parameter warning in optimized mode

  return Point<1>(p * 2.0 - unit_cell_vertex(child_index));
}



template <>
inline Point<2>
GeometryInfo<2>::cell_to_child_coordinates(const Point<2> &        p,
                                           const unsigned int      child_index,
                                           const RefinementCase<2> refine_case)

{
  AssertIndexRange(child_index, GeometryInfo<2>::n_children(refine_case));

  Point<2> point = p;
  switch (refine_case)
    {
      case RefinementCase<2>::cut_x:
        point[0] *= 2.0;
        if (child_index == 1)
          point[0] -= 1.0;
        break;
      case RefinementCase<2>::cut_y:
        point[1] *= 2.0;
        if (child_index == 1)
          point[1] -= 1.0;
        break;
      case RefinementCase<2>::cut_xy:
        point *= 2.0;
        point -= unit_cell_vertex(child_index);
        break;
      default:
        Assert(false, ExcInternalError());
    }

  return point;
}



template <>
inline Point<3>
GeometryInfo<3>::cell_to_child_coordinates(const Point<3> &        p,
                                           const unsigned int      child_index,
                                           const RefinementCase<3> refine_case)

{
  AssertIndexRange(child_index, GeometryInfo<3>::n_children(refine_case));

  Point<3> point = p;
  // there might be a cleverer way to do
  // this, but since this function is called
  // in very few places for initialization
  // purposes only, I don't care at the
  // moment
  switch (refine_case)
    {
      case RefinementCase<3>::cut_x:
        point[0] *= 2.0;
        if (child_index == 1)
          point[0] -= 1.0;
        break;
      case RefinementCase<3>::cut_y:
        point[1] *= 2.0;
        if (child_index == 1)
          point[1] -= 1.0;
        break;
      case RefinementCase<3>::cut_z:
        point[2] *= 2.0;
        if (child_index == 1)
          point[2] -= 1.0;
        break;
      case RefinementCase<3>::cut_xy:
        point[0] *= 2.0;
        point[1] *= 2.0;
        if (child_index % 2 == 1)
          point[0] -= 1.0;
        if (child_index / 2 == 1)
          point[1] -= 1.0;
        break;
      case RefinementCase<3>::cut_xz:
        // careful, this is slightly
        // different from xy and yz due to
        // different internal numbering of
        // children!
        point[0] *= 2.0;
        point[2] *= 2.0;
        if (child_index / 2 == 1)
          point[0] -= 1.0;
        if (child_index % 2 == 1)
          point[2] -= 1.0;
        break;
      case RefinementCase<3>::cut_yz:
        point[1] *= 2.0;
        point[2] *= 2.0;
        if (child_index % 2 == 1)
          point[1] -= 1.0;
        if (child_index / 2 == 1)
          point[2] -= 1.0;
        break;
      case RefinementCase<3>::cut_xyz:
        point *= 2.0;
        point -= unit_cell_vertex(child_index);
        break;
      default:
        Assert(false, ExcInternalError());
    }

  return point;
}



template <int dim>
inline Point<dim>
GeometryInfo<dim>::cell_to_child_coordinates(
  const Point<dim> &  /*p*/ ,
  const unsigned int  /*child_index*/ ,
  const RefinementCase<dim>  /*refine_case*/ )

{
  Assert(false, ExcNotImplemented());
  return {};
}



template <>
inline Point<1>
GeometryInfo<1>::child_to_cell_coordinates(const Point<1> &        p,
                                           const unsigned int      child_index,
                                           const RefinementCase<1> refine_case)

{
  AssertIndexRange(child_index, 2);
  Assert(refine_case == RefinementCase<1>::cut_x, ExcInternalError());
  (void)refine_case; // removes -Wunused-parameter warning in optimized mode

  return (p + unit_cell_vertex(child_index)) * 0.5;
}



template <>
inline Point<3>
GeometryInfo<3>::child_to_cell_coordinates(const Point<3> &        p,
                                           const unsigned int      child_index,
                                           const RefinementCase<3> refine_case)

{
  AssertIndexRange(child_index, GeometryInfo<3>::n_children(refine_case));

  Point<3> point = p;
  // there might be a cleverer way to do
  // this, but since this function is called
  // in very few places for initialization
  // purposes only, I don't care at the
  // moment
  switch (refine_case)
    {
      case RefinementCase<3>::cut_x:
        if (child_index == 1)
          point[0] += 1.0;
        point[0] *= 0.5;
        break;
      case RefinementCase<3>::cut_y:
        if (child_index == 1)
          point[1] += 1.0;
        point[1] *= 0.5;
        break;
      case RefinementCase<3>::cut_z:
        if (child_index == 1)
          point[2] += 1.0;
        point[2] *= 0.5;
        break;
      case RefinementCase<3>::cut_xy:
        if (child_index % 2 == 1)
          point[0] += 1.0;
        if (child_index / 2 == 1)
          point[1] += 1.0;
        point[0] *= 0.5;
        point[1] *= 0.5;
        break;
      case RefinementCase<3>::cut_xz:
        // careful, this is slightly
        // different from xy and yz due to
        // different internal numbering of
        // children!
        if (child_index / 2 == 1)
          point[0] += 1.0;
        if (child_index % 2 == 1)
          point[2] += 1.0;
        point[0] *= 0.5;
        point[2] *= 0.5;
        break;
      case RefinementCase<3>::cut_yz:
        if (child_index % 2 == 1)
          point[1] += 1.0;
        if (child_index / 2 == 1)
          point[2] += 1.0;
        point[1] *= 0.5;
        point[2] *= 0.5;
        break;
      case RefinementCase<3>::cut_xyz:
        point += unit_cell_vertex(child_index);
        point *= 0.5;
        break;
      default:
        Assert(false, ExcInternalError());
    }

  return point;
}



template <>
inline Point<2>
GeometryInfo<2>::child_to_cell_coordinates(const Point<2> &        p,
                                           const unsigned int      child_index,
                                           const RefinementCase<2> refine_case)
{
  AssertIndexRange(child_index, GeometryInfo<2>::n_children(refine_case));

  Point<2> point = p;
  switch (refine_case)
    {
      case RefinementCase<2>::cut_x:
        if (child_index == 1)
          point[0] += 1.0;
        point[0] *= 0.5;
        break;
      case RefinementCase<2>::cut_y:
        if (child_index == 1)
          point[1] += 1.0;
        point[1] *= 0.5;
        break;
      case RefinementCase<2>::cut_xy:
        point += unit_cell_vertex(child_index);
        point *= 0.5;
        break;
      default:
        Assert(false, ExcInternalError());
    }

  return point;
}



template <int dim>
inline Point<dim>
GeometryInfo<dim>::child_to_cell_coordinates(
  const Point<dim> &  /*p*/ ,
  const unsigned int  /*child_index*/ ,
  const RefinementCase<dim>  /*refine_case*/ )
{
  Assert(false, ExcNotImplemented());
  return {};
}



template <int dim>
inline bool
GeometryInfo<dim>::is_inside_unit_cell(const Point<dim> &)
{
  Assert(false, ExcNotImplemented());
  return false;
}

template <>
inline bool
GeometryInfo<1>::is_inside_unit_cell(const Point<1> &p)
{
  return (p[0] >= 0.) && (p[0] <= 1.);
}



template <>
inline bool
GeometryInfo<2>::is_inside_unit_cell(const Point<2> &p)
{
  return (p[0] >= 0.) && (p[0] <= 1.) && (p[1] >= 0.) && (p[1] <= 1.);
}



template <>
inline bool
GeometryInfo<3>::is_inside_unit_cell(const Point<3> &p)
{
  return (p[0] >= 0.) && (p[0] <= 1.) && (p[1] >= 0.) && (p[1] <= 1.) &&
         (p[2] >= 0.) && (p[2] <= 1.);
}



template <int dim>
inline bool
GeometryInfo<dim>::is_inside_unit_cell(const Point<dim> &, const double)
{
  Assert(false, ExcNotImplemented());
  return false;
}

template <>
inline bool
GeometryInfo<1>::is_inside_unit_cell(const Point<1> &p, const double eps)
{
  return (p[0] >= -eps) && (p[0] <= 1. + eps);
}



template <>
inline bool
GeometryInfo<2>::is_inside_unit_cell(const Point<2> &p, const double eps)
{
  const double l = -eps, u = 1 + eps;
  return (p[0] >= l) && (p[0] <= u) && (p[1] >= l) && (p[1] <= u);
}



template <>
inline bool
GeometryInfo<3>::is_inside_unit_cell(const Point<3> &p, const double eps)
{
  const double l = -eps, u = 1.0 + eps;
  return (p[0] >= l) && (p[0] <= u) && (p[1] >= l) && (p[1] <= u) &&
         (p[2] >= l) && (p[2] <= u);
}



template <>
inline unsigned int
GeometryInfo<1>::line_to_cell_vertices(const unsigned int line,
                                       const unsigned int vertex)
{
  (void)line;
  AssertIndexRange(line, lines_per_cell);
  AssertIndexRange(vertex, 2);

  return vertex;
}


template <>
inline unsigned int
GeometryInfo<2>::line_to_cell_vertices(const unsigned int line,
                                       const unsigned int vertex)
{
  constexpr unsigned int cell_vertices[4][2] = {{0, 2}, {1, 3}, {0, 1}, {2, 3}};
  return cell_vertices[line][vertex];
}



template <>
inline unsigned int
GeometryInfo<3>::line_to_cell_vertices(const unsigned int line,
                                       const unsigned int vertex)
{
  AssertIndexRange(line, lines_per_cell);
  AssertIndexRange(vertex, 2);

  constexpr unsigned vertices[lines_per_cell][2] = {{0, 2}, // bottom face
                                                    {1, 3},
                                                    {0, 1},
                                                    {2, 3},
                                                    {4, 6}, // top face
                                                    {5, 7},
                                                    {4, 5},
                                                    {6, 7},
                                                    {0,
                                                     4}, // connects of bottom
                                                    {1, 5}, //   top face
                                                    {2, 6},
                                                    {3, 7}};
  return vertices[line][vertex];
}


template <>
inline unsigned int
GeometryInfo<4>::line_to_cell_vertices(const unsigned int, const unsigned int)
{
  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}

template <>
inline unsigned int
GeometryInfo<3>::standard_to_real_face_vertex(const unsigned int vertex,
                                              const bool face_orientation,
                                              const bool face_flip,
                                              const bool face_rotation)
{
  AssertIndexRange(vertex, GeometryInfo<3>::vertices_per_face);

  // set up a table to make sure that
  // we handle non-standard faces correctly
  //
  // so set up a table that for each vertex (of
  // a quad in standard position) describes
  // which vertex to take
  //
  // first index: four vertices 0...3
  //
  // second index: face_orientation; 0:
  // opposite normal, 1: standard
  //
  // third index: face_flip; 0: standard, 1:
  // face rotated by 180 degrees
  //
  // forth index: face_rotation: 0: standard,
  // 1: face rotated by 90 degrees

  constexpr unsigned int vertex_translation[4][2][2][2] = {
    {{{0, 2},   // vertex 0, face_orientation=false, face_flip=false,
                // face_rotation=false and true
      {3, 1}},  // vertex 0, face_orientation=false, face_flip=true,
                // face_rotation=false and true
     {{0, 2},   // vertex 0, face_orientation=true, face_flip=false,
                // face_rotation=false and true
      {3, 1}}}, // vertex 0, face_orientation=true, face_flip=true,
                // face_rotation=false and true

    {{{2, 3}, // vertex 1 ...
      {1, 0}},
     {{1, 0}, {2, 3}}},

    {{{1, 0}, // vertex 2 ...
      {2, 3}},
     {{2, 3}, {1, 0}}},

    {{{3, 1}, // vertex 3 ...
      {0, 2}},
     {{3, 1}, {0, 2}}}};

  return vertex_translation[vertex][face_orientation][face_flip][face_rotation];
}



template <int dim>
inline unsigned int
GeometryInfo<dim>::standard_to_real_face_vertex(const unsigned int vertex,
                                                const bool,
                                                const bool,
                                                const bool)
{
  Assert(dim > 1, ExcImpossibleInDim(dim));
  AssertIndexRange(vertex, GeometryInfo<dim>::vertices_per_face);
  return vertex;
}

template <int dim>
inline unsigned int
GeometryInfo<dim>::n_children(const RefinementCase<dim> &ref_case)
{
  constexpr unsigned int n_children[RefinementCase<3>::cut_xyz + 1] = {
    0, 2, 2, 4, 2, 4, 4, 8};

  return n_children[ref_case];
}



template <int dim>
inline unsigned int
GeometryInfo<dim>::n_subfaces(const internal::SubfaceCase<dim> &)
{
  Assert(false, ExcNotImplemented());
  return 0;
}

template <>
inline unsigned int
GeometryInfo<1>::n_subfaces(const internal::SubfaceCase<1> &)
{
  Assert(false, ExcImpossibleInDim(1));
  return 0;
}

template <>
inline unsigned int
GeometryInfo<2>::n_subfaces(const internal::SubfaceCase<2> &subface_case)
{
  return (subface_case == internal::SubfaceCase<2>::case_x) ? 2 : 0;
}



template <>
inline unsigned int
GeometryInfo<3>::n_subfaces(const internal::SubfaceCase<3> &subface_case)
{
  const unsigned int nsubs[internal::SubfaceCase<3>::case_isotropic + 1] = {
    0, 2, 3, 3, 4, 2, 3, 3, 4, 4};
  return nsubs[subface_case];
}



template <int dim>
inline double
GeometryInfo<dim>::subface_ratio(const internal::SubfaceCase<dim> &,
                                 const unsigned int)
{
  Assert(false, ExcNotImplemented());
  return 0.;
}

template <>
inline double
GeometryInfo<1>::subface_ratio(const internal::SubfaceCase<1> &,
                               const unsigned int)
{
  return 1;
}


template <>
inline double
GeometryInfo<2>::subface_ratio(const internal::SubfaceCase<2> &subface_case,
                               const unsigned int)
{
  double ratio = 1;
  switch (subface_case)
    {
      case internal::SubfaceCase<2>::case_none:
        // Here, an
        // Assert(false,ExcInternalError())
        // would be the right
        // choice, but
        // unfortunately the
        // current function is
        // also called for faces
        // without children (see
        // tests/fe/mapping.cc).
        //          Assert(false, ExcMessage("Face has no subfaces."));
        // Furthermore, assign
        // following value as
        // otherwise the
        // bits/volume_x tests
        // break
        ratio = 1. / GeometryInfo<2>::max_children_per_face;
        break;
      case internal::SubfaceCase<2>::case_x:
        ratio = 0.5;
        break;
      default:
        // there should be no
        // cases left
        Assert(false, ExcInternalError());
        break;
    }

  return ratio;
}


template <>
inline double
GeometryInfo<3>::subface_ratio(const internal::SubfaceCase<3> &subface_case,
                               const unsigned int              subface_no)
{
  double ratio = 1;
  switch (subface_case)
    {
      case internal::SubfaceCase<3>::case_none:
        // Here, an
        // Assert(false,ExcInternalError())
        // would be the right
        // choice, but
        // unfortunately the
        // current function is
        // also called for faces
        // without children (see
        // tests/bits/mesh_3d_16.cc). Add
        // following switch to
        // avoid diffs in
        // tests/bits/mesh_3d_16
        ratio = 1. / GeometryInfo<3>::max_children_per_face;
        break;
      case internal::SubfaceCase<3>::case_x:
      case internal::SubfaceCase<3>::case_y:
        ratio = 0.5;
        break;
      case internal::SubfaceCase<3>::case_xy:
      case internal::SubfaceCase<3>::case_x1y2y:
      case internal::SubfaceCase<3>::case_y1x2x:
        ratio = 0.25;
        break;
      case internal::SubfaceCase<3>::case_x1y:
      case internal::SubfaceCase<3>::case_y1x:
        if (subface_no < 2)
          ratio = 0.25;
        else
          ratio = 0.5;
        break;
      case internal::SubfaceCase<3>::case_x2y:
      case internal::SubfaceCase<3>::case_y2x:
        if (subface_no == 0)
          ratio = 0.5;
        else
          ratio = 0.25;
        break;
      default:
        // there should be no
        // cases left
        Assert(false, ExcInternalError());
        break;
    }

  return ratio;
}



template <int dim>
RefinementCase<dim - 1> inline GeometryInfo<dim>::face_refinement_case(
  const RefinementCase<dim> &,
  const unsigned int,
  const bool,
  const bool,
  const bool)
{
  Assert(false, ExcNotImplemented());
  return RefinementCase<dim - 1>::no_refinement;
}

template <>
RefinementCase<0> inline GeometryInfo<1>::face_refinement_case(
  const RefinementCase<1> &,
  const unsigned int,
  const bool,
  const bool,
  const bool)
{
  Assert(false, ExcImpossibleInDim(1));

  return RefinementCase<0>::no_refinement;
}


template <>
inline RefinementCase<1>
GeometryInfo<2>::face_refinement_case(
  const RefinementCase<2> &cell_refinement_case,
  const unsigned int       face_no,
  const bool,
  const bool,
  const bool)
{
  const unsigned int dim = 2;
  AssertIndexRange(cell_refinement_case,
                   RefinementCase<dim>::isotropic_refinement + 1);
  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);

  const RefinementCase<dim - 1>
    ref_cases[RefinementCase<dim>::isotropic_refinement +
              1][GeometryInfo<dim>::faces_per_cell / 2] = {
      {RefinementCase<dim - 1>::no_refinement, // no_refinement
       RefinementCase<dim - 1>::no_refinement},

      {RefinementCase<dim - 1>::no_refinement, RefinementCase<dim - 1>::cut_x},

      {RefinementCase<dim - 1>::cut_x, RefinementCase<dim - 1>::no_refinement},

      {RefinementCase<dim - 1>::cut_x, // cut_xy
       RefinementCase<dim - 1>::cut_x}};

  return ref_cases[cell_refinement_case][face_no / 2];
}


template <>
inline RefinementCase<2>
GeometryInfo<3>::face_refinement_case(
  const RefinementCase<3> &cell_refinement_case,
  const unsigned int       face_no,
  const bool               face_orientation,
  const bool  /*face_flip*/ ,
  const bool face_rotation)
{
  const unsigned int dim = 3;
  AssertIndexRange(cell_refinement_case,
                   RefinementCase<dim>::isotropic_refinement + 1);
  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);

  const RefinementCase<dim - 1>
    ref_cases[RefinementCase<dim>::isotropic_refinement + 1]
             [GeometryInfo<dim>::faces_per_cell / 2] = {
               {RefinementCase<dim - 1>::no_refinement, // no_refinement
                RefinementCase<dim - 1>::no_refinement,
                RefinementCase<dim - 1>::no_refinement},

               {RefinementCase<dim - 1>::no_refinement, // cut_x
                RefinementCase<dim - 1>::cut_y,
                RefinementCase<dim - 1>::cut_x},

               {RefinementCase<dim - 1>::cut_x, // cut_y
                RefinementCase<dim - 1>::no_refinement,
                RefinementCase<dim - 1>::cut_y},

               {RefinementCase<dim - 1>::cut_x, // cut_xy
                RefinementCase<dim - 1>::cut_y,
                RefinementCase<dim - 1>::cut_xy},

               {RefinementCase<dim - 1>::cut_y, // cut_z
                RefinementCase<dim - 1>::cut_x,
                RefinementCase<dim - 1>::no_refinement},

               {RefinementCase<dim - 1>::cut_y, // cut_xz
                RefinementCase<dim - 1>::cut_xy,
                RefinementCase<dim - 1>::cut_x},

               {RefinementCase<dim - 1>::cut_xy, // cut_yz
                RefinementCase<dim - 1>::cut_x,
                RefinementCase<dim - 1>::cut_y},

               {RefinementCase<dim - 1>::cut_xy, // cut_xyz
                RefinementCase<dim - 1>::cut_xy,
                RefinementCase<dim - 1>::cut_xy},
             };

  const RefinementCase<dim - 1> ref_case =
    ref_cases[cell_refinement_case][face_no / 2];

  const RefinementCase<dim - 1> flip[4] = {
    RefinementCase<dim - 1>::no_refinement,
    RefinementCase<dim - 1>::cut_y,
    RefinementCase<dim - 1>::cut_x,
    RefinementCase<dim - 1>::cut_xy};

  // correct the ref_case for face_orientation
  // and face_rotation. for face_orientation,
  // 'true' is the default value whereas for
  // face_rotation, 'false' is standard. If
  // <tt>face_rotation==face_orientation</tt>,
  // then one of them is non-standard and we
  // have to swap cut_x and cut_y, otherwise no
  // change is necessary.  face_flip has no
  // influence. however, in order to keep the
  // interface consistent with other functions,
  // we still include it as an argument to this
  // function
  return (face_orientation == face_rotation) ? flip[ref_case] : ref_case;
}



template <int dim>
inline RefinementCase<1>
GeometryInfo<dim>::line_refinement_case(const RefinementCase<dim> &,
                                        const unsigned int)
{
  Assert(false, ExcNotImplemented());
  return RefinementCase<1>::no_refinement;
}

template <>
inline RefinementCase<1>
GeometryInfo<1>::line_refinement_case(
  const RefinementCase<1> &cell_refinement_case,
  const unsigned int       line_no)
{
  (void)line_no;
  const unsigned int dim = 1;
  (void)dim;
  AssertIndexRange(cell_refinement_case,
                   RefinementCase<dim>::isotropic_refinement + 1);
  AssertIndexRange(line_no, GeometryInfo<dim>::lines_per_cell);

  return cell_refinement_case;
}


template <>
inline RefinementCase<1>
GeometryInfo<2>::line_refinement_case(
  const RefinementCase<2> &cell_refinement_case,
  const unsigned int       line_no)
{
  // Assertions are in face_refinement_case()
  return face_refinement_case(cell_refinement_case, line_no);
}


template <>
inline RefinementCase<1>
GeometryInfo<3>::line_refinement_case(
  const RefinementCase<3> &cell_refinement_case,
  const unsigned int       line_no)
{
  const unsigned int dim = 3;
  AssertIndexRange(cell_refinement_case,
                   RefinementCase<dim>::isotropic_refinement + 1);
  AssertIndexRange(line_no, GeometryInfo<dim>::lines_per_cell);

  // array indicating, which simple refine
  // case cuts a line in direction x, y or
  // z. For example, cut_y and everything
  // containing cut_y (cut_xy, cut_yz,
  // cut_xyz) cuts lines, which are in y
  // direction.
  const RefinementCase<dim> cut_one[dim] = {RefinementCase<dim>::cut_x,
                                            RefinementCase<dim>::cut_y,
                                            RefinementCase<dim>::cut_z};

  // order the direction of lines
  // 0->x, 1->y, 2->z
  const unsigned int direction[lines_per_cell] = {
    1, 1, 0, 0, 1, 1, 0, 0, 2, 2, 2, 2};

  return ((cell_refinement_case & cut_one[direction[line_no]]) ?
            RefinementCase<1>::cut_x :
            RefinementCase<1>::no_refinement);
}



template <int dim>
inline RefinementCase<dim>
GeometryInfo<dim>::min_cell_refinement_case_for_face_refinement(
  const RefinementCase<dim - 1> &,
  const unsigned int,
  const bool,
  const bool,
  const bool)
{
  Assert(false, ExcNotImplemented());

  return RefinementCase<dim>::no_refinement;
}

template <>
inline RefinementCase<1>
GeometryInfo<1>::min_cell_refinement_case_for_face_refinement(
  const RefinementCase<0> &,
  const unsigned int,
  const bool,
  const bool,
  const bool)
{
  const unsigned int dim = 1;
  Assert(false, ExcImpossibleInDim(dim));

  return RefinementCase<dim>::no_refinement;
}


template <>
inline RefinementCase<2>
GeometryInfo<2>::min_cell_refinement_case_for_face_refinement(
  const RefinementCase<1> &face_refinement_case,
  const unsigned int       face_no,
  const bool,
  const bool,
  const bool)
{
  const unsigned int dim = 2;
  AssertIndexRange(face_refinement_case,
                   RefinementCase<dim - 1>::isotropic_refinement + 1);
  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);

  if (face_refinement_case == RefinementCase<dim>::cut_x)
    return (face_no / 2) ? RefinementCase<dim>::cut_x :
                           RefinementCase<dim>::cut_y;
  else
    return RefinementCase<dim>::no_refinement;
}


template <>
inline RefinementCase<3>
GeometryInfo<3>::min_cell_refinement_case_for_face_refinement(
  const RefinementCase<2> &face_refinement_case,
  const unsigned int       face_no,
  const bool               face_orientation,
  const bool  /*face_flip*/ ,
  const bool face_rotation)
{
  const unsigned int dim = 3;
  AssertIndexRange(face_refinement_case,
                   RefinementCase<dim - 1>::isotropic_refinement + 1);
  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);

  const RefinementCase<2> flip[4] = {RefinementCase<2>::no_refinement,
                                     RefinementCase<2>::cut_y,
                                     RefinementCase<2>::cut_x,
                                     RefinementCase<2>::cut_xy};

  // correct the face_refinement_case for
  // face_orientation and face_rotation. for
  // face_orientation, 'true' is the default
  // value whereas for face_rotation, 'false'
  // is standard. If
  // <tt>face_rotation==face_orientation</tt>,
  // then one of them is non-standard and we
  // have to swap cut_x and cut_y, otherwise no
  // change is necessary.  face_flip has no
  // influence. however, in order to keep the
  // interface consistent with other functions,
  // we still include it as an argument to this
  // function
  const RefinementCase<dim - 1> std_face_ref =
    (face_orientation == face_rotation) ? flip[face_refinement_case] :
                                          face_refinement_case;

  const RefinementCase<dim> face_to_cell[3][4] = {
    {RefinementCase<dim>::no_refinement, // faces 0 and 1
     RefinementCase<dim>::cut_y, // cut_x in face 0 means cut_y for the cell
     RefinementCase<dim>::cut_z,
     RefinementCase<dim>::cut_yz},

    {RefinementCase<dim>::no_refinement, // faces 2 and 3 (note that x and y are
                                         // "exchanged on faces 2 and 3")
     RefinementCase<dim>::cut_z,
     RefinementCase<dim>::cut_x,
     RefinementCase<dim>::cut_xz},

    {RefinementCase<dim>::no_refinement, // faces 4 and 5
     RefinementCase<dim>::cut_x,
     RefinementCase<dim>::cut_y,
     RefinementCase<dim>::cut_xy}};

  return face_to_cell[face_no / 2][std_face_ref];
}



template <int dim>
inline RefinementCase<dim>
GeometryInfo<dim>::min_cell_refinement_case_for_line_refinement(
  const unsigned int)
{
  Assert(false, ExcNotImplemented());

  return RefinementCase<dim>::no_refinement;
}

template <>
inline RefinementCase<1>
GeometryInfo<1>::min_cell_refinement_case_for_line_refinement(
  const unsigned int line_no)
{
  (void)line_no;
  AssertIndexRange(line_no, 1);

  return RefinementCase<1>::cut_x;
}


template <>
inline RefinementCase<2>
GeometryInfo<2>::min_cell_refinement_case_for_line_refinement(
  const unsigned int line_no)
{
  const unsigned int dim = 2;
  (void)dim;
  AssertIndexRange(line_no, GeometryInfo<dim>::lines_per_cell);

  return (line_no / 2) ? RefinementCase<2>::cut_x : RefinementCase<2>::cut_y;
}


template <>
inline RefinementCase<3>
GeometryInfo<3>::min_cell_refinement_case_for_line_refinement(
  const unsigned int line_no)
{
  const unsigned int dim = 3;
  AssertIndexRange(line_no, GeometryInfo<dim>::lines_per_cell);

  const RefinementCase<dim> ref_cases[6] = {
    RefinementCase<dim>::cut_y,  // lines  0 and  1
    RefinementCase<dim>::cut_x,  // lines  2 and  3
    RefinementCase<dim>::cut_y,  // lines  4 and  5
    RefinementCase<dim>::cut_x,  // lines  6 and  7
    RefinementCase<dim>::cut_z,  // lines  8 and  9
    RefinementCase<dim>::cut_z}; // lines 10 and 11

  return ref_cases[line_no / 2];
}



template <>
inline unsigned int
GeometryInfo<3>::real_to_standard_face_vertex(const unsigned int vertex,
                                              const bool face_orientation,
                                              const bool face_flip,
                                              const bool face_rotation)
{
  AssertIndexRange(vertex, GeometryInfo<3>::vertices_per_face);

  // set up a table to make sure that
  // we handle non-standard faces correctly
  //
  // so set up a table that for each vertex (of
  // a quad in standard position) describes
  // which vertex to take
  //
  // first index: four vertices 0...3
  //
  // second index: face_orientation; 0:
  // opposite normal, 1: standard
  //
  // third index: face_flip; 0: standard, 1:
  // face rotated by 180 degrees
  //
  // forth index: face_rotation: 0: standard,
  // 1: face rotated by 90 degrees

  const unsigned int vertex_translation[4][2][2][2] = {
    {{{0, 2},   // vertex 0, face_orientation=false, face_flip=false,
                // face_rotation=false and true
      {3, 1}},  // vertex 0, face_orientation=false, face_flip=true,
                // face_rotation=false and true
     {{0, 1},   // vertex 0, face_orientation=true, face_flip=false,
                // face_rotation=false and true
      {3, 2}}}, // vertex 0, face_orientation=true, face_flip=true,
                // face_rotation=false and true

    {{{2, 3}, // vertex 1 ...
      {1, 0}},
     {{1, 3}, {2, 0}}},

    {{{1, 0}, // vertex 2 ...
      {2, 3}},
     {{2, 0}, {1, 3}}},

    {{{3, 1}, // vertex 3 ...
      {0, 2}},
     {{3, 2}, {0, 1}}}};

  return vertex_translation[vertex][face_orientation][face_flip][face_rotation];
}



template <int dim>
inline unsigned int
GeometryInfo<dim>::real_to_standard_face_vertex(const unsigned int vertex,
                                                const bool,
                                                const bool,
                                                const bool)
{
  Assert(dim > 1, ExcImpossibleInDim(dim));
  AssertIndexRange(vertex, GeometryInfo<dim>::vertices_per_face);
  return vertex;
}



template <>
inline unsigned int
GeometryInfo<3>::standard_to_real_face_line(const unsigned int line,
                                            const bool         face_orientation,
                                            const bool         face_flip,
                                            const bool         face_rotation)
{
  AssertIndexRange(line, GeometryInfo<3>::lines_per_face);


  // make sure we handle
  // non-standard faces correctly
  //
  // so set up a table that for each line (of a
  // quad) describes which line to take
  //
  // first index: four lines 0...3
  //
  // second index: face_orientation; 0:
  // opposite normal, 1: standard
  //
  // third index: face_flip; 0: standard, 1:
  // face rotated by 180 degrees
  //
  // forth index: face_rotation: 0: standard,
  // 1: face rotated by 90 degrees

  const unsigned int line_translation[4][2][2][2] = {
    {{{2, 0},   // line 0, face_orientation=false, face_flip=false,
                // face_rotation=false and true
      {3, 1}},  // line 0, face_orientation=false, face_flip=true,
                // face_rotation=false and true
     {{0, 3},   // line 0, face_orientation=true, face_flip=false,
                // face_rotation=false and true
      {1, 2}}}, // line 0, face_orientation=true, face_flip=true,
                // face_rotation=false and true

    {{{3, 1}, // line 1 ...
      {2, 0}},
     {{1, 2}, {0, 3}}},

    {{{0, 3}, // line 2 ...
      {1, 2}},
     {{2, 0}, {3, 1}}},

    {{{1, 2}, // line 3 ...
      {0, 3}},
     {{3, 1}, {2, 0}}}};

  return line_translation[line][face_orientation][face_flip][face_rotation];
}



template <int dim>
inline unsigned int
GeometryInfo<dim>::standard_to_real_face_line(const unsigned int line,
                                              const bool,
                                              const bool,
                                              const bool)
{
  Assert(false, ExcNotImplemented());
  return line;
}



template <>
inline unsigned int
GeometryInfo<2>::standard_to_real_line_vertex(const unsigned int vertex,
                                              const bool line_orientation)
{
  return line_orientation ? vertex : (1 - vertex);
}



template <int dim>
inline unsigned int
GeometryInfo<dim>::standard_to_real_line_vertex(const unsigned int vertex,
                                                const bool)
{
  Assert(false, ExcNotImplemented());
  return vertex;
}



template <>
inline std::array<unsigned int, 2>
GeometryInfo<2>::standard_quad_vertex_to_line_vertex_index(
  const unsigned int vertex)
{
  return {{vertex % 2, vertex / 2}};
}



template <int dim>
inline std::array<unsigned int, 2>
GeometryInfo<dim>::standard_quad_vertex_to_line_vertex_index(
  const unsigned int vertex)
{
  Assert(false, ExcNotImplemented());
  (void)vertex;
  return {{0, 0}};
}



template <>
inline std::array<unsigned int, 2>
GeometryInfo<3>::standard_hex_line_to_quad_line_index(const unsigned int i)
{
  // set up a table that for each
  // line describes a) from which
  // quad to take it, b) which line
  // therein it is if the face is
  // oriented correctly
  static const unsigned int lookup_table[GeometryInfo<3>::lines_per_cell][2] = {
    {4, 0}, // take first four lines from bottom face
    {4, 1},
    {4, 2},
    {4, 3},

    {5, 0}, // second four lines from top face
    {5, 1},
    {5, 2},
    {5, 3},

    {0, 0}, // the rest randomly
    {1, 0},
    {0, 1},
    {1, 1}};

  return {{lookup_table[i][0], lookup_table[i][1]}};
}



template <int dim>
inline std::array<unsigned int, 2>
GeometryInfo<dim>::standard_hex_line_to_quad_line_index(const unsigned int line)
{
  Assert(false, ExcNotImplemented());
  (void)line;
  return {{0, 0}};
}



template <>
inline std::array<unsigned int, 2>
GeometryInfo<3>::standard_hex_vertex_to_quad_vertex_index(
  const unsigned int vertex)
{
  // get the corner indices by asking either the bottom or the top face for its
  // vertices. handle non-standard faces by calling the vertex reordering
  // function from GeometryInfo

  // bottom face (4) for first four vertices, top face (5) for the rest
  return {{4 + vertex / 4, vertex % 4}};
}



template <int dim>
inline std::array<unsigned int, 2>
GeometryInfo<dim>::standard_hex_vertex_to_quad_vertex_index(
  const unsigned int vertex)
{
  Assert(false, ExcNotImplemented());
  (void)vertex;
  return {{0, 0}};
}



template <>
inline unsigned int
GeometryInfo<3>::real_to_standard_face_line(const unsigned int line,
                                            const bool         face_orientation,
                                            const bool         face_flip,
                                            const bool         face_rotation)
{
  AssertIndexRange(line, GeometryInfo<3>::lines_per_face);


  // make sure we handle
  // non-standard faces correctly
  //
  // so set up a table that for each line (of a
  // quad) describes which line to take
  //
  // first index: four lines 0...3
  //
  // second index: face_orientation; 0:
  // opposite normal, 1: standard
  //
  // third index: face_flip; 0: standard, 1:
  // face rotated by 180 degrees
  //
  // forth index: face_rotation: 0: standard,
  // 1: face rotated by 90 degrees

  const unsigned int line_translation[4][2][2][2] = {
    {{{2, 0},   // line 0, face_orientation=false, face_flip=false,
                // face_rotation=false and true
      {3, 1}},  // line 0, face_orientation=false, face_flip=true,
                // face_rotation=false and true
     {{0, 2},   // line 0, face_orientation=true, face_flip=false,
                // face_rotation=false and true
      {1, 3}}}, // line 0, face_orientation=true, face_flip=true,
                // face_rotation=false and true

    {{{3, 1}, // line 1 ...
      {2, 0}},
     {{1, 3}, {0, 2}}},

    {{{0, 3}, // line 2 ...
      {1, 2}},
     {{2, 1}, {3, 0}}},

    {{{1, 2}, // line 3 ...
      {0, 3}},
     {{3, 0}, {2, 1}}}};

  return line_translation[line][face_orientation][face_flip][face_rotation];
}



template <int dim>
inline unsigned int
GeometryInfo<dim>::real_to_standard_face_line(const unsigned int line,
                                              const bool,
                                              const bool,
                                              const bool)
{
  Assert(false, ExcNotImplemented());
  return line;
}



template <>
inline unsigned int
GeometryInfo<1>::child_cell_on_face(const RefinementCase<1> &,
                                    const unsigned int face,
                                    const unsigned int subface,
                                    const bool,
                                    const bool,
                                    const bool,
                                    const RefinementCase<0> &)
{
  (void)subface;
  AssertIndexRange(face, faces_per_cell);
  AssertIndexRange(subface, max_children_per_face);

  return face;
}



template <>
inline unsigned int
GeometryInfo<2>::child_cell_on_face(const RefinementCase<2> &ref_case,
                                    const unsigned int       face,
                                    const unsigned int       subface,
                                    const bool  /*face_orientation*/ ,
                                    const bool face_flip,
                                    const bool  /*face_rotation*/ ,
                                    const RefinementCase<1> &)
{
  AssertIndexRange(face, faces_per_cell);
  AssertIndexRange(subface, max_children_per_face);

  // always return the child adjacent to the specified
  // subface. if the face of a cell is not refined, don't
  // throw an assertion but deliver the child adjacent to
  // the face nevertheless, i.e. deliver the child of
  // this cell adjacent to the subface of a possibly
  // refined neighbor. this simplifies setting neighbor
  // information in execute_refinement.
  const unsigned int
    subcells[2][RefinementCase<2>::isotropic_refinement][faces_per_cell]
            [max_children_per_face] = {
              {
                // Normal orientation (face_flip = false)
                {{0, 0}, {1, 1}, {0, 1}, {0, 1}}, // cut_x
                {{0, 1}, {0, 1}, {0, 0}, {1, 1}}, // cut_y
                {{0, 2}, {1, 3}, {0, 1}, {2, 3}}  // cut_xy, i.e., isotropic
              },
              {
                // Flipped orientation (face_flip = true)
                {{0, 0}, {1, 1}, {1, 0}, {1, 0}}, // cut_x
                {{1, 0}, {1, 0}, {0, 0}, {1, 1}}, // cut_y
                {{2, 0}, {3, 1}, {1, 0}, {3, 2}}  // cut_xy, i.e., isotropic
              }};

  return subcells[face_flip][ref_case - 1][face][subface];
}



template <>
inline unsigned int
GeometryInfo<3>::child_cell_on_face(const RefinementCase<3> &ref_case,
                                    const unsigned int       face,
                                    const unsigned int       subface,
                                    const bool               face_orientation,
                                    const bool               face_flip,
                                    const bool               face_rotation,
                                    const RefinementCase<2> &face_ref_case)
{
  const unsigned int dim = 3;

  Assert(ref_case > RefinementCase<dim - 1>::no_refinement,
         ExcMessage("Cell has no children."));
  AssertIndexRange(face, faces_per_cell);
  if (!(subface == 0 &&
        face_ref_case == RefinementCase<dim - 1>::no_refinement))
    {
      AssertIndexRange(subface,
                       GeometryInfo<dim - 1>::n_children(face_ref_case));
    }

  // invalid number used for invalid cases,
  // e.g. when the children are more refined at
  // a given face than the face itself
  const unsigned int e = numbers::invalid_unsigned_int;

  // the whole process of finding a child cell
  // at a given subface considering the
  // possibly anisotropic refinement cases of
  // the cell and the face as well as
  // orientation, flip and rotation of the face
  // is quite complicated. thus, we break it
  // down into several steps.

  // first step: convert the given face refine
  // case to a face refine case concerning the
  // face in standard orientation (, flip and
  // rotation). This only affects cut_x and
  // cut_y
  const RefinementCase<dim - 1> flip[4] = {
    RefinementCase<dim - 1>::no_refinement,
    RefinementCase<dim - 1>::cut_y,
    RefinementCase<dim - 1>::cut_x,
    RefinementCase<dim - 1>::cut_xy};
  // for face_orientation, 'true' is the
  // default value whereas for face_rotation,
  // 'false' is standard. If
  // <tt>face_rotation==face_orientation</tt>,
  // then one of them is non-standard and we
  // have to swap cut_x and cut_y, otherwise no
  // change is necessary.
  const RefinementCase<dim - 1> std_face_ref =
    (face_orientation == face_rotation) ? flip[face_ref_case] : face_ref_case;

  // second step: convert the given subface
  // index to the one for a standard face
  // respecting face_orientation, face_flip and
  // face_rotation

  // first index:  face_ref_case
  // second index: face_orientation
  // third index:  face_flip
  // forth index:  face_rotation
  // fifth index:  subface index
  const unsigned int subface_exchange[4][2][2][2][4] = {
    // no_refinement (subface 0 stays 0,
    // all others are invalid)
    {{{{0, e, e, e}, {0, e, e, e}}, {{0, e, e, e}, {0, e, e, e}}},
     {{{0, e, e, e}, {0, e, e, e}}, {{0, e, e, e}, {0, e, e, e}}}},
    // cut_x (here, if the face is only
    // rotated OR only falsely oriented,
    // then subface 0 of the non-standard
    // face does NOT correspond to one of
    // the subfaces of a standard
    // face. Thus we indicate the subface
    // which is located at the lower left
    // corner (the origin of the face's
    // local coordinate system) with
    // '0'. The rest of this issue is
    // taken care of using the above
    // conversion to a 'standard face
    // refine case')
    {{{{0, 1, e, e}, {0, 1, e, e}}, {{1, 0, e, e}, {1, 0, e, e}}},
     {{{0, 1, e, e}, {0, 1, e, e}}, {{1, 0, e, e}, {1, 0, e, e}}}},
    // cut_y (the same applies as for
    // cut_x)
    {{{{0, 1, e, e}, {1, 0, e, e}}, {{1, 0, e, e}, {0, 1, e, e}}},
     {{{0, 1, e, e}, {1, 0, e, e}}, {{1, 0, e, e}, {0, 1, e, e}}}},
    // cut_xyz: this information is
    // identical to the information
    // returned by
    // GeometryInfo<3>::real_to_standard_face_vertex()
    {{{{0, 2, 1, 3},     // face_orientation=false, face_flip=false,
                         // face_rotation=false, subfaces 0,1,2,3
       {2, 3, 0, 1}},    // face_orientation=false, face_flip=false,
                         // face_rotation=true,  subfaces 0,1,2,3
      {{3, 1, 2, 0},     // face_orientation=false, face_flip=true,
                         // face_rotation=false, subfaces 0,1,2,3
       {1, 0, 3, 2}}},   // face_orientation=false, face_flip=true,
                         // face_rotation=true,  subfaces 0,1,2,3
     {{{0, 1, 2, 3},     // face_orientation=true,  face_flip=false,
                         // face_rotation=false, subfaces 0,1,2,3
       {1, 3, 0, 2}},    // face_orientation=true,  face_flip=false,
                         // face_rotation=true,  subfaces 0,1,2,3
      {{3, 2, 1, 0},     // face_orientation=true,  face_flip=true,
                         // face_rotation=false, subfaces 0,1,2,3
       {2, 0, 3, 1}}}}}; // face_orientation=true,  face_flip=true,
                         // face_rotation=true,  subfaces 0,1,2,3

  const unsigned int std_subface =
    subface_exchange[face_ref_case][face_orientation][face_flip][face_rotation]
                    [subface];
  Assert(std_subface != e, ExcInternalError());

  // third step: these are the children, which
  // can be found at the given subfaces of an
  // isotropically refined (standard) face
  //
  // first index:  (refinement_case-1)
  // second index: face_index
  // third index:  subface_index (isotropic refinement)
  const unsigned int iso_children[RefinementCase<dim>::cut_xyz][faces_per_cell]
                                 [max_children_per_face] = {
                                   // cut_x
                                   {{0, 0, 0, 0},  // face 0, subfaces 0,1,2,3
                                    {1, 1, 1, 1},  // face 1, subfaces 0,1,2,3
                                    {0, 0, 1, 1},  // face 2, subfaces 0,1,2,3
                                    {0, 0, 1, 1},  // face 3, subfaces 0,1,2,3
                                    {0, 1, 0, 1},  // face 4, subfaces 0,1,2,3
                                    {0, 1, 0, 1}}, // face 5, subfaces 0,1,2,3
                                                   // cut_y
                                   {{0, 1, 0, 1},
                                    {0, 1, 0, 1},
                                    {0, 0, 0, 0},
                                    {1, 1, 1, 1},
                                    {0, 0, 1, 1},
                                    {0, 0, 1, 1}},
                                   // cut_xy
                                   {{0, 2, 0, 2},
                                    {1, 3, 1, 3},
                                    {0, 0, 1, 1},
                                    {2, 2, 3, 3},
                                    {0, 1, 2, 3},
                                    {0, 1, 2, 3}},
                                   // cut_z
                                   {{0, 0, 1, 1},
                                    {0, 0, 1, 1},
                                    {0, 1, 0, 1},
                                    {0, 1, 0, 1},
                                    {0, 0, 0, 0},
                                    {1, 1, 1, 1}},
                                   // cut_xz
                                   {{0, 0, 1, 1},
                                    {2, 2, 3, 3},
                                    {0, 1, 2, 3},
                                    {0, 1, 2, 3},
                                    {0, 2, 0, 2},
                                    {1, 3, 1, 3}},
                                   // cut_yz
                                   {{0, 1, 2, 3},
                                    {0, 1, 2, 3},
                                    {0, 2, 0, 2},
                                    {1, 3, 1, 3},
                                    {0, 0, 1, 1},
                                    {2, 2, 3, 3}},
                                   // cut_xyz
                                   {{0, 2, 4, 6},
                                    {1, 3, 5, 7},
                                    {0, 4, 1, 5},
                                    {2, 6, 3, 7},
                                    {0, 1, 2, 3},
                                    {4, 5, 6, 7}}};

  // forth step: check, whether the given face
  // refine case is valid for the given cell
  // refine case. this is the case, if the
  // given face refine case is at least as
  // refined as the face is for the given cell
  // refine case

  // note, that we are considering standard
  // face refinement cases here and thus must
  // not pass the given orientation, flip and
  // rotation flags
  if ((std_face_ref & face_refinement_case(ref_case, face)) ==
      face_refinement_case(ref_case, face))
    {
      // all is fine. for anisotropic face
      // refine cases, select one of the
      // isotropic subfaces which neighbors the
      // same child

      // first index: (standard) face refine case
      // second index: subface index
      const unsigned int equivalent_iso_subface[4][4] = {
        {0, e, e, e},  // no_refinement
        {0, 3, e, e},  // cut_x
        {0, 3, e, e},  // cut_y
        {0, 1, 2, 3}}; // cut_xy

      const unsigned int equ_std_subface =
        equivalent_iso_subface[std_face_ref][std_subface];
      Assert(equ_std_subface != e, ExcInternalError());

      return iso_children[ref_case - 1][face][equ_std_subface];
    }
  else
    {
      // the face_ref_case was too coarse,
      // throw an error
      Assert(false,
             ExcMessage("The face RefineCase is too coarse "
                        "for the given cell RefineCase."));
    }
  // we only get here in case of an error
  return e;
}



template <>
inline unsigned int
GeometryInfo<4>::child_cell_on_face(const RefinementCase<4> &,
                                    const unsigned int,
                                    const unsigned int,
                                    const bool,
                                    const bool,
                                    const bool,
                                    const RefinementCase<3> &)
{
  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



template <>
inline unsigned int
GeometryInfo<2>::face_to_cell_lines(const unsigned int face,
                                    const unsigned int line,
                                    const bool,
                                    const bool,
                                    const bool)
{
  (void)line;
  AssertIndexRange(face, faces_per_cell);
  AssertIndexRange(line, lines_per_face);

  // The face is a line itself.
  return face;
}



template <>
inline unsigned int
GeometryInfo<3>::face_to_cell_lines(const unsigned int face,
                                    const unsigned int line,
                                    const bool         face_orientation,
                                    const bool         face_flip,
                                    const bool         face_rotation)
{
  AssertIndexRange(face, faces_per_cell);
  AssertIndexRange(line, lines_per_face);

  const unsigned lines[faces_per_cell][lines_per_face] = {
    {8, 10, 0, 4},  // left face
    {9, 11, 1, 5},  // right face
    {2, 6, 8, 9},   // front face
    {3, 7, 10, 11}, // back face
    {0, 1, 2, 3},   // bottom face
    {4, 5, 6, 7}};  // top face
  return lines[face][real_to_standard_face_line(
    line, face_orientation, face_flip, face_rotation)];
}



inline unsigned int
GeometryInfo<0>::face_to_cell_lines(const unsigned int,
                                    const unsigned int,
                                    const bool,
                                    const bool,
                                    const bool)
{
  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



template <int dim>
inline unsigned int
GeometryInfo<dim>::face_to_cell_lines(const unsigned int,
                                      const unsigned int,
                                      const bool,
                                      const bool,
                                      const bool)
{
  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



template <int dim>
inline unsigned int
GeometryInfo<dim>::face_to_cell_vertices(const unsigned int face,
                                         const unsigned int vertex,
                                         const bool         face_orientation,
                                         const bool         face_flip,
                                         const bool         face_rotation)
{
  return child_cell_on_face(RefinementCase<dim>::isotropic_refinement,
                            face,
                            vertex,
                            face_orientation,
                            face_flip,
                            face_rotation);
}



inline unsigned int
GeometryInfo<0>::face_to_cell_vertices(const unsigned int,
                                       const unsigned int,
                                       const bool,
                                       const bool,
                                       const bool)
{
  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



template <int dim>
template <typename Number>
inline Point<dim, Number>
GeometryInfo<dim>::project_to_unit_cell(const Point<dim, Number> &q)
{
  Point<dim, Number> p;
  for (unsigned int i = 0; i < dim; i++)
    p[i] = std::min(std::max(q[i], Number(0.)), Number(1.));

  return p;
}



template <int dim>
inline double
GeometryInfo<dim>::distance_to_unit_cell(const Point<dim> &p)
{
  double result = 0.0;

  for (unsigned int i = 0; i < dim; i++)
    {
      result = std::max(result, -p[i]);
      result = std::max(result, p[i] - 1.);
    }

  return result;
}



template <int dim>
inline double
GeometryInfo<dim>::d_linear_shape_function(const Point<dim> & xi,
                                           const unsigned int i)
{
  AssertIndexRange(i, GeometryInfo<dim>::vertices_per_cell);

  switch (dim)
    {
      case 1:
        {
          const double x = xi[0];
          switch (i)
            {
              case 0:
                return 1 - x;
              case 1:
                return x;
            }
          break;
        }

      case 2:
        {
          const double x = xi[0];
          const double y = xi[1];
          switch (i)
            {
              case 0:
                return (1 - x) * (1 - y);
              case 1:
                return x * (1 - y);
              case 2:
                return (1 - x) * y;
              case 3:
                return x * y;
            }
          break;
        }

      case 3:
        {
          const double x = xi[0];
          const double y = xi[1];
          const double z = xi[2];
          switch (i)
            {
              case 0:
                return (1 - x) * (1 - y) * (1 - z);
              case 1:
                return x * (1 - y) * (1 - z);
              case 2:
                return (1 - x) * y * (1 - z);
              case 3:
                return x * y * (1 - z);
              case 4:
                return (1 - x) * (1 - y) * z;
              case 5:
                return x * (1 - y) * z;
              case 6:
                return (1 - x) * y * z;
              case 7:
                return x * y * z;
            }
          break;
        }

      default:
        Assert(false, ExcNotImplemented());
    }
  return -1e9;
}



template <>
Tensor<1, 1> inline GeometryInfo<1>::d_linear_shape_function_gradient(
  const Point<1> &,
  const unsigned int i)
{
  AssertIndexRange(i, GeometryInfo<1>::vertices_per_cell);

  switch (i)
    {
      case 0:
        return Point<1>(-1.);
      case 1:
        return Point<1>(1.);
    }

  return Point<1>(-1e9);
}



template <>
Tensor<1, 2> inline GeometryInfo<2>::d_linear_shape_function_gradient(
  const Point<2> &   xi,
  const unsigned int i)
{
  AssertIndexRange(i, GeometryInfo<2>::vertices_per_cell);

  const double x = xi[0];
  const double y = xi[1];
  switch (i)
    {
      case 0:
        return Point<2>(-(1 - y), -(1 - x));
      case 1:
        return Point<2>(1 - y, -x);
      case 2:
        return Point<2>(-y, 1 - x);
      case 3:
        return Point<2>(y, x);
    }
  return Point<2>(-1e9, -1e9);
}



template <>
Tensor<1, 3> inline GeometryInfo<3>::d_linear_shape_function_gradient(
  const Point<3> &   xi,
  const unsigned int i)
{
  AssertIndexRange(i, GeometryInfo<3>::vertices_per_cell);

  const double x = xi[0];
  const double y = xi[1];
  const double z = xi[2];
  switch (i)
    {
      case 0:
        return Point<3>(-(1 - y) * (1 - z),
                        -(1 - x) * (1 - z),
                        -(1 - x) * (1 - y));
      case 1:
        return Point<3>((1 - y) * (1 - z), -x * (1 - z), -x * (1 - y));
      case 2:
        return Point<3>(-y * (1 - z), (1 - x) * (1 - z), -(1 - x) * y);
      case 3:
        return Point<3>(y * (1 - z), x * (1 - z), -x * y);
      case 4:
        return Point<3>(-(1 - y) * z, -(1 - x) * z, (1 - x) * (1 - y));
      case 5:
        return Point<3>((1 - y) * z, -x * z, x * (1 - y));
      case 6:
        return Point<3>(-y * z, (1 - x) * z, (1 - x) * y);
      case 7:
        return Point<3>(y * z, x * z, x * y);
    }

  return Point<3>(-1e9, -1e9, -1e9);
}



template <int dim>
inline Tensor<1, dim>
GeometryInfo<dim>::d_linear_shape_function_gradient(const Point<dim> &,
                                                    const unsigned int)
{
  Assert(false, ExcNotImplemented());
  return Tensor<1, dim>();
}


unsigned int inline GeometryInfo<0>::n_children(const RefinementCase<0> &)
{
  return 0;
}


namespace internal
{
  namespace GeometryInfoHelper
  {
    // wedge product of a single
    // vector in 2d: we just have to
    // rotate it by 90 degrees to the
    // right
    inline Tensor<1, 2>
    wedge_product(const Tensor<1, 2> (&derivative)[1])
    {
      Tensor<1, 2> result;
      result[0] = derivative[0][1];
      result[1] = -derivative[0][0];

      return result;
    }


    // wedge product of 2 vectors in
    // 3d is the cross product
    inline Tensor<1, 3>
    wedge_product(const Tensor<1, 3> (&derivative)[2])
    {
      return cross_product_3d(derivative[0], derivative[1]);
    }


    // wedge product of dim vectors
    // in dim-d: that's the
    // determinant of the matrix
    template <int dim>
    inline Tensor<0, dim>
    wedge_product(const Tensor<1, dim> (&derivative)[dim])
    {
      Tensor<2, dim> jacobian;
      for (unsigned int i = 0; i < dim; ++i)
        jacobian[i] = derivative[i];

      return determinant(jacobian);
    }
  } // namespace GeometryInfoHelper
} // namespace internal


template <int dim>
template <int spacedim>
inline void
GeometryInfo<dim>::alternating_form_at_vertices
#  ifndef DEAL_II_CXX14_CONSTEXPR_BUG
  (const Point<spacedim> (&vertices)[vertices_per_cell],
   Tensor<spacedim - dim, spacedim> (&forms)[vertices_per_cell])
#  else
  (const Point<spacedim> *vertices, Tensor<spacedim - dim, spacedim> *forms)
#  endif
{
  // for each of the vertices,
  // compute the alternating form
  // of the mapped unit
  // vectors. consider for
  // example the case of a quad
  // in spacedim==3: to do so, we
  // need to see how the
  // infinitesimal vectors
  // (d\xi_1,0) and (0,d\xi_2)
  // are transformed into
  // spacedim-dimensional space
  // and then form their cross
  // product (i.e. the wedge product
  // of two vectors). to this end, note
  // that
  //    \vec x = sum_i \vec v_i phi_i(\vec xi)
  // so the transformed vectors are
  //   [x(\xi+(d\xi_1,0))-x(\xi)]/d\xi_1
  // and
  //   [x(\xi+(0,d\xi_2))-x(\xi)]/d\xi_2
  // which boils down to the columns
  // of the 3x2 matrix \grad_\xi x(\xi)
  //
  // a similar reasoning would
  // hold for all dim,spacedim
  // pairs -- we only have to
  // compute the wedge product of
  // the columns of the
  // derivatives
  for (unsigned int i = 0; i < vertices_per_cell; ++i)
    {
      Tensor<1, spacedim> derivatives[dim];

      for (unsigned int j = 0; j < vertices_per_cell; ++j)
        {
          const Tensor<1, dim> grad_phi_j =
            d_linear_shape_function_gradient(unit_cell_vertex(i), j);
          for (unsigned int l = 0; l < dim; ++l)
            derivatives[l] += vertices[j] * grad_phi_j[l];
        }

      forms[i] = internal::GeometryInfoHelper::wedge_product(derivatives);
    }
}

#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif


