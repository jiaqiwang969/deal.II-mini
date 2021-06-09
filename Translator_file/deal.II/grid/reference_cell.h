//include/deal.II-translator/grid/reference_cell_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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

#ifndef dealii_tria_reference_cell_h
#define dealii_tria_reference_cell_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/ndarray.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class Mapping;

template <int dim>
class Quadrature;

class ReferenceCell;
#endif


namespace internal
{
  namespace ReferenceCell
  {
    /**
     * 一个辅助函数，用于从一个整数中创建一个ReferenceCell对象。ReferenceCell对象是
     * "单子"（实际上是 "多子"
     *
     * - 有多个，但它们只是少数，这些是所有可以使用的）。) 那么需要的是有一种方法来创建这些对象，用它们的内部id来区分存在的少数可能的对象。我们可以通过ReferenceCell的公共构造函数来做到这一点，但这将允许用户在我们设想的范围之外创建这些单元，而我们不想这样做。相反，接受一个整数的构造函数是 "私有 "的，但我们在内部命名空间中有这样一个函数，它是这个类的朋友，可以用来创建对象。
     *
     */
    DEAL_II_CONSTEXPR dealii::ReferenceCell
                      make_reference_cell_from_int(const std::uint8_t kind);

  } // namespace ReferenceCell
} // namespace internal



/**
 * 一个描述可以使用的参考单元的种类的类型。这包括四边形和六面体（即
 * "超立方体"）、三角形和四面体（单纯体），以及使用混合三维网格时必须的金字塔和楔形。
 * 这种类型的对象不应该在用户代码中创建，因此，除了默认的构造函数（会创建一个无效的对象），该类没有一个用户可以访问的构造函数。相反，在ReferenceCells命名空间中定义了有限数量的特定参考单元对象，完全列举了所有可能的值。因此，用户代码应该完全依赖于从这些特殊对象中分配ReferenceCell对象，并与这些特殊对象进行比较。
 * 该类的目的和意图在 @ref GlossReferenceCell "参考单元 "
 * 词汇表条目中描述。
 *
 *
 * @ingroup grid geomprimitives aniso
 *
 *
 */
class ReferenceCell
{
public:
  /**
   * 对于给定的结构尺寸和顶点数量，返回正确的ReferenceCell。例如，如果`dim==2`和`n_vertices==4`，这个函数将返回
   * ReferenceCells::Quadrilateral.
   * 但如果`dim==3`和`n_vertices==4`，它将返回
   * ReferenceCells::Tetrahedron.  。
   *
   */
  static ReferenceCell
  n_vertices_to_type(const int dim, const unsigned int n_vertices);

  /**
   * 默认构造函数。将此对象初始化为一个无效的对象。最终的结果是当前对象等于
   * ReferenceCells::Invalid.
   * 一般来说，ReferenceCell对象是通过从命名空间ReferenceCells中的特殊对象赋值创建的，这是获得有效对象的唯一方法。
   *
   */
  DEAL_II_CONSTEXPR
  ReferenceCell();

  /**
   * @name  查询有关参考单元格种类的信息  
     * @{ 
   *
   */

  /**
   * 如果对象是 ReferenceCells::Vertex,   ReferenceCells::Line,
   * ReferenceCells::Quadrilateral, 或 ReferenceCells::Hexahedron.
   * ，返回`true`。
   *
   */
  bool
  is_hyper_cube() const;

  /**
   * 如果对象是一个顶点、直线、三角形或四面体，则返回真。
   *
   */
  bool
  is_simplex() const;

  /**
   * 返回当前对象所代表的参考单元的尺寸。
   *
   */
  unsigned int
  get_dimension() const;

  /**
   * @}
   *
   */

  /**
   * @name  在参考单元上定义的形状函数、映射、四边形 
     * @{ 
   *
   */

  /**
   * 计算当前参考单元类型在位置 $\xi$ 的第 $i$
   * 个线性形状函数的值。
   *
   */
  template <int dim>
  double
  d_linear_shape_function(const Point<dim> &xi, const unsigned int i) const;

  /**
   * 计算当前参考单元类型在位置 $\xi$ 的 $i$
   * -th线性形状函数的梯度。
   *
   */
  template <int dim>
  Tensor<1, dim>
  d_linear_shape_function_gradient(const Point<dim> & xi,
                                   const unsigned int i) const;

  /**
   * 返回一个与当前参考单元相匹配的默认映射度 @p degree
   * 。如果这个参考单元是一个超立方体，那么返回的映射是一个MappingQGeneric；否则，它是一个用FE_SimplexP（如果参考单元是一个三角形或四面体）、用FE_PyramidP（如果参考单元是一个金字塔）或用FE_WedgeP（如果参考单元是一个楔子）初始化的MappingFE类型对象。
   *
   */
  template <int dim, int spacedim = dim>
  std::unique_ptr<Mapping<dim, spacedim>>
  get_default_mapping(const unsigned int degree) const;

  /**
   * 返回一个与当前参考单元相匹配的默认线性映射。
   * 如果这个参考单元是一个超立方体，那么返回的映射是一个MappingQ1；否则，它是一个用FE_SimplexP（如果参考单元是一个三角形或四面体）、FE_PyramidP（如果参考单元是一个金字塔）或FE_WedgeP（如果参考单元是一个楔子）初始化的MappingFE类型对象。
   * 换句话说，函数名称中的 "线性
   * "一词必须理解为对某些坐标方向的 $d$
   * -线性（即双线性或三线性）。
   *
   */
  template <int dim, int spacedim = dim>
  const Mapping<dim, spacedim> &
  get_default_linear_mapping() const;

  /**
   * 返回一个与给定参考单元相匹配的高斯型正交（QGauss,
   * QGaussSimplex, QGaussPyramid, QGaussWedge）。      @param[in]
   * n_points_1D
   * 每个方向上的正交点的数量（QGauss），或者对于其他类型的正交点，需要准确地集成什么多项式的指示。
   *
   */
  template <int dim>
  Quadrature<dim>
  get_gauss_type_quadrature(const unsigned n_points_1D) const;

  /**
   * 返回一个具有给定参考单元的支持点的正交规则。
   * @note 正交对象的权重是不填的。
   *
   */
  template <int dim>
  const Quadrature<dim> &
  get_nodal_type_quadrature() const;

  /**
   * @}
   *
   */

  /**
   * @name  查询参考单元的构件数量 
     * @{ 
   *
   */

  /**
   * 返回构成有关参考单元的顶点的数量。一个顶点是参考单元的一个
   * "角"（一个零维物体）。
   *
   */
  unsigned int
  n_vertices() const;

  /**
   * 返回一个对象，可以看作是一个包含从零到n_vertices()所有索引的数组。
   *
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  vertex_indices() const;

  /**
   * 返回构成有关参考单元格的线的数量。一条线是参考单元格的一个
   * "边"（一个一维对象）。
   *
   */
  unsigned int
  n_lines() const;

  /**
   * 返回一个对象，可以看作是一个包含从零到n_lines()所有索引的数组。
   *
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  line_indices() const;

  /**
   * 返回构成相关参考单元的面的数量。一个面是包围参考单元的`(dim-1)`维对象。
   *
   */
  unsigned int
  n_faces() const;

  /**
   * 返回一个对象，它可以被认为是一个包含从零到n_faces()所有索引的数组。
   *
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  face_indices() const;

  /**
   * 返回当前对象的面 @p face_no
   * 的参考单元格类型。例如，如果当前对象是
   * ReferenceCells::Tetrahedron, ，那么`face_no`必须介于 $[0,4)$
   * 之间，函数将总是返回 ReferenceCells::Triangle.
   * 如果当前对象是 ReferenceCells::Hexahedron,
   * ，那么`face_no`必须介于 $[0,6)$ 之间，函数将总是返回
   * ReferenceCells::Quadrilateral.
   * 对于楔形和金字塔，返回的对象可能是
   * ReferenceCells::Triangle  或 ReferenceCells::Quadrilateral,
   * ，取决于给定的索引。
   *
   */
  ReferenceCell
  face_reference_cell(const unsigned int face_no) const;

  /**
   * @}
   *
   */

  /**
   * @name  单元内和面上的对象之间的关系  
     * @{ 
   *
   */

  /**
   * 返回哪些子单元与母单元的某个面相邻。
   * 例如，在2D中，一个四边形单元的布局如下。
   * @verbatim
   *    3
   * 2-->--3
   * |     |
   * 0 ^     ^ 1
   * |     |
   * 0-->--1
   *    2
   * @endverbatim
   * 顶点和面用其数字表示，面也用其方向表示。
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
   * ---*---*
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
   * @endverbatim
   * 因此，面0上的子单元是（按面的方向排序）0和2，面3上是2和3，等等。    对于三个空间维度，子单元的确切顺序在这个类的一般文档中规定。    <tt>面的方向</tt>参数目前只针对四边形和六面体。它决定了这个函数如何处理以标准和非标准方向为导向的面。它代表了整个<tt>面的方向</tt>、<tt>面的翻转</tt>和<tt>面的旋转</tt>的位代码，并默认为标准方向。面部方向的概念在这个 @ref GlossFaceOrientation "词汇表 "
   * 条目中得到了解释。
   *
   */
  unsigned int
  child_cell_on_face(const unsigned int  face_n,
                     const unsigned int  subface_n,
                     const unsigned char face_orientation = 1) const;

  /**
   * 对于一个单元格中的给定顶点，返回一对面的索引和该面的顶点索引。
   * @note
   * 在实践中，一个顶点当然通常属于一个以上的面，我们可以返回不同的面和其中相应的索引。这个函数选择哪个面通常并不重要（这个函数也不会故意暴露）。
   *
   */
  std::array<unsigned int, 2>
  standard_vertex_to_face_and_vertex_index(const unsigned int vertex) const;

  /**
   * 对于一个单元格中的某一行，返回一对面的索引和该面中的行的索引。
   * @note
   * 在实践中，一条线通常是一个以上的面的一部分，我们可以返回不同的面和其中相应的索引。这个函数选择哪个面通常并不重要（这个函数也不会故意暴露）。
   *
   */
  std::array<unsigned int, 2>
  standard_line_to_face_and_line_index(const unsigned int line) const;

  /**
   * 将面的行数映射到单元格的行数。
   *
   */
  unsigned int
  face_to_cell_lines(const unsigned int  face,
                     const unsigned int  line,
                     const unsigned char face_orientation) const;

  /**
   * 将面的顶点编号映射到单元的顶点编号。
   *
   */
  unsigned int
  face_to_cell_vertices(const unsigned int  face,
                        const unsigned int  vertex,
                        const unsigned char face_orientation) const;

  /**
   * 根据面的方向，纠正顶点索引。
   *
   */
  unsigned int
  standard_to_real_face_vertex(const unsigned int  vertex,
                               const unsigned int  face,
                               const unsigned char face_orientation) const;

  /**
   * 根据面的方向，纠正线的索引。
   *
   */
  unsigned int
  standard_to_real_face_line(const unsigned int  line,
                             const unsigned int  face,
                             const unsigned char face_orientation) const;

  /**
   * 返回索引为 @p line
   * 的线在单元格内是否为标准方向，给定当前单元格内面的
   * @p face_orientation ，以及该面内线的 @p line_orientation
   * 标志。  @p true
   * 表示线的方向是从顶点0到顶点1，否则就是相反的方向。在1d和2d中，这总是
   * @p true,
   * ，但在3d中，它可能是不同的，见GeometryInfo类的文档中的相应讨论。
   *
   */
  bool
  standard_vs_true_line_orientation(const unsigned int  line,
                                    const unsigned char face_orientation,
                                    const unsigned char line_orientation) const;

  /**
   * @}
   *
   */

  /**
   * @name  参考单元的几何属性  @name  查询一个参考单元的构件数量  @{
   *
   */

  /* 返回  $i$  参考单元的一个面的第-个单位切向量。  这些矢量的排列方式是，两个矢量之间的交积返回单位法向量。      @pre   $i$  必须在零和`dim-1`之间。 
*
*/
  template <int dim>
  Tensor<1, dim>
  unit_tangential_vectors(const unsigned int face_no,
                          const unsigned int i) const;

  /**
   * 返回参考单元格的一个面的单位法向量。
   *
   */
  template <int dim>
  Tensor<1, dim>
  unit_normal_vectors(const unsigned int face_no) const;

  /**
   * 确定由其顶点 @p var_1 描述的当前实体相对于由 @p var_0.
   * 描述的实体的方向。
   *
   */
  template <typename T, std::size_t N>
  unsigned char
  compute_orientation(const std::array<T, N> &vertices_0,
                      const std::array<T, N> &vertices_1) const;

  /**
   * compute_orientation()的逆向函数。
   *
   */
  template <typename T, std::size_t N>
  std::array<T, N>
  permute_according_orientation(const std::array<T, N> &vertices,
                                const unsigned int      orientation) const;

  /**
   * 返回一个给定的 @p vertex_index 所属的面的向量。
   *
   */
  ArrayView<const unsigned int>
  faces_for_given_vertex(const unsigned int vertex_index) const;

  /**
   * @}
   *
   */

  /**
   * @name  在deal.II索引和其他程序使用的格式之间进行转换  
     * @{ 
   *
   */

  /**
   * 将ExodusII的顶点编号映射到deal.II的顶点编号。
   *
   */
  unsigned int
  exodusii_vertex_to_deal_vertex(const unsigned int vertex_n) const;

  /**
   * 将一个ExodusII的面数映射到一个deal.II的面数。
   *
   */
  unsigned int
  exodusii_face_to_deal_face(const unsigned int face_n) const;

  /**
   * 将一个UNV顶点编号映射到一个deal.II顶点编号。
   *
   */
  unsigned int
  unv_vertex_to_deal_vertex(const unsigned int vertex_n) const;

  /**
   * 返回一个VTK线性形状常数，与参考单元对应。
   *
   */
  unsigned int
  vtk_linear_type() const;

  /**
   * 返回对应于参考单元的VTK二次方形状常数。
   *
   */
  unsigned int
  vtk_quadratic_type() const;

  /**
   * 返回对应于参考单元格的VTK拉格朗日形状常数。
   *
   */
  unsigned int
  vtk_lagrange_type() const;

  /**
   * @}
   *
   */

  /**
   * @name  其他函数  @{
   *
   */

  /**
   * 返回当前对象所代表的参考单元格的文本表示。
   *
   */
  std::string
  to_string() const;

  /**
   * 对整数的转换操作。
   *
   */
  constexpr operator std::uint8_t() const;

  /**
   * 用于平等比较的运算符。
   *
   */
  constexpr bool
  operator==(const ReferenceCell &type) const;

  /**
   * 用于不等式比较的运算符。
   *
   */
  constexpr bool
  operator!=(const ReferenceCell &type) const;

  /**
   * 为了使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)进行序列化，从一个流中写入和读取此对象的数据。
   *
   */
  template <class Archive>
  void
  serialize(Archive &archive, const unsigned int  /*version*/ );

  /**
   * @}
   *
   */

private:
  /**
   * 存储这个对象实际对应的变量。
   *
   */
  std::uint8_t kind;

  /**
   * 构造函数。这是用来创建这个类的不同的`静态`成员变量的构造函数。它是
   * "私有 "的，但是可以被作为该类 "朋友
   * "的内部命名空间中的函数调用。
   *
   */
  constexpr ReferenceCell(const std::uint8_t kind);

  /**
   * 一种构造函数
   *
   * - 不完全是私有的，因为它可以被任何人调用，但至少是隐藏在一个内部命名空间中。
   *
   */
  friend DEAL_II_CONSTEXPR ReferenceCell
                           internal::ReferenceCell::make_reference_cell_from_int(const std::uint8_t);
};



inline constexpr ReferenceCell::ReferenceCell(const std::uint8_t kind)
  : kind(kind)
{}



inline constexpr ReferenceCell::operator std::uint8_t() const
{
  return kind;
}



inline constexpr bool
ReferenceCell::operator==(const ReferenceCell &type) const
{
  return kind == type.kind;
}



inline constexpr bool
ReferenceCell::operator!=(const ReferenceCell &type) const
{
  return kind != type.kind;
}



namespace internal
{
  namespace ReferenceCell
  {
    inline DEAL_II_CONSTEXPR dealii::ReferenceCell
                             make_reference_cell_from_int(const std::uint8_t kind)
    {
      // Make sure these are the only indices from which objects can be
      // created.
      Assert((kind == static_cast<std::uint8_t>(-1)) || (kind < 8),
             ExcInternalError());

      // Call the private constructor, which we can from here because this
      // function is a 'friend'.
      return {kind};
    }
  } // namespace ReferenceCell
} // namespace internal



/**
 * 一个命名空间，我们在其中定义对应于特定参考单元的对象。这里定义的对象是对所有可能的参考单元的完整列举，可以在deal.II中使用。
 * @relates  ReferenceCell
 *
 *
 */
namespace ReferenceCells
{
  DEAL_II_CONSTEXPR const ReferenceCell Vertex =
    internal::ReferenceCell::make_reference_cell_from_int(0);
  DEAL_II_CONSTEXPR const ReferenceCell Line =
    internal::ReferenceCell::make_reference_cell_from_int(1);
  DEAL_II_CONSTEXPR const ReferenceCell Triangle =
    internal::ReferenceCell::make_reference_cell_from_int(2);
  DEAL_II_CONSTEXPR const ReferenceCell Quadrilateral =
    internal::ReferenceCell::make_reference_cell_from_int(3);
  DEAL_II_CONSTEXPR const ReferenceCell Tetrahedron =
    internal::ReferenceCell::make_reference_cell_from_int(4);
  DEAL_II_CONSTEXPR const ReferenceCell Pyramid =
    internal::ReferenceCell::make_reference_cell_from_int(5);
  DEAL_II_CONSTEXPR const ReferenceCell Wedge =
    internal::ReferenceCell::make_reference_cell_from_int(6);
  DEAL_II_CONSTEXPR const ReferenceCell Hexahedron =
    internal::ReferenceCell::make_reference_cell_from_int(7);
  DEAL_II_CONSTEXPR const ReferenceCell Invalid =
    internal::ReferenceCell::make_reference_cell_from_int(
      static_cast<std::uint8_t>(-1));

  /**
   * 返回给定维度`dim`的正确单线参考单元类型。根据模板参数`dim`，该函数返回对Vertex、Triangle或Tetrahedron的引用。
   *
   */
  template <int dim>
  constexpr const ReferenceCell &
  get_simplex();

  /**
   * 返回给定维度`dim`的正确超立方体参考单元类型。根据模板参数`dim`，该函数返回对顶点、四边形或六面体的引用。
   *
   */
  template <int dim>
  constexpr const ReferenceCell &
  get_hypercube();
} // namespace ReferenceCells



inline DEAL_II_CONSTEXPR
ReferenceCell::ReferenceCell()
  : ReferenceCell(ReferenceCells::Invalid)
{}



template <class Archive>
inline void
ReferenceCell::serialize(Archive &archive, const unsigned int  /*version*/ )
{
  archive &kind;
}



inline ArrayView<const unsigned int>
ReferenceCell::faces_for_given_vertex(const unsigned int vertex) const
{
  if (*this == ReferenceCells::Line)
    {
      AssertIndexRange(vertex, GeometryInfo<1>::vertices_per_cell);
      return {&GeometryInfo<2>::vertex_to_face[vertex][0], 1};
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      AssertIndexRange(vertex, GeometryInfo<2>::vertices_per_cell);
      return {&GeometryInfo<2>::vertex_to_face[vertex][0], 2};
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      AssertIndexRange(vertex, GeometryInfo<3>::vertices_per_cell);
      return {&GeometryInfo<3>::vertex_to_face[vertex][0], 3};
    }
  else if (*this == ReferenceCells::Triangle)
    {
      AssertIndexRange(vertex, 3);
      static const ndarray<unsigned int, 3, 2> table = {
        {{{0, 2}}, {{0, 1}}, {{1, 2}}}};

      return table[vertex];
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      AssertIndexRange(vertex, 4);
      static const ndarray<unsigned int, 4, 3> table = {
        {{{0, 1, 2}}, {{0, 1, 3}}, {{0, 2, 3}}, {{1, 2, 3}}}};

      return table[vertex];
    }
  else if (*this == ReferenceCells::Wedge)
    {
      AssertIndexRange(vertex, 6);
      static const ndarray<unsigned int, 6, 3> table = {{{{0, 2, 4}},
                                                         {{0, 2, 3}},
                                                         {{0, 3, 4}},
                                                         {{1, 2, 4}},
                                                         {{1, 2, 3}},
                                                         {{1, 3, 4}}}};

      return table[vertex];
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      AssertIndexRange(vertex, 5);
      static const unsigned int X = numbers::invalid_unsigned_int;
      static const ndarray<unsigned int, 5, 4> table = {{{{0, 1, 3, X}},
                                                         {{0, 2, 3, X}},
                                                         {{0, 1, 4, X}},
                                                         {{0, 2, 4, X}},
                                                         {{1, 2, 3, 4}}}};

      return {&table[vertex][0], vertex == 4 ? 4u : 3u};
    }

  Assert(false, ExcNotImplemented());

  return {};
}



inline bool
ReferenceCell::is_hyper_cube() const
{
  return (*this == ReferenceCells::Vertex || *this == ReferenceCells::Line ||
          *this == ReferenceCells::Quadrilateral ||
          *this == ReferenceCells::Hexahedron);
}



inline bool
ReferenceCell::is_simplex() const
{
  return (*this == ReferenceCells::Vertex || *this == ReferenceCells::Line ||
          *this == ReferenceCells::Triangle ||
          *this == ReferenceCells::Tetrahedron);
}



inline unsigned int
ReferenceCell::get_dimension() const
{
  if (*this == ReferenceCells::Vertex)
    return 0;
  else if (*this == ReferenceCells::Line)
    return 1;
  else if ((*this == ReferenceCells::Triangle) ||
           (*this == ReferenceCells::Quadrilateral))
    return 2;
  else if ((*this == ReferenceCells::Tetrahedron) ||
           (*this == ReferenceCells::Pyramid) ||
           (*this == ReferenceCells::Wedge) ||
           (*this == ReferenceCells::Hexahedron))
    return 3;

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



inline unsigned int
ReferenceCell::n_vertices() const
{
  if (*this == ReferenceCells::Vertex)
    return 1;
  else if (*this == ReferenceCells::Line)
    return 2;
  else if (*this == ReferenceCells::Triangle)
    return 3;
  else if (*this == ReferenceCells::Quadrilateral)
    return 4;
  else if (*this == ReferenceCells::Tetrahedron)
    return 4;
  else if (*this == ReferenceCells::Pyramid)
    return 5;
  else if (*this == ReferenceCells::Wedge)
    return 6;
  else if (*this == ReferenceCells::Hexahedron)
    return 8;

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



inline unsigned int
ReferenceCell::n_lines() const
{
  if (*this == ReferenceCells::Vertex)
    return 0;
  else if (*this == ReferenceCells::Line)
    return 1;
  else if (*this == ReferenceCells::Triangle)
    return 3;
  else if (*this == ReferenceCells::Quadrilateral)
    return 4;
  else if (*this == ReferenceCells::Tetrahedron)
    return 6;
  else if (*this == ReferenceCells::Pyramid)
    return 7;
  else if (*this == ReferenceCells::Wedge)
    return 9;
  else if (*this == ReferenceCells::Hexahedron)
    return 12;

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



inline unsigned int
ReferenceCell::n_faces() const
{
  if (*this == ReferenceCells::Vertex)
    return 0;
  else if (*this == ReferenceCells::Line)
    return 2;
  else if (*this == ReferenceCells::Triangle)
    return 3;
  else if (*this == ReferenceCells::Quadrilateral)
    return 4;
  else if (*this == ReferenceCells::Tetrahedron)
    return 4;
  else if (*this == ReferenceCells::Pyramid)
    return 5;
  else if (*this == ReferenceCells::Wedge)
    return 5;
  else if (*this == ReferenceCells::Hexahedron)
    return 6;

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
ReferenceCell::vertex_indices() const
{
  return {0U, n_vertices()};
}



inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
ReferenceCell::line_indices() const
{
  return {0U, n_lines()};
}



inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
ReferenceCell::face_indices() const
{
  return {0U, n_faces()};
}



inline ReferenceCell
ReferenceCell::face_reference_cell(const unsigned int face_no) const
{
  AssertIndexRange(face_no, n_faces());

  if (*this == ReferenceCells::Vertex)
    return ReferenceCells::Invalid;
  else if (*this == ReferenceCells::Line)
    return ReferenceCells::Vertex;
  else if (*this == ReferenceCells::Triangle)
    return ReferenceCells::Line;
  else if (*this == ReferenceCells::Quadrilateral)
    return ReferenceCells::Line;
  else if (*this == ReferenceCells::Tetrahedron)
    return ReferenceCells::Triangle;
  else if (*this == ReferenceCells::Pyramid)
    {
      if (face_no == 0)
        return ReferenceCells::Quadrilateral;
      else
        return ReferenceCells::Triangle;
    }
  else if (*this == ReferenceCells::Wedge)
    {
      if (face_no > 1)
        return ReferenceCells::Quadrilateral;
      else
        return ReferenceCells::Triangle;
    }
  else if (*this == ReferenceCells::Hexahedron)
    return ReferenceCells::Quadrilateral;

  Assert(false, ExcNotImplemented());
  return ReferenceCells::Invalid;
}



inline unsigned int
ReferenceCell::child_cell_on_face(
  const unsigned int  face,
  const unsigned int  subface,
  const unsigned char face_orientation_raw) const
{
  AssertIndexRange(face, n_faces());

  if (*this == ReferenceCells::Vertex)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Line)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Triangle)
    {
      static const ndarray<unsigned int, 3, 2> subcells = {
        {{{0, 1}}, {{1, 2}}, {{2, 0}}}};

      return subcells[face][subface];
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      const bool face_orientation = Utilities::get_bit(face_orientation_raw, 0);
      const bool face_flip        = Utilities::get_bit(face_orientation_raw, 2);
      const bool face_rotation    = Utilities::get_bit(face_orientation_raw, 1);

      return GeometryInfo<2>::child_cell_on_face(
        RefinementCase<2>(RefinementPossibilities<2>::no_refinement),
        face,
        subface,
        face_orientation,
        face_flip,
        face_rotation);
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Wedge)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      const bool face_orientation = Utilities::get_bit(face_orientation_raw, 0);
      const bool face_flip        = Utilities::get_bit(face_orientation_raw, 2);
      const bool face_rotation    = Utilities::get_bit(face_orientation_raw, 1);

      return GeometryInfo<3>::child_cell_on_face(
        RefinementCase<3>(RefinementPossibilities<3>::no_refinement),
        face,
        subface,
        face_orientation,
        face_flip,
        face_rotation);
    }

  Assert(false, ExcNotImplemented());
  return {};
}



inline std::array<unsigned int, 2>
ReferenceCell::standard_vertex_to_face_and_vertex_index(
  const unsigned int vertex) const
{
  AssertIndexRange(vertex, n_vertices());

  if (*this == ReferenceCells::Vertex)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Line)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Triangle)
    {
      static const ndarray<unsigned int, 3, 2> table = {
        {{{0, 0}}, {{0, 1}}, {{1, 1}}}};

      return table[vertex];
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      return GeometryInfo<2>::standard_quad_vertex_to_line_vertex_index(vertex);
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      static const ndarray<unsigned int, 4, 2> table = {
        {{{0, 0}}, {{0, 1}}, {{0, 2}}, {{1, 2}}}};

      return table[vertex];
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      static const ndarray<unsigned int, 5, 2> table = {
        {{{0, 0}}, {{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 2}}}};

      return table[vertex];
    }
  else if (*this == ReferenceCells::Wedge)
    {
      static const ndarray<unsigned int, 6, 2> table = {
        {{{0, 1}}, {{0, 0}}, {{0, 2}}, {{1, 0}}, {{1, 1}}, {{1, 2}}}};

      return table[vertex];
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      return GeometryInfo<3>::standard_hex_vertex_to_quad_vertex_index(vertex);
    }

  Assert(false, ExcNotImplemented());
  return {};
}



inline std::array<unsigned int, 2>
ReferenceCell::standard_line_to_face_and_line_index(
  const unsigned int line) const
{
  AssertIndexRange(line, n_lines());

  if (*this == ReferenceCells::Vertex)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Line)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Triangle)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      static const std::array<unsigned int, 2> table[6] = {
        {{0, 0}}, {{0, 1}}, {{0, 2}}, {{1, 1}}, {{1, 2}}, {{2, 1}}};

      return table[line];
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      static const std::array<unsigned int, 2> table[8] = {{{0, 0}},
                                                           {{0, 1}},
                                                           {{0, 2}},
                                                           {{0, 3}},
                                                           {{1, 2}},
                                                           {{2, 1}},
                                                           {{1, 1}},
                                                           {{2, 2}}};

      return table[line];
    }
  else if (*this == ReferenceCells::Wedge)
    {
      static const std::array<unsigned int, 2> table[9] = {{{0, 0}},
                                                           {{0, 2}},
                                                           {{0, 1}},
                                                           {{1, 0}},
                                                           {{1, 1}},
                                                           {{1, 2}},
                                                           {{2, 0}},
                                                           {{2, 1}},
                                                           {{3, 1}}};

      return table[line];
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      return GeometryInfo<3>::standard_hex_line_to_quad_line_index(line);
    }

  Assert(false, ExcNotImplemented());
  return {};
}



inline unsigned int
ReferenceCell::face_to_cell_lines(const unsigned int  face,
                                  const unsigned int  line,
                                  const unsigned char face_orientation) const
{
  AssertIndexRange(face, n_faces());
  AssertIndexRange(line, face_reference_cell(face).n_lines());

  if (*this == ReferenceCells::Vertex)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Line)
    {
      return GeometryInfo<1>::face_to_cell_lines(
        face,
        line,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }
  else if (*this == ReferenceCells::Triangle)
    {
      return face;
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      return GeometryInfo<2>::face_to_cell_lines(
        face,
        line,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      const static ndarray<unsigned int, 4, 3> table = {
        {{{0, 1, 2}}, {{0, 3, 4}}, {{2, 5, 3}}, {{1, 4, 5}}}};

      return table[face]
                  [standard_to_real_face_line(line, face, face_orientation)];
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Wedge)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      return GeometryInfo<3>::face_to_cell_lines(
        face,
        line,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



inline unsigned int
ReferenceCell::face_to_cell_vertices(const unsigned int  face,
                                     const unsigned int  vertex,
                                     const unsigned char face_orientation) const
{
  AssertIndexRange(face, n_faces());
  AssertIndexRange(vertex, face_reference_cell(face).n_vertices());

  if (*this == ReferenceCells::Vertex)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Line)
    {
      return GeometryInfo<1>::face_to_cell_vertices(
        face,
        vertex,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }
  else if (*this == ReferenceCells::Triangle)
    {
      static const ndarray<unsigned int, 3, 2> table = {
        {{{0, 1}}, {{1, 2}}, {{2, 0}}}};

      return table[face][face_orientation ? vertex : (1 - vertex)];
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      return GeometryInfo<2>::face_to_cell_vertices(
        face,
        vertex,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      static const ndarray<unsigned int, 4, 3> table = {
        {{{0, 1, 2}}, {{1, 0, 3}}, {{0, 2, 3}}, {{2, 1, 3}}}};

      return table[face][standard_to_real_face_vertex(
        vertex, face, face_orientation)];
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      constexpr auto X = numbers::invalid_unsigned_int;
      static const ndarray<unsigned int, 5, 4> table = {{{{0, 1, 2, 3}},
                                                         {{0, 2, 4, X}},
                                                         {{3, 1, 4, X}},
                                                         {{1, 0, 4, X}},
                                                         {{2, 3, 4, X}}}};

      return table[face][standard_to_real_face_vertex(
        vertex, face, face_orientation)];
    }
  else if (*this == ReferenceCells::Wedge)
    {
      constexpr auto X = numbers::invalid_unsigned_int;
      static const ndarray<unsigned int, 6, 4> table = {{{{1, 0, 2, X}},
                                                         {{3, 4, 5, X}},
                                                         {{0, 1, 3, 4}},
                                                         {{1, 2, 4, 5}},
                                                         {{2, 0, 5, 3}}}};

      return table[face][standard_to_real_face_vertex(
        vertex, face, face_orientation)];
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      return GeometryInfo<3>::face_to_cell_vertices(
        face,
        vertex,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



inline unsigned int
ReferenceCell::standard_to_real_face_vertex(
  const unsigned int  vertex,
  const unsigned int  face,
  const unsigned char face_orientation) const
{
  AssertIndexRange(face, n_faces());
  AssertIndexRange(vertex, face_reference_cell(face).n_vertices());

  if (*this == ReferenceCells::Vertex)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Line)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Triangle)
    {
      static const ndarray<unsigned int, 2, 2> table = {{{{1, 0}}, {{0, 1}}}};

      return table[face_orientation][vertex];
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      return GeometryInfo<2>::standard_to_real_line_vertex(vertex,
                                                           face_orientation);
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      static const ndarray<unsigned int, 6, 3> table = {{{{0, 2, 1}},
                                                         {{0, 1, 2}},
                                                         {{2, 1, 0}},
                                                         {{1, 2, 0}},
                                                         {{1, 0, 2}},
                                                         {{2, 0, 1}}}};

      return table[face_orientation][vertex];
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      if (face == 0) // The quadrilateral face
        {
          return GeometryInfo<3>::standard_to_real_face_vertex(
            vertex,
            Utilities::get_bit(face_orientation, 0),
            Utilities::get_bit(face_orientation, 2),
            Utilities::get_bit(face_orientation, 1));
        }
      else // One of the triangular faces
        {
          static const ndarray<unsigned int, 6, 3> table = {{{{0, 2, 1}},
                                                             {{0, 1, 2}},
                                                             {{2, 1, 0}},
                                                             {{1, 2, 0}},
                                                             {{1, 0, 2}},
                                                             {{2, 0, 1}}}};

          return table[face_orientation][vertex];
        }
    }
  else if (*this == ReferenceCells::Wedge)
    {
      if (face > 1) // One of the quadrilateral faces
        {
          return GeometryInfo<3>::standard_to_real_face_vertex(
            vertex,
            Utilities::get_bit(face_orientation, 0),
            Utilities::get_bit(face_orientation, 2),
            Utilities::get_bit(face_orientation, 1));
        }
      else // One of the triangular faces
        {
          static const ndarray<unsigned int, 6, 3> table = {{{{0, 2, 1}},
                                                             {{0, 1, 2}},
                                                             {{2, 1, 0}},
                                                             {{1, 2, 0}},
                                                             {{1, 0, 2}},
                                                             {{2, 0, 1}}}};

          return table[face_orientation][vertex];
        }
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      return GeometryInfo<3>::standard_to_real_face_vertex(
        vertex,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



inline unsigned int
ReferenceCell::standard_to_real_face_line(
  const unsigned int  line,
  const unsigned int  face,
  const unsigned char face_orientation) const
{
  AssertIndexRange(face, n_faces());
  AssertIndexRange(line, face_reference_cell(face).n_lines());

  if (*this == ReferenceCells::Vertex)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Line)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Triangle)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      Assert(false, ExcNotImplemented());
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      static const ndarray<unsigned int, 6, 3> table = {{{{2, 1, 0}},
                                                         {{0, 1, 2}},
                                                         {{1, 0, 2}},
                                                         {{1, 2, 0}},
                                                         {{0, 2, 1}},
                                                         {{2, 0, 1}}}};

      return table[face_orientation][line];
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      if (face == 0) // The quadrilateral face
        {
          return GeometryInfo<3>::standard_to_real_face_line(
            line,
            Utilities::get_bit(face_orientation, 0),
            Utilities::get_bit(face_orientation, 2),
            Utilities::get_bit(face_orientation, 1));
        }
      else // One of the triangular faces
        {
          static const ndarray<unsigned int, 6, 3> table = {{{{2, 1, 0}},
                                                             {{0, 1, 2}},
                                                             {{1, 0, 2}},
                                                             {{1, 2, 0}},
                                                             {{0, 2, 1}},
                                                             {{2, 0, 1}}}};

          return table[face_orientation][line];
        }
    }
  else if (*this == ReferenceCells::Wedge)
    {
      if (face > 1) // One of the quadrilateral faces
        {
          return GeometryInfo<3>::standard_to_real_face_line(
            line,
            Utilities::get_bit(face_orientation, 0),
            Utilities::get_bit(face_orientation, 2),
            Utilities::get_bit(face_orientation, 1));
        }
      else // One of the triangular faces
        {
          static const ndarray<unsigned int, 6, 3> table = {{{{2, 1, 0}},
                                                             {{0, 1, 2}},
                                                             {{1, 0, 2}},
                                                             {{1, 2, 0}},
                                                             {{0, 2, 1}},
                                                             {{2, 0, 1}}}};

          return table[face_orientation][line];
        }
    }
  else if (*this == ReferenceCells::Hexahedron)
    {
      return GeometryInfo<3>::standard_to_real_face_line(
        line,
        Utilities::get_bit(face_orientation, 0),
        Utilities::get_bit(face_orientation, 2),
        Utilities::get_bit(face_orientation, 1));
    }

  Assert(false, ExcNotImplemented());
  return numbers::invalid_unsigned_int;
}



namespace ReferenceCells
{
  template <int dim>
  inline constexpr const ReferenceCell &
  get_simplex()
  {
    switch (dim)
      {
        case 0:
          return ReferenceCells::Vertex;
        case 1:
          return ReferenceCells::Line;
        case 2:
          return ReferenceCells::Triangle;
        case 3:
          return ReferenceCells::Tetrahedron;
        default:
          Assert(false, ExcNotImplemented());
          return ReferenceCells::Invalid;
      }
  }



  template <int dim>
  inline constexpr const ReferenceCell &
  get_hypercube()
  {
    switch (dim)
      {
        case 0:
          return ReferenceCells::Vertex;
        case 1:
          return ReferenceCells::Line;
        case 2:
          return ReferenceCells::Quadrilateral;
        case 3:
          return ReferenceCells::Hexahedron;
        default:
          Assert(false, ExcNotImplemented());
          return ReferenceCells::Invalid;
      }
  }
} // namespace ReferenceCells


inline ReferenceCell
ReferenceCell::n_vertices_to_type(const int dim, const unsigned int n_vertices)
{
  AssertIndexRange(dim, 4);
  AssertIndexRange(n_vertices, 9);

  const auto                                X     = ReferenceCells::Invalid;
  static const ndarray<ReferenceCell, 4, 9> table = {
    {// dim 0
     {{X, ReferenceCells::Vertex, X, X, X, X, X, X, X}},
     // dim 1
     {{X, X, ReferenceCells::Line, X, X, X, X, X, X}},
     // dim 2
     {{X,
       X,
       X,
       ReferenceCells::Triangle,
       ReferenceCells::Quadrilateral,
       X,
       X,
       X,
       X}},
     // dim 3
     {{X,
       X,
       X,
       X,
       ReferenceCells::Tetrahedron,
       ReferenceCells::Pyramid,
       ReferenceCells::Wedge,
       X,
       ReferenceCells::Hexahedron}}}};
  Assert(table[dim][n_vertices] != ReferenceCells::Invalid,
         ExcMessage("The combination of dim = " + std::to_string(dim) +
                    " and n_vertices = " + std::to_string(n_vertices) +
                    " does not correspond to a known reference cell type."));
  return table[dim][n_vertices];
}



template <int dim>
inline double
ReferenceCell::d_linear_shape_function(const Point<dim> & xi,
                                       const unsigned int i) const
{
  AssertDimension(dim, get_dimension());
  if (*this == ReferenceCells::get_hypercube<dim>())
    return GeometryInfo<dim>::d_linear_shape_function(xi, i);

  if (*this ==
      ReferenceCells::Triangle) // see also
                                // BarycentricPolynomials<2>::compute_value
    {
      switch (i)
        {
          case 0:
            return 1.0 - xi[std::min(0, dim - 1)] - xi[std::min(1, dim - 1)];
          case 1:
            return xi[std::min(0, dim - 1)];
          case 2:
            return xi[std::min(1, dim - 1)];
        }
    }

  if (*this ==
      ReferenceCells::Tetrahedron) // see also
                                   // BarycentricPolynomials<3>::compute_value
    {
      switch (i)
        {
          case 0:
            return 1.0 - xi[std::min(0, dim - 1)] - xi[std::min(1, dim - 1)] -
                   xi[std::min(2, dim - 1)];
          case 1:
            return xi[std::min(0, dim - 1)];
          case 2:
            return xi[std::min(1, dim - 1)];
          case 3:
            return xi[std::min(2, dim - 1)];
        }
    }

  if (*this ==
      ReferenceCells::Wedge) // see also
                             // ScalarLagrangePolynomialWedge::compute_value
    {
      return ReferenceCell(ReferenceCells::Triangle)
               .d_linear_shape_function<2>(Point<2>(xi[std::min(0, dim - 1)],
                                                    xi[std::min(1, dim - 1)]),
                                           i % 3) *
             ReferenceCell(ReferenceCells::Line)
               .d_linear_shape_function<1>(Point<1>(xi[std::min(2, dim - 1)]),
                                           i / 3);
    }

  if (*this ==
      ReferenceCells::Pyramid) // see also
                               // ScalarLagrangePolynomialPyramid::compute_value
    {
      const double Q14 = 0.25;
      double       ration;

      const double r = xi[std::min(0, dim - 1)];
      const double s = xi[std::min(1, dim - 1)];
      const double t = xi[std::min(2, dim - 1)];

      if (fabs(t - 1.0) > 1.0e-14)
        {
          ration = (r * s * t) / (1.0 - t);
        }
      else
        {
          ration = 0.0;
        }

      if (i == 0)
        return Q14 * ((1.0 - r) * (1.0 - s) - t + ration);
      if (i == 1)
        return Q14 * ((1.0 + r) * (1.0 - s) - t - ration);
      if (i == 2)
        return Q14 * ((1.0 - r) * (1.0 + s) - t - ration);
      if (i == 3)
        return Q14 * ((1.0 + r) * (1.0 + s) - t + ration);
      else
        return t;
    }

  Assert(false, ExcNotImplemented());

  return 0.0;
}



template <int dim>
inline Tensor<1, dim>
ReferenceCell::d_linear_shape_function_gradient(const Point<dim> & xi,
                                                const unsigned int i) const
{
  AssertDimension(dim, get_dimension());
  if (*this == ReferenceCells::get_hypercube<dim>())
    return GeometryInfo<dim>::d_linear_shape_function_gradient(xi, i);

  if (*this ==
      ReferenceCells::Triangle) // see also
                                // BarycentricPolynomials<2>::compute_grad
    {
      switch (i)
        {
          case 0:
            return Point<dim>(-1.0, -1.0);
          case 1:
            return Point<dim>(+1.0, +0.0);
          case 2:
            return Point<dim>(+0.0, +1.0);
        }
    }

  Assert(false, ExcNotImplemented());

  return Point<dim>(+0.0, +0.0, +0.0);
}


template <int dim>
inline Tensor<1, dim>
ReferenceCell::unit_tangential_vectors(const unsigned int face_no,
                                       const unsigned int i) const
{
  AssertDimension(dim, get_dimension());
  AssertIndexRange(i, dim - 1);

  if (*this == ReferenceCells::get_hypercube<dim>())
    {
      AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);
      return GeometryInfo<dim>::unit_tangential_vectors[face_no][i];
    }
  else if (*this == ReferenceCells::Triangle)
    {
      AssertIndexRange(face_no, 3);
      static const std::array<Tensor<1, dim>, 3> table = {
        {Point<dim>(1, 0),
         Point<dim>(-std::sqrt(0.5), +std::sqrt(0.5)),
         Point<dim>(0, -1)}};

      return table[face_no];
    }
  else if (*this == ReferenceCells::Tetrahedron)
    {
      AssertIndexRange(face_no, 4);
      static const ndarray<Tensor<1, dim>, 4, 2> table = {
        {{{Point<dim>(0, 1, 0), Point<dim>(1, 0, 0)}},
         {{Point<dim>(1, 0, 0), Point<dim>(0, 0, 1)}},
         {{Point<dim>(0, 0, 1), Point<dim>(0, 1, 0)}},
         {{Point<dim>(-std::pow(1.0 / 3.0, 1.0 / 4.0),
                      +std::pow(1.0 / 3.0, 1.0 / 4.0),
                      0),
           Point<dim>(-std::pow(1.0 / 3.0, 1.0 / 4.0),
                      0,
                      +std::pow(1.0 / 3.0, 1.0 / 4.0))}}}};

      return table[face_no][i];
    }
  else if (*this == ReferenceCells::Wedge)
    {
      AssertIndexRange(face_no, 5);
      static const ndarray<Tensor<1, dim>, 5, 2> table = {
        {{{Point<dim>(0, 1, 0), Point<dim>(1, 0, 0)}},
         {{Point<dim>(1, 0, 0), Point<dim>(0, 1, 0)}},
         {{Point<dim>(1, 0, 0), Point<dim>(0, 0, 1)}},
         {{Point<dim>(-1 / std::sqrt(2.0), +1 / std::sqrt(2.0), 0),
           Point<dim>(0, 0, 1)}},
         {{Point<dim>(0, 0, 1), Point<dim>(0, 1, 0)}}}};

      return table[face_no][i];
    }
  else if (*this == ReferenceCells::Pyramid)
    {
      AssertIndexRange(face_no, 5);
      static const ndarray<Tensor<1, dim>, 5, 2> table = {
        {{{Point<dim>(0, 1, 0), Point<dim>(1, 0, 0)}},
         {{Point<dim>(+1.0 / sqrt(2.0), 0, +1.0 / sqrt(2.0)),
           Point<dim>(0, 1, 0)}},
         {{Point<dim>(+1.0 / sqrt(2.0), 0, -1.0 / sqrt(2.0)),
           Point<dim>(0, 1, 0)}},
         {{Point<dim>(1, 0, 0),
           Point<dim>(0, +1.0 / sqrt(2.0), +1.0 / sqrt(2.0))}},
         {{Point<dim>(1, 0, 0),
           Point<dim>(0, +1.0 / sqrt(2.0), -1.0 / sqrt(2.0))}}}};

      return table[face_no][i];
    }

  Assert(false, ExcNotImplemented());

  return {};
}



template <int dim>
inline Tensor<1, dim>
ReferenceCell::unit_normal_vectors(const unsigned int face_no) const
{
  AssertDimension(dim, this->get_dimension());

  if (is_hyper_cube())
    {
      AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);
      return GeometryInfo<dim>::unit_normal_vector[face_no];
    }
  else if (dim == 2)
    {
      Assert(*this == ReferenceCells::Triangle, ExcInternalError());

      // Return the rotated vector
      return cross_product_2d(unit_tangential_vectors<dim>(face_no, 0));
    }
  else if (dim == 3)
    {
      return cross_product_3d(unit_tangential_vectors<dim>(face_no, 0),
                              unit_tangential_vectors<dim>(face_no, 1));
    }

  Assert(false, ExcNotImplemented());

  return {};
}



inline bool
ReferenceCell::standard_vs_true_line_orientation(
  const unsigned int  line,
  const unsigned char face_orientation_raw,
  const unsigned char line_orientation) const
{
  if (*this == ReferenceCells::Hexahedron)
    {
      static const bool bool_table[2][2][2][2] = {
        {{{true, false},    // lines 0/1, face_orientation=false,
                            // face_flip=false, face_rotation=false and true
          {false, true}},   // lines 0/1, face_orientation=false,
                            // face_flip=true, face_rotation=false and true
         {{true, true},     // lines 0/1, face_orientation=true,
                            // face_flip=false, face_rotation=false and true
          {false, false}}}, // lines 0/1, face_orientation=true,
                            // face_flip=true, face_rotation=false and true

        {{{true, true}, // lines 2/3 ...
          {false, false}},
         {{true, false}, {false, true}}}};

      const bool face_orientation = Utilities::get_bit(face_orientation_raw, 0);
      const bool face_flip        = Utilities::get_bit(face_orientation_raw, 2);
      const bool face_rotation    = Utilities::get_bit(face_orientation_raw, 1);

      return (static_cast<bool>(line_orientation) ==
              bool_table[line / 2][face_orientation][face_flip][face_rotation]);
    }
  else
    // TODO: This might actually be wrong for some of the other
    // kinds of objects. We should check this
    return true;
}



namespace internal
{
  template <typename T, std::size_t N>
  class NoPermutation : public dealii::ExceptionBase
  {
  public:
    /**
     * 构造函数。
     *
     */
    NoPermutation(const dealii::ReferenceCell &entity_type,
                  const std::array<T, N> &     vertices_0,
                  const std::array<T, N> &     vertices_1)
      : entity_type(entity_type)
      , vertices_0(vertices_0)
      , vertices_1(vertices_1)
    {}

    /**
     * 解除函数。
     *
     */
    virtual ~NoPermutation() noexcept override = default;

    /**
     * 打印错误信息到 @p out. 。
     *
     */
    virtual void
    print_info(std::ostream &out) const override
    {
      out << "[";

      const unsigned int n_vertices = entity_type.n_vertices();

      for (unsigned int i = 0; i < n_vertices; ++i)
        {
          out << vertices_0[i];
          if (i + 1 != n_vertices)
            out << ",";
        }

      out << "] is not a permutation of [";

      for (unsigned int i = 0; i < n_vertices; ++i)
        {
          out << vertices_1[i];
          if (i + 1 != n_vertices)
            out << ",";
        }

      out << "]." << std::endl;
    }

    /**
     * 实体类型。
     *
     */
    const dealii::ReferenceCell entity_type;

    /**
     * 第一组数值。
     *
     */
    const std::array<T, N> vertices_0;

    /**
     * 第二组数值。
     *
     */
    const std::array<T, N> vertices_1;
  };
} // namespace internal



template <typename T, std::size_t N>
inline unsigned char
ReferenceCell::compute_orientation(const std::array<T, N> &vertices_0,
                                   const std::array<T, N> &vertices_1) const
{
  AssertIndexRange(n_vertices(), N + 1);
  if (*this == ReferenceCells::Line)
    {
      const std::array<T, 2> i{{vertices_0[0], vertices_0[1]}};
      const std::array<T, 2> j{{vertices_1[0], vertices_1[1]}};

      // line_orientation=true
      if (i == std::array<T, 2>{{j[0], j[1]}})
        return 1;

      // line_orientation=false
      if (i == std::array<T, 2>{{j[1], j[0]}})
        return 0;
    }
  else if (*this == ReferenceCells::Triangle)
    {
      const std::array<T, 3> i{{vertices_0[0], vertices_0[1], vertices_0[2]}};
      const std::array<T, 3> j{{vertices_1[0], vertices_1[1], vertices_1[2]}};

      // face_orientation=true, face_rotation=false, face_flip=false
      if (i == std::array<T, 3>{{j[0], j[1], j[2]}})
        return 1;

      // face_orientation=true, face_rotation=true, face_flip=false
      if (i == std::array<T, 3>{{j[1], j[2], j[0]}})
        return 3;

      // face_orientation=true, face_rotation=false, face_flip=true
      if (i == std::array<T, 3>{{j[2], j[0], j[1]}})
        return 5;

      // face_orientation=false, face_rotation=false, face_flip=false
      if (i == std::array<T, 3>{{j[0], j[2], j[1]}})
        return 0;

      // face_orientation=false, face_rotation=true, face_flip=false
      if (i == std::array<T, 3>{{j[2], j[1], j[0]}})
        return 2;

      // face_orientation=false, face_rotation=false, face_flip=true
      if (i == std::array<T, 3>{{j[1], j[0], j[2]}})
        return 4;
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      const std::array<T, 4> i{
        {vertices_0[0], vertices_0[1], vertices_0[2], vertices_0[3]}};
      const std::array<T, 4> j{
        {vertices_1[0], vertices_1[1], vertices_1[2], vertices_1[3]}};

      // face_orientation=true, face_rotation=false, face_flip=false
      if (i == std::array<T, 4>{{j[0], j[1], j[2], j[3]}})
        return 1;

      // face_orientation=true, face_rotation=true, face_flip=false
      if (i == std::array<T, 4>{{j[2], j[0], j[3], j[1]}})
        return 3;

      // face_orientation=true, face_rotation=false, face_flip=true
      if (i == std::array<T, 4>{{j[3], j[2], j[1], j[0]}})
        return 5;

      // face_orientation=true, face_rotation=true, face_flip=true
      if (i == std::array<T, 4>{{j[1], j[3], j[0], j[2]}})
        return 7;

      // face_orientation=false, face_rotation=false, face_flip=false
      if (i == std::array<T, 4>{{j[0], j[2], j[1], j[3]}})
        return 0;

      // face_orientation=false, face_rotation=true, face_flip=false
      if (i == std::array<T, 4>{{j[2], j[3], j[0], j[1]}})
        return 2;

      // face_orientation=false, face_rotation=false, face_flip=true
      if (i == std::array<T, 4>{{j[3], j[1], j[2], j[0]}})
        return 4;

      // face_orientation=false, face_rotation=true, face_flip=true
      if (i == std::array<T, 4>{{j[1], j[0], j[3], j[2]}})
        return 6;
    }

  Assert(false, (internal::NoPermutation<T, N>(*this, vertices_0, vertices_1)));

  return -1;
}



template <typename T, std::size_t N>
inline std::array<T, N>
ReferenceCell::permute_according_orientation(
  const std::array<T, N> &vertices,
  const unsigned int      orientation) const
{
  std::array<T, 4> temp;

  if (*this == ReferenceCells::Line)
    {
      switch (orientation)
        {
          case 1:
            temp = {{vertices[0], vertices[1]}};
            break;
          case 0:
            temp = {{vertices[1], vertices[0]}};
            break;
          default:
            Assert(false, ExcNotImplemented());
        }
    }
  else if (*this == ReferenceCells::Triangle)
    {
      switch (orientation)
        {
          case 1:
            temp = {{vertices[0], vertices[1], vertices[2]}};
            break;
          case 3:
            temp = {{vertices[1], vertices[2], vertices[0]}};
            break;
          case 5:
            temp = {{vertices[2], vertices[0], vertices[1]}};
            break;
          case 0:
            temp = {{vertices[0], vertices[2], vertices[1]}};
            break;
          case 2:
            temp = {{vertices[2], vertices[1], vertices[0]}};
            break;
          case 4:
            temp = {{vertices[1], vertices[0], vertices[2]}};
            break;
          default:
            Assert(false, ExcNotImplemented());
        }
    }
  else if (*this == ReferenceCells::Quadrilateral)
    {
      switch (orientation)
        {
          case 1:
            temp = {{vertices[0], vertices[1], vertices[2], vertices[3]}};
            break;
          case 3:
            temp = {{vertices[2], vertices[0], vertices[3], vertices[1]}};
            break;
          case 5:
            temp = {{vertices[3], vertices[2], vertices[1], vertices[0]}};
            break;
          case 7:
            temp = {{vertices[1], vertices[3], vertices[0], vertices[2]}};
            break;
          case 0:
            temp = {{vertices[0], vertices[2], vertices[1], vertices[3]}};
            break;
          case 2:
            temp = {{vertices[2], vertices[3], vertices[0], vertices[1]}};
            break;
          case 4:
            temp = {{vertices[3], vertices[1], vertices[2], vertices[0]}};
            break;
          case 6:
            temp = {{vertices[1], vertices[0], vertices[3], vertices[2]}};
            break;
          default:
            Assert(false, ExcNotImplemented());
        }
    }
  else
    {
      AssertThrow(false, ExcNotImplemented());
    }

  std::array<T, N> temp_;
  std::copy_n(temp.begin(), N, temp_.begin());

  return temp_;
}


DEAL_II_NAMESPACE_CLOSE

#endif


