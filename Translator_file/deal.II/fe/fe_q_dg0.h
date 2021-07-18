//include/deal.II-translator/fe/fe_q_dg0_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2021 by the deal.II authors
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


#ifndef dealii_fe_q_dg0_h
#define dealii_fe_q_dg0_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor_product_polynomials_const.h>

#include <deal.II/fe/fe_q_base.h>

DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * 实现标量拉格朗日有限元 @p Qp+DG0
 * ，得到每个坐标方向上连续的、度数为 @p
 * p的分片多项式的有限元空间加上局部常数函数的空间。这个类使用基于等距或给定支持点的张量积多项式来实现。
 * 该类的标准构造函数取该有限元的度数 @p p
 * 。或者，它可以取一个正交公式 @p points
 * ，定义一个坐标方向上拉格朗日插值的支持点。
 * 关于<tt>spacedim</tt>模板参数的更多信息，请查阅FiniteElement或Triangulation的文档。
 * 关于这个元素的更多信息，请看。Boffi, D., et al. "Stokes
 * Finite
 * Elements的局部质量保护"。科学计算杂志（2012）：1-18。
 * <h3>Implementation</h3>
 * 构造函数创建一个TensorProductPolynomials对象，其中包括 @p
 * LagrangeEquidistant 度 @p p
 * 的多项式的张量乘积加上局部常数函数。这个 @p
 * TensorProductPolynomialsConst
 * 对象提供形状函数的所有值和导数。在给出正交规则的情况下，构造函数创建一个TensorProductPolynomialsConst对象，其中包括
 * @p 拉格朗日多项式与 @p points
 * 的支持点的张量乘积和一个局部常数函数。
 * 此外，构造函数还填充了 @p interface_constrains, 、 @p
 * 的延长（嵌入）和 @p restriction 矩阵。 <h3>Numbering of the
 * degrees of freedom (DoFs)</h3>
 * TensorProductPolynomialsConst所代表的形状函数的原始排序是张量乘法的编号。然而，单元格上的形状函数被重新编号，从支持点在顶点的形状函数开始，然后是在线上，在四边形上，最后（对于三维）在六边形上。最后，在单元格的中间有一个不连续形状函数的支持点。为了明确起见，这些编号列在下面。
 * <h4>Q1 elements</h4>  <ul>   <li>  1D情况。
 * @verbatim
 *    0---2---1
 * @endverbatim
 * <li>  二维情况。
 * @verbatim
 *    2-------3
 *    |       |
 *    |   5   |
 *    |       |
 *    0-------1
 * @endverbatim
 * <li>  3D情况。
 * @verbatim
 *       6-------7        6-------7
 *      /|       |       /       /|
 *     / |       |      /       / |
 *    /  |       |     /       /  |
 *   4   |  8    |    4-------5   |
 *   |   2-------3    |       |   3
 *   |  /       /     |       |  /
 *   | /       /      |       | /
 *   |/       /       |       |/
 *   0-------1        0-------1
 * @endverbatim
 *
 * 自由度的支持点的各自坐标值如下。  <ul>   <li>  索引0：<tt>[ 0, 0, 0]/tt>；  <li>  索引1：<tt>[ 1, 0, 0]/tt>；  <li>  索引2：<tt>[ 0, 1, 0]/tt>；  <li>  索引3：<tt>[ 1, 1, 0]/tt>；  <li>  索引4。<tt>[ 0, 0, 1]</tt>;  <li>  索引5：<tt>[ 1, 0, 1]</tt>;  <li>  索引6。<tt>[ 0, 1, 1]</tt>;  <li>  索引7：<tt>[ 1, 1, 1]</tt>;  <li>  索引8：<tt>[ 1/2, 1/2, 1/2]</tt>;  </ul>   </ul>  <h4>Q2 elements</h4>  <ul>   <li>  1D情况。
 * @verbatim
 *    0---2---1
 * @endverbatim
 * 索引3与索引2具有相同的坐标
 * <li>  二维情况。
 * @verbatim
 *    2---7---3
 *    |       |
 *    4   8   5
 *    |       |
 *    0---6---1
 * @endverbatim
 * 索引9的坐标与索引2的坐标相同
 * <li>  3D情况。
 * @verbatim
 *       6--15---7        6--15---7
 *      /|       |       /       /|
 *    12 |       19     12      1319
 *    /  18      |     /       /  |
 *   4   |       |    4---14--5   |
 *   |   2---11--3    |       |   3
 *   |  /       /     |      17  /
 *  16 8       9     16       | 9
 *   |/       /       |       |/
 *   0---10--1        0---8---1
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
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
 *      /|       |       /       /|
 *     / |  23   |      /  25   / |
 *    /  |       |     /       /  |
 *     |       |
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
 *   |20-------*    |       |21
 *   |  /       /     |   22  |  /
 *   | /  24   /      |       | /
 *   |/       /       |       |/
 *
 *
 *
 *
 *
 *
 *
 *
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
 * 中心顶点的编号为26和27。
 * 自由度的支持点的各自坐标值如下。  <ul>   <li>  索引0：<tt>[0, 0, 0]</tt>;  <li>  索引1：<tt>[1, 0, 0]</tt>;  <li>  索引2：<tt>[0, 1, 0]</tt>;  <li>  索引3：<tt>[1, 1, 0]/tt>;  <li>  索引4: <tt>[0, 0, 1]</tt>;  <li>  索引5：<tt>[1, 0, 1]</tt>;  <li>  索引6。<tt>[0, 1, 1]</tt>;  <li>  索引7：<tt>[1, 1, 1]</tt>;  <li>  索引8：<tt>[0, 1/2, 0]/tt>;  <li>  ] 索引9：<tt>[1，1/2，0]/tt>； <li>  索引10：<tt>[1/2，0，0]/tt>； <li>  索引11：<tt>[1/2，1，0]/tt>； <li>  ] 索引12：<tt>[0，1/2，1]</tt>； <li>  索引13：<tt>[1，1/2，1]</tt>； <li>  索引14：<tt>[1/2，0，1]</tt>； <li>  ] 索引15：<tt>[1/2, 1, 1]/tt>;  <li>  索引16：<tt>[0, 0, 1/2]/tt>;  <li>  索引17：<tt>[1, 0, 1/2]/tt>;  <li>  ] 索引18：<tt>[0, 1, 1/2]</tt>;  <li>  索引19：<tt>[1, 1, 1/2]</tt>;  <li>  索引20：<tt>[0, 1/2, 1/2]/tt>;  <li>  ] 索引21：<tt>[1，1/2，1/2]</tt>;  <li>  索引22：<tt>[1/2，0，1/2]</tt>;  <li>  索引23：<tt>[1/2，1/2]</tt>;  <li>  ] 索引24：<tt>[1/2, 1/2, 0]</tt>;  <li>  索引25：<tt>[1/2, 1/2, 1]</tt>;  <li>  索引26：<tt>[1/2, 1/2]</tt>;  <li>  ] 索引27：<tt>[1/2, 1/2, 1/2]</tt>;  </ul>   </ul>  <h4>Q3 elements</h4>  <ul>   <li>  1D情况。
 * @verbatim
 *    0--2-4-3--1
 * @endverbatim
 * <li>  二维情况。
 * @verbatim
 *    2--10-11-3
 *    |        |
 *    5  14 15 7
 *    |    16  |
 *    4  12 13 6
 *    |        |
 *    0--8--9--1
 * @endverbatim
 * </ul>  <h4>Q4 elements</h4>  <ul>   <li>  一维情况。
 * @verbatim
 *    0--2--3--4--1
 * @endverbatim
 * 索引5与索引3具有相同的坐标
 * <li>  二维情况。
 * @verbatim
 *    2--13-14-15-3
 *    |           |
 *    6  22 23 24 9
 *    |           |
 *    5  19 20 21 8
 *    |           |
 *    4  16 17 18 7
 *    |           |
 *    0--10-11-12-1
 * @endverbatim
 * 索引21与索引20具有相同的坐标 *  </ul>
 *
 */
template <int dim, int spacedim = dim>
class FE_Q_DG0 : public FE_Q_Base<dim, spacedim>
{
public:
  /**
   * 度数为 @p p
   * 的张量积多项式加上局部常数函数的构造器。
   *
   */
  FE_Q_DG0(const unsigned int p);

  /**
   * 用于支持点为 @p points
   * 的张量积多项式的构造器，加上基于一维正交公式的局部常数函数。有限元的程度是<tt>points.size()-1</tt>。
   * 注意，第一个点必须是0，最后一个是1。
   *
   */
  FE_Q_DG0(const Quadrature<1> &points);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_Q_DG0<dim>(degree)</tt>，
   * @p dim 和 @p degree 用适当的值替换。
   *
   */
  virtual std::string
  get_name() const override;

  // documentation inherited from the base class
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

  /**
   * 返回从给定的有限元插值到现在的矩阵。
   * 然后矩阵的大小是 @p dofs_per_cell
   * 乘以<tt>source.n_dofs_per_cell()</tt>。
   * 这些矩阵只有在源元素也是 @p
   * FE_Q_DG0元素时才能使用。否则，会抛出一个
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented
   * 类型的异常。
   *
   */
  virtual void
  get_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                           FullMatrix<double> &matrix) const override;


  /**
   * 如果形状函数 @p shape_index 在面 @p face_index.
   * 的某处有非零函数值，该函数返回 @p true, 。
   *
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  /**
   * 返回该元素的常数模式列表。对于这个元素，尽管该元素是标量的，但有两个常量模式。第一个常数模式是通常FE_Q基础的所有1，第二个常数模式只使用不连续部分。
   *
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * @copydoc   FiniteElement::compare_for_domination() .
   *
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

private:
  /**
   * 返回 restriction_is_additive
   * 标志。只有最后一个成分是真的。
   *
   */
  static std::vector<bool>
  get_riaf_vector(const unsigned int degree);

  /**
   * 仅供内部使用。它的全称是 @p get_dofs_per_object_vector
   * 函数，它创建了 @p dofs_per_object
   * 向量，在构造函数内需要传递给 @p
   * FiniteElementData的构造函数。
   *
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree);
};



 /*@}*/ 


DEAL_II_NAMESPACE_CLOSE

#endif


