//include/deal.II-translator/fe/fe_q_bubbles_0.txt
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


#ifndef dealii_fe_q_bubbles_h
#define dealii_fe_q_bubbles_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor_product_polynomials_bubbles.h>

#include <deal.II/fe/fe_q_base.h>

DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * 实现标量拉格朗日有限元 $Q_p^+$
 * ，得到每个坐标方向上连续的、程度为 @p p
 * 的分片多项式的有限元空间，加上一些由附加形状函数
 * $\varphi_j(\mathbf x) = 2^{p-1}\left(x_j-frac 12\right)^{p-1}
 * \left[\prod_{i=0}^{dim-1}(x_i(1-x_i))\right]$
 * 跨越的（非归一化）气泡富集空间。 对于 $j=0,\ldots,dim-1$
 * 。 如果 $p$
 * 为1，那么第一个因素就会消失，人们会得到以单元格中点为中心的通常的气泡函数。因为这些最后的形状函数的多项式程度是
 * $p+1$
 * ，这个类所描述的空间中形状函数的总体多项式程度是
 * $p+1$ 。
 * 该类使用基于等距或给定支持点的张量积多项式来实现，与人们向FE_Q类的构造函数提供支持点的方式相同。
 * 关于<tt>spacedim</tt>模板参数的更多信息，请查阅FiniteElement类的文档，或者Triangulation的文档。
 * 由于大的 $p$
 * 的富集度几乎都很小，质量和刚度矩阵的条件数随着 $p$
 * 的增加而迅速增加。下面你可以看到与FE_Q(QGaussLobatto(p+1))在dim=1时的比较。
 * <p ALIGN="center">
 @image html fe_q_bubbles_conditioning.png
 * </p> 因此，对于  $p>3$  应谨慎使用该元素。
 *
 *  <h3>Implementation</h3>
 * 构造函数创建一个TensorProductPolynomials对象，其中包括 @p
 * LagrangeEquidistant 度数为 @p p
 * 的多项式的张量乘积加上气泡富集。这个 @p
 * TensorProductPolynomialsBubbles
 * 对象提供形状函数的所有值和导数。在给出正交规则的情况下，构造函数创建一个TensorProductPolynomialsBubbles对象，其中包括
 * @p Lagrange 多项式与 @p points
 * 的支持点的张量乘积以及上面定义的气泡丰富度。
 * 此外，构造函数还填充了 @p interface_constrains, 、 @p
 * 的延长（嵌入）和 @p restriction 矩阵。
 *
 *  <h3>Numbering of the degrees of freedom (DoFs)</h3>
 * TensorProductPolynomialsBubbles所代表的形状函数的原始排序是张量积的编号。然而，单元格上的形状函数被重新编号，从支持点在顶点的形状函数开始，然后是在线上，在四边形上，最后（对于三维）在六边形上。最后，在单元格的中间有支持点的气泡富集。
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_Q_Bubbles : public FE_Q_Base<dim, spacedim>
{
public:
  /**
   * 度数为 @p p
   * 的张量乘积多项式的构造器，加上气泡富集点
   *
   */
  FE_Q_Bubbles(const unsigned int p);

  /**
   * 用于支持点为 @p points
   * 的张量积多项式的构造器，加上基于一维正交公式的气泡富集。
   * 然后有限元的度数是<tt>points.size()</tt>，与FE_Q类的相应情况相比，其加数来自于额外的气泡函数。更多信息请参见FE_Q构造函数的文档。
   * 注意，第一个点必须是0，最后一个是1。
   *
   */
  FE_Q_Bubbles(const Quadrature<1> &points);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_Q_Bubbles<dim>(degree)</tt>，其中
   * @p dim 和 @p degree 由适当的值替换。
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
   * FE_Q_Bubbles元素时才能使用。否则，会抛出一个
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented
   * 类型的异常。
   *
   */
  virtual void
  get_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                           FullMatrix<double> &matrix) const override;

  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case) const override;

  virtual const FullMatrix<double> &
  get_restriction_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case) const override;

  /**
   * 检查一个面的非零值。    如果形状函数 @p shape_index
   * 在面的非零值，该函数返回 @p true,   @p face_index.
   * 在FiniteElement中实现接口。
   *
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

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
   * 返回restriction_is_additive标志。只有气泡富集的最后一个成分是真的。
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

  /**
   * 额外气泡函数的数量
   *
   */
  const unsigned int n_bubbles;
};



 /*@}*/ 


DEAL_II_NAMESPACE_CLOSE

#endif


