//include/deal.II-translator/fe/fe_bernstein_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2021 by the deal.II authors
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

#ifndef dealii_fe_bernstein_h
#define dealii_fe_bernstein_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_q_base.h>

DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * 实现标量伯恩斯坦有限元 @p that
 * ，我们称之为FE_Bernstein，与FE_Q相类似，得到每个坐标方向上的连续、分片伯恩斯坦多项式
 * @p p
 * 度的有限元空间。这类空间是通过伯恩斯坦基础多项式的张量积多项式实现的。
 *
 *  该类的标准构造函数取该有限元的度数 @p p 。
 * 关于<tt>spacedim</tt>模板参数的更多信息，请查阅FiniteElement或Triangulation的文档。
 * <h3>Implementation</h3>
 * 构造函数创建一个TensorProductPolynomials对象，其中包括度数为
 * @p Bernstein 的多项式的张量积 @p p.  这个 @p
 * TensorProductPolynomials对象提供形状函数的所有值和导数。
 * <h3>Numbering of the degrees of freedom (DoFs)</h3>
 * TensorProductPolynomials所代表的形状函数的原始排序是张量乘法的编号。然而，单元格上的形状函数被重新编号，从支持点在顶点的形状函数开始，然后是在直线上，在四边形上，最后（对于三维）在六边形上。更多细节请参见FE_Q的文档。
 *
 *
 */

template <int dim, int spacedim = dim>
class FE_Bernstein : public FE_Q_Base<dim, spacedim>
{
public:
  /**
   * 度数为 @p p. 的张量乘积多项式的构造函数。
   *
   */
  FE_Bernstein(const unsigned int p);

  /**
   * FE_Bernstein在元素内部不是插值的，这使得这个元素不能定义插值矩阵。将会抛出一个异常。
   * 这个函数覆盖了来自FE_Q_Base的实现。
   *
   */
  virtual void
  get_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                           FullMatrix<double> &matrix) const override;

  /**
   * FE_Bernstein在元素内部没有插值，这使得这个元素无法定义限制矩阵。将会抛出一个异常。
   * 这个函数重写了来自FE_Q_Base的实现。
   *
   */
  virtual const FullMatrix<double> &
  get_restriction_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * FE_Bernstein在元素内部没有插值，这使得这个元素不能定义一个延长矩阵。将会抛出一个异常。
   * 这个函数重写了来自FE_Q_Base的实现。
   *
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * 返回从一个元素的一个面插值到相邻元素的面的矩阵。
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。FE_Bernstein元素家族只为相同类型的元素、有支持点的元素和FE_Nothing提供插值矩阵。对于所有其他元素，会抛出一个
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented
   * 类型的异常。
   *
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double> &                matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * 返回从一个元素的一个面插值到相邻元素的面的矩阵。
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。FE_Bernstein元素家族只为相同类型的元素、有支持点的元素和FE_Nothing提供插值矩阵。对于所有其他元素，会抛出一个
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented
   * 类型的异常。
   *
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source,
    const unsigned int                  subface,
    FullMatrix<double> &                matrix,
    const unsigned int                  face_no = 0) const override;

  /**
   * 返回这个元素是否以新的方式实现了它的悬挂节点约束，这必须被用来使元素
   * "hp-compatible"。
   *
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * 如果在一个顶点上，有几个有限元处于活动状态，hp代码首先为这些FEs中的每个自由度分配不同的全局索引。然后调用这个函数来找出其中哪些应该得到相同的值，从而可以得到相同的全局自由度指数。
   * 因此，该函数返回当前有限元对象的自由度与 @p fe_other,
   * 的自由度之间的相同性列表，该列表是对代表该特定顶点上活动的其他有限元之一的有限元对象的引用。该函数计算两个有限元对象的哪些自由度是相等的，这两个自由度的编号都在零和两个有限元的n_dofs_per_vertex()的相应值之间。每一对的第一个索引表示本元素的一个顶点自由度，而第二个是另一个有限元素的相应索引。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理线上自由度。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理四边形上的自由度。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  /**
   * @copydoc   FiniteElement::compare_for_domination() .
   *
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_Bernstein<dim>(degree)</tt>，其中
   * @p dim 和 @p degree 由适当的值代替。
   *
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

protected:
  /**
   * 仅供内部使用。它的全称是 @p get_dofs_per_object_vector
   * 函数，它创建了 @p dofs_per_object
   * 向量，在构造函数中需要传递给 @p
   * FiniteElementData的构造函数。
   *
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree);

  /**
   * 该函数将伯恩斯坦基函数的编号从分层编号改为列举式编号。
   *
   */
  TensorProductPolynomials<dim>
  renumber_bases(const unsigned int degree);
};



 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


