//include/deal.II-translator/fe/fe_trace_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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

#ifndef dealii_fe_trace_h
#define dealii_fe_trace_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_poly_face.h>
#include <deal.II/fe/fe_q.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个有限元，它是FE_Q元素的轨迹，它是面的多项式的张量乘积，在单元的内部没有定义，而且是连续的。面上的基函数是由一维拉格朗日多项式的张量乘积形成的，张量乘积上有2度以下的等距点和3度开始的高斯-洛巴托点。
 * 这个有限元是FE_Q在面的跟踪空间。
 *
 *
 * @note
 * 由于这些只是面的有限元，只有FEFaceValues和FESubfaceValues能够从任何面的多项式中提取合理的数值。为了使FESystem的使用更加简单，FEValues对象使用这个有限元空间不会失败，但所有提取的形状函数值将等于零。
 *
 *
 */

template <int dim, int spacedim = dim>
class FE_TraceQ
  : public FE_PolyFace<TensorProductPolynomials<dim - 1>, dim, spacedim>
{
public:
  /**
   * 度数为<tt>p</tt>的张量乘积多项式的构造函数。使用此构造函数创建的形状函数对应于每个坐标方向的Legendre多项式。
   *
   */
  FE_TraceQ(unsigned int p);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_DGQ<dim>(degree)</tt>，其中<tt>dim</tt>和<tt>degree</tt>由适当的值替换。
   *
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * 在FiniteElement类中实现相应的函数。
   * 由于当前元素是插值的，所以节点值正好是支持点的值。此外，由于当前元素是标量的，支持点的值需要是长度为1的向量。
   *
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

  /**
   * 如果形状函数 @p shape_index 在面 @p face_index.
   * 的某处有非零函数值，该函数返回 @p true, 。
   *
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  /**
   * 返回该元素的常数模式列表。对于这个元素，它只是简单地返回一行，所有条目都设置为真。
   *
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  /**
   * 返回这个元素是否以新的方式实现了它的悬挂节点约束，这必须用于使元素
   * "hp-compatible"。
   *
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * 返回从一个元素的一个面插值到相邻元素的面的矩阵。
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。这个元素只为同一类型和FE_Nothing的元素提供插值矩阵。对于所有其他元素，会抛出一个
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
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。这个元素只为同一类型和FE_Nothing的元素提供插值矩阵。对于所有其他元素，会抛出一个
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
   * @copydoc   FiniteElement::compare_for_domination() 。
   *
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

private:
  /**
   * 存储FE_Q的副本，用于委托hp-constraints的功能。
   *
   */
  FE_Q<dim, spacedim> fe_q;

  /**
   * 返回每个顶点、线条、四边形、六边形的道夫向量。
   *
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int deg);
};



/**
 * 一维的FE_TraceQ，即元素顶点上的自由度。
 *
 *
 */
template <int spacedim>
class FE_TraceQ<1, spacedim> : public FE_FaceQ<1, spacedim>
{
public:
  /**
   * 构造函数。
   *
   */
  FE_TraceQ(const unsigned int p);

  /**
   * 返回元素的名称
   *
   */
  std::string
  get_name() const override;
};


DEAL_II_NAMESPACE_CLOSE

#endif


