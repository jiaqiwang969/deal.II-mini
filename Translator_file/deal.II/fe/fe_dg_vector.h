//include/deal.II-translator/fe/fe_dg_vector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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

#ifndef dealii_fe_dg_vector_h
#define dealii_fe_dg_vector_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_bdm.h>
#include <deal.II/base/polynomials_nedelec.h>
#include <deal.II/base/polynomials_raviart_thomas.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 基于矢量值多项式的DG元素。
 * 这些元素使用矢量值多项式空间，因为它们已经被引入H<sup>div</sup>和H<sup>curl</sup>符合的有限元，但不使用这些元素的通常连续性。因此，它们适用于涉及这些函数空间的DG和混合公式。
 * 模板参数<tt>PolynomialType</tt>指的是一个矢量值的多项式空间，比如PolynomialsRaviartThomas或者PolynomialsNedelec。注意，多项式空间的维度和参数<tt>dim</tt>必须重合。
 *
 *
 * @ingroup febase
 *
 *
 */
template <class PolynomialType, int dim, int spacedim = dim>
class FE_DGVector : public FE_PolyTensor<dim, spacedim>
{
public:
  /**
   * 程度为 @p p. 的向量元素的构造函数
   *
   */
  FE_DGVector(const unsigned int p, MappingKind m);

  /**
   * 返回一个唯一标识有限元的字符串。这个类返回`FE_DGVector_`加上一块取自多项式对象返回的名称，再加上`<dim>(degree)`，其中
   * @p dim 和 @p degree 用适当的值代替。
   *
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * 如果形状函数 @p shape_index
   * 在面的某处有非零函数值，这个函数就会返回 @p true,
   * 对于这个元素，我们总是返回 @p true.  。
   *
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  virtual std::size_t
  memory_consumption() const override;

private:
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
   * 与单元无关的数据字段。
   * 关于这个类的一般用途的信息，请看基类的文档。
   *
   */
  class InternalData : public FiniteElement<dim>::InternalDataBase
  {
  public:
    /**
     * 具有正交点的形状函数值的数组。每个形状函数都有一行，包含每个正交点的值。
     * 由于形状函数是矢量值（有多少分量就有多少空间维度），所以值是一个张量。
     * 在这个数组中，我们将形状函数的值存储在单元格上的正交点。然后，向实空间单元的转换只需与映射的雅各布系数相乘即可完成。
     *
     */
    std::vector<std::vector<Tensor<1, dim>>> shape_values;

    /**
     * 包含正交点的形状函数梯度的数组。每个形状函数都有一行，包含每个正交点的值。
     * 我们将梯度存储在单元格的正交点上。然后我们只需要在访问实际单元格时应用转换（这是一个矩阵-向量乘法）。
     *
     */
    std::vector<std::vector<Tensor<2, dim>>> shape_gradients;
  };
  Table<3, double> interior_weights;
};



/**
 * 一个基于FE_Nedelec的多项式空间的矢量值DG元素。这个类实现了一个
 * "破碎的
 * "有限元空间，它在单元之间是不连续的，在每个单元上的形状函数等于Nedelec元素的形状函数。
 * 相关的类FE_DGRT用于  step-61  。
 *
 * @ingroup fe
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_DGNedelec : public FE_DGVector<PolynomialsNedelec<dim>, dim, spacedim>
{
public:
  /**
   * 度的不连续N&eacute;d&eacute;lec元素的构造函数  @p p.  。
   *
   */
  FE_DGNedelec(const unsigned int p);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_DGNedelec<dim>(degree)</tt>，其中
   * @p dim 和 @p degree 被适当的值替换。
   *
   */
  virtual std::string
  get_name() const override;
};



/**
 * 一个基于FE_RaviartThomas的多项式空间的矢量值DG元素。该类实现了一个
 * "破碎的
 * "有限元空间，在单元之间是不连续的，在每个单元上的形状函数与Raviart-Thomas元的形状函数相同。
 * 该类在  step-61  中使用。
 *
 *
 * @ingroup fe
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_DGRaviartThomas
  : public FE_DGVector<PolynomialsRaviartThomas<dim>, dim, spacedim>
{
public:
  /**
   * 度的拉维奥特-托马斯元素的构造函数  @p p.  。
   *
   */
  FE_DGRaviartThomas(const unsigned int p);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_DGRaviartThomas<dim>(degree)</tt>，其中
   * @p dim 和 @p 度被适当的值取代。
   *
   */
  virtual std::string
  get_name() const override;
};



/**
 * 一个基于FE_BDM的多项式空间的矢量值DG元素。该类实现了一个
 * "破碎的
 * "有限元空间，在单元之间是不连续的，在每个单元上的形状函数与BDM元素相同。
 * 相关的类FE_DGRT用于  step-61  。
 *
 *
 * @ingroup fe
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_DGBDM : public FE_DGVector<PolynomialsBDM<dim>, dim, spacedim>
{
public:
  /**
   * 度的不连续BDM元素的构造函数  @p p.  。
   *
   */
  FE_DGBDM(const unsigned int p);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_DGBDM<dim>(degree)</tt>，其中
   * @p dim 和 @p degree 由适当的值代替。
   *
   */
  virtual std::string
  get_name() const override;
};


DEAL_II_NAMESPACE_CLOSE

#endif


