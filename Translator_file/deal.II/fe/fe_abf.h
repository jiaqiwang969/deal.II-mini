//include/deal.II-translator/fe/fe_abf_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2021 by the deal.II authors
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

#ifndef dealii_fe_abf_h
#define dealii_fe_abf_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_abf.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * Arnold-Boffi-Falk（ABF）元素的实现，符合空间H<sup>div</sup>。这些元素在网格单元之间产生连续的法向分量的矢量场。
 * 这些元素是基于Arnold, Boffi and
 * Falk的文章：四边形H(div)有限元素，SIAM J. Numer.
 * 分析。第42卷，第6期，第2429-2451页
 * 在这篇文章中，作者证明了通常的RT元素以及BDM和其他提议的H(div)的有限维度子空间在任意的FE网格上不能正常工作。也就是说，在这些网格上收敛率会恶化。作为一个解决方案，作者提出了ABF元素，在本模块中实现。
 * 这类元素没有实现在一维的情况下（<tt>spacedim !=
 * dim</tt>）。
 * @todo
 * 即使这个元素是为二维和三维空间实现的，节点值的定义也依赖于三维中一致方向的面。因此，在复杂的网格上应该注意。
 * <h3>Interpolation</h3>
 * @ref GlossInterpolation  与RT元素相关的 "插值 "
 * 算子的构造是这样的：插值和计算发散是互换的操作。我们从插值任意函数以及#限制矩阵中要求这一点。
 * 这可以通过两种插值方案实现，FE_RaviartThomasNodal中的简化方案和这里的原始方案。
 * <h4>Node values on edges/faces</h4>
 * 在边缘或面上， @ref GlossNodes "节点值 "
 * 是内插函数的法线分量相对于RT多项式轨迹的矩。由于在一个边缘/面的度数为<i>k</i>的RT空间的法线轨迹是<i>Q<sub>k</sub></i>的空间，因此相对于这个空间的矩被取走。
 * <h4>Interior node values</h4>
 * 高阶RT空间有内部节点。这些是相对于<i>Q<sub>k</sub></i>中的函数在单元上的梯度所取的矩（这个空间是混合表述中RT<sub>k</sub>的匹配空间）。
 * <h4>Generalized support points</h4>
 * 上面的节点值依赖于积分，这些积分将由正交规则本身计算出来。广义支持点是一组点，使这种正交能够以足够的精度进行。需要的点是每个面上的QGauss<sub>k+1</sub>以及单元内部的QGauss<sub>k</sub>（或者对于RT<sub>0</sub>来说没有）。更多信息请参见 @ref GlossGeneralizedSupport "广义支持点的术语条目"
 * 。
 *
 *
 */
template <int dim>
class FE_ABF : public FE_PolyTensor<dim>
{
public:
  /**
   * 程度为 @p p. 的ABF元素的构造函数。
   *
   */
  FE_ABF(const unsigned int p);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_ABF<dim>(degree)</tt>，其中
   * @p dim 和 @p degree 由适当的值代替。
   *
   */
  virtual std::string
  get_name() const override;

  /**
   * 如果形状函数 @p shape_index
   * 在面的某处有非零函数值，该函数返回 @p true, 。
   * 现在，这只在一维的RT0中实现。否则，总是返回  @p true.
   * 。
   *
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  // documentation inherited from the base class
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

  virtual std::size_t
  memory_consumption() const override;

  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

private:
  /**
   * ABF元素的顺序。最低阶的元素通常被称为RT0，尽管它们的形状函数是片状二次方程。
   *
   */
  const unsigned int rt_order;

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
   * 初始化FiniteElement类的 @p generalized_support_points 字段，用插值权重（#boundary_weights和#interior_weights）填充表格。从构造函数中调用。    更多信息请参见 @ref GlossGeneralizedSupport "广义支持点的词汇表条目"
   * 。
   *
   */
  void
  initialize_support_points(const unsigned int rt_degree);

  /**
   * 初始化从细化网格单元上的函数到父单元的插值。根据Raviart-Thomas元素的理念，这个限制算子弱地保留了一个函数的发散性。
   *
   */
  void
  initialize_restriction();

  /**
   * 与单元格无关的数据字段。
   * 关于这个类的一般用途的信息，请参见基类的文档。
   *
   */
  class InternalData : public FiniteElement<dim>::InternalDataBase
  {
  public:
    /**
     * 带有正交点的形状函数值的数组。每个形状函数都有一行，包含每个正交点的值。
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

  /**
   * 这些是计算积分时乘以#generalized_face_support_points中的一个函数的系数。它们的组织方式是，每个广义面支持点有一行，面的每个自由度有一列。
   *
   */
  Table<2, double> boundary_weights;
  /**
   * 用于内部自由度内插的预计算系数。这个表的原理与#boundary_weights的原理相同。只是，这个表有第三个坐标，用于评估组件的空间方向。
   *
   */
  Table<3, double> interior_weights;



  /**
   * 这些是计算积分时乘以#generalized_face_support_points中的一个函数的系数。它们的组织方式是，每个广义面支持点有一行，面的每个自由度有一列。
   *
   */
  Table<2, double> boundary_weights_abf;
  /**
   * 用于内部自由度内插的预计算系数。这个表的原理与#boundary_weights的原理相同。只是，这个表有第三个坐标，用于评估组件的空间方向。
   *
   */
  Table<3, double> interior_weights_abf;

  /**
   * 初始化包络模式和符号变化模式。
   * @note
   * 这个函数还没有完全填充正确的实现。它需要在未来的版本中统一实现，以便在包含有翻转面的单元格的网格上工作。
   *
   */
  void
  initialize_quad_dof_index_permutation_and_sign_change();

  // Allow access from other dimensions.
  template <int dim1>
  friend class FE_ABF;
};



 /*@}*/ 


DEAL_II_NAMESPACE_CLOSE

#endif


