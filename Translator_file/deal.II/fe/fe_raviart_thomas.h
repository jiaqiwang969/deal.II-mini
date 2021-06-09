//include/deal.II-translator/fe/fe_raviart_thomas_0.txt
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

#ifndef dealii_fe_raviart_thomas_h
#define dealii_fe_raviart_thomas_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_raviart_thomas.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * Raviart-Thomas（RT）元素的实现。Raviart-Thomas空间的设计是为了解决解只存在于空间
 * $H^\text{div}=\{ {\mathbf u} \in L_2: \text{div}\, {\mathbf u} \in L_2\}$
 * 的问题，而不是更常用的空间 $H^1=\{ u \in L_2: \nabla u \in
 * L_2\}$
 * 。换句话说，解决方案必须是一个矢量场，其发散是可平方整除的，但其梯度可能不是可平方整除的。这个空间（和这些元素）的典型应用是拉普拉斯方程的混合表述和相关情况，例如，见
 * step-20  。 $H^\text{div}$
 * 中的函数的决定性特征是它们一般是不连续的
 *
 * - 但如果你在2D中画一条线（或在3D中画一个面），那么矢量场的<i>normal</i>分量必须在线（或面）上连续，尽管切向分量可能不是。因此，Raviart-Thomas元素的构造是这样的：（i）它是 @ref vector_valued "矢量值"，（ii）形状函数是不连续的，但（iii）每个形状函数所代表的矢量场的法向分量在单元面上是连续的。
 * Raviart-Thomas元素的其他特性是：(i)它是 @ref GlossPrimitive "非原始元素"
 * ；(ii)形状函数的定义使某些面的积分为零或一，而不是常见的某些点值为零或一的情况。然而，有一个FE_RaviartThomasNodal元素使用点值）。
 * 我们遵循通常使用的
 *
 * - 尽管很混乱
 *
 * - RT元素的 "度 "的定义。具体来说，元素的 "度 "表示有限元空间中包含的<i>largest complete polynomial subspace</i>的多项式度，即使该空间可能包含更高多项式度的形状函数。因此，最低阶元素是FE_RaviartThomas(0)，即 "零度 "的Raviart-Thomas元素，尽管这个空间的函数一般是每个变量的一度多项式。这种 "度 "的选择意味着函数本身的近似顺序是<i>degree+1</i>，就像通常的多项式空间一样。如此选择的编号意味着序列为@f[
 * Q_{k+1} \stackrel{\text{grad}}{\rightarrow} \text{Nedelec}_k
 * \stackrel{\text{curl}}{\rightarrow} \text{RaviartThomas}_k
 * \stackrel{\text{div}}{\rightarrow} DGQ_{k} @f]。
 * 该类没有在二维一的情况下实现（<tt>spacedim != dim</tt>）。
 *
 *  <h3>Interpolation</h3>
 * 与RT元素相关的 @ref GlossInterpolation "插值 "
 * 算子的构造是，插值和计算发散是互换的操作。我们从插值任意函数以及#限制矩阵中要求这一点。
 * 这可以通过两种插值方案实现，FE_RaviartThomasNodal中的简化方案和这里的原始方案。
 * <h4>Node values on edges/faces</h4>
 * 在边缘或面上， @ref GlossNodes "节点值 "
 * 是内插函数的法线分量相对于RT多项式轨迹的矩。由于在一个边缘/面的度数为<i>k</i>的RT空间的法线轨迹是空间<i>Q<sub>k</sub></i>，因此相对于这个空间的矩。
 * <h4>Interior node values</h4>
 * 高阶RT空间有内部节点。这些是相对于<i>Q<sub>k</sub></i>中的函数在单元上的梯度所取的矩（这个空间是混合表述中RT<sub>k</sub>的匹配空间）。
 * <h4>Generalized support points</h4>
 * 上面的节点值依赖于积分，这些积分将由正交规则本身计算出来。广义支持点是一组点，使这种正交能够以足够的精度进行。需要的点是每个面上的QGauss<sub>k+1</sub>以及单元内部的QGauss<sub>k+1</sub>（或者对于RT<sub>0</sub>来说没有）。
 *
 *
 */
template <int dim>
class FE_RaviartThomas : public FE_PolyTensor<dim>
{
public:
  /**
   * 度数为 @p p. 的Raviart-Thomas元素的构造函数。
   *
   */
  FE_RaviartThomas(const unsigned int p);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_RaviartThomas<dim>(degree)</tt>，其中
   * @p dim 和 @p degree 用适当的值替换。
   *
   */
  virtual std::string
  get_name() const override;

  // documentation inherited from the base class
  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

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

  /**
   * 返回元素的恒定模式列表。这个方法目前没有正确实现，因为它对所有元件都返回1。
   *
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  virtual std::size_t
  memory_consumption() const override;

private:
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
   * 初始化FiniteElement类的 @p generalized_support_points
   * 字段，用插值权重（#boundary_weights和#interior_weights）填充表格。从构造函数中调用。
   *
   */
  void
  initialize_support_points(const unsigned int rt_degree);

  /**
   * 初始化从细化网格单元上的函数到父单元的插值。根据Raviart-Thomas元素的理念，这个限制算子弱地保留了函数的发散性。
   *
   */
  void
  initialize_restriction();

  /**
   * 这些是计算积分时乘以#generalized_face_support_points中的一个函数的系数。它们的组织方式是，每个广义面支持点有一行，面的每个自由度有一列。    更多信息请参见 @ref GlossGeneralizedSupport "广义支持点词汇表条目"
   * 。
   *
   */
  Table<2, double> boundary_weights;

  /**
   * 用于内部自由度内插的预计算系数。此表的原理与#boundary_weights的原理相同。只是，这个表有第三个坐标，用于评估组件的空间方向。
   *
   */
  Table<3, double> interior_weights;

  /**
   * 填充基类中定义的必要表格，如fe.cc中声明的 <code>adjust_quad_dof_index_for_face_orientation_table</code> 。在非标准面、翻转（旋转+180度）或旋转（旋转+90度）的情况下，我们需要用正确的值来填充它。这些是以三个标志的形式给出的（face_orientation, face_flip, face_rotation），见GeometryInfo<dim>中的文档和这个 @ref GlossFaceOrientation "关于面的方向的词汇条"
   * 。    <h3>Example: Raviart-Thomas Elements of order 2 (tensor polynomial
   * degree 3)</h3> 一个面的道夫被连接到一个  $n\times n$
   * 矩阵，这里  <code>n=3</code>
   * 。在我们的例子中，我们可以想象在一个四边形（面）上有以下的道夫。
   * @verbatim
   * ___________
   * |           |
   * |  6  7  8  |
   * |           |
   * |  3  4  5  |
   * |           |
   * |  0  1  2  |
   * |___________|
   * @endverbatim
   * 我们有一个局部 <code>face_dof_index=i+n*j</code> ，索引
   * <code>i</code> in x-direction and index <code>j</code>
   * 在y方向上从0到 <code>n-1</code>.  To extract <code>i</code> 和
   * <code>j</code> we can use <code>i = face_dof_index % n</code> ，<code>j
   * = dof_index / n</code>（整数分割）。然后指数 <code>i</code>
   * 和 <code>j</code> 可以用来计算偏移量。
   * 对于我们的Raviart-Thomas元素的例子，这意味着如果开关是
   * <code>(true | true | true)</code>
   * ，意味着我们首先将面旋转+90度（逆时针），然后再旋转+180度，但我们不翻转它，因为面有标准方向。翻转轴是指从面的左下角到右上角的对角线。有了这些标志，上面的配置就变成了。
   * @verbatim
   * ___________
   * |           |
   * |  2  5  8  |
   * |           |
   * |  1  4  7  |
   * |           |
   * |  0  3  6  |
   * |___________|
   * @endverbatim
   * 请注意，排列组合的必要性取决于这三个标志的组合。
   * 还有一种模式是，被烫的形状函数的符号变化取决于开关的组合。在上面的例子中，它将是
   * @verbatim
   * ___________
   * |           |
   * |  +
   *
   *
   *
   *
   *
   * -  +  |
   * |           |
   * |  +
   *
   *
   *
   *
   *
   * -  +  |
   * |           |
   * |  +
   *
   *
   *
   *
   *
   * -  +  |
   * |___________|
   * @endverbatim
   * 符号变化的相关表格在FE_PolyTensor中声明。
   *
   */
  void
  initialize_quad_dof_index_permutation_and_sign_change();

  // Allow access from other dimensions.
  template <int dim1>
  friend class FE_RaviartThomas;
};



/**
 * Raviart-Thomas元素的节点函数定义为高斯点的点值。
 * <h3>Description of node values</h3>
 * 对于这个Raviart-Thomas元素，节点值不是相对于某些多项式的单元和面矩，而是正交点中的值。按照自由度编号的一般方案，根据单元格边缘的自然排序，边缘上的节点值排在前面，逐个边缘排列。内部自由度排在最后。
 * 对于一个度数为<i>k</i>的RT元素，我们在每个面上选择<i>(k+1)<sup>d-1</sup></i>个高斯点。这些点在面的方向上是按字母顺序排列的。这样一来，处于<i>Q<sub>k</sub></i>的法线分量就被唯一地确定了。此外，由于这个高斯公式在<i>Q<sub>2k+1</sub></i>上是精确的，这些节点值对应于RT-空间的精确积分矩。
 * 在单元内部，矩是相对于各向异性的<i>Q<sub>k</sub></i>空间而言的，其中测试函数在对应于所考虑的矢量分量的方向上低一度。这是通过使用各向异性的高斯公式进行积分来模仿的。
 * @todo
 * 目前的实现只针对笛卡尔网格。你必须使用MappingCartesian。
 * @todo
 * 即使这个元素是为二维和三维空间实现的，节点值的定义也依赖于三维中一致方向的面。因此，在复杂的网格上应该注意。
 *
 *
 * @note  存储在成员变量 FiniteElementData<dim>::degree
 * 中的度数比构造函数参数高一!
 *
 *
 */
template <int dim>
class FE_RaviartThomasNodal : public FE_PolyTensor<dim>
{
public:
  /**
   * 程度为 @p p. 的Raviart-Thomas元素的构造函数。
   *
   */
  FE_RaviartThomasNodal(const unsigned int p);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_RaviartThomasNodal<dim>(degree)</tt>，其中
   * @p dim 和 @p 度被适当的值取代。
   *
   */
  virtual std::string
  get_name() const override;

  // documentation inherited from the base class
  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim> &source,
                                FullMatrix<double> &      matrix,
                                const unsigned int face_no = 0) const override;

  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim> &source,
    const unsigned int        subface,
    FullMatrix<double> &      matrix,
    const unsigned int        face_no = 0) const override;
  virtual bool
  hp_constraints_are_implemented() const override;

  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(const FiniteElement<dim> &fe_other) const override;

  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(const FiniteElement<dim> &fe_other) const override;

  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim> &fe_other,
                         const unsigned int        face_no = 0) const override;

  /**
   * @copydoc   FiniteElement::compare_for_domination() 。
   *
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim> &fe_other,
                         const unsigned int codim = 0) const override final;

private:
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
   * 计算用于传递给基类构造函数的 @p restriction_is_additive
   * 字段的向量。
   *
   */
  static std::vector<bool>
  get_ria_vector(const unsigned int degree);

  /**
   * 如果形状函数 @p shape_index
   * 在面的某处有非零的函数值，这个函数返回 @p true,
   * 。现在，这只在一维的RT0中实现。否则，总是返回  @p
   * true.  。
   *
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  /**
   * 初始化 FiniteElement<dim>::generalized_support_points 和 FiniteElement<dim>::generalized_face_support_points 字段。从构造函数中调用。    更多信息请参见 @ref GlossGeneralizedSupport "关于广义支持点的词汇表条目"
   * 。
   *
   */
  void
  initialize_support_points(const unsigned int rt_degree);

  /**
   * 初始化包络模式和符号变化模式。
   *
   */
  void
  initialize_quad_dof_index_permutation_and_sign_change();
};


 /*@}*/ 

 /* -------------- declaration of explicit specializations ------------- */ 

#ifndef DOXYGEN

template <>
void
FE_RaviartThomas<1>::initialize_restriction();

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


