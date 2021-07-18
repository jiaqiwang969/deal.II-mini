//include/deal.II-translator/fe/fe_dgq_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
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

#ifndef dealii_fe_dgq_h
#define dealii_fe_dgq_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class MappingQ;
template <int dim>
class Quadrature;
#endif

 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * 基于等距支持点的标量、不连续张量积元素的实现。
 * 这是一个基于拉格朗日多项式的张量乘积的不连续有限元。形状函数是单元格上点的等距网格的拉格朗日内插。这些点按列序编号，<i>x</i>最快，然后是<i>y</i>，然后是<i>z</i>（如果这些坐标在一个给定的空间维度上都存在的话）。例如，这些是3D中<tt>FE_DGQ(1)</tt>的节点排序。
 * @verbatim
 *       6-------7        6-------7
 *      /|       |       /       /|
 *     / |       |      /       / |
 *    /  |       |     /       /  |
 *   4   |       |    4-------5   |
 *   |   2-------3    |       |   3
 *   |  /       /     |       |  /
 *   | /       /      |       | /
 *   |/       /       |       |/
 *   0-------1        0-------1
 * @endverbatim
 * 和<tt>FE_DGQ(2)</tt>。
 * @verbatim
 *       24--25--26       24--25--26
 *      /|       |       /       /|
 *    21 |       |     21  22  23 |
 *    /  15  16  17    /       /  17
 *  18   |       |   18--19--20   |
 *   |12 6---7---8    |       |14 8
 *   9  /       /     9  10  11  /
 *   | 3   4   5      |       | 5
 *   |/       /       |       |/
 *   0---1---2        0---1---2
 * @endverbatim
 * 节点13被放置在六角的内部。
 * 然而，请注意，这些只是形状函数的拉格朗日插值点。即使它们在物理上可能在单元的边界上，但在逻辑上它们是在内部的，因为这些形状函数没有跨越单元边界的连续性要求。虽然是不连续的，但当限制在一个单元时，这个元素的形状函数与FE_Q元素的形状函数完全一样，它们在视觉上显示出来。
 * <h3>Unit support point distribution and conditioning of interpolation</h3>
 * 当以多项式一度或二度构造FE_DGQ元素时，在0和1(线性情况)或0、0.5和1(二次情况)的等距支持点被使用。单位支持点或节点点<i>x<sub>i</sub></i>是那些<i>j</i>个拉格朗日多项式满足
 * $\delta_{ij}$ 属性的点，即一个多项式为1，其他都是0。
 * 对于更高的多项式度数，支持点默认为非流动性的，并选择为<tt>（度数+1）</tt>阶Gauss-Lobatto正交规则的支持点。这种点分布在任意多项式度数下产生条件良好的Lagrange插值。相比之下，基于等距点的多项式随着多项式度数的增加而变得越来越没有条件。在内插法中，这种效应被称为Runge现象。对于Galerkin方法，Runge现象通常在解的质量上不明显，而是在相关系统矩阵的条件数上。例如，10度的等距点的元素质量矩阵的条件数为2.6e6，而Gauss-Lobatto点的条件数约为400。
 * 一维的Gauss-Lobatto点包括单位区间的端点0和+1。内部点被移向端点，这使得靠近元素边界的点分布更加密集。
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_DGQ : public FE_Poly<dim, spacedim>
{
public:
  /**
   * 度数为<tt>p</tt>的张量乘积多项式的构造器。使用该构造器创建的形状函数对应于Gauss-Lobatto支持（节点）点在每个坐标方向的Lagrange插值多项式。
   *
   */
  FE_DGQ(const unsigned int p);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_DGQ<dim>(degree)</tt>，其中<tt>dim</tt>和<tt>degree</tt>由适当的值替换。
   *
   */
  virtual std::string
  get_name() const override;

  /**
   * 返回从给定的有限元插值到现在的矩阵。然后矩阵的大小是
   * @p dofs_per_cell 乘以<tt>source.n_dofs_per_cell()</tt>。
   * 这些矩阵只有在源元素也是 @p
   * FE_DGQ元素时才能使用。否则，会抛出一个
   * FiniteElement<dim>::ExcInterpolationNotImplemented 类型的异常。
   *
   */
  virtual void
  get_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                           FullMatrix<double> &matrix) const override;

  /**
   * 返回从一个元素的一个面插值到相邻元素的面的矩阵。矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。
   * 衍生元素将不得不实现这个函数。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现从一个给定元素的内插，那么他们必须抛出一个类型为
   * FiniteElement<dim>::ExcInterpolationNotImplemented. 的异常。
   *
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double> &                matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * 返回从一个元素的一个面插值到邻近元素的面的矩阵。矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。
   * 衍生元素将不得不实现这个函数。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现从一个给定元素的内插，那么他们必须抛出一个类型为
   * FiniteElement<dim>::ExcInterpolationNotImplemented. 的异常。
   *
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source,
    const unsigned int                  subface,
    FullMatrix<double> &                matrix,
    const unsigned int                  face_no = 0) const override;

  /**
   * 从精细网格空间投射到粗略网格空间。重写FiniteElement中的相应方法，实现懒人评估（在请求时初始化）。
   * 如果这个投影运算符与一个矩阵 @p P,
   * 相关联，那么这里将返回这个矩阵 @p P_i
   * 对单个子单元的限制。    矩阵 @p P 是单元格矩阵 @p
   * P_i的串联或相加，取决于#restriction_is_additive_flags。这区分了插值（连接）和标量积（求和）方面的投影。
   * 行和列指数分别与粗网格和细网格空间有关，与相关运算符的定义一致。
   *
   */
  virtual const FullMatrix<double> &
  get_restriction_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * 网格间的嵌入矩阵。重写FiniteElement中的相应方法，实现懒人评估（查询时初始化）。
   * 从粗网格空间到细网格空间的身份运算符与一个矩阵 @p
   * P. 相关联，该矩阵 @p P_i
   * 对单个子单元的限制在此返回。    矩阵 @p P
   * 是串联的，而不是单元格矩阵 @p
   * P_i的总和。也就是说，如果同一个非零条目<tt>j,k</tt>存在于两个不同的子矩阵
   * @p P_i,
   * 中，其值在两个矩阵中应该是相同的，它只被复制到矩阵
   * @p P 中一次。
   * 行和列指数分别与细格和粗格空间相关，与相关运算符的定义一致。
   * 这些矩阵被组装多层次方法的延长矩阵的程序所使用。
   * 在使用这个矩阵阵列组装单元格之间的转移矩阵时，延长矩阵中的零元素被丢弃，不会填满转移矩阵。
   *
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * @name  支持hp的函数  
     * @{ 
   *
   */

  /**
   * 如果在一个顶点上，有几个有限元处于活动状态，hp代码首先为这些FEs的每个自由度分配不同的全局索引。然后调用这个函数来找出其中哪些应该得到相同的值，从而可以得到相同的全局自由度指数。
   * 因此，该函数返回当前有限元对象的自由度与 @p fe_other,
   * 的自由度之间的相同性列表，该列表是对代表该特定顶点上活动的其他有限元之一的有限元对象的引用。该函数计算两个有限元对象的哪些自由度是等价的，这两个自由度的编号都在0和两个有限元的n_dofs_per_vertex()的相应值之间。每一对的第一个索引表示本元素的一个顶点自由度，而第二个是另一个有限元素的相应索引。
   * 作为一个不连续的元素，这种约束的集合当然是空的。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理的是线上的自由度。
   * 这是一个不连续的元素，这种约束的集合当然是空的。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理的是四边形上的自由度。
   * 这是一个不连续的元素，这种约束的集合当然是空的。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  /**
   * 返回该元素是否以新的方式实现其悬挂节点约束，这必须用于使元素
   * "hp-compatible"。
   * 对于FE_DGQ类，结果总是真（与元素的程度无关），因为它没有悬挂节点（是一个不连续的元素）。
   *
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * @copydoc   FiniteElement::compare_for_domination() .
   *
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

  /**
   * @}
   *
   */

  /**
   * 如果形状函数 @p shape_index
   * 在面的某处有非零函数值，该函数返回 @p true, 。
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
   * 实现FiniteElement类中的相应函数。
   * 因为当前元素是插值的，所以节点值正好是支持点值。此外，由于当前元素是标量的，支持点的值需要是长度为1的向量。
   *
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

protected:
  /**
   * 张量积多项式的构造函数，基于一个任意的多项式向量。这个构造函数在派生类中用于构造例如具有任意结点的元素或基于Legendre多项式的元素。
   * 这些多项式的程度是<tt>polynomials.size()-1</tt>。
   *
   */
  FE_DGQ(const std::vector<Polynomials::Polynomial<double>> &polynomials);

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
   * 计算自由度旋转的重新编号。
   * 这个函数将自由度的张量积编号旋转90度。
   * 它用于计算子体的转移矩阵，只使用第一个子体的矩阵。
   * 方向参数决定了旋转的类型。它是 @p xXyYzZ.
   * 的一个字符，该字符决定了旋转的轴，大小写决定了方向。小写的是在轴的方向上逆时针看。
   * 由于绕Y轴的旋转不被使用，所以也不被实现。
   *
   */
  void
  rotate_indices(std::vector<unsigned int> &indices,
                 const char                 direction) const;

  /* 用于保护限制和嵌入矩阵的初始化的互斥。 
*
*/
  mutable Threads::Mutex mutex;

  // Allow access from other dimensions.
  template <int dim1, int spacedim1>
  friend class FE_DGQ;

  // Allow @p MappingQ class to access to build_renumbering function.
  template <int dim1, int spacedim1>
  friend class MappingQ;
};



/**
 * 实现基于任意节点的拉格朗日多项式的标量、不连续张量乘积元素。这个类的主要目的是提供一个元素，对于这个元素，质量矩阵可以通过选择基函数使其成为对角线，这些基函数在单元的顶点不是零就是一，而是在一组给定的正交点是零或一。如果这组正交点也被用于质量矩阵的积分，那么它将是对角线的。正交点的数量自动决定了这个元素所选择的多项式程度。典型的应用是高斯正交或高斯-洛巴托正交（通过基类提供）。
 * 详见FE_DGQ中的基类文档。
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_DGQArbitraryNodes : public FE_DGQ<dim, spacedim>
{
public:
  /**
   * 基于 Polynomials::Lagrange
   * 正交规则<tt>点</tt>中支持点插值的张量乘积多项式的构造器。这些多项式的程度是<tt>points.size()-1</tt>。
   *
   */
  FE_DGQArbitraryNodes(const Quadrature<1> &points);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_DGQArbitraryNodes<dim>(degree)</tt>，其中<tt>dim</tt>和<tt>degree</tt>被适当的值替换。
   *
   */
  virtual std::string
  get_name() const override;

  /**
   * 在FiniteElement类中实现相应的函数。
   * 由于当前元素是插值的，所以节点值正好是支持点的值。此外，由于当前元素是标量的，支持点的值需要是长度为1的向量。
   *
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;
  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;
};



/**
 * 基于Legendre多项式的标量、不连续张量积元素的实现，由多项式空间的张量积描述
 * Polynomials::Legendre.
 * 张量积是使用TensorProductPolynomials实现的，形状函数的排序与TensorProductPolynomials中一样，是lexicographic。例如，当<tt>度=n</tt>时，2d中的排序是
 * $P_0(x)P_0(y),\ P_1(x)P_0(y),\ \ldots,\ P_n(x)P_0(y),\ P_0(x)P_1(y), \
 * \ldots,\ P_n(x)P_1(y),\ \ldots,\ P_0(x)P_n(y),\ \ldots,\ P_n(x)P_n(y)$
 * ，其中 $\{P_i\}_{i=0}^{n}$ 是定义在 $[0,1]$
 * 上的一维Legendre多项式。与基本的FE_DGQ元素相反，这些元素不是插值的，没有定义支持点。
 * 详见FE_DGQ中的基类文档。
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_DGQLegendre : public FE_DGQ<dim, spacedim>
{
public:
  /**
   * 基于 Polynomials::Legendre
   * 插值的张量乘积多项式的构造函数。
   *
   */
  FE_DGQLegendre(const unsigned int degree);

  /**
   * 返回一个元素的常数模式列表。对于Legendre基础，它返回一行，其中第一个元素（对应于常数模式）被设置为真，所有其他元素被设置为假。
   *
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_DGQLegendre<dim>(degree)</tt>，其中<tt>dim</tt>和<tt>degree</tt>分别由模板参数和传递给构造函数的参数给出的值替换。
   *
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;
};



/**
 * 实现基于Hermite-like多项式的标量、不连续张量乘积元素，由多项式空间
 * Polynomials::HermiteLikeInterpolation. 描述
 * 与基本的FE_DGQ元素相反，这些元素不是插值的，没有定义支持点。
 * 请注意，Hermite多项式只适用于大于或等于3度的情况，因此
 * Polynomials::HermiteLikeInterpolation
 * 的有益特性，即每个维度上只有两个基函数在一个面上有非三值或导数，只适用于更高的度。为了方便使用零度到二度，本类构建了一个通常的拉格朗日基。
 * 详见FE_DGQ中的基类文档。
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_DGQHermite : public FE_DGQ<dim, spacedim>
{
public:
  /**
   * 基于 Polynomials::HermiteLikeInterpolation.
   * 的张量乘积多项式的构造器。
   *
   */
  FE_DGQHermite(const unsigned int degree);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_DGQHermite<dim>(degree)</tt>，其中<tt>dim</tt>和<tt>degree</tt>分别由模板参数和传递给构造函数的参数给出的值替换。
   *
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;
};


 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


