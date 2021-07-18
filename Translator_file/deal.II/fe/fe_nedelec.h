//include/deal.II-translator/fe/fe_nedelec_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2021 by the deal.II authors
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

#ifndef dealii_fe_nedelec_h
#define dealii_fe_nedelec_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomials_nedelec.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_poly_tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * @warning
 * 实现的几个方面是实验性的。目前，在全局细化的网格上使用该元素是安全的，其面的方向是一致的。更详细的注意事项见下面的todo条目。
 * N&eacute;d&eacute;lec元素的实现。N&eacute;d&eacute;lec空间的设计是为了解决解只存在于空间
 * $H^\text{curl}=\{ {\mathbf u} \in L_2: \text{curl}\, {\mathbf u} \in L_2\}$
 * 的问题，而不是更常用的空间 $H^1=\{ u \in L_2: \nabla u \in
 * L_2\}$
 * 。换句话说，解决方案必须是一个矢量场，其卷曲是可平方整除的，但其梯度可能不是可平方整除的。这个空间（和这些元素）的典型应用是麦克斯韦方程和相应的简化，例如麦克斯韦方程的简化版本，只涉及电场
 * $\mathbf E$ ，在没有电流时，必须满足方程 $\text{curl}\,
 * \text{curl}\, {\mathbf E} = 0$ ，或者磁矢量势 $\mathbf A$
 * 在时间独立的情况下必须满足方程
 * $\text{curl}\,\text{curl}\,{\mathbf A} = 4\pi{\mathbf j}$ 。
 * $H^\text{curl}$
 * 中的函数的决定性特征是，它们通常是不连续的。
 *
 * 但如果你在2D中画一条线（或在3D中画一个表面），那么矢量场的<i>tangential</i>分量必须在线（或表面）上连续，尽管法线分量可能不是。因此，N&eacute;d&eacute;lec元素的构造是这样的：（i）它是 @ref vector_valued "矢量值"
 * ，（ii）形状函数是不连续的，但（iii）每个形状函数所代表的矢量场的切向分量在单元面上是连续的。
 * N&eacute;d&eacute;lec元素的其他属性是：（i）它是 @ref GlossPrimitive "非原始元素"
 * ；（ii）形状函数被定义为使某些面的积分为零或一，而不是常见的某些点值为零或一的情况。
 * 我们遵循通常使用的
 *
 * --虽然很混乱
 *
 * - 对N&eacute;d&eacute;lec元素的 "度 "的定义。具体来说，元素的 "度 "表示有限元空间中包含的<i>largest complete polynomial
 * subspace</i>的多项式度，即使该空间可能包含更高多项式度的形状函数。因此，最低阶元素是FE_Nedelec(0)，即 "零度 "的Raviart-Thomas元素，尽管这个空间的函数一般是每个变量的一度多项式。这种 "度 "的选择意味着函数本身的近似顺序是<i>degree+1</i>，就像通常的多项式空间一样。如此选择的编号意味着序列@f[
 * Q_{k+1}
 * \stackrel{\text{grad}}{\rightarrow}
 * \text{Nedelec}_k
 * \stackrel{\text{curl}}{\rightarrow}
 * \text{RaviartThomas}_k
 * \stackrel{\text{div}}{\rightarrow}
 * DGQ_{k}
 * @f]注意，这遵循Brezzi和Raviart的惯例，尽管不是N&eacute;d&eacute;lec的原始论文中使用的。
 * 该类没有在二维一的情况下实现（<tt>spacedim != dim</tt>）。
 * @todo
 * 即使这个元素是为二维和三维空间实现的，节点值的定义也依赖于三维中一致方向的面。因此，在复杂的网格上应该注意。
 *
 *  <h3>Interpolation</h3>
 * @ref GlossInterpolation  与N&eacute;d&eacute;lec元素相关的 "插值 "
 * 算子的构造是这样的：插值和计算卷曲是矩形网格单元上的换算操作。我们从插值任意函数以及#限制性矩阵中要求这一点。
 * <h4>Node values</h4>
 * 参考单元上度数为<i>k</i>的元素的 @ref GlossNodes "节点值 "是。  <ol>   <li>  在边上：切向分量相对于<i>k</i>度的多项式的矩。  <li>  在面：切向分量相对于<tt>dim</tt>-1度的FE_Nedelec多项式的矩值<i>k</i>-1。  <li>  在单元格中：相对于度数为<i>k</i>的FE_Q的多项式的梯度的矩。  </ol>
 * <h4>Generalized support points</h4>
 * 上面的节点值依赖于积分，这些积分将由正交规则本身计算出来。广义支持点是一组点，使这种正交能够以足够的精度进行。需要的点是每个边上的QGauss<sub>k+1</sub>和每个面上的QGauss<sub>k+2</sub>以及单元内部的那些点（或者对于N<sub>1</sub>来说没有）。
 *
 *
 */
template <int dim>
class FE_Nedelec : public FE_PolyTensor<dim>
{
public:
  /**
   *
   */
  FE_Nedelec(const unsigned int order);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_Nedelec<dim>(degree)</tt>，
   * @p dim 和 @p degree 用适当的值替换。
   *
   */
  virtual std::string
  get_name() const override;


  /**
   * 如果形状函数 @p shape_index 在面 @p face_index.
   * 的某处有非零函数值，该函数返回 @p true, 。
   *
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  /**
   * 返回这个元素是否以新的方式实现了它的悬挂节点约束，这必须用于使元素
   * "hp-compatible"。
   * 对于<tt>FE_Nedelec</tt>类，结果总是真（与元素的程度无关），因为它实现了hp-capability所需的完整功能集。
   *
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * @copydoc   FiniteElement::compare_for_domination() .
   *
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim> &fe_other,
                         const unsigned int codim = 0) const override final;

  /**
   * 如果在一个顶点上，有几个有限元处于活动状态，hp代码首先为这些FEs的每个自由度分配不同的全局索引。然后调用这个函数来找出其中哪些应该得到相同的值，从而可以得到相同的全局自由度指数。
   * 因此，该函数返回当前有限元对象的自由度与 @p fe_other,
   * 的自由度之间的相同性列表，后者是对代表在该特定顶点上活动的其他有限元之一的有限元对象的引用。该函数计算两个有限元对象的哪些自由度是等价的，这两个自由度的编号都在0和两个有限元的n_dofs_per_vertex()的相应值之间。每一对的第一个索引表示本元素的一个顶点自由度，而第二个是另一个有限元素的相应索引。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(const FiniteElement<dim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理线上自由度。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(const FiniteElement<dim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理线上的自由度。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim> &fe_other,
                         const unsigned int        face_no = 0) const override;

  /**
   * 返回从一个元素的面插值到邻近元素的面的矩阵。矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。
   * 衍生元素将不得不实现这个函数。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现从一个给定元素的内插，那么他们必须抛出一个类型为
   * <tt>FiniteElement<dim>::ExcInterpolationNotImplemented</tt>. 的异常。
   *
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim> &source,
                                FullMatrix<double> &      matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * 返回从一个元素的面内插到邻近元素的子面的矩阵。矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。
   * 衍生元素将不得不实现这个函数。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现给定元素的插值，那么他们必须抛出一个<tt>ExcInterpolationNotImplemented</tt>类型的异常。
   *
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim> &source,
    const unsigned int        subface,
    FullMatrix<double> &      matrix,
    const unsigned int        face_no = 0) const override;

  /**
   * 从精细网格空间投射到粗略网格空间。如果这个投影运算符与一个矩阵
   * @p P, 相关联，那么这里将返回这个矩阵 @p P_i
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
   * 网格间的嵌入矩阵。
   * 从粗网格空间到细网格空间的同一运算符与一个矩阵 @p
   * P. 相关，该矩阵 @p P_i 对单个子单元的限制在此返回。
   * 矩阵 @p P 是串联的，而不是单元格矩阵 @p
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

  // documentation inherited from the base class
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

  /**
   * 返回一个元素的常数模式列表。
   *
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  virtual std::size_t
  memory_consumption() const override;

  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

private:
  /**
   * 仅供内部使用。它的全称是 @p get_dofs_per_object_vector
   * 函数，它创建了 @p dofs_per_object
   * 向量，在构造函数内需要传递给 @p
   * FiniteElementData的构造函数。
   * 如果可选参数<tt>dg</tt>为真，返回的向量将有分配给单元的所有自由度，面和边上没有。
   *
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree, bool dg = false);

  /**
   * 初始化FiniteElement类的 @p generalized_support_points
   * 字段，用插值权重（#boundary_weights 和
   * interior_weights）填充表格。从构造函数中调用。
   *
   */
  void
  initialize_support_points(const unsigned int order);

  /**
   * 初始化从细化网格单元上的函数到父单元的插值。根据Nédélec元素的理念，这个限制算子弱化地保留了一个函数的卷曲。
   *
   */
  void
  initialize_restriction();

  /**
   * 这些是计算积分时乘以#generalized_face_support_points中的一个函数的系数。    更多信息请参见 @ref GlossGeneralizedSupport "广义支持点的术语条目"
   * 。
   *
   */
  Table<2, double> boundary_weights;

  /**
   * 用于保护限制和嵌入矩阵的初始化的Mutex。
   *
   */
  mutable Threads::Mutex mutex;

  /**
   * 初始化置换模式和符号变化模式。
   * @note
   * 这个函数还没有完全充满正确的实现。它需要在未来的版本中统一实现，以便在包含有翻转面的单元格的网格上工作。
   *
   */
  void
  initialize_quad_dof_index_permutation_and_sign_change();

  // Allow access from other dimensions.
  template <int dim1>
  friend class FE_Nedelec;
};

 /* -------------- declaration of explicit specializations ------------- */ 

#ifndef DOXYGEN

template <>
void
FE_Nedelec<1>::initialize_restriction();

#endif // DOXYGEN

 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


