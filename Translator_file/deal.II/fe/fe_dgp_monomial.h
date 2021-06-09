//include/deal.II-translator/fe/fe_dgp_monomial_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2020 by the deal.II authors
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

#ifndef dealii_fe_dgp_monomial_h
#define dealii_fe_dgp_monomial_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_p.h>

#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * 基于单项式的不连续有限元。
 * 这种有限元实现了完整的多项式空间，即p度的dim-维多项式。例如，在2d中，元素FE_DGP(1)将代表函数
 * $\{1,\hat x,\hat y\}$ 的跨度，这与由 $\{1,\hat x,\hat y,\hat x\hat
 * y\}$
 * 的跨度形成的元素FE_DGQ(1)相反。由于DGP空间中每个四边形只有三个未知数，因此立即可以看出这个元素不可能是连续的。
 * 这个元素的基函数被选择为上面列出的单项式。请注意，这是与FE_DGP类的主要区别，FE_DGP类使用一组完整度数
 * <code>p</code>
 * 的多项式，在单位平方上形成Legendre基。因此，在那里，如果网格单元是平行四边形，质量矩阵是对角线的。这里的基不具有这种特性；然而，它的计算更简单。另一方面，这个元素有一个额外的缺点，即局部单元矩阵通常比源自FE_DGP元素的矩阵有更差的条件数。
 * 这类元素没有在二维一的情况下实现（<tt>spacedim !=
 * dim</tt>）。 <h3>Transformation properties</h3>
 * 值得注意的是，在（双，三）线性映射下，这个元素描述的空间不包含
 * $P(k)$  ，即使我们使用度数为 $k$
 * 的多项式基。因此，例如，在具有非affine单元的网格上，线性函数不能由FE_DGP(1)或FE_DGPMonomial(1)类型的元素准确表示。
* 这可以通过下面的二维例子来理解：考虑顶点在 $(0,0),(1,0),(0,1),(s,s)$ 的单元： @image html dgp_doesnt_contain_p.png  的单元。
 * 对于这个单元，双线性变换 $F$ 产生的关系 $x=\hat x+\hat
 * x\hat y$ 和 $y=\hat y+\hat x\hat y$ 将参考坐标 $\hat x,\hat y$
 * 和实空间坐标 $x,y$
 * 关联起来。在这种映射下，常数函数显然被映射到它自己，但
 * $P_1$ 空间的另外两个形状函数，即 $\phi_1(\hat x,\hat y)=\hat
 * x$ 和 $\phi_2(\hat x,\hat y)=\hat y$ 被映射到
 * $\phi_1(x,y)=\frac{x-t}{t(s-1)},\phi_2(x,y)=t$ ，其中
 * $t=\frac{y}{s-x+sx+y-sy}$  。 对于 $s=1$
 * 这种简单情况，即如果实心单元是单位平方，表达式可以简化为
 * $t=y$ 和 $\phi_1(x,y)=x,\phi_2(x,y)=y$
 * 。然而，对于所有其他情况，函数 $\phi_1(x,y),\phi_2(x,y)$
 * 不再是线性的，也不是它们的任何线性组合。因此，线性函数不在映射的
 * $P_1$ 多项式的范围内。
 *
 *  <h3>Visualization of shape functions</h3>
 * 在2d中，这个元素的形状函数看起来如下。 <h4>$P_0$
 * element</h4> <table> <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P1/P1_DGPMonomial_shape0000.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $P_0$ element,
 * shape function 0 </td>
 *
 * <td align="center"></tr> </table> <h4>$P_1$ element</h4> <table> <tr> <td
 * align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P1/P1_DGPMonomial_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P1/P1_DGPMonomial_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_1$ element, shape function 0 </td>
 *
 * <td align="center"> $P_1$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P1/P1_DGPMonomial_shape0002.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $P_1$ element,
 * shape function 2 </td>
 *
 * <td align="center"></td> </tr> </table>
 *
 *  <h4>$P_2$ element</h4> <table> <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P2/P2_DGPMonomial_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P2/P2_DGPMonomial_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_2$ element, shape function 0 </td>
 *
 * <td align="center"> $P_2$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P2/P2_DGPMonomial_shape0002.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P2/P2_DGPMonomial_shape0003.png
 * </td> </tr> <tr> <td align="center"> $P_2$ element, shape function 2 </td>
 *
 * <td align="center"> $P_2$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P2/P2_DGPMonomial_shape0004.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P2/P2_DGPMonomial_shape0005.png
 * </td> </tr> <tr> <td align="center"> $P_2$ element, shape function 4 </td>
 *
 * <td align="center"> $P_2$ element, shape function 5 </td> </tr> </table>
 * <h4>$P_3$ element</h4> <table> <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 0 </td>
 *
 * <td align="center"> $P_3$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0002.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0003.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 2 </td>
 *
 * <td align="center"> $P_3$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0004.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0005.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 4 </td>
 *
 * <td align="center"> $P_3$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0006.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0007.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 6 </td>
 *
 * <td align="center"> $P_3$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0008.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P3/P3_DGPMonomial_shape0009.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 8 </td>
 *
 * <td align="center"> $P_3$ element, shape function 9 </td> </tr> </table>
 * <h4>$P_4$ element</h4>  <table> <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 0 </td>
 *
 * <td align="center"> $P_4$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0002.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0003.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 2 </td>
 *
 * <td align="center"> $P_4$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0004.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0005.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 4 </td>
 *
 * <td align="center"> $P_4$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0006.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0007.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 6 </td>
 *
 * <td align="center"> $P_4$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0008.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0009.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 8 </td>
 *
 * <td align="center"> $P_4$ element, shape function 9 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0010.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0011.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 10 </td>
 *
 * <td align="center"> $P_4$ element, shape function 11 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0012.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0013.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 12 </td>
 *
 * <td align="center"> $P_4$ element, shape function 13 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPMonomial/P4/P4_DGPMonomial_shape0014.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $P_4$ element,
 * shape function 14 </td>
 *
 * <td align="center"></td> </tr> </table> .
 *
 */
template <int dim>
class FE_DGPMonomial : public FE_Poly<dim>
{
public:
  /**
   * 度数为<tt>p</tt>的多项式空间的构造函数。
   *
   */
  FE_DGPMonomial(const unsigned int p);

  /**
   * 返回一个唯一标识有限元的字符串。这个类返回<tt>FE_DGPMonomial<dim>(degree)</tt>，其中<tt>dim</tt>和<tt>p</tt>被适当的值替换。
   *
   */
  virtual std::string
  get_name() const override;

  /**
   * @name  支持hp的函数  
     * @{ 
   *
   */

  /**
   * 如果在一个顶点上，有几个有限元处于活动状态，hp代码首先为这些FEs的每个自由度分配不同的全局索引。然后调用这个函数来找出其中哪些应该得到相同的值，从而可以得到相同的全局自由度指数。
   * 因此，该函数返回当前有限元对象的自由度与 @p fe_other,
   * 的自由度之间的相同性列表，后者是对代表在该特定顶点上活动的其他有限元之一的有限元对象的引用。该函数计算两个有限元对象的哪些自由度是等价的，两个自由度的编号都在零和两个有限元的n_dofs_per_vertex()的相应值之间。每一对的第一个索引表示本元素的一个顶点自由度，而第二个是另一个有限元素的相应索引。
   * 作为一个不连续的元素，这种约束的集合当然是空的。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(const FiniteElement<dim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理的是线上的自由度。
   * 这是一个不连续的元素，这种约束的集合当然是空的。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(const FiniteElement<dim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理的是四边形上的自由度。
   * 这是一个不连续的元素，这种约束的集合当然是空的。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim> &fe_other,
                         const unsigned int        face_no = 0) const override;

  /**
   * 返回该元素是否以新的方式实现其悬挂节点约束，这必须用于使元素
   * "hp-compatible"。
   * 对于FE_DGPMonomial类，结果总是真（与元素的程度无关），因为它没有悬挂节点（是一个不连续的元素）。
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
   * @}
   *
   */

  /**
   * 返回从给定的有限元插值到现在的矩阵。然后矩阵的大小是
   * @p dofs_per_cell 乘以<tt>source.n_dofs_per_cell()</tt>。
   * 这些矩阵只有在源元素也是 @p FE_Q
   * 元素时才可用。否则，会抛出一个
   * FiniteElement<dim>::ExcInterpolationNotImplemented 类型的异常。
   *
   */
  virtual void
  get_interpolation_matrix(const FiniteElement<dim> &source,
                           FullMatrix<double> &      matrix) const override;

  /**
   * 返回从一个元素的一个面插值到邻近元素的面的矩阵。矩阵的大小是
   * @p  dofs_per_face乘以<tt>source.dofs_per_face</tt>。
   * 衍生元素将不得不实现这个功能。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现从一个给定元素的内插，那么他们必须抛出一个类型为
   * FiniteElement<dim>::ExcInterpolationNotImplemented. 的异常。
   *
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim> &source,
                                FullMatrix<double> &      matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * 返回从一个元素的一个面插值到邻近元素的面的矩阵。矩阵的大小是
   * @p  dofs_per_face乘以<tt>source.dofs_per_face</tt>。
   * 衍生元素将不得不实现这个功能。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现从一个给定元素的内插，那么他们必须抛出一个类型为
   * FiniteElement<dim>::ExcInterpolationNotImplemented. 的异常。
   *
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim> &source,
    const unsigned int        subface,
    FullMatrix<double> &      matrix,
    const unsigned int        face_no = 0) const override;

  /**
   * 如果形状函数 @p shape_index 在面 @p face_index.
   * 的某处有非零函数值，该函数返回 @p true, 。
   *
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   * 这个函数是虚拟的，因为有限元对象通常是通过指向其基类的指针来访问的，而不是类本身。
   *
   */
  virtual std::size_t
  memory_consumption() const override;

  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

private:
  /**
   * 仅供内部使用。其全称是 @p get_dofs_per_object_vector
   * 函数，它创建了 @p dofs_per_object
   * 向量，在构造函数内需要传递给 @p
   * FiniteElementData的构造函数。
   *
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree);

  /**
   * 初始化限制矩阵。从构造函数中调用。
   *
   */
  void
  initialize_restriction();
};

 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


