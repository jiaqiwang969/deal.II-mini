//include/deal.II-translator/fe/fe_dgp_nonparametric_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2020 by the deal.II authors
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

#ifndef dealii_fe_dgp_nonparametric_h
#define dealii_fe_dgp_nonparametric_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * 在映射的正交点上评估的不连续有限元。
 * 警告：这个类还不能正常工作。不要使用它!
 * 这个有限元实现了完整的多项式空间，即 $d$
 * -维多项式的阶数 $k$  。
 * 多项式是没有映射的。因此，它们在任何网格单元上都是常数、线性、二次等。
 * 由于多项式是在实际网格单元的正交点上评估的，所以没有网格转移和内插矩阵。
 * 这个类的目的是实验性的，因此实现起来将是不完整的。
 * 此外，该类没有对一维的情况（<tt>spacedim !=
 * dim</tt>）进行实现。
 *
 *  <h3>Visualization of shape functions</h3>
 * 在2d中，这个元素的形状函数看起来如下。 <h4>$P_0$
 * element</h4> <table> <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P1/P1_DGPNonparametric_shape0000.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $P_0$ element,
 * shape function 0 </td>
 *
 * <td align="center"></tr> </table> <h4>$P_1$ element</h4> <table> <tr> <td
 * align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P1/P1_DGPNonparametric_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P1/P1_DGPNonparametric_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_1$ element, shape function 0 </td>
 *
 * <td align="center"> $P_1$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P1/P1_DGPNonparametric_shape0002.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $P_1$ element,
 * shape function 2 </td>
 *
 * <td align="center"></td> </tr> </table>   <h4>$P_2$ element</h4> <table>
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P2/P2_DGPNonparametric_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P2/P2_DGPNonparametric_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_2$ element, shape function 0 </td>
 *
 * <td align="center"> $P_2$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P2/P2_DGPNonparametric_shape0002.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P2/P2_DGPNonparametric_shape0003.png
 * </td> </tr> <tr> <td align="center"> $P_2$ element, shape function 2 </td>
 *
 * <td align="center"> $P_2$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P2/P2_DGPNonparametric_shape0004.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P2/P2_DGPNonparametric_shape0005.png
 * </td> </tr> <tr> <td align="center"> $P_2$ element, shape function 4 </td>
 *
 * <td align="center"> $P_2$ element, shape function 5 </td> </tr> </table>
 *
 *  <h4>$P_3$ element</h4> <table> <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 0 </td>
 *
 * <td align="center"> $P_3$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0002.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0003.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 2 </td>
 *
 * <td align="center"> $P_3$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0004.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0005.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 4 </td>
 *
 * <td align="center"> $P_3$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0006.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0007.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 6 </td>
 *
 * <td align="center"> $P_3$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0008.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P3/P3_DGPNonparametric_shape0009.png
 * </td> </tr> <tr> <td align="center"> $P_3$ element, shape function 8 </td>
 *
 * <td align="center"> $P_3$ element, shape function 9 </td> </tr> </table>
 * <h4>$P_4$ element</h4>  <table> <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0001.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 0 </td>
 *
 * <td align="center"> $P_4$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0002.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0003.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 2 </td>
 *
 * <td align="center"> $P_4$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0004.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0005.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 4 </td>
 *
 * <td align="center"> $P_4$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0006.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0007.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 6 </td>
 *
 * <td align="center"> $P_4$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0008.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0009.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 8 </td>
 *
 * <td align="center"> $P_4$ element, shape function 9 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0010.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0011.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 10 </td>
 *
 * <td align="center"> $P_4$ element, shape function 11 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0012.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0013.png
 * </td> </tr> <tr> <td align="center"> $P_4$ element, shape function 12 </td>
 *
 * <td align="center"> $P_4$ element, shape function 13 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/DGPNonparametric/P4/P4_DGPNonparametric_shape0014.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $P_4$ element,
 * shape function 14 </td>
 *
 * <td align="center"></td> </tr> </table> 。   <h3> Implementation details
 * </h3>
 * 这个元素没有InternalData类，与所有其他元素不同，因为InternalData类是用来存储那些可以计算一次并多次重复使用的东西（比如参考单元上正交点的形状函数值）。然而，由于该元素没有被映射，这个元素没有任何可以在参考单元上计算的东西
 *
 * - 一切都需要在真实单元上计算
 *
 * 因此，在这样一个对象中没有我们想要存储的东西。因此，我们可以简单地使用
 * FiniteElement::InternalDataBase
 * 已经提供的成员，而不用在这个类的派生类中添加任何东西。
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_DGPNonparametric : public FiniteElement<dim, spacedim>
{
public:
  /**
   * 度数为 @p k. 的张量乘积多项式的构造函数。
   *
   */
  FE_DGPNonparametric(const unsigned int k);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_DGPNonparametric<dim>(degree)</tt>，其中
   * @p dim 和 @p 度数由适当的值代替。
   *
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  // for documentation, see the FiniteElement base class
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  /**
   * 这个函数旨在返回参考单元格上某一点的形状函数的值。然而，由于当前元素没有通过参考单元的映射来实现形状函数，参考单元上不存在形状函数。
   * 因此，正如基类中相应的函数所讨论的，
   * FiniteElement::shape_value(), 这个函数抛出了一个
   * FiniteElement::ExcUnitShapeValuesDoNotExist. 类型的异常。
   *
   */
  virtual double
  shape_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 这个函数旨在返回参考单元格上某一点的形状函数的值。然而，由于当前元素没有通过参考单元的映射来实现形状函数，参考单元上不存在形状函数。
   * 因此，正如基类中相应的函数所讨论的，
   * FiniteElement::shape_value_component(), 这个函数抛出了一个
   * FiniteElement::ExcUnitShapeValuesDoNotExist. 类型的异常。
   *
   */
  virtual double
  shape_value_component(const unsigned int i,
                        const Point<dim> & p,
                        const unsigned int component) const override;

  /**
   * 这个函数的目的是返回一个形状函数在参考单元格上某一点的梯度。然而，由于当前元素没有通过参考单元的映射来实现形状函数，参考单元上不存在形状函数。
   * 因此，正如基类中相应的函数所讨论的，
   * FiniteElement::shape_grad(), 这个函数抛出了一个
   * FiniteElement::ExcUnitShapeValuesDoNotExist. 类型的异常。
   *
   */
  virtual Tensor<1, dim>
  shape_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 这个函数的目的是返回一个形状函数在参考单元格上某一点的梯度。然而，由于当前元素没有通过参考单元的映射来实现形状函数，参考单元上不存在形状函数。
   * 因此，正如基类中相应的函数所讨论的，
   * FiniteElement::shape_grad_component(), 这个函数抛出了一个
   * FiniteElement::ExcUnitShapeValuesDoNotExist. 类型的异常。
   *
   */
  virtual Tensor<1, dim>
  shape_grad_component(const unsigned int i,
                       const Point<dim> & p,
                       const unsigned int component) const override;

  /**
   * 这个函数的目的是返回参考单元上某一点的形状函数的Hessian。然而，由于当前元素没有通过参考单元的映射来实现形状函数，参考单元上不存在形状函数。
   * 因此，正如基类中相应的函数所讨论的，
   * FiniteElement::shape_grad_grad(), 这个函数抛出了一个
   * FiniteElement::ExcUnitShapeValuesDoNotExist. 类型的异常。
   *
   */
  virtual Tensor<2, dim>
  shape_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 这个函数的目的是返回参考单元上某一点的形状函数的Hessian。然而，由于当前元素没有通过参考单元的映射来实现形状函数，参考单元上不存在形状函数。
   * 因此，正如基类中相应的函数所讨论的，
   * FiniteElement::shape_grad_grad_component(), 这个函数抛出了一个
   * FiniteElement::ExcUnitShapeValuesDoNotExist. 类型的异常。
   *
   */
  virtual Tensor<2, dim>
  shape_grad_grad_component(const unsigned int i,
                            const Point<dim> & p,
                            const unsigned int component) const override;

  /**
   * 返回该有限元的多项式程度，即传递给构造函数的值。
   *
   */
  unsigned int
  get_degree() const;

  /**
   * 返回从一个元素的一个面插值到相邻元素的面的矩阵。矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。
   * 衍生元素将不得不实现这个函数。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现从一个给定元素的内插，那么他们必须抛出一个类型为
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented. 的异常。
   *
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double> &                matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * 返回从一个元素的一个面插值到邻近元素的面的矩阵。矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。
   * 衍生元素将不得不实现这个函数。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现从一个给定元素的内插，那么他们必须抛出一个类型为
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented. 的异常。
   *
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source,
    const unsigned int                  subface,
    FullMatrix<double> &                matrix,
    const unsigned int                  face_no = 0) const override;

  /**
   * @name  支持 hp  
     * @{  的函数。
   *
   */

  /**
   * 如果在一个顶点上，有几个有限元处于活动状态，hp代码首先为这些FEs的每个自由度分配不同的全局索引。然后调用这个函数来找出其中哪些应该得到相同的值，从而可以得到相同的全局自由度指数。
   * 因此，该函数返回当前有限元对象的自由度与 @p fe_other,
   * 的自由度之间的相同性列表，后者是对代表在该特定顶点上活动的其他有限元之一的有限元对象的引用。该函数计算两个有限元对象的哪些自由度是等价的，这两个自由度的编号都在0和两个有限元的n_dofs_per_vertex()的相应值之间。每一对的第一个索引表示本元素的一个顶点自由度，而第二个是另一个有限元素的相应索引。
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
   * 对于FE_DGPNonparametric类，结果总是真（与元素的程度无关），因为它没有悬挂节点（是一个不连续的元素）。
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
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   * 这个函数是虚拟的，因为有限元对象通常是通过指向其基类的指针来访问的，而不是类本身。
   *
   */
  virtual std::size_t
  memory_consumption() const override;

protected:
  /**
   * 准备内部数据结构并填入独立于单元的值。
   *
   */
  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim> &       quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  using FiniteElement<dim, spacedim>::fill_fe_face_values;

  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const hp::QCollection<dim - 1> &                            quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          sub_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

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
   * 指向代表这里使用的多项式空间的对象的指针。
   *
   */
  const PolynomialSpace<dim> polynomial_space;

  // Allow access from other dimensions.
  template <int, int>
  friend class FE_DGPNonparametric;
};

 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


