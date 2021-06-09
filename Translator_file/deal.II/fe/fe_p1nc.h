//include/deal.II-translator/fe/fe_p1nc_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2021 by the deal.II authors
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

#ifndef dealii_fe_p1nc_h
#define dealii_fe_p1nc_h

#include <deal.II/base/config.h>

#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/fe/fe.h>


DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * 实现P1非压缩有限元的标量版本，这是一个二维四边形上的片状线性元素。这个实现只适用于二维空间中的二维单元（即维数为0）。
 * 与通常的连续、 $H^1$
 * 符合性有限元不同，P1非符合性元素并不强制要求跨边的连续性。然而，它要求在积分意义上的连续性：空间中的任何函数在两个相邻元素共享的公共边的两边应该有相同的积分值。
 * 因此，不连续元素空间中的每个函数都可以是不连续的，因此不包括在
 * $H^1_0$
 * 中，就像不连续加尔金（DG）有限元空间的基函数一样。另一方面，DG空间中的基函数是完全不连续的，跨越边缘，两边的值之间没有任何关系。
 * 这就是为什么通常的DG方案的弱式计算包含额外的跨边跳跃惩罚项来控制不连续性的原因。
 * 然而，不符合要求的元素通常不需要在它们的弱式计算中加入额外的条款，因为它们沿边的积分从两边都是一样的，也就是说，有<i>some
 * level</i>的连续性。 <h3>Dice Rule</h3>
 * 由于P1非线性空间中的任何函数在每个元素上都是片状线性的，所以每个边的中点的函数值与边上的平均值相同。因此，在这种情况下，跨越每条边的积分值的连续性等同于每条边的中点值的连续性。
 * 因此，对于P1不合格元素，单元格边缘上中点的函数值是很重要的。在四边形上定义（局部）自由度（DoF）的首次尝试是通过使用函数的中点值。
 * 然而，这4个函数并不是线性独立的，因为二维上的线性函数只由3个独立的值唯一决定。一个简单的观察结果是，四边形上的任何线性函数都应该满足'骰子规则'：一个单元格对边的中点上的两个函数值之和等于另一个对边中点上的函数值之和。这被称为'骰子规则'，因为骰子对边的点数加起来也总是相同的数字（在骰子的情况下，为7）。
 * 在公式中，骰子规则被写成 $\phi(m_0) + \phi(m_1) = \phi(m_2) +
 * \phi(m_3)$ ，用于函数空间中的所有 $\phi$ ，其中 $m_j$ 是边
 * $e_j$
 * 的中点。在这里，我们假设在deal.II中使用的、在类GeometryInfo中描述的边的标准编号惯例。
 * 反之，如果给定了4个满足骰子规则的中点值，那么总是存在唯一的线性函数，它与4个中点值重合。
 * 由于骰子规则，任何三个中点的三个值都可以决定最后一个中点的最后一个值。这意味着一个单元格上的独立局部函数的数量是3，这也是二维单元格上的线性多项式空间的维度。
 * <h3>Shape functions</h3>
 * 在介绍自由度之前，我们提出一个单元上的4个局部形状函数。由于骰子规则的存在，我们需要对形状函数进行特殊的构造。尽管以下4个形状函数在一个单元内不是线性独立的，但它们有助于定义在整个域上线性独立的全局基础函数。同样，我们假设在deal.II中使用的顶点的标准编号。
 *
 * @verbatim
 * 2---------|---------3
 * |                   |
 * |                   |
 * |                   |
 * |                   |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -
 * |                   |
 * |                   |
 * |                   |
 * |                   |
 * 0---------|---------1
 * @endverbatim
 *
 * 对于给定单元的每个顶点 $v_j$ ，有两条边，其中 $v_j$
 * 是其中的一个端点。考虑一个线性函数，使其在两条相邻边的中点上的值为0.5，而在其他边的两个中点上的值为0.0。请注意，这些值的集合满足上文所述的骰子规则。我们用
 * $\phi_j$ 表示与顶点 $v_j$
 * 相关的这样一个函数。那么4个形状函数的集合就是一个单元格上的统一分割。
 * $\sum_{j=0}^{3} \phi_j = 1$
 * 。（这很容易看出：在每个边的中点，四个函数的总和加起来是1，因为两个函数的值是0.5，另一个是0.0。因为该函数是全局线性的，唯一能在四个点上有值1的函数也必须是全局等于1的）。)
 * 以下数字表示 $\phi_j$ 对 $j=0,\cdots,3$ 的中点值。
 * <ul>   <li>  形状函数  $\phi_0$  。
 * @verbatim
 * +--------0.0--------+
 * |                   |
 * |                   |
 * |                   |
 * |                   |
 * 0.5                 0.0
 * |                   |
 * |                   |
 * |                   |
 * |                   |
 * +--------0.5--------+
 * @endverbatim
 *
 * <li>  形状函数  $\phi_1$  。
 * @verbatim
 * +--------0.0--------+
 * |                   |
 * |                   |
 * |                   |
 * |                   |
 * 0.0                 0.5
 * |                   |
 * |                   |
 * |                   |
 * |                   |
 * +--------0.5--------+
 * @endverbatim
 *
 * <li>  形状函数  $\phi_2$  。
 * @verbatim
 * +--------0.5--------+
 * |                   |
 * |                   |
 * |                   |
 * |                   |
 * 0.5                 0.0
 * |                   |
 * |                   |
 * |                   |
 * |                   |
 * +--------0.0--------+
 * @endverbatim
 *
 * <li>  形状函数  $\phi_3$  。
 * @verbatim
 * +--------0.5--------+
 * |                   |
 * |                   |
 * |                   |
 * |                   |
 * 0.0                 0.5
 * |                   |
 * |                   |
 * |                   |
 * |                   |
 * +--------0.0--------+
 * @endverbatim
 *
 * </ul>
 * 本地DoFs分别由与顶点相关的形状函数的系数定义。尽管这4个局部DoFs在单个单元内不是线性独立的，但这个定义对于全局DoFs的定义来说是一个很好的出发点。
 * 我们想强调的是，形状函数是在每个单元上构建的，而不是只在参考单元上。通常的有限元是基于
 * "参数化
 * "的概念来定义的。这意味着有限元的函数空间是在一个参考单元上定义的，它通过参考单元的映射被转换到每个单元。然而，P1不合格元素并不遵循这种概念。它在每个单元上定义了一个带有线性形状函数的函数空间，而没有借助于参考单元上的函数空间。换句话说，该元素是在现实空间中定义的，而不是通过参考单元的映射。在这一点上，它类似于FE_DGPN非参数化元素。
 * 因此，这个实现不需要在参考单元上计算形状值。相反，形状值是通过在每个单元上独立构建形状函数来计算的。
 * <h3>Degrees of freedom</h3>
 * 我们接下来要考虑元素的<i>global</i>基函数，因为我们最终要解决的方程组是一个全局系统，而不是局部。与一个节点相关的全局基函数是由每个元素上与该节点相关的局部形状函数的单元逐一组成来定义的。
 * 有一个关于全局基础函数的线性独立性的理论结果，它取决于我们考虑的边界条件的类型。
 * 当给出同质Dirichlet边界条件时，与内部节点相关的全局基础函数是线性独立的。那么，DoFs的数量等于内部节点的数量，因此与标准双线性
 * $Q_1$ 有限元的DoFs的数量相同。
 * 当给出诺伊曼边界条件时，与所有节点（包括边界节点）相关的全局基函数实际上不是线性独立的。存在着一个冗余。因此在这种情况下，DoFs的数量等于所有节点的数量减去1。这还是和普通
 * $Q_1$ 元素一样。 <h3>Unit support points</h3>
 * 对于平滑函数，我们通过使用其节点值作为DoF值，构造一个属于元素空间的分片线性函数。
 * 请注意，对于P1非线性元素，平滑函数的两个节点值和它的内插值一般不重合，这与普通的拉格朗日有限元素不同。当然，提到
 * "节点值
 * "是没有意义的，因为元素空间有不符合性。但是，即使与一个节点相关的单一全局基函数被认为是该节点上唯一的'节点值'，也是如此。例如，考虑与一个节点相关的基函数。考虑通过连接两个中点，分别代表值为0.5和0的水平集的两条线。然后我们通过沿着这两条线的对角线将四边形切成两个子三角形。它给出了另一个值为0.25的水平集，与切割对角线重合。因此，这三个水平集在四边形中都是平行的，它在底层节点上给出的是0.75的值，而不是1的值。无论四边形是矩形、平行四边形还是其他形状，都是如此。
 * <h3>Reference</h3>
 * Park和Sheen关于P1不合格元素的原始论文可在https://doi.org/10.1137/S0036142902404923，见
 * @cite park2003p  。
 *
 *
 */
class FE_P1NC : public FiniteElement<2, 2>
{
public:
  /**
   * P1不顺应元的构造器。  它只适用于二维和二维=0。
   *
   */
  FE_P1NC();

  virtual std::string
  get_name() const override;

  virtual UpdateFlags
  requires_update_flags(const UpdateFlags flags) const override;

  virtual std::unique_ptr<FiniteElement<2, 2>>
  clone() const override;

  /**
   * 解构器。
   *
   */
  virtual ~FE_P1NC() override = default;



private:
  /**
   * 返回由每个物体的自由度数量组成的向量。
   *
   */
  static std::vector<unsigned int>
  get_dpo_vector();

  /**
   * 返回给定单元上的4个局部线性形状函数 $\phi_j(x,y) = a x +
   * b y + c$
   * 的系数。对于每个局部形状函数，数组由三个系数组成，依次为a、b和c。
   *
   */
  static ndarray<double, 4, 3>
  get_linear_shape_coefficients(const Triangulation<2, 2>::cell_iterator &cell);

  /**
   * 在进行单元格数据计算之前，要做一些必要的工作。
   * 由于形状函数是在每个单元格上独立构建的，所以参考单元格上的数据是不需要的。
   * 它返回一个@ InternalDataBase的空变量类型，并更新@
   * update_flags，如果需要的话，为每个单元计算琐碎的零Hessian。
   *
   */
  virtual std::unique_ptr<FiniteElement<2, 2>::InternalDataBase>
  get_data(
    const UpdateFlags update_flags,
    const Mapping<2, 2> &,
    const Quadrature<2> &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<2, 2>
      &output_data) const override;

  using FiniteElement<2, 2>::get_face_data;

  virtual std::unique_ptr<FiniteElement<2, 2>::InternalDataBase>
  get_face_data(
    const UpdateFlags update_flags,
    const Mapping<2, 2> &,
    const hp::QCollection<1> &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<2, 2>
      &output_data) const override;

  virtual std::unique_ptr<FiniteElement<2, 2>::InternalDataBase>
  get_subface_data(
    const UpdateFlags update_flags,
    const Mapping<2, 2> &,
    const Quadrature<1> &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<2, 2>
      &output_data) const override;

  /**
   * 计算当前单元格的数据。
   *
   */
  virtual void
  fill_fe_values(
    const Triangulation<2, 2>::cell_iterator &cell,
    const CellSimilarity::Similarity          cell_similarity,
    const Quadrature<2> &                     quadrature,
    const Mapping<2, 2> &                     mapping,
    const Mapping<2, 2>::InternalDataBase &   mapping_internal,
    const internal::FEValuesImplementation::MappingRelatedData<2, 2>
      &                                          mapping_data,
    const FiniteElement<2, 2>::InternalDataBase &fe_internal,
    internal::FEValuesImplementation::FiniteElementRelatedData<2, 2>
      &output_data) const override;

  using FiniteElement<2, 2>::fill_fe_face_values;

  /**
   * 计算当前单元格面上的数据。
   *
   */
  virtual void
  fill_fe_face_values(
    const Triangulation<2, 2>::cell_iterator &cell,
    const unsigned int                        face_no,
    const hp::QCollection<1> &                quadrature,
    const Mapping<2, 2> &                     mapping,
    const Mapping<2, 2>::InternalDataBase &   mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<2, 2>
      &                     mapping_data,
    const InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<2, 2>
      &output_data) const override;

  /**
   * 计算当前单元格的子面的数据。
   *
   */
  virtual void
  fill_fe_subface_values(
    const Triangulation<2, 2>::cell_iterator &cell,
    const unsigned int                        face_no,
    const unsigned int                        sub_no,
    const Quadrature<1> &                     quadrature,
    const Mapping<2, 2> &                     mapping,
    const Mapping<2, 2>::InternalDataBase &   mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<2, 2>
      &                     mapping_data,
    const InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<2, 2>
      &output_data) const override;

  /**
   * 创建挂边的约束矩阵。
   *
   */
  void
  initialize_constraints();
};



 /*@}*/ 


DEAL_II_NAMESPACE_CLOSE

#endif


