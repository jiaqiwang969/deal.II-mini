//include/deal.II-translator/fe/fe_q_iso_q1_0.txt
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

#ifndef dealii_fe_q_iso_q1_h
#define dealii_fe_q_iso_q1_h

#include <deal.II/base/config.h>

#include <deal.II/base/polynomials_piecewise.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_q_base.h>

DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * 标量拉格朗日有限元 @p Qp-iso-Q1
 * 的实现，它定义了在每个坐标方向上具有 @p p
 * 细分的连续、分片线性元素的有限元空间。它产生的元素具有与
 * @p Qp
 * 元素相同的自由度数，但使用线性插值而不是高阶元素。这种类型的元素在文献中也被称为宏观元素，因为它实际上由几个较小的元素组成，即<i>p</i><tt><sup>dim</sup></tt>这样的子单元。
 * 自由度的编号方式与 @p p.
 * 度的FE_Q完全相同，关于自由度在一个元素内如何编号的详细描述见那里。
 * 这个元素代表了一个缩小了的网格尺寸的Q-线性有限元空间<i>h/p</i>。如果使用等效的正交法，其效果相当于在更细的网格上使用1度的FE_Q，其系数为
 * @p p
 * 。然而，这个元素减少了选择（自适应）网格大小的灵活性，正好是这个系数
 * @p p,
 * ，通常会降低效率。另一方面，在同一网格上，将该元素与
 * @p p 细分的FE_Q度元素 @p p
 * 相比较，可以看出对于光滑问题，收敛性通常要差得多。特别是，
 * @p Qp
 * 元素在L2准则中达到了<i>h<sup>p+1</sup></i>的插值阶数，而这些元素只达到<i>(h/p)<sup>2</sup></i>。由于这两个原因，这个元素作为一个独立的元素通常不是很有用。此外，用这个元素对元素内部边界上的面项进行任何评估都是不可能的，因为deal.II对单元内部的低维积分没有FEFaceValues的等价物。
 * 尽管如此，在一些用例中，这个元素实际上是有用的。  <ol>
 * <li>
 * 在PDEs系统中，某些变量需要比其他变量更高的分辨率，额外的自由度应该用在提高线段的分辨率而不是高阶多项式上，而且你不想为不同的组件使用两个不同的网格。当解中出现不规则现象（冲击）时，就会出现这种情况，而稳定技术对线型元素有效，但对高阶元素无效。
 * </li>
 * <li>  Stokes/Navier Stokes系统，如 step-22
 * 中讨论的系统，可以用Q2-iso-Q1元素代替Q2元素来解决速度问题。与Q1压力相结合，它们给出了一个稳定的混合元对。然而，在大多数情况下，它们的表现比标准（Taylor-Hood
 * $Q_2\times Q_1$ ）方法差。   </li>
 * <li>  用基于 @p Qp-iso-Q1 元素的预处理器对高阶的FE_Q系统 @p p 进行预处理。一些预处理程序，如代数多重网格，使用线性元素比使用高阶元素的效果要好得多，因为它们通常隐含地假设条目之间有稀疏的连接性。然后，在这些元素的基础上创建一个预处理矩阵，产生相同数量的自由度（和一个频谱等效的线性系统），它可以在像CG这样的迭代求解器中与（高阶）系统矩阵相结合。   </li>   </ol> 。
 * <h3>Appropriate integration</h3>
 * 由于这些元素的性质是线段的串联，在为这个元素选择正交公式时必须注意。对于
 * @p p 子元素的标准选择是公式<tt>QIterated<dim>(QGauss<1>(2),
 * p)</tt>，它对应于在更细的网格上用于积分函数的公式。这与FE_Q(p)形成对比，其中QGauss<dim>(p+1)是默认选择。特别是，必须注意不要使用在子元素边界上评估基础函数（及其导数）的正交公式，因为内部边界上的片断函数梯度被设置为零。内部没有进行任何检查以确保这种情况不会发生
 *
 * - 避免这些情况是用户的责任。
 * 还需要注意的是，与FE_Q相比，常规的deal.II设置稀疏性模式和装配矩阵的程序并没有利用这个元素中增加的稀疏性。这是因为
 * DoFTools::make_sparsity_pattern
 * 假设元素内所有自由度之间存在耦合，而FE_Q_iso_Q1有多个细分，耦合度确实较低。
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_Q_iso_Q1 : public FE_Q_Base<dim, spacedim>
{
public:
  /**
   * 构建一个具有给定细分数的FE_Q_iso_Q1元素。细分的数量与FE_Q中的程度相似，即两个元素产生相同的自由度数量。
   *
   */
  FE_Q_iso_Q1(const unsigned int n_subdivisions);

  /**
   * 返回一个能唯一识别有限元的字符串。该类返回<tt>FE_Q_iso_q1<dim>(equivalent_degree)</tt>，其中
   * @p dim 和 @p equivalent_degree由适当的值代替。
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
   * @name 支持hp的函数  
     * @{ 
   *
   */

  /**
   * @copydoc   FiniteElement::compare_for_domination()  支持HP的函数
   *
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

  //@}
};



 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


