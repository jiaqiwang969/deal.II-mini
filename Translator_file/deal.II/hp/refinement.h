//include/deal.II-translator/hp/refinement_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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

#ifndef dealii_hp_refinement_h
#define dealii_hp_refinement_h


#include <deal.II/base/config.h>

#include <deal.II/base/template_constraints.h>

#include <functional>
#include <vector>


DEAL_II_NAMESPACE_OPEN


// forward declarations
#ifndef DOXYGEN
template <typename Number>
class Vector;

template <int dim, int spacedim>
class DoFHandler;
#endif


namespace hp
{
  /**
   * 我们提供自适应方法，使计算资源与数值解的复杂性相一致。误差估计是确定哪里需要调整的适当手段。    然而通过hp-adaptivity，我们有两种方法来实现这些调整。  对于不规则解，动态分配单元尺寸的h-adaptive方法倾向于减少近似误差，而对于平滑解，p-adaptive方法则更适合于动态选择函数空间。这个命名空间收集了决定应用哪种类型的自适应方法的工具。    <h3>Usage</h3> 为了成功应用hp-adaptive方法，我们建议采用以下工作流程。    <ol>   <li>  一个合适的误差估计是任何一种自适应方法的基础。  与纯网格细化类似，我们将以通常的方式确定误差估计值（即KellyErrorEstimator），并对单元格进行细化或粗化（即GridRefinement）。    在这个阶段调用 Triangulation::execute_coarsening_and_refinement() 将按预期执行纯网格细化。      <li>  一旦所有的细化和粗化标志分布在网格上，我们可以确定这些是否符合p-adaptive方法的要求。  相应的函数将在细化和粗化标志之上设置 @p future_fe_indices ，如果它们满足一定的标准。    在细化的情况下，基础 hp::FECollection 的上位元素将被指定为未来的有限元素。  相应地，下级元素将被选择为粗化。      Triangulation::execute_coarsening_and_refinement()  现在将独立提供h-和p-自适应方法。      <li>  现在，可能会有单元同时安排h-适应和p-适应。  如果我们不想同时施加这两种方法，我们需要决定为每个单元单独和明确地选择哪一种。由于网格细化将被默认施加，我们只在上面确定p适应性的资格，所以我们将总是决定支持p适应性的方法。    现在，调用 Triangulation::execute_coarsening_and_refinement() 将在每个单元上唯一地执行h-或p-adaptive方法。      <li>  到此为止，每个单元都知道自己在适应性方面的命运。  我们现在可以继续准备所有的数据结构，以便在网格变化时进行传输。之前设置的细化和粗化标志以及 @p future_fe_indices 将被用于相应地更新数据。    </ol>  作为一个例子，纯p-adaptive方法的实现将看起来像这样。
   * @code
   * // step 1: flag cells for refinement or coarsening
   * Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
   * KellyErrorEstimator<dim>::estimate(
   *   hp_dof_handler,
   *   QGauss<dim-1> (quadrature_points),
   *   std::map<types::boundary_id, const Function<dim, Number>>(),
   *   solution,
   *   estimated_error_per_cell);
   * GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
   *                                                 estimated_error_per_cell,
   *                                                 top_fraction,
   *                                                 bottom_fraction);
   *
   * // step 2: set future finite element indices on flagged cells
   * hp::Refinement::full_p_adaptivity(hp_dof_handler);
   *
   * // step 3: decide whether h- or p-adaptive methods will be supplied
   * hp::Refinement::force_p_over_h(hp_dof_handler);
   *
   * // step 4: prepare solutions to be transferred
   * ...
   *
   * triangulation.execute_coarsening_and_refinement();
   * @endcode
   * @ingroup hp
   *
   */
  namespace Refinement
  {
    /**
     * 一个别名，它定义了一个函数的特性，可以作为决定是进行h-适应还是p-适应的比较标准。
     * 这样的函数需要两个数字作为参数。第一个对应于所提供的标准，而另一个则符合参考。
     * 比较的结果将以布尔值的形式返回。
     *
     */
    template <typename Number>
    using ComparisonFunction =
      std::function<bool(const Number &, const Number &)>;

    /**
     * @name  设置p-adaptivity标志  
     * @{ 
     *
     */

    /**
     * 每个被标记为h-精简的单元也将被标记为p-精简。
     * 这同样适用于粗化。
     * @note   Triangulation::prepare_coarsening_and_refinement()  和
     * hp::Refinement::limit_p_level_difference()
     * 可能改变细化和粗化标志以及未来的有限元指数。
     * 避免在这个特殊函数之前调用它们。
     *
     */
    template <int dim, int spacedim>
    void
    full_p_adaptivity(const dealii::DoFHandler<dim, spacedim> &dof_handler);

    /**
     * 通过参数 @p p_flags.
     * 对那些被特别标记为p-adaptation的单元进行调整，只有当单元事先被标记为细化和粗化时，才会分配未来的有限元素。
     * 参数 @p p_flags
     * 的每个条目都需要对应于一个活动单元。
     * @note   Triangulation::prepare_coarsening_and_refinement()  和
     * hp::Refinement::limit_p_level_difference()
     * 可能改变细化和粗化标志以及未来的有限元指数。
     * 避免在这个特定函数之前调用它们。
     *
     */
    template <int dim, int spacedim>
    void
    p_adaptivity_from_flags(
      const dealii::DoFHandler<dim, spacedim> &dof_handler,
      const std::vector<bool> &                p_flags);

    /**
     * 适应在标准满足一定绝对阈值的单元上使用哪种有限元。
     * 对于p-精化和p-粗化，需要通过参数 @p p_refine_threshold
     * 和 @p p_coarsen_threshold. 提供两个单独的阈值。
     * 如果一个单元当前被标记为精化或粗化，并且其标准成功地与相应的阈值进行比较，我们就认为它是p-adaptivity。让我们对默认情况进行更具体的说明。如果一个单元被标记为细化，并且其准则大于或等于相应的阈值，我们就认为它是p-细化的。这同样适用于p-粗化，但单元格的准则必须低于或等于阈值。然而，可以通过参数
     * @p compare_refine 和 @p compare_coarsen
     * 提供不同的比较函数对象，以施加不同的决策策略。
     * 参数 @p criteria
     * 的每个条目都需要对应于一个活动单元。
     * @note   Triangulation::prepare_coarsening_and_refinement()  和
     * hp::Refinement::limit_p_level_difference()
     * 可能改变细化和粗化标志以及未来的有限元指数。
     * 避免在这个特定函数之前调用它们。
     *
     */
    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_absolute_threshold(
      const dealii::DoFHandler<dim, spacedim> &dof_handler,
      const Vector<Number> &                   criteria,
      const Number                             p_refine_threshold,
      const Number                             p_coarsen_threshold,
      const ComparisonFunction<typename identity<Number>::type>
        &compare_refine = std::greater_equal<Number>(),
      const ComparisonFunction<typename identity<Number>::type>
        &compare_coarsen = std::less_equal<Number>());

    /**
     * 适应在单元格上使用哪种有限元，这些单元格的标准相对于标准值的整体范围满足一定的阈值。
     * 阈值将根据当前设置的细化标记为细化和粗化单元分别确定。对于每一类单元格，我们确定所有标准的最大和最小值，并通过这些极限之间的线性内插来确定阈值。
     * 参数 @p p_refine_fraction 和 @p p_refine_coarsen
     * 被用作插值因子，其中`0`对应于最小值，`1`对应于最大值。默认情况下，均值被认为是阈值。
     * 如果一个单元目前被标记为细化或粗化，并且其标准成功地与相应的阈值相比较，我们就认为它是p-adaptivity。让我们对默认情况进行更具体的说明。如果一个单元被标记为细化，并且其准则大于或等于相应的阈值，我们就认为它是p-细化的。这同样适用于p-粗化，但单元格的准则必须低于或等于阈值。然而，可以通过参数
     * @p compare_refine 和 @p compare_coarsen
     * 提供不同的比较函数对象，以施加不同的决策策略。
     * 参数 @p criteria
     * 的每个条目都需要对应于一个活动单元。参数 @p
     * p_refine_fraction 和 @p p_coarsen_fraction 需要在区间 $[0,1]$  .
     * @note   Triangulation::prepare_coarsening_and_refinement()  和
     * hp::Refinement::limit_p_level_difference()
     * 可能会改变细化和粗化标志，以及未来的有限元索引。
     * 避免在这个特定函数之前调用它们。
     *
     */
    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_relative_threshold(
      const dealii::DoFHandler<dim, spacedim> &dof_handler,
      const Vector<Number> &                   criteria,
      const double                             p_refine_fraction  = 0.5,
      const double                             p_coarsen_fraction = 0.5,
      const ComparisonFunction<typename identity<Number>::type>
        &compare_refine = std::greater_equal<Number>(),
      const ComparisonFunction<typename identity<Number>::type>
        &compare_coarsen = std::less_equal<Number>());

    /**
     * 适应在给定的一部分单元上使用哪一个有限元。
     * 在所有被标记为某种适应类型的单元中，无论是细化还是粗化，我们将在这个子集中确定固定数量的单元，这些单元将被标记为相应的p-适应变量。
     * 对于每个细化和粗化子集，我们将根据提供的参数 @p
     * criteria
     * 确定一个阈值，其中包含每个活动单元的指标。在细化的默认情况下，所有指标大于或等于相应阈值的单元将被考虑为p-细化，而对于粗化，所有指标小于或等于匹配阈值的单元将被考虑。然而，可以通过参数
     * @p compare_refine 和 @p compare_coarsen
     * 提供不同的比较函数对象来实施不同的决策策略。
     * 对于细化，阈值将与具有 @p p_refine_fraction 乘以
     * Triangulation::n_active_cells()
     * 最大指标的单元相关联，而对于粗化，它是具有 @p
     * p_refine_coarsen 乘以 Triangulation::n_active_cells()
     * 最低指标的单元。        参数 @p criteria
     * 的每个条目需要对应于一个活动单元。参数 @p
     * p_refine_fraction 和 @p p_coarsen_fraction 需要在区间 $[0,1]$  .
     * @note   Triangulation::prepare_coarsening_and_refinement()  和
     * hp::Refinement::limit_p_level_difference()
     * 可能会改变细化和粗化标志，以及未来的有限元索引。
     * 避免在这个特定函数之前调用它们。
     *
     */
    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_fixed_number(
      const dealii::DoFHandler<dim, spacedim> &dof_handler,
      const Vector<Number> &                   criteria,
      const double                             p_refine_fraction  = 0.5,
      const double                             p_coarsen_fraction = 0.5,
      const ComparisonFunction<typename identity<Number>::type>
        &compare_refine = std::greater_equal<Number>(),
      const ComparisonFunction<typename identity<Number>::type>
        &compare_coarsen = std::less_equal<Number>());

    /**
     * 根据（未知的）分析解的规则性来调整单元格上使用的有限元。
     * 通过局部Sobolev正则性指数的近似值 $k_K$
     * ，我们可以评估单元格上的局部解属于哪个有限元空间
     * $K$
     * 。由于正则性指数只是一个估计值，我们不会用它来直接分配有限元空间，而是将其视为一个适应性指标。如果一个单元被标记为细化，一旦它满足
     * $k_K > p_{K,\text{super}}$ ，我们将执行p-细化，其中
     * $p_{K,\text{super}}$
     * 是单元上当前活动元素的有限元上级的多项式程度  $K$
     * 。在粗化的情况下，必须满足 $k_K < p_{K,\text{sub}}$
     * 的标准， $p_{K,\text{sub}}$ 是下级元素的度数。
     * 参数 @p sobolev_indices
     * 的每个条目需要对应于一个活动单元。
     * 更多理论细节见  @cite ainsworth1998hp  。
     * @note   Triangulation::prepare_coarsening_and_refinement()  和
     * hp::Refinement::limit_p_level_difference()
     * 可能会改变细化和粗化标志，以及未来的有限元索引。
     * 避免在这个特定函数之前调用它们。
     *
     */
    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_regularity(
      const dealii::DoFHandler<dim, spacedim> &dof_handler,
      const Vector<Number> &                   sobolev_indices);

    /**
     * 根据每个单元的标准与参考值的关系，调整在每个单元上使用的有限元。
     * 如果一个单元当前被标记为细化或粗化，并且其准则成功地与相应的参考比较，我们就考虑对其进行p-adaptivity。除了函数p_adaptivity_from_absolute_threshold()和p_adaptivity_from_relative_threshold()，比较函数对象必须通过参数
     * @p compare_refine 和 @p compare_coarsen. 明确提供，参数 @p
     * criteria 和 @p references
     * 的每个条目需要对应于一个活动单元。
     * @note   Triangulation::prepare_coarsening_and_refinement()  和
     * hp::Refinement::limit_p_level_difference()
     * 可能会改变细化和粗化标志，以及未来的有限元索引。
     * 避免在这个特定函数之前调用它们。
     *
     */
    template <int dim, typename Number, int spacedim>
    void
    p_adaptivity_from_reference(
      const dealii::DoFHandler<dim, spacedim> &                  dof_handler,
      const Vector<Number> &                                     criteria,
      const Vector<Number> &                                     references,
      const ComparisonFunction<typename identity<Number>::type> &compare_refine,
      const ComparisonFunction<typename identity<Number>::type>
        &compare_coarsen);

    /**
     * @}
     *
     */

    /**
     * @name  误差预测  
     * @{ 
     *
     */

    /**
     * 预测当前 @p error_indicators 在所提供的 @p dof_handler,
     * 上发生细化和粗化后将如何适应，并将其结果写入 @p
     * predicted_errors.   @p error_indicators 和 @p predicted_errors
     * 的每个条目对应于底层三角图上的一个活动单元，因此每个容器的大小必须是
     * Triangulation::n_active_cells().
     * 误差被解释为以能量准则测量；这一假设进入了预测中使用的收敛率。输出参数的
     * $l_2$ -norm @p predicted_errors
     * 对应于适应后预测的全局误差。
     * 对于p-适应，预计局部误差会随着分配的有限元的多项式程度呈指数级收敛。因此，程度的每一次增加或减少将通过用户定义的控制参数
     * @p gamma_p. 改变其值。 对于h-适应，我们期望单元 $K$
     * 上的局部误差 $\eta_K$ 与能量准则中的 $(h_K)^{p_K}$
     * 成正比，其中 $h_K$ 表示单元直径， $p_K$
     * 表示当前在单元 $K$ 上分配有限元素的多项式程度。
     * 在h粗化过程中，兄弟姐妹上的有限元可能是不同的，他们的父单元将被分配给属于其最一般的子单元的最小支配有限元。因此，我们将总是在一个封闭的有限元空间上进行插值。
     * 此外，假设要粗化的单元上的有限元足以正确表示解决方案（例如，对于二次解来说，至少是二次基函数），我们有信心说，仅在较大的有限元空间上插值，误差不会改变。
     * 对于p-adaptation，预计局部误差会随着分配的有限元的多项式程度呈指数级收敛。因此，度数的每一次增加或减少将通过用户定义的控制参数改变其值
     * @p gamma_p.
     * 指数收敛的假设只有在h-适应和p-适应方法结合在一起的情况下才有效，即它们都在整个网格中被利用，但不必在一个单元上同时应用。
     * 预测算法制定如下，控制参数 @p gamma_p,  @p gamma_h 和 @p
     * gamma_n
     * 可用于单独影响每种适应类型的预测结果。每个单独单元的结果都储存在
     * @p predicted_errors 输出参数中。      <table> <tr><th>Adaptation
     * type <th colspan="2">Prediction formula <tr><td>no adaptation
     * <td>$\eta_{K,\text{pred}} = \eta_{K} \, \gamma_\text{n}$
     * <td>$\gamma_\text{n} \in (0,\infty)$ <tr><td>p-adaptation
     * <td>$\eta_{K,\text{pred}} = \eta_{K} \,
     * \gamma_\text{p}^{(p_{K,\text{future}}
     *
     * - p_K)}$ <td>$\gamma_\text{p} \in (0,1)$ <tr><td>hp-refinement
     * <td>$\eta_{K,\text{pred}} = \eta_{K} \, \gamma_\text{h} \,
     * 0.5^{p_{K,\text{future}}} \, \gamma_\text{p}^{(p_{K,\text{future}}
     *
     * - p_{K})}$ <td rowspan="2">$\gamma_\text{h} \in (0,\infty)$
     * <tr><td>hp-coarsening <td>$\eta_{K,\text{pred}} = \eta_{K} \,
     * (\gamma_\text{h} \, 0.5^{p_{K,\text{future}}})^{-1} \,
     * \gamma_\text{p}^{(p_{K,\text{future}}
     *
     * - p_{K})}$ </table>
     * 在细化历史的基础上，我们使用预测的误差估计值来决定在下一个适应步骤中如何适应单元。
     * 将上一适应步骤的预测误差与当前步骤的误差估计值相比较，可以证明我们之前选择的适应是否合理，并让我们决定在下一步骤中如何适应。
     * 因此，我们必须将预测的误差从旧的网格转移到适应的网格。当把预测误差转移到适应的网格时，确保配置你的CellDataTransfer对象，将
     * AdaptationStrategies::Refinement::l2_norm() 作为细化策略，
     * AdaptationStrategies::Coarsening::l2_norm() 作为粗化策略。
     * 这可以确保在两个网格上都保留预测误差的 $l_2$
     * -norm。
     * 在这种情况下，我们假设被h精化的单元上的局部误差将在其所有的
     * $n_{K_c}$
     * 子单元上平均分配，而在h粗化的情况下，兄弟姐妹的局部误差将在父单元上加总。这个假设在实践中经常不被满足。例如，如果一个单元处于角部奇点，那么最终最接近奇点的一个子单元将继承大部分的剩余误差
     *
     * - 但这个函数不可能知道奇点在哪里，因此假定是平均分配。        结合从旧网格到适应网格的转移，完整的误差预测算法如下。      <table>
     * <tr><th>Adaptation type <th colspan="2">Prediction formula <tr><td>no
     * adaptation <td>$\eta_{K,\text{pred}} = \eta_{K} \, \gamma_\text{n}$
     * <td>$\gamma_\text{n} \in (0,\infty)$ <tr><td>p-adaptation
     * <td>$\eta_{K,\text{pred}} = \eta_{K} \,
     * \gamma_\text{p}^{(p_{K,\text{future}}
     *
     * - p_K)}$ <td>$\gamma_\text{p} \in (0,1)$ <tr><td>hp-refinement
     * <td>$\left( \eta_{K_c,\text{pred}} \right)^2 = n_{K_c}^{-1} \left(
     * \eta_{K_p} \, \gamma_\text{h} \, 0.5^{p_{K_c,\text{future}}} \,
     * \gamma_\text{p}^{(p_{K_c,\text{future}}
     *
     * - p_{K_p})} \right)^2 \quad \forall K_c \text{ children of } K_p$ <td
     * rowspan="2">$\gamma_\text{h} \in (0,\infty)$ <tr><td>hp-coarsening
     * <td>$\left( \eta_{K_p,\text{pred}} \right)^2 = \sum\limits_{K_c} \left(
     * \eta_{K_c} \, (\gamma_\text{h} \, 0.5^{p_{K_p,\text{future}}})^{-1} \,
     * \gamma_\text{p}^{(p_{K_p,\text{future}}
     *
     * - p_{K_c})} \right)^2 \quad \forall K_c \text{ children of } K_p$
     * </table>
     * 有了这些预测的误差估计值，我们就能够根据单元的细化历史或者说是预测的误差估计值的变化来调整单元的有限元。
     * 如果一个单元被标记为适应，一旦单元 $K$
     * 上的相关误差指标 $\eta_{K}$ 满足 $\eta_{K} <
     * \eta_{K,\text{pred}}$ ，我们要进行p-适应，其中下标
     * $\text{pred}$
     * 表示预测的误差。这对应于我们对平滑性的假设是正确的，否则就采用h-适应性。我们用函数
     * hp::Refinement::p_adaptivity_from_reference() 和一个函数对象
     * `std::less<Number>()` 来实现这一点，这两个比较器参数。
     * 另外，通过另一种策略，我们可以确定所有要适应的细胞中h-适应和p-适应的比例。为此，使用
     * hp::Refinement::p_adaptivity_fixed_number() 和标准
     * $(\eta_{K,\text{pred}}
     *
     * - \eta_{K})$  。
     * 在这两种情况下的第一个适应步骤中，用户需要决定是要进行h-适应还是p-适应。使用
     * $\eta_{K,\text{pred}} = 0$ 时将采用h步，而 $\eta_{K,\text{pred}}
     * = \infty$ 则确保p步。后者可以通过
     * `std::numeric_limits::infinity()`. 实现。
     * 下面的代码片段演示了如何在应用中根据细化历史施加hp-adaptivity。
     * @code
     * // [initialisation...]
     * Vector<float> predicted_error_per_cell(triangulation.n_active_cells());
     * for(unsigned int i = 0; i < triangulation.n_active_cells(); ++i)
     * predicted_error_per_cell[i] = std::numeric_limits<float>::infinity();
     *
     * // [during each refinement step...]
     * // set h-adaptivity flags
     * Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
     * KellyErrorEstimator::estimate(...);
     * GridRefinemet::refine_and_coarsen_fixed_{number|fraction}(...);
     *
     * // set p-adaptivity flags
     * hp::Refinement::p_adaptivity_from_reference(
     * hp_dof_handler,
     * estimated_error_per_cell,
     * predicted_error_per_cell,
     * std::less<float>(),
     * std::less<float>());
     * hp::Refinement::{choose|force}_p_over_h(hp_dof_handler);
     *
     * // predict error for the subsequent adaptation
     * triangulation.prepare_coarsening_and_refinement();
     * hp::Refinement::predict_error(
     * hp_dof_handler,
     * estimated_error_per_cell,
     * predicted_error_per_cell);
     *
     * // perform adaptation
     * CellDataTransfer<dim, spacedim, Vector<float>> cell_data_transfer(
     * triangulation,
     * false,
     * &AdaptationStrategies::Refinement::l2_norm<dim, spacedim, float>,
     * &AdaptationStrategies::Coarsening::l2_norm<dim, spacedim, float>);
     * cell_data_transfer.prepare_coarsening_and_refinement();
     *
     * triangulation.execute_coarsening_and_refinement();
     *
     * Vector<float> transferred_errors(triangulation.n_active_cells());
     * cell_data_transfer.unpack(predicted_error_per_cell, transferred_errors);
     * predicted_error_per_cell = std::move(transferred_errors);
     * @endcode
     * 更多的理论细节见  @cite melenk2001hp
     * ，这个函数的默认参数也来自这里，即  $\gamma_\text{p}^2
     * = 0.4$  、  $\gamma_\text{h}^2 = 4$  、  $\gamma_\text{n}^2 = 1$  。
     * 如果你正在处理  parallel::distributed::Triangulation
     * 对象，你需要特别注意。在这里，p4est决定了网格细化的细节，因此，当我们在适应过程中确定预测的误差时，它能产生更可靠和可信的结果。我们可以通过将这个函数附加到信号
     * Triangulation::Signals::post_p4est_refinement,
     * 中来实现，该信号在p4est得到细化后，但在数据准备传输前被触发。Triangulation对象的细化和粗化标志需要使用
     * internal::parallel::distributed::TemporarilyMatchRefineFlags.
     * 与已经细化的p4est神谕相匹配。因此，像下面这样的结构对于正确预测并行分布式应用中的错误是必要的。
     * @code
     * Vector<float> predicted_errors;
     * triangulation.signals.post_p4est_refinement.connect([&]() {
     * const parallel::distributed::TemporarilyMatchRefineFlags<dim>
     *   refine_modifier(triangulation);
     * predicted_errors.reinit(triangulation.n_active_cells());
     * hp::Refinement::predict_error(dof_handler,
     *                               error_indicators,
     *                               predicted_errors);
     * });
     * @endcode
     * 容器 <code>predicted_errors</code> 然后需要遵循通常的
     * parallel::distributed::CellDataTransfer 工作流程。
     * @note
     * 我们想通过适应性的方式来预测错误的实际发生。
     * 因此，这个函数需要在
     * Triangulation::prepare_coarsening_and_refinement() 和
     * hp::Refinement::limit_p_level_difference(). 之后调用。
     *
     */
    template <int dim, typename Number, int spacedim>
    void
    predict_error(const dealii::DoFHandler<dim, spacedim> &dof_handler,
                  const Vector<Number> &                   error_indicators,
                  Vector<Number> &                         predicted_errors,
                  const double gamma_p = std::sqrt(0.4),
                  const double gamma_h = 2.,
                  const double gamma_n = 1.);

    /**
     * @}
     *
     */

    /**
     * @name  在h-adaptivity和p-adaptivity之间做出决定  
     * @{ 
     *
     */

    /**
     * 在任何情况下都选择p-adaptivity而不是h-adaptivity。
     * 移除所有分配有 @p future_fe_index
     * 的单元格的细化和粗化标志。
     * @note   Triangulation::prepare_coarsening_and_refinement()  和
     * hp::Refinement::limit_p_level_difference()
     * 可能改变细化和粗化标志以及未来的有限元指数。
     * 避免在这个特定函数之前调用它们。
     *
     */
    template <int dim, int spacedim>
    void
    force_p_over_h(const dealii::DoFHandler<dim, spacedim> &dof_handler);

    /**
     * 只要在所有相关单元上调用，就选择p-adaptivity而不是h-adaptivity。        在细化的情况下，关于有限元的信息将被继承。因此，只要有需要，我们就会选择p-细化而不是h-细化，即清除细化标志并提供相应的 @p future_fe_index. 然而对于粗化，我们遵循不同的方法。将一个单元标记为h-粗化并不意味着它最终会被粗化。只有当一个单元和它的所有兄弟姐妹被标记时，它们才会被合并到它们的父单元中。如果我们在上面考虑p粗化，我们必须对所有的兄弟姐妹一起决定他们将如何被粗化。我们区分了三种不同的情况。      <ol>   <li>  不是所有兄弟姐妹都被标记为粗化：p-粗化。        <br>  我们保留  @p future_fe_indices  并清除所有兄弟姐妹的粗化标志。      <li>  所有的兄弟姐妹都被标记为粗化，但不是所有的p-adaptation：h-coarsening。        <br>  我们保留粗化标志，并清除所有兄弟姐妹上的所有 @p future_fe_indices 。      <li>  所有的兄弟姐妹都被标记为粗化和p-适应：p-粗化。        <br>  我们保留 @p future_fe_indices 并清除所有兄弟姐妹上的粗化标志。      </ol>
     * @note  函数 Triangulation::prepare_coarsening_and_refinement()
     * 将清理所有h-粗化标志，如果它们不在所有兄弟姐妹中共享的话。在hp情况下，我们需要把这个决定提前。
     * 如果单元格不会被粗化，但符合p-adaptivity的条件，我们必须相应地设置所有标志。所以这个函数预计了
     * Triangulation::prepare_coarsening_and_refinement()
     * 以后会做出的决定。
     * @note   Triangulation::prepare_coarsening_and_refinement()  和
     * hp::Refinement::limit_p_level_difference()
     * 可能会改变细化和粗化标志，以及未来的有限元指数。
     * 避免在这个特定函数之前调用它们。
     *
     */
    template <int dim, int spacedim>
    void
    choose_p_over_h(const dealii::DoFHandler<dim, spacedim> &dof_handler);

    /**
     * @}
     *
     */

    /**
     * @name  优化p级分布  
     * @{ 
     *
     */

    /**
     * 限制相邻单元之间的p级差。
     * 本质上，对未来FE指数的作用与
     * Triangulation::prepare_coarsening_and_refinement()
     * 对细化标志的作用相同。
     * 详细来说，这个函数限制了相邻单元的水平差异，从而平滑了整个函数空间。未来的FE指数将被提高（而不是降低），因此与相邻单元的水平差永远不会大于
     * @p max_difference.  多个FE层次可能已经通过
     * hp::FECollection::set_hierarchy().
     * 注册，该函数只对一个层次操作，即包含FE指数的那个层次
     * @p contains_fe_index.
     * 未来FE指数不属于相应层次的单元将被忽略掉。
     * 该函数可以在进行适应性调整之前选择性地调用
     * Triangulation::execute_coarsening_and_refinement().
     * 没有必要调用该函数，在库的任何部分也不会自动调用该函数（与三角法的对应函数相反）。
     * 在将被h-coarsened的单元格上，我们强制执行差异标准，就像它已经是一个父单元格一样。也就是说，我们将所有兄弟姐妹的级别设置为其中最高的一个。在这种情况下，所有的同级单元都需要事先通过
     * Triangulation::prepare_coarsening_and_refinement()
     * 终端设置h-coarsenening标志。否则将触发一个断言。
     * 返回是否有任何未来的FE指数被此函数改变。
     *
     */
    template <int dim, int spacedim>
    bool
    limit_p_level_difference(
      const dealii::DoFHandler<dim, spacedim> &dof_handler,
      const unsigned int                       max_difference    = 1,
      const unsigned int                       contains_fe_index = 0);

    /**
     * @}
     *
     */
  } // namespace Refinement
} // namespace hp


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_hp_refinement_h


