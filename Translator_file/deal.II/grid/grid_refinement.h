//include/deal.II-translator/grid/grid_refinement_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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

#ifndef dealii_grid_refinement_h
#define dealii_grid_refinement_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/numerics/vector_tools_common.h>

#include <limits>

DEAL_II_NAMESPACE_OPEN

// forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class Triangulation;
template <typename Number>
class Vector;
#endif

/**
 * 这个命名空间提供了一个帮助细化和粗化三角形的函数集合。尽管命名空间的名称，这些函数实际上并不<i>refine</i>三角化，而只是<i>mark
 * cells for refinement or
 * coarsening</i>。换句话说，它们执行自适应有限元循环中典型的
 * "求解-估计-标记-细化 "循环中的 "标记 "部分。
 * 这个命名空间中的函数形成两类。有辅助函数refine()和coarsen()。对用户来说更重要的是其他函数，它们实现了细化策略，在自适应有限元方法的文献中可以找到。关于这些方法的数学讨论，可以参考D&ouml;rfler,
 * Morin, Nochetto, Rannacher, Stevenson等人的作品。
 *
 *
 * @ingroup grid
 *
 */
namespace GridRefinement
{
  /**
   * 返回一对双倍值，其中第一个是调整后的细化单元分数，第二个是调整后的粗化单元分数。
   * @param[in]  current_n_cells 当前单元格数量。      @param[in]
   * max_n_cells 最大的细胞数。如果当前细胞数 @p current_n_cells
   * 已经超过最大细胞数 @p
   * max_n_cells，细胞的细化分数将被设置为零，细胞的粗化分数将被调整以减少细胞数至@
   * max_n_cells。如果细胞数仅在精炼时才会被超过，那么精炼和粗化分数将以相同的比例调整，以试图达到最大的细胞数。但请注意，由于
   * Triangulation::MeshSmoothing,
   * 的细化扩散，这个数字只是一个指标。这个参数的默认值是对单元格的数量没有限制。
   * @param[in]  top_fraction_of_cells 要求精炼的有效单元的比例。
   * @param[in]  bottom_fraction_of_cells
   * 要求被粗化的活动单元的比例。
   * @note  通常情况下，你不需要明确调用这个函数。将 @p
   * max_n_cells传递给函数refine_and_coarsen_fixed_number()或函数refine_and_coarsen_fixed_fraction()，如果有必要，它们会调用这个函数。
   *
   */
  template <int dim>
  std::pair<double, double>
  adjust_refine_and_coarsen_number_fraction(
    const types::global_cell_index current_n_cells,
    const types::global_cell_index max_n_cells,
    const double                   top_fraction_of_cells,
    const double                   bottom_fraction_of_cells);

  /**
   * 这个函数提供了一个策略来标记单元进行细化和粗化，目的是通过细化所有单元中的一个给定的分数来提供可预测的网格尺寸增长。    该函数接收一个细化向量 @p criteria 和两个介于0和1之间的值，表示要细化和粗化的单元的比例。它根据以下贪婪算法对单元进行标记，以便由 Triangulation::execute_coarsening_and_refinement() 进一步处理。      <ol>   <li>  根据 @p criteria.   <li> 的降序值对单元进行排序，将具有最大细化标准的 @p top_fraction_of_cells 倍 Triangulation::n_active_cells() 活动单元标记为细化。      <li>  将 @p bottom_fraction_of_cells 乘以 Triangulation::n_active_cells() 具有最小细化标准的活动单元标记为粗化。      </ol>  作为一个例子，在没有粗化的情况下，将 @p top_fraction_of_cells 设置为1/3将导致二维的单元数大约翻倍。这是因为这1/3的单元格将被其四个子单元格所取代，从而形成 $4\times \frac 13 N$ 单元格，而其余2/3的单元格则保持不变
   *
   * --因此产生了总共 $4\times \frac 13 N + \frac 23 N = 2N$
   * 个单元。
   * 在三维空间中，通过细化1/7的单元格也能达到同样的效果。因此，这些值经常被使用，因为它们确保后续网格的计算成本足够快地变得昂贵，以至于花费在粗网格上的时间部分不会太大。另一方面，这部分时间也足够小，使网格适应性在每一步中不会细化太多的单元。
   * @note  这个函数只设置粗化和细化的标志。直到你调用
   * Triangulation::execute_coarsening_and_refinement().   @param[in,out]
   * triangulation
   * 这个函数应该标记其单元格进行粗化和细化的三角形。
   * @param[in]  criteria
   * 每个网格单元的细化标准。输入不能是负数。
   * @param[in]  top_fraction_of_cells
   * 要精简的单元的比例。如果这个数字是零，没有单元会被精简。如果它等于1，结果将被标记为全局细化。
   * @param[in]  bottom_fraction_of_cells
   * 要被粗化的单元格的比例。如果这个数字为0，则没有单元会被粗化。
   * @param[in]  max_n_cells
   * 这个参数可以用来指定一个最大的细胞数。如果细化时超过这个数字，那么细化和粗化的比例将被调整，以达到最大的细胞数。但是要注意，由于
   * Triangulation::MeshSmoothing,
   * 的细化扩散，这个数字只是一个指标。这个参数的默认值是对单元格的数量没有限制。
   *
   */
  template <int dim, typename Number, int spacedim>
  void
  refine_and_coarsen_fixed_number(
    Triangulation<dim, spacedim> &triangulation,
    const Vector<Number> &        criteria,
    const double                  top_fraction_of_cells,
    const double                  bottom_fraction_of_cells,
    const unsigned int max_n_cells = std::numeric_limits<unsigned int>::max());

  /**
   * 这个函数提供了一种策略，用于标记细化和粗化的单元，目的是控制误差估计的减少。    也被称为<b>bulk criterion</b>或D&ouml;rfler标记，这个函数计算细化和粗化的阈值，使被标记为细化的 @p criteria 个单元占到总误差的一定比例。我们解释它的细化操作，粗化的工作原理与此类似。    让<i>c<sub>K</sub></i>成为<i>K</i>单元的标准。然后通过公式@f[
   * E = \sum_{K\in \cal T} c_K.
   * @f]计算总误差估计值 如果<i> 0 &lt; a &lt; 1</i>是 @p top_fraction, ，那么我们细化三角形 $\cal T$ 的最小子集 $\cal M$ ，使@f[
   * a E \le \sum_{K\in \cal M} c_K @f]
   * 该算法由refine_and_coarsen_fixed_number()中描述的贪婪算法执行。
   * @note
   * 经常使用的左右两边都有方块的公式，通过实际将<i>c<sub>K</sub></i>的方块存储在向量
   * @p criteria. 中来恢复。
   * 从实现的角度来看，这次我们确实需要对标准数组进行排序。
   * 就像上面描述的其他策略一样，这个函数只计算阈值，然后传递给refine()和coarsen()。
   * @param[in,out]  tria 三角形
   * 这个函数应该对其单元进行粗化和细化标记。
   * @param[in]  criteria 对每个网格单元计算的细化准则。
   * 输入值不能为负值。      @param[in]  top_fraction
   * 应该被细化的总估算值的分数。如果这个数字是零，没有单元会被细化。如果它等于1，结果将被标记为全局细化。
   * @param[in]  bottom_fraction
   * 粗化的估计值的一部分。如果这个数字是0，没有单元将被粗化。
   * @param[in]  max_n_cells
   * 这个参数可以用来指定一个最大的单元数。如果细化时超过这个数字，那么细化和粗化的比例将被调整，以达到最大的细胞数。但是要注意，由于
   * Triangulation::MeshSmoothing,
   * 的细化扩散，这个数字只是一个指标。这个参数的默认值是对单元格的数量没有限制。
   * @param[in]  norm_type
   * 为了确定阈值，单元格子集上的综合误差被计算为这些单元格上的准则的规范。不同类型的准则可用于此目的，目前支持其中的
   * VectorTools::NormType::L1_norm 和 VectorTools::NormType::L2_norm 。
   *
   */
  template <int dim, typename Number, int spacedim>
  void
  refine_and_coarsen_fixed_fraction(
    Triangulation<dim, spacedim> &tria,
    const Vector<Number> &        criteria,
    const double                  top_fraction,
    const double                  bottom_fraction,
    const unsigned int max_n_cells = std::numeric_limits<unsigned int>::max(),
    const VectorTools::NormType norm_type = VectorTools::NormType::L1_norm);



  /**
   * 这个函数对三角网格的单元进行标记，以达到一个最佳的网格，这个目标函数试图在网格细化时平衡减少误差和增加数值成本。具体来说，这个函数的假设是，如果你细化一个单元
   * $K$ ，其误差指标 $\eta_K$
   * 由这个函数的第二个参数提供，那么子单元上的误差（所有子单元一起）将只是
   * $2^{-\text{order}}\eta_K$ ，其中 <code>order</code>
   * 是这个函数的第三个参数。这就假设了误差只是网格上的一个局部属性，可以通过局部细化来减少误差
   *
   * *-这个假设对插值算子来说是真实的，但对通常的Galerkin投影来说不是，尽管它对椭圆问题近似是真实的，在那里Greens函数快速衰减，这里的误差不会受到其他地方太粗的网格的影响。    有了这个，我们可以定义这个函数试图优化的目标函数。让我们假设目前的网格有 $N_0$ 个单元。那么，如果我们细化误差最大的 $m$ 单元，我们期望得到（在 $d$ 空间维度上）@f[
   * N(m) = (N_0-m) + 2^d m = N_0 + (2^d-1)m
   * @f]单元（ $N_0-m$ 没有被细化，而我们细化的每个 $m$ 单元产生 $2^d$ 子单元。另一方面，在精炼 $m$ 单元时，使用上面的假设，我们预计误差将是@f[
   * \eta^\text{exp}(m) = \sum_{K, K\; \text{will not be refined}} \eta_K +
   * \sum_{K, K\; \text{will be refined}} 2^{-\text{order}}\eta_K
   * @f]，其中第一个和延伸到 $N_0-m$ 单元，第二个延伸到将被精炼的 $m$ 单元。请注意， $N(m)$ 是 $m$ 的增函数，而 $\eta^\text{exp}(m)$ 是减函数。    然后，该函数试图找到目标函数@f[
   * J(m) = N(m)^{\text{order}/d} \eta^\text{exp}(m)
   * @f]最小的单元格数量 $m$ 来标记精炼。
   * 这个函数的原理有两个方面。首先，与refine_and_coarsen_fixed_fraction()和refine_and_coarsen_fixed_number()函数相比，这个函数的特性是：如果所有的细化指标都相同（即我们已经达到了每个单元的误差平衡的网格），那么整个网格都被细化。这是基于这样的观察：具有平衡误差指标的网格是所有具有相同单元数的网格中的最优网格（即，具有最小的整体误差）。(关于这一点的证明，见R.
   * Becker, M. Braack, R. Rannacher: "Numerical simulation of laminar flames
   * at low Mach number with adaptive finite elements", Combustion Theory and
   * Modelling, Vol. 3, Nr. 3, p. 503-534 1999; and W. Bangerth, R. Rannacher:
   * "Adaptive Finite Element Methods for Differential Equations", Birkhauser,
   * 2003.)
   * 其次，该函数使用了这样的观察结果：理想情况下，误差表现为
   * $e \approx c N^{-\alpha}$ 与一些常数 $\alpha$
   * 的关系，这些常数取决于维度和有限元程度。它应该
   *
   * - 给出最佳的网格细化
   *
   * - 不太依赖于解的规则性，因为它是基于这样的想法，即所有的奇异点都可以通过细化解决。网格细化是基于我们要使 $c=e N^\alpha$ 变小的想法。这与上面的函数 $J(m)$ 相对应。
   * @note  这个函数最初是由Thomas
   * Richter实现的。它遵循T.Richter, "Parallel Multigrid Method for
   * Adaptive Finite Elements with Application to 3D Flow Problems", PhD
   * thesis, University of Heidelberg,
   * 2005中描述的策略。特别是见第4.3节，第42-43页。
   *
   */
  template <int dim, typename Number, int spacedim>
  void
  refine_and_coarsen_optimize(Triangulation<dim, spacedim> &tria,
                              const Vector<Number> &        criteria,
                              const unsigned int            order = 2);

  /**
   * 标记所有 @p criteria 中数值超过 @p
   * 阈值的网格单元进行细化，但只标记到 @p max_to_mark
   * 单元。    向量 @p criteria
   * 包含每个活动单元的非负值，按照
   * Triangulation::active_cell_iterator. 的规范顺序排序。
   * 这些单元只被标记为细化，它们实际上没有被细化。
   * 要做到这一点，你必须调用
   * Triangulation::execute_coarsening_and_refinement().
   * 这个函数不实现细化策略，它更像是实际策略的一个辅助函数。
   *
   */
  template <int dim, typename Number, int spacedim>
  void
  refine(Triangulation<dim, spacedim> &tria,
         const Vector<Number> &        criteria,
         const double                  threshold,
         const unsigned int max_to_mark = numbers::invalid_unsigned_int);

  /**
   * 标记所有 @p criteria 中的值小于 @p
   * 阈值的网格单元进行粗化。    向量 @p criteria
   * 包含每个活动单元的非负值，按照
   * Triangulation::active_cell_iterator. 的规范顺序排序。
   * 这些单元只被标记为粗化，它们实际上没有被粗化。要做到这一点，你必须调用
   * Triangulation::execute_coarsening_and_refinement().
   * 这个函数并不实现细化策略，它更像是实际策略的一个辅助函数。
   *
   */
  template <int dim, typename Number, int spacedim>
  void
  coarsen(Triangulation<dim, spacedim> &tria,
          const Vector<Number> &        criteria,
          const double                  threshold);

  /**
   * 如果带有单元格标准的向量包含负值，则抛出一个异常。
   *
   */
  DeclException0(ExcNegativeCriteria);

  /**
   * 其中一个阈值参数引起麻烦。或者细化和粗化的阈值重叠了。
   *
   */
  DeclException0(ExcInvalidParameterValue);
} // namespace GridRefinement



DEAL_II_NAMESPACE_CLOSE

#endif // dealii_grid_refinement_h


