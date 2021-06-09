//include/deal.II-translator/numerics/smoothness_estimator_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2021 by the deal.II authors
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

#ifndef dealii_smoothness_estimator_h
#define dealii_smoothness_estimator_h


#include <deal.II/base/config.h>

#include <deal.II/fe/component_mask.h>

#include <deal.II/numerics/vector_tools.h>

#include <functional>
#include <vector>


DEAL_II_NAMESPACE_OPEN


// forward declarations
#ifndef DOXYGEN
template <typename Number>
class Vector;

template <int dim, int spacedim>
class DoFHandler;

namespace FESeries
{
  template <int dim, int spacedim>
  class Fourier;
  template <int dim, int spacedim>
  class Legendre;
} // namespace FESeries

namespace hp
{
  template <int dim, int spacedim>
  class FECollection;
} // namespace hp
#endif


/**
 * 一个用于hp-adaptive FEM的各种平滑度估计策略的命名空间。
 * 平滑度估计是决定一个具有大误差估计的单元是否应该进行h-或p-精简的一种策略。典型的策略是在解特别光滑的情况下决定增加单元上的多项式度数，而如果单元上的解是奇异的，在某些导数上有扭结，或在其他方面不是特别光滑，则会细化网格。所有这些策略都依赖于识别一个函数在一个给定单元上的
 * "平滑 "程度的方法。
 *
 *
 */
namespace SmoothnessEstimator
{
  /**
   * 基于Legendre扩展系数衰减的平滑度估计策略。    在一维中，单元格 $K$ 上的多项式程度 $p$ 的有限元解可以写成@f{eqnarray*}
   *  u_h(x) &=& \sum_j u_j \varphi_j (x) \\
   *  u_{h, k}(x) &=& \sum_{k=0}^{p} a_k \widetilde P_k (x),
   *  \quad a_k = \sum_j {\cal L}_{k,j} u_j
   * @f} 。
   * *其中 $u_j$ 为自由度， $\varphi_j$ 为相应的形状函数。  $\{\widetilde P_k(x)\}$ 是单元格 $K$ 上的 Legendre多项式。  $a_k$ 和 ${\cal L}_{k,j}$ 是每个形状函数的Legendre展开的系数和变换矩阵。由于实际原因，我们将只对参考单元格 $\hat K$ 进行这些矩阵和系数的计算。我们只需要通过这种方式计算一次变换矩阵。然而，结果只适用于从参考单元到实际单元的映射是仿射的情况。我们使用类  FESeries::Legendre  来确定所有系数  $a_k$  。    当且仅当一个函数的Legendre扩展系数以（见 @cite eibner2007hp ）@f[
   * |a_k| \sim c \, \exp(-\sigma k)
   * @f]的方式衰减时，该函数是可分析的，即可由幂级数表示，我们通过对 $k=0,\ldots,p$ 进行@f[
   * \ln |a_k| \sim C
   *
   * - \sigma k @f]的线性回归拟合来确定其衰减率 $\sigma$  ，
   * $p$  是有限元素的多项式程度。  衰减的速度 $\sigma$
   * 可以用来估计平滑度。例如，实施hp-精简准则的一个策略是，如果
   * $\sigma>1$ （见 @cite mavriplis1994hp ），则进行p-精简。
   *
   */
  namespace Legendre
  {
    /**
     * 在这个针对更高维度的估计策略的变体中，我们将考虑描述Legendre多项式 $\widetilde P_{\bf k}$ 的所有模式向量 $\bf k$ ，并对所有系数一次性进行一次最小二乘拟合。如果有多个系数对应相同的模式绝对值 $\|{\bf k}\|_1$ ，我们取其中的最大值。因此，对@f{eqnarray*}
     * \widetilde P_{\bf k}({\bf x}) &=&
     *   \widetilde P_{k_1} (x_1) \ldots \widetilde P_{k_d} (x_d) \\
     * \ln \left( \max\limits_{\|{\bf k}\|_1} |a_{\bf k}| \right) &\sim&
     *   C
     *
     * - \sigma \|{\bf k}\|_1
     * @f}进行最小二乘拟合。
     * 对于 ${\bf k}=(k_1,\ldots,k_d)$ 和 $k_i=0,\ldots,p$ ， $p$
     * 是有限元的多项式程度。        对于有限元近似 @p
     * solution, ，该函数将每个单元的衰减率写入输出向量 @p
     * smoothness_indicators.   @param  [in] fe_legendre  FESeries::Legendre
     * 对象，以计算系数。
     * 这个对象需要被初始化，以使集合中的每个有限元在每个方向上至少有
     * $p+1$ 个系数，其中 $p$ 是其多项式程度。      @param  [in]
     * dof_handler 一个DoFHandler。      @param  [in] solution
     * 一个解决方案向量。      @param  [out] smoothness_indicators
     * 一个平滑度指标的向量。      @param  [in] regression_strategy
     * 决定对绝对值相同的系数子集 $\mathbf{k}$ 使用哪种准则
     * $\|{\bf k}\|_1$  。默认是 VectorTools::Linfty_norm
     * ，用于最大程度的近似。      @param  [in]
     * smallest_abs_coefficient
     * 用于线性回归的最小绝对值的系数。注意，某些函数的Legendre系数可能有零系数的重复模式（即对于在任何坐标方向上围绕元素中点的局部对称或反对称的函数）。因此，这个参数允许忽略线性回归拟合中的小（绝对值）系数。如果有少于两个非零系数，该单元的返回值将是
     * $\sigma=\infty$  。      @param  [in] only_flagged_cells
     * 平稳性指标通常用于决定是进行h-适应还是p-适应。所以在大多数情况下，我们只需要计算那些被标记为细化或粗化的单元的指标。这个参数控制是否考虑这个特定的子集或所有单元。默认情况下，所有细胞都将被考虑。当只考虑被标记的单元时，平滑度指标将只设置在那些被标记的单元的向量条目上；其他的将被设置为一个信号NaN。
     * 更多理论细节见  @cite mavriplis1994hp   @cite houston2005hp
     * @cite eibner2007hp  。
     *
     */
    template <int dim, int spacedim, typename VectorType>
    void
    coefficient_decay(FESeries::Legendre<dim, spacedim> &fe_legendre,
                      const DoFHandler<dim, spacedim> &  dof_handler,
                      const VectorType &                 solution,
                      Vector<float> &                    smoothness_indicators,
                      const VectorTools::NormType        regression_strategy =
                        VectorTools::Linfty_norm,
                      const double smallest_abs_coefficient = 1e-10,
                      const bool   only_flagged_cells       = false);

    /**
     * 在这个针对更高维度的估计策略的变体中，我们只考虑每个坐标方向的模式，即只考虑有一个非零项的模式向量
     * $\bf k$
     * 。我们分别对每个坐标方向进行最小二乘拟合，并取其中最低的衰减率
     * $\sigma$ 。        对于有限元近似 @p solution,
     * ，该函数将每个单元的衰减率写入输出向量 @p
     * smoothness_indicators.   @param  [in] fe_legendre FESeries::Legendre
     * 对象以计算系数。
     * 这个对象需要被初始化为在每个方向上至少有 $p+1$
     * 个系数，其中 $p$ 是要使用的最大多项式程度。
     * @param  [in] dof_handler 一个DoFHandler  @param  [in] solution
     * 一个解决方案向量  @param  [out] smoothness_indicators
     * 一个平滑度指标向量  @param  [in] coefficients_predicate
     * 一个用于在每个坐标方向选择线性回归的 Legendre 系数
     * $a_j$  ,  $j=0,\ldots,p$
     * 的谓语。用户负责更新提供给这个函数的 "标志
     * "向量。注意，它的大小是 $p+1$  ，其中 $p$
     * 是FE基础在给定元素上的多项式程度。默认实现将使用每个坐标方向上的所有Legendre系数，即把向量的所有元素设置为`true'。
     * @param  [in] smallest_abs_coefficient
     * 在每个坐标方向上用于线性回归的最小系数的绝对值。
     * 注意，某些函数的Legendre系数可能有零系数的重复模式（即对于在任何坐标方向上围绕元素中点的局部对称或反对称的函数）。因此，这个参数允许忽略线性回归拟合中的小（绝对值）系数。如果一个坐标方向上的非零系数少于两个，这个方向将被跳过。如果所有的系数都是零，这个单元的返回值将是
     * $\sigma=\infty$  。      @param  [in] only_flagged_cells
     * 平滑性指标通常用于决定是进行h-适应还是p-适应。所以在大多数情况下，我们只需要计算那些被标记为细化或粗化的单元的指标。这个参数控制是否考虑这个特定的子集或所有单元。默认情况下，所有细胞都将被考虑。当只考虑被标记的单元时，平滑度指标将只设置在那些被标记的单元的向量条目上；其他的将被设置为NaN。
     * 更多的理论细节和在deal.II库中的应用见  @cite
     * davydov2017hp  。
     *
     */
    template <int dim, int spacedim, typename VectorType>
    void
    coefficient_decay_per_direction(
      FESeries::Legendre<dim, spacedim> &fe_legendre,
      const DoFHandler<dim, spacedim> &  dof_handler,
      const VectorType &                 solution,
      Vector<float> &                    smoothness_indicators,
      const ComponentMask &coefficients_predicate   = ComponentMask(),
      const double         smallest_abs_coefficient = 1e-10,
      const bool           only_flagged_cells       = false);

    /**
     * 返回一个 FESeries::Legendre
     * 对象，用于Legendre级数展开，其默认配置用于平滑度估计目的。
     * 对于所提供的 @p fe_collection,
     * 的每个有限元，我们使用与其多项式度数相同的模式加2。这包括第一个Legendre多项式，它只是一个常数。此外，对于每个元素，我们使用高斯正交，旨在为所使用的最高阶Legendre多项式产生精确的结果。
     * 由于Legendre扩展只能在标量场上进行，该类不对矢量值的有限元进行操作，因此会抛出一个断言。然而，有限元场的每个分量可以分别被视为标量场，对其进行Legendre扩展也是可能的。为此，可选的参数
     * @p component 定义了每个有限元的哪个分量将被使用。
     * 默认值 @p component
     * 只适用于标量FEs，在这种情况下，它表示唯一的分量将被分解。对于矢量FE，必须明确提供一个非默认值。
     *
     */
    template <int dim, int spacedim>
    FESeries::Legendre<dim, spacedim>
    default_fe_series(
      const hp::FECollection<dim, spacedim> &fe_collection,
      const unsigned int component = numbers::invalid_unsigned_int);
  } // namespace Legendre



  /**
   * 基于傅里叶扩展系数衰减的平滑度估计策略。    根据定义，我们可以将单元格 $K$ 上具有多项式程度 $p$ 的有限元解的傅里叶级数展开 $a_{\bf k}$ 写成矩阵积@f{eqnarray*}
   *  u_h({\bf x}) &=& \sum_j u_j \varphi_j ({\bf x}) \\
   *  u_{h, {\bf k}}({\bf x}) &=&
   *    \sum_{{\bf k}, \|{\bf k}\|\le p} a_{\bf k} \phi_{\bf k}({\bf x}),
   *    \quad a_{\bf k} = \sum_j {\cal F}_{{\bf k},j} u_j
   * @f} 。
   * 与 $u_j$ 的自由度和 $\varphi_j$ 的相应形状函数。
   * $\{\phi_{\bf k}({\bf x}) = \exp(i \, 2 \pi \, {\bf k} \cdot {\bf x}) \}$
   * 是单元格上的指数函数  $K$  。  $a_{\bf k}$ 和 ${\cal
   * F}_{{\bf k},j}$ 是每个形状函数的傅里叶展开的系数和变换矩阵。出于实际考虑，我们将只在参考单元 $\hat K$ 上进行这些矩阵和系数的计算。我们只需要通过这种方式计算一次变换矩阵。然而，结果只适用于从参考单元到实际单元的映射是线性的。我们使用类  FESeries::Fourier  来确定所有系数  $a_{\bf k}$  。    如果单元 $K$ 上的有限元近似是希尔伯特空间 $H^s(K)$ 的一部分，那么对于我们的解决方案@f{eqnarray*}
   * \| \nabla^s u_h({\bf x}) \|_{L^2(K)}^2 &=&
   *   \int\limits_K \left| \nabla^s u_h({\bf x}) \right|^2 d{\bf x} <
   *   \infty \\
   * \| \nabla^s u_{h, {\bf k}}({\bf x}) \|_{L^2(K)}^2 &=&
   *   \int\limits_K \left| \sum\limits_{\bf k} (-i \, 2 \pi {\bf k})^s \,
   *   a_{\bf k} \, \phi_{\bf k}({\bf x}) \right|^2 d{\bf x} =
   *   (2 \pi)^{2s} \sum\limits_{\bf k} \left| a_{\bf k} \right|^2
   *   \|{\bf k}\|_2^{2s} < \infty
   * @f}的有限元和谱系表示，以下积分必须存在
   * 只有当和数在所有 $\epsilon > 0$ 中至少以@f[
   * |a_{\bf k}|^2 \|{\bf k}\|_2^{2s} \|{\bf k}\|_2^{d
   *
   * - 1} =
   *   {\cal O}\left( \|{\bf k}\|_2^{-1-\epsilon} \right)
   * @f]的阶数衰减时，该和数才是有限的。额外的因素源于这样一个事实，即由于我们对位于二维球体上的所有多指数 ${\bf k}$ 求和，我们实际上有位于每个增量 $\|{\bf k}\|_2 +
   * d\|{\bf k}\|_2$ 中的 $\|{\bf k}\|_2^{d-1}$ 模式，需要考虑到的常数。通过指数的比较，我们可以把这个条件改写为 @f[
   * |a_{\bf k}| = {\cal O}\left(\|{\bf k}\|_2^
   *   {-\left(s + \frac d2 + \epsilon \right)} \right)
   * @f] 下一步是估计这些系数随  $\|{\bf k}\|_2$  衰减的速度。因此，我们用回归系数 $\alpha$ 和 $\sigma$ 进行最小二乘拟合@f[
   * \min_{\alpha,\sigma} \frac 12 \sum_{{\bf k}, \|{\bf k}\|_2 \le p} \left(
   * |a_{\bf k}|
   *
   * - \alpha \|{\bf k}\|_2^{-\sigma}\right)^2
   * @f] 。为了简化，我们对我们的最小化问题@f[
   * \min_{\beta,\sigma} Q(\beta,\sigma) = \frac 12 \sum_{{\bf k}, \|{\bf
   * k}\|_2 \le p} \left( \ln |a_{\bf k}|
   *
   * - \beta + \sigma \ln \|{\bf k}\|_2 \right)^2, @f]应用对数，其中
   * $\beta=\ln \alpha$
   * 。现在这是一个问题，对于这个问题，最优性条件
   * $\frac{\partial Q}{\partial\beta}=0,
   * \frac{\partial Q}{\partial\sigma}=0$  , 在  $\beta,\sigma$  中是线性的。我们可以把这些条件写成如下。  @f[
   *  \left(\begin{array}{cc}
   *  \sum_{{\bf k}, \|{\bf k}\|_2 \le p} 1 &
   *  \sum_{{\bf k}, \|{\bf k}\|_2 \le p} \ln \|{\bf k}\|_2
   *  \\
   *  \sum_{{\bf k}, \|{\bf k}\|_2 \le p} \ln \|{\bf k}\|_2 &
   *  \sum_{{\bf k}, \|{\bf k}\|_2 \le p} (\ln \|{\bf k}\|_2)^2
   *  \end{array}\right)
   *  \left(\begin{array}{c}
   *  \beta \\
   *
   * -\sigma
   *  \end{array}\right)
   *  =
   *  \left(\begin{array}{c}
   *  \sum_{{\bf k}, \|{\bf k}\|_2\le p} \ln |a_{{\bf k}}|
   *  \\
   *  \sum_{{\bf k}, \|{\bf k}\|_2\le p} \ln |a_{{\bf k}}| \ln \|{\bf
   * k}\|_2 \end{array}\right)
   * @f] 解决 $\beta$ 和 $\sigma$ 只是一个线性回归拟合，为此我们将使用 FESeries::linear_regression().  虽然我们对 $\beta$ 的实际值不是特别感兴趣，但上面的公式给了我们一种方法来计算指数 $\sigma$ 的值，然后我们可以用它来确定 $u(\hat{\bf x})$ 在 $H^s(K)$  与 $s=\sigma-\frac d2$  中。衰减率 $\sigma$ 将足以作为我们的平滑度指标，并将对任何提供的解决方案的每个单元进行计算。
   * @note  在 step-27 中提供了关于这些函数使用的广泛演示。
   *
   */
  namespace Fourier
  {
    /**
     * 在这个针对更高维度的估计策略的变体中，我们将考虑描述傅里叶多项式 $P_{\bf k}$ 的所有模式向量 $\bf k$ ，并对所有系数一次性进行一次最小二乘拟合。如果有多个系数对应相同的模态绝对值 $\|\bf k\|_2$ ，我们取其中的最大值。    因此，在@f[
     * \ln \left( \max\limits_{\|{\bf k}\|_2} |a_{\bf k}| \right) \sim
     *   C
     *
     * - \sigma \ln \|{\bf k}\|_2
     * @f]上对 ${\bf k}=(k_1,\ldots,k_d)$ 和 $k_i=0,\ldots,p$ 进行最小二乘拟合， $p$ 是有限元的多项式程度。我们排除了 $\|{\bf k}\|_2=0$ 模式以避免对数的奇异性。         @p regression_strategy  参数决定了在绝对值相同的系数子集 $\bf k$ 上使用哪种准则  $\|{\bf k}\|_2$  。默认是 VectorTools::Linfty_norm ，用于最大程度的逼近。        对于提供的定义在DoFHandler  @p dof_handler, 上的求解向量 @p solution ，该函数返回一个有多少个元素的向量 @p smoothness_indicators ，其中每个元素包含估计的规则性 $\sigma$  。        必须提供一个序列扩展对象 @p fe_fourier ，它需要用与 @p dof_handler. 相同的FECollection对象构建。 参数 @p smallest_abs_coefficient 允许忽略线性回归拟合中的小（绝对值）系数。如果一个坐标方向上的非零系数少于两个，这个方向将被跳过。如果所有的系数都是零，这个单元的返回值将是  $\sigma=\infty$  。        平滑度指标通常被用来决定是进行h-适应还是p-适应。因此，在大多数情况下，我们只需要对标记为细化或粗化的单元格计算这些指标。参数 @p only_flagged_cells 控制是否考虑这个特定的子集或所有单元。默认情况下，所有单元格都将被考虑。    当只考虑被标记的单元时，平滑度指标将只设置在那些被标记的单元的向量条目上；其他的将被设置为一个信号NaN。
     *
     */
    template <int dim, int spacedim, typename VectorType>
    void
    coefficient_decay(FESeries::Fourier<dim, spacedim> &fe_fourier,
                      const DoFHandler<dim, spacedim> & dof_handler,
                      const VectorType &                solution,
                      Vector<float> &                   smoothness_indicators,
                      const VectorTools::NormType       regression_strategy =
                        VectorTools::Linfty_norm,
                      const double smallest_abs_coefficient = 1e-10,
                      const bool   only_flagged_cells       = false);

    /**
     * 在这个针对更高维度的估计策略的变体中，我们只考虑每个坐标方向的模式，即只考虑有一个非零项的模式向量
     * $\bf k$
     * 。我们分别对每个坐标方向进行最小二乘拟合，并取其中最低的衰变率
     * $\sigma$ 。         @p coefficients_predicate
     * 参数选择傅里叶系数 $a_j$  ,  $j=0,\ldots,p$
     * 用于每个坐标方向的线性回归。用户负责更新提供给这个函数的
     * "标志 "向量。注意它的大小是 $p+1$  ，其中 $p$
     * 是给定元素上的FE基础的多项式程度。默认实现将使用每个坐标方向上的所有傅里叶系数，即把向量的所有元素设置为`true'。
     * 对于提供的定义在DoFHandler  @p dof_handler, 上的求解向量
     * @p solution ，该函数返回一个有多少个元素的向量 @p
     * smoothness_indicators ，其中每个元素包含估计的规则性
     * $\sigma$  。        必须提供一个序列扩展对象 @p fe_fourier
     * ，它需要用与 @p dof_handler. 相同的FECollection对象构建。
     * 参数 @p smallest_abs_coefficient
     * 允许忽略线性回归拟合中的小（绝对值）系数。如果一个坐标方向上的非零系数少于两个，这个方向将被跳过。如果所有的系数都是零，这个单元的返回值将是
     * $\sigma=\infty$  。
     * 平滑度指标通常被用来决定是进行h-适应还是p-适应。因此在大多数情况下，我们只需要计算那些标记为细化或粗化的单元格的指标。参数
     * @p only_flagged_cells
     * 控制是否考虑这个特定的子集或所有单元。默认情况下，所有单元都将被考虑。
     * 当只考虑被标记的单元时，平滑度指标将只设置在那些被标记的单元的向量条目上；其他的将被设置为一个信号NaN。
     *
     */
    template <int dim, int spacedim, typename VectorType>
    void
    coefficient_decay_per_direction(
      FESeries::Fourier<dim, spacedim> &fe_fourier,
      const DoFHandler<dim, spacedim> & dof_handler,
      const VectorType &                solution,
      Vector<float> &                   smoothness_indicators,
      const ComponentMask &coefficients_predicate   = ComponentMask(),
      const double         smallest_abs_coefficient = 1e-10,
      const bool           only_flagged_cells       = false);

    /**
     * 返回一个用于傅里叶级数展开的 FESeries::Fourier
     * 对象，其默认配置用于平滑度估计目的。
     * 对于所提供的 @p fe_collection,
     * 的每个有限元素，我们使用与它的多项式度数相同的模式加2。此外，对于每个元素，我们使用5点高斯夸父，在每个维度上迭代最大波数，这是由我们从
     * $k = 0$ 开始的模式数量减少了一个。
     * 由于傅里叶展开只能在标量场上进行，这个类不能对矢量值的有限元进行操作，因此会抛出一个断言。然而，有限元场的每个分量都可以分别被视为标量场，对其进行傅里叶扩展也是可能的。为此，可选的参数
     * @p component 定义了每个有限元的哪个分量将被使用。
     * @p component
     * 的默认值只适用于标量FEs，在这种情况下，它表示唯一的分量将被分解。对于矢量FE，必须明确提供一个非默认值。
     *
     */
    template <int dim, int spacedim>
    FESeries::Fourier<dim, spacedim>
    default_fe_series(
      const hp::FECollection<dim, spacedim> &fe_collection,
      const unsigned int component = numbers::invalid_unsigned_int);
  } // namespace Fourier
} // namespace SmoothnessEstimator


DEAL_II_NAMESPACE_CLOSE

#endif


