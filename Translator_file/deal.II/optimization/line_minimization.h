//include/deal.II-translator/optimization/line_minimization_0.txt
//-----------------------------------------------------------
//
//    Copyright (C) 2018 - 2020 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//---------------------------------------------------------------

#ifndef dealii_line_minimization_h
#define dealii_line_minimization_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/std_cxx17/optional.h>
#include <deal.II/base/utilities.h>

#include <deal.II/numerics/history.h>

#include <fstream>
#include <string>


DEAL_II_NAMESPACE_OPEN

/**
 * 一个命名空间，用于与线上最小化有关的各种算法。
 *
 *
 */
namespace LineMinimization
{
  /**
   * 给出  $x\_low$  和  $x\_hi$  连同函数值  $f(x\_low)$  和
   * $f(x\_hi)$  以及梯度  $g(x\_low)$
   * ，返回二次插值函数的局部最小化器。
   * 返回类型是可选的，以便与可能没有给定参数的解的类似函数相适应。
   *
   */
  template <typename NumberType>
  std_cxx17::optional<NumberType>
  quadratic_fit(const NumberType x_low,
                const NumberType f_low,
                const NumberType g_low,
                const NumberType x_hi,
                const NumberType f_hi);

  /**
   * 给出  $x\_low$  和  $x\_hi$  以及这些点的函数值  $f(x\_low)$
   * 和  $f(x\_hi)$  及其梯度 (  $g(x\_low)*g(x\_hi) < 0$  )
   * ，返回三次插值函数的局部最小值（即三次插值函数达到其最小值的位置）。
   * 返回类型是可选的，因为实值解可能不存在。
   *
   */
  template <typename NumberType>
  std_cxx17::optional<NumberType>
  cubic_fit(const NumberType x_low,
            const NumberType f_low,
            const NumberType g_low,
            const NumberType x_hi,
            const NumberType f_hi,
            const NumberType g_hi);

  /**
   * 找出经过  $f\_low=f(x\_low)$  ,  $f\_hi=f(x\_hi)$  和
   * $f\_rec(x\_rec)$  等点，并在  $x\_low$  处有导数  $g\_low$
   * 的三次方多项式的最小化。
   * 返回类型是可选的，因为实值的解决方案可能不存在。
   *
   */
  template <typename NumberType>
  std_cxx17::optional<NumberType>
  cubic_fit_three_points(const NumberType x_low,
                         const NumberType f_low,
                         const NumberType g_low,
                         const NumberType x_hi,
                         const NumberType f_hi,
                         const NumberType x_rec,
                         const NumberType f_rec);

  /**
   * 返回使用函数值  @p f_low  ,  @p f_hi  , 和  @p f_rec[0]
   * 在三个点  @p x_low  ,  @p x_hi  , 和  @p x_rec[0]
   * 上的多项式的最小化，以及在两个点  @p g_low  和  @p g_hi.
   * 上的导数。
   * 这个函数将首先尝试执行cubic_fit()。如果它不成功，或者如果最小值不在提供的
   * @p bounds,
   * 范围内，将执行二次拟合（）。如果quadratic_fit()失败，该函数将退回到一个二分法。
   *
   */
  template <typename NumberType>
  NumberType
  poly_fit(const NumberType                        x_low,
           const NumberType                        f_low,
           const NumberType                        g_low,
           const NumberType                        x_hi,
           const NumberType                        f_hi,
           const NumberType                        g_hi,
           const FiniteSizeHistory<NumberType> &   x_rec,
           const FiniteSizeHistory<NumberType> &   f_rec,
           const FiniteSizeHistory<NumberType> &   g_rec,
           const std::pair<NumberType, NumberType> bounds);

  /**
   * 与poly_fit()相同，但是进行三点的三次拟合（见cubic_fit_three_points()
   * ）。
   *
   */
  template <typename NumberType>
  NumberType
  poly_fit_three_points(const NumberType                        x_low,
                        const NumberType                        f_low,
                        const NumberType                        g_low,
                        const NumberType                        x_hi,
                        const NumberType                        f_hi,
                        const NumberType                        g_hi,
                        const FiniteSizeHistory<NumberType> &   x_rec,
                        const FiniteSizeHistory<NumberType> &   f_rec,
                        const FiniteSizeHistory<NumberType> &   g_rec,
                        const std::pair<NumberType, NumberType> bounds);


  /**
   * 在 $(0,max]$ 中进行具有强沃尔夫条件\f[ f(\alpha) \le f(0) +
   * \alpha \mu f'(0) \\ |f'(\alpha)| \le \eta |f'(0)|
   * \f]的线搜索，使用一维函数 @p func 与函数 @p interpolate
   * 结合，根据其两端的函数值和导数从区间中选择一个新点。
   * 参数 @p a1 是对第一步的试估计。
   * 插值可以使用poly_fit()或poly_fit_three_points()，或任何其他具有类似签名的函数来完成。
   * 该函数实现了第34-35页中的算法2.6.2和2.6.4。
   * @code{.bib}
   * @book{Fletcher2013,
   * title     = {Practical methods of optimization},
   * publisher = {John Wiley \& Sons},
   * year      = {2013},
   * author    = {Fletcher, Roger},
   * isbn      = {978-0-471-49463-8},
   * doi       = {10.1002/9781118723203},
   * }
   * @endcode
   * 这些都是*@code{.bib}
   * @book{Fletcher2013,
   * title     = {Practical methods of optimization},
   * publisher = {John Wiley \& Sons},
   * year      = {2013},
   * author    = {Fletcher, Roger},
   * isbn      = {978-0-471-49463-8},
   * doi       = {10.1002/9781118723203},
   * }
   * @endcode中第60-61页的算法3.5和3.6的微小变化。
   * @code{.bib}
   * @book{Nocedal2006,
   * title     = {Numerical Optimization},
   * publisher = {Springer New York},
   * year      = {2006},
   * author    = {Jorge Nocedal and S. Wright},
   * address   = {233 Spring Street, New York, NY 10013, USA},
   * doi       = {10.1007/978-0-387-40065-5},
   * }
   * @endcode
   * 它包括一个包围阶段和一个缩放阶段，其中使用了 @p
   * interpolate 。    两个使用的例子可能如下。
   * 在第一个例子中，我们希望找到函数  $100 x^4 + (1-x)^2$
   * 的最小值。为了找到近似的解决方案，我们可以使用直线搜索和多项式拟合曲线的方法来执行以下步骤。
   * @code
   * auto func = [](const double x)
   * {
   *   const double f = 100. std::pow(x, 4) + std::pow(1.
   *
   * - x, 2); // Value
   *   const double g = 400. std::pow(x, 3)
   *
   * - 2. (1.
   *
   * - x); // Gradient
   *   return std::make_pair(f, g);
   * };
   *
   * const auto fg0 = func(0);
   * const auto res = LineMinimization::line_search<double>(
   *   func,
   *   fg0.first, fg0.second,
   *   LineMinimization::poly_fit<double>,
   *   0.1, 0.1, 0.01, 100, 20);
   *
   * const double approx_solution = res.first;
   * @endcode
   * 在第二个例子中，我们希望在一个非线性有限元问题的背景下进行直线搜索。下面是一个非优化的反向跟踪算法的实现，当负载步长过大时，它可能是有用的。下面说明了在全局非线性求解器范围内利用该方案的必要基本步骤。
   * @code
   * // Solve some incremental linear system
   * const Vector<double> newton_update = solver_linear_system(...);
   *
   * // Now we check to see if the suggested Newton update is a good one.
   * // First we define what it means to perform linesearch in the context of
   * // this incremental nonlinear finite element problem.
   * auto ls_minimization_function = [&](const double step_size)
   * {
   *   // Scale the full Newton update by the proposed line search step size.
   *   Vector<double> newton_update_trial(newton_update);
   *   newton_update_trial= step_size;
   *   // Ensure that the Dirichlet constraints are correctly applied,
   *   // irrespective of the step size
   *   constraints.distribute(newton_update_trial);
   *   // Now add the constribution from the previously accepted solution
   *   // history.
   *   const Vector<double> solution_total_trial =
   *     get_solution_total(newton_update_trial);
   *
   *   // Recompute the linear system based on the trial newton update
   *   Vector<double> system_rhs (...);
   *   SparseMatrix<double> tangent_matrix (...);
   *   assemble_linear_system(
   *     tangent_matrix, system_rhs, solution_total_trial);
   *   Vector<double> residual_trial (system_rhs);
   *   residual_trial=
   *
   * -1.0; // Residual =
   *
   * -RHS
   *
   *   // Negelect the constrained entries in the consideration
   *   // of the function (value and gradient) to be minimized.
   *   constraints.set_zero(residual_trial);
   *
   *   // Here we compute the function value according to the text given in
   *   // section 5.1.4 of Wriggers, P., "Nonlinear finite element methods",
   *   // 2008.
   *   // The function value correspeonds to equ. 5.11 on p159.
   *   const double f = 0.5 (residual_trial residual_trial); // Value
   *
   *   // However, the corresponding gradient given in eq 5.14 is wrong. The
   *   // suggested result
   *   // const double g =
   *
   * -(residual_0*residual_trial);
   *   // should actually be
   *   // g = G(V + alpha*delta)*[ K(V + alpha*delta)*delta.
   *   Vector<double> tmp;
   *   tmp.reinit(newton_update);
   *   tangent_matrix.vmult(tmp, newton_update);
   *   const double g = tmp residual_trial; // Gradient
   *
   *   return std::make_pair(f, g);
   * };
   *
   * // Next we can write a function to determine if taking the full Newton
   * // step is a good idea or not (i.e. if it offers good convergence
   * // characterisics). This function calls the one we defined above,
   * // and actually only performs the line search if an early exit
   * // criterion is not met.
   * auto perform_linesearch = [&]()
   * {
   *   const auto res_0 = ls_minimization_function(0.0);
   *   Assert(res_0.second < 0.0,
   *          ExcMessage("Gradient should be negative. Current value: " +
   *                      std::to_string(res_0.second)));
   *   const auto res_1 = ls_minimization_function(1.0);
   *
   *   // Check to see if the minimum lies in the interval [0,1] through the
   *   // values of the gradients at the limit points.
   *   // If it does not, then the full step is accepted. This is discussed by
   *   // Wriggers in the paragraph after equ. 5.14.
   *   if (res_0.second res_1.second > 0.0)
   *     return 1.0;
   *
   *   // The values for eta, mu are chosen such that more strict convergence
   *   // conditions are enforced.
   *   // They should be adjusted according to the problem requirements.
   *   const double a1        = 1.0;
   *   const double eta       = 0.5;
   *   const double mu        = 0.49;
   *   const double a_max     = 1.25;
   *   const double max_evals = 20;
   *   const auto   res = LineMinimization::line_search<double>(
   *     ls_minimization_function,
   *     res_0.first, res_0.second,
   *     LineMinimization::poly_fit<double>,
   *     a1, eta, mu, a_max, max_evals));
   *
   *   return res.first; // Final stepsize
   * };
   *
   * // Finally, we can perform the line search and adjust the Newton update
   * // accordingly.
   * const double linesearch_step_size = perform_linesearch();
   * if (linesearch_step_size != 1.0)
   * {
   *   newton_update= linesearch_step_size;
   *   constraints.distribute(newton_update);
   * }
   * @endcode
   * @param  func 一个一维函数，返回给定点的值和导数。
   * @param  f0 原点的函数值。    @param  g0
   * 在原点的函数导数。    @param  interpolate
   * 一个决定在缩放阶段如何进行内插的函数。它接受当前区间/括号（
   * $f\_low$  ,  $f\_hi$
   * ）的值和导数，以及最多五个先前步骤的值和导数。返回的值应在给定的范围内提供。
   * @param  a1 包围阶段的初始试验步骤。    @param  eta
   * 第二个沃尔夫条件（曲率条件）的一个参数。    @param
   * mu 第一个沃尔夫条件中的一个参数（充分减少）。
   * @param  a_max 允许的最大步长。    @param  max_evaluations
   * 允许的最大函数求值数。    @param  debug_output
   * 将额外的调试信息输出到 <code>deallog</code>
   * 静态对象的标志。    @return  该函数返回步长和函数 @p
   * func 被调用的次数。
   *
   */
  template <typename NumberType>
  std::pair<NumberType, unsigned int>
  line_search(
    const std::function<std::pair<NumberType, NumberType>(const NumberType x)>
      &              func,
    const NumberType f0,
    const NumberType g0,
    const std::function<
      NumberType(const NumberType                        x_low,
                 const NumberType                        f_low,
                 const NumberType                        g_low,
                 const NumberType                        x_hi,
                 const NumberType                        f_hi,
                 const NumberType                        g_hi,
                 const FiniteSizeHistory<NumberType> &   x_rec,
                 const FiniteSizeHistory<NumberType> &   f_rec,
                 const FiniteSizeHistory<NumberType> &   g_rec,
                 const std::pair<NumberType, NumberType> bounds)> &interpolate,
    const NumberType                                               a1,
    const NumberType                                               eta = 0.9,
    const NumberType                                               mu  = 0.01,
    const NumberType   a_max           = std::numeric_limits<NumberType>::max(),
    const unsigned int max_evaluations = 20,
    const bool         debug_output    = false);


  // -------------------  inline and template functions ----------------


#ifndef DOXYGEN


  template <typename NumberType>
  std_cxx17::optional<NumberType>
  quadratic_fit(const NumberType x1,
                const NumberType f1,
                const NumberType g1,
                const NumberType x2,
                const NumberType f2)
  {
    Assert(x1 != x2, ExcMessage("Point are the same"));
    const NumberType denom = (2. * g1 * x2 - 2. * g1 * x1 - 2. * f2 + 2. * f1);
    if (denom == 0)
      return {};
    else
      return (g1 * (x2 * x2 - x1 * x1) + 2. * (f1 - f2) * x1) / denom;
  }



  template <typename NumberType>
  std_cxx17::optional<NumberType>
  cubic_fit(const NumberType x1,
            const NumberType f1,
            const NumberType g1,
            const NumberType x2,
            const NumberType f2,
            const NumberType g2)
  {
    Assert(x1 != x2, ExcMessage("Points are the same"));
    const NumberType beta1 = g1 + g2 - 3. * (f1 - f2) / (x1 - x2);
    const NumberType s     = beta1 * beta1 - g1 * g2;
    if (s < 0)
      return {};

    const NumberType beta2 = std::sqrt(s);
    const NumberType denom =
      x1 < x2 ? g2 - g1 + 2. * beta2 : g1 - g2 + 2. * beta2;
    if (denom == 0.)
      return {};

    return x1 < x2 ? x2 - (x2 - x1) * (g2 + beta2 - beta1) / denom :
                     x1 - (x1 - x2) * (g1 + beta2 - beta1) / denom;
  }



  template <typename NumberType>
  std_cxx17::optional<NumberType>
  cubic_fit_three_points(const NumberType x1,
                         const NumberType f1,
                         const NumberType g1,
                         const NumberType x2,
                         const NumberType f2,
                         const NumberType x3,
                         const NumberType f3)
  {
    Assert(x1 != x2, ExcMessage("Points are the same"));
    Assert(x1 != x3, ExcMessage("Points are the same"));
    // f(x) = A *(x-x1)^3 + B*(x-x1)^2 + C*(x-x1) + D
    // =>
    // D = f1
    // C = g1

    // the rest is a system of 2 equations:

    const NumberType x2_shift = x2 - x1;
    const NumberType x3_shift = x3 - x1;
    const NumberType r1       = f2 - f1 - g1 * x2_shift;
    const NumberType r2       = f3 - f1 - g1 * x3_shift;
    const NumberType denom =
      std::pow(x2_shift * x3_shift, 2) * (x2_shift - x3_shift);
    if (denom == 0.)
      return {};

    const NumberType A =
      (r1 * std::pow(x3_shift, 2) - r2 * std::pow(x2_shift, 2)) / denom;
    const NumberType B =
      (r2 * std::pow(x2_shift, 3) - r1 * std::pow(x3_shift, 3)) / denom;
    const NumberType &C = g1;

    // now get the minimizer:
    const NumberType radical = B * B - A * C * 3;
    if (radical < 0)
      return {};

    return x1 + (-B + std::sqrt(radical)) / (A * 3);
  }



  template <typename NumberType>
  NumberType
  poly_fit(const NumberType x1,
           const NumberType f1,
           const NumberType g1,
           const NumberType x2,
           const NumberType f2,
           const NumberType g2,
           const FiniteSizeHistory<NumberType> &,
           const FiniteSizeHistory<NumberType> &,
           const FiniteSizeHistory<NumberType> &,
           const std::pair<NumberType, NumberType> bounds)
  {
    Assert(bounds.first < bounds.second, ExcMessage("Incorrect bounds"));

    // Similar to scipy implementation but we fit based on two points
    // with their gradients and do bisection on bounds.
    // https://github.com/scipy/scipy/blob/v1.0.0/scipy/optimize/linesearch.py#L555-L563

    // First try cubic interpolation
    std_cxx17::optional<NumberType> res = cubic_fit(x1, f1, g1, x2, f2, g2);
    if (res && *res >= bounds.first && *res <= bounds.second)
      return *res;

    // cubic either fails or outside of safe region, do quadratic:
    res = quadratic_fit(x1, f1, g1, x2, f2);
    if (res && *res >= bounds.first && *res <= bounds.second)
      return *res;

    // quadratic either failed or outside of safe region. Do bisection
    // on safe region
    return (bounds.first + bounds.second) * 0.5;
  }



  template <typename NumberType>
  NumberType
  poly_fit_three_points(const NumberType                     x1,
                        const NumberType                     f1,
                        const NumberType                     g1,
                        const NumberType                     x2,
                        const NumberType                     f2,
                        const NumberType                     g2,
                        const FiniteSizeHistory<NumberType> &x_rec,
                        const FiniteSizeHistory<NumberType> &f_rec,
                        const FiniteSizeHistory<NumberType> &  /*g_rec*/ ,
                        const std::pair<NumberType, NumberType> bounds)
  {
    Assert(bounds.first < bounds.second, ExcMessage("Incorrect bounds"));
    AssertDimension(x_rec.size(), f_rec.size());

    // Same as scipy implementation where cubic fit is using 3 points
    // https://github.com/scipy/scipy/blob/v1.0.0/scipy/optimize/linesearch.py#L555-L563

    // First try cubic interpolation after first iteration
    std_cxx17::optional<NumberType> res =
      x_rec.size() > 0 ?
        cubic_fit_three_points(x1, f1, g1, x2, f2, x_rec[0], f_rec[0]) :
        std_cxx17::optional<NumberType>{};
    if (res && *res >= bounds.first && *res <= bounds.second)
      return *res;

    // cubic either fails or outside of safe region, do quadratic:
    res = quadratic_fit(x1, f1, g1, x2, f2);
    if (res && *res >= bounds.first && *res <= bounds.second)
      return *res;

    // quadratic either failed or outside of safe region. Do bisection
    // on safe region
    return (bounds.first + bounds.second) * 0.5;
  }



  template <typename NumberType>
  std::pair<NumberType, unsigned int>
  line_search(
    const std::function<std::pair<NumberType, NumberType>(const NumberType x)>
      &              func,
    const NumberType f0,
    const NumberType g0,
    const std::function<
      NumberType(const NumberType                        x_low,
                 const NumberType                        f_low,
                 const NumberType                        g_low,
                 const NumberType                        x_hi,
                 const NumberType                        f_hi,
                 const NumberType                        g_hi,
                 const FiniteSizeHistory<NumberType> &   x_rec,
                 const FiniteSizeHistory<NumberType> &   f_rec,
                 const FiniteSizeHistory<NumberType> &   g_rec,
                 const std::pair<NumberType, NumberType> bounds)> &choose,
    const NumberType                                               a1,
    const NumberType                                               eta,
    const NumberType                                               mu,
    const NumberType                                               a_max,
    const unsigned int max_evaluations,
    const bool         debug_output)
  {
    // Note that scipy use dcsrch() from Minpack2 Fortran lib for line search
    Assert(mu < 0.5 && mu > 0, ExcMessage("mu is not in (0,1/2)."));
    Assert(eta < 1. && eta > mu, ExcMessage("eta is not in (mu,1)."));
    Assert(a_max > 0, ExcMessage("max is not positive."));
    Assert(a1 > 0 && a1 <= a_max, ExcMessage("a1 is not in (0,max]."));
    Assert(g0 < 0, ExcMessage("Initial slope is not negative"));

    // Growth parameter for bracketing phase:
    // 1 < tau1
    const NumberType tau1 = 9.;
    // shrink parameters for sectioning phase to prevent ai from being
    // arbitrary close to the extremes of the interval.
    // 0 < tau2 < tau3 <= 1/2
    // tau2 <= eta is advisable
    const NumberType tau2 = 0.1; // bound for closeness to a_lo
    const NumberType tau3 = 0.5; // bound for closeness to a_hi

    const NumberType g0_abs = std::abs(g0);
    const NumberType f_min  = f0 + a_max * mu * g0;

    // return True if the first Wolfe condition (sufficient decrease) is
    // satisfied
    const auto w1 = [&](const NumberType a, const NumberType f) {
      return f <= f0 + a * mu * g0;
    };

    // return True if the second Wolfe condition (curvature condition) is
    // satisfied
    const auto w2 = [&](const NumberType g) {
      return std::abs(g) <= eta * g0_abs;
    };

    // Bracketing phase (Algorithm 2.6.2): look for a non-trivial interval
    // which is known to contain an interval of acceptable points.
    // We adopt notation of Noceal.
    const NumberType x    = std::numeric_limits<NumberType>::signaling_NaN();
    NumberType       a_lo = x, f_lo = x, g_lo = x;
    NumberType       a_hi = x, f_hi = x, g_hi = x;
    NumberType       ai = x, fi = x, gi = x;

    // count function calls in i:
    unsigned int i = 0;
    {
      NumberType f_prev, g_prev, a_prev;
      ai     = a1;
      f_prev = f0;
      g_prev = g0;
      a_prev = 0;

      while (i < max_evaluations)
        {
          const auto fgi = func(ai);
          fi             = fgi.first;
          gi             = fgi.second;
          i++;

          if (debug_output)
            deallog << "Bracketing phase: " << i << std::endl
                    << ai << " " << fi << " " << gi << " " << w1(ai, fi) << " "
                    << w2(gi) << " " << f_min << std::endl;

          // first check if we can stop bracketing or the whole line search:
          if (fi <= f_min || ai == a_max)
            {
              if (debug_output)
                deallog << "Reached the maximum step size." << std::endl;
              return std::make_pair(ai, i);
            }

          if (!w1(ai, fi) ||
              (fi >= f_prev && i > 1)) // violate first Wolfe or not descending
            {
              a_lo = a_prev;
              f_lo = f_prev;
              g_lo = g_prev;

              a_hi = ai;
              f_hi = fi;
              g_hi = gi;
              break; // end bracketing
            }

          if (w2(gi)) // satisfies both Wolfe conditions
            {
              if (debug_output)
                deallog << "Satisfied both Wolfe conditions during Bracketing."
                        << std::endl;

              Assert(w1(ai, fi), ExcInternalError());
              return std::make_pair(ai, i);
            }

          if (gi >= 0) // not descending
            {
              a_lo = ai;
              f_lo = fi;
              g_lo = gi;

              a_hi = a_prev;
              f_hi = f_prev;
              g_hi = g_prev;
              break; // end bracketing
            }

          // extrapolation step with the bounds
          const auto bounds =
            std::make_pair(2. * ai - a_prev,
                           std::min(a_max, ai + tau1 * (ai - a_prev)));

          a_prev = ai;
          f_prev = fi;
          g_prev = gi;

          // NOTE: Fletcher's 2.6.2 includes optional extrapolation, we
          // simply take the upper bound
          // Scipy increases by factor of two:
          // https://github.com/scipy/scipy/blob/v1.0.0/scipy/optimize/linesearch.py#L447
          ai = bounds.second;
        }
    }

    AssertThrow(
      i < max_evaluations,
      ExcMessage(
        "Could not find the initial bracket within the given number of iterations."));

    // Check properties of the bracket (Theorem 3.2 in More and Thuente, 94
    // and Eq. 2.6.3 in Fletcher 2013

    // FIXME: these conditions are actually violated for Fig3 and a1=10^3 in
    // More and Thorenton, 94.

    /* Assert((f_lo < f_hi) && w1(a_lo, f_lo), ExcInternalError()) ; Assert((a_hi
* 
* - a_lo) g_lo < 0) && !w2(g_lo), ExcInternalError()); Assert((w1(a_hi, f_hi) || f_hi >= f_lo), ExcInternalError()。  
*
*/

    // keep short history of last points to improve interpolation
    FiniteSizeHistory<NumberType> a_rec(5), f_rec(5), g_rec(5);
    // if neither a_lo nor a_hi are zero:
    if (std::abs(a_lo) > std::numeric_limits<NumberType>::epsilon() &&
        std::abs(a_hi) > std::numeric_limits<NumberType>::epsilon())
      {
        a_rec.add(0);
        f_rec.add(f0);
        g_rec.add(g0);
      }

    // Now sectioning phase: we allow both [a_lo, a_hi] and [a_hi, a_lo]
    while (i < max_evaluations)
      {
        const NumberType a_lo_safe = a_lo + tau2 * (a_hi - a_lo);
        const NumberType a_hi_safe = a_hi - tau3 * (a_hi - a_lo);
        const auto       bounds    = std::minmax(a_lo_safe, a_hi_safe);

        ai = choose(
          a_lo, f_lo, g_lo, a_hi, f_hi, g_hi, a_rec, f_rec, g_rec, bounds);

        const std::pair<NumberType, NumberType> fgi = func(ai);
        fi                                          = fgi.first;
        gi                                          = fgi.second;
        i++;

        if (debug_output)
          deallog << "Sectioning phase: " << i << std::endl
                  << a_lo << " " << f_lo << " " << g_lo << " " << w1(a_lo, f_lo)
                  << " " << w2(g_lo) << std::endl
                  << a_hi << " " << f_hi << " " << g_hi << " " << w1(a_hi, f_hi)
                  << " " << w2(g_hi) << std::endl
                  << ai << " " << fi << " " << gi << " " << w1(ai, fi) << " "
                  << w2(gi) << std::endl;

        if (!w1(ai, fi) || fi >= f_lo)
          // take [a_lo, ai]
          {
            a_rec.add(a_hi);
            f_rec.add(f_hi);
            g_rec.add(g_hi);

            a_hi = ai;
            f_hi = fi;
            g_hi = gi;
          }
        else
          {
            if (w2(gi)) // satisfies both wolf
              {
                if (debug_output)
                  deallog << "Satisfied both Wolfe conditions." << std::endl;
                Assert(w1(ai, fi), ExcInternalError());
                return std::make_pair(ai, i);
              }

            if (gi * (a_hi - a_lo) >= 0)
              // take [ai, a_lo]
              {
                a_rec.add(a_hi);
                f_rec.add(f_hi);
                g_rec.add(g_hi);

                a_hi = a_lo;
                f_hi = f_lo;
                g_hi = g_lo;
              }
            else
              // take [ai, a_hi]
              {
                a_rec.add(a_lo);
                f_rec.add(f_lo);
                g_rec.add(g_lo);
              }

            a_lo = ai;
            f_lo = fi;
            g_lo = gi;
          }
      }

    // if we got here, we could not find the solution
    AssertThrow(
      false,
      ExcMessage(
        "Could not could complete the sectioning phase within the given number of iterations."));
    return std::make_pair(std::numeric_limits<NumberType>::signaling_NaN(), i);
  }

#endif // DOXYGEN

} // namespace LineMinimization

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_line_minimization_h


