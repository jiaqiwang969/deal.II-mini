//include/deal.II-translator/numerics/vector_tools_common_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#ifndef dealii_vector_tools_common_h
#define dealii_vector_tools_common_h


#include <deal.II/base/config.h>

#include <deal.II/base/patterns.h>

DEAL_II_NAMESPACE_OPEN

namespace VectorTools
{
  /**
   * 表示哪个规范/积分将由 integrate_difference()
   * 函数在每个单元格上计算，而 compute_global_error()
   * 则是针对整个域。  让 $f:\Omega \rightarrow \mathbb{R}^c$
   * 是一个具有 $c$ 分量的有限元函数，其中 $c$ 分量由 $f_c$
   * 表示， $\hat{f}$ 是参考函数（ integrate_difference()的 @p
   * fe_function 和 @p exact_solution  参数）。让 $e_c = \hat{f}_c
   *
   * - f_c$ 为二者之间的差异或误差。此外，让 $w:\Omega \rightarrow \mathbb{R}^c$ 为 integrate_difference()的 @p weight 函数，如果没有提供，则假定它等于1。最后，让 $p$ 为 @p exponent 的参数（对于 $L_p$ -norms）。    在下文中，我们用 $E_K$ 表示在单元格 $K$ 上通过 integrate_difference() 计算的局部误差，而 $E$ 是通过 compute_global_error() 计算的全局误差。请注意，积分是以通常的方式用正交法进行近似的。  @f[
   * \int_A f(x) dx \approx \sum_q f(x_q) \omega_q.
   * @f] 同样，对于一个单元上的上位数 $T$  : @f[
   * \sup_{x\in T} |f(x)| dx \approx \max_q |f(x_q)|. @f] 。
   *
   */
  enum NormType
  {
    /**
     * 函数或函数之差在每个单元格上进行积分  $K$  : @f[
     * E_K
     * = \int_K \sum_c (\hat{f}_c
     *
     * - f_c) \, w_c
     * = \int_K \sum_c e_c \, w_c
     * @f] 并求和得到 @f[
     * E = \sum_K E_K = \int_\Omega \sum_c (\hat{f}_c
     *
     * - f_c) \, w_c
     * @f] 或对于  $w \equiv 1$  : @f[
     * E = \int_\Omega (\hat{f}
     *
     * - f) = \int_\Omega e. @f]
     * 注意：这与通常所说的函数的平均值不同，是
     * $\frac{1}{|\Omega|}$
     * 的一个因素。要计算平均值，你也可以使用compute_mean_value()。最后，注意符号：如果
     * $\hat{f}=0$  ，这将计算  $f$  的平均值的负数。
     *
     */
    mean,

    /**
     * 函数的绝对值被整合。    @f[
     * E_K = \int_K \sum_c |e_c| \, w_c
     * @f]和@f[
     * E = \sum_K E_K = \int_\Omega \sum_c |e_c| w_c,
     * @f]或者，对于 $w \equiv 1$ ：@f[
     * E  = \| e \|_{L^1}. @f] 。
     *
     */
    L1_norm,

    /**
     * 函数的平方被整合，结果的平方根在每个单元格上被计算出来。    @f[
     * E_K = \sqrt{ \int_K \sum_c e_c^2 \, w_c }
     * @f]和@f[
     * E = \sqrt{\sum_K E_K^2} = \sqrt{ \int_\Omega  \sum_c e_c^2 \, w_c }
     * @f]或者，对于 $w \equiv 1$ ：@f[
     * E = \sqrt{ \int_\Omega e^2 } = \| e \|_{L^2} @f]。
     *
     */
    L2_norm,

    /**
     * 对 $p$ 的绝对值进行整合，在每个单元格上计算 $p$ 的根。指数 $p$ 是 integrate_difference() 和 compute_mean_value() 的 @p  指数参数。    @f[
     * E_K = \left( \int_K \sum_c |e_c|^p \, w_c \right)^{1/p}
     * @f]和@f[
     * E = \left( \sum_K E_K^p \right)^{1/p}
     * @f]或者，对于 $w \equiv 1$ ：@f[
     * E = \| e \|_{L^p}. @f] 。
     *
     */
    Lp_norm,

    /**
     * 的最大绝对值的函数。    @f[
     * E_K = \sup_K \max_c |e_c| \, w_c
     * @f] 和 @f[
     * E = \max_K E_K = \sup_\Omega \max_c |e_c| \, w_c
     * @f] 或者，对于  $w \equiv 1$  : @f[
     * E  = \sup_\Omega \|e\|_\infty = \| e \|_{L^\infty}. @f] 。
     *
     */
    Linfty_norm,

    /**
     * #L2_norm的梯度。    @f[
     * E_K = \sqrt{ \int_K \sum_c (\nabla e_c)^2 \, w_c }
     * @f]和@f[
     * E = \sqrt{\sum_K E_K^2} = \sqrt{ \int_\Omega \sum_c (\nabla e_c)^2 \,
     * w_c }
     * @f]或者，对于 $w \equiv 1$ ：@f[
     * E = \| \nabla e \|_{L^2}. @f] 。
     *
     */
    H1_seminorm,

    /**
     * #L2_norm的向量场的发散。函数 $f$ 预计有 $c \geq \text{dim}$ 个分量，第一个 @p dim 将被用来计算发散。    @f[
     * E_K = \sqrt{ \int_K \left( \sum_c \frac{\partial e_c}{\partial x_c} \,
     * \sqrt{w_c} \right)^2 }
     * @f]和@f[
     * E = \sqrt{\sum_K E_K^2} = \sqrt{ \int_\Omega \left( \sum_c
     * \frac{\partial e_c}{\partial x_c} \, \sqrt{w_c} \right)^2  }
     * @f]或者，对于 $w \equiv 1$ ：@f[
     * E = \| \nabla \cdot e \|_{L^2}. @f] 。
     *
     */
    Hdiv_seminorm,

    /**
     * 这个规范的平方是#L2_norm的平方加上#H1_seminorm的平方。    @f[
     * E_K = \sqrt{ \int_K \sum_c (e_c^2 + (\nabla e_c)^2) \, w_c }
     * @f]和@f[
     * E = \sqrt{\sum_K E_K^2} = \sqrt{ \int_\Omega \sum_c (e_c^2 + (\nabla
     * e_c)^2) \, w_c }
     * @f]或者，对于 $w \equiv 1$ ：@f[
     * E = \left( \| e \|_{L^2}^2 + \| \nabla e \|_{L^2}^2 \right)^{1/2}. @f]
     * 。
     *
     */
    H1_norm,

    /**
     * #Lp_norm的梯度。    @f[
     * E_K = \left( \int_K \sum_c |\nabla e_c|^p \, w_c \right)^{1/p}
     * @f]和@f[
     * E = \left( \sum_K E_K^p \right)^{1/p} = \left( \int_\Omega \sum_c
     * |\nabla e_c|^p \, w_c \right)^{1/p}
     * @f]或者，对于 $w \equiv 1$ ：@f[
     * E = \| \nabla e \|_{L^p}. @f] 。
     *
     */
    W1p_seminorm,

    /**
     * 与#H1_norm相同，但使用<i>L<sup>p</sup></i>。    @f[
     * E_K = \left( \int_K \sum_c (|e_c|^p + |\nabla e_c|^p) \, w_c
     * \right)^{1/p}
     * @f]和@f[
     * E = \left( \sum_K E_K^p \right)^{1/p} = \left( \int_\Omega \sum_c
     * (|e_c|^p + |\nabla e_c|^p) \, w_c \right)^{1/p}
     * @f]，或者，对于 $w \equiv 1$  : @f[
     * E = \left( \| e \|_{L^p}^p + \| \nabla e \|_{L^p}^p \right)^{1/p}. @f]
     * 。
     *
     */
    W1p_norm,

    /**
     * #梯度的Linfty_norm。    @f[
     * E_K = \sup_K \max_c |\nabla e_c| \, w_c
     * @f]和@f[
     * E = \max_K E_K = \sup_\Omega \max_c |\nabla e_c| \, w_c
     * @f]或者，对于 $w \equiv 1$ ：@f[
     * E = \| \nabla e \|_{L^\infty}. @f] 。
     *
     */
    W1infty_seminorm,

    /**
     * ＃Linfty_norm和＃W1infty_seminorm之和。    @f[
     * E_K = \sup_K \max_c |e_c| \, w_c + \sup_K \max_c |\nabla e_c| \, w_c.
     * @f] 全局规范没有在compute_global_error()中实现，因为不可能从值中计算出全局规范的总和  $E_K$  。作为一种变通方法，你可以分别计算全局的#Linfty_norm和#W1infty_seminorm，然后将它们相加，得到（与  $w \equiv 1$  ）。    @f[
     * E = \| e \|_{L^\infty} + \| \nabla e \|_{L^\infty}. @f]
     *
     */
    W1infty_norm
  };

  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcPointNotAvailableHere,
                   "The given point is inside a cell of a "
                   "parallel::distributed::Triangulation that is not "
                   "locally owned by this processor.");
} // namespace VectorTools

// Make sure we can use NormType with Patterns.
namespace Patterns
{
  namespace Tools
  {
    template <>
    struct Convert<VectorTools::NormType, void>
    {
      /**
       * 返回NormType的正确模式。
       *
       */
      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        return std::make_unique<Patterns::Selection>(
          "mean|L1_norm|L2_norm|Lp_norm|"
          "Linfty_norm|H1_seminorm|Hdiv_seminorm|"
          "H1_norm|W1p_seminorm|W1p_norm|"
          "W1infty_seminorm|W1infty_norm");
      }



      /**
       * 将NormType转换为一个字符串。
       *
       */
      static std::string
      to_string(const VectorTools::NormType &s,
                const Patterns::PatternBase &p =
                  *Convert<VectorTools::NormType>::to_pattern())
      {
        std::string str;
        if (s == VectorTools::mean)
          str = "mean";
        else if (s == VectorTools::L1_norm)
          str = "L1_norm";
        else if (s == VectorTools::L2_norm)
          str = "L2_norm";
        else if (s == VectorTools::Lp_norm)
          str = "Lp_norm";
        else if (s == VectorTools::Linfty_norm)
          str = "Linfty_norm";
        else if (s == VectorTools::H1_seminorm)
          str = "H1_seminorm";
        else if (s == VectorTools::Hdiv_seminorm)
          str = "Hdiv_seminorm";
        else if (s == VectorTools::H1_norm)
          str = "H1_norm";
        else if (s == VectorTools::W1p_seminorm)
          str = "W1p_seminorm";
        else if (s == VectorTools::W1infty_seminorm)
          str = "W1infty_seminorm";
        else if (s == VectorTools::W1infty_norm)
          str = "W1infty_norm";
        else if (s == VectorTools::W1p_norm)
          str = "W1p_norm";
        else
          {
            AssertThrow(false, ExcMessage("Didn't recognize a norm type."));
          }
        AssertThrow(p.match(str), ExcInternalError());
        return str;
      }


      /**
       * 将一个字符串转换为一个NormType。
       *
       */
      static VectorTools::NormType
      to_value(const std::string &          str,
               const Patterns::PatternBase &p =
                 *Convert<VectorTools::NormType>::to_pattern())
      {
        VectorTools::NormType norm = VectorTools::mean;
        AssertThrow(p.match(str),
                    ExcMessage(
                      "String " + str +
                      " cannot be converted to VectorTools::NormType"));

        if (str == "mean")
          norm = VectorTools::mean;
        else if (str == "L1_norm")
          norm = VectorTools::L1_norm;
        else if (str == "L2_norm")
          norm = VectorTools::L2_norm;
        else if (str == "Lp_norm")
          norm = VectorTools::Lp_norm;
        else if (str == "Linfty_norm")
          norm = VectorTools::Linfty_norm;
        else if (str == "H1_seminorm")
          norm = VectorTools::H1_seminorm;
        else if (str == "Hdiv_seminorm")
          norm = VectorTools::Hdiv_seminorm;
        else if (str == "H1_norm")
          norm = VectorTools::H1_norm;
        else if (str == "W1p_seminorm")
          norm = VectorTools::W1p_seminorm;
        else if (str == "W1infty_seminorm")
          norm = VectorTools::W1infty_seminorm;
        else if (str == "W1infty_norm")
          norm = VectorTools::W1infty_norm;
        else if (str == "W1p_norm")
          norm = VectorTools::W1p_norm;
        else
          {
            AssertThrow(false, ExcMessage("Didn't recognize a norm type."));
          }
        return norm;
      }
    };
  } // namespace Tools
} // namespace Patterns

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_common_h


