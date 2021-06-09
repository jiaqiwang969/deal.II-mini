//include/deal.II-translator/base/polynomials_bdm_0.txt
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

#ifndef dealii_polynomials_BDM_h
#define dealii_polynomials_BDM_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/polynomial_space.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_polynomials_base.h>
#include <deal.II/base/thread_management.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 这个类实现了Brezzi和Fortin的<i>Mixed and Hybrid Finite Element
 * Methods</i>中描述的<i>H<sup>div</sup></i>--符合要求的、矢量值的Brezzi-Douglas-Marini（<i>
 * BDM </i>）多项式（参考119页
 *
 * - 124).
 * <i> BDM </i>多项式空间包含整个 $(P_{k})^{n}$
 * 空间（用PolynomialSpace Legendre多项式构建）以及
 * $(P_{k+1})^{n}$ 的一部分（即 $(P_{k})^{n} \subset BDM_{k} \subset
 * (P_{k+1})^{n}$  ）。 此外， $BDM_{k}$ 元素的设计使 $\nabla \cdot
 * q \in P_{k-1} (K)$ 和 $q \cdot n |_{e_{i}} \in P_{k}(e_{i})$
 * 。下面给出二维和三维 $BDM_{k}$ 元素的更多细节。 <dl>
 * <dt> 在二维：  <dd>   $ BDM_{k} = \{\mathbf{q} | \mathbf{q} = p_{k}
 * (x,y) + r \; \text{curl} (x^{k+1}y) + s \; \text{curl} (xy^{k+1}), p_{k}
 * \in (P_{k})^{2} \}$  。 注意：标量函数的卷曲是由
 * $\text{curl}(f(x,y)) = \begin{pmatrix} f_{y}(x,y) \\
 *
 * -f_{x}(x,y) \end{pmatrix}$ 给出的。 用来构建 $BDM_{1}$
 * 形状函数的基础是
 * @f{align*}{
 *   \phi_0 = \begin{pmatrix} 1 \\ 0 \end{pmatrix},
 *   \phi_1 = \begin{pmatrix}
 *
 * -\sqrt{3}+2\sqrt{3}x \\ 0 \end{pmatrix},
 *   \phi_2 = \begin{pmatrix}
 *
 * -\sqrt{3}+2\sqrt{3}y \\ 0 \end{pmatrix},
 *   \phi_3 = \begin{pmatrix} 0 \\ 1 \end{pmatrix},
 *   \phi_4 = \begin{pmatrix} 0 \\
 *
 * -\sqrt{3}+2\sqrt{3}x \end{pmatrix},
 *   \phi_5 = \begin{pmatrix} 0 \\
 *
 * -\sqrt{3}+2\sqrt{3}y \end{pmatrix},
 *   \phi_6 = \begin{pmatrix} x^2 \\
 *
 * -2xy \end{pmatrix},
 *   \phi_7 = \begin{pmatrix} 2xy \\
 *
 * -y^2 \end{pmatrix}.
 * @f}
 *
 * $BDM_{k}$ 空间的维度是 $(k+1)(k+2)+2$  ，每条边有 $k+1$
 * 个未知数，内部有 $k(k-1)$ 个未知数。 <dt> 在三维中：
 * <dd>   $ BDM_{k} = \{\mathbf{q} | \mathbf{q} = p_{k} (x,y,z) +
 * \sum_{i=0}^{k} ( r_{i} \; \text{curl} \begin{pmatrix} 0\\0\\xy^{i+1}z^{k-i}
 * \end{pmatrix} + s_{i} \; \text{curl} \begin{pmatrix} yz^{i+1}x^{k-i}\\0\\0
 * \end{pmatrix} + t_{i} \; \text{curl}
 * \begin{pmatrix}0\\zx^{i+1}y^{k-i}\\0\end{pmatrix}) , p_{k} \in (P_{k})^{3}
 * \}$  。 注意： $BDM_{k}$ 的三维描述不是唯一的。
 * 参见<i>Mixed and Hybrid Finite Element
 * Methods</i>第122页的替代定义。 $BDM_{k}$ 空间的尺寸是
 * $\dfrac{(k+1)(k+2)(k+3)}{2}+3(k+1)$  ，每个面有
 * $\dfrac{(k+1)(k+2)}{2}$ 个未知数， $\dfrac{(k-1)k(k+1)}{2}$
 * 个内部未知数。 </dl>
 *
 *
 *
 * @ingroup Polynomials
 *
 *
 */
template <int dim>
class PolynomialsBDM : public TensorPolynomialsBase<dim>
{
public:
  /**
   * 构造函数。创建给定度数的BDM多项式的所有基函数。
   * @arg
   * k：BDM-空间的度数，即BDM-空间中包含的最大完整多项式空间<i>P<sub>k</sub></i>的度数。
   *
   */
  PolynomialsBDM(const unsigned int k);

  /**
   * 计算每个BDM多项式在 @p unit_point. 的值和一、二次导数
   * 向量的大小必须是零或等于<tt>n()</tt>。
   * 在第一种情况下，该函数将不计算这些值。
   * 如果你需要所有张量积多项式的值或导数，那么使用这个函数，而不是使用任何<tt>compute_value</tt>,
   * <tt>compute_grad</tt>或<tt>compute_grad_grad</tt>函数，见下文，在所有张量积多项式上循环。
   *
   */
  void
  evaluate(const Point<dim> &           unit_point,
           std::vector<Tensor<1, dim>> &values,
           std::vector<Tensor<2, dim>> &grads,
           std::vector<Tensor<3, dim>> &grad_grads,
           std::vector<Tensor<4, dim>> &third_derivatives,
           std::vector<Tensor<5, dim>> &fourth_derivatives) const override;

  /**
   * 返回空间的名称，即<tt>BDM</tt>。
   *
   */
  std::string
  name() const override;

  /**
   * 返回空间<tt>BDM(degree)</tt>中的多项式的数目，而不需要建立PolynomialsBDM的对象。这是由FiniteElement类所要求的。
   *
   */
  static unsigned int
  n_polynomials(const unsigned int degree);

  /**
   * @copydoc   TensorPolynomialsBase::clone() .
   *
   */
  virtual std::unique_ptr<TensorPolynomialsBase<dim>>
  clone() const override;

private:
  /**
   * 一个代表这里使用的多项式空间的对象。构造函数用单项式基础来填充它。
   *
   */
  const PolynomialSpace<dim> polynomial_space;

  /**
   * 单项式的存储。在2D中，这只是阶数为<i>k</i>的多项式。在三维中，我们需要从0度到<i>k</i>的所有多项式。
   *
   */
  std::vector<Polynomials::Polynomial<double>> monomials;

  /**
   * 一个突扰器，可以守护以下的抓取数组。
   *
   */
  mutable Threads::Mutex mutex;

  /**
   * 辅助内存。
   *
   */
  mutable std::vector<double> p_values;

  /**
   * 辅助内存。
   *
   */
  mutable std::vector<Tensor<1, dim>> p_grads;

  /**
   * 辅助存储器。
   *
   */
  mutable std::vector<Tensor<2, dim>> p_grad_grads;

  /**
   * 辅助存储器。
   *
   */
  mutable std::vector<Tensor<3, dim>> p_third_derivatives;

  /**
   * 辅助存储器。
   *
   */
  mutable std::vector<Tensor<4, dim>> p_fourth_derivatives;
};


template <int dim>
inline std::string
PolynomialsBDM<dim>::name() const
{
  return "BDM";
}


DEAL_II_NAMESPACE_CLOSE

#endif


