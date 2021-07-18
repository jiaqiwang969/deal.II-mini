//include/deal.II-translator/base/derivative_form_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2020 by the deal.II authors
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

#ifndef dealii_derivative_form_h
#define dealii_derivative_form_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 这个类表示一个函数的（切向）导数  $ \mathbf F: {\mathbb
 * R}^{\text{dim}} \rightarrow {\mathbb R}^{\text{spacedim}}$
 * 。这样的函数总是用来将参考的dim-dimensional单元映射到spacedim-dimensional空间。对于这样的对象，函数的第一个导数是一个从
 * ${\mathbb R}^{\text{dim}}$ 到 ${\mathbb R}^{\text{spacedim}}$
 * 的线性映射，也就是说，它可以被表示为 ${\mathbb
 * R}^{\text{spacedim}\times \text{dim}}$
 * 的矩阵。这是有道理的，因为人们会用 $\mathbf x\in {\mathbb
 * R}^{\text{dim}}$ 来表示第一个导数， $\nabla \mathbf F(\mathbf x)$
 * ，这样，在方向 $\mathbf d\in {\mathbb R}^{\text{dim}}$
 * 上的定向导数，所以
 *
 * @f{align*}{
 * \nabla \mathbf F(\mathbf x) \mathbf d
 * = \lim_{\varepsilon\rightarrow 0}
 *   \frac{\mathbf F(\mathbf x + \varepsilon \mathbf d)
 *
 * - \mathbf F(\mathbf
 * x)}{\varepsilon},
 * @f}
 * 也就是说，人们需要能够用矩阵 $\nabla \mathbf F(\mathbf x)$
 * 乘以 ${\mathbb R}^{\text{dim}}$
 * 中的一个向量，其结果是函数值的差，这些函数值在
 * ${\mathbb R}^{\text{spacedim}}$ 。因此，矩阵的大小必须是
 * $\text{spacedim}\times\text{dim}$  。 同样，第二个导数是一个从
 * ${\mathbb R}^{\text{dim}} \times  {\mathbb R}^{\text{dim}}$ 到 ${\mathbb
 * R}^{\text{spacedim}}$
 * 的双线性映射，我们可以认为它是一个大小为
 * $\text{spacedim}\times\text{dim}\times\text{dim}$ 的等级3对象。
 * 在deal.II中，我们使用DerivativeForm  @<1,dim,spacedim,Number@>,
 * DerivativeForm  @<2,dim,spacedim,Number@>
 * 等类型的对象来表示这些导数。
 *
 *
 */
template <int order, int dim, int spacedim, typename Number = double>
class DerivativeForm
{
public:
  /**
   * 构造函数。将所有条目初始化为零。
   *
   */
  DerivativeForm() = default;

  /**
   * 来自张量的构造函数。
   *
   */
  DerivativeForm(const Tensor<order + 1, dim, Number> &);

  /**
   * 来自张量的构造函数。
   *
   */
  DerivativeForm(const Tensor<order, spacedim, Tensor<1, dim, Number>> &);

  /**
   * 读写访问操作符。
   *
   */
  Tensor<order, dim, Number> &operator[](const unsigned int i);

  /**
   * 只读访问操作符。
   *
   */
  const Tensor<order, dim, Number> &operator[](const unsigned int i) const;

  /**
   * 赋值运算符。
   *
   */
  DerivativeForm &
  operator=(const Tensor<order + 1, dim, Number> &);

  /**
   * 赋值运算符。
   *
   */
  DerivativeForm &
  operator=(const Tensor<order, spacedim, Tensor<1, dim, Number>> &);

  /**
   * 赋值运算符。
   *
   */
  DerivativeForm &
  operator=(const Tensor<1, dim, Number> &);

  /**
   * 将DerivativeForm <order, dim, dim, Number>转换为Tensor<order+1, dim,
   * Number>。特别是，如果order == 1并且导数是 $\mathbf F(\mathbf
   * x)$ 的Jacobian，那么Tensor[i] =  $\nabla F_i(\mathbf x)$  。
   *
   */
  operator Tensor<order + 1, dim, Number>() const;

  /**
   * 将DerivativeForm<1, dim, 1, Number>转换为Tensor<1, dim, Number>。
   *
   */
  operator Tensor<1, dim, Number>() const;

  /**
   * 返回一个矩形DerivativeForm的转置，被视为一个二维矩阵。
   *
   */
  DerivativeForm<1, spacedim, dim, Number>
  transpose() const;

  /**
   * 计算这个表格的Frobenius准则，即表达式 $\sqrt{\sum_{ij}
   * |DF_{ij}|^2} = \sqrt{\sum_{ij} |\frac{\partial F_i}{\partial x_j}|^2}$
   * 。
   *
   */
  typename numbers::NumberTraits<Number>::real_type
  norm() const;

  /**
   * 计算与变换 $\mathbf F$
   * 的贾可宾相关的体积元素。也就是说，如果 $DF$
   * 是正方形的，它计算 $\det(DF)$
   * ，如果DF不是正方形的，则返回 $\sqrt{\det(DF^T \,DF)}$  。
   *
   */
  Number
  determinant() const;

  /**
   * 假设当前对象存储了映射的雅各布系数  $\mathbf F$
   * ，那么当前函数计算导数的 <i>covariant</i> 形式，即
   * $(\nabla \mathbf F) {\mathbf G}^{-1}$  ，其中  $\mathbf G = (\nabla
   * \mathbf F)^{T}(\nabla \mathbf F)$  。如果 $\nabla \mathbf F$
   * 是一个方形矩阵（即 $\mathbf F: {\mathbb R}^n \mapsto {\mathbb
   * R}^n$ ），那么这个函数就简化为计算 $\nabla {\mathbf
   * F}^{-T}$  。
   *
   */
  DerivativeForm<1, dim, spacedim, Number>
  covariant_form() const;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  static std::size_t
  memory_consumption();

  /**
   * 异常情况。
   *
   */
  DeclException1(ExcInvalidTensorIndex,
                 int,
                 << "Invalid DerivativeForm index " << arg1);

private:
  /**
   * 计算 $A T^{T}$ 的辅助函数，其中A代表当前对象。
   *
   */
  DerivativeForm<1, dim, spacedim, Number>
  times_T_t(const Tensor<2, dim, Number> &T) const;


  /**
   * 持有子元素的张量数组。
   *
   */
  Tensor<order, dim, Number> tensor[spacedim];
};


 /*--------------------------- Inline functions -----------------------------*/ 

#ifndef DOXYGEN

template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<order, dim, spacedim, Number>::DerivativeForm(
  const Tensor<order + 1, dim, Number> &T)
{
  Assert((dim == spacedim),
         ExcMessage("Only allowed for forms with dim==spacedim."));
  if (dim == spacedim)
    for (unsigned int j = 0; j < dim; ++j)
      (*this)[j] = T[j];
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<order, dim, spacedim, Number>::DerivativeForm(
  const Tensor<order, spacedim, Tensor<1, dim, Number>> &T)
{
  for (unsigned int j = 0; j < spacedim; ++j)
    (*this)[j] = T[j];
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<order, dim, spacedim, Number> &
DerivativeForm<order, dim, spacedim, Number>::
operator=(const Tensor<order + 1, dim, Number> &ta)
{
  Assert((dim == spacedim), ExcMessage("Only allowed when dim==spacedim."));

  if (dim == spacedim)
    for (unsigned int j = 0; j < dim; ++j)
      (*this)[j] = ta[j];
  return *this;
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<order, dim, spacedim, Number> &
DerivativeForm<order, dim, spacedim, Number>::
operator=(const Tensor<order, spacedim, Tensor<1, dim, Number>> &T)
{
  for (unsigned int j = 0; j < spacedim; ++j)
    (*this)[j] = T[j];
  return *this;
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<order, dim, spacedim, Number> &
DerivativeForm<order, dim, spacedim, Number>::
operator=(const Tensor<1, dim, Number> &T)
{
  Assert((1 == spacedim) && (order == 1),
         ExcMessage("Only allowed for spacedim==1 and order==1."));

  (*this)[0] = T;

  return *this;
}



template <int order, int dim, int spacedim, typename Number>
inline Tensor<order, dim, Number> &
  DerivativeForm<order, dim, spacedim, Number>::operator[](const unsigned int i)
{
  AssertIndexRange(i, spacedim);

  return tensor[i];
}



template <int order, int dim, int spacedim, typename Number>
inline const Tensor<order, dim, Number> &
  DerivativeForm<order, dim, spacedim, Number>::
  operator[](const unsigned int i) const
{
  AssertIndexRange(i, spacedim);

  return tensor[i];
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<order, dim, spacedim, Number>::
operator Tensor<1, dim, Number>() const
{
  Assert((1 == spacedim) && (order == 1),
         ExcMessage("Only allowed for spacedim==1."));

  return (*this)[0];
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<order, dim, spacedim, Number>::
operator Tensor<order + 1, dim, Number>() const
{
  Assert((dim == spacedim), ExcMessage("Only allowed when dim==spacedim."));

  Tensor<order + 1, dim, Number> t;

  if (dim == spacedim)
    for (unsigned int j = 0; j < dim; ++j)
      t[j] = (*this)[j];

  return t;
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<1, spacedim, dim, Number>
DerivativeForm<order, dim, spacedim, Number>::transpose() const
{
  Assert(order == 1, ExcMessage("Only for rectangular DerivativeForm."));
  DerivativeForm<1, spacedim, dim, Number> tt;

  for (unsigned int i = 0; i < spacedim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      tt[j][i] = (*this)[i][j];

  return tt;
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<1, dim, spacedim, Number>
DerivativeForm<order, dim, spacedim, Number>::times_T_t(
  const Tensor<2, dim, Number> &T) const
{
  Assert(order == 1, ExcMessage("Only for order == 1."));
  DerivativeForm<1, dim, spacedim, Number> dest;
  for (unsigned int i = 0; i < spacedim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      dest[i][j] = (*this)[i] * T[j];

  return dest;
}



template <int order, int dim, int spacedim, typename Number>
inline typename numbers::NumberTraits<Number>::real_type
DerivativeForm<order, dim, spacedim, Number>::norm() const
{
  typename numbers::NumberTraits<Number>::real_type sum_of_squares = 0;
  for (unsigned int i = 0; i < spacedim; ++i)
    sum_of_squares += tensor[i].norm_square();
  return std::sqrt(sum_of_squares);
}



template <int order, int dim, int spacedim, typename Number>
inline Number
DerivativeForm<order, dim, spacedim, Number>::determinant() const
{
  Assert(order == 1, ExcMessage("Only for order == 1."));
  if (dim == spacedim)
    {
      const Tensor<2, dim, Number> T =
        static_cast<Tensor<2, dim, Number>>(*this);
      return dealii::determinant(T);
    }
  else
    {
      Assert(spacedim > dim, ExcMessage("Only for spacedim>dim."));
      const DerivativeForm<1, spacedim, dim, Number> DF_t = this->transpose();
      Tensor<2, dim, Number> G; // First fundamental form
      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
          G[i][j] = DF_t[i] * DF_t[j];

      return (std::sqrt(dealii::determinant(G)));
    }
}



template <int order, int dim, int spacedim, typename Number>
inline DerivativeForm<1, dim, spacedim, Number>
DerivativeForm<order, dim, spacedim, Number>::covariant_form() const
{
  if (dim == spacedim)
    {
      const Tensor<2, dim, Number> DF_t =
        dealii::transpose(invert(static_cast<Tensor<2, dim, Number>>(*this)));
      return DerivativeForm<1, dim, spacedim, Number>(DF_t);
    }
  else
    {
      const DerivativeForm<1, spacedim, dim, Number> DF_t = this->transpose();
      Tensor<2, dim, Number> G; // First fundamental form
      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
          G[i][j] = DF_t[i] * DF_t[j];

      return (this->times_T_t(invert(G)));
    }
}


template <int order, int dim, int spacedim, typename Number>
inline std::size_t
DerivativeForm<order, dim, spacedim, Number>::memory_consumption()
{
  return sizeof(DerivativeForm<order, dim, spacedim, Number>);
}

#endif // DOXYGEN



/**
 * DerivativeForm的一个用途是将其作为一个线性变换来应用。这个函数返回 $\nabla \mathbf F(\mathbf x) \Delta \mathbf x$ ，它近似于当 $\mathbf x$ 被改变的数量 $\Delta \mathbf x$  @f[
 * \nabla \mathbf F(\mathbf x) \; \Delta \mathbf x
 * \approx
 * \mathbf F(\mathbf x + \Delta \mathbf x)
 *
 * - \mathbf F(\mathbf x).
 * @f]时， $\mathbf F(\mathbf x)$ 的变化。这个变换对应于索引符号的@f[
 * [\text{result}]_{i_1,\dots,i_k} = i\sum_{j} \left[\nabla \mathbf F(\mathbf
 * x)\right]_{i_1,\dots,i_k, j} \Delta x_j @f]，对应于矩阵符号的
 * $[\Delta \mathbf x] [\nabla \mathbf F(\mathbf x)]^T$ 。
 * @relatesalso DerivativeForm
 *
 *
 */
template <int spacedim, int dim, typename Number1, typename Number2>
inline Tensor<1, spacedim, typename ProductType<Number1, Number2>::type>
apply_transformation(const DerivativeForm<1, dim, spacedim, Number1> &grad_F,
                     const Tensor<1, dim, Number2> &                  d_x)
{
  Tensor<1, spacedim, typename ProductType<Number1, Number2>::type> dest;
  for (unsigned int i = 0; i < spacedim; ++i)
    dest[i] = grad_F[i] * d_x;
  return dest;
}



/**
 * 与之前的apply_transformation()类似。结果的每一行都对应于由
 * @p grad_F, 转换的 @p D_X 中的一行，在矩阵符号中相当于
 * $\mathrm{D\_X} \, \mathrm{grad\_F}^T$ 。
 * @relatesalso DerivativeForm
 *
 *
 */
// rank=2
template <int spacedim, int dim, typename Number1, typename Number2>
inline DerivativeForm<1,
                      spacedim,
                      dim,
                      typename ProductType<Number1, Number2>::type>
apply_transformation(const DerivativeForm<1, dim, spacedim, Number1> &grad_F,
                     const Tensor<2, dim, Number2> &                  D_X)
{
  DerivativeForm<1, spacedim, dim, typename ProductType<Number1, Number2>::type>
    dest;
  for (unsigned int i = 0; i < dim; ++i)
    dest[i] = apply_transformation(grad_F, D_X[i]);

  return dest;
}



/**
 * 与之前的apply_transformation()类似。结果的每一行都对应于由
 * @p grad_F. 变换的 @p D_X 中的一行。
 * @relatesalso  DerivativeForm
 *
 */
template <int spacedim,
          int dim,
          int n_components,
          typename Number1,
          typename Number2>
inline Tensor<1,
              n_components,
              Tensor<1, spacedim, typename ProductType<Number1, Number2>::type>>
apply_transformation(
  const DerivativeForm<1, dim, spacedim, Number1> &       grad_F,
  const Tensor<1, n_components, Tensor<1, dim, Number2>> &D_X)
{
  Tensor<1,
         n_components,
         Tensor<1, spacedim, typename ProductType<Number1, Number2>::type>>
    dest;
  for (unsigned int i = 0; i < n_components; ++i)
    dest[i] = apply_transformation(grad_F, D_X[i]);

  return dest;
}



/**
 * 与之前的apply_transformation()类似。在矩阵符号中，它计算的是  $DF2 \, DF1^{T}$  。此外，这个操作的结果 $\mathbf A$ 可以解释为 ${\mathbb R}^\text{spacedim}$ 中的度量张量，对应于 ${\mathbb R}^\text{dim}$ 中的欧氏度量张量。对于每一对矢量 $\mathbf u, \mathbf v \in {\mathbb R}^\text{spacedim}$ ，我们有。@f[
 * \mathbf u \cdot \mathbf A \mathbf v =
 * \text{DF2}^{-1}(\mathbf u) \cdot \text{DF1}^{-1}(\mathbf v)
 * @f]
 * @relatesalso DerivativeForm
 *
 *
 */
template <int spacedim, int dim, typename Number1, typename Number2>
inline Tensor<2, spacedim, typename ProductType<Number1, Number2>::type>
apply_transformation(const DerivativeForm<1, dim, spacedim, Number1> &DF1,
                     const DerivativeForm<1, dim, spacedim, Number2> &DF2)
{
  Tensor<2, spacedim, typename ProductType<Number1, Number2>::type> dest;

  for (unsigned int i = 0; i < spacedim; ++i)
    dest[i] = apply_transformation(DF1, DF2[i]);

  return dest;
}



/**
 * 矩形DerivativeForm DF的转置，主要是为了兼容的原因。
 * @relatesalso  DerivativeForm
 *
 *
 */
template <int dim, int spacedim, typename Number>
inline DerivativeForm<1, spacedim, dim, Number>
transpose(const DerivativeForm<1, dim, spacedim, Number> &DF)
{
  DerivativeForm<1, spacedim, dim, Number> tt;
  tt = DF.transpose();
  return tt;
}


DEAL_II_NAMESPACE_CLOSE

#endif


