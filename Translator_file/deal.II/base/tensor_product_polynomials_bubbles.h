//include/deal.II-translator/base/tensor_product_polynomials_bubbles_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2021 by the deal.II authors
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

#ifndef dealii_tensor_product_polynomials_bubbles_h
#define dealii_tensor_product_polynomials_bubbles_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/utilities.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Polynomials   
     * @{ 
 *
 */

/**
 * 一个代表张量积多项式空间的类，由 $dim$
 * （非归一化）形式的 $\varphi_j(\mathbf x) =
 * 2^{\text{degree}-1}\left(x_j-frac 12\right)^{\text{degree}-1}
 * \left[\prod_{i=0}^{dim-1}(x_i(1-x_i))\right]$ 气泡函数对
 * $j=0,\ldots,dim-1$ 进行增强。如果 "度
 * "为1，那么第一个因子就会消失，人们会得到以单元格中点为中心的通常的气泡函数。
 * 该类从TensorProductPolynomials继承了大部分功能。气泡的丰富性是为最后一个索引添加的。
 *
 *
 */
template <int dim>
class TensorProductPolynomialsBubbles : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * 访问此对象的维度，用于检查和自动设置其他类中的维度。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 构造函数。<tt>pols</tt>是一个对象的向量，应该是派生的或以其他方式转换为一维多项式对象。它将被逐个元素复制到一个私有变量中。
   *
   */
  template <class Pol>
  TensorProductPolynomialsBubbles(const std::vector<Pol> &pols);

  /**
   * 打印<tt>tensor_polys</tt>的索引列表到<tt>out</tt>。
   *
   */
  void
  output_indices(std::ostream &out) const;

  /**
   * 设置多项式的排序。需要<tt>renumber.size()==tensor_polys.n()</tt>。
   * 存储一个<tt>renumber</tt>的副本。
   *
   */
  void
  set_numbering(const std::vector<unsigned int> &renumber);

  /**
   * 给予对renumber向量的读取权限。
   *
   */
  const std::vector<unsigned int> &
  get_numbering() const;

  /**
   * 给予对逆向renumber向量的读取权限。
   *
   */
  const std::vector<unsigned int> &
  get_numbering_inverse() const;

  /**
   * 计算每个张量积多项式在<tt>unit_point</tt>的值和一、二导数。
   * 向量的大小必须等于0或等于n()。在第一种情况下，该函数将不计算这些值。
   * 如果你需要所有张量积多项式的值或导数，那么使用这个函数，而不是使用任何一个compute_value(),
   * compute_grad() 或 compute_grad_grad()
   * 函数，见下文，在所有张量积多项式上循环。
   *
   */
  void
  evaluate(const Point<dim> &           unit_point,
           std::vector<double> &        values,
           std::vector<Tensor<1, dim>> &grads,
           std::vector<Tensor<2, dim>> &grad_grads,
           std::vector<Tensor<3, dim>> &third_derivatives,
           std::vector<Tensor<4, dim>> &fourth_derivatives) const override;

  /**
   * 计算<tt>i</tt>第张量积多项式在<tt>unit_point</tt>的值。这里<tt>i</tt>是用张量积的编号给出的。
   * 注意，在所有张量积多项式的循环中使用这个函数并不高效，因为这样底层（一维）多项式的每个点值都要（不必要地）计算多次。
   * 相反，使用evaluate()函数和<tt>values.size()==</tt>n()来一次性获得所有张量多项式的点值，而且效率更高。
   *
   */
  double
  compute_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 计算<tt>i</tt>第1个张量积多项式在<tt>unit_point</tt>处的阶
   * @p order 导数。这里<tt>i</tt>是用张量积的编号给出的。
   * 注意，在所有张量积多项式的循环中使用这个函数并不高效，因为这样一来，底层（一维）多项式的每个导数值都要（不必要地）计算多次。
   * 相反，使用evaluate()函数，见上文，将适当的参数大小设置为n()，以一次获得所有张量多项式的点值，而且效率更高。
   *
   */
  template <int order>
  Tensor<order, dim>
  compute_derivative(const unsigned int i, const Point<dim> &p) const;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_1st_derivative() .
   *
   */
  virtual Tensor<1, dim>
  compute_1st_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_2nd_derivative()
   *
   */
  virtual Tensor<2, dim>
  compute_2nd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_3rd_derivative()
   *
   */
  virtual Tensor<3, dim>
  compute_3rd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_4th_derivative()
   *
   */
  virtual Tensor<4, dim>
  compute_4th_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * 计算<tt>i</tt>第张量积多项式在<tt>unit_point</tt>的梯度。这里<tt>i</tt>是用张量积的编号给出的。
   * 注意，在所有张量积多项式的循环中使用这个函数并不高效，因为这样一来，底层（一维）多项式的每个导数值都要（不必要地）计算多次。
   * 相反，使用evaluate()函数，见上文，用<tt>grads.size()==</tt>n()来一次性获得所有张量多项式的点值，而且效率更高。
   *
   */
  Tensor<1, dim>
  compute_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 计算<tt>i</tt>第1个张量积多项式在<tt>unit_point</tt>的二阶导数（grad_grad）。这里<tt>i</tt>是用张量积的编号给出的。
   * 注意，在所有张量积多项式的循环中使用这个函数并不高效，因为这样一来，底层（一维）多项式的每个导数值都要（不必要地）计算多次。
   * 相反，使用evaluate()函数，见上文，用<tt>grad_grads.size()==</tt>n()来一次性获得所有张量多项式的点值，而且效率更高。
   *
   */
  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 返回张量积多项式的数量，加上气泡的丰富程度。对于<i>n</i>1d多项式，如果多项式的最大度数是1，这就是<i>n<sup>dim</sup>+1</i>，否则就是<i>n<sup>dim</sup>+dim</i>。
   *
   */
  unsigned int
  n() const;

  /**
   * 返回空间的名称，即<tt>TensorProductPolynomialsBubbles</tt>。
   *
   */
  std::string
  name() const override;

  /**
   * @copydoc   ScalarPolynomialsBase::clone() .
   *
   */
  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override;

private:
  /**
   * TensorProductPolynomials对象
   *
   */
  TensorProductPolynomials<dim> tensor_polys;

  /**
   * 用于重新排序多项式的索引图。
   *
   */
  std::vector<unsigned int> index_map;

  /**
   * 用于对多项式进行重新排序的索引图。
   *
   */
  std::vector<unsigned int> index_map_inverse;
};

 /** @} */ 


 /* ---------------- template and inline functions ---------- */ 

#ifndef DOXYGEN

template <int dim>
template <class Pol>
inline TensorProductPolynomialsBubbles<dim>::TensorProductPolynomialsBubbles(
  const std::vector<Pol> &pols)
  : ScalarPolynomialsBase<dim>(1,
                               Utilities::fixed_power<dim>(pols.size()) + dim)
  , tensor_polys(pols)
  , index_map(tensor_polys.n() +
              ((tensor_polys.polynomials.size() <= 2) ? 1 : dim))
  , index_map_inverse(tensor_polys.n() +
                      ((tensor_polys.polynomials.size() <= 2) ? 1 : dim))
{
  const unsigned int q_degree  = tensor_polys.polynomials.size() - 1;
  const unsigned int n_bubbles = ((q_degree <= 1) ? 1 : dim);
  // append index for renumbering
  for (unsigned int i = 0; i < tensor_polys.n() + n_bubbles; ++i)
    {
      index_map[i]         = i;
      index_map_inverse[i] = i;
    }
}


template <int dim>
inline unsigned int
TensorProductPolynomialsBubbles<dim>::n() const
{
  return tensor_polys.n() + dim;
}


template <>
inline unsigned int
TensorProductPolynomialsBubbles<0>::n() const
{
  return numbers::invalid_unsigned_int;
}


template <int dim>
inline const std::vector<unsigned int> &
TensorProductPolynomialsBubbles<dim>::get_numbering() const
{
  return index_map;
}


template <int dim>
inline const std::vector<unsigned int> &
TensorProductPolynomialsBubbles<dim>::get_numbering_inverse() const
{
  return index_map_inverse;
}


template <int dim>
inline std::string
TensorProductPolynomialsBubbles<dim>::name() const
{
  return "TensorProductPolynomialsBubbles";
}


template <int dim>
template <int order>
Tensor<order, dim>
TensorProductPolynomialsBubbles<dim>::compute_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  const unsigned int q_degree      = tensor_polys.polynomials.size() - 1;
  const unsigned int max_q_indices = tensor_polys.n();
  Assert(i < max_q_indices +  /* n_bubbles= */  ((q_degree <= 1) ? 1 : dim),
         ExcInternalError());

  // treat the regular basis functions
  if (i < max_q_indices)
    return tensor_polys.template compute_derivative<order>(i, p);

  const unsigned int comp = i - tensor_polys.n();

  Tensor<order, dim> derivative;
  switch (order)
    {
      case 1:
        {
          Tensor<1, dim> &derivative_1 =
            *reinterpret_cast<Tensor<1, dim> *>(&derivative);

          for (unsigned int d = 0; d < dim; ++d)
            {
              derivative_1[d] = 1.;
              // compute grad(4*\prod_{i=1}^d (x_i(1-x_i)))(p)
              for (unsigned j = 0; j < dim; ++j)
                derivative_1[d] *=
                  (d == j ? 4 * (1 - 2 * p(j)) : 4 * p(j) * (1 - p(j)));
              // and multiply with (2*x_i-1)^{r-1}
              for (unsigned int i = 0; i < q_degree - 1; ++i)
                derivative_1[d] *= 2 * p(comp) - 1;
            }

          if (q_degree >= 2)
            {
              // add \prod_{i=1}^d 4*(x_i(1-x_i))(p)
              double value = 1.;
              for (unsigned int j = 0; j < dim; ++j)
                value *= 4 * p(j) * (1 - p(j));
              // and multiply with grad(2*x_i-1)^{r-1}
              double tmp = value * 2 * (q_degree - 1);
              for (unsigned int i = 0; i < q_degree - 2; ++i)
                tmp *= 2 * p(comp) - 1;
              derivative_1[comp] += tmp;
            }

          return derivative;
        }
      case 2:
        {
          Tensor<2, dim> &derivative_2 =
            *reinterpret_cast<Tensor<2, dim> *>(&derivative);

          double v[dim + 1][3];
          {
            for (unsigned int c = 0; c < dim; ++c)
              {
                v[c][0] = 4 * p(c) * (1 - p(c));
                v[c][1] = 4 * (1 - 2 * p(c));
                v[c][2] = -8;
              }

            double tmp = 1.;
            for (unsigned int i = 0; i < q_degree - 1; ++i)
              tmp *= 2 * p(comp) - 1;
            v[dim][0] = tmp;

            if (q_degree >= 2)
              {
                double tmp = 2 * (q_degree - 1);
                for (unsigned int i = 0; i < q_degree - 2; ++i)
                  tmp *= 2 * p(comp) - 1;
                v[dim][1] = tmp;
              }
            else
              v[dim][1] = 0.;

            if (q_degree >= 3)
              {
                double tmp = 4 * (q_degree - 2) * (q_degree - 1);
                for (unsigned int i = 0; i < q_degree - 3; ++i)
                  tmp *= 2 * p(comp) - 1;
                v[dim][2] = tmp;
              }
            else
              v[dim][2] = 0.;
          }

          // calculate (\partial_j \partial_k \psi) * monomial
          Tensor<2, dim> grad_grad_1;
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            for (unsigned int d2 = 0; d2 < dim; ++d2)
              {
                grad_grad_1[d1][d2] = v[dim][0];
                for (unsigned int x = 0; x < dim; ++x)
                  {
                    unsigned int derivative = 0;
                    if (d1 == x || d2 == x)
                      {
                        if (d1 == d2)
                          derivative = 2;
                        else
                          derivative = 1;
                      }
                    grad_grad_1[d1][d2] *= v[x][derivative];
                  }
              }

          // calculate (\partial_j  \psi) *(\partial_k monomial)
          // and (\partial_k  \psi) *(\partial_j monomial)
          Tensor<2, dim> grad_grad_2;
          Tensor<2, dim> grad_grad_3;
          for (unsigned int d = 0; d < dim; ++d)
            {
              grad_grad_2[d][comp] = v[dim][1];
              grad_grad_3[comp][d] = v[dim][1];
              for (unsigned int x = 0; x < dim; ++x)
                {
                  grad_grad_2[d][comp] *= v[x][d == x];
                  grad_grad_3[comp][d] *= v[x][d == x];
                }
            }

          // calculate \psi *(\partial j \partial_k monomial) and sum
          double psi_value = 1.;
          for (unsigned int x = 0; x < dim; ++x)
            psi_value *= v[x][0];

          for (unsigned int d1 = 0; d1 < dim; ++d1)
            for (unsigned int d2 = 0; d2 < dim; ++d2)
              derivative_2[d1][d2] =
                grad_grad_1[d1][d2] + grad_grad_2[d1][d2] + grad_grad_3[d1][d2];
          derivative_2[comp][comp] += psi_value * v[dim][2];

          return derivative;
        }
      default:
        {
          Assert(false, ExcNotImplemented());
          return derivative;
        }
    }
}



template <int dim>
inline Tensor<1, dim>
TensorProductPolynomialsBubbles<dim>::compute_1st_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<1>(i, p);
}



template <int dim>
inline Tensor<2, dim>
TensorProductPolynomialsBubbles<dim>::compute_2nd_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<2>(i, p);
}



template <int dim>
inline Tensor<3, dim>
TensorProductPolynomialsBubbles<dim>::compute_3rd_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<3>(i, p);
}



template <int dim>
inline Tensor<4, dim>
TensorProductPolynomialsBubbles<dim>::compute_4th_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<4>(i, p);
}

#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif


