//include/deal.II-translator/base/tensor_product_polynomials_const_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2020 by the deal.II authors
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

#ifndef dealii_tensor_product_polynomials_const_h
#define dealii_tensor_product_polynomials_const_h


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
 * 给定的多项式和一个局部常数函数的张量乘积。这个类继承了TensorProductPolynomials的大部分功能。它的工作原理与该类类似，但为最后一个索引增加了一个常数函数。
 *
 *
 */
template <int dim>
class TensorProductPolynomialsConst : public ScalarPolynomialsBase<dim>
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
  TensorProductPolynomialsConst(const std::vector<Pol> &pols);

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
   * 计算<tt>i</tt>第张量积多项式在<tt>unit_point</tt>的<tt>阶</tt>次导数。这里<tt>i</tt>是用张量积的编号给出的。
   * 注意，在所有张量积多项式的循环中使用这个函数并不高效，因为这样一来，底层（一维）多项式的每个导数值都要（不必要地）计算多次。
   * 相反，使用evaluate()函数，见上文，将适当的参数大小设置为n()，可以一次性得到所有张量多项式的点值，而且效率更高。
   * @tparam 顺序 导数的顺序。
   *
   */
  template <int order>
  Tensor<order, dim>
  compute_derivative(const unsigned int i, const Point<dim> &p) const;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_1st_derivative()  阶数
   *
   */
  virtual Tensor<1, dim>
  compute_1st_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_2nd_derivative() 。
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
   * 返回张量积多项式的数量加上常数函数。对于<i>n</i>1d多项式，这就是<i>n<sup>dim</sup>+1</i>。
   *
   */
  unsigned int
  n() const;

  /**
   * 返回空间的名称，即<tt>TensorProductPolynomialsConst</tt>。
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
inline TensorProductPolynomialsConst<dim>::TensorProductPolynomialsConst(
  const std::vector<Pol> &pols)
  : ScalarPolynomialsBase<dim>(1, Utilities::fixed_power<dim>(pols.size()) + 1)
  , tensor_polys(pols)
  , index_map(tensor_polys.n() + 1)
  , index_map_inverse(tensor_polys.n() + 1)
{}



template <int dim>
inline unsigned int
TensorProductPolynomialsConst<dim>::n() const
{
  return tensor_polys.n() + 1;
}



template <int dim>
inline const std::vector<unsigned int> &
TensorProductPolynomialsConst<dim>::get_numbering() const
{
  return index_map;
}


template <int dim>
inline const std::vector<unsigned int> &
TensorProductPolynomialsConst<dim>::get_numbering_inverse() const
{
  return index_map_inverse;
}


template <int dim>
inline std::string
TensorProductPolynomialsConst<dim>::name() const
{
  return "TensorProductPolynomialsConst";
}


template <>
inline unsigned int
TensorProductPolynomialsConst<0>::n() const
{
  return numbers::invalid_unsigned_int;
}


template <int dim>
template <int order>
Tensor<order, dim>
TensorProductPolynomialsConst<dim>::compute_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  const unsigned int max_indices = tensor_polys.n();
  Assert(i <= max_indices, ExcInternalError());

  // treat the regular basis functions
  if (i < max_indices)
    return tensor_polys.template compute_derivative<order>(i, p);
  else
    // this is for the constant function
    return Tensor<order, dim>();
}



template <int dim>
inline Tensor<1, dim>
TensorProductPolynomialsConst<dim>::compute_1st_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<1>(i, p);
}



template <int dim>
inline Tensor<2, dim>
TensorProductPolynomialsConst<dim>::compute_2nd_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<2>(i, p);
}



template <int dim>
inline Tensor<3, dim>
TensorProductPolynomialsConst<dim>::compute_3rd_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<3>(i, p);
}



template <int dim>
inline Tensor<4, dim>
TensorProductPolynomialsConst<dim>::compute_4th_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<4>(i, p);
}

#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif


