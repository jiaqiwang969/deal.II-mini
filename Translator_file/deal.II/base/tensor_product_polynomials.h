//include/deal.II-translator/base/tensor_product_polynomials_0.txt
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

#ifndef dealii_tensor_product_polynomials_h
#define dealii_tensor_product_polynomials_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/ndarray.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/scalar_polynomials_base.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations for friends
// TODO: We may be able to modify these classes so they aren't
// required to be friends
template <int dim>
class TensorProductPolynomialsBubbles;
template <int dim>
class TensorProductPolynomialsConst;

/**
 * @addtogroup Polynomials   
     * @{ 
 *
 */

/**
 * 给定多项式的张量乘积。
 * 给定一个<i>n</i>一维多项式<i>P<sub>1</sub></i>到<i>P<sub>n</sub></i>的向量，这个类产生<i>n<sup>dim</sup></i>形式的多项式<i>Q<sub>ijk</sub>(x,y,z)
 * =
 * P<sub>i</sub>(x)P<sub>j</sub>(y)P<sub>k</sub>(z)</i>。如果基数多项式在区间[-1,1]或[0,1]上是相互正交的，那么张量积多项式在[-1,1]<sup>dim</sup>或[0,1]<sup>dim</sup>上分别是正交的。
 * 索引如下：dim-dimensional多项式的顺序是x-坐标跑得最快，然后是y-坐标，等等。因此，前几个多项式是<i>P<sub>1</sub>(x)P<sub>1</sub>(y),
 * P<sub>2</sub>(x)P<sub>1</sub>(y), P<sub>3</sub>(x)P<sub>1</sub>(y), ...,
 * P<sub>1</sub>(x)P<sub>2</sub>(y), P<sub>2</sub>(x)P<sub>2</sub>(y),
 * P<sub>3</sub>(x)P<sub>2</sub>(y),
 * ...</i>，同样，在三维中也是如此。
 * output_indices()函数打印出dim-dimensional多项式的排序，即对于多项式空间中的每个多项式，它给出了x、y和z方向的一维多项式的指数i,j,k。通过使用set_numbering()函数，可以改变二维多项式的排序。
 * @tparam  PolynomialType
 * 一个满足计算张量积所需接口的类。这个模板参数的典型选择是
 * Polynomials::Polynomial  和  Polynomials::PiecewisePolynomial.  。
 *
 */
template <int dim, typename PolynomialType = Polynomials::Polynomial<double>>
class TensorProductPolynomials : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * 访问此对象的维度，用于检查和自动设置其他类中的维度。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 构造函数。<tt>pols</tt>是一个对象的向量，应该是派生或以其他方式转换为`PolynomialType`类型（类的模板参数）的一维多项式对象。它将被逐个元素复制到一个受保护的成员变量中。
   *
   */
  template <class Pol>
  TensorProductPolynomials(const std::vector<Pol> &pols);

  /**
   * 打印索引的列表到<tt>out</tt>。
   *
   */
  void
  output_indices(std::ostream &out) const;

  /**
   * 设置多项式的排序。要求<tt>renumber.size()==n()</tt>。
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
   * 返回空间的名称，即<tt>TensorProductPolynomials</tt>。
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

  /**
   * 返回这个对象的内存消耗估计值（以字节为单位）。
   *
   */
  virtual std::size_t
  memory_consumption() const override;

  /**
   * 返回给这个类的构造函数的底层一维多项式的副本。
   *
   */
  std::vector<PolynomialType>
  get_underlying_polynomials() const;

protected:
  /**
   * 给予构造函数的多义词向量<tt>pols</tt>的拷贝。
   *
   */
  std::vector<PolynomialType> polynomials;

  /**
   * 用于重新排序多项式的索引图。
   *
   */
  std::vector<unsigned int> index_map;

  /**
   * 用于重新排序多项式的索引图。
   *
   */
  std::vector<unsigned int> index_map_inverse;

  /**
   * 每个张量积多项式<i>i</i>是每个空间方向上的一维多项式的乘积。给出指数<i>i</i>，计算每个空间方向的这些一维多项式的指数。
   *
   */
  void
  compute_index(const unsigned int             i,
                std::array<unsigned int, dim> &indices) const;

  /**
   * TensorProductPolynomialsBubbles有一个TensorProductPolynomials类，所以我们声明它是一个朋友类。
   *
   */
  friend class TensorProductPolynomialsBubbles<dim>;

  /**
   * TensorProductPolynomialsConst有一个TensorProductPolynomials类，所以我们声明它是一个朋友类。
   *
   */
  friend class TensorProductPolynomialsConst<dim>;
};



/**
 * 给定多项式的各向异性张量乘积。 给定一维多项式
 * $P^x_1(x), P^x_2(x), \ldots$  在  $x$  -方向，  $P^y_1(y), P^y_2(y),
 * \ldots$  在  $y$  -方向，以此类推，该类生成形式为
 * $Q_{ijk}(x,y,z) = P^x_i(x)P^y_j(y)P^z_k(z)$  的多项式。（如果  @p
 * dim  实际上只有2，则有明显的概括性。如果  @p dim
 * 实际上只有1，则结果只是传递给构造函数的同一组一维多项式）。
 * 如果每组基数多项式的元素在区间 $[-1,1]$ 或 $[0,1]$
 * 上相互正交，那么张量积多项式分别在 $[-1,1]^d$ 或
 * $[0,1]^d$ 上正交。 得到的 @p dim-dimensional
 * 张量乘积多项式的排序如下。我们在 $x$
 * 坐标上迭代运行最快，然后是 $y$
 * 坐标，等等。例如，对于 @p dim==2,
 * ，前几个多项式是这样的  $P^x_1(x)P^y_1(y)$  ,
 * $P^x_2(x)P^y_1(y)$  ,  $P^x_3(x)P^y_1(y)$  , ...,  $P^x_1(x)P^y_2(y)$  ,
 * $P^x_2(x)P^y_2(y)$  ,  $P^x_3(x)P^y_2(y)$  , 等等。
 *
 */
template <int dim>
class AnisotropicPolynomials : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * 构建器。  @p base_polynomials
   * 是一个一维多项式的表格。该表的行数（索引到 @p
   * base_polynomials)
   * 时的第一个索引需要等于空间维度，每一行的元素（即第二个索引）给出应在这个特定坐标方向使用的多项式。
   * 由于我们要建立<i>anisotropic</i>多项式，作为参数传入的
   * @p dim
   * 项的多项式集合当然可能不同，而且数量也可能不同。
   * 张量积多项式的数量是<tt>Nx*Ny*Nz</tt>，如果空间维数小于3，则去掉条款。
   *
   */
  AnisotropicPolynomials(
    const std::vector<std::vector<Polynomials::Polynomial<double>>>
      &base_polynomials);

  /**
   * 计算每个张量积多项式在<tt>unit_point</tt>的值和一、二次导数。
   * 向量的大小必须等于<tt>0</tt>或者等于<tt>this->n()</tt>。
   * 在第一种情况下，该函数将不计算这些值。
   * 如果你需要所有张量积多项式的值或导数，那么使用这个函数，而不是使用任何<tt>compute_value</tt>,
   * <tt>compute_grad</tt>或<tt>compute_grad_grad</tt>函数，见下文，在所有张量积多项式上循环。
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
   * 计算<tt>i</tt>第1个张量积多项式在<tt>unit_point</tt>的值。这里<tt>i</tt>是用张量积的编号给出的。
   * 注意，在所有张量积多项式的循环中使用这个函数并不高效，因为这样一来，底层（一维）多项式的每个点值都要（不必要地）计算多次。
   * 相反，使用<tt>compute</tt>函数，见上文，用<tt>values.size()==this->n()</tt>来一次性获得所有张量多项式的点值，而且效率更高。
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
   * 相反，使用<tt>compute</tt>函数，见上文，用<tt>grads.size()==this->n()</tt>来一次性获得所有张量多项式的点值，而且效率更高。
   *
   */
  Tensor<1, dim>
  compute_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 计算<tt>i</tt>第1个张量积多项式在<tt>unit_point</tt>的二阶导数（grad_grad）。这里<tt>i</tt>是用张量积的编号给出的。
   * 注意，在所有张量积多项式的循环中使用这个函数并不高效，因为这样一来，底层（一维）多项式的每个导数值都要（不必要地）计算多次。
   * 相反，使用<tt>compute</tt>函数，见上文，用<tt>grad_grads.size()==this->n()</tt>来一次性获得所有张量多项式的点值，而且效率更高。
   *
   */
  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 返回空间的名称，即<tt>AnisotropicPolynomials</tt>。
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
   * 给予构造函数的多义词向量<tt>pols</tt>的副本。
   *
   */
  const std::vector<std::vector<Polynomials::Polynomial<double>>> polynomials;

  /**
   * 每个张量积多项式 $p_i$
   * 是每个空间方向上的一维多项式的乘积。在给定指数<tt>i</tt>的情况下，计算每个空间方向的这些一维多项式的指数。
   *
   */
  void
  compute_index(const unsigned int             i,
                std::array<unsigned int, dim> &indices) const;

  /**
   * 给出构造函数的输入，计算<tt>n_pols</tt>。
   *
   */
  static unsigned int
  get_n_tensor_pols(
    const std::vector<std::vector<Polynomials::Polynomial<double>>> &pols);
};

 /** @} */ 

#ifndef DOXYGEN


 /* ---------------- template and inline functions ---------- */ 


template <int dim, typename PolynomialType>
template <class Pol>
inline TensorProductPolynomials<dim, PolynomialType>::TensorProductPolynomials(
  const std::vector<Pol> &pols)
  : ScalarPolynomialsBase<dim>(1, Utilities::fixed_power<dim>(pols.size()))
  , polynomials(pols.begin(), pols.end())
  , index_map(this->n())
  , index_map_inverse(this->n())
{
  // per default set this index map to identity. This map can be changed by
  // the user through the set_numbering() function
  for (unsigned int i = 0; i < this->n(); ++i)
    {
      index_map[i]         = i;
      index_map_inverse[i] = i;
    }
}


template <int dim, typename PolynomialType>
inline const std::vector<unsigned int> &
TensorProductPolynomials<dim, PolynomialType>::get_numbering() const
{
  return index_map;
}


template <int dim, typename PolynomialType>
inline const std::vector<unsigned int> &
TensorProductPolynomials<dim, PolynomialType>::get_numbering_inverse() const
{
  return index_map_inverse;
}


template <int dim, typename PolynomialType>
inline std::string
TensorProductPolynomials<dim, PolynomialType>::name() const
{
  return "TensorProductPolynomials";
}


template <int dim, typename PolynomialType>
template <int order>
Tensor<order, dim>
TensorProductPolynomials<dim, PolynomialType>::compute_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  std::array<unsigned int, dim> indices;
  compute_index(i, indices);

  ndarray<double, dim, 5> v;
  {
    std::vector<double> tmp(5);
    for (unsigned int d = 0; d < dim; ++d)
      {
        polynomials[indices[d]].value(p(d), tmp);
        v[d][0] = tmp[0];
        v[d][1] = tmp[1];
        v[d][2] = tmp[2];
        v[d][3] = tmp[3];
        v[d][4] = tmp[4];
      }
  }

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
              for (unsigned int x = 0; x < dim; ++x)
                {
                  unsigned int x_order = 0;
                  if (d == x)
                    ++x_order;

                  derivative_1[d] *= v[x][x_order];
                }
            }

          return derivative;
        }
      case 2:
        {
          Tensor<2, dim> &derivative_2 =
            *reinterpret_cast<Tensor<2, dim> *>(&derivative);
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            for (unsigned int d2 = 0; d2 < dim; ++d2)
              {
                derivative_2[d1][d2] = 1.;
                for (unsigned int x = 0; x < dim; ++x)
                  {
                    unsigned int x_order = 0;
                    if (d1 == x)
                      ++x_order;
                    if (d2 == x)
                      ++x_order;

                    derivative_2[d1][d2] *= v[x][x_order];
                  }
              }

          return derivative;
        }
      case 3:
        {
          Tensor<3, dim> &derivative_3 =
            *reinterpret_cast<Tensor<3, dim> *>(&derivative);
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            for (unsigned int d2 = 0; d2 < dim; ++d2)
              for (unsigned int d3 = 0; d3 < dim; ++d3)
                {
                  derivative_3[d1][d2][d3] = 1.;
                  for (unsigned int x = 0; x < dim; ++x)
                    {
                      unsigned int x_order = 0;
                      if (d1 == x)
                        ++x_order;
                      if (d2 == x)
                        ++x_order;
                      if (d3 == x)
                        ++x_order;

                      derivative_3[d1][d2][d3] *= v[x][x_order];
                    }
                }

          return derivative;
        }
      case 4:
        {
          Tensor<4, dim> &derivative_4 =
            *reinterpret_cast<Tensor<4, dim> *>(&derivative);
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            for (unsigned int d2 = 0; d2 < dim; ++d2)
              for (unsigned int d3 = 0; d3 < dim; ++d3)
                for (unsigned int d4 = 0; d4 < dim; ++d4)
                  {
                    derivative_4[d1][d2][d3][d4] = 1.;
                    for (unsigned int x = 0; x < dim; ++x)
                      {
                        unsigned int x_order = 0;
                        if (d1 == x)
                          ++x_order;
                        if (d2 == x)
                          ++x_order;
                        if (d3 == x)
                          ++x_order;
                        if (d4 == x)
                          ++x_order;

                        derivative_4[d1][d2][d3][d4] *= v[x][x_order];
                      }
                  }

          return derivative;
        }
      default:
        {
          Assert(false, ExcNotImplemented());
          return derivative;
        }
    }
}



template <>
template <int order>
Tensor<order, 0>
TensorProductPolynomials<0, Polynomials::Polynomial<double>>::
  compute_derivative(const unsigned int, const Point<0> &) const
{
  AssertThrow(false, ExcNotImplemented());

  return {};
}



template <int dim, typename PolynomialType>
inline Tensor<1, dim>
TensorProductPolynomials<dim, PolynomialType>::compute_1st_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<1>(i, p);
}



template <int dim, typename PolynomialType>
inline Tensor<2, dim>
TensorProductPolynomials<dim, PolynomialType>::compute_2nd_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<2>(i, p);
}



template <int dim, typename PolynomialType>
inline Tensor<3, dim>
TensorProductPolynomials<dim, PolynomialType>::compute_3rd_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<3>(i, p);
}



template <int dim, typename PolynomialType>
inline Tensor<4, dim>
TensorProductPolynomials<dim, PolynomialType>::compute_4th_derivative(
  const unsigned int i,
  const Point<dim> & p) const
{
  return compute_derivative<4>(i, p);
}



template <int dim>
template <int order>
Tensor<order, dim>
AnisotropicPolynomials<dim>::compute_derivative(const unsigned int i,
                                                const Point<dim> & p) const
{
  std::array<unsigned int, dim> indices;
  compute_index(i, indices);

  std::vector<std::vector<double>> v(dim, std::vector<double>(order + 1));
  for (unsigned int d = 0; d < dim; ++d)
    polynomials[d][indices[d]].value(p(d), v[d]);

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
              for (unsigned int x = 0; x < dim; ++x)
                {
                  unsigned int x_order = 0;
                  if (d == x)
                    ++x_order;

                  derivative_1[d] *= v[x][x_order];
                }
            }

          return derivative;
        }
      case 2:
        {
          Tensor<2, dim> &derivative_2 =
            *reinterpret_cast<Tensor<2, dim> *>(&derivative);
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            for (unsigned int d2 = 0; d2 < dim; ++d2)
              {
                derivative_2[d1][d2] = 1.;
                for (unsigned int x = 0; x < dim; ++x)
                  {
                    unsigned int x_order = 0;
                    if (d1 == x)
                      ++x_order;
                    if (d2 == x)
                      ++x_order;

                    derivative_2[d1][d2] *= v[x][x_order];
                  }
              }

          return derivative;
        }
      case 3:
        {
          Tensor<3, dim> &derivative_3 =
            *reinterpret_cast<Tensor<3, dim> *>(&derivative);
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            for (unsigned int d2 = 0; d2 < dim; ++d2)
              for (unsigned int d3 = 0; d3 < dim; ++d3)
                {
                  derivative_3[d1][d2][d3] = 1.;
                  for (unsigned int x = 0; x < dim; ++x)
                    {
                      unsigned int x_order = 0;
                      if (d1 == x)
                        ++x_order;
                      if (d2 == x)
                        ++x_order;
                      if (d3 == x)
                        ++x_order;

                      derivative_3[d1][d2][d3] *= v[x][x_order];
                    }
                }

          return derivative;
        }
      case 4:
        {
          Tensor<4, dim> &derivative_4 =
            *reinterpret_cast<Tensor<4, dim> *>(&derivative);
          for (unsigned int d1 = 0; d1 < dim; ++d1)
            for (unsigned int d2 = 0; d2 < dim; ++d2)
              for (unsigned int d3 = 0; d3 < dim; ++d3)
                for (unsigned int d4 = 0; d4 < dim; ++d4)
                  {
                    derivative_4[d1][d2][d3][d4] = 1.;
                    for (unsigned int x = 0; x < dim; ++x)
                      {
                        unsigned int x_order = 0;
                        if (d1 == x)
                          ++x_order;
                        if (d2 == x)
                          ++x_order;
                        if (d3 == x)
                          ++x_order;
                        if (d4 == x)
                          ++x_order;

                        derivative_4[d1][d2][d3][d4] *= v[x][x_order];
                      }
                  }

          return derivative;
        }
      default:
        {
          Assert(false, ExcNotImplemented());
          return derivative;
        }
    }
}



template <>
template <int order>
Tensor<order, 0>
AnisotropicPolynomials<0>::compute_derivative(const unsigned int,
                                              const Point<0> &) const
{
  AssertThrow(false, ExcNotImplemented());

  return {};
}



template <int dim>
inline Tensor<1, dim>
AnisotropicPolynomials<dim>::compute_1st_derivative(const unsigned int i,
                                                    const Point<dim> & p) const
{
  return compute_derivative<1>(i, p);
}



template <int dim>
inline Tensor<2, dim>
AnisotropicPolynomials<dim>::compute_2nd_derivative(const unsigned int i,
                                                    const Point<dim> & p) const
{
  return compute_derivative<2>(i, p);
}



template <int dim>
inline Tensor<3, dim>
AnisotropicPolynomials<dim>::compute_3rd_derivative(const unsigned int i,
                                                    const Point<dim> & p) const
{
  return compute_derivative<3>(i, p);
}



template <int dim>
inline Tensor<4, dim>
AnisotropicPolynomials<dim>::compute_4th_derivative(const unsigned int i,
                                                    const Point<dim> & p) const
{
  return compute_derivative<4>(i, p);
}



template <int dim>
inline std::string
AnisotropicPolynomials<dim>::name() const
{
  return "AnisotropicPolynomials";
}



#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif


