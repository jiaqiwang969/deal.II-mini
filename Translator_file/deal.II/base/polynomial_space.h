//include/deal.II-translator/base/polynomial_space_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2021 by the deal.II authors
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

#ifndef dealii_polynomial_space_h
#define dealii_polynomial_space_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/ndarray.h>
#include <deal.II/base/point.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/scalar_polynomials_base.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/tensor.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 度数最多为n的多项式的空间在高维中的表示。
 * 给定一个<i>n</i>一维多项式<i>P<sub>0</sub></i>到<i>P<sub>n</sub></i>的向量，其中<i>P<sub>i</sub></i>的度数为<i>i</i>，这个类产生所有形式为<i>
 * P<sub>ijk</sub>(x,y,z) =
 * P<sub>i</sub>(x)P<sub>j</sub>(y)P<sub>k</sub>(z)</i>的二维多项式，其中<i>i</i>、<i>j</i>和<i>k</i>之和小于或等于<i>n</i>。
 * output_indices()函数打印了多项式的排序，即对于多项式空间中的每个dim-dimensional多项式，它给出了x、y和z方向的一维多项式的指数i,j,k。dim-dimensional多项式的排序可以通过set_numbering()函数来改变。
 * 多项式的标准排序是，第一个空间维度的指数变化最快，最后一个空间维度的指数变化最慢。特别是，如果我们为简单起见，取单项式的矢量<i>x<sup>0</sup>,
 * x<sup>1</sup>, x<sup>2</sup>,..., x<sup>n</sup></i>，我们得到 <dl>
 * <dt> 1D  <dd>  <i> x<sup>0</sup>, x<sup>1</sup>,...,x<sup>n</sup></i> <dt>
 * 2D:  <dd>  <i> x<sup>0</sup>y<sup>0</sup>, x<sup>1</sup>y<sup>0</sup>,...,
 * x<sup>n</sup>y<sup>0</sup>,
 * <br>
 * x<sup>0</sup>y<sup>1</sup>, x<sup>1</sup>y<sup>1</sup>,...,
 * x<sup>n-1</sup>y<sup>1</sup>,
 * <br>
 * x<sup>0</sup>y<sup>2</sup>,... x<sup>n-2</sup>y<sup>2</sup>,
 * <br>
 * ...
 * <br>
 * x<sup>0</sup>y<sup>n-1</sup>, x<sup>1</sup>y<sup>n-1</sup>,
 * <br>
 * x<sup>0</sup>y<sup>n</sup> </i> <dt> 3D:  <dd>  <i>
 * x<sup>0</sup>y<sup>0</sup>z<sup>0</sup>,...,
 * x<sup>n</sup>y<sup>0</sup>z<sup>0</sup>,
 * <br>
 * x<sup>0</sup>y<sup>1</sup>z<sup>0</sup>,...,
 * x<sup>n-1</sup>y<sup>1</sup>z<sup>0</sup>,
 * <br>
 * ...
 * <br>
 * x<sup>0</sup>y<sup>n</sup>z<sup>0</sup>,
 * <br>
 * x<sup>0</sup>y<sup>0</sup>z<sup>1</sup>,...
 * x<sup>n-1</sup>y<sup>0</sup>z<sup>1</sup>,
 * <br>
 * ...
 * <br>
 * x<sup>0</sup>y<sup>n-1</sup>z<sup>1</sup>,
 * <br>
 * x<sup>0</sup>y<sup>0</sup>z<sup>2</sup>,...
 * x<sup>n-2</sup>y<sup>0</sup>z<sup>2</sup>,
 * <br>
 * ...
 * <br>
 * x<sup>0</sup>y<sup>0</sup>z<sup>n</sup> </i>  </dl>
 *
 *
 * @ingroup Polynomials
 *
 *
 */
template <int dim>
class PolynomialSpace : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * 访问此对象的尺寸，用于检查和自动设置其他类中的尺寸。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 构造函数。<tt>pols</tt>是一个指向一维多项式的向量，将被复制到一个私有成员变量中。模板参数<tt>pols</tt>的静态类型需要可以转换为
   * Polynomials::Polynomial@<double@>, ，即通常应该是
   * Polynomials::Polynomial@<double@>. 的一个派生类。
   *
   */
  template <class Pol>
  PolynomialSpace(const std::vector<Pol> &pols);

  /**
   * 打印索引列表到<tt>out</tt>。
   *
   */
  template <class StreamType>
  void
  output_indices(StreamType &out) const;

  /**
   * 设置多项式的排序。要求<tt>renumber.size()==n()</tt>。存储一个<tt>renumber</tt>的副本。
   *
   */
  void
  set_numbering(const std::vector<unsigned int> &renumber);

  /**
   * 计算每个多项式在<tt>unit_point</tt>的值和一、二次导数。
   * 向量的大小必须等于0或等于n()。在第一种情况下，函数不会计算这些值，也就是说，你要通过调整那些你想要填充的向量的大小来表明你想要计算什么。
   * 如果你需要所有多项式的值或导数，那么使用这个函数，而不是使用任何一个compute_value(),
   * compute_grad()或compute_grad_grad()函数，见下文，在所有多项式上循环。
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
   * 计算<tt>i</tt>第1个多项式在单位点<tt>p</tt>的值。
   * 可以考虑用evaluate()代替。
   *
   */
  double
  compute_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 计算<tt>i</tt>次多项式在单位点<tt>p</tt>的<tt>阶</tt>次导数。
   * 可以考虑用evaluate()代替。      @tparam  order 导数的阶数。
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
   * @copydoc   ScalarPolynomialsBase::compute_2nd_derivative()
   * ScalarPolynomialsBase::compute_2nd_derivative() 。
   *
   */
  virtual Tensor<2, dim>
  compute_2nd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_3rd_derivative()
   * ScalarPolynomialsBase::compute_3rd_derivative() .
   *
   */
  virtual Tensor<3, dim>
  compute_3rd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_4th_derivative()
   * ScalarPolynomialsBase::compute_4th_derivative() .
   *
   */
  virtual Tensor<4, dim>
  compute_4th_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * 计算<tt>i</tt>第1个多项式在单位点<tt>p</tt>的梯度。
   * 可以考虑用evaluate()代替。
   *
   */
  Tensor<1, dim>
  compute_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 计算<tt>i</tt>次多项式在单位点<tt>p</tt>的二阶导数（grad_grad）。
   * 可以考虑用evaluate()代替。
   *
   */
  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 返回跨越这个类所代表的空间的多项式的数量。这里，如果<tt>N</tt>是给出的一维多项式的数量，那么这个函数的结果是1d中的<i>N</i>，2d中的<i>N(N+1)/2</i>，以及3d中的<i>N(N+1)(N+2)/6</i>。
   *
   */
  static unsigned int
  n_polynomials(const unsigned int n);

  /**
   * 返回空间的名称，即<tt>PolynomialSpace</tt>。
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

protected:
  /**
   * 计算x、y和z方向的数字。给出d维多项式空间中的一个索引<tt>n</tt>，返回索引i,j,k，以便<i>p<sub>n</sub>(x,y,z)
   * = p<sub>i</sub>(x)p<sub>j</sub>(y)p<sub>k</sub>(z)</i>。
   * 在1d和2d中，显然只有i和i,j被返回。
   *
   */
  std::array<unsigned int, dim>
  compute_index(const unsigned int n) const;

private:
  /**
   * 复制给构造函数的多项式的向量<tt>pols</tt>。
   *
   */
  const std::vector<Polynomials::Polynomial<double>> polynomials;

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
};


 /* -------------- declaration of explicit specializations --- */ 

template <>
std::array<unsigned int, 1>
PolynomialSpace<1>::compute_index(const unsigned int n) const;
template <>
std::array<unsigned int, 2>
PolynomialSpace<2>::compute_index(const unsigned int n) const;
template <>
std::array<unsigned int, 3>
PolynomialSpace<3>::compute_index(const unsigned int n) const;



 /* -------------- inline and template functions ------------- */ 

template <int dim>
template <class Pol>
PolynomialSpace<dim>::PolynomialSpace(const std::vector<Pol> &pols)
  : ScalarPolynomialsBase<dim>(pols.size(), n_polynomials(pols.size()))
  , polynomials(pols.begin(), pols.end())
  , index_map(n_polynomials(pols.size()))
  , index_map_inverse(n_polynomials(pols.size()))
{
  // per default set this index map
  // to identity. This map can be
  // changed by the user through the
  // set_numbering function
  for (unsigned int i = 0; i < this->n(); ++i)
    {
      index_map[i]         = i;
      index_map_inverse[i] = i;
    }
}



template <int dim>
inline std::string
PolynomialSpace<dim>::name() const
{
  return "PolynomialSpace";
}


template <int dim>
template <class StreamType>
void
PolynomialSpace<dim>::output_indices(StreamType &out) const
{
  for (unsigned int i = 0; i < this->n(); ++i)
    {
      const std::array<unsigned int, dim> ix = compute_index(i);
      out << i << "\t";
      for (unsigned int d = 0; d < dim; ++d)
        out << ix[d] << " ";
      out << std::endl;
    }
}

template <int dim>
template <int order>
Tensor<order, dim>
PolynomialSpace<dim>::compute_derivative(const unsigned int i,
                                         const Point<dim> & p) const
{
  const std::array<unsigned int, dim> indices = compute_index(i);

  ndarray<double, dim, order + 1> v;
  {
    std::vector<double> tmp(order + 1);
    for (unsigned int d = 0; d < dim; ++d)
      {
        polynomials[indices[d]].value(p(d), tmp);
        for (unsigned int j = 0; j < order + 1; ++j)
          v[d][j] = tmp[j];
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



template <int dim>
inline Tensor<1, dim>
PolynomialSpace<dim>::compute_1st_derivative(const unsigned int i,
                                             const Point<dim> & p) const
{
  return compute_derivative<1>(i, p);
}



template <int dim>
inline Tensor<2, dim>
PolynomialSpace<dim>::compute_2nd_derivative(const unsigned int i,
                                             const Point<dim> & p) const
{
  return compute_derivative<2>(i, p);
}



template <int dim>
inline Tensor<3, dim>
PolynomialSpace<dim>::compute_3rd_derivative(const unsigned int i,
                                             const Point<dim> & p) const
{
  return compute_derivative<3>(i, p);
}



template <int dim>
inline Tensor<4, dim>
PolynomialSpace<dim>::compute_4th_derivative(const unsigned int i,
                                             const Point<dim> & p) const
{
  return compute_derivative<4>(i, p);
}

DEAL_II_NAMESPACE_CLOSE

#endif


