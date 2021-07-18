//include/deal.II-translator/base/polynomials_barycentric_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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


#ifndef dealii_simplex_barycentric_polynomials_h
#define dealii_simplex_barycentric_polynomials_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/scalar_polynomials_base.h>
#include <deal.II/base/table.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 以arycentric坐标实现的多项式。
 * arycentric坐标是一个定义在单片机上的坐标系统，它特别容易操作，因为它把单片机中的坐标表达为顶点的凸组合。例如，一个三角形中的任何一点都可以写成
 * @f[
 * (x, y) = c_0 (x_0, y_0) + c_1 (x_1, y_1) + c_2 (x_2, y_2).
 * @f]
 * 其中每个值 $c_i$
 * 是每个顶点的相对权重（所以中心点在二维中是每个 $c_i
 * = 1/3$
 * ）。由于我们只考虑凸形组合，我们可以把这个方程改写为
 * @f[
 * (x, y) = (1
 *
 * - c_1
 *
 * - c_2) (x_0, y_0) + c_1 (x_1, y_1) + c_2 (x_2, y_2).
 * @f]
 * 这导致三个多项式在二维中等同于 $P^1$
 * 。更确切地说，这个类实现了一个以二维的基础定义的多项式空间，即
 *
 * @f{align*}{
 * t_0(x, y) &= 1
 *
 * - x
 *
 * - y \\
 * t_1(x, y) &= x \\
 * t_2(x, y) &= y
 * @f}
 * 而在三维中。
 *
 * @f{align*}{
 * t_0(x, y) &= 1
 *
 * - x
 *
 * - y
 *
 * - z \\
 * t_1(x, y) &= x             \\
 * t_2(x, y) &= y             \\
 * t_2(x, y) &= z
 * @f}
 *
 * 在实践中，这是一个非常方便的定义单项多项式的基础：例如，TRI6元素的第四个基础函数是
 * @f[
 * 4 t_1(x, y) t_2(x, y).
 * @f]
 * 在 <code>dim</code> -维空间的单项多项式有 <code>dim + 1</code>
 * variables in since <code>t_0</code> 可以用其他单项式来写。
 * 单项式可以方便地用 BarycentricPolynomial::monomial(). 构造。
 *
 *
 * @ingroup Polynomials
 *
 *
 */
template <int dim, typename Number = double>
class BarycentricPolynomial
{
public:
  /**
   * 零点多项式的构造函数。
   *
   */
  BarycentricPolynomial();

  /**
   * 单项式的构造函数。
   *
   */
  BarycentricPolynomial(const TableIndices<dim + 1> &powers,
                        const Number                 coefficient);

  /**
   * 返回指定的单项式。
   *
   */
  static BarycentricPolynomial<dim, Number>
  monomial(const unsigned int d);

  /**
   * 将多项式打印到输出流中，先打印最低阶的项。
   * 例如，第一个P6基函数被打印为 <code>-1 t0^1 + 2 t0^2</code>,
   * where <code>t0</code> 是第一个arycentric变量， <code>t1</code>
   * 是第二个，等等。
   *
   */
  void
  print(std::ostream &out) const;

  /**
   * 每个重心多项式的度数。
   *
   */
  TableIndices<dim + 1>
  degrees() const;

  /**
   * 单数减去。
   *
   */
  BarycentricPolynomial<dim, Number>
  operator-() const;

  /**
   * 增加一个标量。
   *
   */
  template <typename Number2>
  BarycentricPolynomial<dim, Number>
  operator+(const Number2 &a) const;

  /**
   * 减去一个标量。
   *
   */
  template <typename Number2>
  BarycentricPolynomial<dim, Number>
  operator-(const Number2 &a) const;

  /**
   * 乘以一个标量。
   *
   */
  template <typename Number2>
  BarycentricPolynomial<dim, Number> operator*(const Number2 &a) const;

  /**
   * 除以一个标量。
   *
   */
  template <typename Number2>
  BarycentricPolynomial<dim, Number>
  operator/(const Number2 &a) const;

  /**
   * 添加另一个巴里中心多项式。
   *
   */
  BarycentricPolynomial<dim, Number>
  operator+(const BarycentricPolynomial<dim, Number> &augend) const;

  /**
   * 减去另一个arycentric多项式。
   *
   */
  BarycentricPolynomial<dim, Number>
  operator-(const BarycentricPolynomial<dim, Number> &augend) const;

  /**
   * 乘以另一个arycentric多项式。
   *
   */
  BarycentricPolynomial<dim, Number>
  operator*(const BarycentricPolynomial<dim, Number> &multiplicand) const;

  /**
   * 在重心坐标中进行微分。
   *
   */
  BarycentricPolynomial<dim, Number>
  barycentric_derivative(const unsigned int coordinate) const;

  /**
   * 在直角坐标中进行微分。
   *
   */
  BarycentricPolynomial<dim, Number>
  derivative(const unsigned int coordinate) const;

  /**
   * 评估多项式。
   *
   */
  Number
  value(const Point<dim> &point) const;

  /**
   * 以字节为单位，返回该对象的内存使用量的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

protected:
  /**
   * 多项式的系数。指数是整数的索引。
   *
   */
  Table<dim + 1, Number> coefficients;

  /**
   * 用于巴里中心多项式的实用函数
   *
   * - 其方便之处在于以独立于维度的方式一次性循环所有的索引，但我们也需要访问底层表对象的实际索引。这个实用函数将一个积分索引转换为等价的TableIndices数组（也是隐含存储的多项式指数）。
   *
   */
  static TableIndices<dim + 1>
  index_to_indices(const std::size_t &          index,
                   const TableIndices<dim + 1> &extent);
};

/**
 * 基于arycentric多项式空间的标量多项式。
 *
 *
 */
template <int dim>
class BarycentricPolynomials : public ScalarPolynomialsBase<dim>
{
public:
  /**
   * 使得维度可以向外延伸。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 获取指定度数的标准拉格朗日基础。
   *
   */
  static BarycentricPolynomials<dim>
  get_fe_p_basis(const unsigned int degree);

  /**
   * 构造函数将多项式 @p degree 作为输入。
   *
   */
  BarycentricPolynomials(
    const std::vector<BarycentricPolynomial<dim>> &polynomials);

  /**
   * 访问操作符。
   *
   */
  const BarycentricPolynomial<dim> &operator[](const std::size_t i) const;

  /**
   * @copydoc   ScalarPolynomialsBase::evaluate() .
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
   * @copydoc   ScalarPolynomialsBase::compute_value()
   *
   */
  double
  compute_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_1st_derivative() 
   */
  Tensor<1, dim>
  compute_1st_derivative(const unsigned int i,
                         const Point<dim> & p) const override;


  /**
   *
   */
  Tensor<2, dim>
  compute_2nd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_3rd_derivative()
   *
   */
  Tensor<3, dim>
  compute_3rd_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_4th_derivative()
   * ScalarPolynomialsBase::compute_4th_derivative() .
   *
   */
  Tensor<4, dim>
  compute_4th_derivative(const unsigned int i,
                         const Point<dim> & p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_grad()
   *
   */
  Tensor<1, dim>
  compute_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::compute_grad_grad()
   *
   */
  Tensor<2, dim>
  compute_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * @copydoc   ScalarPolynomialsBase::memory_consumption()
   *
   */
  virtual std::size_t
  memory_consumption() const override;

  /**
   * @copydoc   ScalarPolynomialsBase::name()
   *
   */
  std::string
  name() const override;

  /**
   * @copydoc   ScalarPolynomialsBase::clone()
   *
   */
  virtual std::unique_ptr<ScalarPolynomialsBase<dim>>
  clone() const override;

protected:
  std::vector<BarycentricPolynomial<dim>> polys;

  Table<2, BarycentricPolynomial<dim>> poly_grads;

  Table<3, BarycentricPolynomial<dim>> poly_hessians;

  Table<4, BarycentricPolynomial<dim>> poly_third_derivatives;

  Table<5, BarycentricPolynomial<dim>> poly_fourth_derivatives;
};

// non-member template functions for algebra

/**
 * BarycentricPolynomial乘以一个常数。
 *
 *
 */
template <int dim, typename Number1, typename Number2>
BarycentricPolynomial<dim, Number1>
operator*(const Number2 &a, const BarycentricPolynomial<dim, Number1> &bp)
{
  return bp * Number1(a);
}

/**
 * 将一个常数添加到一个BarycentricPolynomial中。
 *
 *
 */
template <int dim, typename Number1, typename Number2>
BarycentricPolynomial<dim, Number1>
operator+(const Number2 &a, const BarycentricPolynomial<dim, Number1> &bp)
{
  return bp + Number1(a);
}

/**
 * 从一个常数中减去一个BarycentricPolynomial。
 *
 *
 */
template <int dim, typename Number1, typename Number2>
BarycentricPolynomial<dim, Number1>
operator-(const Number2 &a, const BarycentricPolynomial<dim, Number1> &bp)
{
  return bp - Number1(a);
}

/**
 * 将一个BarycentricPolynomial写到提供的输出流中。
 *
 *
 */
template <int dim, typename Number>
std::ostream &
operator<<(std::ostream &out, const BarycentricPolynomial<dim, Number> &bp)
{
  bp.print(out);
  return out;
}

// Template function definitions

// BarycentricPolynomial:
template <int dim, typename Number>
BarycentricPolynomial<dim, Number>::BarycentricPolynomial()
{
  TableIndices<dim + 1> extents;
  for (unsigned int d = 0; d < dim + 1; ++d)
    extents[d] = 1;
  coefficients.reinit(extents);

  coefficients(TableIndices<dim + 1>{}) = Number();
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number>::BarycentricPolynomial(
  const TableIndices<dim + 1> &powers,
  const Number                 coefficient)
{
  TableIndices<dim + 1> extents;
  for (unsigned int d = 0; d < dim + 1; ++d)
    extents[d] = powers[d] + 1;
  coefficients.reinit(extents);

  coefficients(powers) = coefficient;
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::monomial(const unsigned int d)
{
  AssertIndexRange(d, dim + 1);
  TableIndices<dim + 1> indices;
  indices[d] = 1;
  return BarycentricPolynomial<dim, Number>(indices, Number(1));
}



template <int dim, typename Number>
void
BarycentricPolynomial<dim, Number>::print(std::ostream &out) const
{
  const auto &coeffs     = this->coefficients;
  auto        first      = index_to_indices(0, coeffs.size());
  bool        print_plus = false;
  if (coeffs(first) != Number())
    {
      out << coeffs(first);
      print_plus = true;
    }
  for (std::size_t i = 1; i < coeffs.n_elements(); ++i)
    {
      const auto indices = index_to_indices(i, coeffs.size());
      if (coeffs(indices) == Number())
        continue;
      if (print_plus)
        out << " + ";
      out << coeffs(indices);
      for (unsigned int d = 0; d < dim + 1; ++d)
        {
          if (indices[d] != 0)
            out << " * t" << d << '^' << indices[d];
        }
      print_plus = true;
    }

  if (!print_plus)
    out << Number();
}



template <int dim, typename Number>
TableIndices<dim + 1>
BarycentricPolynomial<dim, Number>::degrees() const
{
  auto deg = coefficients.size();
  for (unsigned int d = 0; d < dim + 1; ++d)
    deg[d] -= 1;
  return deg;
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::operator-() const
{
  return *this * Number(-1);
}



template <int dim, typename Number>
template <typename Number2>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::operator+(const Number2 &a) const
{
  BarycentricPolynomial<dim, Number> result(*this);
  result.coefficients(index_to_indices(0, result.coefficients.size())) += a;

  return result;
}



template <int dim, typename Number>
template <typename Number2>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::operator-(const Number2 &a) const
{
  return *this + (-a);
}



template <int dim, typename Number>
template <typename Number2>
BarycentricPolynomial<dim, Number> BarycentricPolynomial<dim, Number>::
                                   operator*(const Number2 &a) const
{
  if (a == Number2())
    {
      return BarycentricPolynomial<dim, Number>();
    }

  BarycentricPolynomial<dim, Number> result(*this);
  for (std::size_t i = 0; i < result.coefficients.n_elements(); ++i)
    {
      const auto index = index_to_indices(i, result.coefficients.size());
      result.coefficients(index) *= a;
    }

  return result;
}



template <int dim, typename Number>
template <typename Number2>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::operator/(const Number2 &a) const
{
  Assert(a != Number2(), ExcDivideByZero());
  return *this * (Number(1) / Number(a));
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::
operator+(const BarycentricPolynomial<dim, Number> &augend) const
{
  TableIndices<dim + 1> deg;
  for (unsigned int d = 0; d < dim + 1; ++d)
    {
      deg[d] = std::max(degrees()[d], augend.degrees()[d]);
    }

  BarycentricPolynomial<dim, Number> result(deg, Number());

  auto add_coefficients = [&](const Table<dim + 1, Number> &in) {
    for (std::size_t i = 0; i < in.n_elements(); ++i)
      {
        const auto index = index_to_indices(i, in.size());
        result.coefficients(index) += in(index);
      }
  };

  add_coefficients(this->coefficients);
  add_coefficients(augend.coefficients);
  return result;
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::
operator-(const BarycentricPolynomial<dim, Number> &augend) const
{
  return *this + (-augend);
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number> BarycentricPolynomial<dim, Number>::
                                   operator*(const BarycentricPolynomial<dim, Number> &multiplicand) const
{
  TableIndices<dim + 1> deg;
  for (unsigned int d = 0; d < dim + 1; ++d)
    {
      deg[d] = multiplicand.degrees()[d] + degrees()[d];
    }

  BarycentricPolynomial<dim, Number> result(deg, Number());

  const auto &coef_1   = this->coefficients;
  const auto &coef_2   = multiplicand.coefficients;
  auto &      coef_out = result.coefficients;

  for (std::size_t i1 = 0; i1 < coef_1.n_elements(); ++i1)
    {
      const auto index_1 = index_to_indices(i1, coef_1.size());
      for (std::size_t i2 = 0; i2 < coef_2.n_elements(); ++i2)
        {
          const auto index_2 = index_to_indices(i2, coef_2.size());

          TableIndices<dim + 1> index_out;
          for (unsigned int d = 0; d < dim + 1; ++d)
            index_out[d] = index_1[d] + index_2[d];
          coef_out(index_out) += coef_1(index_1) * coef_2(index_2);
        }
    }

  return result;
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::barycentric_derivative(
  const unsigned int coordinate) const
{
  AssertIndexRange(coordinate, dim + 1);

  if (degrees()[coordinate] == 0)
    return BarycentricPolynomial<dim, Number>();

  auto deg = degrees();
  deg[coordinate] -= 1;
  BarycentricPolynomial<dim, Number> result(deg,
                                            std::numeric_limits<Number>::max());
  const auto &                       coeffs_in  = coefficients;
  auto &                             coeffs_out = result.coefficients;
  for (std::size_t i = 0; i < coeffs_out.n_elements(); ++i)
    {
      const auto out_index   = index_to_indices(i, coeffs_out.size());
      auto       input_index = out_index;
      input_index[coordinate] += 1;

      coeffs_out(out_index) = coeffs_in(input_index) * input_index[coordinate];
    }

  return result;
}



template <int dim, typename Number>
BarycentricPolynomial<dim, Number>
BarycentricPolynomial<dim, Number>::derivative(
  const unsigned int coordinate) const
{
  AssertIndexRange(coordinate, dim);
  return -barycentric_derivative(0) + barycentric_derivative(coordinate + 1);
}



template <int dim, typename Number>
Number
BarycentricPolynomial<dim, Number>::value(const Point<dim> &point) const
{
  // TODO: this is probably not numerically stable for higher order.
  // We really need some version of Horner's method.
  Number result = {};

  // Begin by converting point (which is in Cartesian coordinates) to
  // barycentric coordinates:
  std::array<Number, dim + 1> b_point;
  b_point[0] = 1.0;
  for (unsigned int d = 0; d < dim; ++d)
    {
      b_point[0] -= point[d];
      b_point[d + 1] = point[d];
    }

  // Now evaluate the polynomial at the computed barycentric point:
  for (std::size_t i = 0; i < coefficients.n_elements(); ++i)
    {
      const auto indices = index_to_indices(i, coefficients.size());
      const auto coef    = coefficients(indices);
      if (coef == Number())
        continue;

      auto temp = Number(1);
      for (unsigned int d = 0; d < dim + 1; ++d)
        temp *= std::pow(b_point[d], indices[d]);
      result += coef * temp;
    }

  return result;
}

template <int dim, typename Number>
std::size_t
BarycentricPolynomial<dim, Number>::memory_consumption() const
{
  return coefficients.memory_consumption();
}

template <int dim, typename Number>
TableIndices<dim + 1>
BarycentricPolynomial<dim, Number>::index_to_indices(
  const std::size_t &          index,
  const TableIndices<dim + 1> &extent)
{
  TableIndices<dim + 1> result;
  auto                  temp = index;

  for (unsigned int n = 0; n < dim + 1; ++n)
    {
      std::size_t slice_size = 1;
      for (unsigned int n2 = n + 1; n2 < dim + 1; ++n2)
        slice_size *= extent[n2];
      result[n] = temp / slice_size;
      temp %= slice_size;
    }
  return result;
}

template <int dim>
const BarycentricPolynomial<dim> &BarycentricPolynomials<dim>::
                                  operator[](const std::size_t i) const
{
  AssertIndexRange(i, polys.size());
  return polys[i];
}

DEAL_II_NAMESPACE_CLOSE

#endif


