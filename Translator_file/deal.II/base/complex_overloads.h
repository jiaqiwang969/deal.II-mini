//include/deal.II-translator/base/complex_overloads_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2019 by the deal.II authors
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

#ifndef dealii_complex_overloads_h
#define dealii_complex_overloads_h

#include <deal.II/base/config.h>

#include <deal.II/base/template_constraints.h>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename T, typename U>
struct ProductType;
#endif

#ifndef DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS
/**
 * 提供一个<tt>operator*</tt>，对混合的复杂浮点类型进行操作。令人讨厌的是，标准库并没有提供这样的操作符...
 * @relatesalso  ProductType
 *
 *
 */
template <typename T, typename U>
inline typename std::enable_if<
  std::is_floating_point<T>::value && std::is_floating_point<U>::value,
  typename ProductType<std::complex<T>, std::complex<U>>::type>::type
operator*(const std::complex<T> &left, const std::complex<U> &right)
{
  using result_type =
    typename ProductType<std::complex<T>, std::complex<U>>::type;
  return static_cast<result_type>(left) * static_cast<result_type>(right);
}

/**
 * 提供一个<tt>operator/</tt>，对混合的复杂浮点类型进行操作。令人讨厌的是，标准库并没有提供这样的操作符...
 * @relatesalso  ProductType
 *
 *
 */
template <typename T, typename U>
inline typename std::enable_if<
  std::is_floating_point<T>::value && std::is_floating_point<U>::value,
  typename ProductType<std::complex<T>, std::complex<U>>::type>::type
operator/(const std::complex<T> &left, const std::complex<U> &right)
{
  using result_type =
    typename ProductType<std::complex<T>, std::complex<U>>::type;
  return static_cast<result_type>(left) / static_cast<result_type>(right);
}

/**
 * 提供一个<tt>运算符*</tt>，用于一个复数浮点类型与一个不同的实数浮点类型的标量乘法。令人讨厌的是，标准库并没有提供这样的操作符...
 * @relatesalso  ProductType
 *
 *
 */
template <typename T, typename U>
inline
  typename std::enable_if<std::is_floating_point<T>::value &&
                            std::is_floating_point<U>::value,
                          typename ProductType<std::complex<T>, U>::type>::type
  operator*(const std::complex<T> &left, const U &right)
{
  using result_type = typename ProductType<std::complex<T>, U>::type;
  return static_cast<result_type>(left) * static_cast<result_type>(right);
}

/**
 * 提供一个<tt>运算符/</tt>，用于复数浮点类型与不同实数浮点类型的标量除法。令人讨厌的是，标准库并没有提供这样的操作符...
 * @relatesalso  ProductType
 *
 *
 */
template <typename T, typename U>
inline
  typename std::enable_if<std::is_floating_point<T>::value &&
                            std::is_floating_point<U>::value,
                          typename ProductType<std::complex<T>, U>::type>::type
  operator/(const std::complex<T> &left, const U &right)
{
  using result_type = typename ProductType<std::complex<T>, U>::type;
  return static_cast<result_type>(left) / static_cast<result_type>(right);
}

/**
 * 提供一个<tt>运算符*</tt>，用于一个实数浮点类型与一个不同的复数浮点类型的标量乘法。令人讨厌的是，标准库并没有提供这样的操作符...
 * @relatesalso  ProductType
 *
 *
 */
template <typename T, typename U>
inline typename

  std::enable_if<std::is_floating_point<T>::value &&
                   std::is_floating_point<U>::value,
                 typename ProductType<T, std::complex<U>>::type>::type
  operator*(const T &left, const std::complex<U> &right)
{
  using result_type = typename ProductType<T, std::complex<U>>::type;
  return static_cast<result_type>(left) * static_cast<result_type>(right);
}

/**
 * 提供一个<tt>运算符/</tt>，用于实数浮点类型与不同复数浮点类型的标量除法。令人讨厌的是，标准库并没有提供这样的操作符...
 * @relatesalso  ProductType
 *
 *
 */
template <typename T, typename U>
inline
  typename std::enable_if<std::is_floating_point<T>::value &&
                            std::is_floating_point<U>::value,
                          typename ProductType<T, std::complex<U>>::type>::type
  operator/(const T &left, const std::complex<U> &right)
{
  using result_type = typename ProductType<T, std::complex<U>>::type;
  return static_cast<result_type>(left) / static_cast<result_type>(right);
}
#endif  /* DEAL_II_HAVE_COMPLEX_OPERATOR_OVERLOADS */ 

DEAL_II_NAMESPACE_CLOSE

#endif


