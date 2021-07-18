//include/deal.II-translator/lac/lapack_support_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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

#ifndef dealii_lapack_support_h
#define dealii_lapack_support_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

DEAL_II_NAMESPACE_OPEN

namespace types
{
#ifdef LAPACK_WITH_64BIT_BLAS_INDICES
  /**
   * BLAS中的整数类型。
   *
   */
  using blas_int = long long;
#else
  /**
   * BLAS中的整数类型。
   *
   */
  using blas_int = int;
#endif
} // namespace types

/**
 * 一个包含常量、异常、枚举和其他由deal.II
 * LAPACK绑定使用的效用的命名空间。
 *
 *
 */
namespace LAPACKSupport
{
  /**
   * 大多数可以应用于矩阵的LAPACK函数（例如，通过调用这个类的成员函数）以某种方式改变其内容。例如，它们可以反转矩阵，或者用一个矩阵来代替它，该矩阵的列代表矩阵原始内容的特征向量。
   * 因此，这个枚举的元素被用来跟踪这个对象目前正在存储的内容。
   *
   */
  enum State
  {
    /// Contents is actually a matrix.
    matrix,
    /// Contents is the inverse of a matrix.
    inverse_matrix,
    /// Contents is an LU decomposition.
    lu,
    /// Contents is a Cholesky decomposition.
    cholesky,
    /// Eigenvalue vector is filled
    eigenvalues,
    /// Matrix contains singular value decomposition,
    svd,
    /// Matrix is the inverse of a singular value decomposition
    inverse_svd,
    /// Contents is something useless.
    unusable = 0x8000
  };

  /**
   * %函数打印一个国家的名称。
   *
   */
  inline const char *
  state_name(State s)
  {
    switch (s)
      {
        case matrix:
          return "matrix";
        case inverse_matrix:
          return "inverse matrix";
        case lu:
          return "lu decomposition";
        case cholesky:
          return "cholesky decomposition";
        case eigenvalues:
          return "eigenvalues";
        case svd:
          return "svd";
        case inverse_svd:
          return "inverse_svd";
        case unusable:
          return "unusable";
        default:
          return "unknown";
      }
  }

  /**
   * 一个矩阵可以有某些允许优化的特征，但很难测试。这些在此列出。
   *
   */
  enum Property
  {
    /// No special properties
    general = 0,
    /// Matrix is symmetric
    symmetric = 1,
    /// Matrix is upper triangular
    upper_triangular = 2,
    /// Matrix is lower triangular
    lower_triangular = 4,
    /// Matrix is diagonal
    diagonal = 6,
    /// Matrix is in upper Hessenberg form
    hessenberg = 8
  };

  /**
   * %函数打印一个属性的名称。
   *
   */
  inline const char *
  property_name(const Property s)
  {
    switch (s)
      {
        case general:
          return "general";
        case symmetric:
          return "symmetric";
        case upper_triangular:
          return "upper triangular";
        case lower_triangular:
          return "lower triangular";
        case diagonal:
          return "diagonal";
        case hessenberg:
          return "Hessenberg";
      }

    Assert(false, ExcNotImplemented());
    return "invalid";
  }

  /**
   * 字符常数。
   *
   */
  static const char A = 'A';
  /**
   * 字符常数。
   *
   */
  static const char N = 'N';
  /**
   * 字符常数。
   *
   */
  static const char O = 'O';
  /**
   * 字符常数。
   *
   */
  static const char T = 'T';
  /**
   * 字符常数。
   *
   */
  static const char U = 'U';
  /**
   * 字符常数。
   *
   */
  static const char L = 'L';
  /**
   * 字符常数。
   *
   */
  static const char V = 'V';
  /**
   * 整数常数。
   *
   */
  static const types::blas_int zero = 0;
  /**
   * 整数常数。
   *
   */
  static const types::blas_int one = 1;

  /**
   * 一个LAPACK函数返回一个错误代码。
   *
   */
  DeclException2(ExcErrorCode,
                 std::string,
                 types::blas_int,
                 << "The function " << arg1 << " returned with an error code "
                 << arg2);

  /**
   * 当一个矩阵不在适合操作的状态下时抛出的异常。例如，一个LAPACK程序可能使矩阵处于无法使用的状态，那么vmult就没有意义了。
   *
   */
  DeclException1(
    ExcState,
    State,
    << "The function cannot be called while the matrix is in state "
    << state_name(arg1));

  /**
   * 当矩阵没有合适的操作属性时抛出的异常。
   *
   */
  DeclException1(ExcProperty,
                 Property,
                 << "The function cannot be called with a "
                 << property_name(arg1) << " matrix.");

  /**
   * 如果某个LAPACK函数不可用，因为在配置过程中没有检测到LAPACK安装，就会抛出这个异常。
   *
   */
  DeclException1(
    ExcMissing,
    std::string,
    << "When you ran 'cmake' during installation of deal.II, "
    << "no suitable installation of the BLAS or LAPACK library could "
    << "be found. Consequently, the function <" << arg1
    << "> can not be called. Refer to the doc/readme.html "
    << "file for information on how to ensure that deal.II "
    << "picks up an existing BLAS and LAPACK installation at "
    << "configuration time.");
} // namespace LAPACKSupport


DEAL_II_NAMESPACE_CLOSE

#endif


