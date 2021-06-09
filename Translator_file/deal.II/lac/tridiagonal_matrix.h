//include/deal.II-translator/lac/tridiagonal_matrix_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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

#ifndef dealii_tridiagonal_matrix_h
#define dealii_tridiagonal_matrix_h

#include <deal.II/base/config.h>

#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/lapack_support.h>

#include <iomanip>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// forward declarations
#ifndef DOXYGEN
template <typename number>
class Vector;
#endif

/*!   @addtogroup  Matrix1  
     * @{ 

* 
*
*/


/**
 * 一个四次方三对角矩阵。也就是说，除了对角线和它的左右两边的条目外，所有条目都是零的矩阵。
 * 矩阵有一个额外的对称模式，在这种情况下，只有矩阵的上三角被存储并镜像到下三角，用于矩阵向量操作。
 *
 *
 * @ingroup Matrix1
 *
 *
 */
template <typename number>
class TridiagonalMatrix
{
public:
  ///@name Constructors
  //@{
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * @name  构造函数和初始化
   *
   */
  /**
   * 构造函数生成一个维度为<tt>n</tt>的空矩阵。
   *
   */
  TridiagonalMatrix(size_type n = 0, bool symmetric = false);

  /**
   * 重新初始化矩阵到一个新的尺寸，并将所有条目重置为零。对称性属性也可以设置。
   *
   */
  void
  reinit(size_type n, bool symmetric = false);


  //@}

  ///@name Non-modifying operators
  //@{

  /**
   * 这个矩阵的行数。注意，该矩阵是一个<i>m x m</i>矩阵。
   *
   */
  size_type
  m() const;

  /**
   * 该矩阵的列数。请注意，该矩阵是一个<i>n x n</i>矩阵。
   *
   */
  size_type
  n() const;

  /**
   * 返回该矩阵是否只包含数值为0的元素。这个函数主要用于内部一致性检查，在非调试模式下应该很少使用，因为它需要花费不少时间。
   *
   */
  bool
  all_zero() const;

  //@}

  ///@name Element access
  //@{
  /**
   * 对一个值的只读访问。这只限于<i>|i-j| <= 1</i>的情况。
   *
   */
  number
  operator()(size_type i, size_type j) const;

  /**
   * 对一个值的读写访问。这只限于<i>|i-j| <= 1</i>的情况。
   * @note
   * 在对称存储技术的情况下，条目<i>(i,j)</i>和<i>(j,i)</i>被识别，<b>both</b>存在。如果使用加法进行矩阵组合，必须考虑到这一点，以避免获得双倍的条目。
   *
   */
  number &
  operator()(size_type i, size_type j);

  //@}

  ///@name Multiplications with vectors
  //@{

  /**
   * 矩阵-向量-乘法。从右边乘以<tt>v</tt>，并将结果存入<tt>w</tt>。
   * 如果可选的参数<tt>adding</tt>是<tt>true</tt>，结果将被添加到<tt>w</tt>中。
   * 来源和目的地不能是同一个矢量。
   *
   */
  void
  vmult(Vector<number> &      w,
        const Vector<number> &v,
        const bool            adding = false) const;

  /**
   * 添加矩阵-向量-乘法。与vmult()相同，参数<tt>adding=true</tt>，但广泛用于<tt>deal.II</tt>类。
   * 源和目的不能是同一个向量。
   *
   */
  void
  vmult_add(Vector<number> &w, const Vector<number> &v) const;

  /**
   * 转置矩阵-向量-乘法。从左边乘以<tt>v<sup>T</sup></tt>，并将结果存入<tt>w</tt>。
   * 如果可选的参数<tt>adding</tt>是<tt>true</tt>，结果将被加到<tt>w</tt>中。
   * 来源和目的地不能是同一个矢量。
   *
   */
  void
  Tvmult(Vector<number> &      w,
         const Vector<number> &v,
         const bool            adding = false) const;

  /**
   * 添加转置的矩阵-向量-乘法。与Tvmult()相同，参数为<tt>adding=true</tt>，但广泛用于<tt>deal.II</tt>类。
   * 源和目的不能是同一个向量。
   *
   */
  void
  Tvmult_add(Vector<number> &w, const Vector<number> &v) const;

  /**
   * 建立矩阵标量乘积<tt>u^T M
   * v</tt>。这个函数在有限元背景下建立两个函数的单元标量乘积时大多有用。
   *
   */
  number
  matrix_scalar_product(const Vector<number> &u, const Vector<number> &v) const;

  /**
   * 返回向量<tt>v</tt>相对于该矩阵引起的规范的平方，即<i>(v,Mv)</i>。这很有用，例如在有限元背景下，一个函数的<i>L<sup>2</sup></i>规范等于相对于代表有限元函数节点值的向量的质量矩阵的矩阵规范。
   * 很明显，对于这个操作，矩阵需要是二次的。
   *
   */
  number
  matrix_norm_square(const Vector<number> &v) const;

  //@}

  ///@name LAPACK operations
  //@{
  /**
   * 计算对称三对角矩阵的特征值。
   * @note
   * 这个函数需要配置支持LAPACK的deal.II。此外，矩阵必须使用对称存储技术。
   *
   */
  void
  compute_eigenvalues();
  /**
   * 在调用compute_eigenvalues()后，你可以在这里访问每个特征值。
   *
   */
  number
  eigenvalue(const size_type i) const;
  //@}

  ///@name Miscellanea
  //@{
  /**
   * 以用户定义的格式输出矩阵。
   *
   */
  template <class OutputStream>
  void
  print(OutputStream &     s,
        const unsigned int width     = 5,
        const unsigned int precision = 2) const;
  //@}

private:
  /**
   * 对角线条目。
   *
   */
  std::vector<number> diagonal;
  /**
   * 对角线左边的条目。索引为0的条目总是0，因为第一行没有对角线左边的条目。因此，这个向量的长度与#对角线的长度相同。
   * 对于对称存储，这个向量的长度为零。在这种情况下，#left的第二个元素与#right的第一个元素是一致的。
   *
   */
  std::vector<number> left;
  /**
   * 对角线右边的条目。最后一个条目总是零，因为最后一行在对角线右边没有条目。因此，这个向量的长度和#对角线的长度是一样的。
   *
   */
  std::vector<number> right;

  /**
   * 如果这个标志为真，则只存储对角线右边的条目，并假定矩阵是对称的。
   *
   */
  bool is_symmetric;

  /**
   * 矩阵的状态。通常情况下，该状态为 LAPACKSupport::matrix,
   * ，表示该对象可以用于常规的矩阵操作。
   * 详见该数据类型的解释。
   *
   */
  LAPACKSupport::State state;
};

 /**@}*/ 

//---------------------------------------------------------------------------
#ifndef DOXYGEN

template <typename number>
types::global_dof_index
TridiagonalMatrix<number>::m() const
{
  return diagonal.size();
}



template <typename number>
types::global_dof_index
TridiagonalMatrix<number>::n() const
{
  return diagonal.size();
}


template <typename number>
inline number
TridiagonalMatrix<number>::operator()(size_type i, size_type j) const
{
  AssertIndexRange(i, n());
  AssertIndexRange(j, n());
  Assert(i <= j + 1, ExcIndexRange(i, j - 1, j + 2));
  Assert(j <= i + 1, ExcIndexRange(j, i - 1, i + 2));

  if (j == i)
    return diagonal[i];
  if (j == i - 1)
    {
      if (is_symmetric)
        return right[i - 1];
      else
        return left[i];
    }

  if (j == i + 1)
    return right[i];

  Assert(false, ExcInternalError());
  return 0;
}


template <typename number>
inline number &
TridiagonalMatrix<number>::operator()(size_type i, size_type j)
{
  AssertIndexRange(i, n());
  AssertIndexRange(j, n());
  Assert(i <= j + 1, ExcIndexRange(i, j - 1, j + 2));
  Assert(j <= i + 1, ExcIndexRange(j, i - 1, i + 2));

  if (j == i)
    return diagonal[i];
  if (j == i - 1)
    {
      if (is_symmetric)
        return right[i - 1];
      else
        return left[i];
    }

  if (j == i + 1)
    return right[i];

  Assert(false, ExcInternalError());
  return diagonal[0];
}


template <typename number>
template <class OutputStream>
void
TridiagonalMatrix<number>::print(OutputStream &     s,
                                 const unsigned int width,
                                 const unsigned int) const
{
  for (size_type i = 0; i < n(); ++i)
    {
      if (i > 0)
        s << std::setw(width) << (*this)(i, i - 1);
      else
        s << std::setw(width) << "";

      s << ' ' << (*this)(i, i) << ' ';

      if (i < n() - 1)
        s << std::setw(width) << (*this)(i, i + 1);

      s << std::endl;
    }
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


