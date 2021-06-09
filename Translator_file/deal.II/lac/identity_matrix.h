//include/deal.II-translator/lac/identity_matrix_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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

#ifndef dealii_identity_matrix_h
#define dealii_identity_matrix_h


#include <deal.II/base/config.h>

#include <deal.II/lac/exceptions.h>

DEAL_II_NAMESPACE_OPEN

/*!   @addtogroup  Matrix1  
     * @{ 

 
*
*/


/**
 * 实现一个简单的类，代表给定大小的身份矩阵，即一个条目为
 * $A_{ij}=\delta_{ij}$
 * 的矩阵。虽然它具有矩阵的最重要的成分，特别是可以询问它的大小并对它进行矩阵-向量乘积，但这种类型的矩阵实际上只在两种情况下有用：预处理和初始化其他矩阵。
 * <h4>Initialization</h4>
 * 这个类的主要用处在于它能够初始化其他矩阵，就像这样。
 *
 * @code
 * FullMatrix<double> identity (IdentityMatrix(10));
 * @endcode
 *
 * 这将创建一个 $10\times 10$
 * 矩阵，对角线上是1，其他地方是0。大多数矩阵类型，特别是FullMatrix和SparseMatrix，都有IdentityMatrix的转换构造函数和赋值运算符，因此可以相当容易地用身份矩阵来填充。
 *
 *  <h4>Preconditioning</h4>
 * deal.II有一个专门的类用于此目的，即PreconditionIdentity，比可以在该类的文档中所示的背景下使用。本类可以用大致相同的方式使用，尽管没有任何额外的好处。
 *
 * @code
 * SolverControl           solver_control (1000, 1e-12);
 * SolverCG<>              cg (solver_control);
 * cg.solve (system_matrix, solution, system_rhs,
 *         IdentityMatrix(solution.size()));
 * @endcode
 *
 *
 *
 */
class IdentityMatrix
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 默认构造函数。创建一个零大小的矩阵，以后应该使用reinit()函数来调整其大小。
   *
   */
  IdentityMatrix();

  /**
   * 构造函数。创建一个大小为#n的身份矩阵。
   *
   */
  explicit IdentityMatrix(const size_type n);

  /**
   * 调整矩阵的大小为#n乘#n。
   *
   */
  void
  reinit(const size_type n);

  /**
   * 该矩阵的行数。对于本矩阵，行数和列数当然是相等的。
   *
   */
  size_type
  m() const;

  /**
   * 该矩阵的列数。对于本矩阵，行数和列数当然是相等的。
   *
   */
  size_type
  n() const;

  /**
   * 矩阵-向量的乘法。对于本例来说，这当然相当于简单地将输入向量复制到输出向量。
   *
   */
  template <typename OutVectorType, typename InVectorType>
  void
  vmult(OutVectorType &out, const InVectorType &in) const;

  /**
   * 矩阵向量乘以输出向量的加法。在本例中，这当然相当于简单地将输入向量加到输出向量上。
   *
   */
  template <typename OutVectorType, typename InVectorType>
  void
  vmult_add(OutVectorType &out, const InVectorType &in) const;

  /**
   * 矩阵-向量乘以转置矩阵。在本例中，这当然相当于简单地将输入向量复制到输出向量。
   *
   */
  template <typename OutVectorType, typename InVectorType>
  void
  Tvmult(OutVectorType &out, const InVectorType &in) const;


  /**
   * 与转置矩阵的矩阵向量乘法，并在输出向量上做加法。对于目前的情况，这当然相当于简单地将输入向量加到输出向量。
   *
   */
  template <typename OutVectorType, typename InVectorType>
  void
  Tvmult_add(OutVectorType &out, const InVectorType &in) const;

private:
  /**
   * 这个矩阵的行和列的数量。
   *
   */
  size_type size;
};



// ------------------------- inline and template functions -------------
#ifndef DOXYGEN


inline IdentityMatrix::IdentityMatrix()
  : size(0)
{}



inline IdentityMatrix::IdentityMatrix(const size_type n)
  : size(n)
{}



inline void
IdentityMatrix::reinit(const size_type n)
{
  size = n;
}



inline IdentityMatrix::size_type
IdentityMatrix::m() const
{
  return size;
}



inline IdentityMatrix::size_type
IdentityMatrix::n() const
{
  return size;
}



template <typename OutVectorType, typename InVectorType>
inline void
IdentityMatrix::vmult(OutVectorType &out, const InVectorType &in) const
{
  Assert(out.size() == size, ExcDimensionMismatch(out.size(), size));
  Assert(in.size() == size, ExcDimensionMismatch(in.size(), size));

  out = in;
}



template <typename OutVectorType, typename InVectorType>
inline void
IdentityMatrix::vmult_add(OutVectorType &out, const InVectorType &in) const
{
  Assert(out.size() == size, ExcDimensionMismatch(out.size(), size));
  Assert(in.size() == size, ExcDimensionMismatch(in.size(), size));

  out += in;
}



template <typename OutVectorType, typename InVectorType>
inline void
IdentityMatrix::Tvmult(OutVectorType &out, const InVectorType &in) const
{
  Assert(out.size() == size, ExcDimensionMismatch(out.size(), size));
  Assert(in.size() == size, ExcDimensionMismatch(in.size(), size));

  out = in;
}



template <typename OutVectorType, typename InVectorType>
inline void
IdentityMatrix::Tvmult_add(OutVectorType &out, const InVectorType &in) const
{
  Assert(out.size() == size, ExcDimensionMismatch(out.size(), size));
  Assert(in.size() == size, ExcDimensionMismatch(in.size(), size));

  out += in;
}


#endif

 /**@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


