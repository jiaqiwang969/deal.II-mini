//include/deal.II-translator/lac/exceptions_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2019 by the deal.II authors
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

#ifndef dealii_lac_exceptions_h
#define dealii_lac_exceptions_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

DEAL_II_NAMESPACE_OPEN

namespace LACExceptions
{
  /**
   * @addtogroup Exceptions
   *
   */
  //@{

  /**
   * 这个函数只对二次元矩阵有效。
   *
   */
  DeclExceptionMsg(ExcNotQuadratic,
                   "This function only works for quadratic objects!");

  /**
   * 由于矩阵是奇异的，所以操作无法完成。
   *
   */
  DeclException0(ExcSingular);

  /**
   * 两个块对象的块索引不同。
   *
   */
  DeclException0(ExcDifferentBlockIndices);

  /**
   * 该操作需要一个稀疏的模式。
   *
   */
  DeclExceptionMsg(
    ExcNeedsSparsityPattern,
    "This function requires that the current object have a "
    "sparsity pattern attached to it, but no sparsity pattern "
    "is available. This usually means that there is a missing "
    "reinit() call which would have added the sparsity pattern.");

  /**
   * 当一个PETSc函数报告错误时抛出的异常。如果可能的话，这个异常使用
   * <code>PetscErrorMessage</code> 提供的信息来打印错误的描述。
   * @note
   * 为了向后兼容，无论deal.II是否与PETSc一起编译，这都被定义。
   *
   */
  class ExcPETScError : public dealii::ExceptionBase
  {
  public:
    ExcPETScError(const int error_code);

    virtual void
    print_info(std::ostream &out) const override;

    const int error_code;
  };

  /**
   * 遇到了一个Trilinos函数的错误。请查阅Trilinos文档以了解详情。
   *
   */
  DeclException1(ExcTrilinosError,
                 int,
                 << "An error with error number " << arg1
                 << " occurred while calling a Trilinos function");

  //@}
} // namespace LACExceptions


using namespace LACExceptions;


DEAL_II_NAMESPACE_CLOSE

#endif


