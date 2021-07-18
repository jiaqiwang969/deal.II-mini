//include/deal.II-translator/lac/precondition_selector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#ifndef dealii_precondition_selector_h
#define dealii_precondition_selector_h


#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>

#include <string>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <class number>
class Vector;
template <class number>
class SparseMatrix;
#endif


/*!   @addtogroup Preconditioners  
     * @{ 

 
*
*/

/**
 * 选择预处理程序。这个类的构造函数接收预处理的名称和预处理的阻尼参数
 * @p omega ， @p use_matrix
 * 函数接收矩阵构建预处理函数所使用的矩阵。每次，<tt>operator()</tt>函数被调用时，这个预选的预处理程序、这个矩阵和这个
 * @p omega 被用于预处理。该类被设计为作为 @p 解算器的 @p
 * solve
 * 函数的参数，它涵盖了所有矩阵构建的预处理函数的选择。其他预处理函数的选择，如BlockSOR或ILU，应该由用户在派生类中处理。
 * <h3>Usage</h3> 这个类的最简单的用途是如下。
 *
 * @code
 * // generate a @p SolverControl and a @p VectorMemory
 * SolverControl control;
 * VectorMemory<Vector<double> > memory;
 *
 * // generate a solver
 * SolverCG<SparseMatrix<double>, Vector<double> > solver(control, memory);
 *
 * // generate a @p PreconditionSelector
 * PreconditionSelector<SparseMatrix<double>, Vector<double> >
 * preconditioning("jacobi", 1.);
 *
 * // give a matrix whose diagonal entries are to be used for the
 * // preconditioning. Generally the matrix of the linear equation system Ax=b.
 * preconditioning.use_matrix(A);
 *
 * // call the @p solve function with this preconditioning as last argument
 * solver.solve(A,x,b,preconditioning);
 * @endcode
 * 同样的例子，也使用了 @p SolverSelector 类，内容如下
 *
 * @code
 * // generate a @p SolverControl and a @p VectorMemory
 * SolverControl control;
 * VectorMemory<Vector<double> > memory;
 *
 * // generate a @p SolverSelector that calls the @p SolverCG
 * SolverSelector<SparseMatrix<double>, Vector<double> >
 * solver_selector("cg", control, memory);
 *
 * // generate a @p PreconditionSelector
 * PreconditionSelector<SparseMatrix<double>, Vector<double> >
 * preconditioning("jacobi", 1.);
 *
 * preconditioning.use_matrix(A);
 *
 * solver_selector.solve(A,x,b,preconditioning);
 * @endcode
 * 现在， @p SolverSelector 与 @p
 * PreconditionSelector的结合使用允许用户在程序开始时选择求解器和预处理器，每次求解器启动时（例如在非线性迭代中多次）都会调用预选的求解器和预处理器。
 *
 *
 */
template <typename MatrixType = SparseMatrix<double>,
          typename VectorType = dealii::Vector<double>>
class PreconditionSelector : public Subscriptor
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = typename MatrixType::size_type;

  /**
   * 构造函数。  @p omega 表示预调理的阻尼参数。
   *
   */
  PreconditionSelector(const std::string &                    preconditioning,
                       const typename VectorType::value_type &omega = 1.);

  /**
   * 解构器。
   *
   */
  virtual ~PreconditionSelector() override;

  /**
   * 取出涉及矩阵的预处理所需的矩阵。例如，对于 @p
   * precondition_jacobi,  <tt>~_sor</tt>, <tt>~_ssor</tt>。
   *
   */
  void
  use_matrix(const MatrixType &M);

  /**
   * 返回共域（或范围）空间的维数。注意，矩阵的维度是
   * $m \times n$  。
   *
   */
  size_type
  m() const;

  /**
   * 返回域空间的维度。请注意，矩阵的维度是 $m \times n$  .
   *
   */
  size_type
  n() const;

  /**
   * 预处理程序。调用构造函数中指定的预处理。
   *
   */
  virtual void
  vmult(VectorType &dst, const VectorType &src) const;

  /**
   * 转置前提条件程序。调用构造函数中指定的预设条件。
   *
   */
  virtual void
  Tvmult(VectorType &dst, const VectorType &src) const;

  /**
   * 获取所有实现的预设条件的名称。可能的选项列表包括。    <ul>   <li>  "无"  </li>   <li>  "jacobi"  </li>   <li>  "sor"  </li>   <li>  "ssor"  </li>   </ul>
   *
   */
  static std::string
  get_precondition_names();

  /**
   * @addtogroup  Exceptions  @{
   *
   */


  /**
   * 异常情况。
   *
   */
  DeclException0(ExcNoMatrixGivenToUse);

  //@}
protected:
  /**
   * 存储预处理的名称。
   *
   */
  std::string preconditioning;

private:
  /**
   * 用于矩阵构建预处理功能的矩阵。参见 @p
   * PreconditionUseMatrix. 。
   *
   */
  SmartPointer<const MatrixType, PreconditionSelector<MatrixType, VectorType>>
    A;

  /**
   * 存储预处理程序的阻尼参数。
   *
   */
  const typename VectorType::value_type omega;
};

 /*@}*/ 
 /* --------------------- Inline and template functions ------------------- */ 


template <typename MatrixType, typename VectorType>
PreconditionSelector<MatrixType, VectorType>::PreconditionSelector(
  const std::string &                    preconditioning,
  const typename VectorType::value_type &omega)
  : preconditioning(preconditioning)
  , omega(omega)
{}


template <typename MatrixType, typename VectorType>
PreconditionSelector<MatrixType, VectorType>::~PreconditionSelector()
{
  // release the matrix A
  A = nullptr;
}


template <typename MatrixType, typename VectorType>
void
PreconditionSelector<MatrixType, VectorType>::use_matrix(const MatrixType &M)
{
  A = &M;
}


template <typename MatrixType, typename VectorType>
inline typename PreconditionSelector<MatrixType, VectorType>::size_type
PreconditionSelector<MatrixType, VectorType>::m() const
{
  Assert(A != nullptr, ExcNoMatrixGivenToUse());
  return A->m();
}


template <typename MatrixType, typename VectorType>
inline typename PreconditionSelector<MatrixType, VectorType>::size_type
PreconditionSelector<MatrixType, VectorType>::n() const
{
  Assert(A != nullptr, ExcNoMatrixGivenToUse());
  return A->n();
}



template <typename MatrixType, typename VectorType>
void
PreconditionSelector<MatrixType, VectorType>::vmult(VectorType &      dst,
                                                    const VectorType &src) const
{
  if (preconditioning == "none")
    {
      dst = src;
    }
  else
    {
      Assert(A != nullptr, ExcNoMatrixGivenToUse());

      if (preconditioning == "jacobi")
        {
          A->precondition_Jacobi(dst, src, omega);
        }
      else if (preconditioning == "sor")
        {
          A->precondition_SOR(dst, src, omega);
        }
      else if (preconditioning == "ssor")
        {
          A->precondition_SSOR(dst, src, omega);
        }
      else
        Assert(false, ExcNotImplemented());
    }
}


template <typename MatrixType, typename VectorType>
void
PreconditionSelector<MatrixType, VectorType>::Tvmult(
  VectorType &      dst,
  const VectorType &src) const
{
  if (preconditioning == "none")
    {
      dst = src;
    }
  else
    {
      Assert(A != nullptr, ExcNoMatrixGivenToUse());

      if (preconditioning == "jacobi")
        {
          A->precondition_Jacobi(dst, src, omega); // Symmetric operation
        }
      else if (preconditioning == "sor")
        {
          A->precondition_TSOR(dst, src, omega);
        }
      else if (preconditioning == "ssor")
        {
          A->precondition_SSOR(dst, src, omega); // Symmetric operation
        }
      else
        Assert(false, ExcNotImplemented());
    }
}


template <typename MatrixType, typename VectorType>
std::string
PreconditionSelector<MatrixType, VectorType>::get_precondition_names()
{
  return "none|jacobi|sor|ssor";
}


DEAL_II_NAMESPACE_CLOSE

#endif


