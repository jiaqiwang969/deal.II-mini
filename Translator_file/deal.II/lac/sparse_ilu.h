//include/deal.II-translator/lac/sparse_ilu_0.txt
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


#ifndef dealii_sparse_ilu_h
#define dealii_sparse_ilu_h


#include <deal.II/base/config.h>

#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/sparse_decomposition.h>
#include <deal.II/lac/sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN

/*!   @addtogroup Preconditioners  
     * @{ 

 
*
*/

/**
 * 该类计算稀疏矩阵的不完全LU（ILU）分解，使用相同的稀疏模式或不同的模式。我们所说的不完全是指，与精确的分解不同，不完全的分解也是使用稀疏因子计算的，分解中不适合给定的稀疏结构的条目被丢弃。
 * 本课使用的算法基本上是Y.Saad书中给出的算法的副本："稀疏线性系统的迭代方法"，第二版，在第10.3.2节。
 *
 *  <h3>Usage and state management</h3>
 * 请参考SparseLUDecomposition文档中的建议用法和状态管理。这个类在 @ref step_22  "  step-22  "
 * 教程中使用。
 *
 *
 * @note
 * 这个模板的实例化提供给<tt>  @<float@>  和  @<double@></tt>;
 * 其他可以在应用程序中生成（见手册中 @ref Instantiations
 * 部分）。
 *
 *
 */
template <typename number>
class SparseILU : public SparseLUDecomposition<number>
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = typename SparseLUDecomposition<number>::size_type;

  /**
   * 构造函数。什么都不做。
   * 在使用此对象作为预处理程序之前，请调用 @p initialize
   * 函数。
   *
   */
  SparseILU() = default;

  /**
   * 使 SparseLUDecomposition::AdditionalData 也能被这个类访问。
   *
   */
  using AdditionalData = typename SparseLUDecomposition<number>::AdditionalData;

  /**
   * 对给定的矩阵进行不完全LU因子化。
   * 这个函数需要在这个类的对象被用作预调节器之前被调用。
   * 关于可能的参数的更多细节，请参阅SparseLUDecomposition的类文件和
   * @p  SparseLUDecomposition::AdditionalData 类的文件。    根据 @p
   * parameters,
   * ，这个函数创建一个新的SparsityPattern，或者保持以前的稀疏度，或者采用用户给定的稀疏度，以
   * @p data. 然后，这个函数执行LU分解。
   * 在这个函数被调用后，预处理程序就可以使用了。
   *
   */
  template <typename somenumber>
  void
  initialize(const SparseMatrix<somenumber> &matrix,
             const AdditionalData &          parameters = AdditionalData());

  /**
   * 应用不完全分解，即做一个前向-后向步骤
   * $dst=(LU)^{-1}src$  。    初始化（）函数需要在之前调用。
   *
   */
  template <typename somenumber>
  void
  vmult(Vector<somenumber> &dst, const Vector<somenumber> &src) const;


  /**
   * 应用不完全分解的转置，即做一个向前向后的步骤
   * $dst=(LU)^{-T}src$  。    初始化（）函数需要在之前调用。
   *
   */
  template <typename somenumber>
  void
  Tvmult(Vector<somenumber> &dst, const Vector<somenumber> &src) const;


  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const override;

  /**
   * @addtogroup  Exceptions  @{
   *
   */

  /**
   * 异常情况
   *
   */
  DeclException1(ExcInvalidStrengthening,
                 double,
                 << "The strengthening parameter " << arg1
                 << " is not greater or equal than zero!");
  /**
   * 异常情况
   *
   */
  DeclException1(ExcZeroPivot,
                 size_type,
                 << "While computing the ILU decomposition, the algorithm "
                    "found a zero pivot on the diagonal of row "
                 << arg1
                 << ". This must stop the ILU algorithm because it means "
                    "that the matrix for which you try to compute a "
                    "decomposition is singular.");
  //@}
};

 /*@}*/ 
//---------------------------------------------------------------------------


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_sparse_ilu_h


