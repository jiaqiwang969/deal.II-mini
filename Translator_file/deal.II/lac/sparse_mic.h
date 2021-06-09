//include/deal.II-translator/lac/sparse_mic_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2020 by the deal.II authors
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

#ifndef dealii_sparse_mic_h
#define dealii_sparse_mic_h

#include <deal.II/base/config.h>

#include <deal.II/lac/sparse_decomposition.h>
#include <deal.II/lac/sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN

/*!   @addtogroup Preconditioners  
     * @{ 

 
*
*/

/**
 * 实现对称矩阵的修正不完全Cholesky（MIC(0)）预处理。该类符合SparseLUDecomposition中的状态和使用规范。
 *
 *  <h3>The decomposition</h3>
 * 让一个对称的、正不定的、稀疏的矩阵 $A$ 的形式为 $A = D
 *
 *
 *
 * - L
 *
 * - L^T$ ，其中 $D$ 是 $A$ 的对角线部分， $-L$
 * 是一个严格的下三角形矩阵。矩阵 $A$ 的MIC(0)分解由 $B =
 * (X-L)X^{-1}(X-L^T)$ 定义，其中 $X$ 是一个由条件
 * $\text{rowsum}(A) = \text{rowsum}(B)$ 定义的对角线矩阵。
 *
 *
 */
template <typename number>
class SparseMIC : public SparseLUDecomposition<number>
{
public:
  /**
   * 申报容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 构造函数。什么都不做，所以你必须在事后有时调用 @p
   * decompose 。
   *
   */
  SparseMIC();

  /**
   * 解构器。
   *
   */
  virtual ~SparseMIC() override;

  /**
   * 删除所有成员变量。将类保留在调用构造函数后的直接状态。
   *
   */
  virtual void
  clear() override;

  /**
   * 使基类中的 @p AdditionalData 类型也能被这个类访问。
   *
   */
  using AdditionalData = typename SparseLUDecomposition<number>::AdditionalData;

  /**
   * 对给定的矩阵进行不完全的LU因子化。
   * 这个函数需要在这个类的对象被用作预调节器之前被调用。
   * 关于可能的参数的更多细节，请参阅SparseLUDecomposition的类文件和
   * @p  SparseLUDecomposition::AdditionalData 类的文件。    根据 @p
   * parameters,
   * ，该函数创建一个新的SparsityPattern，或保持以前的稀疏度，或采用用户给定的稀疏度，以
   * @p data. 然后，该函数执行MIC分解。
   * 在这个函数被调用后，预处理程序就可以使用了。
   *
   */
  template <typename somenumber>
  void
  initialize(const SparseMatrix<somenumber> &matrix,
             const AdditionalData &          parameters = AdditionalData());

  /**
   * 应用不完全分解，即做一个前向-后向步骤
   * $dst=(LU)^{-1}src$  。    在调用此函数之前，先调用 @p
   * initialize 。
   *
   */
  template <typename somenumber>
  void
  vmult(Vector<somenumber> &dst, const Vector<somenumber> &src) const;

  /**
   * 应用不完全分解的转置，即做一个向前-向后的步骤
   * $dst=(LU)^{-1}src$  。    在调用此函数之前调用 @p initialize
   * 。
   * @note  这个函数还没有被实现
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
  DeclException0(ExcStrengthenDiagonalTooSmall);
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
  DeclException2(ExcDecompositionNotStable,
                 int,
                 double,
                 << "The diagonal element (" << arg1 << "," << arg1 << ") is "
                 << arg2 << ", but must be positive");

  //@}
private:
  /**
   * 计算的对角线的值。
   *
   */
  std::vector<number> diag;

  /**
   * 对角线的倒数：预先计算以加快vmult。
   *
   */
  std::vector<number> inv_diag;

  /**
   * 计算的 "内部和
   * "的值，即位于对角线右侧的元素的每行之和。
   *
   */
  std::vector<number> inner_sums;

  /**
   * 计算每行的 "内部和"。
   *
   */
  number
  get_rowsum(const size_type row) const;
};

 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_


