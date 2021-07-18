//include/deal.II-translator/lac/sparse_decomposition_0.txt
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

#ifndef dealii_sparse_decomposition_h
#define dealii_sparse_decomposition_h

#include <deal.II/base/config.h>

#include <deal.II/lac/sparse_matrix.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

/*! @addtogroup Preconditioners
 *@{
 */

/**
 * 用于将稀疏矩阵不完全分解为稀疏因子的抽象基类。这个类本身不能使用，只能作为实际实现特定分解的派生类的基类，如SparseILU或SparseMIC。
 * 分解结果以稀疏矩阵的形式存储，这就是为什么这个类派生于SparseMatrix。因为它不是通常意义上的矩阵（存储的条目不是矩阵的条目，而是原始矩阵的两个因子），所以派生是<tt>保护的</tt>而不是<tt>公开的</tt>。
 *
 *  <h3>Fill-in</h3>
 * 稀疏分解经常被用于额外的填充，即分解的稀疏结构比要分解的矩阵更密集。该类的initialize()函数允许通过AdditionalData对象进行填充，只要原始矩阵中的所有条目在分解中也存在，即分解的稀疏模式是原始矩阵中稀疏模式的超集。
 * 这种填充可以通过各种方式完成，其中之一是SparsityPattern类的复制构造器，它允许在给定的稀疏结构中增加边对角线。
 *
 *  <h3>Unified use of preconditioners</h3>
 * 虽然这个类的对象不能直接使用（这个类只是其他实现实际分解的基类），但派生类如SparseILU和SparseMIC可以以通常的形式作为预处理器使用。例如，这样就可以了。
 *
 * @code
 * SparseILU<double> ilu;
 * ilu.initialize(matrix, SparseILU<double>::AdditionalData(...));
 *
 * somesolver.solve (A, x, f, ilu);
 * @endcode
 *
 * 通过AdditionalData对象，可以指定LU分解的额外参数。 1/
 * 矩阵的对角线可以通过添加 <code>strengthen_diagonal</code>
 * 倍的每一行的绝对行项之和来加强各自的对角线项。默认情况下不进行强化。
 * 2/
 * 默认情况下，每个initialize()函数调用都会创建自己的稀疏度。为此，它复制了
 * <code>matrix</code> 的稀疏性，并增加了
 * <code>extra_off_diagonals</code>
 * 所指定的特定数量的额外对角线条目。 3/ 通过设置
 * <code>use_previous_sparsity=true</code>
 * ，稀疏度不会被重新创建，但之前初始化()调用的稀疏度被重新使用（回收）。当需要解决几个相同稀疏度的线性问题时，这可能是有用的，例如，在同一个三角形上的几个牛顿迭代步骤。默认值是
 * <code>false</code>  。 4/ 可以给用户定义的稀疏度为
 * <code>use_this_sparsity</code>  。然后，不创建稀疏度，但
 * <code>*use_this_sparsity</code>
 * 被用来存储分解后的矩阵。关于稀疏度的限制见上面的
 * "填充 "部分）。)
 *
 *  <h3>Particular implementations</h3>
 * 覆盖initialize()和vmult()方法来实现特定的LU分解就足够了，比如真正的LU，或者Cholesky分解。此外，如果该分解需要在每一行的基础上微调对角线的强度，它可以覆盖get_strengthen_diagonal()方法。
 *
 *
 */
template <typename number>
class SparseLUDecomposition : protected SparseMatrix<number>,
                              public virtual Subscriptor
{
protected:
  /**
   * 构造函数。什么都不做。
   * 在使用此对象作为预处理（vmult()）之前，调用初始化()函数。
   *
   */
  SparseLUDecomposition();

public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = typename SparseMatrix<number>::size_type;

  /**
   * 销毁。将析构器标记为纯的，以确保这个类不被直接使用，而只是其派生类。
   *
   */
  virtual ~SparseLUDecomposition() override = 0;

  /**
   * 删除所有成员变量。将类留在调用构造函数后的直接状态中。
   *
   */
  virtual void
  clear() override;

  /**
   * SparseDecomposition的参数。
   *
   */
  class AdditionalData
  {
  public:
    /**
     * 构造函数。关于参数的描述，见下文。
     *
     */
    AdditionalData(const double           strengthen_diagonal   = 0,
                   const unsigned int     extra_off_diagonals   = 0,
                   const bool             use_previous_sparsity = false,
                   const SparsityPattern *use_this_sparsity     = nullptr);

    /**
     * <code>strengthen_diag</code>
     * 倍的绝对行条目之和被添加到对角线条目中。
     * 默认情况下，该值为零，即对角线不被加强。
     *
     */
    double strengthen_diagonal;

    /**
     * 默认情况下， <code>initialize(matrix, data)</code>
     * 函数创建自己的稀疏度。这种稀疏性具有与
     * <code>matrix</code>
     * 相同的SparsityPattern，有一些额外的对角线，其数量由
     * <code>extra_off_diagonals</code> 指定。        用户可以给
     * <code>use_this_sparsity</code> 一个SparsityPattern。
     * 然后使用这个稀疏度， <code>extra_off_diagonals</code>
     * 参数被忽略。
     *
     */
    unsigned int extra_off_diagonals;

    /**
     * 如果这个标志为真，initialize()函数就会使用与之前initialize()调用时相同的稀疏度。
     * 当需要解决几个相同稀疏度的线性问题时，这可能是有用的，例如，在同一个三角形上的几个牛顿迭代步骤。
     *
     */
    bool use_previous_sparsity;

    /**
     * 当一个SparsityPattern被赋予这个参数时，initialize()函数会调用
     * <code>reinit(*use_this_sparsity)</code>
     * 导致这个稀疏度被使用。        注意，
     * <code>*use_this_sparsity</code>
     * 的稀疏性结构和传递给initialize函数的矩阵不需要相等。
     * 允许填入，也允许过滤掉矩阵中的一些元素。
     *
     */
    const SparsityPattern *use_this_sparsity;
  };

  /**
   * 这个函数需要在这个类的对象被用作预处理器之前被调用。
   * 关于可能的参数的更多细节，请参阅类的文档和
   * SparseLUDecomposition::AdditionalData 类的文档。    根据
   * <code>parameters</code>
   * ，这个函数创建一个新的SparsityPattern，或者保持以前的稀疏度，或者采用用户给定的稀疏度
   * <code>data</code>  。然后，这个函数进行LU分解。
   * 这个函数被调用后，预处理程序就可以使用了（使用派生类的
   * <code>vmult</code> 函数）。
   *
   */
  template <typename somenumber>
  void
  initialize(const SparseMatrix<somenumber> &matrix,
             const AdditionalData            parameters);

  /**
   * 返回对象是否为空。它调用继承的 SparseMatrix::empty()
   * 函数。
   *
   */
  bool
  empty() const;

  /**
   * 返回共域（或范围）空间的维度。它调用继承的
   * SparseMatrix::m() 函数。注意，矩阵的维数是 $m \times n$  。
   *
   */
  size_type
  m() const;

  /**
   * 返回域空间的维度。它调用继承的 SparseMatrix::n()
   * 函数。注意，矩阵的维度是 $m \times n$  .
   *
   */
  size_type
  n() const;

  /**
   * 加法
   * 矩阵-向量乘法。在<i>dst</i>上添加<i>M*src</i>，<i>M</i>为该矩阵。
   * 源和目的不能是同一个向量。
   *
   */
  template <class OutVector, class InVector>
  void
  vmult_add(OutVector &dst, const InVector &src) const;

  /**
   * 添加矩阵-向量乘法。将<i>M<sup>T</sup>*src</i>加到<i>dst</i>，<i>M</i>是这个矩阵。这个函数的作用与vmult_add()相同，但取的是转置的矩阵。
   * 来源和目的地不能是同一个向量。
   *
   */
  template <class OutVector, class InVector>
  void
  Tvmult_add(OutVector &dst, const InVector &src) const;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  virtual std::size_t
  memory_consumption() const;

  /**
   * @addtogroup  Exceptions  
     * @{ 
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
  //@}
protected:
  /**
   * 将传递的SparseMatrix复制到这个对象上。这个对象的稀疏度模式保持不变。
   *
   */
  template <typename somenumber>
  void
  copy_from(const SparseMatrix<somenumber> &matrix);

  /**
   * 执行强化循环。对于每一行计算其元素的绝对值之和，确定加强因子（通过get_strengthen_diagonal()），并将对角线条目乘以
   * <code>sf+1</code>  。
   *
   */
  virtual void
  strengthen_diagonal_impl();

  /**
   * 在分解阶段，为第 <code>row</code>
   * 行的对角线条目计算一个加强系数，其元素的绝对值之和为
   * <code>rowsum</code>  。
   * @note  SparseLUDecomposition中的默认实现返回
   * <code>strengthen_diagonal</code>
   * 的值。这个变量在几个派生类中被设置为非零值。
   *
   */
  virtual number
  get_strengthen_diagonal(const number rowsum, const size_type row) const;

  /**
   * 默认的强化值，由get_strengthen_diagonal()返回。
   *
   */
  double strengthen_diagonal;

  /**
   * 对于底层SparsityPattern中的每一行，这个数组包含一个指向该行的第一个对角线条目的指针。在调用prebuild_lower_bound()后变得可用。
   *
   */
  std::vector<const size_type *> prebuilt_lower_bound;

  /**
   * 填充#prebuilt_lower_bound数组。
   *
   */
  void
  prebuild_lower_bound();

private:
  /**
   * 一般来说，这个指针是零，除了没有给这个类提供SparsityPattern的情况。然后，一个SparsityPattern被创建，并被传递给SparseMatrix基类。
   * 尽管如此，SparseLUDecomposition需要保留这个稀疏度的所有权。它保留这个指针，以便在销毁时删除这个稀疏度。
   *
   */
  SparsityPattern *own_sparsity;
};

 /*@}*/ 
//---------------------------------------------------------------------------

#ifndef DOXYGEN

template <typename number>
inline number
SparseLUDecomposition<number>::get_strengthen_diagonal(
  const number  /*rowsum*/ ,
  const size_type  /*row*/ ) const
{
  return strengthen_diagonal;
}



template <typename number>
inline bool
SparseLUDecomposition<number>::empty() const
{
  return SparseMatrix<number>::empty();
}


template <typename number>
inline typename SparseLUDecomposition<number>::size_type
SparseLUDecomposition<number>::m() const
{
  return SparseMatrix<number>::m();
}


template <typename number>
inline typename SparseLUDecomposition<number>::size_type
SparseLUDecomposition<number>::n() const
{
  return SparseMatrix<number>::n();
}

// Note: This function is required for full compatibility with
// the LinearOperator class. ::MatrixInterfaceWithVmultAdd
// picks up the vmult_add function in the protected SparseMatrix
// base class.
template <typename number>
template <class OutVector, class InVector>
inline void
SparseLUDecomposition<number>::vmult_add(OutVector &     dst,
                                         const InVector &src) const
{
  OutVector tmp;
  tmp.reinit(dst);
  this->vmult(tmp, src);
  dst += tmp;
}

// Note: This function is required for full compatibility with
// the LinearOperator class. ::MatrixInterfaceWithVmultAdd
// picks up the vmult_add function in the protected SparseMatrix
// base class.
template <typename number>
template <class OutVector, class InVector>
inline void
SparseLUDecomposition<number>::Tvmult_add(OutVector &     dst,
                                          const InVector &src) const
{
  OutVector tmp;
  tmp.reinit(dst);
  this->Tvmult(tmp, src);
  dst += tmp;
}

//---------------------------------------------------------------------------


template <typename number>
SparseLUDecomposition<number>::AdditionalData::AdditionalData(
  const double           strengthen_diag,
  const unsigned int     extra_off_diag,
  const bool             use_prev_sparsity,
  const SparsityPattern *use_this_spars)
  : strengthen_diagonal(strengthen_diag)
  , extra_off_diagonals(extra_off_diag)
  , use_previous_sparsity(use_prev_sparsity)
  , use_this_sparsity(use_this_spars)
{}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_sparse_decomposition_h


