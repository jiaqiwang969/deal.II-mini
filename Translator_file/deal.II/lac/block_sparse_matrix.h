//include/deal.II-translator/lac/block_sparse_matrix_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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

#ifndef dealii_block_sparse_matrix_h
#define dealii_block_sparse_matrix_h


#include <deal.II/base/config.h>

#include <deal.II/lac/block_matrix_base.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/sparse_matrix.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN


/*!   @addtogroup  Matrix1  
     * @{  

* 
*
*/


/**
 * 基于SparseMatrix类的阻塞式稀疏矩阵。这个类实现了SparseMatrix基对象的特定函数，用于阻塞式稀疏矩阵，并将实际工作中对各个块的大部分调用留给基类中实现的函数。关于这个类何时有用的描述，也请参见这里。
 * @see   @ref GlossBlockLA  "块（线性代数）"
 *
 *
 */
template <typename number>
class BlockSparseMatrix : public BlockMatrixBase<SparseMatrix<number>>
{
public:
  /**
   * 对基类进行类型化定义，以便更简单地访问它自己的别名。
   *
   */
  using BaseClass = BlockMatrixBase<SparseMatrix<number>>;

  /**
   * 对底层矩阵的类型进行类型化定义。
   *
   */
  using BlockType = typename BaseClass::BlockType;

  /**
   * 从基类中导入别名。
   *
   */
  using value_type      = typename BaseClass::value_type;
  using pointer         = typename BaseClass::pointer;
  using const_pointer   = typename BaseClass::const_pointer;
  using reference       = typename BaseClass::reference;
  using const_reference = typename BaseClass::const_reference;
  using size_type       = typename BaseClass::size_type;
  using iterator        = typename BaseClass::iterator;
  using const_iterator  = typename BaseClass::const_iterator;

  /**
   * @name  构造函数和初始化
   *
   */
  //@{
  /**
   * 构造函数；将矩阵初始化为空，没有任何结构，也就是说，矩阵根本无法使用。因此，这个构造函数只对作为类的成员的矩阵有用。所有其他的矩阵都应该在数据流中的一个点上创建，在那里所有必要的信息都是可用的。
   * 你必须在使用前用reinit(BlockSparsityPattern)初始化矩阵。然后每行和每列的块数由该函数决定。
   *
   */
  BlockSparseMatrix() = default;

  /**
   * 构造函数。使用给定的矩阵稀疏度结构来表示该矩阵的稀疏度模式。你可以在以后通过调用reinit()函数来改变稀疏性模式。
   * 这个构造函数用参数中的子疏密模式初始化所有子矩阵。
   * 你必须确保稀疏结构的寿命至少和这个矩阵的寿命一样长，或者只要reinit()没有被调用，就会有新的稀疏结构。
   *
   */
  BlockSparseMatrix(const BlockSparsityPattern &sparsity);

  /**
   * 解构器。
   *
   */
  virtual ~BlockSparseMatrix() override;



  /**
   * 伪拷贝操作符只拷贝空对象。块矩阵的大小需要相同。
   *
   */
  BlockSparseMatrix &
  operator=(const BlockSparseMatrix &);

  /**
   * 这个操作符将一个标量分配给一个矩阵。因为这通常没有什么意义（我们应该把所有的矩阵条目都设置为这个值吗？仅仅是稀疏模式的非零条目？），这个操作只允许在实际要分配的值为零的情况下进行。这个操作符的存在只是为了允许明显的符号<tt>matrix=0</tt>，它将矩阵的所有元素设置为零，但保留之前使用的稀疏模式。
   *
   */
  BlockSparseMatrix &
  operator=(const double d);

  /**
   * 释放所有内存并返回到与调用默认构造函数后相同的状态。它也忘记了之前所绑定的稀疏模式。
   * 这在所有子矩阵上调用 SparseMatrix::clear
   * ，然后将此对象重置为完全没有块。
   *
   */
  void
  clear();

  /**
   * 用给定的稀疏模式重新初始化稀疏矩阵。后者告诉矩阵需要保留多少个非零元素。
   * 基本上，这个函数只调用 SparseMatrix::reinit()
   * 的子矩阵与参数的块状稀疏模式。
   * 你必须确保稀疏结构的寿命至少和这个矩阵的寿命一样长，或者只要reinit(const
   * SparsityPattern &)没有被调用新的稀疏结构。
   * 矩阵的元素被这个函数设置为零。
   *
   */
  virtual void
  reinit(const BlockSparsityPattern &sparsity);
  //@}

  /**
   * @name  矩阵的信息
   *
   */
  //@{
  /**
   * 返回该对象是否为空。如果两个维度都是零或者没有关联BlockSparsityPattern，它就是空的。
   *
   */
  bool
  empty() const;

  /**
   * 返回特定行中的条目数。
   *
   */
  size_type
  get_row_length(const size_type row) const;

  /**
   * 返回这个矩阵的非零元素的数量。实际上，它返回的是稀疏模式中的条目数；如果任何一个条目碰巧是零，无论如何都会被计算在内。
   *
   */
  size_type
  n_nonzero_elements() const;

  /**
   * 返回实际非零元素的数量。只是计算所有块中实际非零元素的数量（绝对值大于阈值）。
   *
   */
  size_type
  n_actually_nonzero_elements(const double threshold = 0.0) const;

  /**
   * 返回一个对该矩阵的底层稀疏模式的（常数）引用。
   * 尽管返回值被声明为<tt>const</tt>，但你应该注意，如果你调用任何对其进行操作的对象的非常量函数，它可能会发生变化。
   *
   */
  const BlockSparsityPattern &
  get_sparsity_pattern() const;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;
  //@}

  /**
   * @name  乘法运算
   *
   */
  //@{
  /**
   * 矩阵-向量乘法：让 $dst = M*src$ 与 $M$ 为这个矩阵。
   *
   */
  template <typename block_number>
  void
  vmult(BlockVector<block_number> &      dst,
        const BlockVector<block_number> &src) const;

  /**
   * 矩阵-向量乘法。就像前面的函数一样，但只适用于矩阵只有一个块列的情况。
   *
   */
  template <typename block_number, typename nonblock_number>
  void
  vmult(BlockVector<block_number> &    dst,
        const Vector<nonblock_number> &src) const;

  /**
   * 矩阵-向量乘法。就像前面的函数，但只适用于矩阵只有一个块行的情况。
   *
   */
  template <typename block_number, typename nonblock_number>
  void
  vmult(Vector<nonblock_number> &        dst,
        const BlockVector<block_number> &src) const;

  /**
   * 矩阵-向量乘法。就像前面的函数，但只适用于矩阵只有一个块的情况。
   *
   */
  template <typename nonblock_number>
  void
  vmult(Vector<nonblock_number> &dst, const Vector<nonblock_number> &src) const;

  /**
   * 矩阵-向量乘法：让 $dst = M^T*src$ 与 $M$
   * 为这个矩阵。这个函数与vmult()的作用相同，但需要转置的矩阵。
   *
   */
  template <typename block_number>
  void
  Tvmult(BlockVector<block_number> &      dst,
         const BlockVector<block_number> &src) const;

  /**
   * 矩阵-向量乘法。就像前面的函数一样，但只适用于矩阵只有一个块行的情况。
   *
   */
  template <typename block_number, typename nonblock_number>
  void
  Tvmult(BlockVector<block_number> &    dst,
         const Vector<nonblock_number> &src) const;

  /**
   * 矩阵-向量乘法。就像前面的函数一样，但只适用于矩阵只有一个块列的情况。
   *
   */
  template <typename block_number, typename nonblock_number>
  void
  Tvmult(Vector<nonblock_number> &        dst,
         const BlockVector<block_number> &src) const;

  /**
   * 矩阵-向量乘法。就像前面的函数，但只适用于矩阵只有一个块的情况。
   *
   */
  template <typename nonblock_number>
  void
  Tvmult(Vector<nonblock_number> &      dst,
         const Vector<nonblock_number> &src) const;
  //@}

  /**
   * @name  预处理方法
   *
   */
  //@{
  /**
   * 应用Jacobi预处理，它将<tt>src</tt>向量的每个元素都乘以各自对角线元素的逆值，并将结果与松弛参数<tt>omega</tt>相乘。
   * 所有的对角线块必须是方形矩阵，才能进行这个操作。
   *
   */
  template <class BlockVectorType>
  void
  precondition_Jacobi(BlockVectorType &      dst,
                      const BlockVectorType &src,
                      const number           omega = 1.) const;

  /**
   * 对一个简单的向量应用雅可比预处理程序。
   * 为此，矩阵必须是单一的正方形块。
   *
   */
  template <typename number2>
  void
  precondition_Jacobi(Vector<number2> &      dst,
                      const Vector<number2> &src,
                      const number           omega = 1.) const;
  //@}

  /**
   * @name  输入/输出
   *
   */
  //@{
  /**
   * 以通常的格式打印矩阵，即作为矩阵而不是作为非零元素的列表。为了提高可读性，不在矩阵中的元素显示为空白，而明确设置为零的矩阵元素则显示为空白。
   * 参数允许对输出格式进行灵活设置。
   * <tt>precision</tt>和<tt>scientific</tt>用于确定数字格式，其中<tt>scientific
   * = false</tt>表示固定点符号。
   * <tt>width</tt>的一个零条目使函数计算出一个宽度，但如果输出粗略的话，可以将其改为一个正值。
   * 此外，还可以指定一个空值的字符。
   * 最后，整个矩阵可以与一个共同的分母相乘，产生更可读的输出，甚至是整数。
   * @attention
   * 如果应用于一个大的矩阵，这个函数可能会产生<b>large</b>量的输出!
   *
   */
  void
  print_formatted(std::ostream &     out,
                  const unsigned int precision   = 3,
                  const bool         scientific  = true,
                  const unsigned int width       = 0,
                  const char *       zero_string = " ",
                  const double       denominator = 1.) const;
  //@}
  /**
   * @addtogroup  Exceptions 
     * @{ 
   *
   */

  /**
   * 异常情况
   *
   */
  DeclException0(ExcBlockDimensionMismatch);
  //@}

private:
  /**
   * 指向用于该矩阵的块状稀疏模式的指针。为了保证它在使用中不被删除，我们使用SmartPointer类来订阅它。
   *
   */
  SmartPointer<const BlockSparsityPattern, BlockSparseMatrix<number>>
    sparsity_pattern;
};



 /*@}*/ 
 /* ------------------------- Template functions ---------------------- */ 



template <typename number>
inline BlockSparseMatrix<number> &
BlockSparseMatrix<number>::operator=(const double d)
{
  Assert(d == 0, ExcScalarAssignmentOnlyForZeroValue());

  for (size_type r = 0; r < this->n_block_rows(); ++r)
    for (size_type c = 0; c < this->n_block_cols(); ++c)
      this->block(r, c) = d;

  return *this;
}



template <typename number>
template <typename block_number>
inline void
BlockSparseMatrix<number>::vmult(BlockVector<block_number> &      dst,
                                 const BlockVector<block_number> &src) const
{
  BaseClass::vmult_block_block(dst, src);
}



template <typename number>
template <typename block_number, typename nonblock_number>
inline void
BlockSparseMatrix<number>::vmult(BlockVector<block_number> &    dst,
                                 const Vector<nonblock_number> &src) const
{
  BaseClass::vmult_block_nonblock(dst, src);
}



template <typename number>
template <typename block_number, typename nonblock_number>
inline void
BlockSparseMatrix<number>::vmult(Vector<nonblock_number> &        dst,
                                 const BlockVector<block_number> &src) const
{
  BaseClass::vmult_nonblock_block(dst, src);
}



template <typename number>
template <typename nonblock_number>
inline void
BlockSparseMatrix<number>::vmult(Vector<nonblock_number> &      dst,
                                 const Vector<nonblock_number> &src) const
{
  BaseClass::vmult_nonblock_nonblock(dst, src);
}



template <typename number>
template <typename block_number>
inline void
BlockSparseMatrix<number>::Tvmult(BlockVector<block_number> &      dst,
                                  const BlockVector<block_number> &src) const
{
  BaseClass::Tvmult_block_block(dst, src);
}



template <typename number>
template <typename block_number, typename nonblock_number>
inline void
BlockSparseMatrix<number>::Tvmult(BlockVector<block_number> &    dst,
                                  const Vector<nonblock_number> &src) const
{
  BaseClass::Tvmult_block_nonblock(dst, src);
}



template <typename number>
template <typename block_number, typename nonblock_number>
inline void
BlockSparseMatrix<number>::Tvmult(Vector<nonblock_number> &        dst,
                                  const BlockVector<block_number> &src) const
{
  BaseClass::Tvmult_nonblock_block(dst, src);
}



template <typename number>
template <typename nonblock_number>
inline void
BlockSparseMatrix<number>::Tvmult(Vector<nonblock_number> &      dst,
                                  const Vector<nonblock_number> &src) const
{
  BaseClass::Tvmult_nonblock_nonblock(dst, src);
}



template <typename number>
template <class BlockVectorType>
inline void
BlockSparseMatrix<number>::precondition_Jacobi(BlockVectorType &      dst,
                                               const BlockVectorType &src,
                                               const number omega) const
{
  Assert(this->n_block_rows() == this->n_block_cols(), ExcNotQuadratic());
  Assert(dst.n_blocks() == this->n_block_rows(),
         ExcDimensionMismatch(dst.n_blocks(), this->n_block_rows()));
  Assert(src.n_blocks() == this->n_block_cols(),
         ExcDimensionMismatch(src.n_blocks(), this->n_block_cols()));

  // do a diagonal preconditioning. uses only
  // the diagonal blocks of the matrix
  for (size_type i = 0; i < this->n_block_rows(); ++i)
    this->block(i, i).precondition_Jacobi(dst.block(i), src.block(i), omega);
}



template <typename number>
template <typename number2>
inline void
BlockSparseMatrix<number>::precondition_Jacobi(Vector<number2> &      dst,
                                               const Vector<number2> &src,
                                               const number omega) const
{
  // check number of blocks. the sizes of the
  // single block is checked in the function
  // we call
  Assert(this->n_block_cols() == 1,
         ExcMessage("This function only works if the matrix has "
                    "a single block"));
  Assert(this->n_block_rows() == 1,
         ExcMessage("This function only works if the matrix has "
                    "a single block"));

  // do a diagonal preconditioning. uses only
  // the diagonal blocks of the matrix
  this->block(0, 0).precondition_Jacobi(dst, src, omega);
}


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_block_sparse_matrix_h


