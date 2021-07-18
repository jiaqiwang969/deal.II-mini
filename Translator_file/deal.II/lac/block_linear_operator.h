//include/deal.II-translator/lac/block_linear_operator_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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

#ifndef dealii_block_linear_operator_h
#define dealii_block_linear_operator_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/linear_operator.h>


DEAL_II_NAMESPACE_OPEN

// Forward declarations:
#ifndef DOXYGEN
namespace internal
{
  namespace BlockLinearOperatorImplementation
  {
    template <typename PayloadBlockType =
                internal::LinearOperatorImplementation::EmptyPayload>
    class EmptyBlockPayload;
  }
} // namespace internal

template <typename Number>
class BlockVector;

template <typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>>
class BlockLinearOperator;
#endif

template <typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>,
          typename BlockMatrixType>
BlockLinearOperator<Range, Domain, BlockPayload>
block_operator(const BlockMatrixType &matrix);

template <std::size_t m,
          std::size_t n,
          typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>>
BlockLinearOperator<Range, Domain, BlockPayload>
block_operator(
  const std::array<std::array<LinearOperator<typename Range::BlockType,
                                             typename Domain::BlockType,
                                             typename BlockPayload::BlockType>,
                              n>,
                   m> &);

template <std::size_t m,
          typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>>
BlockLinearOperator<Range, Domain, BlockPayload>
block_diagonal_operator(
  const std::array<LinearOperator<typename Range::BlockType,
                                  typename Domain::BlockType,
                                  typename BlockPayload::BlockType>,
                   m> &);

template <std::size_t m,
          typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>>
BlockLinearOperator<Range, Domain, BlockPayload>
block_diagonal_operator(
  const LinearOperator<typename Range::BlockType,
                       typename Domain::BlockType,
                       typename BlockPayload::BlockType> &op);



/**
 * 一个用于存储块状线性运算符概念的类。
 * 这个类在LinearOperator（封装了 @p Matrix
 * 接口）的接口上增加了三个额外的功能。
 *
 * @code
 * std::function<unsigned int()> n_block_rows;
 * std::function<unsigned int()> n_block_cols;
 * std::function<BlockType(unsigned int, unsigned int)> block;
 * @endcode
 * 描述（原本不透明的）线性运算符的底层块结构。
 * BlockLinearOperator类型的对象可以通过包装函数与LinearOperator类似地创建。
 *
 * @code
 * dealii::BlockSparseMatrix<double> A;
 * const auto block_op_a = block_operator(A);
 * @endcode
 *
 * 另外，还有几个辅助函数可用于从可能不同类型的多个独立矩阵中创建实例。下面是一个由FullMatrix和SparseMatrixEZ创建的块状对角线矩阵的例子。
 *
 *
 * @code
 * FullMatrix<double> top_left(2, 2);
 * top_left(0, 0) = 2.0;
 * top_left(0, 1) =
 *
 * -1.0;
 * top_left(1, 0) =
 *
 * -1.0;
 * top_left(1, 1) = 2.0;
 *
 * SparseMatrixEZ<double> bottom_right(4, 4, 4);
 * for (std::size_t row_n = 0; row_n < 4; ++row_n)
 * {
 *   bottom_right.add(row_n, row_n, 1.0);
 *   if (row_n < 3)
 *     bottom_right.add(row_n, row_n + 1,
 *
 * -1.0);
 * }
 *
 * auto top_left_op = linear_operator(top_left);
 * auto bottom_right_op = linear_operator(bottom_right);
 * std::array<decltype(top_left_op), 2> operators {{top_left_op,
 *                                                bottom_right_op}};
 * auto block_op = block_diagonal_operator (operators);
 *
 * std::vector<BlockVector<double>::size_type> block_sizes {2, 4};
 * BlockVector<double> src(block_sizes);
 * src = 2.0;
 * BlockVector<double> dst(block_sizes);
 * block_op.vmult(dst, src); // now equal to 2, 2, 0, 0, 0, 2
 * @endcode
 *
 *
 * 一个BlockLinearOperator可以在任何时候被切成一个LinearOperator。这将删除所有关于底层块结构的信息（因为上述
 * <code>std::function</code> 对象不再可用）。
 *
 * - 线性操作符的接口，然而，仍然是完整的。
 *
 *
 * @note  这个类大量使用了 <code>std::function</code>
 * 对象和lambda函数。这种灵活性伴随着运行时间的惩罚。只使用这个对象来封装具有中到大的单个块大小和小块结构的对象（作为经验法则，矩阵块大于
 * $1000\times1000$  ）。
 *
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range, typename Domain, typename BlockPayload>
class BlockLinearOperator
  : public LinearOperator<Range, Domain, typename BlockPayload::BlockType>
{
public:
  using BlockType = LinearOperator<typename Range::BlockType,
                                   typename Domain::BlockType,
                                   typename BlockPayload::BlockType>;

  /**
   * 创建一个空的BlockLinearOperator对象。
   * All<code>std::function</code>
   * 这个类和它的基类LinearOperator的成员对象被初始化为默认变体，在调用时抛出一个异常。
   *
   */
  BlockLinearOperator(const BlockPayload &payload)
    : LinearOperator<Range, Domain, typename BlockPayload::BlockType>(
        typename BlockPayload::BlockType(payload, payload))
  {
    n_block_rows = []() -> unsigned int {
      Assert(
        false,
        ExcMessage(
          "Uninitialized BlockLinearOperator<Range, Domain>::n_block_rows called"));
      return 0;
    };

    n_block_cols = []() -> unsigned int {
      Assert(
        false,
        ExcMessage(
          "Uninitialized BlockLinearOperator<Range, Domain>::n_block_cols called"));
      return 0;
    };

    block = [](unsigned int, unsigned int) -> BlockType {
      Assert(
        false,
        ExcMessage(
          "Uninitialized BlockLinearOperator<Range, Domain>::block called"));
      return BlockType();
    };
  }

  /**
   * 默认的复制构造函数。
   *
   */
  BlockLinearOperator(
    const BlockLinearOperator<Range, Domain, BlockPayload> &) = default;

  /**
   * 模板化的复制构造函数，从一个定义了转换函数
   * <code>block_operator</code> 的对象 @p op
   * 中创建一个BlockLinearOperator对象。
   *
   */
  template <typename Op>
  BlockLinearOperator(const Op &op)
  {
    *this = block_operator<Range, Domain, BlockPayload, Op>(op);
  }

  /**
   * 从一个二维的LinearOperator数组 @p ops
   * 中创建一个BlockLinearOperator。这个构造函数调用相应的block_operator()专用化。
   *
   */
  template <std::size_t m, std::size_t n>
  BlockLinearOperator(const std::array<std::array<BlockType, n>, m> &ops)
  {
    *this = block_operator<m, n, Range, Domain, BlockPayload>(ops);
  }

  /**
   * 从一维数组 @p ops
   * 的LinearOperator创建一个块对角线的BlockLinearOperator。这个构造函数调用相应的block_operator()专用化。
   *
   */
  template <std::size_t m>
  BlockLinearOperator(const std::array<BlockType, m> &ops)
  {
    *this = block_diagonal_operator<m, Range, Domain, BlockPayload>(ops);
  }

  /**
   * 默认的复制赋值运算器。
   *
   */
  BlockLinearOperator<Range, Domain, BlockPayload> &
  operator=(const BlockLinearOperator<Range, Domain, BlockPayload> &) = default;

  /**
   * 为一个定义了转换函数 <code>block_operator</code> 的对象 @p
   * op 模板化的复制赋值运算符。
   *
   */
  template <typename Op>
  BlockLinearOperator<Range, Domain, BlockPayload> &
  operator=(const Op &op)
  {
    *this = block_operator<Range, Domain, BlockPayload, Op>(op);
    return *this;
  }

  /**
   * 从一个二维数组 @p ops 的LinearOperator中复制赋值。
   * 这个赋值运算符调用相应的block_operator()专用化。
   *
   */
  template <std::size_t m, std::size_t n>
  BlockLinearOperator<Range, Domain, BlockPayload> &
  operator=(const std::array<std::array<BlockType, n>, m> &ops)
  {
    *this = block_operator<m, n, Range, Domain, BlockPayload>(ops);
    return *this;
  }

  /**
   * 从一维数组 @p ops
   * 的LinearOperator复制赋值，创建一个块对角线的BlockLinearOperator。这个赋值运算符调用相应的block_operator()专用化。
   *
   */
  template <std::size_t m>
  BlockLinearOperator<Range, Domain, BlockPayload> &
  operator=(const std::array<BlockType, m> &ops)
  {
    *this = block_diagonal_operator<m, Range, Domain, BlockPayload>(ops);
    return *this;
  }

  /**
   * 返回一列中的块数（即 "块行 "的数量，或者数字 $m$
   * ，如果解释为 $m\times n$ 块系统）。
   *
   */
  std::function<unsigned int()> n_block_rows;

  /**
   * 返回一行的块数（即 "块列 "的数量，或者数字 $n$
   * ，如果解释为 $m\times n$ 块系统）。
   *
   */
  std::function<unsigned int()> n_block_cols;

  /**
   * 访问具有给定坐标的块。这个 <code>std::function</code>
   * 对象返回一个LinearOperator，代表BlockLinearOperator的 $(i,j)$
   * -th块。
   *
   */
  std::function<BlockType(unsigned int, unsigned int)> block;
};


namespace internal
{
  namespace BlockLinearOperatorImplementation
  {
    // A helper function to apply a given vmult, or Tvmult to a vector with
    // intermediate storage, similar to the corresponding helper
    // function for LinearOperator. Here, two operators are used.
    // The first one takes care of the first "column" and typically doesn't add.
    // On the other hand, the second operator is normally an adding one.
    template <typename Function1,
              typename Function2,
              typename Range,
              typename Domain>
    void
    apply_with_intermediate_storage(const Function1 &first_op,
                                    const Function2 &loop_op,
                                    Range &          v,
                                    const Domain &   u,
                                    bool             add)
    {
      GrowingVectorMemory<Range> vector_memory;

      typename VectorMemory<Range>::Pointer tmp(vector_memory);
      tmp->reinit(v,  /*bool omit_zeroing_entries =*/ true);

      const unsigned int n = u.n_blocks();
      const unsigned int m = v.n_blocks();

      for (unsigned int i = 0; i < m; ++i)
        {
          first_op(*tmp, u, i, 0);
          for (unsigned int j = 1; j < n; ++j)
            loop_op(*tmp, u, i, j);
        }

      if (add)
        v += *tmp;
      else
        v = *tmp;
    }

    // Populate the LinearOperator interfaces with the help of the
    // BlockLinearOperator functions
    template <typename Range, typename Domain, typename BlockPayload>
    inline void
    populate_linear_operator_functions(
      dealii::BlockLinearOperator<Range, Domain, BlockPayload> &op)
    {
      op.reinit_range_vector = [=](Range &v, bool omit_zeroing_entries) {
        const unsigned int m = op.n_block_rows();

        // Reinitialize the block vector to m blocks:
        v.reinit(m);

        // And reinitialize every individual block with reinit_range_vectors:
        for (unsigned int i = 0; i < m; ++i)
          op.block(i, 0).reinit_range_vector(v.block(i), omit_zeroing_entries);

        v.collect_sizes();
      };

      op.reinit_domain_vector = [=](Domain &v, bool omit_zeroing_entries) {
        const unsigned int n = op.n_block_cols();

        // Reinitialize the block vector to n blocks:
        v.reinit(n);

        // And reinitialize every individual block with reinit_domain_vectors:
        for (unsigned int i = 0; i < n; ++i)
          op.block(0, i).reinit_domain_vector(v.block(i), omit_zeroing_entries);

        v.collect_sizes();
      };

      op.vmult = [&op](Range &v, const Domain &u) {
        const unsigned int m = op.n_block_rows();
        const unsigned int n = op.n_block_cols();
        Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
        Assert(u.n_blocks() == n, ExcDimensionMismatch(u.n_blocks(), n));

        if (PointerComparison::equal(&v, &u))
          {
            const auto first_op = [&op](Range &            v,
                                        const Domain &     u,
                                        const unsigned int i,
                                        const unsigned int j) {
              op.block(i, j).vmult(v.block(i), u.block(j));
            };

            const auto loop_op = [&op](Range &            v,
                                       const Domain &     u,
                                       const unsigned int i,
                                       const unsigned int j) {
              op.block(i, j).vmult_add(v.block(i), u.block(j));
            };

            apply_with_intermediate_storage(first_op, loop_op, v, u, false);
          }
        else
          {
            for (unsigned int i = 0; i < m; ++i)
              {
                op.block(i, 0).vmult(v.block(i), u.block(0));
                for (unsigned int j = 1; j < n; ++j)
                  op.block(i, j).vmult_add(v.block(i), u.block(j));
              }
          }
      };

      op.vmult_add = [&op](Range &v, const Domain &u) {
        const unsigned int m = op.n_block_rows();
        const unsigned int n = op.n_block_cols();
        Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
        Assert(u.n_blocks() == n, ExcDimensionMismatch(u.n_blocks(), n));

        if (PointerComparison::equal(&v, &u))
          {
            const auto first_op = [&op](Range &            v,
                                        const Domain &     u,
                                        const unsigned int i,
                                        const unsigned int j) {
              op.block(i, j).vmult(v.block(i), u.block(j));
            };

            const auto loop_op = [&op](Range &            v,
                                       const Domain &     u,
                                       const unsigned int i,
                                       const unsigned int j) {
              op.block(i, j).vmult_add(v.block(i), u.block(j));
            };

            apply_with_intermediate_storage(first_op, loop_op, v, u, true);
          }
        else
          {
            for (unsigned int i = 0; i < m; ++i)
              for (unsigned int j = 0; j < n; ++j)
                op.block(i, j).vmult_add(v.block(i), u.block(j));
          }
      };

      op.Tvmult = [&op](Domain &v, const Range &u) {
        const unsigned int n = op.n_block_cols();
        const unsigned int m = op.n_block_rows();
        Assert(v.n_blocks() == n, ExcDimensionMismatch(v.n_blocks(), n));
        Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

        if (PointerComparison::equal(&v, &u))
          {
            const auto first_op = [&op](Range &            v,
                                        const Domain &     u,
                                        const unsigned int i,
                                        const unsigned int j) {
              op.block(j, i).Tvmult(v.block(i), u.block(j));
            };

            const auto loop_op = [&op](Range &            v,
                                       const Domain &     u,
                                       const unsigned int i,
                                       const unsigned int j) {
              op.block(j, i).Tvmult_add(v.block(i), u.block(j));
            };

            apply_with_intermediate_storage(first_op, loop_op, v, u, false);
          }
        else
          {
            for (unsigned int i = 0; i < n; ++i)
              {
                op.block(0, i).Tvmult(v.block(i), u.block(0));
                for (unsigned int j = 1; j < m; ++j)
                  op.block(j, i).Tvmult_add(v.block(i), u.block(j));
              }
          }
      };

      op.Tvmult_add = [&op](Domain &v, const Range &u) {
        const unsigned int n = op.n_block_cols();
        const unsigned int m = op.n_block_rows();
        Assert(v.n_blocks() == n, ExcDimensionMismatch(v.n_blocks(), n));
        Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

        if (PointerComparison::equal(&v, &u))
          {
            const auto first_op = [&op](Range &            v,
                                        const Domain &     u,
                                        const unsigned int i,
                                        const unsigned int j) {
              op.block(j, i).Tvmult(v.block(i), u.block(j));
            };

            const auto loop_op = [&op](Range &            v,
                                       const Domain &     u,
                                       const unsigned int i,
                                       const unsigned int j) {
              op.block(j, i).Tvmult_add(v.block(i), u.block(j));
            };

            apply_with_intermediate_storage(first_op, loop_op, v, u, true);
          }
        else
          {
            for (unsigned int i = 0; i < n; ++i)
              for (unsigned int j = 0; j < m; ++j)
                op.block(j, i).Tvmult_add(v.block(i), u.block(j));
          }
      };
    }



    /**
     * 一个用于BlockLinearOperator的假类，不需要任何扩展来促进块矩阵或其子块的操作。
     * 这是通常与deal.II的本地BlockSparseMatrix相关的Payload类。要使用
     * TrilinosWrappers::BlockSparseMatrix 或
     * PETScWrappers::BlockSparseMatrix
     * ，必须用它们相关的BlockPayload初始化BlockLinearOperator。
     * @ingroup LAOperators
     *
     */
    template <typename PayloadBlockType>
    class EmptyBlockPayload
    {
    public:
      /**
       * 每个子块所持有的有效载荷的类型
       *
       */
      using BlockType = PayloadBlockType;

      /**
       * 默认构造函数
       * 由于这个类不做任何特别的事情，不需要特别的配置，我们只有一个通用的构造函数，可以在任何条件下调用。
       *
       */
      template <typename... Args>
      EmptyBlockPayload(const Args &...)
      {}
    };

  } // namespace BlockLinearOperatorImplementation
} // namespace internal



/**
 * @name  创建一个BlockLinearOperator
 *
 *
 */
//@{

/**
 * @relatesalso  BlockLinearOperator 一个将 @p block_matrix
 * 封装为BlockLinearOperator的函数。
 * 在BlockLinearOperator对象创建后，对 @p
 * block_matrix的块结构和单个块所做的所有改变都会被操作者对象所反映。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range,
          typename Domain,
          typename BlockPayload,
          typename BlockMatrixType>
BlockLinearOperator<Range, Domain, BlockPayload>
block_operator(const BlockMatrixType &block_matrix)
{
  using BlockType =
    typename BlockLinearOperator<Range, Domain, BlockPayload>::BlockType;

  BlockLinearOperator<Range, Domain, BlockPayload> return_op{
    BlockPayload(block_matrix, block_matrix)};

  return_op.n_block_rows = [&block_matrix]() -> unsigned int {
    return block_matrix.n_block_rows();
  };

  return_op.n_block_cols = [&block_matrix]() -> unsigned int {
    return block_matrix.n_block_cols();
  };

  return_op.block = [&block_matrix](unsigned int i,
                                    unsigned int j) -> BlockType {
#ifdef DEBUG
    const unsigned int m = block_matrix.n_block_rows();
    const unsigned int n = block_matrix.n_block_cols();
    AssertIndexRange(i, m);
    AssertIndexRange(j, n);
#endif

    return BlockType(block_matrix.block(i, j));
  };

  populate_linear_operator_functions(return_op);
  return return_op;
}



/**
 * @relatesalso  BlockLinearOperator
 * 上述函数的一个变体，它将一个给定的LinearOperator集合 @p
 * ops
 * 封装成一个块结构。这里假定Range和Domain是块状向量，即从
 * @ref BlockVectorBase 中导出。 @p ops
 * 中的各个线性运算符必须作用于块向量的底层向量类型，即在
 * Domain::BlockType 中产生一个 Range::BlockType. 的结果。 列表 @p
 * ops
 * 最好作为初始化器列表传递。例如，考虑一个线性运算块（作用于Vector<double>）。
 *
 * @code
 * op_a00 | op_a01
 *       |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * ---------------
 *       |
 * op_a10 | op_a11
 * @endcode
 * 相应的block_operator调用的形式是
 *
 * @code
 * block_operator<2, 2, BlockVector<double>>({op_a00, op_a01, op_a10, op_a11});
 * @endcode
 *
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <std::size_t m,
          std::size_t n,
          typename Range,
          typename Domain,
          typename BlockPayload>
BlockLinearOperator<Range, Domain, BlockPayload>
block_operator(
  const std::array<std::array<LinearOperator<typename Range::BlockType,
                                             typename Domain::BlockType,
                                             typename BlockPayload::BlockType>,
                              n>,
                   m> &ops)
{
  static_assert(m > 0 && n > 0,
                "a blocked LinearOperator must consist of at least one block");

  using BlockType =
    typename BlockLinearOperator<Range, Domain, BlockPayload>::BlockType;

  // TODO: Create block payload so that this can be initialized correctly
  BlockLinearOperator<Range, Domain, BlockPayload> return_op{BlockPayload()};

  return_op.n_block_rows = []() -> unsigned int { return m; };

  return_op.n_block_cols = []() -> unsigned int { return n; };

  return_op.block = [ops](unsigned int i, unsigned int j) -> BlockType {
    AssertIndexRange(i, m);
    AssertIndexRange(j, n);

    return ops[i][j];
  };

  populate_linear_operator_functions(return_op);
  return return_op;
}



/**
 * @relatesalso  BlockLinearOperator 这个函数提取 @p block_matrix
 * 的对角线块（可以是块状矩阵类型，也可以是BlockLinearOperator），并创建一个具有对角线的BlockLinearOperator。对角线外的元素被初始化为null_operator（有正确的
 * reinit_range_vector 和 reinit_domain_vector 方法）。
 * 在创建BlockLinearOperator对象后，在 @p block_matrix
 * 的各个对角线块上所做的所有改变都会被操作者对象所反映。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>,
          typename BlockMatrixType>
BlockLinearOperator<Range, Domain, BlockPayload>
block_diagonal_operator(const BlockMatrixType &block_matrix)
{
  using BlockType =
    typename BlockLinearOperator<Range, Domain, BlockPayload>::BlockType;

  BlockLinearOperator<Range, Domain, BlockPayload> return_op{
    BlockPayload(block_matrix, block_matrix)};

  return_op.n_block_rows = [&block_matrix]() -> unsigned int {
    return block_matrix.n_block_rows();
  };

  return_op.n_block_cols = [&block_matrix]() -> unsigned int {
    return block_matrix.n_block_cols();
  };

  return_op.block = [&block_matrix](unsigned int i,
                                    unsigned int j) -> BlockType {
#ifdef DEBUG
    const unsigned int m = block_matrix.n_block_rows();
    const unsigned int n = block_matrix.n_block_cols();
    Assert(m == n, ExcDimensionMismatch(m, n));
    AssertIndexRange(i, m);
    AssertIndexRange(j, n);
#endif
    if (i == j)
      return BlockType(block_matrix.block(i, j));
    else
      return null_operator(BlockType(block_matrix.block(i, j)));
  };

  populate_linear_operator_functions(return_op);
  return return_op;
}



/**
 * @relatesalso  BlockLinearOperator
 * 上述函数的一个变体，它从对角线元素的数组 @p ops
 * 中建立起一个块状对角线线性运算器（非对角线块被假定为0）。
 * 列表 @p ops
 * 最好作为一个初始化器列表传递。例如考虑一个线性操作块（作用于Vector<double>）
 * <code>diag(op_a0, op_a1, ...,
 * op_am)</code>。相应的block_operator调用的形式是
 *
 * @code
 * block_diagonal_operator<m, BlockVector<double>>({op_00, op_a1, ..., op_am});
 * @endcode
 *
 * @ingroup LAOperators
 *
 *
 */
template <std::size_t m, typename Range, typename Domain, typename BlockPayload>
BlockLinearOperator<Range, Domain, BlockPayload>
block_diagonal_operator(
  const std::array<LinearOperator<typename Range::BlockType,
                                  typename Domain::BlockType,
                                  typename BlockPayload::BlockType>,
                   m> &ops)
{
  static_assert(
    m > 0, "a blockdiagonal LinearOperator must consist of at least one block");

  using BlockType =
    typename BlockLinearOperator<Range, Domain, BlockPayload>::BlockType;

  std::array<std::array<BlockType, m>, m> new_ops;

  // This is a bit tricky. We have to make sure that the off-diagonal
  // elements of return_op.ops are populated correctly. They must be
  // null_operators, but with correct reinit_domain_vector and
  // reinit_range_vector functions.
  for (unsigned int i = 0; i < m; ++i)
    for (unsigned int j = 0; j < m; ++j)
      if (i == j)
        {
          // diagonal elements are easy:
          new_ops[i][j] = ops[i];
        }
      else
        {
          // create a null-operator...
          new_ops[i][j] = null_operator(ops[i]);
          // ... and fix up reinit_domain_vector:
          new_ops[i][j].reinit_domain_vector = ops[j].reinit_domain_vector;
        }

  return block_operator<m, m, Range, Domain>(new_ops);
}



/**
 * @relatesalso  BlockLinearOperator
 * 上述函数的一个变体，它只接受一个LinearOperator参数 @p op
 * ，并创建一个带有 @p m 副本的块对角线性算子。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <std::size_t m, typename Range, typename Domain, typename BlockPayload>
BlockLinearOperator<Range, Domain, BlockPayload>
block_diagonal_operator(
  const LinearOperator<typename Range::BlockType,
                       typename Domain::BlockType,
                       typename BlockPayload::BlockType> &op)
{
  static_assert(m > 0,
                "a blockdiagonal LinearOperator must consist of at least "
                "one block");

  using BlockType =
    typename BlockLinearOperator<Range, Domain, BlockPayload>::BlockType;
  std::array<BlockType, m> new_ops;
  new_ops.fill(op);

  return block_diagonal_operator(new_ops);
}



//@}
/**
 * @name  对一个BlockLinearOperator的操纵
 *
 */
//@{

/**
 * @relatesalso  LinearOperator  @relatesalso  BlockLinearOperator
 * 这个函数实现了正向替换，以反转一个低级块状三角形矩阵。作为参数，它需要一个BlockLinearOperator
 * @p
 * block_operator，代表一个块状下三角矩阵，以及一个BlockLinearOperator
 * @p diagonal_inverse  ，代表 @p block_operator.  对角线块的反转。
 * 让我们假设我们有一个线性系统，其块结构如下。
 *
 *
 * @code
 * A00 x0 + ...                   = y0
 * A01 x0 + A11 x1 + ...          = y1
 * ...        ...
 * A0n x0 + A1n x1 + ... + Ann xn = yn
 * @endcode
 *
 * 首先， <code>x0 = A00^-1 y0</code>
 * 。然后，我们可以用x0来恢复x1。
 *
 * @code
 *  x1 = A11^-1 ( y1
 *
 * - A01 x0 )
 * @endcode
 * 并因此。
 * @code
 *  xn = Ann^-1 ( yn
 *
 * - A0n x0
 *
 * - ...
 *
 * - A(n-1)n x(n-1) )
 * @endcode
 *
 *
 * @note
 * 我们没有使用BlockLinearOperator参数的所有块。只是使用了
 * @p block_operator 的下三角块矩阵以及 @p diagonal_inverse.
 * 的对角线。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>>
LinearOperator<Domain, Range, typename BlockPayload::BlockType>
block_forward_substitution(
  const BlockLinearOperator<Range, Domain, BlockPayload> &block_operator,
  const BlockLinearOperator<Domain, Range, BlockPayload> &diagonal_inverse)
{
  LinearOperator<Range, Range, typename BlockPayload::BlockType> return_op{
    typename BlockPayload::BlockType(diagonal_inverse)};

  return_op.reinit_range_vector  = diagonal_inverse.reinit_range_vector;
  return_op.reinit_domain_vector = diagonal_inverse.reinit_domain_vector;

  return_op.vmult = [block_operator, diagonal_inverse](Range &      v,
                                                       const Range &u) {
    const unsigned int m = block_operator.n_block_rows();
    Assert(block_operator.n_block_cols() == m,
           ExcDimensionMismatch(block_operator.n_block_cols(), m));
    Assert(diagonal_inverse.n_block_rows() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_rows(), m));
    Assert(diagonal_inverse.n_block_cols() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_cols(), m));
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    if (m == 0)
      return;

    diagonal_inverse.block(0, 0).vmult(v.block(0), u.block(0));
    for (unsigned int i = 1; i < m; ++i)
      {
        auto &dst = v.block(i);
        dst       = u.block(i);
        dst *= -1.;
        for (unsigned int j = 0; j < i; ++j)
          block_operator.block(i, j).vmult_add(dst, v.block(j));
        dst *= -1.;
        diagonal_inverse.block(i, i).vmult(dst,
                                           dst); // uses intermediate storage
      }
  };

  return_op.vmult_add = [block_operator, diagonal_inverse](Range &      v,
                                                           const Range &u) {
    const unsigned int m = block_operator.n_block_rows();
    Assert(block_operator.n_block_cols() == m,
           ExcDimensionMismatch(block_operator.n_block_cols(), m));
    Assert(diagonal_inverse.n_block_rows() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_rows(), m));
    Assert(diagonal_inverse.n_block_cols() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_cols(), m));
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    if (m == 0)
      return;

    GrowingVectorMemory<typename Range::BlockType>            vector_memory;
    typename VectorMemory<typename Range::BlockType>::Pointer tmp(
      vector_memory);

    diagonal_inverse.block(0, 0).vmult_add(v.block(0), u.block(0));

    for (unsigned int i = 1; i < m; ++i)
      {
        diagonal_inverse.block(i, i).reinit_range_vector(
          *tmp,  /*bool omit_zeroing_entries=*/ true);
        *tmp = u.block(i);
        *tmp *= -1.;
        for (unsigned int j = 0; j < i; ++j)
          block_operator.block(i, j).vmult_add(*tmp, v.block(j));
        *tmp *= -1.;
        diagonal_inverse.block(i, i).vmult_add(v.block(i), *tmp);
      }
  };

  return return_op;
}



/**
 * @relatesalso  LinearOperator  @relatesalso  BlockLinearOperator
 * 这个函数实现了反置，以反转一个上块三角矩阵。作为参数，它需要一个BlockLinearOperator
 * @p
 * block_operator，代表一个上块三角矩阵，以及一个BlockLinearOperator
 * @p diagonal_inverse  ，代表 @p block_operator.  的对角块的反转。
 * 让我们假设我们有一个线性系统，其块结构如下。
 *
 *
 * @code
 * A00 x0 + A01 x1 + ... + A0n xn = yn
 *        A11 x1 + ...          = y1
 *                        ...     ..
 *                       Ann xn = yn
 * @endcode
 *  首先， <code>xn = Ann^-1 yn</code>
 * 。然后，我们可以用xn来恢复x(n-1)。
 *
 * @code
 *  x(n-1) = A(n-1)(n-1)^-1 ( y(n-1)
 *
 * - A(n-1)n x(n-1) )
 * @endcode
 * 并因此。
 *
 * @code
 *  x0 = A00^-1 ( y0
 *
 * - A0n xn
 *
 * - ...
 *
 * - A01 x1 )
 * @endcode
 *
 *
 * @note
 * 我们没有使用BlockLinearOperator参数的所有块。只是使用了
 * @p block_operator 的上三角块矩阵以及 @p diagonal_inverse.
 * 的对角线。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range  = BlockVector<double>,
          typename Domain = Range,
          typename BlockPayload =
            internal::BlockLinearOperatorImplementation::EmptyBlockPayload<>>
LinearOperator<Domain, Range, typename BlockPayload::BlockType>
block_back_substitution(
  const BlockLinearOperator<Range, Domain, BlockPayload> &block_operator,
  const BlockLinearOperator<Domain, Range, BlockPayload> &diagonal_inverse)
{
  LinearOperator<Range, Range, typename BlockPayload::BlockType> return_op{
    typename BlockPayload::BlockType(diagonal_inverse)};

  return_op.reinit_range_vector  = diagonal_inverse.reinit_range_vector;
  return_op.reinit_domain_vector = diagonal_inverse.reinit_domain_vector;

  return_op.vmult = [block_operator, diagonal_inverse](Range &      v,
                                                       const Range &u) {
    const unsigned int m = block_operator.n_block_rows();
    Assert(block_operator.n_block_cols() == m,
           ExcDimensionMismatch(block_operator.n_block_cols(), m));
    Assert(diagonal_inverse.n_block_rows() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_rows(), m));
    Assert(diagonal_inverse.n_block_cols() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_cols(), m));
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));

    if (m == 0)
      return;

    diagonal_inverse.block(m - 1, m - 1).vmult(v.block(m - 1), u.block(m - 1));

    for (int i = m - 2; i >= 0; --i)
      {
        auto &dst = v.block(i);
        dst       = u.block(i);
        dst *= -1.;
        for (unsigned int j = i + 1; j < m; ++j)
          block_operator.block(i, j).vmult_add(dst, v.block(j));
        dst *= -1.;
        diagonal_inverse.block(i, i).vmult(dst,
                                           dst); // uses intermediate storage
      }
  };

  return_op.vmult_add = [block_operator, diagonal_inverse](Range &      v,
                                                           const Range &u) {
    const unsigned int m = block_operator.n_block_rows();
    Assert(block_operator.n_block_cols() == m,
           ExcDimensionMismatch(block_operator.n_block_cols(), m));
    Assert(diagonal_inverse.n_block_rows() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_rows(), m));
    Assert(diagonal_inverse.n_block_cols() == m,
           ExcDimensionMismatch(diagonal_inverse.n_block_cols(), m));
    Assert(v.n_blocks() == m, ExcDimensionMismatch(v.n_blocks(), m));
    Assert(u.n_blocks() == m, ExcDimensionMismatch(u.n_blocks(), m));
    GrowingVectorMemory<typename Range::BlockType>            vector_memory;
    typename VectorMemory<typename Range::BlockType>::Pointer tmp(
      vector_memory);

    if (m == 0)
      return;

    diagonal_inverse.block(m - 1, m - 1)
      .vmult_add(v.block(m - 1), u.block(m - 1));

    for (int i = m - 2; i >= 0; --i)
      {
        diagonal_inverse.block(i, i).reinit_range_vector(
          *tmp,  /*bool omit_zeroing_entries=*/ true);
        *tmp = u.block(i);
        *tmp *= -1.;
        for (unsigned int j = i + 1; j < m; ++j)
          block_operator.block(i, j).vmult_add(*tmp, v.block(j));
        *tmp *= -1.;
        diagonal_inverse.block(i, i).vmult_add(v.block(i), *tmp);
      }
  };

  return return_op;
}

//@}

DEAL_II_NAMESPACE_CLOSE

#endif


