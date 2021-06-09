//include/deal.II-translator/lac/trilinos_linear_operator_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

#ifndef dealii_trilinos_linear_operator_h
#define dealii_trilinos_linear_operator_h

#include <deal.II/base/config.h>

#if defined(DEAL_II_WITH_TRILINOS)

#  include <deal.II/lac/block_linear_operator.h>
#  include <deal.II/lac/linear_operator.h>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  // Forward declarations:
#  ifndef DOXYGEN
  class SparseMatrix;
  class PreconditionBase;
  class BlockSparseMatrix;

  namespace internal
  {
    namespace LinearOperatorImplementation
    {
      class TrilinosPayload;
    }

    namespace BlockLinearOperatorImplementation
    {
      template <typename PayloadBlockType>
      class TrilinosBlockPayload;
    }
  } // namespace internal
#  endif

  /**
   * @name  创建一个LinearOperator
   *
   */
  //@{


  /**
   * @relatesalso  LinearOperator 一个封装通用 @p matrix
   * 对象的函数，基于一个 @p operator_exemplar,
   * ，作用于一个兼容的Vector类型，变成一个LinearOperator。
   * 这个函数等同于 dealii::linear_operator,
   * ，但通过预选适当的模板参数，确保与Trilinos操作完全兼容。
   * @ingroup TrilinosWrappers
   *
   */
  template <typename Range, typename Domain = Range, typename Matrix>
  inline LinearOperator<
    Range,
    Domain,
    TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload>
  linear_operator(const TrilinosWrappers::SparseMatrix &operator_exemplar,
                  const Matrix &                        matrix)
  {
    using OperatorExemplar = TrilinosWrappers::SparseMatrix;
    using Payload =
      TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload;
    return dealii::
      linear_operator<Range, Domain, Payload, OperatorExemplar, Matrix>(
        operator_exemplar, matrix);
  }


  /**
   * @relatesalso  LinearOperator
   * 一个将作用于兼容矢量类型的通用 @p matrix
   * 对象封装为LinearOperator的函数。    这个函数等同于
   * dealii::linear_operator,
   * ，但通过预先选择适当的模板参数，确保与Trilinos操作完全兼容。
   * @ingroup TrilinosWrappers
   *
   */
  template <typename Range, typename Domain = Range>
  inline LinearOperator<
    Range,
    Domain,
    TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload>
  linear_operator(const TrilinosWrappers::SparseMatrix &matrix)
  {
    using Matrix = TrilinosWrappers::SparseMatrix;
    using Payload =
      TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload;
    return dealii::linear_operator<Range, Domain, Payload, Matrix, Matrix>(
      matrix, matrix);
  }


  //@}
  /**
   * @name  创建一个BlockLinearOperator
   *
   */
  //@{


  /**
   * @relatesalso  BlockLinearOperator 一个将 @p block_matrix
   * 封装为BlockLinearOperator的函数。    这个函数等同于
   * dealii::block_operator,
   * ，但通过预选适当的模板参数，确保与Trilinos操作完全兼容。
   * @ingroup TrilinosWrappers
   *
   */
  template <typename Range, typename Domain = Range>
  inline BlockLinearOperator<
    Range,
    Domain,
    TrilinosWrappers::internal::BlockLinearOperatorImplementation::
      TrilinosBlockPayload<TrilinosWrappers::internal::
                             LinearOperatorImplementation::TrilinosPayload>>
  block_operator(const TrilinosWrappers::BlockSparseMatrix &block_matrix)
  {
    using BlockMatrix = TrilinosWrappers::BlockSparseMatrix;
    using PayloadBlockType =
      TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload;
    using BlockPayload = TrilinosWrappers::internal::
      BlockLinearOperatorImplementation::TrilinosBlockPayload<PayloadBlockType>;
    return dealii::block_operator<Range, Domain, BlockPayload, BlockMatrix>(
      block_matrix);
  }


  /**
   * @relatesalso  BlockLinearOperator
   * 上述函数的一个变体，从对角线元素的数组 @p ops
   * 中建立一个块状对角线线性算子（非对角线块被假设为0）。
   * 这个函数等同于 dealii::block_operator,
   * ，但通过预选适当的模板参数，确保与Trilinos操作完全兼容。
   * @ingroup TrilinosWrappers
   *
   */
  template <std::size_t m,
            std::size_t n,
            typename Range,
            typename Domain = Range>
  inline BlockLinearOperator<
    Range,
    Domain,
    TrilinosWrappers::internal::BlockLinearOperatorImplementation::
      TrilinosBlockPayload<TrilinosWrappers::internal::
                             LinearOperatorImplementation::TrilinosPayload>>
  block_operator(
    const std::array<
      std::array<
        LinearOperator<typename Range::BlockType,
                       typename Domain::BlockType,
                       TrilinosWrappers::internal::
                         LinearOperatorImplementation::TrilinosPayload>,
        n>,
      m> &ops)
  {
    using PayloadBlockType =
      TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload;
    using BlockPayload = TrilinosWrappers::internal::
      BlockLinearOperatorImplementation::TrilinosBlockPayload<PayloadBlockType>;
    return dealii::block_operator<m, n, Range, Domain, BlockPayload>(ops);
  }


  /**
   * @relatesalso  BlockLinearOperator 这个函数提取 @p block_matrix
   * 的对角线块（可以是块状矩阵类型或BlockLinearOperator），并以对角线创建BlockLinearOperator。对角线外的元素被初始化为null_operator（有正确的
   * reinit_range_vector 和 reinit_domain_vector 方法）。
   * 这个函数等同于 dealii::block_diagonal_operator,
   * ，但通过预选适当的模板参数，确保与Trilinos操作完全兼容。
   * @ingroup TrilinosWrappers
   *
   */
  template <typename Range, typename Domain = Range>
  inline BlockLinearOperator<
    Range,
    Domain,
    TrilinosWrappers::internal::BlockLinearOperatorImplementation::
      TrilinosBlockPayload<TrilinosWrappers::internal::
                             LinearOperatorImplementation::TrilinosPayload>>
  block_diagonal_operator(
    const TrilinosWrappers::BlockSparseMatrix &block_matrix)
  {
    using BlockMatrix = TrilinosWrappers::BlockSparseMatrix;
    using PayloadBlockType =
      TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload;
    using BlockPayload = TrilinosWrappers::internal::
      BlockLinearOperatorImplementation::TrilinosBlockPayload<PayloadBlockType>;
    return dealii::
      block_diagonal_operator<Range, Domain, BlockPayload, BlockMatrix>(
        block_matrix);
  }


  /**
   * @relatesalso  BlockLinearOperator
   * 上述函数的一个变体，从对角线元素的数组 @p ops
   * 中建立一个块状对角线线性算子（非对角线块被假设为0）。
   * 这个函数等同于 dealii::block_diagonal_operator,
   * ，但通过预选适当的模板参数，确保与Trilinos操作完全兼容。
   * @ingroup TrilinosWrappers
   *
   */
  template <std::size_t m, typename Range, typename Domain = Range>
  inline BlockLinearOperator<
    Range,
    Domain,
    TrilinosWrappers::internal::BlockLinearOperatorImplementation::
      TrilinosBlockPayload<TrilinosWrappers::internal::
                             LinearOperatorImplementation::TrilinosPayload>>
  block_diagonal_operator(
    const std::array<
      LinearOperator<typename Range::BlockType,
                     typename Domain::BlockType,
                     TrilinosWrappers::internal::LinearOperatorImplementation::
                       TrilinosPayload>,
      m> &ops)
  {
    using PayloadBlockType =
      TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload;
    using BlockPayload = TrilinosWrappers::internal::
      BlockLinearOperatorImplementation::TrilinosBlockPayload<PayloadBlockType>;
    return dealii::block_diagonal_operator<m, Range, Domain, BlockPayload>(ops);
  }

  //@}

} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_TRILINOS
#endif


