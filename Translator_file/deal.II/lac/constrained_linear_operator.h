//include/deal.II-translator/lac/constrained_linear_operator_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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

#ifndef dealii_constrained_linear_operator_h
#define dealii_constrained_linear_operator_h

#include <deal.II/base/config.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>


DEAL_II_NAMESPACE_OPEN


/**
 * @name  间接对LinearOperator施加约束
 *
 *
 */
//@{


/**
 * 这个函数接收一个AffineConstraints对象 @p constraints
 * 和一个运算器示例 @p exemplar
 * （这个示例通常是一个描述系统矩阵的线性运算器
 *
 * -它只用于创建适当大小的域和范围向量，它的动作<tt>vmult</tt>从不使用）。)
 * 一个与底层AffineConstraints对象的 "同质动作
 * "相关的LinearOperator对象被返回。 在向量 <code>u</code>
 * 上应用LinearOperator对象的结果是一个向量 <code>v</code>
 * ，它存储了在 <code>u</code> 上调用 AffineConstraints::distribute()
 * 的结果。
 *
 * - 有一个重要的区别：不均匀性不被应用，而是始终被视为0。
 * 这个函数创建的LinearOperator对象主要在内部用于constrained_linear_operator()，以建立一个修正的线性方程组。如何用这种方法解决线性方程组，在 @ref
 * constraints 模块中有详细解释。
 *
 *
 *
 * @note
 * 目前，这个函数对于分布式数据结构可能无法正确工作。
 * @relatesalso  LinearOperator
 * @ingroup constraints
 *
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
distribute_constraints_linear_operator(
  const AffineConstraints<typename Range::value_type> &constraints,
  const LinearOperator<Range, Domain, Payload> &       exemplar)
{
  LinearOperator<Range, Domain, Payload> return_op = exemplar;

  return_op.vmult_add = [&constraints](Range &v, const Domain &u) {
    Assert(!dealii::PointerComparison::equal(&v, &u),
           dealii::ExcMessage("The domain and range vectors must be different "
                              "storage locations"));

    // First, add vector u to v unconditionally and clean up constrained
    // degrees of freedom later.
    v += u;

    const auto &locally_owned_elements = v.locally_owned_elements();
    for (const auto &line : constraints.get_lines())
      {
        const auto i = line.index;
        if (locally_owned_elements.is_element(i))
          {
            v(i) -= u(i);
            const auto &entries = line.entries;
            for (types::global_dof_index j = 0; j < entries.size(); ++j)
              {
                const auto pos = entries[j].first;
                v(i) += u(pos) * entries[j].second;
              }
          }
      }

    v.compress(VectorOperation::add);
  };

  return_op.Tvmult_add = [&constraints](Domain &v, const Range &u) {
    Assert(!dealii::PointerComparison::equal(&v, &u),
           dealii::ExcMessage("The domain and range vectors must be different "
                              "storage locations"));

    // First, add vector u to v unconditionally and clean up constrained
    // degrees of freedom later.
    v += u;

    const auto &locally_owned_elements = v.locally_owned_elements();
    for (const auto &line : constraints.get_lines())
      {
        const auto i = line.index;

        if (locally_owned_elements.is_element(i))
          {
            v(i) -= u(i);
          }

        const auto &entries = line.entries;
        for (types::global_dof_index j = 0; j < entries.size(); ++j)
          {
            const auto pos = entries[j].first;
            if (locally_owned_elements.is_element(pos))
              v(pos) += u(i) * entries[j].second;
          }
      }

    v.compress(VectorOperation::add);
  };

  return_op.vmult = [vmult_add = return_op.vmult_add](Range &       v,
                                                      const Domain &u) {
    v = 0.;
    vmult_add(v, u);
  };

  return_op.Tvmult = [Tvmult_add = return_op.Tvmult_add](Domain &     v,
                                                         const Range &u) {
    v = 0.;
    Tvmult_add(v, u);
  };

  return return_op;
}


/**
 * 给定一个AffineConstraints @p constraints 和一个运算符示例 @p
 * 示例，返回一个LinearOperator，该运算符是对受限自由度子空间的投影，即结果向量中对应于非受限自由度的所有条目被设置为零。
 *
 * @relatesalso  LinearOperator
 *
 * @ingroup constraints
 *
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
project_to_constrained_linear_operator(
  const AffineConstraints<typename Range::value_type> &constraints,
  const LinearOperator<Range, Domain, Payload> &       exemplar)
{
  LinearOperator<Range, Domain, Payload> return_op = exemplar;

  return_op.vmult_add = [&constraints](Range &v, const Domain &u) {
    const auto &locally_owned_elements = v.locally_owned_elements();
    for (const auto &line : constraints.get_lines())
      {
        const auto i = line.index;
        if (locally_owned_elements.is_element(i))
          {
            v(i) += u(i);
          }
      }

    v.compress(VectorOperation::add);
  };

  return_op.Tvmult_add = [&constraints](Domain &v, const Range &u) {
    const auto &locally_owned_elements = v.locally_owned_elements();
    for (const auto &line : constraints.get_lines())
      {
        const auto i = line.index;
        if (locally_owned_elements.is_element(i))
          {
            v(i) += u(i);
          }
      }

    v.compress(VectorOperation::add);
  };

  return_op.vmult = [vmult_add = return_op.vmult_add](Range &       v,
                                                      const Domain &u) {
    v = 0.;
    vmult_add(v, u);
  };

  return_op.Tvmult = [Tvmult_add = return_op.Tvmult_add](Domain &     v,
                                                         const Range &u) {
    v = 0.;
    Tvmult_add(v, u);
  };

  return return_op;
}


/**
 * 给定一个AffineConstraints对象 @p constraints 和一个LinearOperator
 * @p linop,
 * ，该函数创建一个LinearOperator对象，由三个操作和一个正则化组成。
 *
 * @code
 * Ct linop C + Id_c;
 * @endcode
 * 与
 *
 * @code
 * C = distribute_constraints_linear_operator(constraints, linop);
 * Ct = transpose_operator(C);
 * Id_c = project_to_constrained_linear_operator(constraints, linop);
 * @endcode
 * 而 <code>Id_c</code>
 * 是对由所有与受限自由度相关的向量条目组成的子空间的投影。
 * 这个LinearOperator对象与constrained_right_hand_side()一起使用，建立了以下修改后的线性方程系统。@f[
 * (C^T A C + Id_c) x = C^T (b
 *
 * - A\,k)
 * @f] 具有给定的（无约束的）系统矩阵  $A$  ，右手边  $b$  ，以及具有不均匀性的线性约束  $C$  。
 * 在 @ref
 * constraints 模块中对这种方法进行了详细解释。
 *
 *
 *
 * @note
 * 目前，这个函数对于分布式数据结构可能无法正确工作。
 * @relatesalso  LinearOperator
 * @ingroup constraints
 *
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
constrained_linear_operator(
  const AffineConstraints<typename Range::value_type> &constraints,
  const LinearOperator<Range, Domain, Payload> &       linop)
{
  const auto C    = distribute_constraints_linear_operator(constraints, linop);
  const auto Ct   = transpose_operator(C);
  const auto Id_c = project_to_constrained_linear_operator(constraints, linop);
  return Ct * linop * C + Id_c;
}


/**
 * 给定一个AffineConstraints对象 @p constraints, 一个LinearOperator @p
 * linop和一个右手边 @p right_hand_side,
 * ，这个函数创建一个PackagedOperation，存储以下计算结果。
 *
 * @code
 * Ct (right_hand_side
 *
 * - linop k)
 * @endcode
 * 与
 *
 * @code
 * C = distribute_constraints_linear_operator(constraints, linop);
 * Ct = transpose_operator(C);
 * @endcode
 *
 * 这个LinearOperator对象与constrained_right_hand_side()一起用于建立以下修改后的线性方程组。@f[
 * (C^T A C + Id_c) x = C^T (b
 *
 * - A\,k)
 * @f] 具有给定的（无约束的）系统矩阵  $A$  ，右手边  $b$  ，以及具有不均匀性的线性约束  $C$  。
 * 在 @ref
 * constraints 模块中对这种方法进行了详细解释。
 *
 *
 *
 * @note
 * 目前，这个函数对于分布式数据结构可能无法正确工作。
 * @relatesalso  LinearOperator
 * @ingroup constraints
 *
 *
 */
template <typename Range, typename Domain, typename Payload>
PackagedOperation<Range>
constrained_right_hand_side(
  const AffineConstraints<typename Range::value_type> &constraints,
  const LinearOperator<Range, Domain, Payload> &       linop,
  const Range &                                        right_hand_side)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = linop.reinit_range_vector;

  return_comp.apply_add = [&constraints, &linop, &right_hand_side](Range &v) {
    const auto C  = distribute_constraints_linear_operator(constraints, linop);
    const auto Ct = transpose_operator(C);

    GrowingVectorMemory<Domain>            vector_memory;
    typename VectorMemory<Domain>::Pointer k(vector_memory);
    linop.reinit_domain_vector(*k,  /*bool fast=*/ false);
    constraints.distribute(*k);

    v += Ct * (right_hand_side - linop * *k);
  };

  return_comp.apply = [apply_add = return_comp.apply_add](Range &v) {
    v = 0.;
    apply_add(v);
  };

  return return_comp;
}

//@}

DEAL_II_NAMESPACE_CLOSE

#endif


