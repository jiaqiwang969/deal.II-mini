//include/deal.II-translator/lac/linear_operator_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2021 by the deal.II authors
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

#ifndef dealii_linear_operator_h
#define dealii_linear_operator_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/vector_memory.h>

#include <array>
#include <functional>
#include <type_traits>

DEAL_II_NAMESPACE_OPEN

// Forward declarations:
#ifndef DOXYGEN
namespace internal
{
  namespace LinearOperatorImplementation
  {
    class EmptyPayload;
  }
} // namespace internal

template <typename Number>
class Vector;

class PreconditionIdentity;

template <typename Range  = Vector<double>,
          typename Domain = Range,
          typename Payload =
            internal::LinearOperatorImplementation::EmptyPayload>
class LinearOperator;
#endif

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload,
  typename OperatorExemplar,
  typename Matrix>
LinearOperator<Range, Domain, Payload>
linear_operator(const OperatorExemplar &, const Matrix &);

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload,
  typename Matrix>
LinearOperator<Range, Domain, Payload>
linear_operator(const Matrix &);

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload>
LinearOperator<Range, Domain, Payload>
null_operator(const LinearOperator<Range, Domain, Payload> &);

template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
identity_operator(const LinearOperator<Range, Domain, Payload> &);


/**
 * 一个用于存储线性算子的抽象概念的类。 该类本质上由
 * <code>std::function</code> 对象组成，这些对象通过实现抽象的
 * @p Matrix 接口来存储如何应用线性运算符的知识。
 *
 * @code
 * std::function<void(Range &, const Domain &)> vmult;
 * std::function<void(Range &, const Domain &)> vmult_add;
 * std::function<void(Domain &, const Range &)> Tvmult;
 * std::function<void(Domain &, const Range &)> Tvmult_add;
 * @endcode
 *
 * 但是，与通常的矩阵对象不同，线性运算符的域和范围也被绑定到类型级别的LinearOperator类上。因为这个原因，
 * <code>LinearOperator <Range, Domain></code> 有两个额外的函数对象
 *
 * @code
 * std::function<void(Range &, bool)> reinit_range_vector;
 * std::function<void(Domain &, bool)> reinit_domain_vector;
 * @endcode
 * 储存了如何初始化（调整大小+内部数据结构） @p Range 和
 * @p Domain 空间的任意矢量的知识。
 * 这个类的主要目的是为复杂的矩阵-向量操作提供语法糖，使用户不必手工创建、设置和处理中间存储位置。
 * 作为一个例子，考虑操作 $(A+k\,B)\,C$ ，其中 $A$ 、 $B$ 和
 * $C$ 表示（可能不同的）矩阵。为了构造一个LinearOperator
 * <code>op</code> 来存储这个操作的知识，我们可以这样写。
 *
 *
 * @code
 * #include <deal.II/lac/linear_operator_tools.h>
 *
 * dealii::SparseMatrix<double> A, B, C;
 * const double k = ...;
 *
 * // Setup and assembly of matrices
 *
 * const auto op_a = linear_operator(A);
 * const auto op_b = linear_operator(B);
 * const auto op_c = linear_operator(C);
 *
 * const auto op = (op_a + k op_b) op_c;
 * @endcode
 *
 *
 *
 * @note  这个类大量使用了 <code>std::function</code>
 * 对象和lambda函数。这种灵活性带来了运行时间的惩罚。只用这个对象来封装大中型的矩阵对象（作为经验法则，稀疏矩阵的大小为
 * $1000\times1000$ ，或更大）。
 *
 *
 * @note
 * 为了将Trilinos或PETSc稀疏矩阵和预处理程序与LinearOperator类结合使用，有必要通过一个额外的Payload来扩展LinearOperator类的功能。
 * 例如。代表矩阵求逆的LinearOperator实例通常需要调用一些线性求解器。这些求解器可能没有与LinearOperator的接口（例如，它可能代表一个复合操作）。因此，
 * TrilinosWrappers::internal::LinearOperatorImplementation::TrilinosPayload
 * 为LinearOperator提供了一个接口扩展，以便它可以被传递给求解器并被求解器使用，就像它是一个Trilinos算子一样。这意味着特定Trilinos算子的所有必要功能都在Payload类中被重载。这包括算子-向量乘法和反算子-向量乘法，其中算子可以是
 * TrilinosWrappers::SparseMatrix 或 TrilinosWrappers::PreconditionBase
 * ，向量是本地Trilinos向量。
 * 另一种情况是，当构建复合运算时（通过运算符重载），有效载荷为LinearOperator类提供了重要补充。在这种情况下，又有必要提供一个接口，以产生与Trilinos求解器使用的Trilinos算子兼容的这种复合操作的结果。
 *
 *
 * @note
 * LinearOperator的许多用例会导致中间表达式需要一个PackagedOperation。为了一次性包含所有必要的头文件，可以考虑使用
 *
 * @code
 * #include <deal.II/lac/linear_operator_tools.h>
 * @endcode
 *
 * 为了使用完整的LinearOperator和PackagedOperation
 *
 *
 * @note
 * 为了确保提供正确的有效载荷，在各自的TrilinosWrappers（以及未来的PETScWrappers）命名空间中提供了线性操作符的封装函数。
 *
 *
 * @note   step-20
 * 教程中有一个LinearOperator类的详细使用例子。
 *
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Range, typename Domain, typename Payload>
class LinearOperator : public Payload
{
public:
  /**
   * 创建一个空的LinearOperator对象。
   * 当一个有效载荷被传递给这个构造函数时，产生的运算器会被构造成一个功能性的有效载荷。
   * 在这两种情况下，这个构造函数产生的对象实际上不能用于任何线性运算符的操作，并且在调用时将抛出一个异常。
   *
   */
  LinearOperator(const Payload &payload = Payload())
    : Payload(payload)
    , is_null_operator(false)
  {
    vmult = [](Range &, const Domain &) {
      Assert(false,
             ExcMessage("Uninitialized LinearOperator<Range, "
                        "Domain>::vmult called"));
    };

    vmult_add = [](Range &, const Domain &) {
      Assert(false,
             ExcMessage("Uninitialized LinearOperator<Range, "
                        "Domain>::vmult_add called"));
    };

    Tvmult = [](Domain &, const Range &) {
      Assert(false,
             ExcMessage("Uninitialized LinearOperator<Range, "
                        "Domain>::Tvmult called"));
    };

    Tvmult_add = [](Domain &, const Range &) {
      Assert(false,
             ExcMessage("Uninitialized LinearOperator<Range, "
                        "Domain>::Tvmult_add called"));
    };

    reinit_range_vector = [](Range &, bool) {
      Assert(false,
             ExcMessage("Uninitialized LinearOperator<Range, "
                        "Domain>::reinit_range_vector method called"));
    };

    reinit_domain_vector = [](Domain &, bool) {
      Assert(false,
             ExcMessage("Uninitialized LinearOperator<Range, "
                        "Domain>::reinit_domain_vector method called"));
    };
  }

  /**
   * 默认的复制构造函数。
   *
   */
  LinearOperator(const LinearOperator<Range, Domain, Payload> &) = default;

  /**
   * 模板化的复制构造函数，从定义了转换函数
   * <code>linear_operator</code> 的对象 @p op
   * 中创建一个LinearOperator对象。
   *
   */
  template <typename Op,
            typename = typename std::enable_if<
              !std::is_base_of<LinearOperator<Range, Domain, Payload>,
                               Op>::value>::type>
  LinearOperator(const Op &op)
  {
    *this = linear_operator<Range, Domain, Payload, Op>(op);
  }

  /**
   * 默认的复制赋值运算符。
   *
   */
  LinearOperator<Range, Domain, Payload> &
  operator=(const LinearOperator<Range, Domain, Payload> &) = default;

  /**
   * 定义了转换函数 <code>linear_operator</code> 的对象 @p op
   * 的模板化的复制赋值运算符。
   *
   */
  template <typename Op,
            typename = typename std::enable_if<
              !std::is_base_of<LinearOperator<Range, Domain, Payload>,
                               Op>::value>::type>
  LinearOperator<Range, Domain, Payload> &
  operator=(const Op &op)
  {
    *this = linear_operator<Range, Domain, Payload, Op>(op);
    return *this;
  }

  /**
   * 将LinearOperator对象应用于 @p Domain 空间的向量u，得到 @p
   * Range 空间的向量v。
   *
   */
  std::function<void(Range &v, const Domain &u)> vmult;

  /**
   * 将LinearOperator对象应用于 @p Domain
   * 空间的向量u。其结果被添加到向量v中。
   *
   */
  std::function<void(Range &v, const Domain &u)> vmult_add;

  /**
   * 对 @p Range 空间的向量u应用转置LinearOperator对象，得到 @p
   * Domain 空间的向量v。
   *
   */
  std::function<void(Domain &v, const Range &u)> Tvmult;

  /**
   * 对 @p Range 空间的向量 @p u
   * 应用转置LinearOperator对象。结果被添加到向量 @p v. 。
   *
   */
  std::function<void(Domain &v, const Range &u)> Tvmult_add;

  /**
   * 初始化Range空间的一个向量v，可以直接作为vmult应用中的目标参数使用。与向量类的reinit函数类似，布尔值决定是否进行快速初始化，也就是说，如果它被设置为false，向量的内容就被设置为0。
   *
   */
  std::function<void(Range &v, bool omit_zeroing_entries)> reinit_range_vector;

  /**
   * 初始化域空间的一个向量，以便在vmult的应用中可以直接作为源参数使用。与向量类的reinit函数类似，布尔值决定是否进行快速初始化，也就是说，如果它被设置为false，向量的内容就被设置为0。
   *
   */
  std::function<void(Domain &v, bool omit_zeroing_entries)>
    reinit_domain_vector;

  /**
   * @name  原地向量空间操作
   *
   */
  //@{

  /**
   * 用LinearOperator @p second_op 进行加法运算，同 @p Domain 和 @p
   * Range. 。
   *
   */
  LinearOperator<Range, Domain, Payload> &
  operator+=(const LinearOperator<Range, Domain, Payload> &second_op)
  {
    *this = *this + second_op;
    return *this;
  }

  /**
   * 用LinearOperator  @p second_op 和 @p Range. 做减法。
   *
   */
  LinearOperator<Range, Domain, Payload> &
  operator-=(const LinearOperator<Range, Domain, Payload> &second_op)
  {
    *this = *this - second_op;
    return *this;
  }

  /**
   * 与 @p Domain 空间的内形态 @p second_op
   * 的LinearOperator的构成。
   *
   */
  LinearOperator<Range, Domain, Payload> &
  operator*=(const LinearOperator<Domain, Domain, Payload> &second_op)
  {
    *this = *this * second_op;
    return *this;
  }

  /**
   * LinearOperator与 @p number 从右边开始的标量乘法。
   *
   */
  LinearOperator<Range, Domain, Payload>
  operator*=(typename Domain::value_type number)
  {
    *this = *this * number;
    return *this;
  }

  /**
   * 这个bool用于确定线性运算符是否是空运算符。在这种情况下，该类能够优化一些操作，如乘法或加法。
   *
   */
  bool is_null_operator;

  //@}
};


/**
 * @name  矢量空间操作
 *
 *
 */
//@{

/**
 * @relatesalso LinearOperator 两个线性运算符 @p first_op 和 @p
 * second_op 的相加，由 $(\mathrm{first\_op}+\mathrm{second\_op})x
 * \dealcoloneq \mathrm{first\_op}(x) + \mathrm{second\_op}(x)$ 给出。
 *
 * @ingroup LAOperators
 *
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
operator+(const LinearOperator<Range, Domain, Payload> &first_op,
          const LinearOperator<Range, Domain, Payload> &second_op)
{
  if (first_op.is_null_operator)
    {
      return second_op;
    }
  else if (second_op.is_null_operator)
    {
      return first_op;
    }
  else
    {
      LinearOperator<Range, Domain, Payload> return_op{
        static_cast<const Payload &>(first_op) +
        static_cast<const Payload &>(second_op)};

      return_op.reinit_range_vector  = first_op.reinit_range_vector;
      return_op.reinit_domain_vector = first_op.reinit_domain_vector;

      // ensure to have valid computation objects by catching first_op and
      // second_op by value

      return_op.vmult = [first_op, second_op](Range &v, const Domain &u) {
        first_op.vmult(v, u);
        second_op.vmult_add(v, u);
      };

      return_op.vmult_add = [first_op, second_op](Range &v, const Domain &u) {
        first_op.vmult_add(v, u);
        second_op.vmult_add(v, u);
      };

      return_op.Tvmult = [first_op, second_op](Domain &v, const Range &u) {
        second_op.Tvmult(v, u);
        first_op.Tvmult_add(v, u);
      };

      return_op.Tvmult_add = [first_op, second_op](Domain &v, const Range &u) {
        second_op.Tvmult_add(v, u);
        first_op.Tvmult_add(v, u);
      };

      return return_op;
    }
}


/**
 * @relatesalso LinearOperator 两个线性运算符 @p first_op 和 @p
 * second_op 的减法，由 $(\mathrm{first\_op}-\mathrm{second\_op})x
 * \dealcoloneq \mathrm{first\_op}(x)
 *
 *
 *
 * - \mathrm{second\_op}(x)$ 给出。
 *
 * @ingroup LAOperators
 *
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
operator-(const LinearOperator<Range, Domain, Payload> &first_op,
          const LinearOperator<Range, Domain, Payload> &second_op)
{
  if (first_op.is_null_operator)
    {
      return -1. * second_op;
    }
  else if (second_op.is_null_operator)
    {
      return first_op;
    }
  else
    {
      // implement with addition and scalar multiplication
      return first_op + (-1. * second_op);
    }
}


/**
 * @relatesalso  LinearOperator ScalarOperator对象 @p op 与 @p number
 * 的标量乘法，从左边开始。 @p Domain 和 @p Range
 * 类型必须实现以下 <code>operator*=</code>
 * 成员函数，接受适当的标量Number类型进行重新缩放。
 *
 *
 * @code
 * Domain & operator=(Domain::value_type);
 * Range & operator=(Range::value_type);
 * @endcode
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
operator*(typename Range::value_type                    number,
          const LinearOperator<Range, Domain, Payload> &op)
{
  static_assert(
    std::is_convertible<typename Range::value_type,
                        typename Domain::value_type>::value,
    "Range and Domain must have implicitly convertible 'value_type's");

  if (op.is_null_operator)
    {
      return op;
    }
  else if (number == 0.)
    {
      return null_operator(op);
    }
  else
    {
      LinearOperator<Range, Domain, Payload> return_op = op;

      // ensure to have valid computation objects by catching number and op by
      // value

      return_op.vmult = [number, op](Range &v, const Domain &u) {
        op.vmult(v, u);
        v *= number;
      };

      return_op.vmult_add = [number, op](Range &v, const Domain &u) {
        v /= number;
        op.vmult_add(v, u);
        v *= number;
      };

      return_op.Tvmult = [number, op](Domain &v, const Range &u) {
        op.Tvmult(v, u);
        v *= number;
      };

      return_op.Tvmult_add = [number, op](Domain &v, const Range &u) {
        v /= number;
        op.Tvmult_add(v, u);
        v *= number;
      };

      return return_op;
    }
}


/**
 * @relatesalso  LinearOperator
 * 从右边开始对一个ScalarOperator对象进行标量乘法。 @p Domain
 * 和 @p Range 类型必须实现以下 <code>operator*=</code>
 * 成员函数，用于重新缩放。
 *
 *
 * @code
 * Domain & operator=(Domain::value_type);
 * Range & operator=(Range::value_type);
 * @endcode
 *
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
operator*(const LinearOperator<Range, Domain, Payload> &op,
          typename Domain::value_type                   number)
{
  static_assert(
    std::is_convertible<typename Range::value_type,
                        typename Domain::value_type>::value,
    "Range and Domain must have implicitly convertible 'value_type's");

  return number * op;
}

//@}


/**
 * @name  LinearOperator的组成和操作
 *
 *
 */
//@{

/**
 * @relatesalso LinearOperator 两个线性运算符 @p first_op 和 @p
 * second_op 的组合，由 $(\mathrm{first\_op}*\mathrm{second\_op})x
 * \dealcoloneq \mathrm{first\_op}(\mathrm{second\_op}(x))$ 给出。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range,
          typename Intermediate,
          typename Domain,
          typename Payload>
LinearOperator<Range, Domain, Payload>
operator*(const LinearOperator<Range, Intermediate, Payload> & first_op,
          const LinearOperator<Intermediate, Domain, Payload> &second_op)
{
  if (first_op.is_null_operator || second_op.is_null_operator)
    {
      LinearOperator<Range, Domain, Payload> return_op;
      return_op.reinit_domain_vector = second_op.reinit_domain_vector;
      return_op.reinit_range_vector  = first_op.reinit_range_vector;
      return null_operator(return_op);
    }
  else
    {
      LinearOperator<Range, Domain, Payload> return_op{
        static_cast<const Payload &>(first_op) *
        static_cast<const Payload &>(second_op)};

      return_op.reinit_domain_vector = second_op.reinit_domain_vector;
      return_op.reinit_range_vector  = first_op.reinit_range_vector;

      // ensure to have valid computation objects by catching first_op and
      // second_op by value

      return_op.vmult = [first_op, second_op](Range &v, const Domain &u) {
        GrowingVectorMemory<Intermediate> vector_memory;

        typename VectorMemory<Intermediate>::Pointer i(vector_memory);
        second_op.reinit_range_vector(*i,  /*bool omit_zeroing_entries =*/ true);
        second_op.vmult(*i, u);
        first_op.vmult(v, *i);
      };

      return_op.vmult_add = [first_op, second_op](Range &v, const Domain &u) {
        GrowingVectorMemory<Intermediate> vector_memory;

        typename VectorMemory<Intermediate>::Pointer i(vector_memory);
        second_op.reinit_range_vector(*i,  /*bool omit_zeroing_entries =*/ true);
        second_op.vmult(*i, u);
        first_op.vmult_add(v, *i);
      };

      return_op.Tvmult = [first_op, second_op](Domain &v, const Range &u) {
        GrowingVectorMemory<Intermediate> vector_memory;

        typename VectorMemory<Intermediate>::Pointer i(vector_memory);
        first_op.reinit_domain_vector(*i,  /*bool omit_zeroing_entries =*/ true);
        first_op.Tvmult(*i, u);
        second_op.Tvmult(v, *i);
      };

      return_op.Tvmult_add = [first_op, second_op](Domain &v, const Range &u) {
        GrowingVectorMemory<Intermediate> vector_memory;

        typename VectorMemory<Intermediate>::Pointer i(vector_memory);
        first_op.reinit_domain_vector(*i,  /*bool omit_zeroing_entries =*/ true);
        first_op.Tvmult(*i, u);
        second_op.Tvmult_add(v, *i);
      };

      return return_op;
    }
}


/**
 * @relatesalso  LinearOperator 返回 @p op. 的转置线性运算。
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Domain, Range, Payload>
transpose_operator(const LinearOperator<Range, Domain, Payload> &op)
{
  LinearOperator<Domain, Range, Payload> return_op{op.transpose_payload()};

  return_op.reinit_range_vector  = op.reinit_domain_vector;
  return_op.reinit_domain_vector = op.reinit_range_vector;

  return_op.vmult      = op.Tvmult;
  return_op.vmult_add  = op.Tvmult_add;
  return_op.Tvmult     = op.vmult;
  return_op.Tvmult_add = op.vmult_add;

  return return_op;
}


/**
 * @relatesalso  LinearOperator
 * 返回一个代表LinearOperator的逆运算的对象  @p op.  。
 * 该函数需要引用 @p solver 和 @p preconditioner
 * 一个迭代求解器和一个预处理器，这些都是LinearOperator对象的
 * <code>vmult</code> and <code>Tvmult</code> 实现中使用的。 创建的
 * LinearOperator 对象存储了对  @p solver  和  @p preconditioner.
 * 的引用。因此，这两个对象必须在 LinearOperator
 * 对象的整个生命周期内保持有效引用。 @p solver
 * 对象的内部数据结构将在调用  <code>vmult</code> or
 * <code>Tvmult</code>  时被修改。
 *
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Payload,
          typename Solver,
          typename Preconditioner,
          typename Range  = typename Solver::vector_type,
          typename Domain = Range>
LinearOperator<Domain, Range, Payload>
inverse_operator(const LinearOperator<Range, Domain, Payload> &op,
                 Solver &                                      solver,
                 const Preconditioner &                        preconditioner)
{
  LinearOperator<Domain, Range, Payload> return_op{
    op.inverse_payload(solver, preconditioner)};

  return_op.reinit_range_vector  = op.reinit_domain_vector;
  return_op.reinit_domain_vector = op.reinit_range_vector;

  return_op.vmult = [op, &solver, &preconditioner](Range &v, const Domain &u) {
    op.reinit_range_vector(v,  /*bool omit_zeroing_entries =*/ false);
    solver.solve(op, v, u, preconditioner);
  };

  return_op.vmult_add = [op, &solver, &preconditioner](Range &       v,
                                                       const Domain &u) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer v2(vector_memory);
    op.reinit_range_vector(*v2,  /*bool omit_zeroing_entries =*/ false);
    solver.solve(op, *v2, u, preconditioner);
    v += *v2;
  };

  return_op.Tvmult = [op, &solver, &preconditioner](Range &v, const Domain &u) {
    op.reinit_range_vector(v,  /*bool omit_zeroing_entries =*/ false);
    solver.solve(transpose_operator(op), v, u, preconditioner);
  };

  return_op.Tvmult_add = [op, &solver, &preconditioner](Range &       v,
                                                        const Domain &u) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer v2(vector_memory);
    op.reinit_range_vector(*v2,  /*bool omit_zeroing_entries =*/ false);
    solver.solve(transpose_operator(op), *v2, u, preconditioner);
    v += *v2;
  };

  return return_op;
}


/**
 * @relatesalso  LinearOperator 上述函数的变体，接受 LinearOperator
 * @p preconditioner  作为预处理参数。
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Payload,
          typename Solver,
          typename Range  = typename Solver::vector_type,
          typename Domain = Range>
LinearOperator<Domain, Range, Payload>
inverse_operator(const LinearOperator<Range, Domain, Payload> &op,
                 Solver &                                      solver,
                 const LinearOperator<Range, Domain, Payload> &preconditioner)
{
  LinearOperator<Domain, Range, Payload> return_op{
    op.inverse_payload(solver, preconditioner)};

  return_op.reinit_range_vector  = op.reinit_domain_vector;
  return_op.reinit_domain_vector = op.reinit_range_vector;

  return_op.vmult = [op, &solver, preconditioner](Range &v, const Domain &u) {
    op.reinit_range_vector(v,  /*bool omit_zeroing_entries =*/ false);
    solver.solve(op, v, u, preconditioner);
  };

  return_op.vmult_add = [op, &solver, preconditioner](Range &       v,
                                                      const Domain &u) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer v2(vector_memory);
    op.reinit_range_vector(*v2,  /*bool omit_zeroing_entries =*/ false);
    solver.solve(op, *v2, u, preconditioner);
    v += *v2;
  };

  return_op.Tvmult = [op, &solver, preconditioner](Range &v, const Domain &u) {
    op.reinit_range_vector(v,  /*bool omit_zeroing_entries =*/ false);
    solver.solve(transpose_operator(op), v, u, preconditioner);
  };

  return_op.Tvmult_add = [op, &solver, preconditioner](Range &       v,
                                                       const Domain &u) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer v2(vector_memory);
    op.reinit_range_vector(*v2,  /*bool omit_zeroing_entries =*/ false);
    solver.solve(transpose_operator(op), *v2, u, preconditioner);
    v += *v2;
  };

  return return_op;
}


/**
 * @relatesalso  LinearOperator
 * 上述函数的变体，没有预处理参数。在这种情况下， @p op
 * 参数的ident_operator()被用作预处理。这等同于使用PreconditionIdentity。
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Payload,
          typename Solver,
          typename Range  = typename Solver::vector_type,
          typename Domain = Range>
LinearOperator<Domain, Range, Payload>
inverse_operator(const LinearOperator<Range, Domain, Payload> &op,
                 Solver &                                      solver)
{
  return inverse_operator(op, solver, identity_operator(op));
}


/**
 * @relatesalso  LinearOperator
 * 上述函数的特殊重载，需要一个PreconditionIdentity参数。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Payload,
          typename Solver,
          typename Range  = typename Solver::vector_type,
          typename Domain = Range>
LinearOperator<Domain, Range, Payload>
inverse_operator(const LinearOperator<Range, Domain, Payload> &op,
                 Solver &                                      solver,
                 const PreconditionIdentity &)
{
  return inverse_operator(op, solver);
}

//@}


/**
 * @name  创建一个LinearOperator
 *
 *
 */
//@{

/**
 * @relatesalso  LinearOperator
 * 返回一个LinearOperator，它是向量空间的标识  @p Range.  。
 * 该函数需要一个 <code>std::function</code> 对象 @p reinit_vector
 * 作为参数，以初始化LinearOperator对象的
 * <code>reinit_range_vector</code> 和 <code>reinit_domain_vector</code>
 * 对象。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <
  typename Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload>
LinearOperator<Range, Range, Payload>
identity_operator(const std::function<void(Range &, bool)> &reinit_vector)
{
  LinearOperator<Range, Range, Payload> return_op{Payload()};

  return_op.reinit_range_vector  = reinit_vector;
  return_op.reinit_domain_vector = reinit_vector;

  return_op.vmult = [](Range &v, const Range &u) { v = u; };

  return_op.vmult_add = [](Range &v, const Range &u) { v += u; };

  return_op.Tvmult = [](Range &v, const Range &u) { v = u; };

  return_op.Tvmult_add = [](Range &v, const Range &u) { v += u; };

  return return_op;
}


/**
 * @relatesalso  LinearOperator
 * 返回一个LinearOperator，它是向量空间的标识  @p Range.  。
 * 该函数接收一个LinearOperator  @p op
 * 并使用其范围初始化器来创建一个身份算子。与上面的函数不同，这个函数还确保底层的Payload与输入的Payload相匹配
 * @p op.
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
identity_operator(const LinearOperator<Range, Domain, Payload> &op)
{
  auto return_op = identity_operator<Range, Payload>(op.reinit_range_vector);
  static_cast<Payload &>(return_op) = op.identity_payload();

  return return_op;
}


/**
 * @relatesalso  LinearOperator 返回LinearOperator  @p op,
 * 的空变体，即具有优化的 LinearOperator::vmult,
 * LinearOperator::vmult_add, 等函数，并将
 * LinearOperator::is_null_operator 设置为真。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
null_operator(const LinearOperator<Range, Domain, Payload> &op)
{
  LinearOperator<Range, Domain, Payload> return_op{op.null_payload()};

  return_op.is_null_operator = true;

  return_op.reinit_range_vector  = op.reinit_range_vector;
  return_op.reinit_domain_vector = op.reinit_domain_vector;

  return_op.vmult = [](Range &v, const Domain &) { v = 0.; };

  return_op.vmult_add = [](Range &, const Domain &) {};

  return_op.Tvmult = [](Domain &v, const Range &) { v = 0.; };

  return_op.Tvmult_add = [](Domain &, const Range &) {};

  return return_op;
}


/**
 * @relatesalso  LinearOperator
 * 返回一个LinearOperator，作为一个均值过滤器。这个矩阵的vmult()函数减去向量的均值。
 * 该函数需要一个 <code>std::function</code> 对象 @p reinit_vector
 * 作为参数来初始化LinearOperator对象的
 * <code>reinit_range_vector</code> 和 <code>reinit_domain_vector</code>
 * 对象。
 *
 *
 * @ingroup LAOperators
 *
 */
template <
  typename Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload>
LinearOperator<Range, Range, Payload>
mean_value_filter(const std::function<void(Range &, bool)> &reinit_vector)
{
  LinearOperator<Range, Range, Payload> return_op{Payload()};

  return_op.reinit_range_vector  = reinit_vector;
  return_op.reinit_domain_vector = reinit_vector;

  return_op.vmult = [](Range &v, const Range &u) {
    const auto mean = u.mean_value();

    v = u;
    v.add(-mean);
  };

  return_op.vmult_add = [](Range &v, const Range &u) {
    const auto mean = u.mean_value();

    v += u;
    v.add(-mean);
  };

  return_op.Tvmult     = return_op.vmult_add;
  return_op.Tvmult_add = return_op.vmult_add;

  return return_op;
}


/**
 * @relatesalso  LinearOperator
 * 返回一个LinearOperator，作为一个均值过滤器。这个矩阵的vmult()函数减去向量的均值。
 * 该函数接收一个LinearOperator  @p op
 * 并使用其范围初始化器来创建一个均值过滤器运算器。该函数还确保底层的Payload与输入的Payload相匹配
 * @p op.  。
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Range, typename Domain, typename Payload>
LinearOperator<Range, Domain, Payload>
mean_value_filter(const LinearOperator<Range, Domain, Payload> &op)
{
  auto return_op = mean_value_filter<Range, Payload>(op.reinit_range_vector);
  static_cast<Payload &>(return_op) = op.identity_payload();

  return return_op;
}


namespace internal
{
  namespace LinearOperatorImplementation
  {
    /**
     * 一个辅助类，负责初始化一个向量，使其可以直接作为目标参数，或者在矩阵的vmult应用中的源参数。
     * 这个类的通用版本只是分别调用
     * <code>Vector::reinit()</code> 和 <code>Matrix::m()</code> or
     * <code>Matrix::n()</code> 的结果。
     * 这个类专门用于更复杂的数据结构，如
     * TrilinosWrappers::MPI::Vector, 等。
     *
     */
    template <typename Vector>
    class ReinitHelper
    {
    public:
      /**
       * 初始化Range空间的一个向量v，以便在vmult的应用中直接作为目标参数使用。与向量类的reinit函数类似，布尔值决定是否进行快速初始化，即如果它被设置为false，向量的内容就被设置为0。这个类的通用版本只是调用
       * <code>Vector::reinit()</code> ，其结果是
       * <code>Matrix::m()</code>  。
       *
       */
      template <typename Matrix>
      static void
      reinit_range_vector(const Matrix &matrix,
                          Vector &      v,
                          bool          omit_zeroing_entries)
      {
        v.reinit(matrix.m(), omit_zeroing_entries);
      }

      /**
       * 初始化域空间的一个向量，在vmult的应用中可直接作为源参数使用。与向量类的reinit函数类似，布尔值决定是否进行快速初始化，也就是说，如果它被设置为false，向量的内容就被设置为0。
       *
       */
      template <typename Matrix>
      static void
      reinit_domain_vector(const Matrix &matrix,
                           Vector &      v,
                           bool          omit_zeroing_entries)
      {
        v.reinit(matrix.n(), omit_zeroing_entries);
      }
    };


    /**
     * 一个假的类，用于不需要任何扩展的LinearOperators，以方便对矩阵进行操作。
     * 这是通常与deal.II的原生SparseMatrix相关的Payload类。为了使用Trilinos和PETSc稀疏矩阵类，有必要用它们相关的Payload初始化LinearOperator。
     * @ingroup LAOperators
     *
     */
    class EmptyPayload
    {
    public:
      /**
       * 默认构造函数
       * 由于这个类不做任何特别的事情，也不需要特别的配置，所以我们只有一个通用构造函数，可以在任何条件下调用。
       *
       */
      template <typename... Args>
      EmptyPayload(const Args &...)
      {}


      /**
       * 返回一个为身份操作配置的有效载荷
       *
       */
      EmptyPayload
      identity_payload() const
      {
        return *this;
      }


      /**
       * 返回一个为空操作而配置的有效载荷
       *
       */
      EmptyPayload
      null_payload() const
      {
        return *this;
      }


      /**
       * 返回一个为转置操作配置的有效载荷
       *
       */
      EmptyPayload
      transpose_payload() const
      {
        return *this;
      }


      /**
       * 返回一个为反转操作配置的有效载荷
       *
       */
      template <typename Solver, typename Preconditioner>
      EmptyPayload
      inverse_payload(Solver &, const Preconditioner &) const
      {
        return *this;
      }
    };

    /**
     * 返回一个配置为支持两个LinearOperators相加的有效载荷的操作符
     *
     */
    inline EmptyPayload
    operator+(const EmptyPayload &, const EmptyPayload &)
    {
      return {};
    }

    /**
     * 返回被配置为支持两个LinearOperators的乘法的有效载荷的操作符。
     *
     */
    inline EmptyPayload operator*(const EmptyPayload &, const EmptyPayload &)
    {
      return {};
    }



    // A trait class that determines whether type T provides public
    // (templated or non-templated) vmult_add member functions
    template <typename Range, typename Domain, typename T>
    class has_vmult_add_and_Tvmult_add
    {
      template <typename C>
      static std::false_type
      test(...);

      template <typename C>
      static auto
      test(Range *r, Domain *d)
        -> decltype(std::declval<C>().vmult_add(*r, *d),
                    std::declval<C>().Tvmult_add(*d, *r),
                    std::true_type());

    public:
      // type is std::true_type if Matrix provides vmult_add and Tvmult_add,
      // otherwise it is std::false_type

      using type = decltype(test<T>(nullptr, nullptr));
    };


    // A helper function to apply a given vmult, or Tvmult to a vector with
    // intermediate storage
    template <typename Function, typename Range, typename Domain>
    void
    apply_with_intermediate_storage(Function      function,
                                    Range &       v,
                                    const Domain &u,
                                    bool          add)
    {
      GrowingVectorMemory<Range> vector_memory;

      typename VectorMemory<Range>::Pointer i(vector_memory);
      i->reinit(v,  /*bool omit_zeroing_entries =*/ true);

      function(*i, u);

      if (add)
        v += *i;
      else
        v = *i;
    }


    // A helper class to add a reduced matrix interface to a LinearOperator
    // (typically provided by Preconditioner classes)
    template <typename Range, typename Domain, typename Payload>
    class MatrixInterfaceWithoutVmultAdd
    {
    public:
      template <typename Matrix>
      void
      operator()(LinearOperator<Range, Domain, Payload> &op,
                 const Matrix &                          matrix)
      {
        op.vmult = [&matrix](Range &v, const Domain &u) {
          if (PointerComparison::equal(&v, &u))
            {
              // If v and u are the same memory location use intermediate
              // storage
              apply_with_intermediate_storage(
                [&matrix](Range &b, const Domain &a) { matrix.vmult(b, a); },
                v,
                u,
                 /*bool add =*/ false);
            }
          else
            {
              matrix.vmult(v, u);
            }
        };

        op.vmult_add = [&matrix](Range &v, const Domain &u) {
          // use intermediate storage to implement vmult_add with vmult
          apply_with_intermediate_storage(
            [&matrix](Range &b, const Domain &a) { matrix.vmult(b, a); },
            v,
            u,
             /*bool add =*/ true);
        };

        op.Tvmult = [&matrix](Domain &v, const Range &u) {
          if (PointerComparison::equal(&v, &u))
            {
              // If v and u are the same memory location use intermediate
              // storage
              apply_with_intermediate_storage(
                [&matrix](Domain &b, const Range &a) { matrix.Tvmult(b, a); },
                v,
                u,
                 /*bool add =*/ false);
            }
          else
            {
              matrix.Tvmult(v, u);
            }
        };

        op.Tvmult_add = [&matrix](Domain &v, const Range &u) {
          // use intermediate storage to implement Tvmult_add with Tvmult
          apply_with_intermediate_storage(
            [&matrix](Domain &b, const Range &a) { matrix.Tvmult(b, a); },
            v,
            u,
             /*bool add =*/ true);
        };
      }
    };


    // A helper class to add the full matrix interface to a LinearOperator
    template <typename Range, typename Domain, typename Payload>
    class MatrixInterfaceWithVmultAdd
    {
    public:
      template <typename Matrix>
      void
      operator()(LinearOperator<Range, Domain, Payload> &op,
                 const Matrix &                          matrix)
      {
        // As above ...

        MatrixInterfaceWithoutVmultAdd<Range, Domain, Payload>().operator()(
          op, matrix);

        // ... but add native vmult_add and Tvmult_add variants:

        op.vmult_add = [&matrix](Range &v, const Domain &u) {
          if (PointerComparison::equal(&v, &u))
            {
              apply_with_intermediate_storage(
                [&matrix](Range &b, const Domain &a) { matrix.vmult(b, a); },
                v,
                u,
                 /*bool add =*/ true);
            }
          else
            {
              matrix.vmult_add(v, u);
            }
        };

        op.Tvmult_add = [&matrix](Domain &v, const Range &u) {
          if (PointerComparison::equal(&v, &u))
            {
              apply_with_intermediate_storage(
                [&matrix](Domain &b, const Range &a) { matrix.Tvmult(b, a); },
                v,
                u,
                 /*bool add =*/ true);
            }
          else
            {
              matrix.Tvmult_add(v, u);
            }
        };
      }
    };
  } // namespace LinearOperatorImplementation
} // namespace internal


/**
 * @relatesalso  LinearOperator 一个将作用于兼容矢量类型的通用
 * @p matrix
 * 对象封装为LinearOperator的函数。被创建的LinearOperator对象存储了一个对矩阵对象的引用。因此，
 * @p matrix
 * 必须在LinearOperator对象的整个生命周期内保持有效引用。
 * 在创建LinearOperator对象之后，对 @p matrix
 * 所做的所有更改都会被操作者对象所反映。例如，首先创建一个LinearOperator，然后调整大小，以后再重新组合矩阵，这是一个有效的程序。
 * 有关的矩阵类必须提供以下最小接口。
 *
 *
 * @code
 * class Matrix
 * {
 * public:
 * // (type specific) information how to create a Range and Domain vector
 * // with appropriate size and internal layout
 *
 * // Application of matrix to vector src, writes the result into dst.
 * vmult(Range &dst, const Domain &src);
 *
 * // Application of the transpose of matrix to vector src, writes the
 * // result into dst. (Depending on the usage of the linear operator
 * // class this can be a dummy implementation throwing an error.)
 * Tvmult(Range &dst, const Domain &src);
 * };
 * @endcode
 *
 * 如果有以下（可选）接口，则使用该接口。
 *
 *
 * @code
 * class Matrix
 * {
 * public:
 * // Application of matrix to vector src, adds the result to dst.
 * vmult_add(Range &dst, const Domain &src);
 *
 * // Application of the transpose of matrix to vector src, adds the
 * // result to dst.
 * Tvmult_add(Range &dst, const Domain &src);
 * };
 * @endcode
 *
 * 如果矩阵不提供 <code>vmult_add</code> 和 <code>Tvmult_add</code>
 * ，则以 <code>vmult</code> and <code>Tvmult</code>
 * 的方式实现（需要中间存储）。
 *
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Range, typename Domain, typename Payload, typename Matrix>
LinearOperator<Range, Domain, Payload>
linear_operator(const Matrix &matrix)
{
  // implement with the more generic variant below...
  return linear_operator<Range, Domain, Payload, Matrix, Matrix>(matrix,
                                                                 matrix);
}


/**
 * @relatesalso  LinearOperator
 * 上述函数的变体，需要一个操作者对象 @p
 * operator_exemplar作为额外参考。这个对象被用来填充
 * reinit_domain_vector 和 reinit_range_vector 函数对象。引用 @p
 * matrix 用于构造vmult、Tvmult等。
 * 这个变体可以，例如，用于封装预处理程序（通常不暴露任何关于底层矩阵的信息）。
 *
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Range,
          typename Domain,
          typename Payload,
          typename OperatorExemplar,
          typename Matrix>
LinearOperator<Range, Domain, Payload>
linear_operator(const OperatorExemplar &operator_exemplar, const Matrix &matrix)
{
  using namespace internal::LinearOperatorImplementation;
  // Initialize the payload based on the input exemplar matrix
  LinearOperator<Range, Domain, Payload> return_op{
    Payload(operator_exemplar, matrix)};

  // Always store a reference to matrix and operator_exemplar in the lambda
  // functions. This ensures that a modification of the matrix after the
  // creation of a LinearOperator wrapper is respected - further a matrix
  // or an operator_exemplar cannot usually be copied...

  return_op.reinit_range_vector =
    [&operator_exemplar](Range &v, bool omit_zeroing_entries) {
      internal::LinearOperatorImplementation::ReinitHelper<
        Range>::reinit_range_vector(operator_exemplar, v, omit_zeroing_entries);
    };

  return_op.reinit_domain_vector = [&operator_exemplar](
                                     Domain &v, bool omit_zeroing_entries) {
    internal::LinearOperatorImplementation::ReinitHelper<
      Domain>::reinit_domain_vector(operator_exemplar, v, omit_zeroing_entries);
  };

  typename std::conditional<
    has_vmult_add_and_Tvmult_add<Range, Domain, Matrix>::type::value,
    MatrixInterfaceWithVmultAdd<Range, Domain, Payload>,
    MatrixInterfaceWithoutVmultAdd<Range, Domain, Payload>>::type()
    .
    operator()(return_op, matrix);

  return return_op;
}



/**
 * @relatesalso  LinearOperator 上述函数的变体，将LinearOperator  @p
 * operator_exemplar作为一个额外的引用。reinit_domain_vector 和
 * reinit_range_vector 函数是从  @p operator_exemplar
 * 对象中复制的。 参考 @p matrix 用于构造vmult、Tvmult等。
 * 这个变体可以，例如，用于封装预处理程序（通常不暴露任何关于底层矩阵的信息）。
 *
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Range, typename Domain, typename Payload, typename Matrix>
LinearOperator<Range, Domain, Payload>
linear_operator(const LinearOperator<Range, Domain, Payload> &operator_exemplar,
                const Matrix &                                matrix)
{
  using namespace internal::LinearOperatorImplementation;
  // Initialize the payload based on the LinearOperator exemplar
  auto return_op = operator_exemplar;

  typename std::conditional<
    has_vmult_add_and_Tvmult_add<Range, Domain, Matrix>::type::value,
    MatrixInterfaceWithVmultAdd<Range, Domain, Payload>,
    MatrixInterfaceWithoutVmultAdd<Range, Domain, Payload>>::type()
    .
    operator()(return_op, matrix);

  return return_op;
}


//@}

#ifndef DOXYGEN

//
// Ensure that we never capture a reference to a temporary by accident.
// to avoid "stack use after free".
//

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload,
  typename OperatorExemplar,
  typename Matrix,
  typename =
    typename std::enable_if<!std::is_lvalue_reference<Matrix>::value>::type>
LinearOperator<Range, Domain, Payload>
linear_operator(const OperatorExemplar &, Matrix &&) = delete;

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload,
  typename OperatorExemplar,
  typename Matrix,
  typename = typename std::enable_if<
    !std::is_lvalue_reference<OperatorExemplar>::value>::type,
  typename = typename std::enable_if<
    !std::is_same<OperatorExemplar,
                  LinearOperator<Range, Domain, Payload>>::value>::type>
LinearOperator<Range, Domain, Payload>
linear_operator(OperatorExemplar &&, const Matrix &) = delete;

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload,
  typename OperatorExemplar,
  typename Matrix,
  typename =
    typename std::enable_if<!std::is_lvalue_reference<Matrix>::value>::type,
  typename = typename std::enable_if<
    !std::is_lvalue_reference<OperatorExemplar>::value>::type,
  typename = typename std::enable_if<
    !std::is_same<OperatorExemplar,
                  LinearOperator<Range, Domain, Payload>>::value>::type>
LinearOperator<Range, Domain, Payload>
linear_operator(OperatorExemplar &&, Matrix &&) = delete;

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload,
  typename Matrix,
  typename =
    typename std::enable_if<!std::is_lvalue_reference<Matrix>::value>::type>
LinearOperator<Range, Domain, Payload>
linear_operator(const LinearOperator<Range, Domain, Payload> &,
                Matrix &&) = delete;

template <
  typename Range   = Vector<double>,
  typename Domain  = Range,
  typename Payload = internal::LinearOperatorImplementation::EmptyPayload,
  typename Matrix,
  typename =
    typename std::enable_if<!std::is_lvalue_reference<Matrix>::value>::type>
LinearOperator<Range, Domain, Payload>
linear_operator(Matrix &&) = delete;

template <typename Payload,
          typename Solver,
          typename Preconditioner,
          typename Range  = typename Solver::vector_type,
          typename Domain = Range,
          typename        = typename std::enable_if<
            !std::is_lvalue_reference<Preconditioner>::value>::type,
          typename = typename std::enable_if<
            !std::is_same<Preconditioner, PreconditionIdentity>::value>::type,
          typename = typename std::enable_if<
            !std::is_same<Preconditioner,
                          LinearOperator<Range, Domain, Payload>>::value>::type>
LinearOperator<Domain, Range, Payload>
inverse_operator(const LinearOperator<Range, Domain, Payload> &,
                 Solver &,
                 Preconditioner &&) = delete;

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


