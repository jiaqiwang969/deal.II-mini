//include/deal.II-translator/lac/packaged_operation_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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

#ifndef dealii_packaged_operation_h
#define dealii_packaged_operation_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/vector_memory.h>

#include <functional>

DEAL_II_NAMESPACE_OPEN

// Forward declarations:
#ifndef DOXYGEN
template <typename Number>
class Vector;
template <typename Range, typename Domain, typename Payload>
class LinearOperator;
template <typename Range = Vector<double>>
class PackagedOperation;
#endif


/**
 * 一个用于存储计算的类。
 * PackagedOperation类允许对涉及向量和线性运算符的表达式进行懒惰评估。这是通过存储计算表达式来实现的，只有当对象被隐含地转换为向量对象，或者
 * <code>apply</code> （或 <code>apply_add</code>
 * ）被手动调用时才执行计算。这就避免了不必要的中间结果的临时存储。
 * 该类本质上由 <code>std::function</code>
 * 对象组成，这些对象存储了如何生成计算结果并将其存储在向量中的知识。
 *
 * @code
 * std::function<void(Range &)> apply;
 * std::function<void(Range &)> apply_add;
 * @endcode
 *
 * 与LinearOperator类相似，它也有关于如何初始化 @p Range
 * 空间的向量的知识。
 *
 * @code
 * std::function<void(Range &, bool)> reinit_vector;
 * @endcode
 *
 * 作为一个例子，考虑多个向量的相加问题
 *
 * @code
 * Vector<double> a, b, c, d;
 * // ..
 * Vector<double> result = a + b
 *
 * - c + d;
 * @endcode
 * 或者计算一个残差  $b-Ax$  。
 *
 * @code
 * SparseMatrix<double> A;
 * Vector<double> b, x;
 * // ..
 * const auto op_a = linear_operator(A);
 *
 * auto residual =  b
 *
 * - op_a x;
 * @endcode
 * 表达式  <code>residual</code>  是
 * <code>PackagedOperation<Vector<double>></code>  的类型。它存储对
 * <code>A</code>, <code>b</code> and <code>x</code>
 * 的引用，并将实际计算推迟到 <code>apply</code>, or
 * <code>apply_add</code> 被明确调用时进行。
 *
 * @code
 * Vector<double> y;
 * residual.reinit_vector(y);
 * residual.apply(y);
 * residual.apply_add(y);
 * @endcode
 * 或者直到 @p PackagedOperation 对象被隐式转换。
 *
 * @code
 * Vector<double> y;
 * y = residual;
 * y += residual;
 * y
 *
 * -= residual;
 * @endcode
 *
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
template <typename Range>
class PackagedOperation
{
public:
  /**
   * 创建一个空的PackagedOperation对象。所有的
   * <code>std::function</code>
   * 成员对象都用默认的变体初始化，在调用时抛出一个异常。
   *
   */
  PackagedOperation()
  {
    apply = [](Range &) {
      Assert(false,
             ExcMessage(
               "Uninitialized PackagedOperation<Range>::apply called"));
    };

    apply_add = [](Range &) {
      Assert(false,
             ExcMessage(
               "Uninitialized PackagedOperation<Range>::apply_add called"));
    };

    reinit_vector = [](Range &, bool) {
      Assert(false,
             ExcMessage("Uninitialized PackagedOperation<Range>::reinit_vector "
                        "method called"));
    };
  }

  /**
   * 默认的复制构造函数。
   *
   */
  PackagedOperation(const PackagedOperation<Range> &) = default;

  /**
   * 从引用向量创建PackagedOperation对象的构造函数  @p u.
   * PackagedOperation返回  @p u.
   * 被创建的PackagedOperation对象存储了对  @p u.
   * 的引用。因此，向量必须在PackagedOperation对象的整个生命周期内保持有效引用。在创建PackagedOperation对象后，在
   * @p u 上所做的所有改变都会被操作者对象所反映。
   *
   */
  PackagedOperation(const Range &u)
  {
    *this = u;
  }

  /**
   * 默认的复制赋值运算符。
   *
   */
  PackagedOperation<Range> &
  operator=(const PackagedOperation<Range> &) = default;

  /**
   * 复制赋值操作符，从一个引用向量创建一个PackagedOperation对象
   * @p u.  PackagedOperation返回  @p u.
   * 被创建的PackagedOperation对象存储了对  @p u.
   * 的引用。因此，向量必须在PackagedOperation对象的整个生命周期中保持有效引用。在PackagedOperation对象创建后，在
   * @p u 上所做的所有改变都会被操作者对象所反映。
   *
   */
  PackagedOperation<Range> &
  operator=(const Range &u)
  {
    apply = [&u](Range &v) { v = u; };

    apply_add = [&u](Range &v) { v += u; };

    reinit_vector = [&u](Range &v, bool omit_zeroing_entries) {
      v.reinit(u, omit_zeroing_entries);
    };

    return *this;
  }

  /**
   * 将一个PackagedOperation转换为其结果。
   * 这个转换操作者创建一个Range空间的向量，并对其调用
   * <code>apply()</code> 。
   *
   */
  operator Range() const
  {
    Range result_vector;

    reinit_vector(result_vector,  /*bool omit_zeroing_entries=*/ true);
    apply(result_vector);

    return result_vector;
  }

  /**
   * @name  就地的向量空间操作
   *
   */
  //@{

  /**
   * 用具有相同 @p Range. 的PackagedOperation @p second_comp 加法。
   *
   */
  PackagedOperation<Range> &
  operator+=(const PackagedOperation<Range> &second_comp)
  {
    *this = *this + second_comp;
    return *this;
  }

  /**
   * 用相同的 @p second_comp 的打包操作 @p 做减法，范围。
   *
   */
  PackagedOperation<Range> &
  operator-=(const PackagedOperation<Range> &second_comp)
  {
    *this = *this - second_comp;
    return *this;
  }

  /**
   * 将一个常数 @p offset （属于 @p Range
   * 空间）添加到一个PackagedOperation的结果中。
   *
   */
  PackagedOperation<Range> &
  operator+=(const Range &offset)
  {
    *this = *this + PackagedOperation<Range>(offset);
    return *this;
  }

  /**
   * 从一个PackagedOperation的结果中减去一个常数 @p offset
   * （属于 @p Range 空间）。
   *
   */
  PackagedOperation<Range> &
  operator-=(const Range &offset)
  {
    *this = *this - PackagedOperation<Range>(offset);
    return *this;
  }

  /**
   * 将PackagedOperation与一个 @p number. 的标量乘法。
   *
   */
  PackagedOperation<Range> &
  operator*=(typename Range::value_type number)
  {
    *this = *this * number;
    return *this;
  }
  //@}

  /**
   * 将PackagedOperation的结果存储在 @p Range
   * 空间的一个向量v中。
   *
   */
  std::function<void(Range &v)> apply;

  /**
   * 将PackagedOperation的结果添加到 @p Range
   * 空间的一个向量v中。
   *
   */
  std::function<void(Range &v)> apply_add;

  /**
   * 初始化Range空间的一个向量v，使其可以直接作为apply或apply_add应用中的目标参数使用。与向量类的reinit函数类似，布尔值决定是否进行快速初始化，即如果它被设置为false，则向量的内容被设置为0。
   *
   */
  std::function<void(Range &v, bool omit_zeroing_entries)> reinit_vector;
};


/**
 * @name  矢量空间操作
 *
 *
 */
//@{

/**
 * @relatesalso PackagedOperation（PackagedOperation
 * 两个PackagedOperation对象 @p first_comp 和 @p second_comp
 * 的加法，通过向量空间加法给出相应结果。
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Range>
PackagedOperation<Range>
operator+(const PackagedOperation<Range> &first_comp,
          const PackagedOperation<Range> &second_comp)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = first_comp.reinit_vector;

  // ensure to have valid PackagedOperation objects by catching first_comp and
  // second_comp by value

  return_comp.apply = [first_comp, second_comp](Range &v) {
    first_comp.apply(v);
    second_comp.apply_add(v);
  };

  return_comp.apply_add = [first_comp, second_comp](Range &v) {
    first_comp.apply_add(v);
    second_comp.apply_add(v);
  };

  return return_comp;
}

/**
 * @relatesalso PackagedOperation（PackagedOperation
 * 两个PackagedOperation对象 @p first_comp 和 @p
 * second_comp的减法，由相应结果的向量空间相加给出。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range>
PackagedOperation<Range>
operator-(const PackagedOperation<Range> &first_comp,
          const PackagedOperation<Range> &second_comp)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = first_comp.reinit_vector;

  // ensure to have valid PackagedOperation objects by catching first_comp and
  // second_comp by value

  return_comp.apply = [first_comp, second_comp](Range &v) {
    second_comp.apply(v);
    v *= -1.;
    first_comp.apply_add(v);
  };

  return_comp.apply_add = [first_comp, second_comp](Range &v) {
    first_comp.apply_add(v);
    v *= -1.;
    second_comp.apply_add(v);
    v *= -1.;
  };

  return return_comp;
}

/**
 * @relatesalso PackagedOperation（PackagedOperation
 * 一个PackagedOperation对象 @p comp 与一个标量 @p number
 * 的标量乘法，这个标量是由PackagedOperation的结果与 @p
 * number. 的标量给出的。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range>
PackagedOperation<Range> operator*(const PackagedOperation<Range> &comp,
                                   typename Range::value_type      number)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = comp.reinit_vector;

  // the trivial case: number is zero
  if (number == 0.)
    {
      return_comp.apply = [](Range &v) { v = 0.; };

      return_comp.apply_add = [](Range &) {};
    }
  else
    {
      return_comp.apply = [comp, number](Range &v) {
        comp.apply(v);
        v *= number;
      };

      return_comp.apply_add = [comp, number](Range &v) {
        v /= number;
        comp.apply_add(v);
        v *= number;
      };
    }

  return return_comp;
}

/**
 * @relatesalso PackagedOperation（PackagedOperation
 * 一个PackagedOperation对象 @p comp 与一个标量 @p number
 * 的标量乘法，这个标量是由PackagedOperation的结果与 @p
 * number. 的标量给定的。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range>
PackagedOperation<Range> operator*(typename Range::value_type      number,
                                   const PackagedOperation<Range> &comp)
{
  return comp * number;
}

/**
 * @relatesalso PackagedOperation（PackagedOperation
 * 在PackagedOperation的结果中添加一个常数 @p offset （属于 @p
 * Range 空间）。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range>
PackagedOperation<Range>
operator+(const PackagedOperation<Range> &comp, const Range &offset)
{
  return comp + PackagedOperation<Range>(offset);
}

/**
 * @relatesalso PackagedOperation（PackagedOperation
 * 在PackagedOperation的结果中添加一个常数 @p offset （属于 @p
 * Range 空间）。
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Range>
PackagedOperation<Range>
operator+(const Range &offset, const PackagedOperation<Range> &comp)
{
  return PackagedOperation<Range>(offset) + comp;
}

/**
 * @relatesalso PackagedOperation（PackagedOperation
 * 从一个PackagedOperation的结果中减去一个常数 @p offset
 * （属于 @p Range 空间）。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range>
PackagedOperation<Range>
operator-(const PackagedOperation<Range> &comp, const Range &offset)
{
  return comp - PackagedOperation<Range>(offset);
}


/**
 * @relatesalso PackagedOperation 从一个常数 @p offset （属于 @p Range
 * 空间）减去一个计算结果。其结果是一个PackagedOperation对象，应用这个计算。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range>
PackagedOperation<Range>
operator-(const Range &offset, const PackagedOperation<Range> &comp)
{
  return PackagedOperation<Range>(offset) - comp;
}

//@}


/**
 * @name  创建一个PackagedOperation对象
 *
 *
 */
//@{

namespace internal
{
  namespace PackagedOperationImplementation
  {
    // Poor man's trait class that determines whether type T is a vector:
    // FIXME: Implement this as a proper type trait - similar to
    // isBlockVector

    template <typename T>
    class has_vector_interface
    {
      template <typename C>
      static std::false_type
      test(...);

      template <typename C>
      static std::true_type
      test(decltype(&C::operator+=),
           decltype(&C::operator-=),
           decltype(&C::l2_norm));

    public:
      // type is std::true_type if Matrix provides vmult_add and Tvmult_add,
      // otherwise it is std::false_type

      using type = decltype(test<T>(nullptr, nullptr, nullptr));
    }; // namespace
  }    // namespace PackagedOperationImplementation
} // namespace internal


/**
 * @relatesalso PackagedOperation（PackagedOperation
 * 创建一个PackagedOperation对象，存储两个向量的相加。
 * 创建的PackagedOperation对象存储了对 @p u 和 @p v.
 * 的引用。因此，向量必须在PackagedOperation对象的整个生命周期内保持有效引用。在创建PackagedOperation对象后，在
 * @p u 或 @p v 上所做的所有更改都会被操作者对象所反映。
 *
 *
 * @ingroup LAOperators
 *
 *
 */

template <typename Range,
          typename = typename std::enable_if<
            internal::PackagedOperationImplementation::has_vector_interface<
              Range>::type::value>::type>
PackagedOperation<Range>
operator+(const Range &u, const Range &v)
{
  PackagedOperation<Range> return_comp;

  // ensure to have valid PackagedOperation objects by catching op by value
  // u is caught by reference

  return_comp.reinit_vector = [&u](Range &x, bool omit_zeroing_entries) {
    x.reinit(u, omit_zeroing_entries);
  };

  return_comp.apply = [&u, &v](Range &x) {
    x = u;
    x += v;
  };

  return_comp.apply_add = [&u, &v](Range &x) {
    x += u;
    x += v;
  };

  return return_comp;
}


/**
 * @relatesalso PackagedOperation（PackagedOperation
 * 创建一个PackagedOperation对象，存储两个向量的减法。
 * 创建的PackagedOperation对象存储了对 @p u 和 @p v.
 * 的引用。因此，这些向量必须在PackagedOperation对象的整个生命周期内保持有效的引用。在创建PackagedOperation对象后，在
 * @p u 或 @p v 上所做的所有改变都会被操作者对象所反映。
 *
 *
 * @ingroup LAOperators
 *
 *
 */

template <typename Range,
          typename = typename std::enable_if<
            internal::PackagedOperationImplementation::has_vector_interface<
              Range>::type::value>::type>
PackagedOperation<Range>
operator-(const Range &u, const Range &v)
{
  PackagedOperation<Range> return_comp;

  // ensure to have valid PackagedOperation objects by catching op by value
  // u is caught by reference

  return_comp.reinit_vector = [&u](Range &x, bool omit_zeroing_entries) {
    x.reinit(u, omit_zeroing_entries);
  };

  return_comp.apply = [&u, &v](Range &x) {
    x = u;
    x -= v;
  };

  return_comp.apply_add = [&u, &v](Range &x) {
    x += u;
    x -= v;
  };

  return return_comp;
}


/**
 * @relatesalso PackagedOperation（PackagedOperation
 * 创建一个PackagedOperation对象，该对象存储一个具有 @p
 * number. 的向量的缩放比例。
 * 创建的PackagedOperation对象存储了对 @p u.
 * 的引用。因此，向量必须在PackagedOperation对象的整个生命周期内保持有效的引用。在创建PackagedOperation对象后，在
 * @p u 或 @p v 上所做的所有改变都会被操作者对象所反映。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range,
          typename = typename std::enable_if<
            internal::PackagedOperationImplementation::has_vector_interface<
              Range>::type::value>::type>
PackagedOperation<Range> operator*(const Range &              u,
                                   typename Range::value_type number)
{
  return PackagedOperation<Range>(u) * number;
}


/**
 * @relatesalso PackagedOperation（PackagedOperation
 * 创建一个PackagedOperation对象，该对象存储了对一个具有 @p
 * number. 的向量的缩放。 创建的PackagedOperation对象存储了对
 * @p u.
 * 的引用。因此，向量必须在PackagedOperation对象的整个生命周期内保持有效的引用。在创建PackagedOperation对象后，在
 * @p u 或 @p v 上所做的所有改变都会被操作者对象所反映。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range,
          typename = typename std::enable_if<
            internal::PackagedOperationImplementation::has_vector_interface<
              Range>::type::value>::type>
PackagedOperation<Range> operator*(typename Range::value_type number,
                                   const Range &              u)
{
  return number * PackagedOperation<Range>(u);
}


/**
 * @relatesalso PackagedOperation（PackagedOperation
 * 从一个LinearOperator和对域空间的向量 @p u
 * 的引用中创建一个PackagedOperation对象。该对象存储PackagedOperation
 * $\text{op} \,u$  （用矩阵表示）。  <code>return</code>  (
 * <code>return_add</code>) are implemented with <code>vmult(__1,u)</code>  (
 * <code>vmult_add(__1,u)</code>  ) 。
 * 创建的PackagedOperation对象存储了对 @p u.
 * 的引用。因此，在PackagedOperation对象的整个生命周期中，该向量必须保持有效的引用。在创建PackagedOperation对象后，在
 * @p u 上所做的所有改变都会被操作者对象所反映。
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Range, typename Domain, typename Payload>
PackagedOperation<Range>
operator*(const LinearOperator<Range, Domain, Payload> &op, const Domain &u)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = op.reinit_range_vector;

  // ensure to have valid PackagedOperation objects by catching op by value
  // u is caught by reference

  return_comp.apply = [op, &u](Range &v) { op.vmult(v, u); };

  return_comp.apply_add = [op, &u](Range &v) { op.vmult_add(v, u); };

  return return_comp;
}


/**
 * @relatesalso PackagedOperation（PackagedOperation
 * 从一个LinearOperator和一个对Range空间的向量 @p u
 * 的引用创建一个PackagedOperation对象。该对象存储PackagedOperation
 * $\text{op}^T \,u$  （用矩阵符号表示）。  <code>return</code>  (
 * <code>return_add</code>) are implemented with <code>Tvmult(__1,u)</code>  (
 * <code>Tvmult_add(__1,u)</code>  ) 。
 * 创建的PackagedOperation对象存储了对 @p u.
 * 的引用。因此，在PackagedOperation对象的整个生命周期内，该向量必须保持有效引用。在创建PackagedOperation对象后，在
 * @p u 上所做的所有改变都会被操作者对象所反映。
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Range, typename Domain, typename Payload>
PackagedOperation<Domain>
operator*(const Range &u, const LinearOperator<Range, Domain, Payload> &op)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = op.reinit_domain_vector;

  // ensure to have valid PackagedOperation objects by catching op by value
  // u is caught by reference

  return_comp.apply = [op, &u](Domain &v) { op.Tvmult(v, u); };

  return_comp.apply_add = [op, &u](Domain &v) { op.Tvmult_add(v, u); };

  return return_comp;
}


/**
 * @relatesalso PackagedOperation(PackagedOperation)
 * 一个PackagedOperation对象与一个LinearOperator的组合。该对象存储了计算
 * $\text{op} \,comp$ （用矩阵符号表示）。
 *
 *
 * @ingroup LAOperators
 *
 */
template <typename Range, typename Domain, typename Payload>
PackagedOperation<Range>
operator*(const LinearOperator<Range, Domain, Payload> &op,
          const PackagedOperation<Domain> &             comp)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = op.reinit_range_vector;

  // ensure to have valid PackagedOperation objects by catching op by value
  // u is caught by reference

  return_comp.apply = [op, comp](Domain &v) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer i(vector_memory);
    op.reinit_domain_vector(*i,  /*bool omit_zeroing_entries =*/ true);

    comp.apply(*i);
    op.vmult(v, *i);
  };

  return_comp.apply_add = [op, comp](Domain &v) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer i(vector_memory);
    op.reinit_range_vector(*i,  /*bool omit_zeroing_entries =*/ true);

    comp.apply(*i);
    op.vmult_add(v, *i);
  };

  return return_comp;
}


/**
 * @relatesalso PackagedOperation(PackagedOperation)
 * 一个PackagedOperation对象与一个LinearOperator的组合。该对象存储了计算
 * $\text{op}^T \,comp$ （用矩阵符号表示）。
 *
 *
 * @ingroup LAOperators
 *
 *
 */
template <typename Range, typename Domain, typename Payload>
PackagedOperation<Domain>
operator*(const PackagedOperation<Range> &              comp,
          const LinearOperator<Range, Domain, Payload> &op)
{
  PackagedOperation<Range> return_comp;

  return_comp.reinit_vector = op.reinit_domain_vector;

  // ensure to have valid PackagedOperation objects by catching op by value
  // u is caught by reference

  return_comp.apply = [op, comp](Domain &v) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer i(vector_memory);
    op.reinit_range_vector(*i,  /*bool omit_zeroing_entries =*/ true);

    comp.apply(*i);
    op.Tvmult(v, *i);
  };

  return_comp.apply_add = [op, comp](Domain &v) {
    GrowingVectorMemory<Range> vector_memory;

    typename VectorMemory<Range>::Pointer i(vector_memory);
    op.reinit_range_vector(*i,  /*bool omit_zeroing_entries =*/ true);

    comp.apply(*i);
    op.Tvmult_add(v, *i);
  };

  return return_comp;
}

//@}

DEAL_II_NAMESPACE_CLOSE

#endif


