//include/deal.II-translator/base/tensor_accessors_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#ifndef dealii_tensor_accessors_h
#define dealii_tensor_accessors_h

#include <deal.II/base/config.h>

#include <deal.II/base/table_indices.h>
#include <deal.II/base/template_constraints.h>


DEAL_II_NAMESPACE_OPEN

/**
 * 这个命名空间是一个在通用张量对象（任意等级）上工作的算法的集合。
 * 在一个单独的命名空间中以通用方式实现这种功能的理由是
 *
 *
 *
 *
 *
 * - 以方便代码的重用，从而避免代码的重复。
 *
 *
 *
 *
 *
 * - 要有一个定义明确的接口，可以交换低级别的实现。
 *
 * 一个张量对象有一个等级的概念，并允许索引操作符的等级-次数递归应用，例如，如果
 * <code>t</code>
 * 是一个等级为4的张量对象，下面的访问是有效的。
 *
 * @code
 * t[1][2][1][4]
 * @endcode
 *
 * deal.II对诸如 dealii::Tensor<rank,  dim, Number>和
 * dealii::SymmetricTensor<rank,  dim,
 * Number>这样的递归对象有自己的实现。
 * 在这个命名空间中实现的方法和算法，是完全通用的。更确切地说，它可以对嵌套的c-style数组进行操作，或者对具有最小接口的类类型
 * <code>T</code> 进行操作，该接口提供了一个本地别名
 * <code>value_type</code> 和一个索引操作符 <code>operator[](unsigned
 * int)</code> ，它返回一个 <code>value_type</code>
 * 的（常量或非常量）引用 。
 *
 * @code
 * template <...>
 * class T
 * {
 *   using value_type = ...;
 *   value_type & operator[](unsigned int);
 *   const value_type & operator[](unsigned int) const;
 * };
 * @endcode
 *
 * 这个命名空间提供了访问、重新排序和收缩此类对象的基元。
 *
 *
 * @ingroup geomprimitives
 *
 */
namespace TensorAccessors
{
  // forward declarations
  namespace internal
  {
    template <int index, int rank, typename T>
    class ReorderedIndexView;
    template <int position, int rank>
    struct ExtractHelper;
    template <int no_contr, int rank_1, int rank_2, int dim>
    class Contract;
    template <int rank_1, int rank_2, int dim>
    class Contract3;
  } // namespace internal


  /**
   * 这个类提供了一个本地别名 @p value_type
   * ，表示用operator[](unsigned
   * int)访问的结果类型。更准确地说， @p  value_type将是
   *
   *
   *
   *
   *
   *
   * -  <code>T::value_type</code>  如果T是一个提供别名的张量类 <code>value_type</code> ，并且没有const限定词。
   *
   *
   *
   *
   *
   * -  <code>const T::value_type</code>  如果T是一个提供别名的张量类 <code>value_type</code> ，并且确实有一个const限定词。
   *
   *
   *
   *
   *
   * -  <code>const T::value_type</code>  如果T是一个提供别名的张量类 <code>value_type</code> ，并且确实有一个const限定词。
   *
   *
   *
   * -  <code>A</code> if T is of array type <code>A[...]</code>
   *
   *
   *
   *
   * -  <code>const A</code> if T is of array type <code>A[...]</code> 并且确实有一个const限定词。
   *
   */
  template <typename T>
  struct ValueType
  {
    using value_type = typename T::value_type;
  };

  template <typename T>
  struct ValueType<const T>
  {
    using value_type = const typename T::value_type;
  };

  template <typename T, std::size_t N>
  struct ValueType<T[N]>
  {
    using value_type = T;
  };

  template <typename T, std::size_t N>
  struct ValueType<const T[N]>
  {
    using value_type = const T;
  };


  /**
   * 这个类提供了一个本地别名 @p value_type
   * ，在通过``operator[](unsigned int)``递归后的 @p deref_steps
   * 别名相等。此外，恒定性通过ValueType类型特征得以保留，即，如果T是恒定的，ReturnType<rank，
   * T>::value_type 也将是恒定的。
   *
   */
  template <int deref_steps, typename T>
  struct ReturnType
  {
    using value_type =
      typename ReturnType<deref_steps - 1,
                          typename ValueType<T>::value_type>::value_type;
  };

  template <typename T>
  struct ReturnType<0, T>
  {
    using value_type = T;
  };


  /**
   * 为等级为 @p rank 的张量对象的引用 @p t 提供一个
   * "张量视图"，其中索引 @p index
   * 被移到最后。作为一个例子，考虑一个dim=5空间维度的五阶张量，可以通过5个递归
   * <code>operator[]()</code> 的调用来访问。
   * @code
   * Tensor<5, dim> tensor;
   * tensor[0][1][2][3][4] = 42.;
   * @endcode
   * 索引1（第2个索引，计数从0开始）现在可以通过以下方式移到末尾
   * @code
   * auto tensor_view = reordered_index_view<1, 5>(tensor);
   * tensor_view[0][2][3][4][1] == 42.; // is true
   * @endcode
   * dealii::Tensor
   * 类型的用法完全是为了举例说明。这个函数实现的机制可用于相当普遍的张量类型
   * @p T.
   * 这个重排设施的目的是能够在两个（或多个）张量的任意索引上收缩。
   *
   *
   *
   *
   *
   * - 将指数重新排序，铭记在心，并将其置于张量的末端
   *
   *
   * - 使用下面的收缩函数，收缩张量的_最后一个元素。
   * @note  这个函数返回一个由数组下标操作符
   * <code>operator[](unsigned int)</code> 和描述其返回值的别名
   * <code>value_type</code> 组成的内部类对象。      @tparam  index
   * 要移到最后的索引。索引从0开始计算，因此有效范围是
   * $0\le\text{index}<\text{rank}$  。    @tparam  秩 递归对象的秩
   * @p t   @tparam  T 秩  @p rank.   @p T
   * 的递归对象必须提供一个本地别名  <code>value_type</code>
   * 和一个索引操作符  <code>operator[]()</code>  ，它返回一个
   * <code>value_type</code>  的（常量或非常量）引用.
   *
   */
  template <int index, int rank, typename T>
  constexpr DEAL_II_ALWAYS_INLINE internal::ReorderedIndexView<index, rank, T>
                                  reordered_index_view(T &t)
  {
    static_assert(0 <= index && index < rank,
                  "The specified index must lie within the range [0,rank)");

    return internal::ReorderedIndexView<index, rank, T>(t);
  }


  /**
   * 返回一个数组类型 @p ArrayType 对象 @p indices. 所描述的 @p
   * t 类型的一个子对象的引用（常量或非常量）。 例如。
   * @code
   * Tensor<5, dim> tensor;
   * TableIndices<5> indices (0, 1, 2, 3, 4);
   * TensorAccessors::extract(tensor, indices) = 42;
   * @endcode
   * 这相当于  <code>tensor[0][1][2][3][4] = 42.</code>  。      @tparam
   * T 一个等级为 @p rank.  @p T
   * 的递归对象必须提供一个本地别名 <code>value_type</code>
   * 和一个索引操作符 <code>operator[]()</code> ，它返回一个
   * <code>value_type</code> 的（常量或非常量）引用
   * 。此外，它的张力等级必须等于或大于  @p rank.   @tparam
   * ArrayType 一个类似数组的对象，例如  std::array,  或
   * dealii::TableIndices  ，至少存储  @p rank  索引，可以通过
   * operator[]() 访问。
   *
   */
  template <int rank, typename T, typename ArrayType>
  constexpr DEAL_II_ALWAYS_INLINE typename ReturnType<rank, T>::value_type &
  extract(T &t, const ArrayType &indices)
  {
    return internal::ExtractHelper<0, rank>::template extract<T, ArrayType>(
      t, indices);
  }


  /**
   * 这个函数收缩两个递归对象 @p left 和 @p right ，并将结果存储在 @p result. 中，收缩是在两个递归对象的_最后_ @p no_contr 索引上完成的。    @f[
   * \text{result}_{i_1,..,i_{r1},j_1,..,j_{r2}}
   * = \sum_{k_1,..,k_{\mathrm{no\_contr}}}
   *   \mathrm{left}_{i_1,..,i_{r1},k_1,..,k_{\mathrm{no\_contr}}}
   *   \mathrm{right}_{j_1,..,j_{r2},k_1,..,k_{\mathrm{no\_contr}}}
   * @f] 调用这个函数相当于编写以下低级代码。
   * @code
   * for(unsigned int i_0 = 0; i_0 < dim; ++i_0)
   *   ...
   *     for(unsigned int i_ = 0; i_ < dim; ++i_)
   *       for(unsigned int j_0 = 0; j_0 < dim; ++j_0)
   *         ...
   *           for(unsigned int j_ = 0; j_ < dim; ++j_)
   *             {
   *               result[i_0]..[i_][j_0]..[j_] = 0.;
   *               for(unsigned int k_0 = 0; k_0 < dim; ++k_0)
   *                 ...
   *                   for(unsigned int k_ = 0; k_ < dim; ++k_)
   *                     result[i_0]..[i_][j_0]..[j_] +=
   *                       left[i_0]..[i_][k_0]..[k_]
   *                         right[j_0]..[j_][k_0]..[k_];
   *             }
   * @endcode
   * 与r = rank_1 + rank_2
   *
   * - 2 no_contr, l = rank_1
   *
   * - no_contr, l1 = rank_1, and c = no_contr.
   * @note  类型 @p T1,   @p T2, 和 @p T3 必须具有rank rank_1 + rank_2
   *
   * - 2 no_contr, rank_1, 或 rank_2，分别。很明显，no_contr必须小于或等于rank_1和rank_2。
   *
   */
  template <int no_contr,
            int rank_1,
            int rank_2,
            int dim,
            typename T1,
            typename T2,
            typename T3>
  constexpr inline DEAL_II_ALWAYS_INLINE void
  contract(T1 &result, const T2 &left, const T3 &right)
  {
    static_assert(rank_1 >= no_contr,
                  "The rank of the left tensor must be "
                  "equal or greater than the number of "
                  "contractions");
    static_assert(rank_2 >= no_contr,
                  "The rank of the right tensor must be "
                  "equal or greater than the number of "
                  "contractions");

    internal::Contract<no_contr, rank_1, rank_2, dim>::
      template contract<T1, T2, T3>(result, left, right);
  }


  /**
   * 三个递归对象的完全收缩。    @f[
   * \sum_{i_1,..,i_{r1},j_1,..,j_{r2}}
   * \text{left}_{i_1,..,i_{r1}}
   * \text{middle}_{i_1,..,i_{r1},j_1,..,j_{r2}}
   * \text{right}_{j_1,..,j_{r2}}
   * @f] 调用此函数相当于编写以下低级代码。
   * @code
   * T1 result = T1();
   * for(unsigned int i_0 = 0; i_0 < dim; ++i_0)
   *   ...
   *     for(unsigned int i_ = 0; i_ < dim; ++i_)
   *       for(unsigned int j_0 = 0; j_0 < dim; ++j_0)
   *         ...
   *           for(unsigned int j_ = 0; j_ < dim; ++j_)
   *             result += left[i_0]..[i_]
   *                         middle[i_0]..[i_][j_0]..[j_]
   *                         right[j_0]..[j_];
   * @endcode
   * @note  类型 @p T2,   @p T3, 和 @p T4 必须分别具有等级rank_1,
   * rank_1 + rank_2, 和 rank_3。  @p T1 必须是一个标量类型。
   *
   */
  template <int rank_1,
            int rank_2,
            int dim,
            typename T1,
            typename T2,
            typename T3,
            typename T4>
  constexpr T1
  contract3(const T2 &left, const T3 &middle, const T4 &right)
  {
    return internal::Contract3<rank_1, rank_2, dim>::
      template contract3<T1, T2, T3, T4>(left, middle, right);
  }


  namespace internal
  {
    // -------------------------------------------------------------------------
    // Forward declarations and type traits
    // -------------------------------------------------------------------------

    template <int rank, typename S>
    class StoreIndex;
    template <typename T>
    class Identity;
    template <int no_contr, int dim>
    class Contract2;

    /**
     * 一个内部使用的类型特质，允许嵌套应用函数reordered_index_view(T
     * &t)。
     * 问题是，在处理实际的张量类型时，我们必须通过引用返回子张量
     *
     * - 但有时，特别是对于返回r值的StoreIndex和ReorderedIndexView，我们必须按值返回。
     *
     */
    template <typename T>
    struct ReferenceType
    {
      using type = T &;
    };

    template <int rank, typename S>
    struct ReferenceType<StoreIndex<rank, S>>
    {
      using type = StoreIndex<rank, S>;
    };

    template <int index, int rank, typename T>
    struct ReferenceType<ReorderedIndexView<index, rank, T>>
    {
      using type = ReorderedIndexView<index, rank, T>;
    };


    // TODO: Is there a possibility to just have the following block of
    // explanation on an internal page in doxygen? If, yes. Doxygen
    // wizards, your call!

    // -------------------------------------------------------------------------
    // Implementation of helper classes for reordered_index_view
    // -------------------------------------------------------------------------

    // OK. This is utterly brutal template magic. Therefore, we will not
    // comment on the individual internal helper classes, because this is
    // of not much value, but explain the general recursion procedure.
    //
    // (In order of appearance)
    //
    // Our task is to reorder access to a tensor object where a specified
    // index is moved to the end. Thus we want to construct an object
    // <code>reordered</code> out of a <code>tensor</code> where the
    // following access patterns are equivalent:
    // @code
    //   tensor    [i_0]...[i_index-1][i_index][i_index+1]...[i_n]
    //   reordered [i_0]...[i_index_1][i_index+1]...[i_n][i_index]
    // @endcode
    //
    // The first task is to get rid of the application of
    // [i_0]...[i_index-1]. This is a classical recursion pattern - relay
    // the task from <index, rank> to <index-1, rank-1> by accessing the
    // subtensor object:

    template <int index, int rank, typename T>
    class ReorderedIndexView
    {
    public:
      constexpr ReorderedIndexView(typename ReferenceType<T>::type t)
        : t_(t)
      {}

      using value_type = ReorderedIndexView<index - 1,
                                            rank - 1,
                                            typename ValueType<T>::value_type>;

      // Recurse by applying index j directly:
      constexpr DEAL_II_ALWAYS_INLINE value_type
                                      operator[](unsigned int j) const
      {
        return value_type(t_[j]);
      }

    private:
      typename ReferenceType<T>::type t_;
    };

    // At some point we hit the condition index == 0 and rank > 1, i.e.,
    // the first index should be reordered to the end.
    //
    // At this point we cannot be lazy any more and have to start storing
    // indices because we get them in the wrong order. The user supplies
    //   [i_0][i_1]...[i_{rank - 1}]
    // but we have to call the subtensor object with
    //   [i_{rank - 1}[i_0][i_1]...[i_{rank-2}]
    //
    // So give up and relay the task to the StoreIndex class:

    template <int rank, typename T>
    class ReorderedIndexView<0, rank, T>
    {
    public:
      constexpr ReorderedIndexView(typename ReferenceType<T>::type t)
        : t_(t)
      {}

      using value_type = StoreIndex<rank - 1, internal::Identity<T>>;

      constexpr DEAL_II_ALWAYS_INLINE value_type
                                      operator[](unsigned int j) const
      {
        return value_type(Identity<T>(t_), j);
      }

    private:
      typename ReferenceType<T>::type t_;
    };

    // Sometimes, we're lucky and don't have to do anything. In this case
    // just return the original tensor.

    template <typename T>
    class ReorderedIndexView<0, 1, T>
    {
    public:
      constexpr ReorderedIndexView(typename ReferenceType<T>::type t)
        : t_(t)
      {}

      using value_type =
        typename ReferenceType<typename ValueType<T>::value_type>::type;

      constexpr DEAL_II_ALWAYS_INLINE value_type
                                      operator[](unsigned int j) const
      {
        return t_[j];
      }

    private:
      typename ReferenceType<T>::type t_;
    };

    // Here, Identity is a helper class to ground the recursion in
    // StoreIndex. Its implementation is easy - we haven't stored any
    // indices yet. So, we just provide a function apply that returns the
    // application of an index j to the stored tensor t_:

    template <typename T>
    class Identity
    {
    public:
      constexpr Identity(typename ReferenceType<T>::type t)
        : t_(t)
      {}

      using return_type = typename ValueType<T>::value_type;

      constexpr DEAL_II_ALWAYS_INLINE typename ReferenceType<return_type>::type
      apply(unsigned int j) const
      {
        return t_[j];
      }

    private:
      typename ReferenceType<T>::type t_;
    };

    // StoreIndex is a class that stores an index recursively with every
    // invocation of operator[](unsigned int j): We do this by recursively
    // creating a new StoreIndex class of lower rank that stores the
    // supplied index j and holds a copy of the current class (with all
    // other stored indices). Again, we provide an apply member function
    // that knows how to apply an index on the highest rank and all
    // subsequently stored indices:

    template <int rank, typename S>
    class StoreIndex
    {
    public:
      constexpr StoreIndex(S s, int i)
        : s_(s)
        , i_(i)
      {}

      using value_type = StoreIndex<rank - 1, StoreIndex<rank, S>>;

      constexpr DEAL_II_ALWAYS_INLINE value_type
                                      operator[](unsigned int j) const
      {
        return value_type(*this, j);
      }

      using return_type =
        typename ValueType<typename S::return_type>::value_type;

      constexpr typename ReferenceType<return_type>::type
      apply(unsigned int j) const
      {
        return s_.apply(j)[i_];
      }

    private:
      const S   s_;
      const int i_;
    };

    // We have to store indices until we hit rank == 1. Then, upon the next
    // invocation of operator[](unsigned int j) we have all necessary
    // information available to return the actual object.

    template <typename S>
    class StoreIndex<1, S>
    {
    public:
      constexpr StoreIndex(S s, int i)
        : s_(s)
        , i_(i)
      {}

      using return_type =
        typename ValueType<typename S::return_type>::value_type;
      using value_type = return_type;

      constexpr DEAL_II_ALWAYS_INLINE return_type &
                                      operator[](unsigned int j) const
      {
        return s_.apply(j)[i_];
      }

    private:
      const S   s_;
      const int i_;
    };


    // -------------------------------------------------------------------------
    // Implementation of helper classes for extract
    // -------------------------------------------------------------------------

    // Straightforward recursion implemented by specializing ExtractHelper
    // for position == rank. We use the type trait ReturnType<rank, T> to
    // have an idea what the final type will be.
    template <int position, int rank>
    struct ExtractHelper
    {
      template <typename T, typename ArrayType>
      constexpr static typename ReturnType<rank - position, T>::value_type &
      extract(T &t, const ArrayType &indices)
      {
        return ExtractHelper<position + 1, rank>::template extract<
          typename ValueType<T>::value_type,
          ArrayType>(t[indices[position]], indices);
      }
    };

    // For position == rank there is nothing to extract, just return the
    // object.
    template <int rank>
    struct ExtractHelper<rank, rank>
    {
      template <typename T, typename ArrayType>
      constexpr static T &
      extract(T &t, const ArrayType &)
      {
        return t;
      }
    };


    // -------------------------------------------------------------------------
    // Implementation of helper classes for contract
    // -------------------------------------------------------------------------

    // Straightforward recursive pattern:
    //
    // As long as rank_1 > no_contr, assign indices from the left tensor to
    // result. This builds up the first part of the nested outer loops:
    //
    // for(unsigned int i_0; i_0 < dim; ++i_0)
    //   ...
    //     for(i_; i_ < dim; ++i_)
    //       [...]
    //         result[i_0]..[i_] ... left[i_0]..[i_] ...

    template <int no_contr, int rank_1, int rank_2, int dim>
    class Contract
    {
    public:
      template <typename T1, typename T2, typename T3>
      constexpr inline DEAL_II_ALWAYS_INLINE static void
      contract(T1 &result, const T2 &left, const T3 &right)
      {
        for (unsigned int i = 0; i < dim; ++i)
          Contract<no_contr, rank_1 - 1, rank_2, dim>::contract(result[i],
                                                                left[i],
                                                                right);
      }
    };

    // If rank_1 == no_contr leave out the remaining no_contr indices for
    // the contraction and assign indices from the right tensor to the
    // result. This builds up the second part of the nested loops:
    //
    //  for(unsigned int i_0 = 0; i_0 < dim; ++i_0)
    //    ...
    //      for(unsigned int i_ = 0; i_ < dim; ++i_)
    //        for(unsigned int j_0 = 0; j_0 < dim; ++j_0)
    //          ...
    //            for(unsigned int j_ = 0; j_ < dim; ++j_)
    //             [...]
    //               result[i_0]..[i_][j_0]..[j_] ... left[i_0]..[i_] ...
    //               right[j_0]..[j_]
    //

    template <int no_contr, int rank_2, int dim>
    class Contract<no_contr, no_contr, rank_2, dim>
    {
    public:
      template <typename T1, typename T2, typename T3>
      constexpr inline DEAL_II_ALWAYS_INLINE static void
      contract(T1 &result, const T2 &left, const T3 &right)
      {
        for (unsigned int i = 0; i < dim; ++i)
          Contract<no_contr, no_contr, rank_2 - 1, dim>::contract(result[i],
                                                                  left,
                                                                  right[i]);
      }
    };

    // If rank_1 == rank_2 == no_contr we have built up all of the outer
    // loop. Now, it is time to do the actual contraction:
    //
    // [...]
    //   {
    //     result[i_0]..[i_][j_0]..[j_] = 0.;
    //     for(unsigned int k_0 = 0; k_0 < dim; ++k_0)
    //       ...
    //         for(unsigned int k_ = 0; k_ < dim; ++k_)
    //           result[i_0]..[i_][j_0]..[j_] += left[i_0]..[i_][k_0]..[k_] *
    //           right[j_0]..[j_][k_0]..[k_];
    //   }
    //
    //  Relay this summation to another helper class.

    template <int no_contr, int dim>
    class Contract<no_contr, no_contr, no_contr, dim>
    {
    public:
      template <typename T1, typename T2, typename T3>
      constexpr inline DEAL_II_ALWAYS_INLINE static void
      contract(T1 &result, const T2 &left, const T3 &right)
      {
        result = Contract2<no_contr, dim>::template contract2<T1>(left, right);
      }
    };

    // Straightforward recursion:
    //
    // Contract leftmost index and recurse one down.

    template <int no_contr, int dim>
    class Contract2
    {
    public:
      template <typename T1, typename T2, typename T3>
      constexpr inline DEAL_II_ALWAYS_INLINE static T1
      contract2(const T2 &left, const T3 &right)
      {
        // Some auto-differentiable numbers need explicit
        // zero initialization.
        if (dim == 0)
          {
            T1 result = dealii::internal::NumberType<T1>::value(0.0);
            return result;
          }
        else
          {
            T1 result =
              Contract2<no_contr - 1, dim>::template contract2<T1>(left[0],
                                                                   right[0]);
            for (unsigned int i = 1; i < dim; ++i)
              result +=
                Contract2<no_contr - 1, dim>::template contract2<T1>(left[i],
                                                                     right[i]);
            return result;
          }
      }
    };

    // A contraction of two objects of order 0 is just a scalar
    // multiplication:

    template <int dim>
    class Contract2<0, dim>
    {
    public:
      template <typename T1, typename T2, typename T3>
      constexpr DEAL_II_ALWAYS_INLINE static T1
      contract2(const T2 &left, const T3 &right)
      {
        return left * right;
      }
    };


    // -------------------------------------------------------------------------
    // Implementation of helper classes for contract3
    // -------------------------------------------------------------------------

    // Fully contract three tensorial objects
    //
    // As long as rank_1 > 0, recurse over left and middle:
    //
    // for(unsigned int i_0; i_0 < dim; ++i_0)
    //   ...
    //     for(i_; i_ < dim; ++i_)
    //       [...]
    //         left[i_0]..[i_] ... middle[i_0]..[i_] ... right

    template <int rank_1, int rank_2, int dim>
    class Contract3
    {
    public:
      template <typename T1, typename T2, typename T3, typename T4>
      constexpr static inline T1
      contract3(const T2 &left, const T3 &middle, const T4 &right)
      {
        // Some auto-differentiable numbers need explicit
        // zero initialization.
        T1 result = dealii::internal::NumberType<T1>::value(0.0);
        for (unsigned int i = 0; i < dim; ++i)
          result += Contract3<rank_1 - 1, rank_2, dim>::template contract3<T1>(
            left[i], middle[i], right);
        return result;
      }
    };

    // If rank_1 ==0, continue to recurse over middle and right:
    //
    // for(unsigned int i_0; i_0 < dim; ++i_0)
    //   ...
    //     for(i_; i_ < dim; ++i_)
    //       for(unsigned int j_0; j_0 < dim; ++j_0)
    //         ...
    //           for(j_; j_ < dim; ++j_)
    //             [...]
    //               left[i_0]..[i_] ... middle[i_0]..[i_][j_0]..[j_] ...
    //               right[j_0]..[j_]

    template <int rank_2, int dim>
    class Contract3<0, rank_2, dim>
    {
    public:
      template <typename T1, typename T2, typename T3, typename T4>
      constexpr static inline T1
      contract3(const T2 &left, const T3 &middle, const T4 &right)
      {
        // Some auto-differentiable numbers need explicit
        // zero initialization.
        T1 result = dealii::internal::NumberType<T1>::value(0.0);
        for (unsigned int i = 0; i < dim; ++i)
          result +=
            Contract3<0, rank_2 - 1, dim>::template contract3<T1>(left,
                                                                  middle[i],
                                                                  right[i]);
        return result;
      }
    };

    // Contraction of three tensorial objects of rank 0 is just a scalar
    // multiplication.

    template <int dim>
    class Contract3<0, 0, dim>
    {
    public:
      template <typename T1, typename T2, typename T3, typename T4>
      constexpr static T1
      contract3(const T2 &left, const T3 &middle, const T4 &right)
      {
        return left * middle * right;
      }
    };

    // -------------------------------------------------------------------------

  }  /* namespace internal */ 
}  /* namespace TensorAccessors */ 

DEAL_II_NAMESPACE_CLOSE

#endif  /* dealii_tensor_accessors_h */ 


