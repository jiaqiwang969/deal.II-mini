//include/deal.II-translator/base/parallel_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2020 by the deal.II authors
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

#ifndef dealii_parallel_h
#define dealii_parallel_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/synchronous_iterator.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/thread_management.h>

#include <cstddef>
#include <functional>
#include <memory>
#include <tuple>

#ifdef DEAL_II_WITH_TBB
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <tbb/blocked_range.h>
#  include <tbb/parallel_for.h>
#  include <tbb/parallel_reduce.h>
#  include <tbb/partitioner.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
#endif


// TODO[WB]: allow calling functions to pass along a tbb::affinity_partitioner
// object to ensure that subsequent calls use the same cache lines

DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace internal
  {
    /**
     * 帮助结构，告诉我们是否可以对给定的 @p Number
     * 类型使用SIMD指令。
     *
     */
    template <typename Number>
    struct EnableOpenMPSimdFor
    {
      static const bool value = true;
    };

#ifdef __INTEL_COMPILER
    // Disable long double SIMD instructions on ICC. This is to work around a
    // bug that generates wrong code at least up to intel 15 (see
    // tests/lac/vector-vector, tests/lac/intel-15-bug, and the discussion at
    // https://github.com/dealii/dealii/issues/598).
    template <>
    struct EnableOpenMPSimdFor<long double>
    {
      static const bool value = false;
    };
#endif



    /**
     * 将F类型的函数对象转换成可以应用于同步迭代器范围的所有元素的对象。
     *
     */
    template <typename F>
    struct Body
    {
      /**
       * 构造器。取出并打包给定的函数对象。
       *
       */
      Body(const F &f)
        : f(f)
      {}

      template <typename Range>
      void
      operator()(const Range &range) const
      {
        for (typename Range::const_iterator p = range.begin(); p != range.end();
             ++p)
          apply(f, *p);
      }

    private:
      /**
       * 存储的函数对象。
       *
       */
      const F f;

      /**
       * 将F应用于一个有两个元素的迭代器集合。
       *
       */
      template <typename I1, typename I2>
      static void
      apply(const F &f, const std::tuple<I1, I2> &p)
      {
        *std::get<1>(p) = f(*std::get<0>(p));
      }

      /**
       * 将F应用于一个有三个元素的迭代器集合。
       *
       */
      template <typename I1, typename I2, typename I3>
      static void
      apply(const F &f, const std::tuple<I1, I2, I3> &p)
      {
        *std::get<2>(p) = f(*std::get<0>(p), *std::get<1>(p));
      }

      /**
       * 将F应用于一个有三个元素的迭代器集合。
       *
       */
      template <typename I1, typename I2, typename I3, typename I4>
      static void
      apply(const F &f, const std::tuple<I1, I2, I3, I4> &p)
      {
        *std::get<3>(p) = f(*std::get<0>(p), *std::get<1>(p), *std::get<2>(p));
      }
    };



    /**
     * 取一个函数对象并从中创建一个Body对象。我们在这个辅助函数中这样做，因为我们必须指定F的实际数据类型。
     *
     * - 对于函数对象来说，这往往是非常复杂的。
     *
     */
    template <typename F>
    Body<F>
    make_body(const F &f)
    {
      return Body<F>(f);
    }



#ifdef DEAL_II_WITH_TBB
    /**
     * 封装  tbb::parallel_for.  。
     *
     */
    template <typename Iterator, typename Functor>
    void
    parallel_for(Iterator           x_begin,
                 Iterator           x_end,
                 const Functor &    functor,
                 const unsigned int grainsize)
    {
      tbb::parallel_for(tbb::blocked_range<Iterator>(x_begin, x_end, grainsize),
                        functor,
                        tbb::auto_partitioner());
    }



    /**
     * 当提供一个affinite_partitioner时，封装 tbb::parallel_for 。
     *
     */
    template <typename Iterator, typename Functor>
    void
    parallel_for(Iterator                                          x_begin,
                 Iterator                                          x_end,
                 const Functor &                                   functor,
                 const unsigned int                                grainsize,
                 const std::shared_ptr<tbb::affinity_partitioner> &partitioner)
    {
      tbb::parallel_for(tbb::blocked_range<Iterator>(x_begin, x_end, grainsize),
                        functor,
                        *partitioner);
    }
#endif
  } // namespace internal

  /**
   * 一个执行<code>*out++ = predicate(*in++)</code>动作的算法，其中 <code>in</code> 迭代器在给定的输入范围内。    这个算法的作用与 std::transform 的作用差不多。不同的是，当deal.II被配置为使用多线程时，该函数可以并行运行。    如果是并行运行，迭代器范围会被分割成若干块，每块都被打包成一个任务，并交给线程积木调度器，在计算资源可用的情况下进行处理。一旦所有的块都被处理完毕，该函数就会返回。最后一个参数表示每个任务的迭代器范围的最小元素数；这个数字必须大到足以摊销新任务的启动成本，小到足以确保任务可以合理地负载平衡。    关于这个函数适用于哪类问题的讨论，见 @ref threads "多处理器的并行计算 "
   * 模块。
   *
   */
  template <typename InputIterator, typename OutputIterator, typename Predicate>
  void
  transform(const InputIterator &begin_in,
            const InputIterator &end_in,
            OutputIterator       out,
            const Predicate &    predicate,
            const unsigned int   grainsize)
  {
#ifndef DEAL_II_WITH_TBB
    // make sure we don't get compiler
    // warnings about unused arguments
    (void)grainsize;

    for (OutputIterator in = begin_in; in != end_in;)
      *out++ = predicate(*in++);
#else
    using Iterators     = std::tuple<InputIterator, OutputIterator>;
    using SyncIterators = SynchronousIterators<Iterators>;
    Iterators x_begin(begin_in, out);
    Iterators x_end(end_in, OutputIterator());
    internal::parallel_for(SyncIterators(x_begin),
                           SyncIterators(x_end),
                           internal::make_body(predicate),
                           grainsize);
#endif
  }



  /**
   * 一个执行<code>*out++ = predicate(*in1++, in2++)</code>动作的算法，其中 <code>in1</code> 迭代器在给定的输入范围内，使用tbb的并行for操作符。    这个算法几乎做了 std::transform 的工作。不同的是，当deal.II被配置为使用多线程时，该函数可以并行运行。    如果是并行运行，迭代器范围会被分割成若干块，每块都被打包成一个任务，并交给线程构件调度器，在计算资源可用的情况下进行处理。一旦所有的块都被处理完毕，该函数就会返回。最后一个参数表示每个任务的迭代器范围的最小元素数；这个数字必须大到足以摊销新任务的启动成本，小到足以确保任务可以合理地负载平衡。    关于这个函数适用于哪类问题的讨论，见 @ref threads "多处理器的并行计算 "
   * 模块。
   *
   */
  template <typename InputIterator1,
            typename InputIterator2,
            typename OutputIterator,
            typename Predicate>
  void
  transform(const InputIterator1 &begin_in1,
            const InputIterator1 &end_in1,
            InputIterator2        in2,
            OutputIterator        out,
            const Predicate &     predicate,
            const unsigned int    grainsize)
  {
#ifndef DEAL_II_WITH_TBB
    // make sure we don't get compiler
    // warnings about unused arguments
    (void)grainsize;

    for (OutputIterator in1 = begin_in1; in1 != end_in1;)
      *out++ = predicate(*in1++, *in2++);
#else
    using Iterators =
      std::tuple<InputIterator1, InputIterator2, OutputIterator>;
    using SyncIterators = SynchronousIterators<Iterators>;
    Iterators x_begin(begin_in1, in2, out);
    Iterators x_end(end_in1, InputIterator2(), OutputIterator());
    internal::parallel_for(SyncIterators(x_begin),
                           SyncIterators(x_end),
                           internal::make_body(predicate),
                           grainsize);
#endif
  }



  /**
   * 一个执行<code>*out++ = predicate(*in1++, in2++,in3++)</code>动作的算法，其中 <code>in1</code> 迭代器在给定的输入范围内。    这种算法的作用与 std::transform 的作用差不多。不同的是，当deal.II被配置为使用多线程时，该函数可以并行运行。    如果是并行运行，迭代器范围会被分割成若干块，每块都被打包成一个任务，并交给线程积木调度器，在计算资源可用时进行处理。一旦所有的块都被处理完毕，该函数就会返回。最后一个参数表示每个任务的迭代器范围的最小元素数；这个数字必须大到足以摊销新任务的启动成本，小到足以确保任务可以合理地负载平衡。    关于这个函数适用于哪类问题的讨论，见 @ref threads "多处理器的并行计算 "
   * 模块。
   *
   */
  template <typename InputIterator1,
            typename InputIterator2,
            typename InputIterator3,
            typename OutputIterator,
            typename Predicate>
  void
  transform(const InputIterator1 &begin_in1,
            const InputIterator1 &end_in1,
            InputIterator2        in2,
            InputIterator3        in3,
            OutputIterator        out,
            const Predicate &     predicate,
            const unsigned int    grainsize)
  {
#ifndef DEAL_II_WITH_TBB
    // make sure we don't get compiler
    // warnings about unused arguments
    (void)grainsize;

    for (OutputIterator in1 = begin_in1; in1 != end_in1;)
      *out++ = predicate(*in1++, *in2++, *in3++);
#else
    using Iterators = std::
      tuple<InputIterator1, InputIterator2, InputIterator3, OutputIterator>;
    using SyncIterators = SynchronousIterators<Iterators>;
    Iterators x_begin(begin_in1, in2, in3, out);
    Iterators x_end(end_in1,
                    InputIterator2(),
                    InputIterator3(),
                    OutputIterator());
    internal::parallel_for(SyncIterators(x_begin),
                           SyncIterators(x_end),
                           internal::make_body(predicate),
                           grainsize);
#endif
  }


  namespace internal
  {
#ifdef DEAL_II_WITH_TBB
    /**
     * 接受一个范围参数，并以其开始和结束来调用给定的函数。
     *
     */
    template <typename RangeType, typename Function>
    void
    apply_to_subranges(const tbb::blocked_range<RangeType> &range,
                       const Function &                     f)
    {
      f(range.begin(), range.end());
    }
#endif
  } // namespace internal


  /**
   * 这个函数将给定的函数参数 @p f 应用于范围
   * <code>[begin,end)</code>
   * 中的所有元素，并可能以并行方式进行。 step-69
   * 中给出了它的一个使用实例。
   * 然而，在许多情况下，在每个元素上调用一个函数并不高效，所以这个函数在子范围上调用给定的函数对象。
   * 换句话说：如果给定的范围  <code>[begin,end)</code>
   * 小于grainsize，或者多线程没有启用，那么我们就调用
   * <code>f(begin,end)</code>
   * ；否则，我们可能会执行，可能是%并行的，一连串的调用
   * <code>f(b,e)</code>  ，其中  <code>[b,e)</code> are subintervals of
   * <code>[begin,end)</code>  和我们对  <code>f(.,.)</code>
   * 的调用集合将发生在不相交的子区间，共同覆盖原始区间
   * <code>[begin,end)</code>  。
   * 很多时候，被调用的函数当然要获得额外的信息，例如，对于迭代器参数的给定值，要对哪个对象进行处理。这可以通过
   * <i>binding</i>
   * 某些参数来实现。例如，这里有一个对全矩阵 $A$
   * 和向量 $x,y$ 的矩阵-向量乘法 $y=Ax$ 的实现。
   * @code
   * void matrix_vector_product (const FullMatrix &A,
   *                             const Vector     &x,
   *                             Vector           &y)
   * {
   *   parallel::apply_to_subranges
   *      (0, A.n_rows(),
   *       [&](const unsigned int begin_row,
   *           const unsigned int end_row)
   *       {
   *         mat_vec_on_subranges(begin_row, end_row, A, x, y);
   *       },
   *       50);
   * }
   *
   * void mat_vec_on_subranges (const unsigned int begin_row,
   *                            const unsigned int end_row,
   *                            const FullMatrix &A,
   *                            const Vector     &x,
   *                            Vector           &y)
   * {
   *   for (unsigned int row=begin_row; row!=end_row; ++row)
   *     for (unsigned int col=0; col<x.size(); ++col)
   *       y(row) += A(row,col) x(col);
   * }
   * @endcode
   * 注意我们是如何使用lambda函数将 <code>mat_vec_on_subranges</code> 从一个需要5个参数的函数转换为一个需要2个参数的函数，并将其余参数绑定。由此产生的函数对象只需要两个参数，`begin_row'和`end_row'，而其他参数都是固定的。    如果在单线程模式下，该代码将在整个范围 <code>[0,n_rows)</code> 上调用 <code>mat_vec_on_subranges</code> 一次。然而，在多线程模式下，它可能会在这个区间的子范围内被多次调用，可能允许一个以上的CPU核心来处理部分工作。     @p grainsize 参数（在上面的例子中为50）确保子区间不会变得太小，以避免在调度子区间的CPU资源上花费更多的时间，而不是在做实际工作。    关于这个函数适用于哪类问题的讨论，也可以参见 @ref threads "多处理器的并行计算 "
   * 模块。
   *
   */
  template <typename RangeType, typename Function>
  void
  apply_to_subranges(const RangeType &                         begin,
                     const typename identity<RangeType>::type &end,
                     const Function &                          f,
                     const unsigned int                        grainsize)
  {
#ifndef DEAL_II_WITH_TBB
    // make sure we don't get compiler
    // warnings about unused arguments
    (void)grainsize;

    f(begin, end);
#else
    internal::parallel_for(begin,
                           end,
                           [&f](const tbb::blocked_range<RangeType> &range) {
                             internal::apply_to_subranges<RangeType, Function>(
                               range, f);
                           },
                           grainsize);
#endif
  }



  /**
   * 这是一个专门用于for循环的类，其固定范围由无符号整数给出。这是一个抽象的基类，一个实际的工作者函数是从它派生出来的。有一个公共函数apply可以并行地发出for循环，只要有足够的工作要做（即元素的数量大于grain_size），就把工作细分到可用的处理器核心上。在这个函数中，一个虚拟函数apply_to_subrange被调用，它指定了两个整数的范围<tt>[lower,
   * upper)</tt>，这需要在一个派生类中定义。
   * 这个类所涵盖的并行化情况是函数apply_to_subranges所能实现的一个子集（它也涵盖了可能不是由整数范围描述的更一般的迭代器的情况）。然而，对于简单的整数范围，人们可能更喜欢这个类，比如当有许多结构相似的循环时，例如对一个指针数组进行一些简单的复制或算术操作。在这种情况下，apply_to_subranges将产生大量的代码（或者说，大量的符号），因为它将由
   * std::bind
   * 产生的长名称传递给TBB中函数的模板化并行。这可以大大增加编译时间和目标代码的大小。同样地，不正确地使用
   * std::bind
   * 往往会导致非常隐晦的错误信息，这一点可以通过这个类来避免（只有一个虚拟函数需要在派生类中定义）。最后，在并行函数的背景下，虚拟函数的额外成本是可以忽略不计的。将工作实际发布到线程上的成本要高得多，而这又应该比for循环中的实际工作要少得多。
   *
   */
  struct ParallelForInteger
  {
    /**
     * 解构器。变成了虚拟的，以确保派生类也有虚拟的析构器。
     *
     */
    virtual ~ParallelForInteger() = default;

    /**
     * 这个函数在给定的范围内运行for循环<tt>[lower,upper)</tt>，当end-begin大于最小并行粒度时，可能是并行的。这个函数被标记为const，因为当几个线程同时处理相同的数据时，它的任何改变派生类数据的操作本质上都不是线程安全的。
     *
     */
    void
    apply_parallel(const std::size_t begin,
                   const std::size_t end,
                   const std::size_t minimum_parallel_grain_size) const;

    /**
     * 在派生类中定义的用于处理子范围的虚拟函数。
     * 这个函数被标记为const，因为当几个线程同时处理相同的数据时，任何改变派生类的数据的操作本质上都不是线程安全的。
     *
     */
    virtual void
    apply_to_subrange(const std::size_t, const std::size_t) const = 0;
  };



  namespace internal
  {
#ifdef DEAL_II_WITH_TBB
    /**
     * 一个符合TBB
     * parallel_reduce函数正文要求的类。第一个模板参数表示要进行还原的类型。第二个模板参数表示对每个子范围应调用的函数对象的类型。
     *
     */
    template <typename ResultType, typename Function>
    struct ReductionOnSubranges
    {
      /**
       * 一个变量，它将保存还原的结果。
       *
       */
      ResultType result;

      /**
       * 构造器。取出对每个子范围要调用的函数对象，以及相对于还原操作的中性元素。
       * 第二个参数表示一个函数对象，它将被用来把两个计算的结果还原成一个数字。如果我们想简单地累积整数结果，一个例子是
       * std::plus<int>().  。
       *
       */
      template <typename Reductor>
      ReductionOnSubranges(const Function & f,
                           const Reductor & reductor,
                           const ResultType neutral_element = ResultType())
        : result(neutral_element)
        , f(f)
        , neutral_element(neutral_element)
        , reductor(reductor)
      {}

      /**
       * 分割构造函数。关于这一点，请看TBB书中的更多细节。
       *
       */
      ReductionOnSubranges(const ReductionOnSubranges &r, tbb::split)
        : result(r.neutral_element)
        , f(r.f)
        , neutral_element(r.neutral_element)
        , reductor(r.reductor)
      {}

      /**
       * 连接操作：合并不同子区间上的计算结果。
       *
       */
      void
      join(const ReductionOnSubranges &r)
      {
        result = reductor(result, r.result);
      }

      /**
       * 在指定的范围内执行给定的函数。
       *
       */
      template <typename RangeType>
      void
      operator()(const tbb::blocked_range<RangeType> &range)
      {
        result = reductor(result, f(range.begin(), range.end()));
      }

    private:
      /**
       * 在每个子区间上调用的函数对象。
       *
       */
      const Function f;

      /**
       * 相对于还原操作的中性元素。在调用分割构造函数时需要这个，因为在这种情况下我们必须重新设置结果变量。
       *
       */
      const ResultType neutral_element;

      /**
       * 用来将两次调用的结果减少为一个数字的函数对象。
       *
       */
      const std::function<ResultType(ResultType, ResultType)> reductor;
    };
#endif
  } // namespace internal


  /**
   * 这个函数的工作原理很像apply_to_subranges()，但是它允许将在每个子范围内计算的数字结果累积为一个数字。
   * 这个数字的类型是由需要明确指定的ResultType模板参数给出的。
   * 使用这个函数的一个例子是计算一个正方形矩阵  $A$
   * 和一个向量  $x$  的表达式  $x^T A x$
   * 的值。对行的求和可以被并行化，整个代码可能看起来像这样。
   * @code
   * void matrix_norm (const FullMatrix &A,
   *                   const Vector     &x)
   * {
   *   return
   *    std::sqrt
   *     (parallel::accumulate_from_subranges<double>
   *      (0, A.n_rows(),
   *       [&](const unsigned int begin_row,
   *           const unsigned int end_row)
   *       {
   *         mat_vec_on_subranges(begin_row, end_row, A, x, y);
   *       },
   *       50);
   * }
   *
   * double
   * mat_norm_sqr_on_subranges (const unsigned int begin_row,
   *                            const unsigned int end_row,
   *                            const FullMatrix &A,
   *                            const Vector     &x)
   * {
   *   double norm_sqr = 0;
   *   for (unsigned int row=begin_row; row!=end_row; ++row)
   *     for (unsigned int col=0; col<x.size(); ++col)
   *       norm_sqr += x(row) A(row,col) x(col);
   *   return norm_sqr;
   * }
   * @endcode
   * 这里，如果 <code>mat_norm_sqr_on_subranges</code> 的范围小于最小粒度（上面选择的是50），或者如果deal.II被配置为不使用多线程，则在整个 <code>[0,A.n_rows())</code> 上调用。否则，它可以在给定范围的子集上调用，并在内部累积各个子集的结果。      @warning  如果ResultType是一个浮点类型，那么累加就不是一个关联操作。换句话说，如果给定的函数对象在三个子范围上被调用三次，返回值  $a,b,c$  ，那么这个函数的返回结果是  $(a+b)+c$  。然而，根据这三个子任务在可用CPU资源上的分布情况，结果也可能是 $(a+c)+b$ 或其他任何排列组合；因为浮点加法不是关联性的（当然，与实数%数的加法相反），多次调用这个函数的结果可能在舍入的顺序上有所不同。    关于这个函数适用于哪类问题的讨论，请参见 @ref threads "多处理器的并行计算 "
   * 模块。
   *
   */
  template <typename ResultType, typename RangeType, typename Function>
  ResultType
  accumulate_from_subranges(const Function &                          f,
                            const RangeType &                         begin,
                            const typename identity<RangeType>::type &end,
                            const unsigned int                        grainsize)
  {
#ifndef DEAL_II_WITH_TBB
    // make sure we don't get compiler
    // warnings about unused arguments
    (void)grainsize;

    return f(begin, end);
#else
    internal::ReductionOnSubranges<ResultType, Function> reductor(
      f, std::plus<ResultType>(), 0);
    tbb::parallel_reduce(tbb::blocked_range<RangeType>(begin, end, grainsize),
                         reductor,
                         tbb::auto_partitioner());
    return reductor.result;
#endif
  }


  // --------------------- for loop affinity partitioner -----------------------

  /**
   * 一个以线程安全的方式包装TBB亲和分区器的类。在Vector中，我们使用一个共享指针，在相同大小的不同向量之间共享一个亲和分区器，以提高数据（和NUMA）的定位性。然而，当一个外部任务进行多个向量操作时，共享指针可能导致竞赛条件。这个类只允许一个实例获得一个分区器。其他对象不能使用该对象，需要创建自己的副本。
   *
   */
  namespace internal
  {
    class TBBPartitioner
    {
    public:
      /**
       * 构造函数。
       *
       */
      TBBPartitioner();

#ifdef DEAL_II_WITH_TBB
      /**
       * 销毁器。检查该对象是否不再被使用，即所有的循环都已完成。
       *
       */
      ~TBBPartitioner();

      /**
       * 返回一个亲和分区器。如果类所拥有的分区器是自由的，它将在此返回。如果另一个线程还没有释放它，就会创建一个新的对象。要再次释放该分区器，请通过release_one_partitioner()调用返回它。
       *
       */
      std::shared_ptr<tbb::affinity_partitioner>
      acquire_one_partitioner();

      /**
       * 在通过acquisition_one_partitioner()在tbb循环中使用分区器后，这个调用使分区器再次可用。
       *
       */
      void
      release_one_partitioner(std::shared_ptr<tbb::affinity_partitioner> &p);

    private:
      /**
       * 存储的分区器，可以在多次运行 tbb::parallel_for
       * 中积累知识。
       *
       */
      std::shared_ptr<tbb::affinity_partitioner> my_partitioner;

      /**
       * 一个标志，表示该分区器是否已经获得但尚未释放，即它在其他地方使用。
       *
       */
      bool in_use;

      /**
       * 一个突变器，用于保护对in_use标志的访问。
       *
       */
      std::mutex mutex;
#endif
    };
  } // namespace internal
} // namespace parallel


namespace internal
{
  namespace VectorImplementation
  {
    /**
     * 如果我们对向量进行并行计算（例如，我们将两个向量相加得到第三个向量，并对所有元素进行并行循环），那么这个变量决定了最小的元素数，对于这个最小的元素数，将一个范围的元素再分割开来分配给不同的线程是有利的。
     * 这个变量可以作为一个全局可写变量，以便让测试组也能测试并行情况。默认情况下，它被设置为几千个元素，这是一个测试工具通常不会遇到的情况。因此，在testuite中，我们将其设置为一个
     *
     * - 一个巨大的无益的值，但肯定会测试并行操作。
     *
     */
    extern unsigned int minimum_parallel_grain_size;
  } // namespace VectorImplementation


  namespace SparseMatrixImplementation
  {
    /**
     * 就像 internal::VectorImplementation::minimum_parallel_grain_size,
     * 一样，但现在表示一个矩阵的行数，应该作为最低限度的工作。
     *
     */
    extern unsigned int minimum_parallel_grain_size;
  } // namespace SparseMatrixImplementation

} // end of namespace internal


 /* --------------------------- inline functions ------------------------- */ 

namespace parallel
{
#ifdef DEAL_II_WITH_TBB

  namespace internal
  {
    /**
     * 这是TBB为ParallelForInteger类实际调用的函数。
     *
     */
    struct ParallelForWrapper
    {
      ParallelForWrapper(const parallel::ParallelForInteger &worker)
        : worker_(worker)
      {}

      void
      operator()(const tbb::blocked_range<std::size_t> &range) const
      {
        worker_.apply_to_subrange(range.begin(), range.end());
      }

      const parallel::ParallelForInteger &worker_;
    };
  } // namespace internal

#endif


  inline void
  ParallelForInteger::apply_parallel(
    const std::size_t begin,
    const std::size_t end,
    const std::size_t minimum_parallel_grain_size) const
  {
#ifndef DEAL_II_WITH_TBB
    // make sure we don't get compiler
    // warnings about unused arguments
    (void)minimum_parallel_grain_size;

    apply_to_subrange(begin, end);
#else
    internal::ParallelForWrapper worker(*this);
    internal::parallel_for(begin, end, worker, minimum_parallel_grain_size);
#endif
  }

} // end of namespace parallel

DEAL_II_NAMESPACE_CLOSE

#endif


