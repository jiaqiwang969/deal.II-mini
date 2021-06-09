//include/deal.II-translator/base/thread_management_0.txt
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

#ifndef dealii_thread_management_h
#  define dealii_thread_management_h


#  include <deal.II/base/config.h>

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/multithread_info.h>
#  include <deal.II/base/std_cxx17/tuple.h>
#  include <deal.II/base/template_constraints.h>

#  include <atomic>
#  include <condition_variable>
#  include <functional>
#  include <future>
#  include <iterator>
#  include <list>
#  include <memory>
#  include <mutex>
#  include <thread>
#  include <tuple>
#  include <utility>
#  include <vector>



DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup threads */ 
 /*@{*/ 


/**
 * 一个用于实现deal.II中线程管理的命名空间。这个命名空间的大部分内容在deal.II的文档页面上链接到的一份报告中详细讨论过。
 *
 *
 * @ingroup threads
 *
 *
 */
namespace Threads
{
  /**
   * 一个实现<a
   * href="https://en.wikipedia.org/wiki/Lock_(computer_science)">mutex</a>的类。
   * Mutexes用于锁定数据结构，以确保在同一时间内只有一个执行线程可以访问它们。
   * 这个类是对 `std::mutex`.
   * 的一个薄包装，唯一的区别是这个类是可复制的，而
   * `std::mutex` 是不可复制的。
   * 事实上，当复制时，接收对象不会从被复制的对象中复制任何状态，也就是说，会创建一个全新的mutex。如果一个突变体被用作成员变量来锁定一个类的其他成员变量，那么这些语义与常见的使用情况是一致的：在这种情况下，被复制到对象的突变体应该只保护被复制到对象的成员，而不是被复制到对象和被复制从对象的成员。因为在类被复制的时候，目的地的成员变量还没有被使用，其对应的mutex也应该保持在原来的状态。
   *
   */
  class Mutex : public std::mutex
  {
  public:
    /**
     * 默认构造函数。
     *
     */
    Mutex() = default;

    /**
     * 复制构造函数。正如在这个类的文档中所讨论的，不会从作为参数给出的对象中复制状态。
     *
     */
    Mutex(const Mutex &)
      : std::mutex()
    {}

    /**
     * 拷贝操作符。正如在这个类的文档中所讨论的，不从作为参数的对象中复制状态。
     *
     */
    Mutex &
    operator=(const Mutex &)
    {
      return *this;
    }
  };
} // namespace Threads


namespace Threads
{
  /**
   * 将范围 <code>[begin,end)</code> into <code>n_intervals</code>
   * 分割为大小相等的子区间。如果整个范围内的元素数不正好被
   * <code>n_intervals</code>
   * 整除，最后一个区间会稍微大一点。迭代器的类型必须满足前向迭代器的要求，即
   * <code>operator++</code> 必须可用，当然也必须是可分配的。
   * 一个子间隔的列表以一对迭代器的向量形式返回，其中每一对迭代器表示的范围是
   * <code>[begin[i],end[i])</code>  .
   * @ingroup threads
   *
   */
  template <typename ForwardIterator>
  std::vector<std::pair<ForwardIterator, ForwardIterator>>
  split_range(const ForwardIterator &begin,
              const ForwardIterator &end,
              const unsigned int     n_intervals);

  /**
   * 将区间 <code>[begin,end)</code>
   * 分割成（几乎）同等大小的子区间。这个函数的工作原理与之前的函数基本相同，不同的是，现在取的是定义整个区间的值，而不是迭代器。
   * @ingroup threads
   *
   */
  std::vector<std::pair<unsigned int, unsigned int>>
  split_interval(const unsigned int begin,
                 const unsigned int end,
                 const unsigned int n_intervals);

  /**
   * @cond 内部
   *
   */

  /**
   * 一个命名空间，用于实现线程子系统的辅助函数和类似功能。这个命名空间的成员不打算公开使用。
   *
   */
  namespace internal
  {
    /**
     * @internal
     * 如果在子线程中抛出了一个异常，它不会被传播到主线程中。因此，由应用程序的主函数或其他一些部分提供的异常处理程序将无法捕捉这些异常。因此，我们必须在每个子线程的顶层函数中提供一个异常处理程序，至少要捕捉到异常并打印出一些信息，而不是让操作系统在没有信息的情况下直接杀死程序。因此，在我们用作新线程入口的每个函数中，我们都安装了一个try-catch块，如果捕获到
     * <code>std::exception</code>
     * 类型的异常，它就会将控制权移交给这个函数，然后它将提供一些输出。
     *
     */
    [[noreturn]] void
    handle_std_exception(const std::exception &exc);

    /**
     * @internal  和上面一样，但是异常的类型不是从
     * <code>std::exception</code>
     * 派生出来的，所以没有什么办法提供更有用的东西。
     *
     */
    [[noreturn]] void
    handle_unknown_exception();
  } // namespace internal

  /**
   * @endcond
   *
   */

} // namespace Threads

 /* ----------- implementation of functions in namespace Threads ---------- */ 
#  ifndef DOXYGEN
namespace Threads
{
  template <typename ForwardIterator>
  std::vector<std::pair<ForwardIterator, ForwardIterator>>
  split_range(const ForwardIterator &begin,
              const ForwardIterator &end,
              const unsigned int     n_intervals)
  {
    using IteratorPair = std::pair<ForwardIterator, ForwardIterator>;

    // in non-multithreaded mode, we often have the case that this
    // function is called with n_intervals==1, so have a shortcut here
    // to handle that case efficiently

    if (n_intervals == 1)
      return (std::vector<IteratorPair>(1, IteratorPair(begin, end)));

    // if more than one interval requested, do the full work
    const unsigned int n_elements              = std::distance(begin, end);
    const unsigned int n_elements_per_interval = n_elements / n_intervals;
    const unsigned int residual                = n_elements % n_intervals;

    std::vector<IteratorPair> return_values(n_intervals);

    return_values[0].first = begin;
    for (unsigned int i = 0; i < n_intervals; ++i)
      {
        if (i != n_intervals - 1)
          {
            return_values[i].second = return_values[i].first;
            // note: the cast is performed to avoid a warning of gcc
            // that in the library `dist>=0' is checked (dist has a
            // template type, which here is unsigned if no cast is
            // performed)
            std::advance(return_values[i].second,
                         static_cast<signed int>(n_elements_per_interval));
            // distribute residual in division equally among the first
            // few subintervals
            if (i < residual)
              ++return_values[i].second;

            return_values[i + 1].first = return_values[i].second;
          }
        else
          return_values[i].second = end;
      }
    return return_values;
  }
} // namespace Threads

#  endif // DOXYGEN

namespace Threads
{
  namespace internal
  {
    /**
     * @internal
     * 给定一个任意类型的RT，存储其中的一个元素，并通过函数get()和set()授予对它的访问。对于引用类型（需要作为指向被引用对象的指针来存储）和void类型，有专门的规定。
     * 这个函数与 `std::promise`/`std::future`
     * 的类的组合没有什么不同。不同的是，一个
     * `std::promise` 只能通过 `std::future::get()`
     * 读取一次（据推测，这种设计是由于 `std::future::get()`
     * 可以抛出一个先前存储在 `std::promise`). 中的异常
     * 另一方面，这个类使结果可用于任意次数。它也不存储任何异常（尽管它们会被使用当前类的类转发）。
     *
     */
    template <typename RT>
    struct return_value
    {
    private:
      RT value;

    public:
      using reference_type = RT &;

      inline return_value()
        : value()
      {}

      inline reference_type
      get()
      {
        return value;
      }

      inline void
      set(RT &&v)
      {
        value = std::move(v);
      }

      inline void
      set_from(std::future<RT> &v)
      {
        value = std::move(v.get());
      }
    };


    /**
     * @internal
     * 给定一个任意类型的RT，存储其中的一个元素，并通过函数get()和set()授予其访问权。这是对引用类型的特殊化：由于引用在构造时间之后不能被设置，我们存储一个指针，而指针持有被引用对象的地址。
     * 这个函数与 `std::promise`/`std::future`
     * 类的组合没有什么不同。不同的是，一个 `std::promise`
     * 只能通过 `std::future::get()`
     * 读取一次（据推测，这种设计是由于 `std::future::get()`
     * 可以抛出一个先前存储在 `std::promise`). 中的异常
     * 另一方面，这个类使结果可用于任意次数。它也不存储任何异常（尽管它们会被使用当前类的类转发）。
     *
     */
    template <typename RT>
    struct return_value<RT &>
    {
    private:
      RT *value;

    public:
      using reference_type = RT &;

      inline return_value()
        : value(nullptr)
      {}

      inline reference_type
      get() const
      {
        return *value;
      }

      inline void
      set(RT &v)
      {
        value = &v;
      }

      inline void
      set_from(std::future<RT &> &v)
      {
        value = &v.get();
      }
    };


    /**
     * @internal
     * 给定一个任意类型的RT，存储其中的一个元素，并通过函数get()和set()授予其访问权。这是对void类型的特殊化：显然没有任何东西可以存储，所以没有函数set()，而函数get()则返回void。
     * 这个函数与 `std::promise`/`std::future`
     * 类的组合没有什么不同。不同的是，一个 `std::promise`
     * 只能通过 `std::future::get()`
     * 读取一次（据推测，这种设计是由于 `std::future::get()`
     * 可以抛出一个先前存储在 `std::promise`). 中的异常
     * 另一方面，这个类使结果可用于任意次数。它也不存储任何异常（尽管它们会被使用当前类的类转发）。
     *
     */
    template <>
    struct return_value<void>
    {
      using reference_type = void;

      static inline void
      get()
      {}


      inline void
      set_from(std::future<void> &)
      {}
    };
  } // namespace internal



  namespace internal
  {
    template <typename RT>
    inline void
    call(const std::function<RT()> & function,
         internal::return_value<RT> &ret_val)
    {
      ret_val.set(function());
    }


    inline void
    call(const std::function<void()> &function, internal::return_value<void> &)
    {
      function();
    }
  } // namespace internal



  namespace internal
  {
    /**
     * 一个代表线程的类。对于每个线程，我们正好创建一个这样的对象
     *
     * - 正是一个，因为它携带着线程上调用的函数的返回值。        虽然我们每个线程只有一个这样的对象，但几个 Threads::Thread 对象可以引用这个描述符。如果所有的Thread对象都超出了范围，那么ThreadDescriptor将在被销毁之前从线程中分离出来。
     *
     */
    template <typename RT>
    struct ThreadDescriptor
    {
      /**
       * 一个代表线程启动的对象。
       *
       */
      std::thread thread;

      /**
       * 一个对象，它将保存线程上调用的函数所返回的值。
       * 返回值被保存在一个shared_ptr中，因为我们可能会放弃ThreadDescriptor。
       * 这可以确保该对象保持活力，直到线程退出。
       *
       */
      std::shared_ptr<return_value<RT>> ret_val;

      /**
       * 一个原子型的bool变量，最初是假的，当一个新的线程开始时被设置为真，一旦join()被调用，就被设置为假。
       * 我们使用这个变量来确保我们可以在同一个线程上调用两次join()。出于某种原因，如果试图调用
       * std::thread::join 两次，C++标准库会抛出一个
       * std::system_error 异常（事实上，在第二次调用之前，
       * std::thread::joinable
       * 返回false），但这是一个有点可取的做法，因为我们不必跟踪
       * join() 是否在之前被调用。
       * 使用这个变量，只要我们之前调用过join()，这个变量就会被设置为true，我们就可以跳过第二次调用
       * std::thread::join()
       * 。对这个变量的访问是由下面的突变来保护的。
       * @note
       * 历史上，我们不需要这个变量的mutex：线程只能从最初创建它的线程中被加入。因此，在不创建线程的函数中发生的一切（比如下面的join()函数）在外界看来都是原子的。由于我们在调用
       * std::thread::join,
       * 的同一个函数中清除和测试thread_is_active，这些动作是原子的，不需要mutex。当然，两个线程可以同时在同一个线程对象上调用join()，但这个动作无论如何都是未定义的，因为它们不能同时加入同一个线程。也就是说，最近的C++标准似乎不再有这样的要求，即唯一可以调用join()的线程是创建该线程的那个。`pthread_join`似乎也没有这个要求。
       * 因此，我们实际上可以从不同的线程加入，我们在base/thread_validity_07中测试了这一点。
       * @note 我们需要使用 std::atomic<bool> 的原因在
       * Task::task_has_finished. 的文档中有详细讨论。
       *
       */
      std::atomic<bool> thread_is_active;

      /**
       * 守护对前一个变量的访问的Mutex。
       *
       */
      Mutex thread_is_active_mutex;

      /**
       * 默认构造函数。
       *
       */
      ThreadDescriptor()
        : thread_is_active(false)
      {}

      ~ThreadDescriptor()
      {
        if (!thread_is_active)
          return;
        thread.detach();
        thread_is_active = false;
      }

      /**
       * 启动线程并让它将其返回值放入ret_val对象中。
       *
       */
      void
      start(const std::function<RT()> &function)
      {
        thread_is_active = true;
        ret_val          = std::make_shared<return_value<RT>>();
        thread           = std::thread(thread_entry_point, function, ret_val);
      }


      /**
       * 等待线程结束。
       *
       */
      void
      join()
      {
        // see if the thread hasn't been joined yet. if it has, then
        // join() is a no-op. use schmidt's double-checking strategy
        // to use the mutex only when necessary
        if (thread_is_active == false)
          return;

        std::lock_guard<std::mutex> lock(thread_is_active_mutex);
        if (thread_is_active == true)
          {
            Assert(thread.joinable(), ExcInternalError());
            thread.join();
            thread_is_active = false;
          }
      }

    private:
      /**
       * 在该线程上运行的函数。
       *
       */
      static void
      thread_entry_point(const std::function<RT()> &       function,
                         std::shared_ptr<return_value<RT>> ret_val)
      {
        // call the function in question. since an exception that is
        // thrown from one of the called functions will not propagate
        // to the main thread, it will kill the program if not treated
        // here before we return to the operating system's thread
        // library
        try
          {
            call(function, *ret_val);
          }
        catch (const std::exception &exc)
          {
            internal::handle_std_exception(exc);
          }
        catch (...)
          {
            internal::handle_unknown_exception();
          }
      }
    };
  } // namespace internal


  /**
   * 一个代表被催生的线程的对象。这个对象可以在用户空间中自由复制，所有的实例将代表同一个线程，并可以要求等待它的终止和访问它的返回值。
   * 线程可以被放弃，也就是说，如果你只是调用
   * Threads::new_thread
   * ，但并不关心返回的对象，或者如果你把返回的
   * Threads::Thread
   * 对象分配给一个随后超出范围的对象，那么之前创建的线程仍然会继续工作。你只是不能再访问它的返回值，而且还可能发生你的程序在该线程完成工作之前就终止了。
   * 模板参数的默认值是  <code>void</code>
   * ，所以如果你在一个新线程上调用的函数没有返回值，你可以省略模板参数。
   * @ingroup threads   @deprecated  使用  std::thread  或  std::jthread
   * 来代替。
   * @note
   * 由于该类在ThreadGroup中使用，它的构造函数，而不是该类本身，已被弃用，以允许用以下方式进行编译
   *
   *
   *
   *
   *
   *
   * - error=deprecated-declarations.
   *
   */
  template <typename RT = void>
  class Thread
  {
  public:
    /**
     * 用一个函数对象构造一个线程对象。
     *
     */
    DEAL_II_DEPRECATED
    Thread(const std::function<RT()> &function)
      : thread_descriptor(new internal::ThreadDescriptor<RT>())
    {
      // in a second step, start the thread.
      thread_descriptor->start(function);
    }

    /**
     * 默认的构造函数。除了给它分配一个持有new_thread()函数创建的数据的线程对象外，你不能对这样构造的线程对象做什么。
     *
     */
    DEAL_II_DEPRECATED
    Thread() = default;

    /**
     * 复制构造函数。
     *
     */
    DEAL_II_DEPRECATED
    Thread(const Thread<RT> &t)
      : thread_descriptor(t.thread_descriptor)
    {}

    /**
     * 加入这个对象所代表的线程，也就是等待它完成。
     * 如果你使用了这个类的默认构造函数，并且没有给它分配一个线程对象，那么这个函数就是一个无用功。
     *
     */
    void
    join() const
    {
      if (thread_descriptor)
        thread_descriptor->join();
    }

    /**
     * 获取线程的函数的返回值。由于只有在线程结束后才能使用，这个函数在内部也调用join()。只要对象指的是同一个任务，你就可以多次调用这个函数，并期望每次都能得到相同的返回值。(返回的对象被移动的情况除外，见下文)。
     * @note  函数返回一个<i>non-@p const
     * reference</i>给返回对象，而不是返回对象。这允许编写诸如以下的代码
     * @code
     * Threads::Thread<int> t = Threads::new_thread (
     *   ...function returning an int...);
     * t.return_value() = 42;      // overwrite returned value
     * int i = t.return_value();   // i is now 42
     * @endcode
     * 你很少会有写这种代码的需要。另一方面，该函数需要返回一个可写的（非
     * @p const) 引用，以支持像这样的代码。
     * @code
     * std::unique_ptr<int> create_int (const std::string &s)
     * {
     *   ...
     * }
     *
     * void f()
     * {
     *   Threads::Thread<std::unique_ptr<int>>
     *     t = Threads::new_thread (&create_int, "42");
     *
     *   std::unique_ptr<int> i = std::move(t.return_value());
     *   ...
     * }
     * @endcode
     * 这里，需要对 `std::move` 返回的对象（即
     * <code>std::unique_ptr</code> 对象），因为
     * <code>std::unique_ptr</code>
     * 对象不能被复制。换句话说，要想从线程返回的对象中得到指针，就需要移动它，而为了移动它，当前函数需要返回一个可写的（非
     * @p const) 引用。
     *
     */
    typename internal::return_value<RT>::reference_type
    return_value()
    {
      join();
      return thread_descriptor->ret_val->get();
    }

    /**
     * 如果这个对象曾经有一个线程与之相关联，无论是通过使用非默认的构造函数还是通过赋值，则返回true。
     *
     */
    bool
    valid() const
    {
      return static_cast<bool>(thread_descriptor);
    }


    /**
     * 检查线程对象的平等性。由于这个类的对象存储了一个隐含的指向对象的指针，而这个指针对每个线程来说只存在一次，所以检查只是比较这些指针。
     *
     */
    bool
    operator==(const Thread &t) const
    {
      return thread_descriptor == t.thread_descriptor;
    }

  private:
    /**
     * 代表线程的对象的共享指针，并抽象出操作系统的函数来对其工作。这也确保了只要有至少一个订阅者，该对象就会存在。
     *
     */
    std::shared_ptr<internal::ThreadDescriptor<RT>> thread_descriptor;
  };


  namespace internal
  {
    /**
     * 一个一般的模板，如果t是引用类型，返回 std::ref(t)
     * ，否则返回t。
     * t是引用类型的情况在下面声明的部分特殊化中处理。
     *
     */
    template <typename T>
    struct maybe_make_ref
    {
      static T
      act(T &t)
      {
        return t;
      }
    };



    /**
     * 一个一般的模板，如果t是引用类型，返回 std::ref(t)
     * ，否则返回t。
     * t是引用类型的情况在这个部分特殊化中处理。
     *
     */
    template <typename T>
    struct maybe_make_ref<T &>
    {
      static std::reference_wrapper<T>
      act(T &t)
      {
        return std::ref(t);
      }
    };
  } // namespace internal



  // ----------- thread starters for functions not taking any parameters

  /**
   * new_thread函数的重载，适用于可以转换为 std::function<RT
   * ()>的对象，即任何可以像函数对象一样被调用而没有参数并返回RT（或void）类型的对象。
   * @ingroup threads
   *
   */
  template <typename RT>
  inline Thread<RT>
  new_thread(const std::function<RT()> &function)
  {
    return Thread<RT>(function);
  }



  /**
   * new_thread()函数的重载，适用于可以像函数对象一样被调用而没有参数的对象。特别是，这个函数允许用使用
   * std::bind, 产生的对象或使用lambda函数来调用
   * Threads::new_thread()
   * 。例如，在编写如下代码时，可以调用这个函数
   * @code
   * Threads::Thread<int>
   * thread = Threads::new_thread ( [] () {
   *                                        do_this();
   *                                        then_do_that();
   *                                        return 42;
   *                                      });
   * @endcode
   * 这里，我们在一个单独的线程上运行函数序列
   * <code>do_this()</code> and <code>then_do_that()</code>
   * ，通过使这里声明的lambda函数成为在该线程上执行的函数。然后lambda函数返回42（在这里有点无意义，但它当然可以是一些计算出来的数字），这将是你后来在线程（即lambda函数的主体）完成后可以通过
   * <code>thread.return_value()</code> 检索到的返回值。
   * @note
   * 每个lambda函数（或者你在这里传递给new_thread()函数的其他东西，例如一个
   * std::bind()
   * 表达式的结果）都有一个返回类型，因此会返回一个该类型的对象。这个类型可以通过在这个函数的声明中使用的C++11
   * <code>decltype</code>
   * 语句来推断，然后它被用作当前函数返回的 Threads::Thread
   * 对象的模板参数。
   * 在上面的例子中，因为lambda函数返回42（在C++中它的数据类型是
   * <code>int</code> ），推断的类型是 <code>int</code>
   * ，任务对象的类型将是 <code>Task@<int@></code>
   * 。换句话说，在用户代码中不能<i>necessary</i>明确指定λ或
   * std::bind
   * 表达式的返回类型，尽管可以通过（完全等价）写出明确的规定
   * @code
   * Threads::Thread<int>
   *   thread = Threads::new_thread ( [] ()
   *
   * -> int {
   *                                                 do_this();
   *                                                 then_do_that();
   *                                                 return 42;
   *                                               });
   * @endcode
   *
   * @note
   * 在实践中，你将传递给new_thread()的lambda函数当然会更加复杂。
   * 特别是，它们可能会从周围的环境中<i>capture</i>变量，并在lambda中使用它们。
   * 更多关于lambda函数的工作原理，请参见https://en.wikipedia.org/wiki/Anonymous_function#C.2B.2B_.28since_C.2B.2B11.29。
   * @note
   * 如果你将一个lambda函数作为参数传递给当前函数，该函数捕获了一个变量<i>by
   * reference</i>，或者如果你使用一个 std::bind
   * 将一个函数参数与一个引用变量绑定，使用 std::ref() 或
   * std::cref(),
   * ，那么显然你只能在你引用或捕获的变量具有至少延伸到线程结束时的生命周期时才能这样做。
   * @ingroup CPP11
   *
   */
  template <typename FunctionObjectType>
  inline auto
  new_thread(FunctionObjectType function_object)
    -> Thread<decltype(function_object())>
  {
    using return_type = decltype(function_object());
    return Thread<return_type>(std::function<return_type()>(function_object));
  }



  /**
   * 对非成员或静态成员函数的new_thread函数进行重载。
   * @ingroup threads
   *
   */
  template <typename RT, typename... Args>
  inline Thread<RT>
  new_thread(RT (*fun_ptr)(Args...), typename identity<Args>::type... args)
  {
    auto dummy = std::make_tuple(internal::maybe_make_ref<Args>::act(args)...);
    return new_thread(
      [dummy, fun_ptr]() -> RT { return std_cxx17::apply(fun_ptr, dummy); });
  }



  /**
   * 对成员函数的非const new_thread函数的重载。
   * @ingroup threads
   *
   */
  template <typename RT, typename C, typename... Args>
  inline Thread<RT>
  new_thread(RT (C::*fun_ptr)(Args...),
             typename identity<C>::type &c,
             typename identity<Args>::type... args)
  {
    // NOLINTNEXTLINE(modernize-avoid-bind) silence clang-tidy
    return new_thread(std::function<RT()>(std::bind(
      fun_ptr, std::ref(c), internal::maybe_make_ref<Args>::act(args)...)));
  }

  /**
   * 对常量成员函数new_thread函数的重载。
   * @ingroup threads
   *
   */
  template <typename RT, typename C, typename... Args>
  inline Thread<RT>
  new_thread(RT (C::*fun_ptr)(Args...) const,
             typename identity<const C>::type &c,
             typename identity<Args>::type... args)
  {
    // NOLINTNEXTLINE(modernize-avoid-bind) silence clang-tidy
    return new_thread(std::function<RT()>(std::bind(
      fun_ptr, std::cref(c), internal::maybe_make_ref<Args>::act(args)...)));
  }

  // ------------------------ ThreadGroup -------------------------------------

  /**
   * 一个线程对象的容器。允许添加新的线程对象并一起等待它们。线程对象需要对被调用的函数有相同的返回值。
   * @ingroup threads   @deprecated  用任务组代替。
   *
   */
  template <typename RT = void>
  class DEAL_II_DEPRECATED ThreadGroup
  {
  public:
    /**
     * 在集合中添加另一个线程对象。
     *
     */
    ThreadGroup &
    operator+=(const Thread<RT> &t)
    {
      threads.push_back(t);
      return *this;
    }

    /**
     * 等待集合中的所有线程都完成。如果有些线程已经被等待过了，这也不是什么问题，也就是说，你可以多次调用这个函数，如果你愿意，你也可以在后续调用这个函数之间添加新的线程对象。
     *
     */
    void
    join_all() const
    {
      for (typename std::list<Thread<RT>>::const_iterator t = threads.begin();
           t != threads.end();
           ++t)
        t->join();
    }

  private:
    /**
     * 线程对象的列表。
     *
     */
    std::list<Thread<RT>> threads;
  };


  namespace internal
  {
    /**
     * 通过评估动作来设置一个 std::promise 对象的值。
     *
     */
    template <typename RT, typename Function>
    void
    evaluate_and_set_promise(Function &function, std::promise<RT> &promise)
    {
      promise.set_value(function());
    }


    /**
     * 通过评估动作来设置一个 std::promise
     * 对象的值。这个函数是前一个函数的特殊化，用于返回类型为`void`的情况。因此，我们不能设置一个值。但我们确实评估了函数对象，并在没有参数的情况下调用
     * `std::promise::set_value()` 。
     *
     */
    template <typename Function>
    void
    evaluate_and_set_promise(Function &function, std::promise<void> &promise)
    {
      function();
      promise.set_value();
    }
  } // namespace internal



  /**
   * 这个类描述了一个任务对象，即通过调用
   * Threads::new_task(). 得到的东西。其想法是，
   * Threads::new_task()
   * 允许人们在C++运行时系统认为方便时运行一个函数。
   *
   * - 通常，当有一个空闲的处理器可用时。这可以用来在没有立即需要结果的情况下在后台运行，或者在有其他事情可以并行进行的情况下。  每当需要该背景任务的结果时，可以调用join()来等待任务完成，或者调用return_value()来获得在该背景任务上运行的函数所返回的值。    这个类在概念上类似于由 [`std::async`](https://en.cppreference.com/w/cpp/thread/async) 返回的 [`std::future`](https://en.cppreference.com/w/cpp/thread/future) 类（其本身类似于 Threads::new_task() 的作用）。主要的概念差异是，人们只能调用 `std::future::get()` 一次，而可以根据需要多次调用 Threads::Task::return_value() 。因此，它可以与 [`std::shared_future`](https://en.cppreference.com/w/cpp/thread/shared_future) 类相媲美。然而， `std::shared_future` 不能用于不能被复制的类型
   *
   * --例如对 `std::unique_ptr`, 的一个特殊限制。
   * @ingroup threads
   *
   */
  template <typename RT = void>
  class Task
  {
  public:
    /**
     * 构建一个任务对象，给定一个要在任务上执行的函数对象，然后安排这个函数的执行。然而，当
     * MultithreadInfo::n_threads()
     * 返回1时，即如果deal.II运行时系统已被配置为只使用一个线程，那么只要执行给定的函数对象即可。
     * @post
     * 使用这个构造函数会自动使任务对象可加入（）。
     *
     */
    Task(const std::function<RT()> &function_object)
    {
      if (MultithreadInfo::n_threads() > 1)
        task_data = std::make_shared<TaskData>(
          std::async(std::launch::async, function_object));
      else
        {
          // Only one thread allowed. So let the task run to completion
          // and just emplace a 'ready' future.
          //
          // The design of std::promise/std::future is unclear, but it
          // seems that the intent is to obtain the std::future before
          // we set the std::promise. So create the TaskData object at
          // the top and then run the task and set the returned
          // value. Since everything here happens sequentially, it
          // really doesn't matter in which order all of this is
          // happening.
          std::promise<RT> promise;
          task_data = std::make_shared<TaskData>(promise.get_future());
          try
            {
              internal::evaluate_and_set_promise(function_object, promise);
            }
          catch (...)
            {
              try
                {
                  // store anything thrown in the promise
                  promise.set_exception(std::current_exception());
                }
              catch (...)
                {}
              // set_exception() may throw too
            }
        }
    }

    /**
     * 默认构造函数。你不能对这样构造的任务对象做很多事情，除了给它分配一个持有由
     * Threads::new_task() 函数创建的数据的任务对象。
     * @post
     * 使用这个构造函数会使对象处于不可连接的状态，即joinable()将返回false。
     *
     */
    Task() = default;

    /**
     * 加入这个对象所代表的任务，即等待它完成。
     * 一个任务可以被多次加入（虽然第一次join()操作可能会阻塞，直到任务完成运行，但所有连续的加入尝试将立即返回）。
     * 如果在这个对象被初始化的任务上执行的操作抛出一个异常，而不是定期返回，那么调用当前的join()函数将首先等待该任务完成，然后反过来抛出该任务操作最初抛出的异常。这允许将异常从独立线程上执行的任务传播到调用线程。
     * (这种行为与
     * [`std::future`](https://en.cppreference.com/w/cpp/thread/future),
     * 不同， `std::future::wait()`
     * 函数只等待操作的完成，而只有当人们调用
     * `std::future::get()`. 时才会传播异常。然而，当把 "void
     * "函数放到独立的任务上时，这种做法就很尴尬了，因为这些任务实际上不会返回任何东西；因此，对这类任务调用
     * `std::task::wait()` 比调用 `std::task::get()`
     * 函数更自然，因为后者实际上不会返回任何可以得到的东西。)
     * @pre
     * 如果你使用了这个类的默认构造函数，并且没有给它分配一个任务对象，你就不能调用这个函数。换句话说，函数joinable()必须返回true。
     *
     */
    void
    join() const
    {
      // Make sure we actually have a task that we can wait for.
      AssertThrow(joinable(), ExcNoTask());

      task_data->wait();
    }

    /**
     * 返回当前对象是否可以被加入。一旦一个任务（通常是用
     * Threads::new_task())
     * 创建的）实际上已经被分配给它，你就可以加入一个任务对象。另一方面，如果该对象已被默认构建，该函数返回false。
     * 一个任务可以被多次加入（虽然第一次join()操作可能会阻塞，直到任务完成运行，但所有连续的加入尝试将立即返回）。因此，如果这个函数返回真，它将继续返回真，直到它所报告的任务对象被从另一个对象分配到。
     *
     */
    bool
    joinable() const
    {
      return (task_data != nullptr);
    }


    /**
     * 获取任务的函数的返回值。由于它只有在线程结束后才能使用，这个函数在内部也调用join()。只要对象指的是同一个任务，你就可以多次调用这个函数，并期望每次都能得到相同的返回值。(返回的对象被移动的情况除外，见下文)。
     * @note  函数返回一个<i>non-@p const
     * reference</i>给返回对象，而不是返回对象。这允许编写诸如以下的代码
     * @code
     * Threads::Task<int> t = Threads::new_task (...function returning an
     * int...); t.return_value() = 42;      // overwrite returned value int i =
     * t.return_value();   // i is now 42
     * @endcode
     * 你很少会有写这种代码的需要。另一方面，该函数需要返回一个可写的（非
     * @p const) 引用，以支持像这样的代码。
     * @code
     * std::unique_ptr<int> create_int (const std::string &s) { ... }
     *
     * void f()
     * {
     *   Threads::Task<std::unique_ptr<int>>
     *     t = Threads::new_task (&create_int, "42");
     *
     *   std::unique_ptr<int> i = std::move(t.return_value());
     *   ...
     * }
     * @endcode
     * 这里，需要对 `std::move` 返回的对象（即
     * <code>std::unique_ptr</code> 对象），因为
     * <code>std::unique_ptr</code>
     * 对象不能被复制。换句话说，要想从任务返回的对象中得到指针，需要移动它，为了移动它，当前函数需要返回一个可写的（非
     * @p const) 引用。        这个函数在内部调用了 join()
     * 成员函数。因此，正如那里所解释的，如果被打包的任务抛出一个异常，那么这个异常会被join()函数重新抛出，如果你之前没有调用join()，那么也会被当前函数抛出。
     * @pre
     * 如果你使用了这个类的默认构造函数，并且没有给它分配一个任务对象，你就不能调用这个函数。换句话说，函数joinable()必须返回true。
     *
     */
    typename internal::return_value<RT>::reference_type
    return_value()
    {
      // Make sure we actually have a task that we can wait for.
      AssertThrow(joinable(), ExcNoTask());

      // Then return the promised object. If necessary, wait for the promise to
      // be set.
      return task_data->get();
    }


    /**
     * @addtogroup  Exceptions  
     * @{ 
     *
     */

    /**
     * 异常情况
     *
     */
    DeclExceptionMsg(ExcNoTask,
                     "The current object is not associated with a task that "
                     "can be joined. It may have been detached, or you "
                     "may have already joined it in the past.");
    //@}
  private:
    /**
     * 一个持有 std::future
     * 的数据结构，任务将其返回值存入其中。由于只能调用
     * std::future::get()
     * 一次，我们在get()成员函数中这样做，然后将返回的对象移到`returned_object`成员变量中，我们可以从那里多次读取它，如果它不能被复制，也可以从那里移开。
     *
     */
    class TaskData
    {
    public:
      /**
       * 构造函数。初始化一个 std::future
       * 对象，并假定这样设置的任务还没有完成。
       *
       */
      TaskData(std::future<RT> &&future)
        : future(std::move(future))
        , task_has_finished(false)
      {}

      /**
       * 等待 std::future 对象准备好，即等待 std::promise
       * 接收其值的时间。如果这已经发生了，这个函数可以遵循一个快速路径。
       *
       */
      void
      wait()
      {
        // If we have previously already moved the result, then we don't
        // need a lock and can just return.
        if (task_has_finished)
          return;

        // Else, we need to go under a lock and try again. A different thread
        // may have waited and finished the task since then, so we have to try
        // a second time. (This is Schmidt's double-checking pattern.)
        std::lock_guard<std::mutex> lock(mutex);
        if (task_has_finished)
          return;
        else
          {
            // Wait for the task to finish and then move its
            // result. (We could have made the set_from() function
            // that we call here wait for the future to be ready --
            // which happens implicitly when it calls future.get() --
            // but that would have required putting an explicit
            // future.wait() into the implementation of
            // internal::return_value<void>::set_from(), which is a
            // bit awkward: that class doesn't actually need to set
            // anything, and so it looks odd to have the explicit call
            // to future.wait() in the set_from() function. Avoid the
            // issue by just explicitly calling future.wait() here.)
            future.wait();
            returned_object.set_from(future);

            // Now we can safely set the flag and return.
            task_has_finished = true;
          }
      }



      typename internal::return_value<RT>::reference_type
      get()
      {
        wait();
        return returned_object.get();
      }

    private:
      /**
       * 一个用于同步访问该类数据结构的mutex。
       *
       */
      std::mutex mutex;

      /**
       * 与当前类所代表的任务相关的承诺。
       *
       */
      std::future<RT> future;

      /**
       * 一个布尔值，表示有关任务是否已经完成。
       * @note 我们在这里使用一个 `std::atomic_bool`
       * ，因为我们必须确保线程之间的并发读取和存储是正确同步的，并且在一个特定线程上的顺序读取不会被重新排序或优化掉。一个
       * std::atomic
       * [1]实现了这一点，因为（如果没有其他注释的话）对布尔的读取和存储都受制于
       * std::memory_order_seq_cst
       * 内存排序[2]。这确保了Schmidt的双重检查确实有效。更多信息（以及可能更有效的实现）请参见[3]。
       * [1] https://en.cppreference.com/w/cpp/atomic/atomic [2]
       * https://en.cppreference.com/w/cpp/atomic/memory_order [3]
       * https://preshing.com/20130930/double-checked-locking-is-fixed-in-cpp11/
       *
       */
      std::atomic<bool> task_has_finished;

      /**
       * 一旦 std::future 交付后，返回值将被移到的地方。
       *
       */
      internal::return_value<RT> returned_object;
    };

    /**
     * 一个指向描述任务及其返回值的对象描述符的指针。
     *
     */
    std::shared_ptr<TaskData> task_data;
  };



  /**
   * new_task函数的重载，用于可以转换为 std::function<RT
   * ()>的对象，即任何可以像函数对象一样被调用而不需要参数并返回RT（或void）类型的对象。
   * @note  当 MultithreadInfo::n_threads()
   * 返回1时，即如果deal.II运行时系统被配置为只使用一个线程，那么这个函数只是立即执行给定的函数对象，并将返回值存储在由该函数返回的任务对象中。
   * @note   Threads::new_task() 本质上等同于调用
   * `std::async(std::launch::async,
   * ...)`，因为它在后台运行给定的任务。更多信息见https://en.cppreference.com/w/cpp/thread/async。
   * @ingroup threads
   *
   */
  template <typename RT>
  inline Task<RT>
  new_task(const std::function<RT()> &function)
  {
    return Task<RT>(function);
  }



  /**
   * new_task函数的重载，对象可以像函数对象一样被调用，没有参数。特别是，这个函数允许用使用
   * std::bind, 或使用lambda函数产生的对象来调用
   * Threads::new_task()
   * 。例如，在编写如下代码时，可以调用这个函数
   * @code
   * Threads::Task<int>
   * task = Threads::new_task ( [] () {
   *                                    do_this();
   *                                    then_do_that();
   *                                    return 42;
   *                                  });
   * @endcode
   * 在这里，我们把对函数序列 <code>do_this()</code> and
   * <code>then_do_that()</code>
   * 的调用安排在一个单独的任务上，把这里声明的lambda函数作为任务上执行的函数。然后lambda函数返回42（这在这里有点无意义，但它当然可以是一些计算出来的数字），这将是你以后在任务（即lambda函数的主体）完成后可以通过
   * <code>task.return_value()</code> 检索到的返回值。
   * @note  当 MultithreadInfo::n_threads()
   * 返回1时，即如果deal.II运行时系统被配置为只使用一个线程，那么这个函数只是立即执行给定的函数对象，并将返回值存储在由该函数返回的任务对象中。
   * @note
   * 每个lambda函数（或者你在这里传递给new_task()函数的其他东西，例如一个
   * std::bind()
   * 表达式的结果）都有一个返回类型，并因此返回一个该类型的对象。这个类型可以通过在这个函数的声明中使用的C++11
   * <code>decltype</code>
   * 语句来推断，然后它被用作当前函数返回的 Threads::Task
   * 对象的模板参数。
   * 在上面的例子中，因为lambda函数返回42（在C++中它的数据类型是
   * <code>int</code> ），推断的类型是 <code>int</code>
   * ，任务对象的类型将是 <code>Task@<int@></code>
   * 。换句话说，在用户代码中不能<i>necessary</i>明确指定lambda或
   * std::bind
   * 表达式的返回类型，尽管可以通过（完全等价）写出明确的方式来做到这一点
   * @code
   * Threads::Task<int>
   *   task = Threads::new_task ( [] ()
   *
   * -> int {
   *                                             do_this();
   *                                             then_do_that();
   *                                             return 42;
   *                                           });
   * @endcode
   *
   * @note
   * 在实践中，你将传递给new_task()的lambda函数当然会更加复杂。
   * 特别是，它们可能会从周围的环境中<i>capture</i>变量并在lambda中使用它们。
   * 关于lambda函数如何工作的更多信息，请参阅https://en.wikipedia.org/wiki/Anonymous_function#C.2B.2B_.28since_C.2B.2B11.29。
   * @note
   * 如果你把一个lambda函数作为参数传递给当前函数，该函数捕获了一个变量<i>by
   * reference</i>，或者如果你使用一个 std::bind ，用 std::ref()
   * 或 std::cref(),
   * 将一个函数参数绑定到一个引用变量，那么显然你只能在你引用或捕获的变量有一个至少延伸到任务结束的时间的寿命时才能这样做。
   * @note   Threads::new_task() 本质上等同于调用
   * `std::async(std::launch::async,
   * ...)`，因为它在后台运行指定的任务。更多信息见https://en.cppreference.com/w/cpp/thread/async。
   * @ingroup CPP11
   *
   */
  template <typename FunctionObjectType>
  inline auto
  new_task(FunctionObjectType function_object)
    -> Task<decltype(function_object())>
  {
    using return_type = decltype(function_object());
    dealii::MultithreadInfo::initialize_multithreading();
    return new_task(std::function<return_type()>(function_object));
  }



  /**
   * new_task函数的重载，用于非成员或静态成员函数。更多信息见同名的其他函数。
   * @ingroup threads
   *
   */
  template <typename RT, typename... Args>
  inline Task<RT>
  new_task(RT (*fun_ptr)(Args...), typename identity<Args>::type... args)
  {
    auto dummy = std::make_tuple(internal::maybe_make_ref<Args>::act(args)...);
    return new_task(
      [dummy, fun_ptr]() -> RT { return std_cxx17::apply(fun_ptr, dummy); });
  }



  /**
   * 非const
   * new_task函数的重载。更多信息请参见同名的其他函数。
   * @ingroup threads
   *
   */
  template <typename RT, typename C, typename... Args>
  inline Task<RT>
  new_task(RT (C::*fun_ptr)(Args...),
           typename identity<C>::type &c,
           typename identity<Args>::type... args)
  {
    // NOLINTNEXTLINE(modernize-avoid-bind) silence clang-tidy
    return new_task(std::function<RT()>(std::bind(
      fun_ptr, std::ref(c), internal::maybe_make_ref<Args>::act(args)...)));
  }

  /**
   * new_task函数的重载。更多信息请参见同名的其他函数。
   * @ingroup threads
   *
   */
  template <typename RT, typename C, typename... Args>
  inline Task<RT>
  new_task(RT (C::*fun_ptr)(Args...) const,
           typename identity<const C>::type &c,
           typename identity<Args>::type... args)
  {
    // NOLINTNEXTLINE(modernize-avoid-bind) silence clang-tidy
    return new_task(std::function<RT()>(std::bind(
      fun_ptr, std::cref(c), internal::maybe_make_ref<Args>::act(args)...)));
  }


  // ------------------------ TaskGroup -------------------------------------

  /**
   * 一个任务对象的容器。允许添加新的任务对象并一起等待它们。任务对象需要对被调用的函数有相同的返回值。
   * 注意，对join_all()的调用必须与添加子任务的调用在同一个线程上执行。否则，可能会出现死锁。换句话说，一个任务对象绝不应该传递给另一个任务来调用join()方法。
   * @ingroup tasks
   *
   */
  template <typename RT = void>
  class TaskGroup
  {
  public:
    /**
     * 将另一个任务对象添加到集合中。
     *
     */
    TaskGroup &
    operator+=(const Task<RT> &t)
    {
      tasks.push_back(t);
      return *this;
    }


    /**
     * 返回有多少个任务被放入这个组。这个函数不区分其中有多少任务已经运行并完成，仍在等待被安排到CPU资源，或目前正在运行。已经加入的任务也仍然被计算在内。
     *
     */
    std::size_t
    size() const
    {
      return tasks.size();
    }


    /**
     * 等待集合中的所有任务完成。如果其中有些任务已经被等待过了，也不是什么问题，也就是说，你可以多次调用这个函数，如果你愿意，你也可以在后续调用这个函数之间添加新任务对象。
     *
     */
    void
    join_all() const
    {
      for (auto &t : tasks)
        t.join();
    }

  private:
    /**
     * 任务对象的列表。
     *
     */
    std::list<Task<RT>> tasks;
  };

} // namespace Threads

/**
 * @}
 *
 *
 */


//---------------------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
// end of #ifndef dealii_thread_management_h
#endif
//---------------------------------------------------------------------------


