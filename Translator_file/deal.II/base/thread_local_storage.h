//include/deal.II-translator/base/thread_local_storage_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2021 by the deal.II authors
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

#ifndef dealii_thread_local_storage_h
#  define dealii_thread_local_storage_h


#  include <deal.II/base/config.h>

#  include <deal.II/base/exceptions.h>

#  include <list>
#  include <map>
#  include <memory>
#  include <shared_mutex>
#  include <thread>
#  include <vector>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup threads */ 
 /*@{*/ 

#  ifndef DOXYGEN
class LogStream;
#  endif

namespace Threads
{
#  ifndef DOXYGEN
  namespace internal
  {
    /* 解决方法。当涉及到STL容器和包含不可复制的对象T时，标准不幸地在 std::is_copy_constructible 的类型特征中存在一个不幸的设计 "缺陷"。通过解压一些常用的容器来解决这个问题。   
*
*/
    template <typename T>
    struct unpack_container
    {
      using type = T;
    };

    template <typename T, typename A>
    struct unpack_container<std::vector<T, A>>
    {
      using type = T;
    };

    template <typename T, typename A>
    struct unpack_container<std::list<T, A>>
    {
      using type = T;
    };
  } // namespace internal
#  endif

  /**
   * @brief A class that provides a separate storage location on each thread
   * ，访问该对象。
   * 这个类提供的方式是，每个访问它的线程都有自己的T类型对象的副本。从本质上讲，在多线程程序中，访问这个对象永远不会导致竞赛条件，因为除了当前的线程，没有其他线程可以访问它。
   * <h3>Construction and destruction</h3>
   * 这个类的对象可以是默认构建的，也可以通过提供一个
   * "典范"，即一个T类型的对象，这样，每次我们需要在一个还没有这样一个对象的线程上创建一个T，它就会从典范中复制出来。
   * 在销毁这个类的对象时，所有对应于访问过这个对象的线程的T对象都会被销毁。请注意，这可能是在线程终止的时间之前。
   * <h3>Access</h3>
   * 这个对象所存储的T对象可以使用get()函数来访问。当从不同的线程访问时，它提供了一个对唯一对象的引用。T类型的对象是懒惰地创建的，也就是说，只有当线程实际调用get()时才会创建它们。
   *
   */
  template <typename T>
  class ThreadLocalStorage
  {
    static_assert(
      std::is_copy_constructible<
        typename internal::unpack_container<T>::type>::value ||
        std::is_default_constructible<T>::value,
      "The stored type must be either copyable, or default constructible");

  public:
    /**
     * 默认构造函数。使用默认构造函数初始化每个线程局部对象。
     *
     */
    ThreadLocalStorage() = default;

    /**
     * 复制构造函数。
     *
     */
    ThreadLocalStorage(const ThreadLocalStorage &);

    /**
     * 移动构造函数。构造函数从参数中移动所有的内部数据结构。
     *
     */
    ThreadLocalStorage(ThreadLocalStorage &&t) noexcept;

    /**
     * 一种复制构造函数。通过给定的对象初始化一个内部典范。该范例反过来被用来初始化每个线程的本地对象，而不是调用默认的构造函数。
     *
     */
    explicit ThreadLocalStorage(const T &t);

    /**
     * 一种移动构造函数。将给定的对象移动到一个内部样板中。该范例反过来被用来初始化每个线程的本地对象，而不是调用默认的构造函数。
     *
     */
    explicit ThreadLocalStorage(T &&t);

    /**
     * 复制赋值运算符。
     *
     */
    ThreadLocalStorage &
    operator=(const ThreadLocalStorage &t);

    /**
     * 移动赋值运算符。
     *
     */
    ThreadLocalStorage &
    operator=(ThreadLocalStorage &&t) noexcept;

    /**
     * 返回这个对象在当前线程中存储的数据的引用，这个函数被调用。
     * 请注意，没有一个成员函数get()是常数，并且像人们期望的那样返回一个常数引用。原因是，如果在一个还没有创建线程本地对象的线程上调用这样一个成员函数，那么就必须先创建这样一个对象，这当然是一个非常量操作。如果你需要从一个常量成员函数中调用一个类的成员变量的get()函数，那么你需要声明成员变量
     * <code>mutable</code> 以允许这种访问。
     *
     */
    T &
    get();

    /**
     * 与上述相同，除了 @p exists
     * 在当前线程已经存在一个元素的情况下被设置为真，否则为假。
     *
     */
    T &
    get(bool &exists);

    /**
     * 转换操作符，简单地将线程-本地对象转换为它所存储的数据类型。这个函数等同于调用get()成员函数；它的目的是使TLS对象看起来更像它所存储的对象。
     *
     */
    operator T &();

    /**
     * 将给定的参数复制到用于表示当前线程的存储空间中。以
     * <code>tls_data = object</code>
     * 的形式调用这个函数，相当于调用 <code>tls_data.get() =
     * object</code>
     * 。这个操作符的意图是使ThreadLocalStorage对象看起来更像它在当前线程上代表的对象。
     * @param  t 要复制到当前线程使用的存储空间的对象。
     * @return  当前的对象，在做出改变之后
     *
     */
    ThreadLocalStorage<T> &
    operator=(const T &t);

    /**
     * 将给定的参数移入用于表示当前线程的存储空间。以<code>tls_data
     * =
     * object</code>的方式调用此函数，相当于调用<code>tls_data.get()
     * =
     * object</code>。这个操作符的意图是使ThreadLocalStorage对象看起来更像它在当前线程中代表的对象。移动赋值运算符。
     * @param  t 要复制到当前线程所用存储空间的对象。
     * @return  当前的对象，在做出改变之后
     *
     */
    ThreadLocalStorage<T> &
    operator=(T &&t);

    /**
     * 移除为所有用此对象创建过的线程所存储的线程本地对象（即，在此线程上至少调用过一次get()。这包括当前的线程。如果你随后在这个线程或任何其他线程上调用get()，新的对象将再次被创建。
     * 如果deal.II被配置为不使用多线程，那么这个函数根本就不做任何事情。请注意，这当然有不同的语义，因为在多线程背景下，对象被删除并再次创建（如果调用了该类的适当构造函数，可以通过从样本对象中复制），而在多线程背景下，对象根本就没有被触及。同时，这个函数的目的是释放其他线程可能已经为他们自己的线程局部对象分配的内存，在这之后，每次使用这个对象都需要某种初始化。在多线程和非多线程的情况下，这都是必要的。
     *
     */
    void
    clear();

  private:
    /**
     * 我们存储的数据元素。
     *
     */
    std::map<std::thread::id, T> data;

    /**
     * 一个突变器，用于保护插入数据对象。
     * 我们在这里使用一个 std::shared_timed_mutex （或
     * std::shared_mutex ，如果有的话），以便能够使用
     * std::unique_lock 和 std::shared_lock
     * 来实现读者-作者锁（https://en.wikipedia.org/wiki/Readers%E2%80%93writer_lock）。
     *
     */
#  ifdef DEAL_II_HAVE_CXX17
    mutable std::shared_mutex insertion_mutex;
#  else
    mutable std::shared_timed_mutex insertion_mutex;
#  endif

    /**
     * 创建一个新的（针对线程的）副本的范例。
     *
     */
    std::shared_ptr<const T> exemplar;

    friend class dealii::LogStream;
  };
} // namespace Threads
/**
 * @}
 *
 */

#  ifndef DOXYGEN
namespace Threads
{
  // ----------------- inline and template functions --------------------------


  template <typename T>
  ThreadLocalStorage<T>::ThreadLocalStorage(const ThreadLocalStorage<T> &t)
    : exemplar(t.exemplar)
  {
    // Raise a reader lock while we are populating our own data in order to
    // avoid copying over an invalid state.
    std::shared_lock<decltype(insertion_mutex)> lock(t.insertion_mutex);
    data = t.data;
  }



  template <typename T>
  ThreadLocalStorage<T>::ThreadLocalStorage(ThreadLocalStorage<T> &&t) noexcept
    : exemplar(std::move(t.exemplar))
  {
    // We are nice and raise the writer lock before copying over internal
    // data structures from the argument.
    //
    // The point is a bit moot, though: Users of ThreadLocalStorage
    // typically obtain their thread's thread-local object through the
    // get() function. That function also acquires the lock, but
    // whether or not we do that here really doesn't make any
    // difference in terms of correctness: If another thread manages
    // to call get() just before we get here, then the result of that
    // get() function immediately becomes invalid; if it manages to
    // call get() at the same time as this function if there were no
    // locking here, it might access undefined state; and if it
    // manages to call get() just after we moved away the state --
    // well, then it just got lucky to escape the race condition, but
    // the race condition is still there.
    //
    // On the other hand, there is no harm in doing at least
    // conceptually the right thing, so ask for that lock:
    std::unique_lock<decltype(insertion_mutex)> lock(t.insertion_mutex);
    data = std::move(t.data);
  }



  template <typename T>
  inline ThreadLocalStorage<T>::ThreadLocalStorage(const T &t)
    : exemplar(std::make_shared<const T>(t))
  {}



  template <typename T>
  inline ThreadLocalStorage<T>::ThreadLocalStorage(T &&t)
    : exemplar(std::make_shared<T>(std::forward<T>(t)))
  {}



  template <typename T>
  inline ThreadLocalStorage<T> &
  ThreadLocalStorage<T>::operator=(const ThreadLocalStorage<T> &t)
  {
    // We need to raise the reader lock of the argument and our writer lock
    // while copying internal data structures.
    std::shared_lock<decltype(insertion_mutex)> reader_lock(t.insertion_mutex);
    std::unique_lock<decltype(insertion_mutex)> writer_lock(insertion_mutex);

    data     = t.data;
    exemplar = t.exemplar;

    return *this;
  }



  template <typename T>
  inline ThreadLocalStorage<T> &
  ThreadLocalStorage<T>::operator=(ThreadLocalStorage<T> &&t) noexcept
  {
    // We need to raise the writer lock of the argument (because we're
    // moving information *away* from that object) and the writer lock
    // of our object while copying internal data structures.
    //
    // That said, the same issue with acquiring the source lock as
    // with the move constructor above applies here as well.
    std::unique_lock<decltype(insertion_mutex)> reader_lock(t.insertion_mutex);
    std::unique_lock<decltype(insertion_mutex)> writer_lock(insertion_mutex);

    data     = std::move(t.data);
    exemplar = std::move(t.exemplar);

    return *this;
  }


#    ifndef DOXYGEN
  namespace internal
  {
    /* 我们必须确保，如果相应的元素不能被复制构建，就不要调用 "data.emplace(id,exemplar)"。我们使用一些SFINAE魔法来解决C++14没有 "if constexpr "的问题。   
*
*/
    template <typename T>
    typename std::enable_if_t<
      std::is_copy_constructible<typename unpack_container<T>::type>::value,
      T &>
    construct_element(std::map<std::thread::id, T> &  data,
                      const std::thread::id &         id,
                      const std::shared_ptr<const T> &exemplar)
    {
      if (exemplar)
        {
          const auto it = data.emplace(id, *exemplar).first;
          return it->second;
        }
      return data[id];
    }

    template <typename T>
    typename std::enable_if_t<
      !std::is_copy_constructible<typename unpack_container<T>::type>::value,
      T &>
    construct_element(std::map<std::thread::id, T> &data,
                      const std::thread::id &       id,
                      const std::shared_ptr<const T> &)
    {
      return data[id];
    }
  } // namespace internal
#    endif


  template <typename T>
  inline T &
  ThreadLocalStorage<T>::get(bool &exists)
  {
    const std::thread::id my_id = std::this_thread::get_id();

    // Note that std::map<..>::emplace guarantees that no iterators or
    // references to stored objects are invalidated. We thus only have to
    // ensure that we do not perform a lookup while writing, and that we
    // do not write concurrently. This is precisely the "reader-writer
    // lock" paradigm supported by C++14 by means of the std::shared_lock
    // and the std::unique_lock.

    {
      // Take a shared ("reader") lock for lookup and record the fact
      // whether we could find an entry in the boolean exists.
      std::shared_lock<decltype(insertion_mutex)> lock(insertion_mutex);

      const auto it = data.find(my_id);
      if (it != data.end())
        {
          exists = true;
          return it->second;
        }
      else
        {
          exists = false;
        }
    }

    {
      // Take a unique ("writer") lock for manipulating the std::map. This
      // lock ensures that no other threat does a lookup at the same time.
      std::unique_lock<decltype(insertion_mutex)> lock(insertion_mutex);

      return internal::construct_element(data, my_id, exemplar);
    }
  }


  template <typename T>
  inline T &
  ThreadLocalStorage<T>::get()
  {
    bool exists;
    return get(exists);
  }


  template <typename T>
  inline ThreadLocalStorage<T>::operator T &()
  {
    return get();
  }


  template <typename T>
  inline ThreadLocalStorage<T> &
  ThreadLocalStorage<T>::operator=(const T &t)
  {
    get() = t;
    return *this;
  }


  template <typename T>
  inline ThreadLocalStorage<T> &
  ThreadLocalStorage<T>::operator=(T &&t)
  {
    get() = std::forward<T>(t);
    return *this;
  }


  template <typename T>
  inline void
  ThreadLocalStorage<T>::clear()
  {
    std::unique_lock<decltype(insertion_mutex)> lock(insertion_mutex);
    data.clear();
  }
} // namespace Threads

#  endif // DOXYGEN

//---------------------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
// end of #ifndef dealii_thread_local_storage_h
#endif
//---------------------------------------------------------------------------


