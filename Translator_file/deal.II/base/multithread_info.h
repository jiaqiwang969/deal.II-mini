//include/deal.II-translator/base/multithread_info_0.txt
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

#ifndef dealii_multithread_info_h
#  define dealii_multithread_info_h
//---------------------------------------------------------------------------


#  include <deal.II/base/config.h>

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/types.h>

#  include <memory>

// forward declaration from <taskflow/taskflow.hpp>
namespace tf
{
  class Executor;
}

DEAL_II_NAMESPACE_OPEN

/**
 * 这个类提供了在多线程程序中可能有用的系统信息。 目前，这只是CPU的数量。如果deal.II在编译时支持多线程，一些函数将使用多个线程进行操作。目前该库同时支持基于线程和基于任务的并行性。
 * @ref threads
 * 描述了每种的不同用途。用于基于任务的并行方法的默认线程数由线程构件库自动选择。有关这方面的更多信息，请参见
 * @ref threads 。
 * 基于线程的并行方法需要明确地创建线程，并且可能希望使用与系统中CPU数量相关的线程数量。推荐的线程数可以用
 * MultithreadInfo::n_threads(), 来查询，而系统中的内核数则由
 * MultithreadInfo::n_cores(). 返回。
 *
 *
 * @ingroup threads
 *
 *
 */
class MultithreadInfo
{
public:
  /**
   * 构造函数。这个构造函数被删除，因为不需要构造这个类的实例（所有成员都是静态的）。
   *
   */
  MultithreadInfo() = delete;

  /**
   * 系统中CPU的数量。    这个内部调用
   * [<code>std::thread::hardware_concurrency</code>](https://en.cppreference.com/w/cpp/thread/thread/hardware_concurrency)
   * ，但如果调用返回错误，则将结果设置为1。
   *
   */
  static unsigned int
  n_cores();

  /**
   * 返回要使用的线程数。这最初被设置为系统的核心数（见n_cores()），但可以通过set_thread_limit()和环境变量DEAL_II_NUM_THREADS进一步限制。
   *
   */
  static unsigned int
  n_threads();

  /**
   * 返回这个对象的内存消耗估计值，单位是字节。
   * 这不是精确的（但通常会很接近），因为计算树的内存使用量（例如，
   * <tt>std::map</tt>) 是困难的。
   *
   */
  static std::size_t
  memory_consumption();

  /**
   * 将要使用的最大线程数设置为环境变量DEAL_II_NUM_THREADS和给定参数（或其默认值）的最小值。这将影响到TBB的初始化。如果两者都没有给出，则使用TBB的默认值（基于系统中的核心数量）。
   * 在你main()中的代码运行之前，这个例程会自动执行，并使用默认参数（使用静态构造函数）。它也被
   * Utilities::MPI::MPI_InitFinalize.
   * 执行。如果你有一个基于MPI的代码，使用
   * Utilities::MPI::MPI_InitFinalize 的构造函数的适当参数。
   *
   */
  static void
  set_thread_limit(
    const unsigned int max_threads = numbers::invalid_unsigned_int);

  /**
   * 如果TBB使用单线程运行，则返回，要么是因为线程亲和性，要么是因为通过调用set_thread_limit来设置。这在PETScWrappers中使用，以避免使用非线程安全的接口。
   *
   */
  static bool
  is_running_single_threaded();

  /**
   * 确保多线程API被初始化。这通常不需要在usercode中调用。
   *
   */
  static void
  initialize_multithreading();


#  ifdef DEAL_II_WITH_TASKFLOW
  /**
   * 从taskflow返回对全局Executor的引用。
   * 该执行器被设置为使用n_threads()工作线程，你可以使用set_thread_limit()和DEAL_II_NUM_THREADS环境变量控制。
   *
   */
  static tf::Executor &
  get_taskflow_executor();
#  endif

private:
  /**
   * 代表最大线程数的变量。
   *
   */
  static unsigned int n_max_threads;

#  ifdef DEAL_II_WITH_TASKFLOW
  /**
   * 存储一个用N个工作者构建的任务流执行器（来自set_thread_limit）。
   *
   */
  static std::unique_ptr<tf::Executor> executor;
#  endif
};



//---------------------------------------------------------------------------
DEAL_II_NAMESPACE_CLOSE
// end of #ifndef dealii_multithread_info_h
#endif
//---------------------------------------------------------------------------


