//include/deal.II-translator/base/logstream_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_logstream_h
#define dealii_logstream_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/thread_local_storage.h>

#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <stack>
#include <string>


DEAL_II_NAMESPACE_OPEN

/**
 * 一个简化了执行日志过程的类。它通过提供  <ul>   <li>  前缀的推送和弹出机制，以及  <li>  将信息分发到文件和控制台的可能性来实现。  </ul>
 * 这个类的通常用法是通过预生成的对象<tt>deallog</tt>。典型的设置步骤是。  <ul>   <li>  <tt>deallog.depth_console(n)</tt>: 将屏幕上的输出限制在外循环。  <li>   <tt>deallog.attach(std::ostream)</tt>:  将日志信息写进文件。  <li>  <tt>deallog.depth_file(n)</tt>: 限制输出到文件的外循环。  </ul>
 * 在进入程序的新阶段之前，例如一个新的循环，可以通过
 * <tt>LogStream::Prefix
 * p("loopname");</tt>设置一个新的前缀。前缀的析构器将从堆栈中弹出前缀文本。
 * 通过<tt>&lt;&lt;</tt>操作符写入，<tt> deallog << "This is a log
 * notice";</tt>将被本地缓冲线程，直到遇到 <tt>std::flush</tt>
 * 或 <tt>std::endl</tt>
 * ，这将触发写出到控制台，如果设置了日志文件。
 * <h3>LogStream and thread safety</h3>
 * 在并发线程的附近，LogStream的行为方式如下。  <ul>   <li>  每一次对Logstream的写操作<tt>&lt;&lt;</tt>（或使用其中一个特殊的成员函数）都会在线程本地存储中进行缓冲。  <li>  一个 <tt>std::flush</tt> 或 <tt>std::endl</tt> 将触发一个写出到控制台和（如果附加）到文件流。这个写出是有顺序的，所以并发的线程的输出不会交错。  <li>  在一个新的线程上，调用写出，以及调用#push或#pop将把创建LogStream实例的 "祝福 "线程的当前前缀复制到线程本地存储。此后，前缀是线程本地的。  </ul>
 *
 *
 * @ingroup textoutput
 *
 */
class LogStream : public Subscriptor
{
public:
  /**
   * 一个允许安全生成和移除前缀的子类。
   * 在一个区块的某个地方，创建一个这样的对象，它将作为前缀出现在LogStream输出中，就像
   * @p deallog.
   * 在区块结束时，前缀将自动被删除，当这个对象被销毁时。
   * 换句话说，这样创建的对象的范围决定了前缀的寿命。使用这样一个对象的好处是，无论你以何种方式退出这个范围，前缀都会被删除
   *
   * --通过 <code>continue</code>, <code>break</code>, <code>return</code>
   * 、 <code>throw</code>
   * ，或者仅仅通过到达关闭括号。在所有这些情况下，没有必要记得使用
   * LogStream::pop(). 来手动弹出前缀
   * 在这一点上，它的工作方式就像更著名的 std::unique_ptr
   * 和 std::lock_guard 类。
   *
   */
  class Prefix
  {
  public:
    /**
     * 为 @p deallog,
     * 设置一个新的前缀，这个前缀将在变量被销毁时被删除。
     *
     */
    Prefix(const std::string &text);

    /**
     * 为给定的流设置一个新的前缀，这个前缀将在变量被销毁时被删除。
     *
     */
    Prefix(const std::string &text, LogStream &stream);

    /**
     * 删除与此变量相关的前缀。
     *
     */
    ~Prefix();

  private:
    /**
     * 一个指向应用前缀的LogStream对象的指针。
     *
     */
    SmartPointer<LogStream, LogStream::Prefix> stream;
  };


  /**
   * 标准构造函数。构造函数将输出流设置为 <tt>std::cout</tt>
   * ，深度设置为零。(使用attach()和depth_console()来改变这个。)
   *
   */
  LogStream();


  /**
   * 解构器。
   *
   */
  ~LogStream() override;


  /**
   * 启用输出到第二个流<tt>o</tt>。      @param  o
   * 附加这个输出流。      @param[in]  print_job_id
   * 当前进程的JobIdentifier是否应该被打印到流中。
   * @param[in]  flags 要在输出流上设置的格式标志  @p o.
   *
   */
  void
  attach(std::ostream &                o,
         const bool                    print_job_id = true,
         const std::ios_base::fmtflags flags        = std::ios::showpoint |
                                               std::ios::left);


  /**
   * 禁用对第二个流的输出。你可能想在以前连接到这个对象的流上调用<tt>close</tt>。
   *
   */
  void
  detach();


  /**
   * 返回默认流（<tt>std_out</tt>）。
   *
   */
  std::ostream &
  get_console();


  /**
   * 返回文件流。
   *
   */
  std::ostream &
  get_file_stream();


  /**
   * 如果文件流已经被连接，返回 @p true ，否则返回 @p false
   * 。
   *
   */
  bool
  has_file() const;


  /**
   * 返回前缀字符串。
   *
   */
  const std::string &
  get_prefix() const;


  /**
   * 在堆栈中推送另一个前缀。前缀自动用冒号隔开，最后一个前缀后面有一个双冒号。
   * 一个更简单的添加前缀的方法（无需手动添加相应的pop()）是使用
   * LogStream::Prefix
   * 类。使用该类的好处是，只要Prefix对象超出了范围，就会发出相应的pop()调用
   *
   * - 无论是在代码块的末尾，还是在最近的 @p return 语句，或者是因为中间函数调用导致的异常没有被立即捕获。
   *
   */
  void
  push(const std::string &text);


  /**
   * 删除用push()添加的最后一个前缀。
   *
   */
  void
  pop();


  /**
   * 在控制台打印的最大级别数。默认为0，不会产生任何输出。这个函数允许人们将控制台输出限制在最高级别的迭代。只有小于<tt>n</tt>前缀的输出被打印。在<tt>n=0</tt>的情况下调用这个函数，将不会写出控制台输出。参见
   * step-3 中关于此方法的使用实例。
   * 该参数的前一个值将被返回。
   *
   */
  unsigned int
  depth_console(const unsigned int n);


  /**
   * 写入日志文件的最大层数。其功能与<tt>depth_console</tt>相同，尽管如此，这个函数应该谨慎使用，因为它可能破坏日志文件的价值。
   * 该参数的前一个值将被返回。
   *
   */
  unsigned int
  depth_file(const unsigned int n);


  /**
   * 记录线程的ID。
   *
   */
  bool
  log_thread_id(const bool flag);


  /**
   * 为底层流设置精度，并返回之前的流精度。这个函数模仿http://www.cplusplus.com/reference/ios/ios_base/precision/
   *
   */
  std::streamsize
  precision(const std::streamsize prec);


  /**
   * 设置底层流的宽度，并返回上一个流的宽度。此函数模仿http://www.cplusplus.com/reference/ios/ios_base/width/
   *
   */
  std::streamsize
  width(const std::streamsize wide);


  /**
   * 设置底层流的标志，并返回之前的流标志。此函数模仿http://www.cplusplus.com/reference/ios/ios_base/flags/
   *
   */
  std::ios::fmtflags
  flags(const std::ios::fmtflags f);


  /**
   * 处理ostream操作器。这将整个事情传递给模板函数，但
   * <tt>std::endl</tt>
   * 操纵器除外，对其进行特殊操作：将临时流缓冲区包括头写入文件和
   * <tt>std::cout</tt> 并清空缓冲区。
   * 反正这个函数的重载是需要的，因为编译器不能像以前的通用模板那样，将
   * @p std::endl 这样的操纵器直接与模板参数 @p T
   * 绑定。这是由于 @p  std::endl 实际上是 @p std::ostream,  @p
   * std::wostream,
   * 以及潜在的更多此类函数的重载集。因此，这个函数有必要从这个重载集合中挑选一个元素。
   *
   */
  LogStream &
  operator<<(std::ostream &(*p)(std::ostream &));


  /**
   * 返回这个对象的内存消耗估计值，单位是字节。
   * 这不是精确的（但通常会很接近），因为计算树（例如，
   * <tt>std::map</tt>) ）的内存用量是很困难的。
   *
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * 围绕线程本地前缀的内部包装器。这个私有函数将返回正确的内部前缀栈。更重要的是，一个新的线程本地堆栈将从创建这个LogStream实例的
   * "受祝福 "线程（通常，在deallog的情况下，是 "主
   * "线程）的当前堆栈中复制出来。
   *
   */
  std::stack<std::string> &
  get_prefixes() const;

  /**
   * 堆栈中的字符串，这些字符串被打印在每一行的开头，以便于识别输出是在哪里产生的。
   *
   */
  mutable Threads::ThreadLocalStorage<std::stack<std::string>> prefixes;

  /**
   * 我们记录创建此对象的线程ID。我们需要这个信息，以便在一个新的线程上第一次使用deallog时，从这个
   * "父 "线程 "偷 "出当前的前缀。
   *
   */
  std::thread::id parent_thread;

  /**
   * 默认的流，即输出要去的地方。这个流默认为
   * <tt>std::cout</tt>, ，但可以通过构造函数设置为其他流。
   *
   */
  std::ostream *std_out;

  /**
   * 指向一个流的指针，输出的副本将被送到那里。通常，这将是一个文件流。
   * 你可以通过<tt>attach</tt>函数设置和重置这个流。
   *
   */
  std::ostream *file;

  /**
   * 表示要打印到标准输出的前缀数量的值。如果超过这个数量的前缀被推到堆栈，那么将不会产生输出，直到前缀的数量缩减到这个数字以下。
   *
   */
  unsigned int std_depth;

  /**
   * 输出到文件的最大前缀深度也一样。
   *
   */
  unsigned int file_depth;

  /**
   * 打印线程ID的标志。
   *
   */
  bool print_thread_id;

  /**
   * 指示输出当前是否在新行的标志
   *
   */
  bool at_newline;

  /**
   * 打印行的头部。
   *
   */
  void
  print_line_head();

  /**
   * 围绕 "线程本地
   * "输出流的内部包装器。这个私有函数将为operator<<返回正确的内部ostringstream缓冲区。
   *
   */
  std::ostringstream &
  get_stream();

  /**
   * 我们使用我们的线程本地存储设施，为每个发送日志信息的线程生成一个字符串流。
   *
   */
  Threads::ThreadLocalStorage<std::shared_ptr<std::ostringstream>> outstreams;

  template <typename T>
  friend LogStream &
  operator<<(LogStream &log, const T &t);
};


 /* ----------------------------- Inline functions and templates ----------------
 */ 


/**
 * 通过LogStream输出一个不变的东西。
 *
 *
 * @note
 * 我们将这个操作符声明为一个非成员函数，这样就有可能在C++11重载解析规则下用更专业的模板化版本来重载它
 *
 *
 */
template <typename T>
inline LogStream &
operator<<(LogStream &log, const T &t)
{
  // print to the internal stringstream
  log.get_stream() << t;
  return log;
}



/**
 * deal.II的标准日志对象。
 *
 *
 */
extern LogStream deallog;



DEAL_II_NAMESPACE_CLOSE

#endif


