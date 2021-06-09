//include/deal.II-translator/base/timer_0.txt
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

#ifndef dealii_timer_h
#define dealii_timer_h

#include <deal.II/base/config.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>

#include <chrono>
#include <list>
#include <map>
#include <string>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个时钟，与 <code>std::chrono</code>
 * 中的时钟概念兼容，其now()方法返回一个时间点，表示当前进程使用的CPU时间量。
 *
 *
 */
struct CPUClock
{
  /**
   * 持续时间类型。Windows默认以1/64秒的倍数来测量CPU时间，而POSIX使用微秒，所以为了统一，使用微秒。
   *
   */
  using duration = std::chrono::microseconds;

  /**
   * 有符号的积分类型，用于存储count()返回的值。
   *
   */
  using rep = duration::rep;

  /**
   * 代表一个周期长度的比率（以秒为单位）。
   *
   */
  using period = duration::period;

  /**
   * 时间点类型。
   *
   */
  using time_point = std::chrono::time_point<CPUClock, duration>;

  /**
   * 表示时钟单调增长的布尔值。
   *
   */
  static const bool is_steady = true;

  /**
   * 返回当前进程所使用的CPU时间的数量。不幸的是，这需要特定平台的调用，所以这个函数在既不是Windows也不是POSIX的平台上返回0。
   *
   */
  static time_point
  now() noexcept;
};

/**
 * 计时器类提供了一种方法来测量挂钟的时间量（即挂钟上经过的时间量）和一个应用程序的某些部分所使用的CPU时间量。该类还提供了在MPI通信器上同步经过的时间的设施。
 * <h3>Usage</h3>
 * 计时器类可以被多次启动和停止。它既存储了最后一次启动-停止周期的时间量，或
 * <em>  圈  </em>
 * ，也存储了所有圈的总时间。下面是一个例子。
 *
 *
 * @code
 * Timer timer; // creating a timer also starts it
 *
 * // do some complicated computations here
 * // ...
 *
 * timer.stop();
 *
 * std::cout << "Elapsed CPU time: " << timer.cpu_time() << " seconds.\n";
 * std::cout << "Elapsed wall time: " << timer.wall_time() << " seconds.\n";
 *
 * // reset timer for the next thing it shall do
 * timer.reset();
 * @endcode
 * 另外，你也可以重新启动定时器，而不是重设它。然后，连续调用start()和stop()之间的时间（也就是圈数）将被累积。这个类的用法在
 * step-28 的教程程序中也有解释。
 *
 *
 * @note  TimerOutput（结合 TimerOutput::Scope)
 * 类提供了一个方便的方法来为多个命名的部分计时并总结输出。
 *
 *
 * @note
 * 这个类的实现取决于系统。特别是，CPU时间是由所有线程的总和累积而成的，通常会超过墙的时间。
 *
 *
 * @ingroup utilities
 *
 *
 */
class Timer
{
public:
  /**
   * 构造函数。将累计次数设置为零，并调用 Timer::start().
   * 。
   *
   */
  Timer();

  /**
   * 构造函数，指定CPU时间应在给定的通信器上进行求和。如果
   * @p sync_lap_times 是 <code>true</code>
   * ，那么定时器将把上一圈的耗时墙和CPU时间设置为它们在所提供的通信器上的最大值。只有当
   * Timer::stop()
   * 在定时器被查询到时间长度值之前被调用时，才会进行这种同步。
   * 该构造函数调用 Timer::start(). 。
   * @note
   * 在通过通信器发生同步之前，定时器被停止；同步的额外费用不被衡量。
   *
   */
  Timer(const MPI_Comm &mpi_communicator, const bool sync_lap_times = false);

  /**
   * 返回一个对数据结构的引用，该结构包含在给定通信器中所有MPI进程中测得的最后一圈的墙体时间的基本统计数据。这个结构在调用
   * Timer::stop() 之前不包含有意义的值。
   *
   */
  const Utilities::MPI::MinMaxAvg &
  get_last_lap_wall_time_data() const;

  /**
   * 返回一个数据结构的引用，该结构包含在给定通信器中所有MPI进程中测得的累积壁挂时间的基本统计数据。在调用
   * Timer::stop() 之前，这个结构不包含有意义的值。
   *
   */
  const Utilities::MPI::MinMaxAvg &
  get_accumulated_wall_time_data() const;

  /**
   * 打印由 Timer::get_last_lap_wall_time_data()
   * 返回的数据到给定的流。
   *
   */
  template <class StreamType>
  void
  print_last_lap_wall_time_data(StreamType &stream) const;

  /**
   * 打印由 Timer::get_accumulated_wall_time_data()
   * 返回的数据到给定的流中。
   *
   */
  template <class StreamType>
  void
  print_accumulated_wall_time_data(StreamType &stream) const;

  /**
   * 开始测量新的一圈。如果 <code>sync_lap_times</code> 是
   * <code>true</code>
   * ，那么就会使用MPI屏障，以确保所有进程在同一墙面时间开始测圈。
   *
   */
  void
  start();

  /**
   * 停止定时器。这将更新单圈时间和累积时间。如果
   * <code>sync_lap_times</code> is <code>true</code>
   * ，那么在通信器中的所有处理器上，单圈时间是同步的（即单圈时间被设置为最大单圈时间）。
   * 返回累积的CPU时间，单位为秒。
   *
   */
  double
  stop();

  /**
   * 停止计时器，如果它正在运行，并将所有测量值重置为默认状态。
   *
   */
  void
  reset();

  /**
   * 相当于调用 Timer::reset() ，然后再调用 Timer::start(). 。
   *
   */
  void
  restart();

  /**
   * 在不停止定时器的情况下，返回当前累积的壁挂时间（包括当前的一圈，如果定时器正在运行），单位为秒。
   *
   */
  double
  wall_time() const;

  /**
   * 返回最后一圈的壁挂时间，单位是秒。计时器不会因为这个函数而停止。
   *
   */
  double
  last_wall_time() const;

  /**
   * 在不停止计时器的情况下，返回累计的CPU时间（包括当前的一圈，如果计时器正在运行），单位为秒。
   * 如果向构造函数提供了一个MPI通信器，那么返回值是通信器中所有处理器的所有累积CPU时间的总和。
   *
   */
  double
  cpu_time() const;

  /**
   * 返回最后一圈的CPU时间，单位是秒。计时器不会因为这个函数而停止。
   *
   */
  double
  last_cpu_time() const;

private:
  /**
   * 计时器类存储了两个不同时钟的计时信息：一个挂钟和一个CPU使用时钟。由于处理这两个时钟的逻辑在大多数地方都是相同的，所以我们把每个时钟的相关测量数据收集到这个
   * <code>struct</code>  。      @tparam  clock_type_
   * 其测量值被存储的时钟的类型。这个类应该符合
   * <code>std::chrono</code>
   * 所期望的通常的时钟接口（即正确的别名和一个静态
   * <code>now()</code> 方法）。
   *
   */
  template <class clock_type_>
  struct ClockMeasurements
  {
    /**
     * 存储时钟类型。
     *
     */
    using clock_type = clock_type_;

    /**
     * 提供时钟的时间点类型。
     *
     */
    using time_point_type = typename clock_type::time_point;

    /**
     * 提供的时钟的持续时间类型。
     *
     */
    using duration_type = typename clock_type::duration;

    /**
     * 对应于当前一圈开始的时间点。这是通过调用
     * <code>clock_type::now()</code> 获得的。
     *
     */
    time_point_type current_lap_start_time;

    /**
     * 几圈的累计时间。
     *
     */
    duration_type accumulated_time;

    /**
     * 最后一圈的时间。
     *
     */
    duration_type last_lap_time;

    /**
     * 构造函数。将 <code>current_lap_start_time</code>
     * 设置为当前的时钟时间，持续时间为零。
     *
     */
    ClockMeasurements();

    /**
     * 通过设置 <code>current_lap_start_time</code>
     * 为当前的时钟时间和持续时间为零来重置时钟。
     *
     */
    void
    reset();
  };

  /**
   * 挂钟的别名。
   *
   */
  using wall_clock_type = std::chrono::steady_clock;

  /**
   * CPU时钟的别名。
   *
   */
  using cpu_clock_type = CPUClock;

  /**
   * 墙面时间的测量集合。
   *
   */
  ClockMeasurements<wall_clock_type> wall_times;

  /**
   * 收集CPU时间测量值。
   *
   */
  ClockMeasurements<cpu_clock_type> cpu_times;

  /**
   * 计时器是否目前正在运行。
   *
   */
  bool running;

  /**
   * 用于同步和组合各种时间值的通信器：更多信息请参见相关构造函数的文档。
   *
   */
  MPI_Comm mpi_communicator;

  /**
   * 在 Timer::start() 和 Timer::stop().
   * 中存储墙面时间和CPU时间是否在通信器上同步进行。
   *
   */
  bool sync_lap_times;

  /**
   * 用于并行壁挂时间测量的结构，包括MPI通信器已知的最后一圈时间的最小、最大和所有处理器的平均值。
   *
   */
  Utilities::MPI::MinMaxAvg last_lap_wall_time_data;

  /**
   * 一个用于平行墙时间测量的结构，包括所有进程中记录的最小时间、最大时间以及平均时间，定义为所有单独时间的总和除以MPI_Comm中的MPI进程数，以获得总运行时间。
   *
   */
  Utilities::MPI::MinMaxAvg accumulated_wall_time_data;
};



// TODO: The following class is not thread-safe
/**
 * 该类可用于从程序中不同子段的时间测量中生成格式化的输出。它可以创建几个部分来执行程序的某些方面。一个部分可以被多次输入。通过改变OutputFrequency和OutputType中的选项，用户可以选择是在每次加入一个部分时产生输出，还是只在程序结束时产生。此外，还可以显示CPU时间、墙壁时间或两者。
 * 该类被用于大量的收集定时数据的教程程序中。  step-77
 * 是一个使用它的相对简单的顺序程序的例子。  step-40
 * 和下面提到的其他几个程序将其用于并行计算。
 *
 *  <h3>Usage</h3> 这个类的使用可以是如下的。
 *
 * @code
 * TimerOutput timer (std::cout, TimerOutput::summary,
 *                    TimerOutput::wall_times);
 *
 * timer.enter_subsection ("Setup dof system");
 * setup_dofs();
 * timer.leave_subsection();
 *
 * timer.enter_subsection ("Assemble");
 * assemble_system_1();
 * timer.leave_subsection();
 *
 * timer.enter_subsection ("Solve");
 * solve_system_1();
 * timer.leave_subsection();
 *
 * timer.enter_subsection ("Assemble");
 * assemble_system_2();
 * timer.leave_subsection();
 *
 * timer.enter_subsection ("Solve");
 * solve_system_2();
 * timer.leave_subsection();
 *
 * // do something else...
 * @endcode
 * 当运行时，这个程序将返回这样的输出。
 *
 * @code
 * +---------------------------------------------+------------+------------+
 * | Total wallclock time elapsed since start    |      88.8s |            |
 * |                                             |            |            |
 * | Section                         | no. calls |  wall time | % of total |
 * +---------------------------------+-----------+------------+------------+
 * | Assemble                        |         2 |      19.7s |        22% |
 * | Solve                           |         2 |      3.03s |       3.4% |
 * | Setup dof system                |         1 |      3.97s |       4.5% |
 * +---------------------------------+-----------+------------+------------+
 * @endcode
 * 输出会看到我们进入了装配和解决部分两次，并报告我们在那里花费了多少时间。此外，这个类还测量了TimerOutput对象从开始到终止所花费的总时间。在这种情况下，我们做了很多其他的事情，所以我们测量的函数的时间比例离100%很远。
 *
 *  <h3>Using scoped timers</h3> 上面的方案，你必须有对
 * TimerOutput::enter_subsection() 和 TimerOutput::leave_subsection()
 * 的调用，如果这些调用之间的部分包含 <code>return</code>
 * 语句或可能抛出异常，那就很尴尬了。在这种情况下，很容易忘记我们还是需要以某种方式在某个地方离开这个部分。一个更简单的方法是使用
 * "范围
 * "节。这是一个变量，当你创建它的时候，它就进入了一个部分，当你销毁它的时候，它就离开了这个部分。如果这是一个特定范围（大括号之间的代码块）的局部变量，而你由于
 * <code>return</code>
 * 语句或异常而离开了这个范围，那么这个变量就会被销毁，定时段也会自动离开。因此，我们可以把上面的代码块写成下面的样子，结果完全一样，但现在是异常安全的。
 *
 * @code
 * TimerOutput timer (std::cout, TimerOutput::summary,
 *                    TimerOutput::wall_times);
 *
 * {
 *   TimerOutput::Scope timer_section(timer, "Setup dof system");
 *   setup_dofs();
 * }
 *
 * {
 *   TimerOutput::Scope timer_section(timer, "Assemble");
 *   assemble_system_1();
 * }
 *
 * {
 *   TimerOutput::Scope timer_section(timer, "Solve");
 *   solve_system_1();
 * }
 *
 * {
 *   TimerOutput::Scope timer_section(timer, "Assemble");
 *   assemble_system_2();
 * }
 *
 * {
 *   TimerOutput::Scope timer_section(timer, "Solve");
 *   solve_system_2();
 * }
 *
 * // do something else...
 * @endcode
 *
 *
 *  <h3>Usage in parallel programs using MPI</h3>
 * 在一个建立在MPI上的并行程序中，以如上所示的方式使用该类会导致这样一种情况：每个进程对相应的部分进行计时，然后在最后输出结果的计时信息。这是很烦人的，因为你会得到大量的输出
 *
 * - 每个处理器的一组时序信息。
 * 这可以通过只让一个处理器产生屏幕输出来避免，典型的做法是使用ConditionalOStream类型的对象而不是
 * <code>std::cout</code> 来写入屏幕（例如，见 step-17 、 step-18
 * 、 step-32 和 step-40 ，它们都使用这种方法）。
 * 这样，只有一个处理器输出计时信息，通常是MPI宇宙中的第一个进程。然而，如果你以上面的代码片段为例，想象一下，如果
 * <code>setup_dofs()</code>
 * 在零号处理器上是快的，而在其他处理器中至少有一个是慢的；并且如果
 * <code>assemble_system_1()</code>
 * 做的第一件事是需要所有处理器进行通信，会发生什么。在这种情况下，在零号处理器上，名称为
 * <code>"Setup dof system"</code>
 * 的计时部分将在零号处理器上产生较短的运行时间，而
 * <code> "Assemble"</code> 部分将花费很长的时间：不是因为
 * <code>assemble_system_1()</code>
 * 需要特别长的时间，而是因为在我们计时的处理器上（或者说，我们产生输出的处理器）刚好需要等待很长时间，直到其他处理器最终完成
 * <code>setup_dofs()</code> 并开始参与 <code>assemble_system_1()</code>
 * 。换句话说，所报告的时间是不可靠的，因为它反映了来自其他处理器的运行时间。此外，本节在零号处理器上的运行时间与本节在其他处理器上的运行时间无关，而是与<i>the
 * previous section</i>在另一个处理器上的运行时间有关。
 * 避免这种情况的第一个方法是，在我们开始和停止计时部分之前，在并行代码中引入一个障碍。这可以确保所有进程都在同一个地方，然后计时信息反映了所有处理器的最大运行时间。为了实现这一点，你需要用一个MPI通信器对象来初始化TimerOutput对象，例如像下面的代码。
 *
 * @code
 * TimerOutput timer (MPI_COMM_WORLD,
 *                    pcout,
 *                    TimerOutput::summary,
 *                    TimerOutput::wall_times);
 * @endcode
 * 这里， <code>pcout</code>
 * 是一个ConditionalOStream类型的对象，确保我们只在单个处理器上产生输出。参见
 * step-32 、 step-40 和 step-42
 * 的教程程序，了解该类的这种用法。
 * 应对这个问题的第二个变种是打印更多关于记录时间的信息，以便能够理解这种不平衡，而不需要实际添加障碍。虽然这种方法仍然受到不同MPI进程之间不平衡的影响，但是它的输出不是等级0的任意时间，而是MPI结果的最小、平均和最大，使用来自
 * Utilities::MPI::MinMaxAvg.
 * 的信息，因为数据也配备了达到最小和最大的等级ID，这种方法允许识别某些减速发生在哪个等级。如果可以容忍MPI等级之间从一个部分到下一个部分的一些不平衡，那么这个策略就比障碍变量更有优势，因为它不会在没有必要的地方同步程序，而是试图显示在各个阶段观察到的不平衡。为了使用这个变体，在没有任何本地打印设置和没有通信器的情况下初始化输出对象。
 *
 * @code
 * TimerOutput timer (pcout,
 *                    TimerOutput::never,
 *                    TimerOutput::wall_times);
 * @endcode
 * 然后调用
 *
 * @code
 * timer.print_wall_time_statistics(MPI_COMM_WORLD);
 * @endcode
 * 在适当的地方。这里，输出被写入传递给构造函数的
 * <code>pcout</code>
 * 类型的ConditionalOStream对象中，确保信息只被打印一次。参见
 * step-67
 * 中关于这个变体的使用实例。除了所有MPI等级的基本最小、平均和最大时间外，
 * TimerOutput::print_wall_time_statistics()
 * 函数还需要第二个参数来指定定量的输出，例如，最慢和最快等级的10/%所花的时间，以获得对统计分布的额外洞察力。
 *
 *
 * @ingroup utilities
 *
 *
 */
class TimerOutput
{
public:
  /**
   * 帮助类，在TimerOutput中进入/退出部分，构建一个简单的基于范围的对象。这个类的目的在TimerOutput的文档中有所解释。
   *
   */
  class Scope
  {
  public:
    /**
     * 在定时器中输入给定的部分。在调用stop()或析构器运行时自动退出。
     *
     */
    Scope(dealii::TimerOutput &timer_, const std::string &section_name);

    /**
     * 解构器调用stop()
     *
     */
    ~Scope();

    /**
     * 如果你想在执行析构器之前退出作用域，请调用这个函数。
     *
     */
    void
    stop();

  private:
    /**
     * 对TimerOutput对象的引用
     *
     */
    dealii::TimerOutput &timer;

    /**
     * 我们需要退出的部分的名称
     *
     */
    const std::string section_name;

    /**
     * 我们是否还需要退出我们所处的部分？
     *
     */
    bool in;
  };

  /**
   * 一个枚举数据类型，它描述了是否在我们每次退出一个部分时产生输出，只是在最后，两者都是，或者永远不产生。
   *
   */
  enum OutputFrequency
  {
    /**
     * 每次调用后生成输出。
     *
     */
    every_call,
    /**
     * 在结尾处生成摘要输出。
     *
     */
    summary,
    /**
     * 在每次调用后和最后的总结中都生成输出。
     *
     */
    every_call_and_summary,
    /**
     * 不产生任何输出。
     *
     */
    never
  };

  /**
   * 一个枚举数据类型，描述从定时器获取数据时要返回的数据类型。
   *
   */
  enum OutputData
  {
    /**
     * 输出CPU时间。
     *
     */
    total_cpu_time,
    /**
     * 输出挂钟时间。
     *
     */
    total_wall_time,
    /**
     * 输出调用次数。
     *
     */
    n_calls
  };

  /**
   * 一个枚举数据类型，描述了每当我们生成输出时，是否显示CPU时间、壁挂时间，或者同时显示CPU和壁挂时间。
   *
   */
  enum OutputType
  {
    /**
     * 输出CPU时间。
     *
     */
    cpu_times,
    /**
     * 输出挂钟时间。
     *
     */
    wall_times,
    /**
     * 在不同的表格中同时输出CPU和挂钟时间。
     *
     */
    cpu_and_wall_times,
    /**
     * 在一个表中输出CPU和墙面时钟时间。
     *
     */
    cpu_and_wall_times_grouped
  };

  /**
   * 构造函数。      @param  stream 流（类型为 std::ostream)
   * ，输出被写入其中。    @param  output_frequency
   * 一个变量，表示何时将输出写入给定的流。    @param
   * output_type
   * 一个变量，表示输出应该代表哪种时间（CPU或墙面时间）。
   *
   */
  TimerOutput(std::ostream &        stream,
              const OutputFrequency output_frequency,
              const OutputType      output_type);

  /**
   * 构造函数。      @param  stream
   * 输出被写入的流（ConditionalOstream类型）。    @param
   * output_frequency 表示何时将输出写入给定流的一个变量。
   * @param  output_type
   * 一个变量，表示输出应该代表哪种时间（CPU或墙面时间）。
   *
   */
  TimerOutput(ConditionalOStream &  stream,
              const OutputFrequency output_frequency,
              const OutputType      output_type);

  /**
   * 构造函数，将MPI通信器作为输入。这样构造的定时器将把MPI网络中所有处理器的CPU时间加起来计算CPU时间，或者取所有处理器的最大值，取决于
   * @p output_type
   * 的值。关于这个构造函数的原理和一个例子，请参见该类的文档。
   * @param  mpi_comm
   * 一个MPI通信器，我们应该在其上积累或以其他方式同步我们在每个MPI进程上产生的计时信息。
   * @param  stream 流（类型为 std::ostream) ，输出被写入其中。
   * @param  output_frequency 表示何时将输出写入给定流的变量。
   * @param  output_type
   * 一个变量，表示输出应该代表哪种时间（CPU或墙面时间）。在这个并行的上下文中，当这个参数选择CPU时间时，那么时间会在参与MPI通信器的所有进程中累积。如果这个参数选择了墙时间，那么报告的时间就是这部分所有处理器运行时间的最大值。后者是通过在启动和停止每个部分的定时器之前放置一个
   * <code>MPI_Barrier</code> 调用来计算的。
   *
   */
  TimerOutput(const MPI_Comm &      mpi_comm,
              std::ostream &        stream,
              const OutputFrequency output_frequency,
              const OutputType      output_type);

  /**
   * 构造函数需要一个MPI通信器作为输入。这样构造的定时器将把MPI网络中所有处理器的CPU时间加起来计算CPU时间，或者取所有处理器的最大值，这取决于
   * @p output_type
   * 的值。关于这个构造函数的原理和一个例子，请参见该类的文档。
   * @param  mpi_comm
   * 一个MPI通信器，我们应该在其上积累或以其他方式同步我们在每个MPI进程上产生的时间信息。
   * @param  stream 输出被写入的流（ConditionalOstream类型）。
   * @param  output_frequency
   * 表示何时将输出写入给定流的一个变量。    @param
   * output_type
   * 一个变量，表示输出应该代表哪种时间（CPU或墙面时间）。在这个并行的上下文中，当这个参数选择CPU时间时，那么时间会在参与MPI通信器的所有进程中累积。如果这个参数选择了墙时间，那么报告的时间就是这部分所有处理器运行时间的最大值。后者是通过在启动和停止每个部分的定时器之前放置一个
   * <code>MPI_Barrier</code> 调用来计算的）。
   *
   */
  TimerOutput(const MPI_Comm &      mpi_comm,
              ConditionalOStream &  stream,
              const OutputFrequency output_frequency,
              const OutputType      output_type);

  /**
   * 销毁器。如果写摘要输出的选项被设置，则调用print_summary()。
   *
   */
  ~TimerOutput();

  /**
   * 通过给定一个字符串的名称来打开一个部分。如果该名称已经存在，则再次进入该部分并累积次数。
   *
   */
  void
  enter_subsection(const std::string &section_name);

  /**
   * 离开一个部分。如果没有给定名称，则留下最后输入的部分。
   *
   */
  void
  leave_subsection(const std::string &section_name = "");

  /**
   * 获得一张地图，上面有每个小节的指定类型的收集的数据。
   *
   */
  std::map<std::string, double>
  get_summary_data(const OutputData kind) const;

  /**
   * 打印一个格式化的表格，总结各部分所消耗的时间。
   *
   */
  void
  print_summary() const;

  /**
   * 打印一个格式化的表格，总结各部分所消耗的墙体时间，使用各部分时间的最小值、平均值和最大值以及达到最小值和最大值的MPI等级的统计数据。注意这个调用只有在构造TimerOutput对象时没有MPI_Comm参数时才能提供有用的信息，以使各个部分的运行不受障碍物的干扰。
   * 可选的参数`quantile`允许在运行时间的分布方面为输出增加两列。如果quantile=0.1，除了最小值和最大值外，还打印10%的最低数据的值和等级，以及分布函数的90%的值和等级。quantile
   * "的值需要在0（除了最小和最大之外不打印任何量值）和0.5（当给出中位数时）之间。
   *
   */
  void
  print_wall_time_statistics(const MPI_Comm &mpi_comm,
                             const double    print_quantile = 0.) const;

  /**
   * 通过调用这个函数，所有的输出都可以被禁用。如果想以灵活的方式控制输出，而又不想在程序中加入大量的<tt>if</tt>条款，这个函数和enable_output()就很有用。
   *
   */
  void
  disable_output();

  /**
   * 如果之前用disable_output()禁用了该类的输出，该函数重新启用。如果想以灵活的方式控制输出，而不在程序中加入大量的<tt>if</tt>子句，这个函数与disable_output()一起使用是很有用的。
   *
   */
  void
  enable_output();

  /**
   * 重置记录的定时信息。
   *
   */
  void
  reset();

private:
  /**
   * 什么时候向输出流输出信息。
   *
   */
  OutputFrequency output_frequency;

  /**
   * 是否显示CPU时间、壁挂时间，或同时显示CPU和壁挂时间。
   *
   */
  OutputType output_type;


  /**
   * 一个用于整体运行时间的计时器对象。如果我们使用的是MPI，这个计时器也会在所有MPI进程中累积。
   *
   */
  Timer timer_all;

  /**
   * 一个结构，将我们收集的关于每个部分的所有信息分组。
   *
   */
  struct Section
  {
    Timer        timer;
    double       total_cpu_time;
    double       total_wall_time;
    unsigned int n_calls;
  };

  /**
   * 一个所有章节及其信息的列表。
   *
   */
  std::map<std::string, Section> sections;

  /**
   * 我们要输出的流对象。
   *
   */
  ConditionalOStream out_stream;

  /**
   * 一个布尔变量，用于设置该类的输出目前是开还是关。
   *
   */
  bool output_is_enabled;

  /**
   * 一个已经进入而没有退出的部分的列表。该列表按照进入章节的顺序保存，但是如果给leave_subsection()函数一个参数，元素可以在中间被移除。
   *
   */
  std::list<std::string> active_sections;

  /**
   * mpi通信器
   *
   */
  MPI_Comm mpi_communicator;

  /**
   * 一个锁，确保这个类在与多个线程一起使用时也能给出合理的结果。
   *
   */
  Threads::Mutex mutex;
};



 /* ---------------- inline functions ----------------- */ 


inline void
Timer::restart()
{
  reset();
  start();
}



inline const Utilities::MPI::MinMaxAvg &
Timer::get_last_lap_wall_time_data() const
{
  return last_lap_wall_time_data;
}



inline const Utilities::MPI::MinMaxAvg &
Timer::get_accumulated_wall_time_data() const
{
  return accumulated_wall_time_data;
}



template <class StreamType>
inline void
Timer::print_last_lap_wall_time_data(StreamType &stream) const
{
  const Utilities::MPI::MinMaxAvg &statistic = get_last_lap_wall_time_data();
  stream << statistic.max << " wall,"
         << " max @" << statistic.max_index << ", min=" << statistic.min << " @"
         << statistic.min_index << ", avg=" << statistic.avg << std::endl;
}



template <class StreamType>
inline void
Timer::print_accumulated_wall_time_data(StreamType &stream) const
{
  const Utilities::MPI::MinMaxAvg &statistic = get_accumulated_wall_time_data();
  stream << statistic.max << " wall,"
         << " max @" << statistic.max_index << ", min=" << statistic.min << " @"
         << statistic.min_index << ", avg=" << statistic.avg << std::endl;
}



inline TimerOutput::Scope::Scope(dealii::TimerOutput &timer_,
                                 const std::string &  section_name_)
  : timer(timer_)
  , section_name(section_name_)
  , in(true)
{
  timer.enter_subsection(section_name);
}



inline void
TimerOutput::Scope::stop()
{
  if (!in)
    return;
  in = false;

  timer.leave_subsection(section_name);
}


DEAL_II_NAMESPACE_CLOSE

#endif


