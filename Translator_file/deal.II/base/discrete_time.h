//include/deal.II-translator/base/discrete_time_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

#ifndef dealii_discrete_time_h
#define dealii_discrete_time_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 这个类提供了一种方法来跟踪随时间变化的模拟的时间。它管理着从开始时间
 * $T_{\text{start}}$ 到结束时间 $T_{\text{end}}$
 * 的向前迈进。它还允许在仿真过程中调整时间步长。该类提供了必要的接口，可以被纳入任何时间依赖性的仿真。这个类的用法在
 * step-19  和  step-21  中演示。
 * 该类提供了一些保证在任何时候都是真实的不变量。
 * 当前的模拟时间是在开始时间和结束时间之间的封闭区间内（
 * $T_{\text{start}} \le t \le T_{\text{end}}$
 * ）。每当时间递增时，步长是正的（ $dt > 0$ ）。
 * 换句话说，时间以严格的升序前进（  $m < n \Leftrightarrow
 * t_m < t_n$  ）。
 * 这个类所遵循的模型是，人们通过构造函数或使用set_desired_next_step_size()函数设置期望的*时间步长。这个步长将用于接下来所有对advance_time()函数的调用，但在模拟结束时可以稍作调整，以确保模拟时间与结束时间完全一致。这种调整对于以下原因是很有用的。
 * 假设你通过使用`for`循环来循环所有的时间步骤
 *
 * @code
 * for (DiscreteTime time(0., 1., 0.3);
 *      time.is_at_end() == false;
 *      time.advance_time())
 * {
 *   // Insert simulation code here
 * }
 * @endcode
 * 或者，如果你更喜欢这种方式，可以使用等效的 "while
 * "循环。
 *
 * @code
 * DiscreteTime time(0., 1., 0.3);
 * while (time.is_at_end() == false)
 * {
 *   // Insert simulation code here
 *
 *   time.advance_time();
 * }
 * @endcode
 *
 * 在上面的例子中，时间从 $T_{\text{start}} = 0$ 开始，直到
 * $T_{\text{end}}=1$  。假设时间步数  $dt = 0.3$
 * 在循环内没有被修改，时间从  $t = 0$  推进到  $t = 0.3$  ,
 * $t = 0.6$  ,  $t = 0.9$  ，最后它在  $t = 1.0$
 * 到达结束时间。在这里，最后的步长需要从它的理想值0.3减少到
 * $dt = 0.1$
 * ，以确保我们正好在指定的结束时间完成模拟。事实上，你应该假设不仅最后的时间步长可以调整，而且以前的时间步长也可以调整。
 *
 * - 例如，这个类可能会擅自将时间步长的减少分散到几个时间步长中，并将时间从 $t=0$ ，增加到 $0.3$ ， $0.6$ ， $0.8$ ，最后到 $t=T_{\text{end}}=1$ ，以避免时间步长从一个步骤到另一个步骤的变化太大。
 * 另一种需要调整时间步长的情况（这次是调整为稍大的数值）是如果时间增量刚好低于最终时间。例如，想象一下，与上述情况类似，但结束时间不同。
 *
 * @code
 * for (DiscreteTime time(0., 1.21, 0.3);
 *      time.is_at_end() == false;
 *      time.advance_time())
 * {
 *   // Insert simulation code here
 * }
 * @endcode
 * 这里，从 $t=0.9$ 到 $t=1.2$ 的时间步长刚好低于最终时间
 * $T_{\text{end}}=1.21$  。而不是用一个长度为 $dt=0.01$
 * 的非常小的步骤来跟进，该类人将最后的时间步骤（或最后的时间步骤）稍微拉长，以达到所需的结束时间。
 * 上面的例子清楚地表明，给这个类的时间步长只是一个期望的*步长。你可以使用get_next_step_size()函数查询实际的时间步长。
 *
 *  # ### 时间步长的细节
 * 因为在我们的模拟中，时间是以离散的方式向前推进的，所以我们需要讨论如何增加时间。在时间步进过程中，我们每一步都会进入两个独立的交替状态。
 * *快照**阶段（*当前**阶段，*一致**阶段 阶段）。)
 * 在这部分算法中，我们处于 $t = t_n$
 * ，所有的模拟量（位移、应变、温度等）都是 $t = t_n$
 * 的最新情况。在这个阶段，当前时间*是指 $t_n$
 * ，下一个时间*是指 $t_{n+1}$ ，上一个时间*是指 $t_{n-1}$
 * 。其他有用的符号量是下一个*时间步长  $t_{n+1}
 *
 * - t_n$  和上一个*时间步长  $t_n
 *
 * - t_{n-1}$
 * 。在这个阶段，使用用户代码中的打印命令生成文本输出是一个完美的场合。此外，还可以在这里准备后处理输出，然后可以通过可视化程序如
 * "Tecplot"、"Paraview "和 "VisIt
 * "查看。此外，在快照阶段，代码可以评估前一个步骤的质量，并决定是否要增加或减少时间步长。在这里可以通过调用set_desired_next_step_size()来修改下一个时间步骤的步骤大小。更新**阶段（*过渡**阶段，*不一致**阶段
 * 阶段）。) 在这部分程序中，模拟的内部状态被从 $t_n$
 * 更新到 $t_{n+1}$
 * 。所有的变量都需要逐一更新，步数被增加，时间被增加
 * $dt = t_{n+1}
 *
 * - t_n$
 * ，时间积分算法被用来更新其他模拟量。在这个阶段的中间，一些变量已经被更新到
 * $t_{n+1}$ ，但其他变量仍然代表它们在 $t_n$
 * 的值。因此，我们把这个阶段称为不一致阶段，要求在这个阶段内不发生与状态变量相关的后处理输出。状态变量，即那些与时间、解场和任何内部变量有关的变量，并不是同步的，然后逐一得到更新。一般来说，更新变量的顺序是任意的，但如果它们之间存在相互依赖关系，就应该注意一些。例如，如果某个变量如
 * $x$ 依赖于另一个变量如 $y$ 的计算，那么 $y$ 必须在 $x$
 * 被更新之前被更新。
 * 问题是，在更新状态量之前，时间是否应该被递增。存在多种可能性，取决于程序和配方的要求，可能还有程序员的偏好。
 * 时间在*其余的更新之前被递增。在这种情况下，即使时间被递增到
 * $t_{n+1}$ ，也不是所有变量都被更新。在这个更新阶段，
 * $dt$
 * 等于previous*时间步长。Previous*意味着它是指之前执行的`advance_time()`命令的
 * $dt$
 * 。在下面的示例代码中，我们假设`a`和`b`是两个需要在这个时间步长中更新的状态变量。
 * @code
 *     time.advance_time();
 *     new_a = update_a(a, b, time.get_previous_step_size());
 *     b = update_b(a, b, time.get_previous_step_size());
 *     a = new_a;
 * @endcode
 * 这里，代码开始时是一致的状态，但是一旦调用advance_time()，时间变量、`a`和`b`就不再相互一致，直到最后一条语句之后。在这一点上，这些变量又都是一致的。
 * 在*所有变量已经更新为 $t_{n+1}$ 之后，时间从 $t_n$
 * 递增到 $t_{n+1}$  。在更新阶段， $dt$
 * 被表示为下一个*时间步长。下一个*意味着 $dt$
 * 的步长对应于随后发生的`advance_time()`命令。
 * @code
 *     new_a = update_a(a, b, time.get_next_step_size());
 *     b = update_b(a, b, time.get_next_step_size());
 *     a = new_a;
 *     time.advance_time();
 * @endcode
 * 时间是在其他更新的中间递增的。在这种情况下， $dt$
 * 将对应xt*或previous*，取决于它是在调用`advance_time()'之前还是之后使用。
 * @code
 *     new_a = update_a(a, b, time.get_next_step_size());
 *     time.advance_time();
 *     b = update_b(a, b, time.get_previous_step_size());
 *     a = new_a;
 * @endcode
 *
 * 需要注意的是，在更新阶段， $dt$
 * 是指*下一个**或*上一个**时间步长，这取决于是否已经调用了advance_time()。当前*时间步长的概念定义不明确。事实上，在更新阶段，每个变量的定义都取决于它是否已经被更新，因此被称为*不一致的阶段**。
 * 下面的代码片段显示了在一个完整的时间依赖性仿真的背景下，快照阶段和更新阶段的代码部分。这段代码遵循了教程实例中的编码惯例。请注意，尽管这个例子是以
 * "for "循环的格式写的，但它也可以等同于写成 "while "或
 * "do while "循环（如 step-21 所示）。
 *
 * @code
 * // pre-processing/setup stage {
 * make_grid();
 * setup_system();
 * for (DiscreteTime time(0., 1., 0.1);  // } end pre-processing/setup stage
 *    time.is_at_end() == false;
 *    time.advance_time())             // part of the update stage, runs at
 *                                     // the end of the loop body
 * {
 * // snapshot stage {
 * const double time_of_simulation = time.get_next_time();
 * const double timestep_size      = time.get_next_step_size();
 *
 * std::cout
 *   << "Timestep: " << time.get_step_number() << "
 *
 * -- "
 *   << "Solving for the solution at "
 *   << "t = " << time_of_simulation << " with "
 *   << "dt = " << timestep_size << "." << std::endl;
 * // } end snapshot stage
 *
 * // update stage {
 * assemble_system(time_of_simulation, timestep_size);
 * solve();
 * update_solutions();
 * // } end update stage
 *
 * // snapshot stage {
 * output_results(time_of_solution);
 *
 * // propose a new timestep size if need be
 * // time.set_desired_next_step_size(...);
 * // } end snapshot stage
 * }
 * @endcode
 * step-19
 * 中的`run()`函数显示了一个非常类似的例子，其中对advance_time()的调用结束了更新阶段，随后用当时的时间生成了图形输出。
 *
 *
 */
class DiscreteTime
{
public:
  /**
   * 构造函数。      @param[in]  start_time 仿真开始时的时间。
   * @param[in]  end_time 仿真结束时的时间。      @param[in]
   * desired_start_step_size
   * 用于第一步时间递增的预期步长。不保证这个值会被实际用作第一步的大小，这一点在介绍中已经讨论过。
   * @pre   @p desired_start_step_size  必须为非负数。
   * @note   @p desired_start_step_size
   * 是一个可选的参数。如果没有提供或指定为零，表明时间步长的所需尺寸将在代码中的不同位置计算。在这种情况下，创建的对象不能增加时间，直到通过调用set_desired_next_step_size()改变步长大小。
   *
   */
  DiscreteTime(const double start_time,
               const double end_time,
               const double desired_start_step_size = 0.);

  /**
   * 返回当前时间。
   *
   */
  double
  get_current_time() const;

  /**
   * 返回如果我们将时间提前一步，将达到的下一个时间。
   * @note 如果模拟到了结束时间，该方法返回结束时间。
   *
   */
  double
  get_next_time() const;

  /**
   * 返回我们上次调用`advance_time()`之前的时间。
   * @note 如果模拟处于开始时间，此方法返回开始时间。
   *
   */
  double
  get_previous_time() const;

  /**
   * 返回开始时间。
   *
   */
  double
  get_start_time() const;

  /**
   * 返回时间区间的结束时间。
   * 最后的时间步骤正好在这一点上结束。这个精确的浮点数是非常重要的，因为它允许我们将当前时间与结束时间进行等价比较，并决定我们是否已经到达了模拟的终点。
   *
   */
  double
  get_end_time() const;

  /**
   * 返回是否还没有发生任何步骤。
   *
   */
  bool
  is_at_start() const;

  /**
   * 返回时间是否已经到达结束时间。
   *
   */
  bool
  is_at_end() const;

  /**
   * 返回从当前时间步长到下一个时间步长的大小。正如在类的介绍中所讨论的，这是实际的*时间步长，可能与在构造函数中或通过set_desired_next_step_size()函数设置的desired*时间步长不同。
   * @note 如果仿真处于结束时间，该方法返回0。
   *
   */
  double
  get_next_step_size() const;

  /**
   * 返回上一步的步长。
   * @note  如果仿真处于开始时间，该方法返回0。
   *
   */
  double
  get_previous_step_size() const;

  /**
   * 返回仿真时间被增加的次数。
   * 当仿真处于开始时间时，返回0。
   *
   */
  unsigned int
  get_step_number() const;

  /**
   * 设置下一个时间步长的期望*值。通过调用这个方法，我们表明下次调用advance_time()时，我们希望用
   * @p time_step_size 来推进仿真时间。
   * 但是，如果步长过大，导致下一次仿真时间超过结束时间，步长就会被截断。
   * 此外，如果步长使下一次模拟时间接近结束时间（但略微低于结束时间），步长将被调整，使下一次模拟时间与结束时间完全一致。
   *
   */
  void
  set_desired_next_step_size(const double time_step_size);

  /**
   * 设置下一个时间步长的实际*值。通过调用这个方法，我们表明下次调用advance_time()时，将使用
   * @p time_step_size 来推进仿真时间。
   * @note
   * set_next_step_size()和set_desired_next_step_size()的区别在于，前者完全使用提供的
   * $dt$ 而不做任何调整，但如果 $dt$
   * 不在可接受的范围内，就会产生一个错误（在调试模式）。
   * 一般来说，set_desired_next_step_size()是首选方法，因为它可以根据
   * $T_{\text{end}}$  智能地调整 $dt$  。    @pre   $0 < dt \le
   * T_{\text{end}}
   *
   * - t$  。
   *
   */
  void
  set_next_step_size(const double time_step_size);

  /**
   * 根据当前步骤的数值推进当前时间。
   * 如果你想调整下一个时间步长，请在调用此方法之前调用set_desired_next_step_size()方法。
   * 如果你重复调用这个函数，时间会以相同的步长增加，直到达到结束时间。参见set_desired_next_step_size()的文档，了解自动调整步长的规则。
   * @pre
   * 当前时间必须小于结束时间。如果对象已经到了结束时间，就不能推进时间。创建这条规则是为了避免在循环内调用
   * advance_time() 时产生无限循环。      @pre
   * 时间步长必须为非零。如果当前的步长为零，请在调用advance_time()之前通过调用set_desired_next_step_size()来改变它。
   *
   */
  void
  advance_time();

  /**
   * 设置当前时间等于开始时间，并将步长设置为初始步长。
   *
   */
  void
  restart();

private:
  /**
   * 时间间隔的起始时间。
   *
   */
  double start_time;

  /**
   * 时间间隔的结束。
   *
   */
  double end_time;

  /**
   * 当前的时间。
   *
   */
  double current_time;

  /**
   * 下一个步骤的时间。
   * @note
   * 在内部，下一个模拟时间被存储，而不是当前的步长。例如，当方法set_desired_next_step_size()被调用时，它计算出适当的下一个模拟时间并存储起来。当advance_time()被调用时，current_time被next_time所取代。这种对内部状态的选择使得代码更加简单，并确保当我们在最后一步调用advance_time()时，时间的浮点值与结束时间完全一致。
   *
   */
  double next_time;

  /**
   * 之前的时间。
   *
   */
  double previous_time;

  /**
   * 第一个步骤的大小。
   *
   */
  double start_step_size;

  /**
   * 步数，即模拟时间被递增的次数。
   *
   */
  unsigned int step_number;
};


 /*---------------------- Inline functions ------------------------------*/ 


inline double
DiscreteTime::get_start_time() const
{
  return start_time;
}



inline double
DiscreteTime::get_end_time() const
{
  return end_time;
}



inline bool
DiscreteTime::is_at_start() const
{
  return step_number == 0;
}



inline bool
DiscreteTime::is_at_end() const
{
  return current_time == end_time;
}



inline double
DiscreteTime::get_next_step_size() const
{
  return next_time - current_time;
}



inline double
DiscreteTime::get_previous_step_size() const
{
  return current_time - previous_time;
}



inline double
DiscreteTime::get_current_time() const
{
  return current_time;
}



inline double
DiscreteTime::get_next_time() const
{
  return next_time;
}



inline double
DiscreteTime::get_previous_time() const
{
  return previous_time;
}



inline unsigned int
DiscreteTime::get_step_number() const
{
  return step_number;
}


DEAL_II_NAMESPACE_CLOSE

#endif


