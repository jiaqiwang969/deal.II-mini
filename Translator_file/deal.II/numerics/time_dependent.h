//include/deal.II-translator/numerics/time_dependent_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#ifndef dealii_time_dependent_h
#  define dealii_time_dependent_h


 /*----------------------------   time-dependent.h ---------------------------*/ 


#  include <deal.II/base/config.h>

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/smartpointer.h>
#  include <deal.II/base/subscriptor.h>

#  include <utility>
#  include <vector>

DEAL_II_NAMESPACE_OPEN

// forward declarations
#  ifndef DOXYGEN
class TimeStepBase;
template <typename number>
class Vector;
template <int dim, int spacedim>
class Triangulation;
#  endif

/**
 * 这个类为时间相关问题提供了一个抽象接口，它解决了这一类问题中最令人讨厌的一些方面：数据管理。这些问题经常需要大量的计算机资源，最明显的是计算时间、主内存和磁盘空间。减少主内存往往是最迫切的需求，不过实现它的方法几乎总是相当混乱，很快就会导致代码在散落在程序各处的地方存储和重新加载数据，而且有时变得无法维护。本类试图提供一个更加结构化的接口，尽管很简单，这是在我搞了几个月的波浪方程模拟之后在我的脑海中出现的。
 * 这个类的设计主要是为解决与时间有关的偏微分方程而定制的，在这些方程中，每两个时间段的计算网格可能不同，而且与这个类的开销相比，每个时间段的计算需要相当长的时间。由于本类中没有提到问题的类别，因此它并不局限于PDEs，尽管如此，大型普通矩阵微分方程的求解器似乎有可能成功地使用相同的设置，因此本类。
 *
 *  <h3>Overview</h3>
 * 使用时间步进方案的时间相关问题求解器的一般结构大约如下：我们有一个时间步进对象的集合，我们随后在其上求解我们的问题。为了做到这一点，我们需要知道之前的零个或几个时间步长的数据（当使用单步或多步方法时，即），也许还需要一些时间步长的数据（例如，这些时间步长的计算网格）。根据有关问题，可以在所有时间步上进行第二次循环，解决一个对偶问题，循环可以向前运行（每个时间步的一个对偶问题）或向后运行（使用一个全局对偶问题）。在这些循环中的一个或使用一个单独的循环，可以计算误差估计值并完善网格。这些循环中的每一个都是由一个为下一个循环准备每个时间步对象的调用开始的，然后才实际开始循环本身。
 * 我们将用术语 "扫频
 * "来表示所有这些循环的完整集合。由于这个库主要是关于自适应方法的，所以在一个扫频中的最后一个循环可能会产生精炼的网格，我们将对这些精炼的网格进行另一次扫频。因此，一个总的运行往往是几个扫频的序列。因此，全局设置看起来像这样。
 *
 * @verbatim
 *  for sweep=0 to n_sweeps-1
 *  {
 *    for i=0 to n_timesteps-1
 *      initialize timestep i for this sweep, e.g. for setting up
 *      data structures, creating temporary files, etc.
 *
 *    for i=0 to n_timesteps-1
 *      prepare timestep i for loop 0
 *    for i=0 to n_timesteps-1
 *      perform loop 0 on timestep i   (e.g. solve primal problem)
 *
 *    for i=0 to n_timesteps-1
 *      prepare timestep i for loop 1
 *    for i=0 to n_timesteps-1
 *      perform loop 1 on timestep i   (e.g. solve dual problem)
 *
 *    for i=0 to n_timesteps-1
 *      prepare timestep i for loop 2
 *    for i=0 to n_timesteps-1
 *      perform loop 2 on timestep i   (e.g. compute error information)
 *
 *    ...
 *
 *    for i=0 to n_timesteps-1
 *      notify timestep i of the end of the sweep, e.g. for cleanups,
 *      deletion of temporary files, etc.
 *  }
 * @endverbatim
 * 用户可以指定一个循环应向前或向后运行（例如，后者是解决全局对偶问题所需要的）。
 * 从全局的角度来看，我们注意到，当一个循环访问一个时间段时（例如，解决原始或对偶问题，或计算误差信息），我们需要关于这个时间段、之前的一个或多个时间段以及未来的零个或多个时间段的信息。然而，往往不需要知道这些时间步骤的所有信息，而且在不再需要数据的时候，在第一个可能的时间删除数据往往是一种计算要求。同样，数据也应该在可能的最晚时间被重新加载。
 * 为了促进这些原则，开发了唤醒和让一个时间步长对象睡眠的概念。假设我们有一个时间步长方案，需要向前看一个时间步长，并且需要最后两个时间步长的数据，下面的伪代码描述了当我们从时间步长
 * @p n-1移到时间步长 @p n:
 * 时，这个类的中心循环函数将做什么
 *
 * @verbatim
 * wake up timestep n+1 with signal 1
 * wake up timestep n with signal 0
 * do computation on timestep n
 * let timestep n sleep with signal 0
 * let timestep n-1 sleep with signal 1
 * let timestep n-2 sleep with signal 2
 *
 * move from n to n+1
 * @endverbatim
 * 这里的信号号表示被发送信号的时间段到进行计算的时间段的距离。对
 * @p  wake_up和 @p sleep
 * 函数的信号0的调用原则上可以被吸收到进行计算的函数中；然而，我们使用这些多余的信号是为了将计算和数据管理相互分离，允许将所有围绕网格管理、数据重载和存储的东西放到一组函数中，而将计算放到另一组中。
 * 在上面的例子中，可能的动作是：时间步数<tt>n+1</tt>重建计算网格（有一个专门的类可以为你做这个）；时间步数
 * @p n
 * 建立矩阵并将解决方案向量设置为合适的大小，也许使用一个初始猜测；然后它进行计算；然后它删除矩阵，因为它们不被后续时间步数需要；时间步数
 * @p n-1
 * ]删除那些只被前面一个时间步骤所需要的数据向量；时间步骤
 * @p n-2
 * 删除剩余的向量并删除计算网格，在某处存储如何最终重建它的信息。
 * 从上面给定的草图中可以看出，每个时间步骤对象都看到以下的事件序列。
 *
 * @verbatim
 * wake up with signal 1
 * wake up with signal 0
 * do computation
 * sleep with signal 0
 * sleep with signal 1
 * sleep with signal 2
 * @endverbatim
 * 这个模式对每个扫频中的每个循环都是重复的。
 * 对于每个扫频中的不同循环，可以分别选择向前看（即对
 * @p wake_up 函数的最大信号数）和向后看（即对 @p sleep
 * 函数的最大信号数）的时间步数。例如，在计算误差估计时，通常只需要往后看一个时间步长（在某些情况下，甚至可以完全不往前看或往后看，在这种情况下，将只发送零号信号），而对于时间步长法，则需要至少往后看一个。
 * 最后，关于前看和后看的方向的说明：前看总是指循环运行的方向，即对于向前运行的循环，
 * @p wake_up
 * 是为时间值大于先前计算的时间段对象而调用的，而 @p
 * sleep
 * 是为时间值较低的时间段调用的。如果循环运行的方向相反，例如在解决一个全局对偶问题时，这个顺序是相反的。
 *
 *  <h3>Implementation</h3>
 * 使用这个类的程序的主循环通常会像下面这样，从一个没有作为库的一部分分发的应用程序中修改而来。
 *
 * @code
 * template <int dim>
 * void TimeDependent_Wave<dim>::run_sweep (const unsigned int sweep_no)
 * {
 *   start_sweep (sweep_no);
 *
 *   solve_primal_problem ();
 *
 *   if (compute_dual_problem)
 *     solve_dual_problem ();
 *
 *   postprocess ();
 *
 *   if (sweep_no != number_of_sweeps-1)
 *     refine_grids ();
 *
 *   write_statistics ();
 *
 *   end_sweep ();
 * }
 *
 *
 *
 *
 *
 * template <int dim>
 * void WaveProblem<dim>::run ()
 * {
 *   for (unsigned int sweep=0; sweep<number_of_sweeps; ++sweep)
 *     timestep_manager.run_sweep (sweep);
 * }
 * @endcode
 * 这里， @p timestep_manager
 * 是一个类型为TimeDependent_Wave<dim>的对象，它是一个从TimeDependent派生的类。
 * @p start_sweep,   @p  solve_primal_problem,  @p solve_dual_problem,   @p
 * postprocess  和  @p  end_sweep
 * 是继承自这个类的函数。它们都在这个对象中的所有时间段上做循环，并在每个对象上调用各自的函数。例如，这里有两个函数是由库实现的。
 *
 * @code
 * void TimeDependent::start_sweep (const unsigned int s)
 * {
 *   sweep_no = s;
 *
 *  // reset the number each time step has, since some time steps might have
 *  // been added since the last time we visited them.
 *  // also set the sweep we will process in the sequel
 *   for (unsigned int step=0; step<timesteps.size(); ++step)
 *     {
 *       timesteps[step]->set_timestep_no (step);
 *       timesteps[step]->set_sweep_no (sweep_no);
 *     }
 *
 *   for (unsigned int step=0; step<timesteps.size(); ++step)
 *     timesteps[step]->start_sweep ();
 * }
 *
 *
 *
 *
 * void
 * TimeDependent::solve_primal_problem ()
 * {
 *   do_loop([](TimeStepBaseconst time_step)
 *             { time_step->init_for_primal_problem(); },
 *           [](TimeStepBaseconst time_step)
 *             { time_step->solve_primal_problem(); },
 *           timestepping_data_primal,
 *           forward);
 * }
 * @endcode
 * 后一个函数相当清楚地显示了大多数循环的调用方式（
 * @p solve_primal_problem,   @p solve_dual_problem,   @p postprocess,   @p
 * refine_grids和 @p write_statistics
 * 都有这种形式，其中后两个函数给出了派生时间段类的函数，而不是来自基类）。函数
 * TimeStepBase::init_for_primal_problem
 * 和该类定义的其他操作的相应函数仅用于存储目前执行的循环将进行的操作类型。
 * 可以看出，大部分工作是由该类的 @p do_loop
 * 函数完成的，它接收两个函数的地址，这两个函数用于初始化循环的所有时间步长对象和实际执行一些动作。下一个参数给出了一些关于前视和后视的信息，最后一个参数表示了循环运行的方向。
 * 使用lambda函数可以做一些巧妙的技巧，就像下面这种情况下的函数
 * @p refine_grids: 。
 *
 * @code
 * ...
 * compute the thresholds for refinement
 * ...
 *
 * do_loop([](TimeStepBase_Tria<dim>const time_step)
 *           { time_step->init_for_refinement(); },
 *         [=](TimeStepBase_Wave<dim>const time_step)
 *           {
 *             time_step->solve_primal_problem(
 *               TimeStepBase_Tria<dim>::RefinementData (
 *                 top_threshold, bottom_threshold)));
 *           },
 *         TimeDependent::TimeSteppingData (0,1),
 *         TimeDependent::forward);
 * @endcode
 * TimeStepBase_Wave<dim>::refine_grid
 * 是一个带参数的函数，与上面在循环中使用的所有其他函数不同。然而，在这种特殊情况下，参数对所有时间段都是一样的，而且在循环开始之前就已经知道了，所以我们把它固定下来，做成一个对外界来说不需要参数的函数对象。
 * 因为它是这个类的中心函数，所以我们最后展示了一个精简版的
 * @p do_loop
 * 方法，展示它的目的是为了让大家更好地了解这个类的内部结构。为了简洁起见，我们省略了处理向后运行循环的部分，以及检查唤醒和睡眠操作是否在<tt>0...n_timesteps-1</tt>以外的时间步长上进行。
 *
 * @code
 * template <typename InitFunctionObject, typename LoopFunctionObject>
 * void TimeDependent::do_loop (InitFunctionObject      init_function,
 *                         LoopFunctionObject      loop_function,
 *                         const TimeSteppingData &timestepping_data,
 *                         const Direction         direction)
 * {
 *   // initialize the time steps for a round of this loop
 *   for (unsigned int step=0; step<n_timesteps; ++step)
 *     init_function (static_cast<typename InitFunctionObject::argument_type>
 *                      (timesteps[step]));
 *
 *   // wake up the first few time levels
 *   for (int step=-timestepping_data.look_ahead; step<0; ++step)
 *     for (int look_ahead=0;
 *          look_ahead<=timestepping_data.look_ahead;
 *          ++look_ahead)
 *       timesteps[step+look_ahead]->wake_up(look_ahead);
 *
 *
 *
 *
 *   for (unsigned int step=0; step<n_timesteps; ++step)
 *     {
 *       // first thing: wake up the timesteps ahead as necessary
 *       for (unsigned int look_ahead=0;
 *            look_ahead<=timestepping_data.look_ahead;
 *            ++look_ahead)
 *         timesteps[step+look_ahead]->wake_up(look_ahead);
 *
 *
 *
 *
 *       // actually do the work
 *       loop_function(
 *         static_cast<typename LoopFunctionObject::argument_type> (
 *           timesteps[step]));
 *
 *       // let the timesteps behind sleep
 *       for (unsigned int look_back=0;
 *            look_back<=timestepping_data.look_back;
 *            ++look_back)
 *         timesteps[step-look_back]->sleep(look_back);
 *     }
 *
 *   // make the last few timesteps sleep
 *   for (int step=n_timesteps;
 *        step<=n_timesteps+timestepping_data.look_back;
 *        ++step)
 *     for (int look_back=0;
 *          look_back<=timestepping_data.look_back;
 *          ++look_back)
 *       timesteps[step-look_back]->sleep(look_back);
 * }
 * @endcode
 *
 *
 */
class TimeDependent
{
public:
  /**
   * 持有两个基本实体的结构，它们控制着所有时间步数的循环：我们应在当前步数的前面多少个时间步数开始唤醒时间步数对象，以及在后面多少个时间步数调用它们的
   * @p sleep 方法。
   *
   */
  struct TimeSteppingData
  {
    /**
     * 构造函数；参见不同的字段，以了解参数的含义。
     *
     */
    TimeSteppingData(const unsigned int look_ahead,
                     const unsigned int look_back);

    /**
     * 这表示分时算法需要向前看的时间步数。通常情况下，这个数字将是0，因为提前看的算法不能作为时隙方案，因为它们不能只从过去的知识中计算它们的数据，因此在时间上是全局的。
     * 然而，在其他情况下，当不想访问下一个时间步骤的数据时，可能有必要向前看，但例如要知道下一个网格，下一个时间层次上的二元问题的解决方案，等等。
     * 请注意，对于一个向后走的二元问题，"向前看
     * "意味着向更小的时间值看。
     * 这个数字的值决定了，时间步长管理器开始在每个时间步长中调用
     * @p wake_up 函数。
     *
     */
    const unsigned int look_ahead;

    /**
     * 这是一个与上述变量相反的变量。它表示在目前的时间步数之后的时间步数，我们需要保留所有的数据，以便在目前的时间水平上进行计算。
     * 对于单步方案（例如欧拉方案或克拉克-尼科尔森方案），这个值是1。
     * 这个数字的值决定了在一个时间层上做完计算后，时间步长管理器将在每个时间步长上调用
     * @p 睡眠函数。
     *
     */
    const unsigned int look_back;
  };

  /**
   * 枚举提供了由 @p
   * do_loop执行的循环可能运行的不同方向。
   *
   */
  enum Direction
  {
    /**
     * 往前走。
     *
     */
    forward,
    /**
     * 在后退方向上走。
     *
     */
    backward
  };

  /**
   * 构造函数。
   *
   */
  TimeDependent(const TimeSteppingData &data_primal,
                const TimeSteppingData &data_dual,
                const TimeSteppingData &data_postprocess);


  /**
   * 销毁器。这将删除给<tt>insert_*</tt>和 @p add_timestep
   * 函数的指针所指向的对象，也就是说，它将删除在每个时间步长进行计算的对象。
   *
   */
  virtual ~TimeDependent();

  /**
   * 在任何位置添加一个时间步长。该位置是一个指向现有时间步长对象的指针，或者一个表示时间步长序列结束的空指针。如果
   * @p position
   * 为非空，新的时间步长将被插入到相应的元素之前。
   * 请注意，通过将一个对象交给这个函数，TimeDependent对象承担了该对象的所有权；因此它也将负责删除其管理的对象。
   * 还有一个函数， @p add_timestep,
   * ，在列表的末尾插入一个时间步骤。
   * 注意，这个函数不会改变存储在其他时间步长对象中的时间步长，也不会设置这个新时间步长的时间步长。这只有在调用
   * @p start_sweep
   * 函数时才会完成。在不改变时间步数的情况下，对时空三角的操作更简单，因为可以一直使用上一次扫频中使用的时间步数。
   *
   */
  void
  insert_timestep(const TimeStepBase *position, TimeStepBase *new_timestep);

  /**
   * 就像 @p insert_timestep, 一样，但在最后插入。
   * 这种机制通常会导致像这样的设置循环
   * @code
   * for (int i=0; i<N; ++i)
   * manager.add_timestep(new MyTimeStep());
   * @endcode
   *
   *
   */
  void
  add_timestep(TimeStepBase *new_timestep);

  /**
   * 删除一个时间段。只有当你想在两次扫描之间删除它时，才有必要调用这个功能；在这个对象的生命周期结束时，会自动注意删除时间步长对象。解构器对对象的删除也是通过这个函数完成的。
   * 注意，这个函数不会改变存储在其他时间步数对象中的时间步数。这只有在调用
   * @p
   * start_sweep函数时才会完成。在不改变时间步数的情况下，对时空三角的操作更加简单，因为可以一直使用上一次扫频中使用的时间步数。
   *
   */
  void
  delete_timestep(const unsigned int position);

  /**
   * 解决原始问题；通过TimeStepBase类的 @p init_for_primal_problem
   * 和 @p solve_primal_problem do_loop函数使用该类的函数。
   * 向前看和向后看是由给构造函数的 @p
   * timestepping_data_primal对象决定的。
   *
   */
  void
  solve_primal_problem();

  /**
   * 解决对偶问题；通过TimeStepBase类的 @p do_loop
   * 函数使用该类的 @p init_for_dual_problem 和 @p solve_dual_problem
   * 函数。    向前看和向后看是由给构造函数的 @p
   * timestepping_data_dual 对象决定的。
   *
   */
  void
  solve_dual_problem();

  /**
   * 做一轮后处理；通过本类的 @p do_loop
   * 函数使用TimeStepBase类的 @p init_for_postprocessing 和 @p
   * postprocess 函数。    向前看和向后看是由给构造函数的 @p
   * timestepping_data_postprocess对象决定的。
   *
   */
  void
  postprocess();

  /**
   * 在所有时间步长上做一个循环，在开始时调用 @p
   * init_function ，在每个时间步长的 @p loop_function 上调用。
   * @p timestepping_data
   * 决定在当前时间段前后有多少个时间段调用 @p  wake_up和
   * @p sleep 函数。
   * 为了了解这个函数的工作原理，请注意，函数 @p
   * solve_primal_problem只包括以下调用。
   * @code
   * do_loop([](TimeStepBaseconst time_step)
   *         { time_step->init_for_primal_problem(); },
   *       [](TimeStepBaseconst time_step)
   *         { time_step->solve_primal_problem(); },
   *       timestepping_data_primal,
   *       forward);
   * @endcode
   * 还请注意，这两个函数所来自的给定类不一定是TimeStepBase，但也可以是一个派生类，即
   * @p static_castable
   * 来自一个TimeStepBase。该函数可以是该类的一个虚函数（甚至是一个纯函数），如果实现该函数的实际类是一个通过虚基类派生的类，从而无法通过
   * @p static_cast 从TimeStepBase类到达，这应该会有帮助。
   * 不使用上述形式，你同样可以使用<tt>[args...](Xconst
   * x){x->unary_function(args...);}</tt>，让 @p do_loop
   * 函数用指定参数调用给定的函数。
   *
   */
  template <typename InitFunctionObject, typename LoopFunctionObject>
  void
  do_loop(InitFunctionObject      init_function,
          LoopFunctionObject      loop_function,
          const TimeSteppingData &timestepping_data,
          const Direction         direction);


  /**
   * 为下一次扫描初始化对象。这个函数具体做了以下工作：给每个时间层分配它目前在数组中的编号（如果时间层被插入或删除，这个编号可能会改变），并将本次扫描的编号传送给这些对象。
   * 在上述数字设定后，它还调用每个时间层对象的 @p
   * start_sweep 函数。
   * 这个函数是虚拟的，所以你可以重载它。然而，你不应该忘记在你的重载版本中也调用这个函数，最好是在你的函数的开头，因为这是某种
   * "类似构造器 "的函数，应该自下而上地调用。
   * 这个函数的默认实现在所有时间步长对象上调用 @p
   * start_sweep 。
   *
   */
  virtual void
  start_sweep(const unsigned int sweep_no);

  /**
   * 与上述函数类似，调用每个时间步长对象的 @p end_sweep
   * 。这个函数的 @p virtualness 与前一个函数的 @p virtualness
   * 同样适用。
   * @note  这个函数并不保证 @p end_sweep
   * 对连续的时间步长连续调用，相反，调用该函数的时间步长对象的顺序是任意的。因此，你不应该认为该函数已经为以前的时间步数调用过了。如果在多线程模式下，几个时间步骤的
   * @p end_sweep
   * 函数可能会被同时调用，所以如果你的程序需要，你应该使用同步机制。
   *
   */
  virtual void
  end_sweep();

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 异常情况。
   *
   */
  DeclExceptionMsg(ExcInvalidPosition,
                   "You cannot insert a time step at the specified position.");

protected:
  /**
   * 持有指向时间级别对象的指针的向量。这是此对象操作的主要数据。请注意，这个对象占有这个集合中的指针所指向的对象。
   *
   */
  std::vector<SmartPointer<TimeStepBase, TimeDependent>> timesteps;

  /**
   * 当前扫频的编号。这将由每次扫描开始时调用的 @p
   * start_sweep 函数重置。
   *
   */
  unsigned int sweep_no;

  /**
   * 一些告诉 @p solve_primal_problem
   * 函数要做什么的标志。更多信息请参见该结构的文档。
   *
   */
  const TimeSteppingData timestepping_data_primal;

  /**
   * 一些告诉 @p solve_dual_problem
   * 函数要做什么的标志。更多信息请参见此结构的文档。
   *
   */
  const TimeSteppingData timestepping_data_dual;

  /**
   * 一些告诉 @p postprocess
   * 函数要做什么的标志。更多信息请参见此结构的文档。
   *
   */
  const TimeSteppingData timestepping_data_postprocess;

private:
  /**
   * 只对某些时间段做<tt>end_sweep()/tt>的工作。这在多线程模式下是很有用的。
   *
   */
  void
  end_sweep(const unsigned int begin_timestep, const unsigned int end_timestep);
};



/**
 * 在时间相关问题中，一个时间步长的基类。这个类只提供了基本的框架，定义了必要的虚拟函数（即
 * @p sleep 和 @p wake_up),
 * 前一个和后一个网格的接口，以及一些在所有时间步长的新循环开始前要调用的函数。
 *
 *
 */
class TimeStepBase : public Subscriptor
{
public:
  /**
   * 表示接下来要解决的问题类型的枚举。
   *
   */
  enum SolutionState
  {
    /**
     * 接下来解决原始问题。
     *
     */
    primal_problem = 0x0,
    /**
     * 接下来解决对偶问题。
     *
     */
    dual_problem = 0x1,
    /**
     * 接下来进行后处理。
     *
     */
    postprocess = 0x2
  };

  /**
   * 构造函数。除了设置时间外，这里不做任何事情。
   *
   */
  TimeStepBase(const double time);

  /**
   * 解构器。目前，它什么都不做。
   *
   */
  virtual ~TimeStepBase() override = default;

  /**
   * 删除了拷贝构造函数，以避免浅层拷贝的意外行为。
   *
   */
  TimeStepBase(const TimeStepBase &) = delete;

  /**
   * 删除了复制赋值操作符，以避免浅层拷贝的意外行为。
   *
   */
  TimeStepBase &
  operator=(const TimeStepBase &) = delete;

  /**
   * 重建这个时间级别工作所需的所有数据。这个函数的作用是在所有的变量和数据结构被送入睡眠状态后，或者在我们第一次访问这个时间层时，让它们重新开始工作。特别是，它被用来重建三角形、自由度处理程序，在数据向量被存储到磁盘的情况下重新加载它们，等等。
   * 这个函数的实际实现没有任何作用。
   * 由于这是一个重要的任务，如果你选择在自己的类中重载这个函数（很可能是这样），你应该从自己的函数中调用这个函数，最好是在开始的时候，这样你的函数就可以对已经存在的三角测量产生影响。
   *
   */
  virtual void
  wake_up(const unsigned int);

  /**
   * 这是一个与 @p wake_up.
   * 相反的函数，它用于删除数据或在当前扫描不再需要它们时将其保存到磁盘。
   * 典型的数据种类是数据向量、自由度处理程序、三角测量对象等，它们占据了大量的内存，因此可能被外部化。
   * 默认情况下，这个函数不做任何事情。
   *
   */
  virtual void
  sleep(const unsigned int);

  /**
   * 每次在开始新的扫频之前都会调用这个函数。你可能想在计算过程中设置一些需要的字段，等等。然而，你应该很小心，不要安装大的对象，这应该推迟到
   * @p wake_up 函数被调用。
   * 这个函数的一个典型动作是整理解题过程中需要的临时文件的名称，等等。
   * 在这个函数被调用时， @p timestep_no,   @p
   * sweep_no的值以及指向上一个和下一个时间步长对象的指针已经具有正确的值。
   * 这个函数的默认实现不做任何事情。
   *
   */
  virtual void
  start_sweep();

  /**
   * 这个函数类似于上面的函数，但它是在扫频结束时调用的。你通常想在这个函数中进行清理，比如删除临时文件之类的。
   *
   */
  virtual void
  end_sweep();

  /**
   * 在每个时间层的原始问题被解决之前，这个函数被调用（也就是在第一个时间层的解决发生之前）。默认情况下，该函数设置该类的
   * @p next_action 变量。
   * 你可以重载这个函数，但你应该在你自己的函数中调用这个函数。
   *
   */
  virtual void
  init_for_primal_problem();

  /**
   * 和上面一样，但在一轮对偶问题解决之前调用。
   *
   */
  virtual void
  init_for_dual_problem();

  /**
   * 同上，但在一轮后处理步骤之前调用。
   *
   */
  virtual void
  init_for_postprocessing();

  /**
   * 当需要解决这个时间层次上的原始问题时，这个函数被管理器对象调用。它在
   * @p wake_up 函数被调用之后，在 @p sleep
   * 函数将被调用之前被调用。
   * 由于明显的原因，没有默认的实现，所以你必须重载这个函数。
   *
   */
  virtual void
  solve_primal_problem() = 0;

  /**
   * 当需要解决这个时间层次上的二元问题时，这个函数会被管理器对象调用。它在
   * @p wake_up 函数被调用后，在 @p sleep
   * 函数将被调用前被调用。
   * 有一个默认的实现是什么都不做，因为有些问题可能不需要解决对偶问题。然而，当被调用时，它将中止程序，因为那时你真的应该重载该函数。
   *
   */
  virtual void
  solve_dual_problem();

  /**
   * 当需要对这个时间级别进行后处理时，这个函数被管理器对象调用。它在
   * @p wake_up 函数被调用后，在 @p sleep
   * 函数将被调用前被调用。有一个默认的实现是什么都不做，因为有些问题可能不需要做后处理步骤，例如，如果在解决原始问题时已经完成了一切。然而，当被调用时，它将中止程序，因为那时你真的应该重载这个函数。
   *
   */
  virtual void
  postprocess_timestep();

  /**
   * 返回这个时间步骤的时间值。
   *
   */
  double
  get_time() const;

  /**
   * 返回这个时间步长的数字。注意，如果增加或删除时间步长，这个数字在不同的扫频之间可能会有所不同。
   *
   */
  unsigned int
  get_timestep_no() const;

  /**
   * 计算与上一个时间步长的时间差。如果这个时间步长是第一个时间步长，这个函数将导致一个异常。虽然这种行为看起来有点激烈，但在大多数情况下是合适的，因为如果没有上一个时间步长，无论如何你都需要特殊处理，这样就不会返回无效的值，这可能导致错误的但未被注意的计算结果。(在这种情况下，唯一合理的返回值不会是零，因为有效的计算可以用它来完成，而是一个变性的值，如
   * @p NaN.
   * 。然而，发现计算结果都是变性的值或得到一个异常并没有太大的区别；在后一种情况下，你至少可以得到你问题所在的确切位置。)
   *
   */
  double
  get_backward_timestep() const;

  /**
   * 返回到下一个时间步骤的时间差。关于没有下一个时间步骤的情况，和上面的函数同样适用。
   *
   */
  double
  get_forward_timestep() const;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   * 你要在派生类中重载这个函数来计算派生类所使用的内存量，并把这个函数的结果加到你的结果中。
   *
   */
  virtual std::size_t
  memory_consumption() const;

protected:
  /**
   * 指向列表中前一个时间步骤对象的指针。
   *
   */
  const TimeStepBase *previous_timestep;

  /**
   * 指向列表中下一个时间步长对象的指针。
   *
   */
  const TimeStepBase *next_timestep;

  /**
   * 我们目前所处的扫频的编号。这个数字在扫频开始前由时间级别管理器重置。
   *
   */
  unsigned int sweep_no;

  /**
   * 时间步骤的编号，从零开始计算。这个数字在每次扫频开始时由时间级别管理器重置，因为一些时间步长可能在上一次扫频后被插入或删除。
   *
   */
  unsigned int timestep_no;

  /**
   * 本级操作的离散时间。
   *
   */
  const double time;

  /**
   * 存储原始问题或对偶问题的解决是否实际的变量，或任何其他指定的行动。这个变量是由<tt>init_for_*</tt>函数设置的。
   *
   */
  unsigned int next_action;

private:
  /**
   * 重置指向上一个时间步长的指针；只应由时间层管理器对象调用。
   * 该函数在管理器对象的设置以及每当插入或删除一个时间步长时被调用。
   *
   */
  void
  set_previous_timestep(const TimeStepBase *previous);

  /**
   * 重置指向下一个时间步长的指针；只应由时间级管理器对象调用。
   * 该函数在管理器对象的设置以及每当插入或删除一个时间步长时被调用。
   *
   */
  void
  set_next_timestep(const TimeStepBase *next);

  /**
   * 设置该时间步长在时间步长列表中的编号。这个函数在每次扫描开始时由时间步长管理对象调用，以更新由于增加或删除时间层而可能改变的信息。
   *
   */
  void
  set_timestep_no(const unsigned int step_no);

  /**
   * 设置我们目前所处的扫频的编号。这个函数是由时间层管理对象在每次扫频的启动时调用的。
   *
   */
  void
  set_sweep_no(const unsigned int sweep_no);

  // make the manager object a friend
  friend class TimeDependent;
};



/**
 * 名称空间，其中声明了一些类，这些类封装了TimeStepBase_Tria()类的标志。这些曾经是该类的本地数据类型，但是一些编译器在某些方面扼杀了它们，所以我们把它们放到了一个自己的命名空间。
 *
 *
 */
namespace TimeStepBase_Tria_Flags
{
  /**
   * 这个结构是用来告诉TimeStepBase_Tria()类应该如何处理网格。它有一些标志，定义了网格应被重新制作以及何时可以删除的时刻。此外，还有一个变量规定了网格是否应该保留在内存中，或者应该在两次使用之间删除以节省内存。
   *
   */
  template <int dim>
  struct Flags
  {
    /**
     * 默认的构造函数；产生一个异常，所以不能真正使用。
     *
     */
    Flags();

    /**
     * 构造函数；请看不同字段对参数含义的描述。
     *
     */
    Flags(const bool         delete_and_rebuild_tria,
          const unsigned int wakeup_level_to_build_grid,
          const unsigned int sleep_level_to_delete_grid);

    /**
     * 这个标志决定了 @p sleep 和 @p wake_up
     * 函数是否要删除和重建三角结构。
     * 虽然对于小问题来说，这不是必要的，但对于大问题来说，这对于节省内存是必不可少的。
     * 原因是内存中可能有几百个时间层，每个时间层都有自己的三角图，如果每个时间层有很多单元，可能需要大量的内存。在所有时间层上总共有100.000.000个单元的情况并不罕见，这使得这个标志可以理解。
     *
     */
    const bool delete_and_rebuild_tria;

    /**
     * 这个数字表示 @p wake_up
     * 函数的参数，它应该在这个参数上重建网格。显然，它应小于或等于传递给时间步长管理对象的
     * @p look_ahead
     * 数字；如果它相等，那么网格将在第一次调用 @p wake_up
     * 函数时重建。如果 @p delete_and_rebuild_tria 是 @p false,
     * 这个数字没有意义。
     *
     */
    const unsigned int wakeup_level_to_build_grid;

    /**
     * 这是一个与上面相反的标志：它决定了在调用 @p sleep
     * 时，网格应被删除。
     *
     */
    const unsigned int sleep_level_to_delete_grid;
  };



  /**
   * 这个结构是用来告诉TimeStepBase_Tria()类应该如何完善网格的。在我们解释所有不同的变量之前，先握紧一些术语。    <ul>   <li> 修正：在对三角形的一些单元进行了遵循某种给定标准的标记后，我们可能想根据另一种标准改变这个网格上的标记单元的数量，这个单元的数量可能只比之前的网格上的单元数量多或少一定的分数。这种细化标志的改变在后面将被称为 "修正"。    <li>  适应：为了使一个网格和下一个网格之间的变化不大，我们可能想在两个网格中的一个上标记一些额外的单元格，这样就不会有太严重的差异。这个过程将被称为 "适应"。    </ul>  <h3>Description of flags</h3>  <ul>   <li>   @p max_refinement_level:  切断某一级别的单元的细化。这个标志并不影响单元的标记，所以在较粗的层次上并没有比平时更多的单元被标记。相反，标志都被设置，但是当涉及到实际细化时，最大细化水平被截断。    这个选项只有在你想比较全局细化和自适应细化时才真正有用，因为你不希望后者的细化程度超过全局细化。      <li>   @p first_sweep_with_correction:  当使用上面定义的单元数校正时，可能值得只在以后的扫描中开始使用，而不是已经在第一次扫描中开始。如果这个变量是零，那么就从第一次扫描开始，否则就从更高的扫描开始。只从后面开始的理由是，我们不想在一开始就阻止网格的发展，而只在我们开始对计算的实际结果感兴趣的扫描中施加限制。      <li>   @p min_cells_for_correction:  如果我们想要一个更自由的网格发展过程，我们可能希望对单元数少的网格也施加较少的规则。这个变量为要进行修正的网格的单元数设定了一个下限。      <li>   @p cell_number_corridor_top:  一个网格的单元数可能高于前一个网格的单元数的分数。常见值为10%（即0.1）。该变量的命名源于为细化后的单元格数量定义一个目标走廊的目标。      <li>   @p cell_number_corridor_bottom:  一个网格的单元数可能低于前一个网格的单元数的分数。常见值为5%（即0.05）。通常这个数字会小于 @p cell_number_corridor_top ，因为单元格数量的增加是无害的（尽管它增加了解决问题所需的数字工作量），而急剧减少可能会降低最终结果的精度，即使在减少之前计算的时间步数是以高精确度计算的。    但是请注意，如果你也计算对偶问题，那么时间方向是相反的，所以定义单元数走廊的两个值应该是差不多的。      <li>   @p correction_relaxations:  这是一个数字对的列表，其含义如下：正如 @p min_cells_for_correction, 一样，如果网格的单元数很少，可能值得减少对网格的要求。  本变量存储了一个单元格编号的列表以及一些数值，这些数值告诉我们，单元格编号的走廊应该被放大一定的系数。例如，如果这个列表是<tt>((100 5) (200 3) (500 2))</tt>，这将意味着对于单元数低于100的网格，<tt>cell_number_corridor_*</tt>变量在应用前要乘以5，对于单元数低于200的网格，要乘以3，以此类推。      @p correction_relaxations 实际上是这种列表的一个向量。这个向量中的每个条目都表示一个扫描的放松规则。最后一个条目定义了所有后续扫描的松弛规则。采用这种方案是为了允许在以后的扫描中进行更严格的修正，而在第一次扫描中的放松可能更宽松。    有一个静态变量 @p default_correction_relaxations ，你可以把它作为一个默认值。它是一个空列表，因此没有定义放宽。      <li>   @p cell_number_correction_steps:  通常情况下，如果你想修正单元格的数量，会计算单元格数量的目标走廊，并对一些额外的单元格进行标记或去除标记。但由于标记和除旗后产生的单元格数不容易计算，通常不会在走廊内。因此，我们需要迭代来达到我们的目标。通常需要三到四次迭代，但使用这个变量，你可以减少允许的迭代次数；在两次迭代后打破循环，经常会产生好的结果。将该变量设置为零将导致完全没有修正步骤。      <li>   @p mirror_flags_to_previous_grid:  如果当前网格上的一个单元被标记为细化，也要标记前一个网格上的相应单元。例如，如果计算时空单元的误差指标，但只存储在第二个网格中，这就很有用。现在，由于第一个网格对指标的贡献与第二个网格相同，如果有必要，对这两个网格进行标记可能是有用的。如果设置了本变量就可以做到这一点。      <li>   @p adapt_grids:  在上面定义的意义上，使现在的网格与以前的网格相适应。这里实际做的是：如果从以前的网格到现在的网格会导致某些单元的双倍细化或双倍粗化，那么我们就尝试标记这些单元进行细化或粗化，从而避免双倍步骤。  显然，超过两次的细化和粗化也会被发现。    网格自适应可以尝试避免两个网格之间的这种变化，但它永远不能保证不发生这种变化。这是因为下一个网格可能会改变现在的网格，但同样地，现在的网格和上一个网格之间的细化水平可能会有跳跃；这只能通过反复循环所有的网格来避免，来回循环，直到没有任何变化为止，如果有很多时间步长的非常大的网格，这显然是不可能的。    </ul>
   *
   */
  template <int dim>
  struct RefinementFlags
  {
    /**
     * 描述校正过程的一些放松的数据类型的类型定义。更多信息见该类的一般描述。
     *
     */
    using CorrectionRelaxations =
      std::vector<std::vector<std::pair<unsigned int, double>>>;

    /**
     * 松弛的默认值：没有松弛。
     *
     */
    static CorrectionRelaxations default_correction_relaxations;

    /**
     * 构造函数。默认值的选择是对网格细化几乎没有限制的。
     *
     */
    RefinementFlags(const unsigned int max_refinement_level        = 0,
                    const unsigned int first_sweep_with_correction = 0,
                    const unsigned int min_cells_for_correction    = 0,
                    const double       cell_number_corridor_top    = (1 << dim),
                    const double       cell_number_corridor_bottom = 1,
                    const CorrectionRelaxations &correction_relaxations =
                      CorrectionRelaxations(),
                    const unsigned int cell_number_correction_steps  = 0,
                    const bool         mirror_flags_to_previous_grid = false,
                    const bool         adapt_grids                   = false);

    /**
     * 一个时间层次的三角形中的单元的最大层次。如果它被设置为零，那么对一个粗网格单元可能进行的细化次数没有限制。通常情况下，如果出于某种原因，你想在自适应过程中限制细化，例如避免单元数量过大，或者与有一定细化次数的网格进行比较，就会使用这个字段。
     *
     */
    const unsigned int max_refinement_level;

    /**
     * 第一次扫面执行单元格数量校正步骤；对于之前的扫面，单元格只被标记，不执行对之前网格的数量校正。
     *
     */
    const unsigned int first_sweep_with_correction;


    /**
     * 仅当有超过此数量的细胞时，才应用前一时间水平的细胞数校正。
     *
     */
    const unsigned int min_cells_for_correction;

    /**
     * 一个时间层上的细胞数可能与前一个时间层上的细胞数不同的分数（第一：上偏差，第二：下偏差）。
     *
     */
    const double cell_number_corridor_top;

    /**
     * @ref cell_number_corridor_top
     *
     */
    const double cell_number_corridor_bottom;

    /**
     * 对校正步骤的放松列表。
     *
     */
    const std::vector<std::vector<std::pair<unsigned int, double>>>
      correction_relaxations;

    /**
     * 为了将某一时间层的单元格数量调整到上一时间层的单元格数量，要进行的迭代次数。零意味着：不做这样的迭代。
     *
     */
    const unsigned int cell_number_correction_steps;

    /**
     * 标记所有在这个时间段被标记的单元格，以便在上一个时间段也进行细化。这在误差指标是通过对时间空间单元的积分计算出来的，但现在与离散时间层面上的网格相关联的情况下很有用。然而，由于误差贡献来自于两个网格，因此对两个网格进行细化是合适的。
     * 由于前一个网格并没有将标志映照到前一个网格上，所以这并不会导致单元数几乎无限增长。然而，你应该只在开启单元格号校正的情况下使用这个标志。
     * 镜像是在细胞数校正完成后，但在网格适应前完成的，所以这个网格上的细胞数不会明显受到前一个网格上额外标记的细胞的影响。
     *
     */
    const bool mirror_flags_to_previous_grid;

    /**
     * 使这个网格与前一个网格相适应。
     *
     */
    const bool adapt_grids;

    /**
     * 异常情况
     *
     */
    DeclException1(ExcInvalidValue,
                   double,
                   << "The value " << arg1
                   << " for the cell number corridor does not fulfill "
                      "its natural requirements.");
  };



  /**
   * 给予实际细化函数的结构，告诉它在粗化和细化时要采取哪些阈值。实际的细化标准是通过调用虚拟函数
   * @p  get_tria_refinement_criteria加载的。
   *
   */
  template <int dim>
  struct RefinementData
  {
    /**
     * 构造函数
     *
     */
    RefinementData(const double refinement_threshold,
                   const double coarsening_threshold = 0);

    /**
     * 细化的阈值：具有较大数值的单元格将被细化（至少在第一轮；细化过程的后续步骤可能也会标记其他单元格，或者从具有高于该阈值的标准的单元格中移除标记）。
     *
     */
    const double refinement_threshold;

    /**
     * 粗化的相同阈值：如果可能的话，阈值较小的单元将被粗化。
     *
     */
    const double coarsening_threshold;

    /**
     * 异常情况
     *
     */
    DeclException1(ExcInvalidValue,
                   double,
                   << "The value " << arg1
                   << " for the cell refinement thresholds does not fulfill "
                      "its natural requirements.");
  };
} // namespace TimeStepBase_Tria_Flags



/**
 * TimeStepBase的特殊化，解决了网格处理的某些方面。特别是，这个类被认为可以处理网格，这些网格可以在每个时间步长上分别进行自适应的细化，或者在时间步长之间进行松散耦合。通过基类中声明的
 * @p sleep 和 @p wake_up
 * 函数，它还负责在内存资源是一个点时删除和重建网格。
 * 除此之外，它还提供了一些函数，对与时间有关的问题做了一些相当棘手的细化规则，试图避免在随后的时间层之间的网格有太多的变化，同时也试图保留单独细化每个网格的自由。有很多控制这个函数的标志和数字，它们可能会极大地改变这个函数的行为
 *
 * - 进一步的信息见标志的描述。
 *
 *
 */
template <int dim>
class TimeStepBase_Tria : public TimeStepBase
{
public:
  /**
   * 将TimeStepBase_Tria_Flags()命名空间的数据类型类型化为本地范围。
   *
   */
  using Flags = typename TimeStepBase_Tria_Flags::Flags<dim>;
  using RefinementFlags =
    typename TimeStepBase_Tria_Flags::RefinementFlags<dim>;
  using RefinementData = typename TimeStepBase_Tria_Flags::RefinementData<dim>;


  /**
   * 扩展基类中的枚举，表示下一个要做的动作。
   *
   */
  enum SolutionState
  {
    /**
     * 接下来执行网格细化。
     *
     */
    grid_refinement = 0x1000
  };


  /**
   * 默认构造函数。什么都不做，只是抛出一个异常。我们需要有这样一个构造函数，以满足派生类的需要，它们把这个类作为一个虚拟基类，不需要调用这个构造函数，因为它们不是终端类。编译器希望知道一个可以调用的构造函数，因为它不能知道这个类不是终端类。
   *
   */
  TimeStepBase_Tria();

  /**
   * 构造函数。接受一个粗略的网格，这个时间层的网格将从这个粗略的网格派生出来，还有一些引导这个对象的行为的标志。
   * 粗略网格的所有权保留在这个对象的创建者那里。
   * 然而，它被锁定，不会被破坏，以保证粗略网格的寿命比这个对象需要的时间长。
   * 你需要给这个函数提供一般的标志结构，因为它无论如何都需要；如果你不打算调用这个类的细化函数，细化标志可以省略。
   *
   */
  TimeStepBase_Tria(
    const double                   time,
    const Triangulation<dim, dim> &coarse_grid,
    const Flags &                  flags,
    const RefinementFlags &        refinement_flags = RefinementFlags());

  /**
   * 解构器。目前，这不过是释放给构造函数的粗略网格三角的锁。
   *
   */
  virtual ~TimeStepBase_Tria() override;

  /**
   * 重建这个时间层所需的所有数据。这个函数的作用是在所有的变量和数据结构被送入睡眠状态后，或者在我们第一次访问这个时间层时，让它们重新开始工作。特别是，它被用来重建三角形、自由度处理程序，在数据向量被存储到磁盘的情况下重新加载它们，等等。默认情况下，如果在
   * @p sleep
   * 函数中设置了相应的标志来破坏三角形，该函数将重建三角形。它也会在我们第一次碰到这个函数并且
   * @p wakeup_level
   * 等于<tt>flags.wakeup_level_to_build_grid</tt>时这样做，与上述标志的值无关。(实际上，只要三角形指针等于Null指针并且
   * @p wakeup_level 的值是正确的，它就会这样做)。
   * 由于这是一项重要的任务，如果你选择在自己的类中重载这个函数（很可能是这样），你应该从自己的函数中调用这个函数，最好是在开始的时候，这样你的函数就可以对已经存在的三角关系生效。
   *
   */
  virtual void
  wake_up(const unsigned int wakeup_level) override;

  /**
   * 这是一个与 @p wake_up.
   * 相反的函数，它用于删除数据或在当前扫描不再需要它们时将其保存到磁盘。
   * 典型的数据种类是数据向量、自由度处理程序、三角测量对象等，它们占据了大量的内存，因此可以被外部化。
   * 默认情况下，如果用户在此对象的标志中这样指定，三角测量将被删除，细化历史将被保存，以便相应的
   * @p wake_up
   * 函数可以重建它。因此，你应该从你的重载版本中调用这个函数，最好是在最后，这样你的函数就可以在你需要的时候使用三角测量。
   *
   */
  virtual void
  sleep(const unsigned int) override;

  /**
   * 根据传递给这个对象的构造函数的标志和传递给这个函数的数据进行细化。关于这个函数的工作描述，请参考这个类的一般文档。
   * 事实上，这个函数实际上并没有细化或粗化三角形，而只是设置各自的标志。这样做的原因是，通常你不会在事后立即需要这个网格，而是在下一次扫描时才需要，所以只需存储这些标志，并在下次需要时重建它。另外，可能下一步会想增加或删除一些标志，所以我们无论如何都要使用这个网格来等待。
   *
   */
  void
  refine_grid(const RefinementData data);

  /**
   * 细化循环的相应初始函数；除了将 @p next_action 设置为 @p
   * grid_refinement外，在默认实现中不做任何事情，但可以被重载。
   *
   */
  virtual void
  init_for_refinement();

  /**
   * 虚拟函数，应该用当前三角形的细化标准填充向量。它在
   * @p refine_grid
   * 函数内使用，以获得当前时间步长的标准，因为在使用时间步长管理对象的循环时，它们不能通过其参数传递。
   *
   */
  virtual void
  get_tria_refinement_criteria(Vector<float> &criteria) const = 0;

  /**
   * 三角形的细化标志被存储在一个局部变量中，因此可以进行恢复。粗化标志也被存储。
   *
   */
  void
  save_refine_flags();

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   * 你要在派生类中重载这个函数来计算派生类所使用的内存量，并将这个函数的结果加到你的结果中。
   *
   */
  virtual std::size_t
  memory_consumption() const override;

  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcGridNotDeleted,
                   "When calling restore_grid(), you must have previously "
                   "deleted the triangulation.");

protected:
  /**
   * 在这个时间层次上使用的三角法。由于这是每一个时间步进方案都需要有的东西，我们可以安全地把它放到基类中。请注意，如果在
   * @p flags 结构中指定了这样的行为，三角形会经常被 @p
   * sleep 和 @p wake_up 函数删除和重建以节省内存。
   *
   */
  SmartPointer<Triangulation<dim, dim>, TimeStepBase_Tria<dim>> tria;

  /**
   * 指向一个网格的指针，该网格将被用作该时间层的粗略网格。
   * 这个指针是通过构造函数设置的；所有权仍然属于这个管理对象的所有者。
   *
   */
  SmartPointer<const Triangulation<dim, dim>, TimeStepBase_Tria<dim>>
    coarse_grid;

  /**
   * 一些关于这个时间级别应如何表现的标志。参见此结构的文档以了解更多信息。
   *
   */
  const Flags flags;

  /**
   * 控制细化过程的标志；更多信息请参见相应结构的文档。
   *
   */
  const RefinementFlags refinement_flags;

private:
  /**
   * 保存该时间层次上不同扫描的细化和粗化标志的向量。因此，这些向量保存了网格的历史。
   *
   */
  std::vector<std::vector<bool>> refine_flags;

  /**
   * @ref refine_flags
   *
   */
  std::vector<std::vector<bool>> coarsen_flags;

  /**
   * 根据保存的数据恢复网格。为此，粗略的网格被复制，并使用保存的标志逐步重建网格。
   *
   */
  void
  restore_grid();
};



 /*--------------------------- template functions ----------------------------*/ 

template <typename InitFunctionObject, typename LoopFunctionObject>
void
TimeDependent::do_loop(InitFunctionObject      init_function,
                       LoopFunctionObject      loop_function,
                       const TimeSteppingData &timestepping_data,
                       const Direction         direction)
{
  // the following functions looks quite
  // disrupted due to the recurring switches
  // for forward and backward running loops.
  //
  // I chose to switch at every place where
  // it is needed, since it is so easy
  // to overlook something when you change
  // some code at one place when it needs
  // to be changed at a second place, here
  // for the other direction, also.

  const unsigned int n_timesteps = timesteps.size();

  // initialize the time steps for
  // a round of this loop
  for (unsigned int step = 0; step < n_timesteps; ++step)
    switch (direction)
      {
        case forward:
          init_function((&*timesteps[step]));
          break;
        case backward:
          init_function((&*timesteps[n_timesteps - step - 1]));
          break;
      };


  // wake up the first few time levels
  for (int step = -static_cast<int>(timestepping_data.look_ahead); step < 0;
       ++step)
    for (int look_ahead = 0;
         look_ahead <= static_cast<int>(timestepping_data.look_ahead);
         ++look_ahead)
      switch (direction)
        {
          case forward:
            if (step + look_ahead >= 0)
              timesteps[step + look_ahead]->wake_up(look_ahead);
            break;
          case backward:
            if (n_timesteps - (step + look_ahead) < n_timesteps)
              timesteps[n_timesteps - (step + look_ahead)]->wake_up(look_ahead);
            break;
        };


  for (unsigned int step = 0; step < n_timesteps; ++step)
    {
      // first thing: wake up the
      // timesteps ahead as necessary
      for (unsigned int look_ahead = 0;
           look_ahead <= timestepping_data.look_ahead;
           ++look_ahead)
        switch (direction)
          {
            case forward:
              if (step + look_ahead < n_timesteps)
                timesteps[step + look_ahead]->wake_up(look_ahead);
              break;
            case backward:
              if (n_timesteps > (step + look_ahead))
                timesteps[n_timesteps - (step + look_ahead) - 1]->wake_up(
                  look_ahead);
              break;
          };


      // actually do the work
      switch (direction)
        {
          case forward:
            loop_function((&*timesteps[step]));
            break;
          case backward:
            loop_function((&*timesteps[n_timesteps - step - 1]));
            break;
        };

      // let the timesteps behind sleep
      for (unsigned int look_back = 0; look_back <= timestepping_data.look_back;
           ++look_back)
        switch (direction)
          {
            case forward:
              if (step >= look_back)
                timesteps[step - look_back]->sleep(look_back);
              break;
            case backward:
              if (n_timesteps - (step - look_back) <= n_timesteps)
                timesteps[n_timesteps - (step - look_back) - 1]->sleep(
                  look_back);
              break;
          }
    }

  // make the last few timesteps sleep
  for (int step = n_timesteps;
       step < static_cast<int>(n_timesteps + timestepping_data.look_back);
       ++step)
    for (int look_back = 0;
         look_back <= static_cast<int>(timestepping_data.look_back);
         ++look_back)
      switch (direction)
        {
          case forward:
            if ((step - look_back >= 0) &&
                (step - look_back < static_cast<int>(n_timesteps)))
              timesteps[step - look_back]->sleep(look_back);
            break;
          case backward:
            if ((step - look_back >= 0) &&
                (step - look_back < static_cast<int>(n_timesteps)))
              timesteps[n_timesteps - (step - look_back) - 1]->sleep(look_back);
            break;
        };
}

DEAL_II_NAMESPACE_CLOSE

 /*----------------------------   time-dependent.h ---------------------------*/ 
#endif
 /*----------------------------   time-dependent.h ---------------------------*/ 


