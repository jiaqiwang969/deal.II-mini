//include/deal.II-translator/A-headers/multithreading_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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


/**
 *    @defgroup threads Parallel computing with multiple processors accessing
 * shared memory
 * @ingroup Parallel
 * @brief A module discussing the use of parallelism on shared memory  机器。参见详细的文档和 @ref MTToC "目录"
 * ，下面是该模块成员的冗长列表。 @dealiiVideoLecture{39,40}
 * 在具有多个处理器（或多核处理器）的机器上，以%并行方式运行计算的几个部分往往是有利的。例如，我们可以有几个线程在%并行运行，每个线程组装三角形的一个子集的单元格矩阵，然后将它们写入全局矩阵。由于组装矩阵通常是一个昂贵的操作，这经常导致在多处理器机器上显著节省计算时间。
 * 通过Threads命名空间中的函数和类，deal.II支持在共享内存（SMP）机器上以%并行方式运行操作。MultithreadInfo类允许查询系统的某些属性，如CPU的数量。这些用于%并行计算的设施将在下文中描述。
 * step-9 、 step-13 、 step-14 、 step-32 、 step-35 和 step-37
 * 的教程程序也展示了它们在实践中的使用，其中从 step-32
 * 开始的程序采用了更现代的做事风格，基本上我们描述<i>what</i>可以用%parallel完成，而旧的教程程序描述<i>how</i>事必须用%parallel完成。
 * 另一方面，在分布式内存机器（即集群）上运行的程序需要一个建立在MPI和PETSc或Trilinos之上的不同编程模型。这在
 * step-17  ,  step-18  和  step-32  示例程序中有所描述。 @anchor
 * MTToC  <table class="tutorial" width="50%"> <tr><th><b>%Table of
 * contents</b></th></tr> <tr><td width="100%" valign="top">
 * <ol>
 * <li> @ref MTTasks "Task-based parallelism"
 * <li> @ref MTUsing "Using tasks from within deal.II"
 * <li> @ref MTHow "How scheduling tasks works and when task-based programming is not efficient"
 * <li> @ref MTSimpleLoops "Abstractions for tasks: Simple loops"
 * <li> @ref MTComplexLoops "Abstractions for tasks: More complex loops"
 * <li> @ref MTWorkStream "Abstractions for tasks: Work streams"
 * <li> @ref MTTaskSynchronization "Tasks and synchronization"
 * <li> @ref MTThreads "Thread-based parallelism"
 * <li> @ref MTTaskThreads "Controlling the number of threads used for tasks"
 * </ol> </td> </tr> </table>  。
 *   @anchor  MTTasks <h3>Task-based parallelism</h3> 。
 * 在共享内存机器上并行的传统观点是将程序分解成<i>threads</i>，即以%并行的方式运行程序的不同部分<i>at the same time</i>（如果你的机器上的线程多于处理器内核，操作系统会在将执行切换到另一个线程之前短暂地轮流运行每个线程，从而模拟线程并发运行）。下面描述了deal.II的线程设施（见 @ref MTThreads  "基于线程的并行性"
 * ），但我们首先想讨论一个通常比线程更合适的抽象。<i>tasks</i>.
 * 任务本质上是一个程序的各个部分。其中一些任务是独立的，而另一些任务则依赖于之前的任务，要先完成。举例来说，考虑一下大多数教程程序所具有的
 * <code>setup_dofs</code> 函数的一部分的典型布局。
 * @code
 * 1  dof_handler.distribute_dofs (fe);
 * 2  DoFTools::make_hanging_node_constraints (dof_handler, hanging_node_constraints);
 * 3  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
 * 4  hanging_node_constraints.condense (sparsity_pattern);
 * @endcode
 * 这里，每个操作都需要大量的计算。但要注意的是，并不是所有的操作都相互依赖：显然我们不能在1之前运行语句2-4，4需要等待语句2和3的完成。但语句2和3是独立的：它们可以以任何顺序运行，或以%并行方式运行。实质上，我们已经确定了四个<i>tasks</i>，其中一些是相互依赖的，而其他是独立的。在目前的例子中，任务被识别为单独的C++语句，但往往它们更普遍地与整个代码块相吻合。
 * 这里的重点是这样的。如果我们想使用线程来利用任务2和3的独立性，我们将启动两个线程，并在自己的线程上运行任务2和3；然后我们将等待这两个线程完成（一个称为
 * "加入线程
 * "的操作），并继续执行语句4。实现这一目标的代码是这样的（实际的语法在下面有更详细的解释）。
 * @code
 * dof_handler.distribute_dofs (fe);
 *
 * Threads::Thread<void>
 *   thread_1 = Threads::new_thread (&DoFTools::make_hanging_node_constraints,
 *                                   dof_handler, hanging_node_constraints);
 * Threads::Thread<void>
 *   thread_2 = Threads::new_thread (&DoFTools::make_sparsity_pattern,
 *                                   dof_handler, sparsity_pattern);
 * thread_1.join();
 * thread_2.join();
 * hanging_node_constraints.condense (sparsity_pattern);
 * @endcode
 *
 * 但如果你的电脑只有一个处理器核心，或者我们有两个处理器核心，但已经有一个不同的程序部分在与上面的代码%并行运行，那该怎么办？在这种情况下，上面的代码仍然会启动新的线程，但程序不会运行得更快，因为没有额外的计算资源可用；相反，程序会运行得更慢，因为线程必须被创建和销毁，而且操作系统必须将线程安排到超额的计算资源上。
 * deal.II本身并没有实现对线程的任务调度。为此，我们使用了<a
 * href="http://www.threadingbuildingblocks.org" target="_top">Threading
 * Building Blocks (TBB)
 * library</a>，并为此提供了简单的封装器。TBB抽象了如何启动或停止线程、在单个线程上启动任务等细节，并提供了可以在许多不同系统中移植的接口。
 *
 *
 * @anchor  MTUsing <h3>Using tasks from within deal.II</h3>。
 * 理想情况下，启动任务（以及类似的线程）的语法，对于上面的例子，应该是这样的。
 * @code
 * Threads::Task<void>
 *   thread
 *   = new_task DoFTools::make_hanging_node_constraints (dof_handler,
 *                                                       hanging_node_constraints);
 * @endcode
 * 换句话说，我们希望通过简单地在调用前加上一个关键字（比如这里的
 * <code>new_task</code> ，线程的关键字类似 <code>new_thread</code>
 * ）来表明函数调用应该在一个单独的任务上运行。前缀的调用将返回一个任务的句柄，我们可以用它来等待任务的完成，也可以用它来查询被调用函数的返回值（除非它是空的，就像这里一样）。
 * 由于C++不支持创建新的关键字，我们必须要有一点创造性。所选择的方式是引入一个函数
 * <code>new_task</code>
 * ，该函数将调用的函数以及调用的参数作为参数。
 * <code>new_task</code>
 * 函数被重载，以适应没有、1个、2个和最多9个参数的函数启动任务。在deal.II中，这些函数生活在Threads命名空间中。因此，我们上面尝试做的实际代码看起来是这样的。
 * @code
 * Threads::Task<void>
 *   thread
 *   = Threads::new_task (&DoFTools::make_hanging_node_constraints,
 *                        dof_handler,
 *                        hanging_node_constraints);
 * @endcode
 *
 * 同样，如果我们想在不同的任务上调用一个成员函数，我们可以通过指定调用函数的对象作为函数指针后的第一个参数来实现。
 * @code
 * class C {
 *   public:
 *     double f(int);
 * };
 *
 * int main () {
 *   C c;
 *
 *   // call f(13) as usual, i.e. using the current processor:
 *   c.f(13);
 *
 *   // call f(42) as a separate task, to be scheduled
 *   // whenever processor resources become available:
 *   Threads::Task<double>
 *     task = Threads::new_task (&C::f, c, 42);
 *
 *   // do something else in between:
 *   ...;
 *
 *   // having finished our other business, wait until the task
 *   // above has terminated and get the value returns by c.f(42):
 *   double result = task.return_value();
 * @endcode
 * 在这里，首先注意我们如何传递对象 <code>c</code> （即
 * <code>this</code> pointer the function <code>C::f</code>
 * 会看到），好像它是函数的第一个参数。其次，注意我们如何通过调用
 * Threads::Task::return_value().
 * 在单独的任务上获取函数返回的值，这个函数意味着等待任务的完成，也就是说，最后一行完全相当于
 * @code
 *   task.join ();
 *   double result = task.return_value();
 * @endcode
 *
 * 还要注意，如果 <code>C::f</code>
 * 想启动自己的任务，也是完全有效的。
 * @code
 * class C {
 *   public:
 *     double f(int);
 *   private:
 *     double f1(int);
 *     double f2(int);
 * };
 *
 * double C::f (int i) {
 *   Threads::Task<double> t1 = Threads::new_task (&C::f1,this, i);
 *   Threads::Task<double> t2 = Threads::new_task (&C::f2,this, i);
 *   return t1.return_value() + t2.return_value();
 * }
 *
 * int main () {
 *   C c;
 *
 *   Threads::Task<double>
 *     task = Threads::new_task (&C::f, c, 42);
 *
 *   // do something else in between:
 *   ...;
 *
 *   double result = task.return_value();
 * @endcode
 * 这里，我们让 <code>C::f</code> 计算其返回值为
 * <code>c.f1(i)+c.f2(i)</code>
 * 。如果有足够的CPU资源，那么加法的两个部分以及
 * <code>main()</code>
 * 中的其他东西都将以%的速度并行运行。如果没有，那么我们最终会在其中一个需要返回值的地方阻塞，从而释放出必要的CPU资源来运行所有这些生成的任务来完成。
 *
 *  在许多情况下，比如上面概述的 <code>setup_dofs</code>
 * 函数的介绍性例子，人们可以确定几个独立的作业，它们可以作为任务运行，但必须等待所有的作业在一个点上完成。我们可以通过存储所有
 * Threads::new_task() 调用的返回对象，并对其中每一个调用
 * Threads::Task::join()
 * 来做到这一点。一个更简单的方法是将所有这些任务对象放入一个
 * Threads::TaskGroup
 * 对象中，并一次性地等待所有的任务。然后，代码会是这样的。
 * @code
 * dof_handler.distribute_dofs (fe);
 *
 * Threads::TaskGroup<void> task_group;
 * task_group += Threads::new_task (&DoFTools::make_hanging_node_constraints,
 *                                  dof_handler, hanging_node_constraints);
 * task_group += Threads::new_task (&DoFTools::make_sparsity_pattern,
 *                                  dof_handler, sparsity_pattern);
 * task_group.join_all ();
 * hanging_node_constraints.condense (sparsity_pattern);
 * @endcode
 *
 *  @anchor  MTHow <h3>How scheduling tasks works and when task-based
 * programming is not efficient</h3> 。
 * 任务如何调度运行的确切细节是deal.II用于任务的Threading
 * Building
 * Blocks（TBB）库的内部%。TBB的文档对任务如何被安排到线程中给出了详细的描述，但对实际使用多少个线程却没有提及。然而，一个合理的猜测是，假设TBB创建的线程数与系统中的处理器核心数一样多。这样，它就能充分利用整个系统，而不会有太多的线程让操作系统不得不定期中断，以便其他线程能在可用的处理器核心上运行。
 * 那么问题来了，TBB调度器接受任务并让线程执行它们。线程完全执行任务，也就是说，TBB调度器不会中途中断一个任务，让另一个任务取得一些中途的进展。这确保了缓存总是热的，例如，避免了抢占式中断的开销。
 * 缺点是，只有当线程真正在做一些事情时，CPU核心才会被充分利用，这意味着（i）必须有足够的任务可用，以及（ii）这些任务真正在做一些事情。请注意，这两个条件都必须满足；特别是，这意味着，如果我们已经确定了足够数量的任务，但如果其中一些任务踌躇不前，例如，因为一个任务正在向磁盘写入数据（这个过程中，CPU经常需要等待磁盘完成交易）或正在等待输入，那么CPU内核就没有得到充分利用。其他情况是，任务在其他外部事件上阻塞，例如通过突扰器与其他任务或线程同步。在这种情况下，调度器会让一个任务在一个线程上运行，但并没有注意到这个线程并没有完全利用CPU核心。
 * 在这样的情况下，<i>does</i>创建一个新的线程（见下文 @ref MTThreads "基于线程的并行性"
 * ）是有意义的，操作系统可以在他们等待外部事物时将其搁置，并让不同的线程（例如运行TBB调度的任务的线程）同时使用CPU。
 *
 *  @anchor  MTSimpleLoops <h3>Abstractions for tasks: Simple loops</h3>。
 * 有些循环在数据上的执行体是完全独立的，因此可以以%并行方式执行。TBB库不是先验地将循环分割成固定数量的块并在任务或线程上执行，而是使用以下概念：循环迭代的范围被分割成一定数量的子范围（例如CPU核数的2或3倍），并平均分配给线程；然后线程执行子范围，如果它们完成了工作，就从其他线程那里偷取整个或部分子范围以保持忙碌。这样一来，即使不是每个循环迭代都需要同样多的工作，或者一些CPU核心因为操作系统中断了其他工作而落后，工作也是平衡的。
 * TBB库的原语有点笨拙，所以deal.II为最经常使用的操作提供了包装例程。最简单的一个类似于
 * std::transform
 * 的做法：它需要一个或多个范围的输入运算符，一个输出迭代器和一个函数对象。一个典型的
 * std::transform 的实现会是这样的。
 * @code
 *   template <typename InputIterator1, typename InputIterator,
 *             typename OutputIterator, typename FunctionObject>
 *   void transform (const InputIterator1 &begin_in_1,
 *                   const InputIterator1 &end_in_1,
 *                   const InputIterator2 &begin_in_2,
 *                   const OutputIterator &begin_out,
 *                   FunctionObject       &function)
 *   {
 *     InputIterator1 in_1 = begin_in_1;
 *     InputIterator2 in_2 = begin_in_2;
 *     OutputIterator out  = begin_out;
 *
 *     for (; in_1 != end_in_1; ++in_1, ++in_2, ++out)
 *      out = function(*in_1,in_2);
 *   }
 * @endcode
 *
 * 在很多情况下， <code>function</code>
 * 没有状态，因此我们可以将这个循环分成几个子范围，如上文所解释的。因此，deal.II有一组函数
 * parallel::transform
 * ，看起来和上面的函数一样，但它们是以%parallel的方式进行工作的（对于接受一个、两个或更多参数的函数对象，有几个版本有一个、两个和更多的输入迭代器）。调用这些函数的唯一区别是，它们需要一个额外的最后一个参数，表示
 * <code>[begin_in_1,end_in_1)</code>
 * 的子范围的最小尺寸；它应该足够大，这样我们就不会在调度子范围到处理器上花费更多的时间，但又足够小，处理器可以有效地进行负载平衡。一个经验法则似乎是，如果执行一个子程序需要少于2000条指令，那么这个子程序就太小了。
 * 如何使用这些函数的一个例子是向量操作，如 $z = x+y$
 * 中的加法，所有三个对象都是Vector<Number>类型。
 * @code
 *   parallel::transform (x.begin(), x.end(),
 *                        y.begin(),
 *                        z.begin(),
 *                        [](const Number first, const Number second)
 *                        {
 *                          return first+second;
 *                        },
 *                        1000);
 * @endcode
 *
 * 在这个例子中，我们用一个<i>lambda
 * expression</i>来构建一个函数对象，它接受两个参数并返回两个参数的和。当我们想把向量
 * $x$ 和 $y$ 的各个元素相加，并把两者之和写入 $z$
 * 的元素中时，这正是我们需要的东西。我们在这里得到的函数对象完全为编译器所知，当它展开
 * parallel::transform
 * 所产生的循环时，就像我们以其明显的形式写出的循环一样。
 * @code
 *     InputIterator1 in_1 = x.begin();
 *     InputIterator2 in_2 = y.begin();
 *     OutputIterator out  = z.begin();
 *
 *     for (; in_1 != x.end(); ++in_1, ++in_2, ++out)
 *      out =in_1 +in_2;
 * @endcode
 *
 * 还要注意的是，我们已经确保没有任何一个CPU能得到整个循环中小于1000次迭代的那块（除非整个范围更小）。
 *
 *  @anchor  MTComplexLoops <h3>Abstractions for tasks: More complex
 * loops</h3>。
 * 如果在每个迭代中进行的操作不需要大量的设置成本，并且可以被编译器内联，那么上一节所示的方案是有效的。<i>Lambda
 * expressions</i>正是这种类型，从而消除了调用外部函数的开销。然而，在有些情况下，在每个迭代中调用一些对象或函数是低效的。
 * 这种情况的一个例子是稀疏矩阵-向量乘法。如果你知道数据是如何以压缩行格式存储的，比如在SparseMatrix类中，那么一个矩阵-向量乘积函数看起来是这样的。
 * @code
 *  void SparseMatrix::vmult (const Vector &src,
 *                            Vector       &dst) const
 *  {
 *    const double      val_ptr    = &values[0];
 *    const unsigned intcolnum_ptr = &colnums[0];
 *    Vector::iterator dst_ptr = dst.begin();
 *
 *    for (unsigned int row=0; row<n_rows; ++row, ++dst_ptr)
 *      {
 *        double s = 0.;
 *        const doubleconst val_end_of_row = &values[rowstart[row+1]];
 *        while (val_ptr != val_end_of_row)
 *          s +=val_ptr++ src(*colnum_ptr++);
 *       dst_ptr = s;
 *      }
 *  }
 * @endcode
 * 在for循环中，我们计算矩阵的某一行与右侧向量
 * <code>src</code> 的点积，并将其写入 <code>dst</code>
 * 向量的相应元素中。通过利用<i>next</i>行的元素跟随当前行<i>immediately</i>的元素，使代码更加有效，也就是说，在循环体的开始，我们不必重新设置指向每一行的值和列数的指针。
 * 使用上面的 parallel::transform
 * 函数，原则上我们可以把这段代码写成如下。
 * @code
 *  void SparseMatrix::vmult_one_row (const Vector     &src,
 *                                    Vector           &dst,
 *                                    Vector::iterator &dst_row) const
 *  {
 *    const unsigned int  row = (dst_row
 *
 * - dst.begin());
 *
 *    const double      val_ptr    = &values[rowstart[row]];
 *    const unsigned intcolnum_ptr = &colnums[rowstart[row]];
 *
 *    double s = 0.;
 *    const doubleconst val_end_of_row = &values[rowstart[row+1]];
 *    while (val_ptr != val_end_of_row)
 *      s +=val_ptr++ src(*colnum_ptr++);
 *   dst_row = s;
 *  }
 *
 *  void SparseMatrix::vmult (const Vector &src,
 *                            Vector       &dst) const
 *  {
 *    parallel::transform (dst.begin(), dst.end(),
 *                         std::bind (&SparseMatrix::vmult_one_row,
 *                                    this,
 *                                    std::cref(src),
 *                                    std::ref(dst),
 *                                    std::_1),
 *                         200);
 *  }
 * @endcode
 * 注意我们如何使用 <code>std::bind</code>
 * 来<i>bind</i>某些参数给 <code>vmult_one_row</code>
 * 函数，留下一个参数，从而使 parallel::transform
 * 函数认为传递的函数参数是单值的。还要注意的是，我们需要把源向量和目的向量作为（const）引用，以防止
 * <code>std::bind</code> 按值传递它们（意味着对 <code>src</code>
 * 的拷贝和把结果写入 <code>dst</code>
 * 的临时拷贝，这都不是我们想要的）。最后，注意到一个矩阵的最小200行的粒度，应该由单个CPU核心来处理。
 * 问题是，虽然这样做是正确的，但效率不高：我们必须在循环的每个迭代中设置
 * <code>row, val_ptr, colnum_ptr</code>
 * 变量。此外，由于现在每一行要调用的函数对象不再是简单的<i>lambda
 * expression</i>，在循环的每一次迭代中都有一个隐含的函数调用，包括参数传递。
 * 一个更有效的方法是让TBB将原始范围分割成子范围，然后不是在循环的每个元素上调用目标函数，而是在整个范围上调用。这一点由
 * parallel::apply_to_subranges 函数提供便利。
 * @code
 *  void
 *  SparseMatrix::vmult_on_subrange (const unsigned int  begin_row,
 *                                   const unsigned int  end_row,
 *                                   const Vector     &src,
 *                                   Vector           &dst)
 *  {
 *    const double      val_ptr    = &values[rowstart[begin_row]];
 *    const unsigned intcolnum_ptr = &colnums[rowstart[begin_row]];
 *    Vector::iterator dst_ptr = dst.begin() + begin_row;
 *
 *    for (unsigned int row=begin_row; row<end_row; ++row, ++dst_ptr)
 *      {
 *        double s = 0.;
 *        const doubleconst val_end_of_row = &values[rowstart[row+1]];
 *        while (val_ptr != val_end_of_row)
 *          s +=val_ptr++ src(*colnum_ptr++);
 *       dst_ptr = s;
 *      }
 *  }
 *
 *  void SparseMatrix::vmult (const Vector &src,
 *                            Vector       &dst) const
 *  {
 *     parallel::apply_to_subranges (0, n_rows(),
 *                                   std::bind (vmult_on_subrange,
 *                                              this,
 *                                              std::_1, std::_2,
 *                                              std::cref(src),
 *                                              std::ref(dst)),
 *                                   200);
 *  }
 * @endcode
 * 这里，我们在每个元素至少200个的子范围上调用
 * <code>vmult_on_subrange</code>
 * 函数，这样初始设置成本就可以摊销。
 * 一个相关的操作是当元素上的循环各自产生一个结果，然后必须累积起来（除了数字的加法之外，其他的减少操作也可以）。一个例子是形成矩阵规范
 * $x^T M x$ （如果 $M$
 * 是正定的，它才是真正的规范，但我们暂时假设它是）。对于稀疏矩阵来说，一个顺序的实现会是这样的。
 * @code
 *  double SparseMatrix::mat_norm (const Vector &x) const
 *  {
 *    const double      val_ptr    = &values[0];
 *    const unsigned intcolnum_ptr = &colnums[0];
 *
 *    double norm_sqr = 0;
 *
 *    for (unsigned int row=0; row<n_rows; ++row, ++dst_ptr)
 *      {
 *        double s = 0.;
 *        const doubleconst val_end_of_row = &values[rowstart[row+1]];
 *        while (val_ptr != val_end_of_row)
 *          s +=val_ptr++ x(*colnum_ptr++);
 *        norm_sqr += x(row) s;
 *      }
 *
 *    return std::sqrt (norm_sqr);
 *  }
 * @endcode
 *
 * 如果我们能把这个操作分成几个子行，每个子行都计算自己的那部分规范的平方，把各个子行的结果加在一起，然后取结果的平方根，那就更好了。这就是
 * parallel::accumulate_from_subranges
 * 函数所做的（注意，你必须指定结果类型作为模板参数，而且像往常一样，最后一个参数是可以在单个CPU核上调度的外循环元素的最小数量）。
 * @code
 *  double
 *  SparseMatrix::mat_norm_sqr_on_subrange (const unsigned int  begin_row,
 *                                          const unsigned int  end_row,
 *                                          const Vector     &x)
 *  {
 *    const double      val_ptr    = &values[rowstart[begin_row]];
 *    const unsigned intcolnum_ptr = &colnums[rowstart[begin_row]];
 *    Vector::iterator dst_ptr = dst.begin() + begin_row;
 *
 *    double norm_sqr = 0;
 *
 *    for (unsigned int row=begin_row; row<end_row; ++row, ++dst_ptr)
 *      {
 *        double s = 0.;
 *        const doubleconst val_end_of_row = &values[rowstart[row+1]];
 *        while (val_ptr != val_end_of_row)
 *          s +=val_ptr++ x(*colnum_ptr++);
 *        norm_sqr += x(row) s;
 *      }
 *
 *    return norm_sqr;
 *  }
 *
 *  double SparseMatrix::mat_norm (const Vector &x) const
 *  {
 *    return
 *      std::sqrt
 *      (parallel::accumulate_from_subranges (0, n_rows(),
 *                                            std::bind (mat_norm_sqr_on_subrange,
 *                                                       this,
 *                                                       std::_1, std::_2,
 *                                                       std::cref(x)),
 *                                            200));
 *  }
 * @endcode
 *
 *  @anchor  MTWorkStream <h3>Abstractions for tasks: Work streams</h3>。
 * 在介绍中所示的例子中，我们已经确定了一些可以作为独立任务运行的函数。理想情况下，这个任务的数量要大于CPU核心的数量（以保持它们的忙碌），但也不要过于庞大（以免数百万的任务淹没调度器，然后不得不分配给2或4个核心，例如）。然而，在有些情况下，我们有几千甚至几百万个相对独立的作业：例如，在网格的每个单元上组装对全局线性系统的局部贡献；在每个单元上评估误差估计器；或者在每个单元上对计算的数据进行后处理，这些都属于这一类。这些情况可以用我们称之为WorkStream的软件设计模式来处理。在下文中，我们将介绍这种模式的原理及其实现；更多的细节以及用它可以实现的加速的例子在 @ref workstream_paper "WorkStream论文 "
 * 中给出。 像这样的代码可以这样写。
 * @code
 * template <int dim>
 * void MyClass<dim>::assemble_on_one_cell (const typename DoFHandler<dim>::active_cell_iterator &cell)
 * { ... }
 *
 * template <int dim>
 * void MyClass<dim>::assemble_system ()
 * {
 *   Threads::TaskGroup<void> task_group;
 *   for (typename DoFHandler<dim>::active_cell_iterator
 *          cell = dof_handler.begin_active();
 *        cell != dof_handler.end(); ++cell)
 *     task_group += Threads::new_task (&MyClass<dim>::assemble_on_one_cell,
 *                                     this,
 *                                      cell);
 *   task_group.join_all ();
 * }
 * @endcode
 * 在一个大的网格上，可能有一百万个单元，这会产生大量的任务；虽然它会让所有的CPU核忙上一阵子，但首先创建这么多的任务，对它们进行调度，然后等待它们，这样的开销可能不会导致高效的代码。一个更好的策略是，如果调度器能够以某种方式表明它有可用的资源，在这一点上，我们将给它提供另一个新创建的任务，我们将这样做，直到我们的任务用完，并且所创建的任务已经被工作。
 * 这基本上就是 WorkStream::run
 * 函数所做的。你给它一个迭代器范围，它可以从中抽取对象进行工作（在上面的例子中，它是由
 * <code>dof_handler.begin_active()</code> 到 <code>dof_handler.end()</code>
 * 给出的区间），以及一个对每个项目进行工作的函数（函数
 * <code>MyClass::assemble_on_one_cell</code>
 * ）和一个对象，如果它是一个成员函数。
 * 在下文中，让我们阐述一下为什么WorkStream命名空间中的函数会以这种方式实现的理由。关于其实现的更多信息可以在 @ref workstream_paper "WorkStream文件 "
 * 中找到。要看到WorkStream类在类似上述任务中的实际应用，请看
 * step-9 、 step-13 、 step-14 、 step-32 、 step-35 或 step-37
 * 的教程程序。 首先，考虑到上面的简单描述，那么
 * <code>MyClass::assemble_system</code>
 * 函数的写法可以是这样的（注意，这并不是很正确的语法，下面会介绍）。
 * @code
 * template <int dim>
 * void MyClass<dim>::assemble_on_one_cell (const typename DoFHandler<dim>::active_cell_iterator &cell)
 * { ... }
 *
 * template <int dim>
 * void MyClass<dim>::assemble_system ()
 * {
 *   WorkStream::run (dof_handler.begin_active(),
 *                    dof_handler.end(),
 *                   this,
 *                    &MyClass<dim>::assemble_on_one_cell);
 * }
 * @endcode
 *
 * 然而，这至少有三个问题。 <ul>   <li>  首先，让我们看一下 <code>MyClass::assemble_on_one_cell</code> 函数可能的样子。
 * @code
 * template <int dim>
 * void MyClass<dim>::assemble_on_one_cell (const typename DoFHandler<dim>::active_cell_iterator &cell)
 * {
 *   FEValues<dim> fe_values (...);
 *   FullMatrix<double> cell_matrix (...);
 *   Vector<double>     cell_rhs (...);
 *
 *   // assemble local contributions
 *   fe_values.reinit (cell);
 *   for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
 *     for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
 *       for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
 *         cell_matrix(i,j) += ...;
 *   ...same for cell_rhs...
 *
 *   // now copy results into global system
 *   std::vector<unsigned int> dof_indices (...);
 *   cell->get_dof_indices (dof_indices);
 *   for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
 *     for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
 *       system_matrix.add (dof_indices[i], dof_indices[j],
 *                          cell_matrix(i,j));
 *   ...same for rhs...
 * }
 * @endcode
 *
 * 这里的问题是，几个任务，每个都在运行
 * <code>MyClass::assemble_on_one_cell</code> ，可能会试图写入对象
 * <code>MyClass::system_matrix</code>  <i>at the same
 * time</i>。这可以通过使用 Threads::Mutex,
 * 的显式同步来避免，例如，看起来像这样。
 * @code
 *   // now copy results into global system
 *   std::vector<unsigned int> dof_indices (...);
 *   cell->get_dof_indices (dof_indices);
 *
 *   static Threads::Mutex mutex;
 *   mutex.acquire ();
 *   for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
 *     for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
 *       system_matrix.add (dof_indices[i], dof_indices[j],
 *                          cell_matrix(i,j));
 *   ...same for rhs...
 *   mutex.release ();
 * }
 * @endcode
 *
 * 通过使mutex成为静态变量，它在全局范围内只存在一次（即对所有可能在%parallel中运行的任务都存在一次），并且只有一个任务可以进入由mutex上的acquisition/release调用保护的区域。顺便说一句，写这段代码的更好方法是这样的，确保即使在抛出异常的情况下也能释放mutex，而且不需要记住写对
 * Threads::Mutex::release(): 的调用。
 * @code
 *   // now copy results into global system
 *   static Threads::Mutex mutex;
 *   Threads::Mutex::ScopedLock lock (mutex);
 *   for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
 *     for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
 *       system_matrix.add (dof_indices[i], dof_indices[j],
 *                          cell_matrix(i,j));
 *   ...same for rhs...
 * }
 * @endcode
 * 在这里，从ScopedLock被创建到它被销毁，在代码块结束时，mutex一直被锁定。
 * 请注意，尽管我们现在避免了多个线程可能向同一个对象写入的竞赛条件，但这段代码的效率并不高：在多核机器上，mutexes是很昂贵的，而且我们还在某些时候阻塞了线程，这对于任务来说是低效的，正如上文 @ref MTHow "调度任务如何工作以及基于任务的编程何时不高效 "
 * 一节中所解释的。
 * <li>
 * 第二个正确性问题是，即使我们使用mutex锁定全局矩阵和右手边的对象，我们也是以一种或多或少的随机顺序这样做的：虽然任务是按照我们正常遍历单元格的顺序创建的，但不能保证当我们到了要把局部复制到全局贡献的时候，顺序仍然是我们按顺序计算的那样。换句话说，我们可能会把单元格1的贡献加在单元格0的贡献之前。这看起来无害，因为加法是交换和关联的，但事实上，如果用浮点运算，就不是这样了。
 * $a+b+c \neq a+c+b$
 *
 * - 以 $a=1, b=-1, c=10^{-20}$ 为例（因为 $1+10^{-20}=1$ 在浮点运算中，使用双精度）。
 * 结果是，最终出现在全局矩阵和右手边的确切数值会很接近，但可能会因任务完成的顺序不同而有接近四舍五入的差异。这不是一个理想的结果，因为这样的结果是不可复制的。
 * 因此，WorkStream类的设计方式是使用两个函数：
 * <code>MyClass::assemble_on_one_cell</code>
 * 计算本地贡献并将其存储在某个地方（我们接下来会讨论这个问题），第二个函数，例如
 * <code>MyClass::copy_local_to_global</code>
 * ，将每个单元的计算结果复制到全局对象。在WorkStream类中实现的技巧是：(i)
 * <code>MyClass::copy_local_to_global</code>
 * 永远不会以%并行方式运行超过一次，所以我们不需要通过互斥来同步执行，(ii)它在单元格上的运行顺序与它们在迭代器范围内出现的顺序完全相同，也就是说，我们以同样的方式将元素加入全局矩阵中<i>every
 * time, independently of when the computation of these element
 * finishes</i>。 我们现在只需要讨论
 * <code>MyClass::assemble_on_one_cell</code> 如何向
 * <code>MyClass::copy_local_to_global</code>
 * 传达它所计算的内容。这样做的方法是使用一个持有所有临时数据的对象。
 * @code
 * struct PerTaskData {
 *   FullMatrix<double>        cell_matrix;
 *   Vector<double>            cell_rhs;
 *   std::vector<unsigned int> dof_indices;
 * }
 *
 * template <int dim>
 * void MyClass<dim>::assemble_on_one_cell (const typename DoFHandler<dim>::active_cell_iterator &cell,
 *                                          PerTaskData &data)
 * {
 *   FEValues<dim> fe_values (...);
 *
 *   data.cell_matrix = 0;
 *   data.cell_rhs    = 0;
 *
 *   // assemble local contributions
 *   fe_values.reinit (cell);
 *   for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
 *     for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
 *       for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
 *         data.cell_matrix(i,j) += ...;
 *   ...same for cell_rhs...
 *
 *   cell->get_dof_indices (data.dof_indices);
 * }
 *
 * template <int dim>
 * void MyClass<dim>::copy_local_to_global (const PerTaskData &data)
 * {
 *   for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
 *     for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
 *       system_matrix.add (data.dof_indices[i], data.dof_indices[j],
 *                          data.cell_matrix(i,j));
 *   ...same for rhs...
 * }
 *
 * template <int dim>
 * void MyClass<dim>::assemble_system ()
 * {
 *   PerTaskData per_task_data;
 *   ...initialize members of per_task_data to the correct sizes...
 *
 *   WorkStream::run (dof_handler.begin_active(),
 *                    dof_handler.end(),
 *                   this,
 *                    &MyClass<dim>::assemble_on_one_cell,
 *                    &MyClass<dim>::copy_local_to_global,
 *                    per_task_data);
 * }
 * @endcode
 *
 * 这种工作方式是，我们创建一个样本 <code>per_task_data</code>
 * 对象，工作流对象将在每个以%并行方式运行的任务中复制一次。对于每个任务，这个对象将首先被传递给以%并行方式运行的
 * <code>MyClass::assemble_on_one_cell</code>
 * 的几个实例中的一个，该实例将在单个单元上获得的数据填入该对象，然后再传递给顺序运行的
 * <code>MyClass::copy_local_to_global</code>
 * ，将数据复制到全局对象中。当然，在实践中，如果我们有数以百万计的单元，我们不会产生数以百万计的
 * <code>per_task_data</code> 对象；相反，我们在这些对象被
 * <code>MyClass::copy_local_to_global</code>
 * 使用后进行回收，并将其送回
 * <code>MyClass::assemble_on_one_cell</code>
 * 的另一个实例；这意味着我们实际创建的此类对象的数量是调度器使用的线程数的一小部分，通常与系统中的CPU核数量相当。
 * <li>  最后一个值得解决的问题是，按照上面
 * <code>MyClass::assemble_on_one_cell</code>
 * 函数的写法，我们在每次调用该函数时都要创建和销毁一个FEValues对象，也就是说，为三角形中的每个单元创建一次。这是一个非常昂贵的操作，因为FEValues类试图在其构造函数中做大量的工作，试图减少我们对每个单元的操作数量（即增加
 * ${\cal O}(1)$
 * 中的常数来初始化这样一个对象，以减少在三角形的 $N$
 * 单元上调用 ${\cal O}(N)$
 * 操作的常数）。在每个单元上创建和销毁一个FEValues对象会使这种努力失效。
 * 避免这种情况的方法是将FEValues对象放入第二个结构中，该结构将保存从头开始的数据，并在构造函数中初始化它。
 * @code
 * struct PerTaskData {
 *   FullMatrix<double>        cell_matrix;
 *   Vector<double>            cell_rhs;
 *   std::vector<unsigned int> dof_indices;
 *
 *   PerTaskData (const FiniteElement<dim> &fe)
 *              :
 *              cell_matrix (fe.dofs_per_cell, fe.dofs_per_cell),
 *              cell_rhs (fe.dofs_per_cell),
 *              dof_indices (fe.dofs_per_cell)
 *     {}
 * }
 *
 * struct ScratchData {
 *   FEValues<dim>             fe_values;
 *
 *   ScratchData (const FiniteElement<dim> &fe,
 *                const Quadrature<dim>    &quadrature,
 *                const UpdateFlags         update_flags)
 *              :
 *              fe_values (fe, quadrature, update_flags)
 *     {}
 *
 *   ScratchData (const ScratchData &scratch)
 *              :
 *              fe_values (scratch.fe_values.get_fe(),
 *                         scratch.fe_values.get_quadrature(),
 *                         scratch.fe_values.get_update_flags())
 *     {}
 * }
 * @endcode
 * 然后在组装函数中使用这个FEValues对象。
 * @code
 * template <int dim>
 * void MyClass<dim>::assemble_on_one_cell (const typename DoFHandler<dim>::active_cell_iterator &cell,
 *                                          ScratchData &scratch,
 *                                          PerTaskData &data)
 * {
 *   scratch.fe_values.reinit (cell);
 *   ...
 * }
 * @endcode
 * 就像 <code>PerTaskData</code> 结构一样，我们将创建一个
 * <code>ScratchData</code>
 * 的样本对象，并将其传递给工作流对象，后者将根据需要多次复制它。为了使其发挥作用，
 * <code>ScratchData</code>
 * 结构需要是可复制的。由于FEValues对象相当复杂，不能隐式复制，我们为
 * <code>ScratchData</code>
 * 结构提供了我们自己的复制构造函数。
 * 同样的方法，把东西放到 <code>ScratchData</code>
 * 数据结构中，应该用于所有昂贵的构造。这尤其适用于所有在构造时需要分配内存的东西；例如，如果一个函数的值需要在正交点进行评估，那么这就很昂贵。
 * @code
 * template <int dim>
 * void MyClass<dim>::assemble_on_one_cell (const typename DoFHandler<dim>::active_cell_iterator &cell,
 *                                          ScratchData &scratch,
 *                                          PerTaskData &data)
 * {
 *   std::vector<double> rhs_values (fe_values.n_quadrature_points);
 *   rhs_function.value_list (data.fe_values.get_quadrature_points,
 *                            rhs_values)
 *   ...
 * }
 * @endcode
 * 而这是一个更便宜的方法。
 * @code
 * struct ScratchData {
 *   std::vector<double>       rhs_values;
 *   FEValues<dim>             fe_values;
 *
 *   ScratchData (const FiniteElement<dim> &fe,
 *                const Quadrature<dim>    &quadrature,
 *                const UpdateFlags         update_flags)
 *              :
 *              rhs_values (quadrature.size()),
 *              fe_values (fe, quadrature, update_flags)
 *     {}
 *
 *    ScratchData (const ScratchData &scratch)
 *              :
 *              rhs_values (scratch.rhs_values),
 *              fe_values (scratch.fe_values.get_fe(),
 *                         scratch.fe_values.get_quadrature(),
 *                         scratch.fe_values.get_update_flags())
 *     {}
 * }
 *
 * template <int dim>
 * void MyClass<dim>::assemble_on_one_cell (const typename DoFHandler<dim>::active_cell_iterator &cell,
 *                                          ScratchData &scratch,
 *                                          PerTaskData &data)
 * {
 *   rhs_function.value_list (scratch.fe_values.get_quadrature_points,
 *                            scratch.rhs_values)
 *   ...
 * }
 * @endcode
 *
 * </ul>
 * 作为最后一点。如果由于某种原因，我的汇编器和复制器函数不符合上述分别有三个和一个参数的签名，怎么办？这也不是什么问题。WorkStream
 * 命名空间提供了  WorkStream::run()
 * 函数的两个版本：一个是接收一个对象和两个成员函数的地址，另一个是简单地接收两个函数对象，可以分别用三个和一个参数调用。所以，换句话说，下面的两个调用是完全相同的。
 * @code
 *   WorkStream::run (dof_handler.begin_active(),
 *                    dof_handler.end(),
 *                   this,
 *                    &MyClass<dim>::assemble_on_one_cell,
 *                    &MyClass<dim>::copy_local_to_global,
 *                    per_task_data);
 *   // ...is the same as:
 *   WorkStream::run (dof_handler.begin_active(),
 *                    dof_handler.end(),
 *                    std::bind(&MyClass<dim>::assemble_on_one_cell,
 *                             this,
 *                              std::_1,
 *                              std::_2,
 *                              std::_3),
 *                    std::bind(&MyClass<dim>::copy_local_to_global,
 *                             this,
 *                              std::_1),
 *                    per_task_data);
 * @endcode
 * 注意 <code>std::bind</code> 是如何通过将成员函数绑定到
 * <code>*this</code>
 * 对象上产生一个需要三个参数的函数对象。   <code>std::_1,
 * std::_2</code> and <code>std::_3</code>
 * 是第一个、第二个和第三个参数的占位符，可以在以后指定。换句话说，例如，如果
 * <code>p</code>  是第一次调用  <code>std::bind</code>
 * 的结果，那么调用 <code>p(cell, scratch_data, per_task_data)</code>
 * 将导致执行  <code>this-@>assemble_on_one_cell (cell, scratch_data,
 * per_task_data)</code>  ，即  <code>std::bind</code>
 * 已经将对象绑定到函数指针上，但为以后留下了三个参数。
 * 同样地，让我们假设 <code>MyClass::assemble_on_one_cell</code>
 * 在一个非线性、时间依赖性问题的求解器中具有如下签名。
 * @code
 * template <int dim>
 * void
 * MyClass<dim>::assemble_on_one_cell (const Vector<double> &linearization_point,
 *                                     const typename DoFHandler<dim>::active_cell_iterator &cell,
 *                                     ScratchData &scratch,
 *                                     PerTaskData &data,
 *                                     const double current_time)
 * { ... }
 * @endcode
 * 因为WorkStream希望能够只用三个参数来调用worker函数，第一个参数是迭代器，第二个和第三个参数是ScratchData和PerTaskData对象，所以我们需要向它传递以下内容。
 * @code
 *   WorkStream::run (dof_handler.begin_active(),
 *                    dof_handler.end(),
 *                    std::bind(&MyClass<dim>::assemble_on_one_cell,
 *                             this,
 *                              current_solution,
 *                              std::_1,
 *                              std::_2,
 *                              std::_3,
 *                              previous_time+time_step),
 *                    std::bind(&MyClass<dim>::copy_local_to_global,
 *                             this,
 *                              std::_1),
 *                    per_task_data);
 * @endcode
 * 这里，我们将对象、线性化点参数和当前时间参数绑定到函数中，然后再交给
 * WorkStream::run().   WorkStream::run()
 * 将简单地用单元格和scratch以及每个任务对象调用函数，这些对象将在
 * <code>std::_1, std::_2</code>  和  <code>std::_3</code>
 * 指示的位置被填入。
 * 上面显示的 WorkStream::run 函数有一些细化。例如，人们可能会意识到，只有在复制-本地到全局的函数比本地装配函数快得多的情况下，上面的基本想法才能扩展，因为前者必须按顺序运行。这种限制只能通过安排更多的并行工作来改善。这就导致了这样一个概念：通过记录哪些写操作相互冲突，为我们工作的单元格（或更广泛的迭代器）的图形着色。因此，有一个 WorkStream::run 的第三个版本，它不只是取一个迭代器的范围，而是取一个由可以同时工作的元素组成的向量。这个概念在 @ref workstream_paper 的 "WorkStream论文 "
 * 中得到了非常详细的解释，同时还有常见例子的性能评估。
 *
 *  @anchor  MTTaskSynchronization <h3>Tasks and synchronization</h3>。
 * 任务是强大的，但它们确实有其局限性：为了使事情变得高效，任务调度器从不自己中断任务。除了调用
 * Threads::Task::join
 * 函数来等待另一个任务完成的情况外，任务调度器总是将一个任务运行到完成。缺点是，调度器看不到一个任务是否真的在空转，例如，如果它在等待其他事情发生（文件IO完成，来自键盘的输入等）。在这样的情况下，任务调度器原则上可以运行一个不同的任务，但由于它不知道任务在做什么，所以它不知道。因此，那些确实在等待外部事件发生的函数不适合做任务，应该使用线程（见下文）。
 * 然而，在有些情况下，任务不仅是对工作的一种糟糕的抽象，而且实际上可以不使用。原则上，任务不能通过使用mutex或条件变量与其他任务同步（参见
 * Threads::Mutex 和 Threads::ConditionVariable
 * 类）。原因是，如果任务A需要等待任务B完成某件事情，那么只有在保证任务B最终能够运行并完成任务的情况下，这才会奏效。现在想象一下，你有2个处理器，任务A1和A2目前正在运行；我们假设他们已经排好了任务B1和B2的队列，现在正用一个突变器等待这些排队的任务完成（部分）工作。由于机器只有两个处理器，任务调度器只有在A1或A2完成后才会启动B1或B2
 *
 * 但这并没有发生，因为它们是使用操作系统资源（mutex）而不是任务调度器资源在等待。结果是一个死锁。
 * 底线是，任务不能使用互斥或条件变量来与其他任务同步。如果任务之间的通信是必要的，你需要使用线程，因为操作系统确保所有线程最终都能运行，与线程的总数无关。然而，请注意，如果你只在每个任务上分别使用一个
 * Thread::Mutex
 * 来保护对任务可能写入的变量的访问，那么情况就不一样了：这种对mutexes的使用是可以的；任务可能只是不想等待另一个任务做什么。
 *
 *  @anchor  MTThreads <h3>Thread-based parallelism</h3> 。
 * 尽管任务是一种更高层次的描述事物的方式，但有些情况是不适合用任务的（关于其中一些情况的讨论，见上文 @ref MTHow "调度任务是如何进行的以及基于任务的编程在什么情况下是无效的"
 * ）。一般来说，不能完全利用CPU的工作不适合任务，而适合线程。
 * 在这样的情况下，你可以采用显式启动线程，而不是任务，使用的语法与上面基本相同。例如，如果你的应用程序中有一个生成图形输出的函数，然后估计误差，为自适应网格方案的下一次迭代细化网格，它可以是这样的。
 * @code
 * template <int dim>
 * void MyClass<dim>::output_and_estimate_error () const
 * {
 *   DataOut<dim> data_out;
 *   data_out.attach_dof_handler (dof_handler);
 *   data_out.add_data_vector (solution, "solution");
 *   data_out.build_patches ();
 *
 *   std::ofstream output ("solution.vtk");
 *
 *   Threads::Thread<void>
 *     thread = Threads::new_thread (&DataOut<dim>::write_vtk, data_out, output);
 *
 *   Vector<float> error_per_cell (triangulation.n_active_cells());
 *   KellyErrorEstimator<dim>::estimate (dof_handler,
 *                                      QGauss<dim-1>(3),
 *                                      typename FunctionMap<dim>::type(),
 *                                      solution,
 *                                      estimated_error_per_cell);
 *   thread.join ();
 * @endcode
 *
 * 在这里， Threads::new_thread
 * 启动了给定的函数，该函数在一个新的线程上向输出文件写入数据，该线程可以与其他一切并行运行：在与
 * KellyErrorEstimator::estimate() 函数并行的%， DataOut::write_vtk()
 * 函数将在一个单独的线程上运行。这种执行是独立于负责任务的调度器的，但这并不是一个问题，因为向文件写入大量数据并不会让CPU非常忙碌。
 * 创建线程的工作方式与任务基本相同，即你可以用
 * Threads::Thread::join(), 来等待一个线程的终止，用
 * Threads::Thread::return_value(),
 * 来查询一个已完成的线程的返回值，你可以将线程分组到一个
 * Threads::ThreadGroup 对象中，等待所有线程的完成。
 *
 *  @anchor  MTTaskThreads <h3>Controlling the number of threads used for
 * tasks</h3>
 * 如前所述，deal.II没有实现向线程调度任务，甚至没有启动线程本身。TBB库在决定使用多少个线程方面做得很好，他们不建议明确设置线程的数量。然而，在大型对称多处理（SMP）机器上，特别是那些有资源/任务管理器的机器，或者在对某些部分内存的访问是可能的，但对远处的处理器来说非常昂贵的系统上（例如，非常大的NUMA
 * SMP机器），可能有必要明确设置线程数，以防止TBB使用过多的CPU。另一个用例是，如果你在一台机器上运行多个MPI作业，并且每个作业只应使用可用处理器内核的一个子集。
 * 明确设置线程数是通过在调用其他可能创建线程的函数之前调用
 * MultithreadInfo::set_thread_limit()
 * 来完成的。在实践中，它应该是你在  <code>main()</code>
 * 中首先调用的函数之一。
 * 如果你用MPI运行你的程序，那么你可以使用MPI_InitFinalize类的构造函数的可选第三个参数来实现相同的目标。
 *
 *
 * @note
 * deal.II内部有一小部分地方也明确使用了基于线程的并行性，例如用于运行需要等待输入或输出发生的后台任务，因此不会消耗太多的CPU时间。这种线程不在TBB任务调度器的控制下运行，因此不受上述程序的影响。在某些情况下，deal.II也会调用BLAS库，它有时也可能启动自己的线程。你将不得不查阅你的BLAS安装文档，以确定如何为这些操作设置线程数。
 *
 */


