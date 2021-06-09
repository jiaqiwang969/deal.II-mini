//include/deal.II-translator/A-headers/exceptions_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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
 *    @defgroup Exceptions Exceptions and assertions
 * 该模块包含在deal.II的异常机制中使用的类。 <h2>Brief
 * overview</h2>
 * 异常有两种不同的使用方式。   <ul>
 * <li>
 * 静态断言。这些是只在调试模式下启用的检查，而不是在发布（或优化，生产）模式下。在deal.II中，静态断言通常用于检查函数的参数是否满足某些属性，内部数据结构是否一致，以及类似的断言。例如，静态断言用于确保两个相加的向量具有相同的组件数量
 *
 * - 其他的一切反正都没有任何意义。
 * 这种检查是由 @p Assert
 * 宏在库中的几千个地方进行的。另外，从 step-5
 * 开始的几个教程程序显示了如何做到这一点。
 * 如果一个静态断言被违反了，异常机制会产生一个异常类型，指出到底是什么出了问题，显示适当的信息，包括检测到问题的确切位置，然后中止程序
 *
 * - 如果你试图添加两个不同长度的向量，在程序中没有什么可以应对的情况，你必须去修正程序代码来代替。一般来说，甚至没有理由使用通常的C++异常机制来 @p throw 一个异常对象，因为在这种情况下，更高一级的函数没有办法纠正这种情况，并以一种有用的方式处理它
 *
 * - 并不是程序收到了坏的数据；程序只是出现了错误，人们无法智能地解决这个问题。
 * （有时将 @p Assert
 * 宏的行为从中止程序改为抛出异常是很有用的。另一方面，异常不允许从类的析构器中传播出去。
 * 为此，有一个叫做  @p AssertNothrow,
 * 的宏的变体，可以在析构器中使用。这些用例将在本页面下面进一步讨论)。
 *
 * <li>
 * 动态断言。这些用于检查依赖于外部事物的条件，这些外部事物可能在一次程序运行中与下一次不同，例如，一个输出文件是否可以被写入。
 * 这些是不应该被静态检查的东西，因为不能保证在调试模式下满足条件的程序，在随后的发布模式下也会满足该条件。
 *
 * - 换句话说，只在调试模式下检查这些情况是不够的。
 * 相反，我们必须在程序的执行过程中每次都要检查这些情况。在deal.II中，这是用 step-9 、 step-13 和以下教程中介绍的 @p AssertThrow 宏来完成的。该宏检查一个条件，如果违反了，就使用C++  <code>throw</code> 机制抛出一个本模块中声明的类型之一的异常。由于这些是运行时异常，这就给了程序捕捉异常的机会，例如，将输出写入一个可写文件。   </ul>
 *
 *  <h2>Detailed description</h2>
 * <tt>deal.II</tt>中的错误处理机制通常以两种方式使用。
 * 第一种是只在调试模式下使用错误检查，对没有经过全面测试的程序很有用。当程序不再显示错误时，可以关闭错误处理，并以此获得更好的性能，因为在库中对错误的检查是相当频繁的（典型的速度提升是4倍！）。这种异常生成模式对于内部一致性检查是最有用的，比如范围检查或函数参数的有效性检查。这种类型的错误通常是编程错误，程序终止时应该有尽可能详细的信息，包括异常产生的位置和原因。
 * 第二种模式是用于错误检查，这种检查应该始终处于开启状态，比如I/O错误、内存请求失败等等。关掉这个模式没有什么意义，因为这种错误同样可能发生在经过测试和未经测试的程序中。这类异常不会终止程序，而是以<tt>C++</tt>的方式抛出异常，允许程序捕捉它们并最终做一些处理。由于在异常不能被正确处理的情况下，打印出一些信息可能是有用的，所以额外的信息会像第一种模式一样被传递。后者使得有必要提供一系列的宏，将这些额外的信息输入到异常类中；原则上，这可以由程序员自己每次手工完成，但由于这些信息可以自动获得，所以为此提供了一个宏。
 * 这两种模式都使用异常类，除了<tt>C++</tt>标准的
 * <tt>std::exception</tt> 类之外，它们还需要有特殊的功能。
 * 这样的类是由以下几行代码声明的。
 * @code
 *   DeclException2 (ExcDomain, int, int,
 *                   << "Index= " << arg1 << "Upper Bound= " << arg2);
 * @endcode
 *
 * 这声明了一个名为<tt>ExcDomain</tt>的异常类，它有两个变量作为附加信息（默认命名为<tt>arg1</tt>和<tt>arg2</tt>），它输出给定的序列（它被附加到一个
 * <tt>std::ostream</tt>
 * 变量的名称上，因此有奇怪的语法）。还有其他<tt>DeclExceptionN</tt>宏，用于有更多或没有参数的异常类。按照惯例，所有异常类的名字都是以<tt>Exc...</tt>开头的，而且大多数都是在本地声明的，以用于它要使用的类（少数非常频繁的异常也是在StandardExceptions命名空间中声明的，在任何地方都可以使用）。在全局范围内声明异常是可能的，但是会污染全局命名空间，可读性较差，而且大多数时候是不必要的。
 * 由于异常类在两种错误检查模式下的声明方式相同，所以可以使用通过<tt>DeclExceptionN(...)</tt>宏系列声明的异常来进行静态和动态检查。
 *
 *  <h3>Use of the debug mode exceptions (static checks)</h3>
 * 要使用异常机制进行调试模式的错误检查，请在你的源代码中写下类似以下的行。
 * @code
 *  Assert (n<dim, ExcDomain(n,dim));
 * @endcode
 * 通过宏扩展，它基本上做了以下工作（尽管实际的代码稍微复杂一些）。
 * @code
 *  #ifdef DEBUG
 *      if (!(cond))
 *        {
 *          // issue error of class ExcDomain(n,dim)
 *        }
 *  #else
 *      // do nothing
 *  #endif
 * @endcode
 * 也就是说，只有当预处理器变量<tt>DEBUG</tt>被设置，并且违反了给定的条件（在这种情况下<tt>n
 * < dim</tt>），它才会发出错误。
 * 如果异常是用<tt>DeclException0
 * (...)</tt>宏来声明的，也就是说，没有任何附加参数，那么它的名字就必须用括号给出。
 * <tt>Assert (i>m, ExcSomewhat());</tt> <h4>How it works internally</h4>
 * 如果设置了<tt>DEBUG</tt>预处理器指令，调用<tt>Assert (cond,
 * exc);</tt>基本上会被预处理器转换为以下序列。
 * @code
 *  if (!(cond))
 *    deal_II_exceptions::internals::issue_error_noreturn
 *           (__FILE__,
 *            __LINE__,
 *            __PRETTY_FUNCTION__,
 *            #cond,
 *            #exc,
 *            exc);
 * @endcode
 *
 * （注意，函数名称和确切的调用序列可能会随着时间的推移而改变，但一般原则是不变的）。也就是说，如果给定的条件被违反，那么发生异常的文件和行以及条件本身和异常对象的调用序列就会被传递给
 * deal_II_exceptions::internals::issue_error_noreturn()
 * 函数。此外，一个由<tt>exc</tt>给出的对象被创建（这通常是一个未命名的对象，如<tt>ExcDomain(n,
 * dim)</tt>类的<tt>ExcDomain</tt>）并被转移到这个函数。
 * <tt>__PRETTY_FUNCTION__</tt>是一些编译器定义的宏，给出了函数的名称。如果使用的是另一个编译器，如果编译器为我们提供了这个函数，我们就尝试将其设置为合理的函数，否则就<tt>"（不可用）"</tt>。
 * 在<tt>issue_error_noreturn</tt>中，通过调用set_fields()函数，将给定的数据转移到<tt>exc</tt>对象中；之后，程序要么被中止（关于异常的信息被打印到deallog），要么抛出异常。<tt>Assert</tt>宏做了第一条路径（打印和中止）；<tt>AssertThrow</tt>做了第二条（抛出）。这种行为与本文前面对静态和动态断言的描述是一致的。如果能从操作系统中获得，输出也可能包含堆栈跟踪，以显示错误发生的位置。 @ref
 * Tutorial 的几个程序显示了一个典型的输出。
 * 如果预处理器变量<tt>DEBUG</tt>没有设置，那么<tt>Assert</tt>宏就会扩展为<tt>{}</tt>。
 * 有时，除了程序流不应该到达某个点之外，没有其他有用的异常条件，例如，<tt>switch</tt>语句的<tt>default</tt>部分。在这种情况下，通过下面的结构引发异常。
 * @code
 *  Assert (false, ExcInternalError());
 * @endcode
 * 参见 step-7 和其他几个教程程序中对该结构的使用。
 * 如上所述，一旦对<tt>Assert</tt>的调用失败，程序就会终止。然而，有一种情况我们不想这样做，即当一个C++异常被激活时。发生这种情况的通常情况是，有人通过<tt>AssertThrow</tt>机制抛出一个异常（见下文），在堆栈被解开的同时，导致堆栈框架上面的其他对象被破坏。如果其他对象引用了被销毁的对象，一些析构器会通过<tt>Assert</tt>引发一个异常。如果我们当时中止程序，我们只会看到一个对象被销毁的消息，而这个对象仍然被某个地方引用，但我们永远不会看到触发这个异常的原始异常。(你可以在调试器中通过在函数<tt>__throw</tt>上设置断点来看到它，但你不能从程序本身看到它。)在这种情况下，我们使用一个C++标准库函数来检测另一个活动异常的存在，并且不终止程序，以允许抛出的异常传播到某个可以显示其信息的地方。
 * 由于一个失败的断言导致一连串的其他断言是很常见的，我们只打印第一个消息。如果程序被中止了，那就没有问题。如果不是这样（因为一个C++的异常是有效的），就只显示第一条，并显示一条关于被抑制的后续信息的信息。
 *
 *  <h3>Use of run-time exceptions (dynamic checks)</h3>
 * C++有一种机制来表明发生了一些特殊情况：可以由<tt>throw</tt>语句触发的异常和由<tt>catch</tt>子句捕获的异常，例如见https://en.wikipedia.org/wiki/C%2B%2B#Exception_handling
 * 和 http://www.cplusplus.com/doc/tutorial/exceptions/ 。
 * 在一些基本的层面上，典型的C++异常是一个对象被放置在一些特殊的地方，然后函数通过一个特殊的返回路径退出当前的范围（例如，当前的函数）。
 * 这往往足以说明是什么问题触发了异常，但更多的时候，如果能得到更多的信息就更好了：例如，问题发生在代码的哪一行，或者代码想写进稀疏矩阵的哪个不存在的条目。
 * 因此，deal.II中的动态断言对这种机制进行了一些扩展。
 * 通常情况下，人们会通过以下代码引发一个异常，例如
 * @code
 *  if (!(cond))
 *    throw ExcSomething();
 * @endcode
 * 然后用语句来捕获它
 * @code
 *  try
 *    {
 *      do_something ();
 *    }
 *  catch (std::exception &e)
 *    {
 *      std::cerr << "Exception occurred:" << std::endl
 *                << e.what ()
 *                << std::endl;
 *      do_something_to_receiver ();
 *    }
 * @endcode
 * <tt>std::exception</tt>
 * 是一个标准的<tt>C++</tt>类，为异常提供了基本的功能，比如虚拟函数<tt>what()</tt>，它返回异常本身的一些信息。如果一个异常不能被正确处理，这些信息是很有用的，在这种情况下，应该尽可能地打印出精确的描述。
 * 这里的问题是，要想从<tt>what()</tt>中获得重要而有用的信息，就必须在我们的异常类中重载这个函数，并在异常类中调用带有额外参数的<tt>throw</tt>操作符。首先，重载<tt>what</tt>函数是使用<tt>DeclExceptionN</tt>宏来完成的，但是把正确的信息，也就是上面解释的<tt>Assert</tt>扩展，如果想每次都写下来，需要做一些工作。
 * @code
 *  if (!(cond))
 *    {
 *      ExcSomething e(additional information);
 *      e.set_fields (__FILE__, __LINE__, __PRETTY_FUNCTION__,
 *                    "condition as a string",
 *                    "name of condition as a string");
 *      throw e;
 *    };
 * @endcode
 *
 * 为此，我们发明了宏<tt>AssertThrow</tt>。它所做的工作主要与<tt>Assert</tt>宏相同，但它并不中止程序；相反，它抛出一个异常，如上所示。使用模式是
 * @code
 *  AssertThrow (cond, ExcSomething(additional information));
 * @endcode
 *
 * 要检查的条件被纳入宏中，以允许将被违反的条件作为一个字符串传递。<tt>AssertThrow</tt>宏的扩展不受<tt>DEBUG</tt>预处理器变量的影响。
 *
 *  <h3>Description of the DeclExceptionN macro family</h3>
 * 有一整个系列的<tt>DeclExceptionX</tt>宏，其中<tt>X</tt>将被附加参数的数量所取代（目前是0到5）。
 * 这些宏被用来以如下方式声明异常类。
 * @code
 *  DeclException2 (ExcDomain,
 *                  int,
 *                  int,
 *                  << " i=" << arg1 << ", m=" << arg2);
 * @endcode
 * 第一个参数表示要创建的异常类的名称。
 * 接下来的参数是参数的类型（这里有两种类型，与<tt>DeclExceptionX</tt>中的<tt>X</tt>相对应），最后是输出序列，你可以用它来打印附加信息。
 * 输出序列的语法有点奇怪，但一旦你看到这个宏是如何定义的，就会明白了（同样是示意性的，实际的函数名称和定义可能会随着时间的推移而改变，并有所不同）。
 * @code
 * class name : public ExceptionBase
 * {
 * public:
 *  name (const type1 a1, const type2 a2) : arg1 (a1), arg2(a2)
 *  {}
 *
 *  virtual void print_info (std::ostream &out) const
 *  {
 *    out << "    " outsequence << std::endl;
 *  }
 * private:
 *  type1 arg1;
 *  type2 arg2;
 * };
 * @endcode
 *
 * 如果按照指定的方式声明，你以后就可以按照以下方式使用这个异常类。
 * @code
 *  int i=5;
 *  int m=3;
 *  Assert (i<m, MyExc2(i,m));
 * @endcode
 * 而如果条件失败的话，输出结果将是
 * @code
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * --------------------------------------------------------
 *  An error occurred in line <301> of file <exc-test.cc>.
 *  The violated condition was:
 *    i<m
 *  The name and call sequence of the exception was:
 *    MyExc2(i,m)
 *  Additional Information:
 *    i=5, m=3
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * --------------------------------------------------------
 * @endcode
 *
 * 显然对于<tt>DeclException0(name)</tt>宏，不允许有任何类型，也不允许有任何输出序列。
 *
 *  <h3>A corner case of @p Assert: The @p AssertNothrow macro</h3> @p Assert
 * 宏的默认实现，如上所述，将关于到底出了什么问题的详细信息打印到屏幕上，然后中止程序。终止程序是很有用的，因为它可以轻松地找到出错的地方
 *
 * - 包括我们如何到达那个地方的所有信息
 *
 * - 通过在调试器中运行该程序。
 * 另一方面，在有些情况下，中止程序可能是不可取的，我们需要以一种更优雅的方式退出程序
 *
 * - 即使在这些情况下，我们真的没有什么可以做的，仍然可以产生一个有意义的结果。一个例子是，如果一个deal.II程序是在一个更大的软件框架中运行一个模块。例如，想想这样一种情况：deal.II程序计算的流场与一些优化程序提供的一组输入变量相对应：如果外部的优化器提供了一个负的密度作为输入（一个可能要通过 @p Assert), 检查的条件，那么这显然是没有意义的，流解器不能产生一个有意义的答案；但它应该很好地告诉优化器，而不是直接中止整个过程（优化器和流解器）。
 * 为此，我们可以调用
 * deal_II_exceptions::disable_abort_on_exception() ，将 @p Assert
 * 所做的事情从中止程序转换为与 @p AssertThrow
 * 基本相同，即使用C++ @p throw
 * 机制来引发一个异常。然后这个异常可以在更高层次上被捕获
 *
 * 例如，在位于流解算器之上的优化例程中，然后可以决定它要对这种情况做什么。
 * 这一切都很好，但是C++不允许在类的析构器中或在当前从调用栈中更高的析构器中调用的函数中抛出异常。为此，有一个单独的宏，
 * @p AssertNothrow, ，可以在析构器中使用。它的作用就像 @p
 * Assert 通常做的那样
 *
 * - 特别是，它只检查调试模式下的条件
 *
 * - 但它对 deal_II_exceptions::disable_abort_on_exception(): 的影响是免疫的，它只会中止程序，而不会抛出异常。
 *
 * @author  Wolfgang Bangerth, 1998-2017年
 *
 */


