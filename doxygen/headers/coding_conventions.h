//include/deal.II-translator/A-headers/coding_conventions_0.txt
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


/**
 * @page CodingConventions Coding conventions used throughout deal.II
 * 在整个deal.II中，我们努力保持我们的编程风格和我们提供的界面种类尽可能的一致。为此，我们采用了一套编码惯例，并尽可能地加以遵守。它们有两个部分：风格问题，以及我们称之为
 * "防御性编程
 * "的东西，后者是试图让我们的代码帮助我们发现错误。
 * 在阅读它们的时候，重要的是要记住，风格并不是神赐予的，也不是比任何其他的惯例都要好；它们的目的只是为了使deal.II尽可能的统一。统一性减少了我们产生的错误的数量，因为我们可以，例如，总是假设输入参数在函数调用的输出参数之前。它们也简化了代码的阅读，因为有些东西通过查看代码的写法就已经很清楚了，而不需要去查找某些东西的确切定义。
 * <h3>Notes on deal.II indentation</h3> <p>  deal.II使用
 * <code>clang-format</code>  6.0来规范缩进。Astyle文件提供在
 * @code
 * \${SOURCE_DIR}/.clang-format
 * @endcode
 *  <p>  在提交之前，你应该运行
 * @code
 * clang-format
 *
 * -i <file>
 * @endcode
 * 在你的每一个文件上。这将确保缩进符合本页面中列出的风格指南。
 * 这是很麻烦的。因此，更容易的是，你可以直接运行
 * @code
 * make indent
 * @endcode
 * 在你设置的编译库的任何目录中，对所有最近被修改的源文件进行缩进。如果你想确保所有提交的缩进都是正确的，你可能想设置一个提交后钩子。一种方法是将
 * <code>\${SOURCE_DIR}/contrib/git-hooks/pre-commit</code> 复制到
 * <code>\${SOURCE_DIR}/.git/hooks/pre-commit</code>
 * ，并确保其可执行。
 * 如果你工作的系统安装了不止一个版本的
 * <code>clang-format</code>
 * （或者它不在路径中），你应该将上述 <code>make
 * indent</code> 命令改为
 * @code
 * make DEAL_II_CLANG_FORMAT=/path/to/clang-6.0/clang-format indent
 * @endcode
 * 以指向正确的可执行文件。 </p> <h3>Style issues</h3>
 * <ol>   <li>  %返回某物数量（单元数、自由度等）的函数应以  <code>n_*</code>  开始。例如。    SparsityPatternBase::n_nonzero_elements().</li>
 * <li>  %设置位或标志的函数应以  <code>set_*</code>
 * 开始；清除位或标志的函数应以  <code>clear_*</code>  命名。
 * 例子。   CellAccessor::set_refine_flag().</li>
 * <li>
 * 应使用传统的逻辑运算符，而不是它们的英文等价物（即，使用
 * <code>&&</code>, <code>||</code>, and <code>!</code> 而不是
 * <code>and</code>, <code>or</code>, and <code>not</code>  ）。
 * <li>
 * 在实现文件中，在每个函数之后，希望有三个空行，以使可读性更好。一个空行出现在函数中，用于对代码块进行分组，因为两个空行不足以明显地区分出代码属于两个不同的函数。
 * </li>
 * <li>
 * 每当一个整数变量只能承担非负值，它就被标记为无符号。这同样适用于只能返回正值或零值的函数。例子。
 * Triangulation::n_active_cells().</li>
 * <li>
 * 每当一个函数的参数不会被改变，它就应该被标记为const，即使是通过值传递。一般来说，我们将输入参数标记为const。这有助于作为一个额外的文档工具来澄清参数的意图（输入、输出或两者），并让编译器在这样的参数被改变时发出警告，这往往是不由自主的或糟糕的风格。
 * </li>
 * <li>
 * 每当一个函数不改变嵌入类/对象的任何成员变量时，它就应该被标记为const。
 * </li>
 * <li>
 * %函数和变量名称不能只由一个或两个字母组成，除非该变量是一个纯粹的计数索引。
 * </li>
 * <li>  类型别名（  <code>using</code>  -declarations）优于
 * <code>typedef</code>  -declarations。 </li>
 * <li>
 * 使用GeometryInfo中的几何信息来获取每个单元的面数、每个单元的子数、与面3相邻的子单元的子指数等，而不是像
 * <code>2*dim</code>  、  <code>(1@<@<dim)</code>  和  <code>{0,3}</code>
 * 那样直接写进代码。这减少了出错的可能性，提高了代码的可读性。
 * </li>
 * <li>  类声明的布局如下：首先是公共函数块，从构造函数开始，然后是解构函数。如果有公共成员变量，这些变量必须出现在构造函数之前。只有当公有变量是常量（特别是静态和常量）或不可避免时，才应使用公有变量。    <br>  在公共成员之后，应列出受保护成员，最后是私有成员。顺序如上：先是变量，然后是函数。    <br>  异常应在公共部分的末尾声明，然后再开始非公共部分。    <br>  对于既不是 <code>static const</code> nor <code>static constexpr</code> 的成员变量，我们不使用C++11风格的类成员初始化；即，代替
 * @code
 * class Foo
 * {
 *  int a = 42;
 *  intb = nullptr;
 * };
 * @endcode
 * 写
 * @code
 * class Foo
 * {
 *  Foo();
 *
 *  int a;
 *  intb;
 * };
 *
 *
 *
 *
 *
 * inline Foo::Foo()
 * : a(42)
 * , b(nullptr)
 * {}
 * @endcode
 * </li>
 * <li>
 * 如果一个函数有输入和输出参数，通常输入参数应在输出参数之前，除非有充分的理由改变这个顺序。最常见的原因是带有默认值的输入参数的尾部。
 * <li>
 * 异常用于%内部参数检查和通过Assert宏进行一致性检查。像C++语言那样的异常处理（
 * <code>try/throw/catch</code>
 * ，并使用AssertThrow宏）用于处理运行时的错误（如I/O故障），在任何情况下都必须开启，而不仅仅是在调试模式下。
 * </li>
 * <li>
 * 有时通过使用几个非成员函数来实现一个类是有意义的，这些非成员函数不是公共接口的一部分，而且只在当前源文件中被调用。这样的自由函数应该放在一个内部命名空间中，其结构如下。
 * @code
 * namespace internal
 * {
 *  namespace ClassNameImplementation
 *  {
 *    // free functions go here
 *  }
 * }
 * @endcode
 * 其中 <code>ClassName</code> 是调用类的名称。
 * <li>
 * 类、命名空间和类型一般使用大写字母来表示词的开头（例如TriaIterator）&mdash；有时也称为<a
 * href="http://en.wikipedia.org/wiki/Camel_case"><i>camel case</i><i>camel
 * case</i></a>&mdash；而函数和变量使用小写字母和下划线来分隔单词。
 * 唯一的例外是Triangulation和DoFHandler中的迭代器别名（命名为cell_iterator、active_line_iterator等），以使其与标准库容器类的联系清晰。
 * </li>
 * <li>
 * 对于有多个模板参数的类，维度通常放在数据类型指定器之前，即我们使用Point<dim,number>而不是Point<number,dim>。
 * <li>
 * 在deal.II中，有几个地方我们在头文件中使用正向声明。这样做的原因是，当我们只需要将某个类型标记为某个函数的参数时，不使用头文件，希望可以提高编译速度。在deal.II中使用的惯例是，如果我们需要的只是一个类型名称，那么这个类型可以在我们需要的头文件中向前声明；如果一个函数（或成员函数）可以返回一个值，那么这个值的类型声明应该可以得到（通过包括必要的头文件）。例如，
 * <code>deal.II/dofs/dof_handler.h</code> 包括
 * <code>deal.II/dofs/dof_accessor.h</code> ，这样就可以写出类似
 * <code>dof_handler.begin_active()->is_active()</code>
 * 的东西，而不需要明确包括声明 <code>begin_active()</code>
 * 所返回对象类型的头文件。
 * <li>  每个类都必须有至少200页的文档；-) </li>
 * </ol>
 *
 *  <h3>Instantiation of templated functions/classes</h3> <p>
 * deal.II中的大多数类和函数是模板化的。这就带来了一个问题：如果有的话，这些对象是如何以及在哪里被实例化的。在整个deal.II中，我们采用了以下惯例。
 * </p>
 * <ol>
 * <li>
 * 如果我们可以列举出所有可能的模板参数（例如，维度只能是1、2或3），那么一个函数模板就会进入
 * <code>.cc</code>
 * 文件，我们明确地将所有可能性实例化。用户不会有任何需要看到这些函数模板，因为他们无论如何都不想为其他模板参数实例化这些函数。
 * </li>
 * <li>
 * 如果我们不能列举所有可能的模板参数（例如，vectortypes
 *
 * 因为用户可能想定义他们自己的向量种类），但至少知道一些常见的使用情况，那么该函数将被放入一个
 * <code>.templates.h</code> file. We \#include it into the <code>.cc</code>
 * 文件中，并为所有常见的参数实例化函数。对于几乎所有的用户来说，这样做就可以了
 *
 * - 他们只使用我们已经实例化的（向量、矩阵......）类型，对他们来说， <code>.templates.h</code> 文件不会有任何意义。这也不会减慢他们的编译速度，因为他们看到的东西都会包括 <code>.templates.h</code> 文件。但是定义了自己的（向量、矩阵......）类型的用户可以通过包括 <code>.templates.h</code> 文件，用自己的用户定义的类型实例化模板函数。
 * <li>
 * 最后，如果我们不能事先假定模板参数将采取哪些值（例如，任何从Subscriptor派生的类都可以作为参数），函数的定义将在头文件的底部提供声明。这些定义应该用<code>\#ifndefDOXYGEN
 * ...\#以防止Doxygen发现它们。 </li>
 * </ol>
 * <p>  对于前两种情况，实例化指令被定义在
 * <code>.inst.in</code>
 * 文件中。它们由一个叫做expand_instantiations的二进制文件处理（由
 * <code>cmake/scripts/expand_instantiations.cc</code>
 * 构建），参数根据你的配置通过cmake动态定义（见构建目录中的
 * <code>cmake/config/template-arguments.in</code> ）。正是这些
 * <code>.inst</code> 文件最终被包含在相应的 <code>.cc</code>
 * 文件中。   </p>
 *
 *  <h3>Defensive programming</h3> <p>
 * 防御性编程是我们经常使用的一个术语，当我们在谈论写代码的时候，要考虑到错误会发生。在这里，错误有两种方式：第一，我自己在写函数时可能犯错；第二，别人在调用我的函数时可能犯错。不管是哪种情况，我都希望我的代码写得(i)尽可能不出错，(ii)编译器已经可以发现一些错误，(iii)剩下的错误相对容易发现，比如说因为程序中止了。因此，防御性编程是一组使这些目标更有可能实现的策略。
 * </p>
 * <p>  随着时间的推移，我们已经学会了一些这方面的技术，其中一些我们在此列出。 <ol>   <li>  <i>Assert preconditions on parameters:</i> 人们总是用错误或无意义的参数调用函数。作为典型的例子，考虑一个微不足道的向量加法的实现。
 * @code
 *  Vector &
 *  operator+=(Vector       &lhs,
 *             const Vector &rhs)
 *  {
 *    for (unsigned int i=0; i<lhs.size(); ++i)
 *      lhs(i) += rhs(i);
 *    return lhs;
 *  }
 * @endcode
 * 虽然正确，但如果两个向量的大小不一样，这个函数就会陷入困境。你认为用不同大小的向量来调用这个函数是愚蠢的？是的，当然是这样的。但这种情况经常发生：人们忘记重新初始化一个向量，或者它在不同的函数中被重置，等等。这种情况时有发生。所以，如果你遇到这种不幸运的情况，可能要花很长时间才能弄清楚发生了什么，因为你很可能只是读取了未初始化的内存，或者也许你正在向
 * <code>lhs</code> 向量实际上并不拥有的内存写入。
 * 这两种情况都不会导致程序的立即终止，但你可能会在以后的时间里出现随机错误。如果程序只是在这里立即停止，那就容易多了。下面的实现正是这样做的。
 * @code
 *  Vector &
 *  operator+=(Vector       &lhs,
 *             const Vector &rhs)
 *  {
 *    Assert (lhs.size() == rhs.size(),
 *            ExcDimensionMismatch(lhs.size(), rhs.size());
 *    for (unsigned int i=0; i<lhs.size(); ++i)
 *      lhs(i) += rhs(i);
 *    return lhs;
 *  }
 * @endcode
 * <code>Assert</code>
 * 宏确保条件在运行时为真，否则会打印一个包含由第二个参数编码的信息的字符串，并中止程序。这样，当你写一个新的程序恰好调用这个函数时，你会马上知道你的错误，并有机会修复它，而不需要认真调试什么。
 * <p>
 * 作为一般准则，每当你实现一个新的函数时，要考虑<i>preconditions</i>上的参数，即该函数对每一个参数或它们的组合期望是什么。
 * 然后为所有这些先决条件写断言。在某些情况下，这可能是半打断言，但请记住，每个断言都是一个潜在的错误，已经通过琐碎的手段发现。
 * <p>
 * 最后，让我们说说断言当然是昂贵的：当你把一个程序与库的调试版本链接时，它们可能会使程序慢3或5倍。但是如果你考虑到你的<i>overall</i>开发时间，快速发现错误的能力可能远远超过你等待程序完成的时间。此外，在优化模式下，对Assert宏的调用被从程序中移除（你可能只在知道在调试模式下一切运行正常后才使用。优化后的库比调试后的库快3-5倍，但代价是更难发现错误。
 * </li>
 * <li>  <i>Assert postconditions:</i>
 * 如果一个函数计算出一些非微不足道的东西，那么代码中可能有一个错误。为了找到这些，可以使用后置条件：就像你对输入参数的有用值有一定的了解一样，你对可能的返回值也有一定的了解。例如，一个计算向量规范的函数希望规范是正的。你可以把它写成这样。
 * @code
 *  double norm(const Vector &v)
 *  {
 *    double s = 0;
 *    for (unsigned int i=0; i<v.size(); ++i)
 *      s += v(i) v(i);
 *
 *    Assert (s >= 0, ExcInternalError());
 *    return std::sqrt(s);
 *  }
 * @endcode
 * 这个函数太简单了，无法真正证明这个断言的正确性，但是想象一下计算的长度，你就可以看到这个断言是如何帮助你确保（或<i>hedge</i>）自己不犯错误的。请注意，人们可以争辩说，一旦我们运行了若干次程序，发现条件从未被触发，就应该删除这个断言。但最好还是把它留在原处：它为未来（和读者）编码了你对该函数的了解；如果有人出现，用更有效的算法取代了该函数的实现，断言可以帮助确保该函数继续做它应该做的事。
 * </li>
 * <li>  <i>Assert internal states:</i>
 * 类似地，如果你有一个复杂的算法，使用断言来确保你对正在发生的事情的心理模型与确实的事实相符。例如，假设你正在编写一个函数，以确保网格尺寸不会在本地发生太大的变化。你最终可能会得到如下的代码。
 * @code
 *  for (const auto &cell = triangulation.active_cell_iterators())
 *    for (unsigned int face=0; ...)
 *      {
 *        if (something)
 *          { ... }
 *        else
 *          {
 *            // we have a cell whose neighbor must
 *            // be at the boundary if we got here
 *          }
 *      }
 * @endcode
 * 导致我们进入else-branch的条件可能很复杂，虽然我们认为我们到达这里的唯一可能性是邻居在边界上，这可能是真的，但在我们的实现中可能存在一个bug。也可能是我们的思维出现了错误，或者有人在同一个函数中改变了上面的代码而忘记了这里的问题，或者在库中一个完全不同的位置的改变使得这个假设站不住脚。在所有这些情况下，我们的断言的明确陈述可以确保这些问题被轻易发现。
 * </li>
 * <li>  <i>Initialize variables at the point of their declaration if they
 * live on the stack:</i>
 * 传统的C语言要求在函数的开头声明变量，即使它们只在下面进一步使用。这就导致了我们可以想象在1d代码中出现这样的代码。
 * @code
 *  template <int dim>
 *  void foo ()
 *  {
 *    Point<dim> cell_center;
 *    ... // something lengthy and complicated
 *    for (const auto &cell = dof_handler.active_cell_iterators())
 *      {
 *        cell_center = (cell->vertex(0) + cell->vertex(1)) / 2;
 *        ...
 *      }
 *    ...
 *  }
 * @endcode
 * 问题是，如果声明和初始化之间的代码又长又复杂，你不可能在一页纸上查到一个变量的类型是什么，它的值可能是什么。事实上，甚至可能不太清楚该变量是否被用于初始化，或者它是否被意外地留在了未初始化状态。
 * <p>  一个更好的方法是这样做的。
 * @code
 *  template <int dim>
 *  void foo ()
 *  {
 *    ... // something lengthy and complicated
 *    for (const auto &cell = dof_handler.active_cell_iterators())
 *      {
 *        Point<dim> cell_center = (cell->vertex(0) + cell->vertex(1)) / 2;
 *        ...
 *      }
 *    ...
 *  }
 * @endcode
 * 这使得变量的类型更加清晰，而且事实上它只在初始化时使用。此外，如果有人想阅读代码，看看这个变量实际上在做什么，在最内层的作用域中声明和初始化它，使这项任务更容易：我们不必在声明之外向上寻找它，也不必在当前作用域的末端向下寻找，因为这是变量死亡的地方。
 * <p>
 * 作为最后的说明，很明显，你只能对完全生活在堆栈中的变量做这种事情，而不需要在堆上分配内存。在deal.II中，这只适用于像
 * <code>int, double, char</code>
 * 等内置类型，以及点和张量类。其他的东西都有像
 * <code>std::vector</code>
 * 这样的成员变量，这需要内存分配&mdash；你不想在循环内声明这些，至少在循环被频繁遍历的情况下。
 * </li>
 * <li>  <i>Make variables const:</i>
 * 继续上面的例子，注意在大多数情况下，我们不会再改变如此初始化的变量。换句话说，如果是这种情况，我们不妨把事情写成如下。
 * @code
 *  template <int dim>
 *  void foo ()
 *  {
 *    ... // something lengthy and complicated
 *    for (const auto &cell = dof_handler.active_cell_iterators())
 *      {
 *        const Point<dim> cell_center = (cell->vertex(0) +
 *                                        cell->vertex(1)) / 2;
 *        ...
 *      }
 *    ...
 *  }
 * @endcode
 * 通过标记变量为常量，我们可以确保我们不会意外地改变它。例如，编译器可以捕获这样的代码。
 * @code
 *      if (cell_center[0] = 0)
 *        ...
 * @endcode
 * 这很可能是指 <code>==</code>
 * 而不是一个赋值。通过将该变量标记为常量，编译器就会告诉我们这个错误。也许同样重要的是，代码的人类读者不需要再往下看，变量的值是否真的在声明和使用之间的某个地方被改变了&mdash；如果它被标记为const，就不可能了。
 * </li>
 * <li>  <i>Make input arguments of functions const:</i>
 * 对于函数参数，本质上也是如此。如果你不打算改变一个变量（这通常是输入参数的情况），那么就把它标记为常量。例如，下面这个函数应该把它的参数作为一个常数。
 * @code
 *   template <int dim>
 *   typename Triangulation<dim>::cell_iterator
 *   CellAccessor<dim>::child(const unsigned int child_no)
 *   {
 *     ...
 *     return something;
 *   }
 * @endcode
 * 这里，用户调用 <code>cell-@>child(3)</code>  ，例如。真的没有理由这个函数要改变 <code>child_no</code> 这个参数的值&mdash；所以把它标记为常量：这既可以帮助代码的读者理解这是一个函数的输入参数，我们不需要在下面搜索它是否被改变过，又可以帮助编译器在我们不小心改变了这个值时帮助我们找到bug。 </ol>
 *
 *
 */


