//include/deal.II-translator/base/parameter_handler_0.txt
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

#ifndef dealii_parameter_handler_h
#define dealii_parameter_handler_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/subscriptor.h>

#include <boost/archive/basic_archive.hpp>
#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/property_tree/ptree_serialization.hpp>
#include <boost/serialization/split_member.hpp>

#include <map>
#include <memory>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// forward declarations for interfaces and friendship
#ifndef DOXYGEN
class LogStream;
class MultipleParameterLoop;
#endif

/**
 * ParameterHandler类为输入文件提供了一个标准接口，该文件在运行时提供了程序参数，如时间步长、几何尺寸、右手边等。程序的输入是以文件、流或内存中的字符串的形式给出的，使用的文本包括
 * @code
 *   set Time step size = 0.3
 *   set Geometry       = [0,1]x[0,3]
 * @endcode
 * 输入可以被排序为分节树，以便给输入一个逻辑结构，输入文件可以包括其他文件。
 * ParameterHandler类在  step-29  ,  step-33  , 和  step-34  中讨论。
 * <h3>Declaring entries</h3>
 * 为了使用ParameterHandler对象的设施，首先必须知道输入文件可能包含或不包含的不同条目。这可以通过以下方式完成。
 * @code
 *   ...
 *   ParameterHandler prm;
 *   prm.declare_entry ("Time step size",
 *                      "0.2",
 *                      Patterns::Double(),
 *                      "Some documentation");
 *   prm.declare_entry ("Geometry",
 *                      "[0,1]x[0,1]",
 *                      Patterns::Anything());
 *   ...
 * @endcode
 * 每个条目都要用函数
 * declare_entry()来声明。第一个参数是条目的名称（简而言之：条目）。第二个参数是在输入文件中没有指定条目时要采取的默认答案。第三个参数是一个正则表达式，输入（和默认答案）必须与之匹配。
 * Patterns中定义了几个这样的正则表达式。这个参数可以省略，在这种情况下，它将默认为
 * Patterns::Anything,
 * ，即一个匹配每个输入字符串的模式。第四个参数可用于记录条目的意图或预期格式；在使用print_parameters()函数编写ParameterHandler对象的所有条目时，其值将作为注释打印出来，以便更容易理解参数文件。它也可以被省略，在这种情况下，将不会打印这样的文件。
 * 条目可以位于子节中，形成一种输入树。例如，线性求解程序的输入参数应该被分类在一个名为<tt>线性求解器</tt>或任何其他合适名称的小节中。这可以通过以下方式实现。
 *
 * @code
 * ...
 * LinEq eq;
 * eq.declare_parameters (prm);
 * ...
 *
 * void LinEq::declare_parameters (ParameterHandler &prm)
 * {
 * prm.enter_subsection("Linear solver");
 * {
 *   prm.declare_entry ("Solver",
 *                      "CG",
 *                      Patterns::Selection("CG|GMRES|GaussElim"),
 *                      "Name of a linear solver for the inner iteration");
 *   prm.declare_entry ("Maximum number of iterations", "20",
 *                      ParameterHandler::RegularExpressions::Integer());
 *   ...
 * }
 * prm.leave_subsection ();
 * }
 * @endcode
 *
 * 分节可以是嵌套的。例如，一个非线性求解器可能有一个线性求解器作为成员对象。那么函数调用树会是这样的（如果类<tt>NonLinEq</tt>有一个成员变量<tt>eq</tt>类型<tt>LinEq</tt>）。
 * @code
 * void NonLinEq::declare_parameters (ParameterHandler &prm)
 * {
 * prm.enter_subsection ("Nonlinear solver");
 * {
 *   prm.declare_entry ("Nonlinear method",
 *                      "Newton-Raphson",
 *                      ParameterHandler::RegularExpressions::Anything());
 *   eq.declare_parameters (prm);
 * }
 * prm.leave_subsection ();
 * }
 * @endcode
 *
 * 对于声明不同条目的类成员函数，我们建议使用通用名称<tt>declare_parameters</tt>。在一般情况下，这个方法可以是<tt>静态的</tt>，因为这些条目将不依赖于任何先前的知识。那些条目在逻辑上应该被分组为子段的类应该自己声明这些子段。如果一个类有两个或更多相同类型的成员变量，它们都应该有自己的参数，这个父类的方法<tt>declare_parameters</tt>负责将它们归入不同的子段。
 *
 * @code
 * void NonLinEq::declare_parameters (ParameterHandler &prm)
 * {
 * prm.enter_subsection ("Nonlinear solver");
 * {
 *   prm.enter_subsection ("Linear solver 1");
 *   {
 *     eq1.declare_parameters (prm);
 *   }
 *   prm.leave_subsection ();
 *
 *   prm.enter_subsection ("Linear solver 2");
 *   {
 *     eq2.declare_parameters (prm);
 *   }
 *   prm.leave_subsection ();
 * }
 * prm.leave_subsection ();
 * }
 * @endcode
 *
 *  <h3>Input files and special characters</h3>
 * 对于上面的第一个例子，输入文件将看起来像下面这样。
 * @code
 *   ...
 *   subsection Nonlinear solver
 *     set Nonlinear method = Gradient
 *     # this is a comment
 *     subsection Linear solver
 *       set Solver                       = CG
 *       set Maximum number of iterations = 30
 *     end
 *   end
 *   ...                       # other stuff
 * @endcode
 * <tt>subsection</tt>,
 * <tt>set</tt>和<tt>end</tt>这些词可以用小写字母或大写字母书写。前面和后面的空白被删除，多个空白被浓缩成一个。由于后者也适用于条目的名称，如果在声明中使用了多个空白，条目名称将无法被识别。
 * 在条目名称和值中，不允许使用以下字符。<tt>\#</tt>,
 * <tt>{</tt>, <tt>}</tt>, <tt>|</tt>.
 * 它们的使用是为MultipleParameterLoop类保留的。
 * 以\#开头的注释被跳过。
 * 通过字符<tt>\\</tt>允许续行，该字符必须是该行的最后一个字符（除了空白处，空白处被忽略）。当一行是续行时（即前一行以<tt>\</tt>结束），那么，与<tt>C</tt>预处理器的默认行为不同，该行开头的所有空白都被忽略了。
 * 我们建议使用以下方案来命名条目：第一个词用大写字母开头，以后用小写字母。这同样适用于<tt>=</tt>符号右边的可能条目值。
 *
 *  <h3>Including other input files</h3>
 * 一个输入文件可以包括其他包含文件，使用的语法是
 * @code
 *   ...
 *   include some_other_file.prm
 *   ...
 * @endcode
 * 这样引用的文件是相对于当前目录搜索的（不是相对于包含参数文件所在的目录搜索的，因为这不是所有三个版本的parse_input()函数都知道的）。
 *
 *  <h3>Reading data from input sources</h3>
 * 为了读取输入，有三种可能性：从一个 <tt>std::istream</tt>
 * 对象中读取，从一个给定名称的文件中读取，以及从内存中的字符串中读取，其中各行由<tt>
 * @\n</tt> 字符分隔。这些可能性的使用情况如下。
 * @code
 *   ParameterHandler prm;
 *   ...
 *   // declaration of entries
 *   ...
 *   prm.parse_input (std::cin); // read input from standard in,
 *   // or
 *   prm.parse_input ("simulation.prm");
 *   // or
 *   charin = "set Time step size = 0.3 \n ...";
 *   prm.parse_input_from_string (in);
 *   ...
 * @endcode
 * 你可以连续使用几个输入源。被改变一次以上的条目在每次使用时都会被覆盖。
 * 你不应该在使用parse_input()后，试图用 declare_entry()和
 * enter_subsection()来声明条目，并使用尚不清楚的分段名称。这种情况下的结果是不明确的。
 * 如果在读取输入时发生错误，错误信息将被写入
 * <tt>std::cerr</tt>  ，并且阅读器函数的返回值为
 * <code>false</code>
 * 。这与deal.II中几乎所有的其他函数相反，如果发生错误，通常会抛出一个异常；这种行为上的差异是一个遗留问题，即这个类在deal.II之前就已经存在，而且之前是为一个不同的项目编写的。
 *
 *  <h3>Using the %ParameterHandler Graphical User Interface</h3>
 * 除了使用上面显示的手写输入文件外，还有一种方法是使用伴随这个类的图形用户界面（GUI）。
 * 更多细节见<a href="https://github.com/dealii/parameter_gui">the
 * parameter_gui github repository</a>。 <h3>Getting entry values out of a
 * %ParameterHandler object</h3>
 * 每个类通过调用get()成员函数，从ParameterHandler对象中获取数据，比如这样。
 * @code
 *   void NonLinEq::get_parameters (ParameterHandler &prm)
 *   {
 *     prm.enter_subsection ("Nonlinear solver");
 *     std::string method = prm.get ("Nonlinear method");
 *     eq.get_parameters (prm);
 *     prm.leave_subsection ();
 *   }
 * @endcode
 * get()返回给定条目的值。如果输入源中没有指定该条目，则返回默认值。你必须完全按照你在声明子段时的做法来输入和离开子段。你可以选择遍历子树的顺序。
 * 可以通过向get()提供一个代表获取值的路径的字符串向量来避免对enter_subsection()和leave_subsection()的调用。例如，以下两个版本的get_parameters()将产生相同的结果。
 * @code
 *   void NonLinEq::get_parameters (ParameterHandler &prm)
 *   {
 *     prm.enter_subsection ("Equation 1 Settings");
 *     prm.enter_subsection ("Linear solver");
 *     solver_ = prm.get ("Solver");
 *     prm.leave_subsection ();
 *     prm.leave_subsection ();
 *   }
 * @endcode
 *
 *
 * @code
 *   void NonLinEq::get_parameters (const ParameterHandler &prm)
 *   {
 *     std::vector<std::string> path =
 *       {"Equation 1 Settings", "Linear solver"};
 *     solver_ = prm.get (path, "Solver");
 *   }
 * @endcode
 *
 * 后一种方法允许ParameterHandler的引用是 @p const.  。
 * 保证只返回与给定的正则表达式相匹配的条目，也就是说，不符合正则表达式的输入条目值不会被存储。
 * 你可以使用get()来检索文本形式的参数，使用get_integer()来获得一个整数，或者使用get_double()来获得一个双数。你也可以使用get_bool()。如果字符串不能被转换为整数、双数或布尔，它将导致一个内部错误。不过，如果你为这个条目正确地指定了正则表达式，这就不应该发生；你不应该试图从一个没有设置相应正则表达式的条目中得到一个整数或双数。内部错误是通过Assert()宏程序提出的，它只在调试模式下工作。
 * 如果你想打印出所有用户可选择的特征，请使用print_parameters()函数。一般来说，在日志文件的开头打印所有的参数是一个好主意，因为这样输入和输出都在一个文件中，这使得以后的匹配更加容易。此外，该函数也会打印那些在输入文件中没有被修改，因而被设置为默认值的条目；由于默认值在程序开发过程中可能会改变，所以你无法知道输入文件中没有指定的参数值。
 *
 *
 * <h3>Adding Actions to Parameters</h3>
 * 在读取一个参数值的时候，通常会有一些事情发生，这很方便。这可以是检查它是否有效
 *
 * 例如，在参数文件中列出的一个文件是否存在
 *
 * 或者启动其他的响应，比如在ParameterHandler之外设置一个变量（如下图所示的例子）。在几乎所有的情况下，这个
 * "动作
 * "也可以在通过parse_input()读取所有参数后启动，但有时<i>convenient</i>会马上进行。
 * add_action()函数可以在通过declare_entry()声明一个参数后调用，这就方便了。"动作
 * "实质上是指向将为具有相关动作的参数调用的函数的指针。这些函数把参数的值作为参数，然后可以对它做任何事情
 *
 * 例如，将其保存在ParameterHandler对象之外的某个地方。(具体什么时候调用动作，在add_action()函数的文档里有描述)。当然，在C++中，我们通常不会传递一个函数的地址，但是一个动作可以是一个类似于函数的对象（接受一个字符串作为参数），这个对象是通过调用诸如<a
 * href="http://en.cppreference.com/w/cpp/language/lambda">lambda
 * function</a>的形式产生的。
 *
 * @code
 * [] (const std::string &value) { ... do something with the value ... }
 * @endcode
 * 并且附在一个特定的参数上。
 * 这种行为的一个典型例子如下：假设你有一个程序，它声明了一个参数，表示它要运行的迭代次数，例如
 *
 * @code
 * class MyAlgorithm
 * {
 *    public:
 *      void run ();
 *    private:
 *      unsigned int n_iterations;
 * };
 * @endcode
 * 那么人们可以使用 @p run()
 * 中的代码片段从参数文件中获得这个参数，如下所示。
 *
 * @code
 * void MyAlgorithm::run ()
 * {
 *   ParameterHandler prm;
 *   prm.declare_entry ("Number of iterations",  // name of parameter
 *                      "10",                    // default value
 *                      Patterns::Integer(1,100),// allowed values: 1...100
 *                      "The number of ...");    // some documentation
 *
 *   // next read the parameter from an input file...
 *   prm.parse_input ("my_algorithm.prm");
 *
 *   // ...and finally get the value for use in the program:
 *   n_iterations = prm.get_integer ("Number of iterations");
 *
 *   ... actual code doing something useful follows here...
 * @endcode
 *
 * 这个两步过程
 *
 * - 首先声明参数，然后再读取它
 *
 * - 有点麻烦，因为必须先声明<i>all</i>参数，以后再从ParameterHandler对象中检索参数。在大型程序中，这两件事也经常发生在不同的函数中。
 * 为了避免这种情况，如果我们能把声明和检索都放在同一个地方就好了。这可以通过动作来实现，然后函数会看起来像这样。
 *
 * @code
 * void MyAlgorithm::run ()
 * {
 *   ParameterHandler prm;
 *   prm.declare_entry ("Number of iterations",  // name of parameter
 *                      "10",                    // default value
 *                      Patterns::Integer(1,100),// allowed values: 1...100
 *                      "The number of ...");    // some documentation
 *   prm.add_action ("Number of iterations",
 *                   [&](const std::string &value)
 *                   {
 *                     this->n_iterations = Utilities::string_to_int(value);
 *                   });
 *
 *   // next read the parameter from an input file...
 *   prm.parse_input ("my_algorithm.prm");
 *
 *   ... actual code doing something useful follows here...
 * @endcode
 * 这里，动作由一个lambda函数组成，它将这个参数的值作为一个字符串，然后将其转换为一个整数，存储在它所属的变量中。这个动作是在调用
 * <code>prm.parse_input()</code>
 * 时执行的，所以现在不再需要在以后提取参数的值。此外，设置成员变量的代码就在实际声明参数的地方旁边，所以我们不再需要在代码库中设置两个单独的部分来处理输入参数。
 * 当然，我们有可能执行比上面所示的设置成员变量更多的动作，尽管这只是一个典型的案例。
 *
 *  <h3>Style guide for data retrieval</h3>
 * 我们建议每个从ParameterHandler对象中获取数据的类都提供一个名为<tt>get_parameters</tt>的函数。这应该被声明为<tt>虚拟的</tt>。派生类中的<tt>get_parameters</tt>函数应该调用
 * <tt>BaseClass::get_parameters</tt> 函数。
 *
 *  <h3>Experience with large parameter lists</h3>
 * 经验表明，在定义大量参数的程序中（例如超过50个），定义一个额外的类来保存这些参数是有利的。这个类更像一个C风格的结构，有大量的变量，通常是公共变量。然后它至少有两个函数，用来声明和解析参数。在主程序中，主类有一个这个参数类的对象，并将参数的声明和解析委托给这个对象。
 * 这种方法的好处是，你可以不把技术细节（声明和解析）放在主类之外，另外也不用用几十个或更多的变量来表示参数，使你的主类变得混乱。
 *
 *
 * <h3>Worked Example</h3> 这就是代码。
 * @code
 * #include <deal.II/base/parameter_handler.h>
 *
 * #include <iostream>
 * #include <string>
 *
 * using namespace dealii;
 * class LinearEquation
 * {
 * public:
 *   static void declare_parameters (ParameterHandler &prm);
 *   void get_parameters (ParameterHandler &prm);
 * private:
 *   std::string method;
 *   int         max_iterations;
 * };
 *
 *
 *
 *
 *
 * class Problem
 * {
 * private:
 *   LinearEquation eq1, eq2;
 *   std::string matrix1, matrix2;
 *   std::string outfile;
 * public:
 *   static void declare_parameters (ParameterHandler &prm);
 *   void get_parameters (ParameterHandler &prm);
 *
 *   void do_something ();
 * };
 *
 *
 *
 *
 *
 * void LinearEquation::declare_parameters (ParameterHandler &prm)
 * {
 *   // declare parameters for the linear solver in a subsection
 *   prm.enter_subsection ("Linear solver");
 *   {
 *     prm.declare_entry ("Solver",
 *                        "CG",
 *                        Patterns::Selection("CG|BiCGStab|GMRES"),
 *                        "Name of a linear solver for the inner iteration");
 *     prm.declare_entry ("Maximum number of iterations",
 *                        "20",
 *                        Patterns::Integer());
 *   }
 *   prm.leave_subsection ();
 * }
 *
 *
 *
 *
 *
 * void LinearEquation::get_parameters (ParameterHandler &prm)
 * {
 *   prm.enter_subsection ("Linear solver");
 *   {
 *     method         = prm.get ("Solver");
 *     max_iterations = prm.get_integer ("Maximum number of iterations");
 *   }
 *   prm.leave_subsection ();
 *   std::cout << "  LinearEquation: method=" << method
 *             << ", max_iterations=" << max_iterations
 *             << std::endl;
 * }
 *
 *
 *
 *
 *
 * void Problem::declare_parameters (ParameterHandler &prm)
 * {
 *   // first some global parameter entries
 *   prm.declare_entry (
 *     "Output file",
 *     "out",
 *     Patterns::Anything(),
 *     "Name of the output file, either relative or absolute");
 *   prm.declare_entry ("Equation 1", "Laplace",
 *                      Patterns::Anything(),
 *                      "String identifying the equation we want to solve");
 *   prm.declare_entry ("Equation 2",
 *                      "Elasticity",
 *                      Patterns::Anything());
 *
 *   // declare parameters for the first equation
 *   prm.enter_subsection ("Equation 1 Settings");
 *   {
 *     prm.declare_entry ("Matrix type",
 *                        "Sparse",
 *                        Patterns::Selection("Full|Sparse|Diagonal"),
 *                        "Type of the matrix to be used, either full, "
 *                        "sparse, or diagonal");
 *     LinearEquation::declare_parameters (prm);  // for eq1
 *   }
 *   prm.leave_subsection ();
 *
 *   // declare parameters for the second equation
 *   prm.enter_subsection ("Equation 2 Settings");
 *   {
 *     prm.declare_entry ("Matrix type",
 *                        "Sparse",
 *                        Patterns::Selection("Full|Sparse|Diagonal"));
 *     LinearEquation::declare_parameters (prm);  // for eq2
 *   }
 *   prm.leave_subsection ();
 * }
 *
 *
 *
 *
 *
 * void Problem::get_parameters (ParameterHandler &prm)
 * {
 *   // entries of the problem class
 *   outfile = prm.get ("Output file");
 *   std::string equation1 = prm.get ("Equation 1"),
 *               equation2 = prm.get ("Equation 2");
 *
 *   // get parameters for the first equation
 *   prm.enter_subsection ("Equation 1 Settings");
 *   {
 *     matrix1 = prm.get ("Matrix type");
 *     eq1.get_parameters (prm); // for eq1
 *   }
 *   prm.leave_subsection ();
 *
 *   // get parameters for the second equation
 *   prm.enter_subsection ("Equation 2 Settings");
 *   {
 *     matrix2 = prm.get ("Matrix type");
 *     eq2.get_parameters (prm); // for eq2
 *   }
 *   prm.leave_subsection ();
 *   std::cout
 *     << "  Problem: outfile=" << outfile << '\n'
 *     << "           eq1="     << equation1 << ", eq2=" << equation2 << '\n'
 *     << "           matrix1=" << matrix1 << ", matrix2=" << matrix2
 *     << std::endl;
 * }
 *
 *
 *
 *
 *
 * void Problem::do_something ()
 * {
 *   // While this example does nothing here, at this point in the program
 *   // all of the parameters are known so we can start doing computations.
 * }
 *
 *
 *
 *
 *
 * int main ()
 * {
 *   ParameterHandler prm;
 *   Problem p;
 *   p.declare_parameters (prm);
 *   // read input from "prmtest.prm"; giving argv[1] would also be a
 *   // good idea
 *   prm.parse_input ("prmtest.prm");
 *   // print parameters to std::cout as ASCII text
 *   std::cout << "\n\n";
 *   prm.print_parameters (std::cout, ParameterHandler::Text);
 *   // get parameters into the program
 *   std::cout << "\n\n" << "Getting parameters:" << std::endl;
 *   p.get_parameters (prm);
 *   // now run the program with these input parameters
 *   p.do_something ();
 * }
 * @endcode
 *
 *  这是输入文件（名为 "prmtest.prm"）。
 * @code
 * # first declare the types of equations
 * set Equation 1 = Poisson
 * set Equation 2 = Stokes
 *
 * subsection Equation 1 Settings
 *   set Matrix type = Sparse
 *   subsection Linear solver # parameters for linear solver 1
 *     set Solver                       = Gauss-Seidel
 *     set Maximum number of iterations = 40
 *   end
 * end
 *
 * subsection Equation 2 Settings
 *   set Matrix type = Full
 *   subsection Linear solver
 *     set Solver                       = CG
 *     set Maximum number of iterations = 100
 *   end
 * end
 * @endcode
 *
 * 这里是程序的输出。
 * @code
 * Line <8> of file <prmtest.prm>:
 *     The entry value
 *         Gauss-Seidel
 *     for the entry named
 *         Solver
 *     does not match the given pattern
 *         [Selection CG|BiCGStab|GMRES ]
 *
 *
 *
 *
 * # Listing of Parameters
 * #
 *
 * ---------------------
 * # String identifying the equation we want to solve
 * set Equation 1  = Poisson # default: Laplace
 * set Equation 2  = Stokes  # default: Elasticity
 *
 * # Name of the output file, either relative to the present path or absolute
 * set Output file = out
 *
 *
 *
 *
 * subsection Equation 1 Settings
 *   # Type of the matrix to be used, either full, sparse, or diagonal
 *   set Matrix type = Sparse
 *
 *
 *
 *
 *   subsection Linear solver
 *     set Maximum number of iterations = 40 # default: 20
 *     # Name of a linear solver for the inner iteration
 *     set Solver                       = CG
 *   end
 *
 * end
 *
 *
 *
 *
 * subsection Equation 2 Settings
 *   set Matrix type = Full # default: Sparse
 *
 *
 *
 *
 *   subsection Linear solver
 *     set Maximum number of iterations = 100 # default: 20
 *     # Name of a linear solver for the inner iteration
 *     set Solver                       = CG
 *   end
 *
 * end
 *
 *
 *
 *
 *
 *
 *
 *
 * Getting parameters:
 *   LinearEquation: method=CG, max_iterations=40
 *   LinearEquation: method=CG, max_iterations=100
 *   Problem: outfile=out
 *            eq1=Poisson, eq2=Stokes
 *            matrix1=Sparse, matrix2=Full
 * @endcode
 *
 *
 *   <h3>Representation of Parameters</h3>
 * 下面是一些关于参数表示法的更多内部信息。
 * 从逻辑上讲，参数和它们所排列的嵌套部分可以被认为是一个分层的目录结构，或者说是一棵树。以下面的代码为例，它声明了一组参数和它们所在的部分。
 * @code
 *   ParameterHandler prm;
 *
 *   prm.declare_entry ("Maximal number of iterations",
 *                      "10",
 *                      Patterns::Integer (1, 1000),
 *                      "A parameter that describes the maximal number of "
 *                      "iterations the CG method is to take before giving "
 *                      "up on a matrix.");
 *   prm.enter_subsection ("Preconditioner");
 *   {
 *     prm.declare_entry(
 *       "Kind",
 *       "SSOR",
 *       Patterns::Selection ("SSOR|Jacobi"),
 *       "A string that describes the kind of preconditioner to use.");
 *     prm.declare_entry(
 *       "Relaxation factor",
 *       "1.0",
 *       Patterns::Double (0, 1),
 *       "The numerical value (between zero and one) for the "
 *       "relaxation factor to use in the preconditioner.");
 *   }
 *   prm.leave_subsection ();
 * @endcode
 *
 * 我们可以把这样安排的参数看作一个文件系统，其中每个参数都是一个目录。这个目录的名称是参数的名称，在这个目录中存在着描述参数的文件。在编写本文档时，这些文件是（其他字段，如表示
 * "行动 "的字段也可能存在于每个目录中）。
 *
 *
 *
 * -  <code>value</code>  : 这个文件的内容是这个参数的当前值；最初，这个文件的内容等于参数的默认值。
 *
 *
 *
 * -  <code>default_value</code>  : 这个文件的内容是该参数的默认值。
 *
 *
 *
 * -  <code>pattern</code>  : 描述参数可能值的模式的文字表述。
 *
 *
 *
 * -  <code>pattern_index</code>  : 一个索引 Patterns::PatternBase 对象的数字，该对象用于描述该参数。
 *
 *
 *
 * -  <code>documentation</code>  : 这个文件的内容是为作为 ParameterHandler::declare_entry 调用的最后一个参数的参数提供的文件。除了 <code>value</code> 文件外，文件的内容在声明一个参数后永远不会改变。
 * 另外，这个文件系统中的一个目录可能没有一个叫做
 * <code>value</code>
 * 的文件在里面。在这种情况下，该目录代表上面声明的一个分节，目录的名称将与分节的名称相对应。然后，它里面根本没有文件，但它可能有更多的目录：其中一些目录将是参数（通过文件的存在表示）或进一步嵌套的子节。
 * 鉴于这种解释，上面的代码将导致数据的分层表示，看起来像这样（文件的内容在右边用不同的字体表示）。
*  @image html parameter_handler.png
 * 一旦参数被读入， <code>value</code> "文件
 * "的内容可能会有所不同，而其他文件则保持不动。 使用
 * ParameterHandler::print_parameters() 函数，将 ParameterHandler::XML
 * 作为第二个参数，我们可以得到这个数据结构在XML中的完整表示。它将看起来像这样。
 * @code
 * <?xml version="1.0" encoding="utf-8"?>
 * <ParameterHandler>
 *   <Maximal_20number_20of_20iterations>
 *     <value>10</value>
 *     <default_value>10</default_value>
 *     <documentation>
 *       A parameter that describes the maximal number of iterations the CG
 *       method is to take before giving up on a matrix.
 *     </documentation>
 *     <pattern>0</pattern>
 *     <pattern_description>
 *       [Integer range 1...1000 (inclusive)]
 *     </pattern_description>
 *   </Maximal_20number_20of_20iterations>
 *   <Preconditioner>
 *     <Kind><value>SSOR</value>
 *       <default_value>SSOR</default_value>
 *       <documentation>
 *         A string that describes the kind of preconditioner to use.
 *       </documentation>
 *       <pattern>1</pattern>
 *       <pattern_description>SSOR|Jacobi</pattern_description>
 *     </Kind>
 *     <Relaxation_20factor>
 *       <value>1.0</value>
 *       <default_value>1.0</default_value>
 *       <documentation>
 *         The numerical value (between zero and one) for the relaxation
 *         factor to use in the preconditioner.
 *       </documentation>
 *       <pattern>2</pattern>
 *       <pattern_description>
 *         [Floating point range 0...1 (inclusive)]
 *       </pattern_description>
 *     </Relaxation_20factor>
 *   </Preconditioner>
 * <ParameterHandler>
 * @endcode
 * 这种表示法与上面讨论的目录/文件结构非常相似。唯一不同的是，目录和文件名被篡改了：因为它们应该只包含字母和数字，所以它们名称中每一个不是字母或数字的字符都被下划线替换，后面是其两位数的十六进制表示。此外，特殊名称
 * "value
 * "在用作参数名称时也会被处理，因为这个名称也用于命名层次结构中的特殊文件。最后，整个树被包装成一个标签
 * <code>%ParameterHandler</code>
 * ，以满足XML的要求，即每个文件中只有一个顶级结构。
 * 树状结构（及其XML表示）是图形用户界面（见上文）用来表示目录/文件集合等参数的。
 *
 *
 *
 * @ingroup input
 *
 */
class ParameterHandler : public Subscriptor
{
public:
  /**
   * 用于 ParameterHandler::print_parameters().
   * 等函数的可能输出格式列表 选项可分为两组。
   *
   *
   *
   *
   *
   * - 格式选项。PRM, LaTeX, Description, XML, JSON
   *
   *
   *
   *
   * - 文体选项。Short, KeepDeclarationOrder 每次只能指定一个格式选项。任何接受OutputStyle作为选项的函数，如果你指定了一个以上的选项，就会抛出。    我们提供了一些常用选项组合的快捷方式。  例如，ShortPRM以PRM格式打印参数，同时跳过文档。
   *
   */
  enum OutputStyle
  {
    /**
     * 默认文体风格：打印文档，并按字母顺序排列所有参数。
     *
     */
    DefaultStyle = 0x0000,

    /**
     * 为ParameterHandler写输入，没有注释或改变默认值。
     *
     */
    Short = 0x0001,

    /**
     * 保持参数的顺序，因为它们已经被声明了。
     *
     */
    KeepDeclarationOrder = 0x0002,

    /**
     * 编写适合由 ParameterHandler::parse_input()
     * 再次读取的人类可读输出。
     *
     */
    PRM = 0x0010,

    /**
     * 写出适合由 ParameterHandler::parse_input()
     * 再次读取的人类可读输出。          @deprecated
     * 使用`PRM`而不是`Text`。
     *
     */
    Text = PRM,

    /**
     * 将参数写成LaTeX表格。
     *
     */
    LaTeX = 0x0020,

    /**
     * 写出声明的参数与描述和可能的值。
     * @note 这种格式不适合再读回来。
     *
     */
    Description = 0x0040,

    /**
     * 将所有内容写成一个<a
     * href="http://en.wikipedia.org/wiki/XML">XML</a>的文件，适合由
     * ParameterHandler::parse_input_from_xml() 再次读取。
     * 关于输出的例子，请看这个类的一般文档。
     *
     */
    XML = 0x0080,

    /**
     * 把所有的东西都写成一个<a
     * href="http://en.wikipedia.org/wiki/JSON">JSON</a>文件，适合由
     * ParameterHandler::parse_input_from_json() 再次读取。
     *
     */
    JSON = 0x0100,

    /**
     * 写出ParameterHandler的内容，没有注释或改变默认值。
     *
     */
    ShortPRM = PRM | Short,

    /**
     * 写下ParameterHandler的内容，不加评论或改变默认值。
     * @deprecated  使用`ShortPRM`而不是`ShortText`。
     *
     */
    ShortText = ShortPRM,

    /**
     * 把ParameterHandler的内容写成XML文件，不加注释或改变默认值。
     *
     */
    ShortXML = XML | Short,

    /**
     * 把ParameterHandler的内容写成JSON文件，不加注释或改变默认值。
     *
     */
    ShortJSON = JSON | Short,

    /**
     * 把ParameterHandler的内容写成LaTeX文件，不加注释或改变默认值。
     *
     */
    ShortLaTeX = LaTeX | Short,
  };



  /**
   * 构造函数。
   *
   */
  ParameterHandler();

  /**
   * 解构器。声明这个只是为了有一个虚拟的析构器，这更安全，因为我们有虚拟函数。
   * 它实际上没有什么了不起的作用。
   *
   */
  virtual ~ParameterHandler() override = default;

  /**
   * 抑制自动的CopyConstructor。
   *
   */
  ParameterHandler(const ParameterHandler &) = delete;

  /**
   * 抑制自动赋值运算符。
   *
   */
  ParameterHandler &
  operator=(const ParameterHandler &) = delete;

  /**
   * 解析流中的每一行，直到流返回<tt>eof</tt>条件或错误，以提供已知参数字段的值。第二个参数可以用来表示我们正在读取的文件的名称（如果那是输入流所代表的）；这只在为异常创建输出时使用。
   * 如果提供了非空的  @p last_line
   * ，ParameterHandler对象将在遇到  @p last_line
   * 后停止解析行。
   * 这在添加应被手动解析的额外数据时很方便。    如果
   * @p skip_undefined  是  <code>true</code>
   * ，参数处理程序将跳过未定义的部分和条目。这对于部分解析参数文件非常有用，例如，只获得问题的空间维度。默认情况下，所有的条目和子节都应该被声明。
   * 该函数将它在输入文件中遇到的所有参数的值设置为提供的值。在输入文件中没有明确列出的参数将保持其先前的值，这将是提供给
   * declare_entry()
   * 的默认值，除非先前读取过不同的输入文件。
   * 每个参数的值都会与提供给 declare_entry()
   * 的该参数的模式相匹配，并对每个参数执行之前可能由
   * add_action()
   * 设置的所有相关动作。如果一个参数不符合它的模式，或者一个相关的动作抛出了一个异常，那么为该参数提供的值就不会被设置，当前对象就会恢复到当前函数被调用之前的小节。对输入流不会有进一步的处理，也就是说，在参数的值不满足其模式之后的所有东西都被忽略了。
   *
   */
  virtual void
  parse_input(std::istream &     input,
              const std::string &filename       = "input file",
              const std::string &last_line      = "",
              const bool         skip_undefined = false);

  /**
   * 解析来自指定参数文件的输入 @p filename
   * ，与正在使用的输入文件类型（prm、xml、json）无关。这个函数选择的代码路径是从文件名的结尾处提取的，所以用户必须确保输入文件的内容与文件名一致。
   * 参数 @p last_line 将只用于.prm类型的参数文件。
   * 参见其他parse_input函数的文档。
   * 用户可以指定输入文件中未加入参数处理程序的参数是否会被
   * @p skip_undefined
   * 跳过（启用部分解析），以及代码是否会断言所有用标志`has_to_be_set=true`声明的参数处理程序的参数确实在输入文件中找到。
   * 如果函数是以`skip_undefined=true`调用的，建议同时设置`assert_mandatory_entries_are_found=true`。例如，这可以确保在输入文件中有错别字的参数不会被跳过，而这种错误在其他情况下仍然无法被识别。
   *
   */
  virtual void
  parse_input(const std::string &filename,
              const std::string &last_line                          = "",
              const bool         skip_undefined                     = false,
              const bool         assert_mandatory_entries_are_found = false);

  /**
   * 从一个字符串中解析输入，以填充已知的参数字段。字符串中的各行必须用<tt>
   * @\n</tt> 字符分开。
   * 该函数实质上是将整个文件读入一个流，然后用该流调用其他的parse_input()函数。更多信息见那里。
   *
   */
  virtual void
  parse_input_from_string(const std::string &s,
                          const std::string &last_line      = "",
                          const bool         skip_undefined = false);

  /**
   * 从一个XML流中解析输入，以填充已知的参数字段。这可能来自一个最初由print_parameters()函数使用XML输出风格编写的文件，然后根据需要手工修改，或者来自一个使用该方法编写的文件，然后由图形化参数GUI修改（见该类的一般文档）。
   *
   */
  virtual void
  parse_input_from_xml(std::istream &input, const bool skip_undefined = false);

  /**
   * 解析来自JSON流的输入，以填充已知的参数字段。这可能来自于最初由print_parameters()函数使用JSON输出风格编写的文件，然后根据需要手工修改，或者来自一个知道如何为ParameterHandler输入编写JSON格式的单独程序。
   *
   */
  virtual void
  parse_input_from_json(std::istream &input, const bool skip_undefined = false);

  /**
   * 清除所有内容。
   *
   */
  void
  clear();


  /**
   * 声明一个新的条目，名称为<tt>entry</tt>，默认为任何输入都必须与<tt>pattern</tt>相匹配（默认：任何模式）。
   * 如果默认值与给定的模式不匹配，该函数会产生一个ExcValueDoesNotMatchPattern类型的异常，使用C++的抛出机制。然而，这个异常只在条目被创建时产生<i>after</i>；如果你的代码中不可能有合理的参数默认值，那么你可以捕捉并忽略这个异常。
   * 参数 @p documentation
   * 默认为一个空字符串，用于为每个条目添加记录文本，当这个类被要求使用print_parameters()函数将所有声明写到一个流中时，该文本将被打印出来作为注释。
   * 可以使用参数 @p has_to_be_set
   * 来声明这个参数，这个参数的默认值必须被这个类所提供的一个方法所覆盖。一个参数是否被成功设置，可以通过函数get_entries_wrongly_not_set()和assert_that_entries_have_been_set()进行查询。
   * @note
   * 一个条目可以被声明多次而不产生错误，例如，为了覆盖一个早期的默认值。
   *
   */
  void
  declare_entry(const std::string &          entry,
                const std::string &          default_value,
                const Patterns::PatternBase &pattern = Patterns::Anything(),
                const std::string &          documentation = "",
                const bool                   has_to_be_set = false);

  /**
   * 在当前章节中为名称为 @p entry
   * 的参数附加一个动作。这个动作需要是一个类似于函数的对象，把参数的值作为一个（字符串）参数。关于动作的更多描述，请参见这个类的一般文档，以及例子。
   * 该动作在三种不同的情况下被执行。
   *
   *
   *
   *
   *
   * - 在当前函数的末尾加上名称为 @p name, 的参数的默认值。这很有用，因为它允许动作对每个参数至少执行一次它需要做的事情，即使是那些实际上没有在输入文件中指定的参数（因此保持其默认值）。
   *
   *
   *
   *
   *
   *
   * - 在 ParameterHandler::set() 函数中，明确地为一个参数设置一个值。
   *
   *
   *
   *
   *
   *
   * - 在parse_input()函数和类似的函数中，例如parse_input_from_string()。在这里，只要从输入中读取与之相关的参数，在确定所读取的值与该参数所对应的模式相匹配后，在实际保存该值之前，就会执行该动作。    为同一个参数添加多个动作是有效的。在这种情况下，它们将按照添加的顺序被执行。
   * @note
   * 动作可以修改其范围内的所有种类的变量。一个动作唯一不应该修改的是它所连接的ParameterHandler对象。换句话说，它不允许进入或离开当前ParameterHandler对象的部分。原则上，在当前部分的其他参数上调用
   * ParameterHandler::get()
   * 和相关函数是可以接受的，但由于不能保证它们从输入文件中读取的顺序，你将不希望依赖这些函数将返回的值。
   * @note
   * 在一个动作中抛出一个异常通常不是一个好主意，但产生的结果与试图从一个文件中读取一个参数的结果基本相同，因为该参数的值不符合与该参数相关的模式。换句话说，刚刚读取的值被丢弃，
   * ParameterHandler::parse_input()
   * 停止从文件中读取任何进一步的内容。更多信息见
   * ParameterHandler::parse_input() 。
   *
   */
  void
  add_action(const std::string &                                  entry,
             const std::function<void(const std::string &value)> &action);

  /**
   * 声明一个新的条目名称 @p entry, 将其默认值设置为变量
   * @p parameter,
   * 的内容，并创建一个动作，当文件被解析，或条目被设置为新值时，将用更新的值填充
   * @p 参数。    默认情况下，要使用的模式是通过调用函数
   * Patterns::Tools::Convert<T>::to_pattern(),
   * 获得的，但也可以使用一个自定义的模式。
   * 可以使用参数 @p has_to_be_set
   * 来声明这个参数，其默认值必须被这个类提供的一个方法所覆盖。一个参数是否被成功设置，可以通过函数get_entries_wrongly_not_set()和assert_that_entries_have_been_set()进行查询。
   *
   */
  template <class ParameterType>
  void
  add_parameter(const std::string &          entry,
                ParameterType &              parameter,
                const std::string &          documentation = "",
                const Patterns::PatternBase &pattern =
                  *Patterns::Tools::Convert<ParameterType>::to_pattern(),
                const bool has_to_be_set = false);

  /**
   * 为一个现有条目创建一个别名。这提供了一种方法，可以用另一个名字来引用输入文件中的一个参数。该别名将在当前部分，而被引用的条目需要是当前部分的一个现有条目。
   * 这个函数的主要目的是允许用一种向后兼容的方式来改变应用程序输入文件中的名称，因为向后兼容很重要。这可以通过在调用
   * declare_entry()
   * 时改变参数的名称，然后创建一个别名将旧名称映射到新名称来实现。这样，旧的输入文件可以继续引用旧名称下的参数，它们将自动被映射到新的参数名称。
   * 在一个输入文件中多次设置同一个参数是有效的。
   * 在这种情况下，最终选择的值仅仅是最后设置的值。这个规则也适用于别名，参数的最终值是通过参数的当前名称或通过其可能的多个别名设置的最后一个值。例如，如果你有一个输入文件，看起来像
   * @code
   * set parm1       = 1
   * set parm1_alias = 2
   * @endcode
   * 其中 <code>parm1_alias</code>
   * 是一个通过以下方式声明的别名
   * @code
   * prm.declare_alias ("parm1", "parm1_alias");
   * @endcode
   * 那么名为 <code>parm1</code>
   * 的参数的最终值将是2，而不是1。      @param
   * existing_entry_name
   * 在当前章节中的一个现有参数的名称，该别名应该指的是这个。
   * @param  alias_name 第一个参数所引用的参数的另一个名称。
   * @param  alias_is_deprecated
   * 如果为真，将该别名标记为已废弃。如果你调用print_parameters()，这将被列在别名的描述中，并且在读取包含这个已废弃的别名的输入文件时，你将在屏幕上得到一个警告。这个参数的目的是能够允许使用一个旧的参数名称（见上文），但是要明确这个旧名称最终会被删除。
   *
   */
  void
  declare_alias(const std::string &existing_entry_name,
                const std::string &alias_name,
                const bool         alias_is_deprecated = false);

  /**
   * 输入一个分节。如果它还不存在，就创建它。
   *
   */
  void
  enter_subsection(const std::string &subsection);

  /**
   * 离开目前的分节。
   *
   */
  void
  leave_subsection();

  /**
   * 检查当前树中是否存在一个分节或分节路径。
   * 输入参数 @p sub_path 被认为是相对于当前选择的路径。
   *
   */
  bool
  subsection_path_exists(const std::vector<std::string> &sub_path) const;

  /**
   * 入口的返回值  @p entry_string.
   * 如果入口被改变，则返回改变后的值，否则返回默认值。如果需要一个未声明的条目的值，
   * @p Assert 将失败。
   *
   */
  std::string
  get(const std::string &entry_string) const;

  /**
   * 条目 @p entry_string. 的返回值
   * 如果条目被改变，那么返回改变后的值，否则返回默认值。如果需要一个未声明的条目的值，
   * @p Assert 将失败。  如果 @p entry_subsection_path
   * 是非空的，那么将从该路径所代表的分节中获取数值，而不是当前分节。
   * @p entry_subsection_path
   * 中的第一个字符串必须是当前分节的名称，接下来的每个字符串必须是前面一个分节的名称。
   *
   */
  std::string
  get(const std::vector<std::string> &entry_subsection_path,
      const std::string &             entry_string) const;

  /**
   * 返回条目 @p entry_string 的值为 <code>long int</code>
   * 。（选择一个长的int，这样即使是非常大的无符号值也可以由这个函数返回）。
   *
   */
  long int
  get_integer(const std::string &entry_string) const;

  /**
   * 条目 @p entry_string 的返回值为 <code>long int</code>
   * 。（选择一个长的int，这样即使是非常大的无符号值也能被这个函数返回）。
   * 如果 @p entry_subsection_path
   * 是非空的，那么该值将从该路径所代表的分节而不是当前分节中得到。
   *
   */
  long int
  get_integer(const std::vector<std::string> &entry_subsection_path,
              const std::string &             entry_string) const;

  /**
   * 返回条目 @p entry_name 的值为 @p double. 。
   *
   */
  double
  get_double(const std::string &entry_name) const;

  /**
   * 如果 @p entry_subsection_path
   * 不是空的，将从该路径所代表的分节而不是当前分节中获取数值。
   *
   */
  double
  get_double(const std::vector<std::string> &entry_subsection_path,
             const std::string &             entry_string) const;
  /**
   * 条目 @p entry_name 的返回值为 @p bool. ，该条目对 @p true,
   * 可以是 "真 "或 "是"，对 @p false 可以是 "假 "或 "否"。
   *
   */
  bool
  get_bool(const std::string &entry_name) const;

  /**
   * 条目 @p entry_name 的返回值为 @p bool.  该条目对于 @p true,
   * 可以分别为 "真 "或 "是"，对于 @p false 可以为 "假 "或
   * "否"。  如果 @p entry_subsection_path
   * 不为空，则将从该路径所代表的分节中获取数值，而不是当前分节。
   *
   */
  bool
  get_bool(const std::vector<std::string> &entry_subsection_path,
           const std::string &             entry_string) const;

  /**
   * 将当前存储的<tt>entry_name</tt>的值改为第二个参数中给出的值。
   * 该参数必须已经存在于当前分节中。
   * 如果新值不符合该条目的模式，该函数会抛出一个ExcValueDoesNotMatchPattern类型的异常。
   *
   */
  void
  set(const std::string &entry_name, const std::string &new_value);

  /**
   * 和上面一样，但有一个重载，第二个参数是一个字符指针。这是必要的，因为否则对
   * <code>set("abc","def")</code>
   * 的调用将被映射为该函数以一个字符串和一个bool作为参数，这当然不是最常见的目的。
   * 如果新值不符合这个条目的模式，该函数会抛出一个ExcValueDoesNotMatchPattern类型的异常。
   *
   */
  void
  set(const std::string &entry_name, const char *new_value);

  /**
   * 将目前为<tt>entry_name</tt>存储的值改为第二个参数中给出的值。
   * 该参数必须已经存在于本小节中。
   * 如果新值不符合该条目的模式，该函数会抛出一个ExcValueDoesNotMatchPattern类型的异常。
   *
   */
  void
  set(const std::string &entry_name, const long int new_value);

  /**
   * 将目前为<tt>entry_name</tt>存储的值改为第二个参数中给出的值。
   * 该参数必须已经存在于本分节中。
   * 为了内部的目的，新的值需要转换为一个字符串。
   * 这是用16位数的精度完成的，所以设定的值和你用get_double()得到的值可能在第16位数上有所不同。
   * 如果新值不符合这个条目的模式，该函数会抛出一个ExcValueDoesNotMatchPattern类型的异常。
   *
   */
  void
  set(const std::string &entry_name, const double new_value);

  /**
   * 将目前为<tt>entry_name</tt>存储的值改为第二个参数中给出的值。
   * 该参数必须已经存在于本小节中。
   * 如果新值不符合该条目的模式，该函数会抛出一个ExcValueDoesNotMatchPattern类型的异常。
   *
   */
  void
  set(const std::string &entry_name, const bool new_value);

  /**
   * 打印所有具有给定 @p style 至 @p out. 的参数
   * 在打印之前，所有当前参数和子部分默认按字母顺序排序。
   * 这种行为可以通过设置可选的参数 @p style
   * 来禁用<tt>KeepDeclarationOrder</tt>：在这种情况下，条目的打印顺序与它们被声明的顺序相同。
   * 在<tt>PRM</tt>、<tt>XML</tt>和<tt>JSON</tt>格式中，输出的格式是可以用于以后再次输入。这对于记录特定运行的参数是最有用的，因为如果你用这个函数把参数输出到一个日志文件中，你总是可以通过简单地把输出复制到输入文件来恢复结果。
   * 除了每个条目的名称和值之外，输出还包括条目的默认值（如果它与实际值不同），以及给
   * declare_entry() 函数的记录字符串（如果有）。
   * 通过使用标志<tt>Short</tt>与<tt>PRM</tt>、<tt>XML</tt>、<tt>JSON</tt>或<tt>LaTeX</tt>相结合（或通过使用快捷键<tt>ShortPRM</tt>。<tt>ShortXML</tt>,
   * <tt>ShortJSON</tt>, 或
   * <tt>ShortLaTeX</tt>），可以生成一个缩小的输出，只包含数值，跳过文档。
   * 在<tt>XML</tt>格式中，输出从一个根元素<tt>ParameterHandler</tt>开始，以便得到一个有效的XML文档和它下面的所有子节。
   * 在<tt>LaTeX</tt>格式中，输出包含同样的信息，但其格式是这样的：所产生的文件可以输入到一个latex文档中，例如这个对象处理运行时参数的代码的手册。然后，参数的各个部分由latex章节和分节命令以及嵌套的枚举来表示。
   * 你可以通过为每个条目自动生成的标签来引用特定的参数部分和单个参数。标签的格式为
   * <code>parameters:section1/subsection1</code>  和
   * <code>parameters:section1/subsection1/someentry</code>
   * 。由于特殊字符可能出现在章节和条目名称中，这些字符将被
   * "混杂"。在这里，除了 <code>[a-zA-Z0-9]</code>
   * 以外的所有字符都被 <code>_XX</code>, where <code>XX</code>
   * 所取代，该字符在十六进制编码中的两位数ascii代码（因此，例如，一个空格变成
   * <code>_20</code> ）。
   * 虽然这个函数在大多数输出中（名称、默认值等）都转义了LaTeX特有的字符（反斜杠、下划线等），但文档字符串是按原样传递的。这意味着你可以在描述中使用数学环境和其他格式，但你需要自己转义引号、反斜线、下划线等。
   * 此外，所有的参数名称都用 <code>@\index</code>
   * 语句列在两个索引中，称为 <code>prmindex</code>
   * （每个参数的名称都列在索引中）和
   * <code>prmindexfull</code>
   * ，参数名称按其存在的章节排序。默认情况下，LaTeX程序会忽略这些
   * <code>@\index</code>
   * 命令，但它们可以通过在latex文件的序言中使用以下命令来生成索引。
   * @code
   * \usepackage{imakeidx}
   * \makeindex[name=prmindex, title=Index of run-time parameter entries]
   * \makeindex[name=prmindexfull,
   *          title=Index of run-time parameters with section names]
   * @endcode
   * 并在文件的结尾这样做。
   * @code
   * \printindex[prmindex]
   * \printindex[prmindexfull]
   * @endcode
   *
   */
  std::ostream &
  print_parameters(std::ostream &out, const OutputStyle style) const;



  /**
   * 打印所有参数到 @p filename
   * 给定的文件中，并使用给定的输出样式 @p style.
   * 这个函数从指定文件名的扩展名推断出输出格式。支持的扩展名是`prm',
   * `xml', `tex', 和
   * `json'。因此，只要在文件名中加入这些扩展名之一，就没有必要通过
   * @p style 参数指定输出格式。不过，如果在 @p style
   * 参数中指定了输出格式，输出格式必须与文件名的扩展名一致。
   * 如果没有指定扩展名或不支持扩展名，输出格式将从 @p
   * style 参数中推断出来。
   * 如果既不支持扩展名，也不支持 @p style
   * 参数中的格式规范，则会抛出一个断言。      @param
   * filename 输出文件名。    @param  style 产生输出的样式。
   *
   */
  void
  print_parameters(const std::string &filename,
                   const OutputStyle  style = DefaultStyle) const;

  /**
   * 打印参数到一个日志流。这个函数允许将所有参数打印到一个日志文件中。各部分将以通常的日志文件风格缩进。
   * 默认情况下，所有当前的参数和子部分都是按字母顺序排列的。
   * 这种行为可以通过设置可选的参数 @p style
   * 来禁用<tt>KeepDeclarationOrder</tt>：在这种情况下，条目的打印顺序与它们被声明的顺序相同。
   * @note  @p style 中所有与排序无关的样式设置都被忽略了。
   *
   */
  void
  log_parameters(LogStream &out, const OutputStyle style = DefaultStyle);

  /**
   * 本分节中的日志参数。该分节由<tt>subsection_path</tt>成员变量决定。这个变量通过enter_subsection()和leave_subsection()函数进入和离开子段来控制。
   * 所有当前的参数和子段默认是按字母顺序排序的。
   * 这种行为可以通过设置可选的参数 @p style
   * 来禁用<tt>KeepDeclarationOrder</tt>：在这种情况下，条目的打印顺序与它们被声明的顺序相同。
   * @note   @p style
   * 中所有与排序无关的样式设置都被忽略了。
   * 在大多数情况下，你不会想直接使用这个函数，而是让它被前面的函数递归调用。
   *
   */
  void
  log_parameters_section(LogStream &       out,
                         const OutputStyle style = DefaultStyle);

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写入一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  save(Archive &ar, const unsigned int version) const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从一个流中读取此对象的数据，以达到序列化的目的。
   *
   */
  template <class Archive>
  void
  load(Archive &ar, const unsigned int version);

#ifdef DOXYGEN
  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从流中写入和读取此对象的数据，以达到序列化的目的。
   *
   */
  template <class Archive>
  void
  serialize(Archive &archive, const unsigned int version);
#else
  // This macro defines the serialize() method that is compatible with
  // the templated save() and load() method that have been implemented.
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

  /**
   * 检验是否相等。
   *
   */
  bool
  operator==(const ParameterHandler &prm2) const;

  /**
   * 返回一组参数名称（包括分段名称），对应于参数处理程序的那些条目，这些条目没有被一个从输入文件解析参数的函数或明确调用set()函数之一来设置，但已被声明为必须设置的强制性参数（通过
   * declare_entry() 函数或 add_parameter() 函数的最后参数）。
   *
   */
  std::set<std::string>
  get_entries_wrongly_not_set() const;

  /**
   * 断言参数处理程序中那些标记为`has_to_be_set =
   * true`的条目已经被设置。如果这些参数中至少有一个没有被设置，就会出现异常。
   *
   */
  void
  assert_that_entries_have_been_set() const;

  /**
   * @addtogroup  Exceptions  
     * @{ 
   *
   */

  /**
   * 异常情况
   *
   */
  DeclException1(ExcEntryAlreadyExists,
                 std::string,
                 << "The following entry already exists: " << arg1 << ".");
  /**
   * 异常情况
   *
   */
  DeclException2(ExcValueDoesNotMatchPattern,
                 std::string,
                 std::string,
                 << "The string <" << arg1
                 << "> does not match the given pattern <" << arg2 << ">.");
  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(
    ExcAlreadyAtTopLevel,
    "You can't leave a subsection if you are already at the top level "
    "of the subsection hierarchy.");
  /**
   * 异常情况
   *
   */
  DeclException1(ExcEntryUndeclared,
                 std::string,
                 << "You can't ask for entry <" << arg1
                 << "> you have not yet declared.");

  /**
   * 当 "分节 "和 "结束
   * "语句的数量不相等时的例外。第一个参数是文件名，第二个参数是进入解析器之前和之后的分节路径的格式化列表。
   *
   */
  DeclException2(ExcUnbalancedSubsections,
                 std::string,
                 std::string,
                 << "There are unequal numbers of 'subsection' and 'end' "
                    "statements in the parameter file <"
                 << arg1 << ">." << (arg2.size() > 0 ? "\n" + arg2 : ""));

  /**
   * 当在解析参数文件时，解析器遇到文件中的分节没有事先声明时的异常情况。
   *
   */
  DeclException3(ExcNoSubsection,
                 int,
                 std::string,
                 std::string,
                 << "Line <" << arg1 << "> of file <" << arg2
                 << ": There is "
                    "no such subsection to be entered: "
                 << arg3);

  /**
   * 无法解析的行的一般异常，以行号、文件名和对该行无法解析的原因的简要描述作为参数。
   *
   */
  DeclException3(ExcCannotParseLine,
                 int,
                 std::string,
                 std::string,
                 << "Line <" << arg1 << "> of file <" << arg2 << ">: " << arg3);

  /**
   * 参数文件中不符合所提供模式的条目的异常。参数依次是：行号、文件名、条目值、条目名和模式描述。
   *
   */
  DeclException5(ExcInvalidEntryForPattern,
                 int,
                 std::string,
                 std::string,
                 std::string,
                 std::string,
                 << "Line <" << arg1 << "> of file <" << arg2
                 << ">:\n"
                    "    The entry value \n"
                 << "        " << arg3 << '\n'
                 << "    for the entry named\n"
                 << "        " << arg4 << '\n'
                 << "    does not match the given pattern:\n"
                 << "        " << arg5);

  /**
   * 当一个XML文件完全不能被读取时的异常。当没有名为
   * "ParameterHandler
   * "的顶层XML元素或有多个顶层元素时就会发生这种情况。
   *
   */
  DeclExceptionMsg(ExcInvalidXMLParameterFile,
                   "The provided file could not be parsed as a "
                   "ParameterHandler description.");

  /**
   * 当一个XML参数文件中的条目与所提供的模式不匹配时的异常。参数依次是条目值、条目名称和模式的描述。
   * @deprecated
   * 使用ExcValueDoesNotMatchPattern代替ExcInvalidEntryForPatternXML。
   *
   */
  using ExcInvalidEntryForPatternXML DEAL_II_DEPRECATED =
    ExcValueDoesNotMatchPattern;

  /**
   * 当include语句中给出的文件不能被打开时的异常。参数依次是include语句的行号，当前参数文件名，以及打算列入的文件名。
   *
   */
  DeclException3(
    ExcCannotOpenIncludeStatementFile,
    int,
    std::string,
    std::string,
    << "Line <" << arg1 << "> of file <" << arg2
    << ">: This line "
       "contains an 'include' or 'INCLUDE' statement, but the given "
       "file to include <"
    << arg3 << "> cannot be opened.");

  //@}

private:
  /**
   * 访问路径中的元素进入参数树时使用的分隔符。
   *
   */
  static const char path_separator = '.';

  /**
   * 目前选择的子部分的路径；空列表表示顶层
   *
   */
  std::vector<std::string> subsection_path;

  /**
   * 完整的章节和条目树。关于数据如何存储在这个变量中的描述，请参见这个类的一般文档。
   * 这个变量是一个指针，这样我们就可以使用一个不完整的类型，而不需要包括boost中的所有property_tree的东西。这可以解决gcc
   * 4.5的一个问题。
   *
   */
  std::unique_ptr<boost::property_tree::ptree> entries;

  /**
   * 一个地图，为每个添加到参数处理程序的条目存储一对布尔变量。第一个布尔变量描述了参数是否必须根据函数
   * declare_entry() 或 add_parameter()
   * 的最后一个参数来设置，第二个布尔变量包含了参数是否被任何解析输入参数的函数或本类的设置函数设置过的信息。
   *
   */
  std::map<std::string, std::pair<bool, bool>> entries_set_status;

  /**
   * 一个用于描述此对象的参数的模式列表。与参数相对应的属性树中的每一个节点都向这个数组存储一个索引。
   *
   */
  std::vector<std::unique_ptr<const Patterns::PatternBase>> patterns;

  /**
   * 一个与参数相关的动作列表。这些是由add_action()函数添加的。与单个参数相对应的属性树中的节点在这个数组中存储索引，以便引用特定的动作。
   *
   */
  std::vector<std::function<void(const std::string &)>> actions;

  /**
   * 返回标识当前进入属性树的路径的字符串。这只是一个路径，也就是说，它没有被path_separator字符终止。
   * 这个函数只是调用collate_path_string()，参数为 @p
   * subsection_path 。
   *
   */
  std::string
  get_current_path() const;

  /**
   * 给出条目名称作为参数，该函数使用当前分节计算出进入参数树的完整路径。
   *
   */
  std::string
  get_current_full_path(const std::string &name) const;

  /**
   * 这个函数给定一个来自当前分节的路径和一个条目的名称，计算一个进入参数树的全路径。
   *
   */
  std::string
  get_current_full_path(const std::vector<std::string> &sub_path,
                        const std::string &             name) const;

  /**
   * 扫描一行的输入。<tt>input_filename</tt>和<tt>current_line_n</tt>是输入文件的名称和当前扫描的行数（这些在异常信息中用来显示哪里发生了解析错误）。如果该行包含一个未声明的分段或条目，如果该行的条目不符合其给定的模式，或者如果该行不能被理解为一个有效的参数文件表达式，这个函数将引发一个异常。
   * 该函数修改了它的参数，但同时也以值的形式接受了它，所以调用者的变量没有被改变。
   * 如果  @p skip_undefined  是  <code>true</code>
   * ，解析器将跳过未定义的部分和条目。这对于部分解析参数文件是很有用的，例如，只获得问题的空间维度。默认情况下，所有的条目和子节都应该被声明。
   *
   */
  void
  scan_line(std::string        line,
            const std::string &input_filename,
            const unsigned int current_line_n,
            const bool         skip_undefined);

  /**
   * 打印出由 @p target_subsection_path
   * 参数给出的分节的参数，以及其中的所有分节的递归。这个函数是由print_parameters()函数调用的，并对所有
   * @p style
   * 参数实现，除了XML和JSON（在那里我们可以通过BOOST函数输出整个参数集）。
   * @p indent_level
   * 参数表示输出应该缩进多少个空格，以便子部分正确嵌套在更高部分的输出中。
   *
   */
  void
  recursively_print_parameters(
    const boost::property_tree::ptree & tree,
    const std::vector<std::string> &    target_subsection_path,
    const ParameterHandler::OutputStyle style,
    const unsigned int                  indent_level,
    std::ostream &                      out) const;

  friend class MultipleParameterLoop;
};

/**
 * 全局操作符，返回一个对象，其中所有的位都被设置，这些位在第一个或第二个参数中都被设置。这个操作符的存在是因为如果它不存在，那么位-或<tt>操作符|</tt>的结果将是一个整数，当我们试图将它赋值给一个类型为
 * ParameterHandler::OutputStyle.
 * 的对象时，会引发一个编译器警告。
 *
 *
 */
inline ParameterHandler::OutputStyle
operator|(const ParameterHandler::OutputStyle f1,
          const ParameterHandler::OutputStyle f2)
{
  return static_cast<ParameterHandler::OutputStyle>(
    static_cast<unsigned int>(f1) | static_cast<unsigned int>(f2));
}

/**
 * 多参数循环（MultipleParameterLoop）类提供了一种简单的可能性，可以在程序的一次运行中测试多个参数集。为此，它使用ParameterHandler类以标准化的形式读入数据，搜索变量入口值并对所有的参数组合进行循环。
 * 变体条目值是这样给出的。
 * @verbatim
 *   set Time step size = { 0.1 | 0.2 | 0.3 }
 * @endverbatim
 * 循环将进行三次程序运行，对<tt>时间步长</tt>的每个值进行一次运行，而所有其他参数都按指定值或默认值运行。如果在输入中存在几个变体条目值，则对每个变体值的组合进行一次循环。
 * @verbatim
 *   set Time step size = { 0.1 | 0.2 }
 *   set Solver         = { CG  | GMRES }
 * @endverbatim
 * 将导致程序的四次运行，两个求解器的时间步长分别为0.1和0.2。
 * 除了变体条目，这个类还支持<i>array
 * entries</i>，看起来像这样。
 * @verbatim
 *   set Output file = ofile.{{ 1 | 2 | 3 | 4 }}
 * @endverbatim
 * 这表明如果有变体条目产生了总共四个不同的运行，那么我们将把它们的结果分别写入文件<tt>ofile.1</tt>、<tt>ofile.2</tt>、<tt>ofile.3</tt>和<tt>ofile.4</tt>。数组条目本身不产生主循环的多次运行，但是如果有变体条目，那么在主循环的<i>n</i>次运行中，也会返回数组的<i>n</i>次值。
 * 由于不同的变体是按照声明的顺序构建的，而不是按照变体条目在输入文件中出现的顺序，所以可能很难猜测不同的变体和数组中适当的条目之间的映射。你将不得不检查声明的顺序，或者只使用一个变体条目。
 * 保证只有在声明一个条目时给出的正则表达式（模式）相匹配的选择才会被反馈给程序。如果一个变量的值与正则表达式不匹配，就会存储默认值并发出错误。在循环的第一次运行之前，所有可能的值都会被检查是否符合，因此错误会在程序的一开始就发出。
 *
 *  <h3>Usage</h3>
 * 这个类的用法与ParameterHandler类相似。首先必须声明条目和分节，然后进行循环，在循环中设置不同的参数集，创建一个新的用户类的实例，然后调用。以ParameterHandler类的例子为例，扩展后的程序会是这样的。
 * @code
 *   class HelperClass : public MultipleParameterLoop::UserClass
 *   {
 *   public:
 *     HelperClass ();
 *
 *     virtual void create_new (const unsigned int run_no);
 *     virtual void run (ParameterHandler &prm);
 *
 *     static void declare_parameters (ParameterHandler &prm);
 *   private:
 *     std::unique_ptr<Problem> p;
 *   };
 *
 *
 *
 *
 *   HelperClass::HelperClass () : p(0) {}
 *
 *
 *
 *
 *   void HelperClass::create_new (const unsigned int run_no)
 *   {
 *     p = std::make_unique<Problem>());
 *   }
 *
 *
 *
 *
 *   void HelperClass::declare_parameters (ParameterHandler &prm)
 *   {
 *     Problem::declare_parameters (prm);
 *   }
 *
 *
 *
 *
 *   void HelperClass::run (ParameterHandler &prm)
 *   {
 *     p->get_parameters (prm);
 *     p->do_useful_work ();
 *   }
 *
 *
 *
 *
 *
 *   int main ()
 *   {
 *     class MultipleParameterLoop prm;
 *     HelperClass h;
 *     HelperClass::declare_parameters (prm);
 *     prm.parse_input ("prmtest.prm");
 *     prm.loop (h);
 *     return 0;
 *   }
 * @endcode
 *
 * 可以看出，首先必须建立一个新的辅助类。这必须包含一个问题类的虚拟构造函数。你也可以从
 * MultipleParameterLoop::UserClass
 * 中派生出你的问题类，并让<tt>create_new</tt>清除所有成员变量。如果你能以某种方式访问所有继承的成员变量，这就是推荐的程序。第三种可能性是使用多重继承，从
 * MultipleParameterLoop::UserClass
 * 和问题类中派生出一个辅助类。在任何情况下，<tt>create_new</tt>都必须提供一个干净的问题对象，这是第二和第三种可能性的问题。
 * 派生类还必须提供声明条目和运行程序的成员函数。运行程序包括从ParameterHandler对象中获取参数。
 * 在定义了这个辅助类的对象和多参数循环类的对象后，必须以与ParameterHandler类相同的方式声明条目。然后要读取输入。最后，循环被调用。这将执行以下步骤。
 * @code
 *   for (each combination)
 *     {
 *       UserObject.create_new (run_no);
 *
 *       // set parameters for this run
 *
 *       UserObject.run (*this);
 *     }
 * @endcode
 * <tt>UserObject</tt>是<tt>loop</tt>函数的参数。<tt>create_new</tt>被赋予运行的编号（从1开始），以使每次运行的输出文件的命名不同。
 *
 *  <h3>Syntax for variant and array entry values</h3>
 * 变量值的指定与<tt>prefix{ v1 | v2 | v3 | ...
 * }postfix</tt>。开头括号<tt>{</tt>右边的空白会被忽略，而结尾括号<tt>}</tt>左边的空白也不会被忽略。中间符号<tt>|</tt>周围的空白也会被忽略。空的选择<tt>prefix{
 * v1 |
 * }postfix</tt>也是允许的，并产生字符串<tt>prefixv1postfix</tt>和<tt>prefixpostfix</tt>。
 * 除了双括号之外，数组值的语法是相同的。<tt>prefix{{ v1 |
 * v2 | v3 }}postfix</tt>。
 *
 *  <h3>Worked example</h3>
 * 给出上述对ParameterHandler示例程序的扩展和以下输入文件
 * @verbatim
 *   set Equation 1 = Poisson
 *   set Equation 2 = Navier-Stokes
 *   set Output file= results.{{ 1 | 2 | 3 | 4 | 5 | 6 }}
 *
 *   subsection Equation 1
 *     set Matrix type = Sparse
 *     subsection Linear solver
 *       set Solver                       = CG
 *       set Maximum number of iterations = { 10 | 20 | 30 }
 *     end
 *   end
 *
 *   subsection Equation 2
 *     set Matrix type = Full
 *     subsection Linear solver
 *       set Solver                       = { BiCGStab | GMRES }
 *       set Maximum number of iterations = 100
 *     end
 *   end
 * @endverbatim
 * 这就是输出结果。
 * @verbatim
 *   LinEq: Method=CG, MaxIterations=10
 *   LinEq: Method=BiCGStab, MaxIterations=100
 *   Problem: outfile=results.1
 *            eq1=Poisson, eq2=Navier-Stokes
 *            Matrix1=Sparse, Matrix2=Full
 *   LinEq: Method=CG, MaxIterations=20
 *   LinEq: Method=BiCGStab, MaxIterations=100
 *   Problem: outfile=results.2
 *            eq1=Poisson, eq2=Navier-Stokes
 *            Matrix1=Sparse, Matrix2=Full
 *   LinEq: Method=CG, MaxIterations=30
 *   LinEq: Method=BiCGStab, MaxIterations=100
 *   Problem: outfile=results.3
 *            eq1=Poisson, eq2=Navier-Stokes
 *            Matrix1=Sparse, Matrix2=Full
 *   LinEq: Method=CG, MaxIterations=10
 *   LinEq: Method=GMRES, MaxIterations=100
 *   Problem: outfile=results.4
 *            eq1=Poisson, eq2=Navier-Stokes
 *            Matrix1=Sparse, Matrix2=Full
 *   LinEq: Method=CG, MaxIterations=20
 *   LinEq: Method=GMRES, MaxIterations=100
 *   Problem: outfile=results.5
 *            eq1=Poisson, eq2=Navier-Stokes
 *            Matrix1=Sparse, Matrix2=Full
 *   LinEq: Method=CG, MaxIterations=30
 *   LinEq: Method=GMRES, MaxIterations=100
 *   Problem: outfile=results.6
 *            eq1=Poisson, eq2=Navier-Stokes
 *            Matrix1=Sparse, Matrix2=Full
 * @endverbatim
 * 由于<tt>create_new</tt>得到了运行的编号，所以也可以输出运行的编号。
 *
 *
 *
 * @ingroup input
 *
 *
 */
class MultipleParameterLoop : public ParameterHandler
{
public:
  /**
   * 这是辅助类或问题类必须派生的类。
   *
   */
  class UserClass
  {
  public:
    /**
     * 解构器。它实际上不做任何事情，但声明它是为了强迫派生类有一个虚拟的析构器。
     *
     */
    virtual ~UserClass() = default;

    /**
     * <tt>create_new</tt>必须提供一个干净的对象，要么创建一个新的，要么清理一个旧的。
     *
     */
    virtual void
    create_new(const unsigned int run_no) = 0;

    /**
     * 获取参数并运行任何必要的动作。
     *
     */
    virtual void
    run(ParameterHandler &prm) = 0;
  };

  /**
   * 构造函数
   *
   */
  MultipleParameterLoop();

  /**
   * 解构器。声明这个只是为了有一个虚拟的析构器，这更安全，因为我们有虚拟函数。它实际上没有什么了不起的作用。
   *
   */
  virtual ~MultipleParameterLoop() override = default;

  /**
   * 从一个流中读取输入，直到该流返回<tt>eof</tt>条件或错误。第二个参数可以用来表示我们正在读取的文件的名称（如果那是输入流所代表的）；这只在为错误信息创建输出时使用。
   * 如果提供了非空的 @p last_line
   * ，ParameterHandler对象将在遇到 @p last_line 后停止解析行。
   * 这在添加应被手动解析的额外数据时很方便。    如果
   * @p skip_undefined  是  <code>true</code>
   * ，参数处理程序将跳过未定义的部分和条目。这对于部分解析参数文件非常有用，例如，只获得问题的空间维度。默认情况下，所有条目和子节都应该被声明。
   * @note
   * 这是ParameterHandler实现的三个<tt>parse_input</tt>函数中唯一的重载，这个类用新的行为重写了这个函数。这是因为其他两个<tt>parse_input</tt>函数只是重新格式化它们的输入，然后调用这个版本。
   *
   */
  virtual void
  parse_input(std::istream &     input,
              const std::string &filename       = "input file",
              const std::string &last_line      = "",
              const bool         skip_undefined = false) override;

  /**
   * 重载虚函数（比如 ParameterHandler::parse_input,
   * ，它有两套不同的输入参数类型）会导致非重载的函数被隐藏。通过明确使用
   * ParameterHandler::parse_input
   * 的两个变体，然后重载我们关心的那个，来解决这个问题。
   *
   */
  using ParameterHandler::parse_input;

  /**
   * 运行中心循环。
   *
   */
  void
  loop(UserClass &uc);

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * 列表中的一个对象有多个值的条目。
   *
   */
  class Entry
  {
  public:
    /**
     * 声明多条目是什么：一个变体条目（在大括号里<tt>{</tt>,
     * <tt>}</tt>）或一个数组（在双大括号里<tt>{{</tt>,
     * <tt>}</tt>）。
     *
     */
    enum MultipleEntryType
    {
      /**
       * 一个变体条目。
       *
       */
      variant,
      /**
       * 一个数组条目。
       *
       */
      array
    };

    /**
     * 构造函数
     *
     */
    Entry()
      : type(array)
    {}

    /**
     * 用给定的分段路径、名称和值构造一个对象。分割成不同的变体是由<tt>split_different_values</tt>稍后完成的。
     *
     */
    Entry(const std::vector<std::string> &Path,
          const std::string &             Name,
          const std::string &             Value);

    /**
     * 将条目值分割成不同的分支。
     *
     */
    void
    split_different_values();

    /**
     * 变体条目的路径。
     *
     */
    std::vector<std::string> subsection_path;

    /**
     * 条目名称。
     *
     */
    std::string entry_name;

    /**
     * 原始变体值。
     *
     */
    std::string entry_value;

    /**
     * 由输入文件中的内容构建的条目值列表。
     *
     */
    std::vector<std::string> different_values;

    /**
     * 存储该条目是一个变体条目还是一个数组。
     *
     */
    MultipleEntryType type;

    /**
     * 确定这个对象的内存消耗（以字节为单位）的估计值。
     *
     */
    std::size_t
    memory_consumption() const;
  };

  /**
   * 变体条目值的列表。
   *
   */
  std::vector<Entry> multiple_choices;

  /**
   * 从变体的不同组合中构建的分支数量。这显然等于要执行的运行次数。
   *
   */
  unsigned int n_branches;

  /**
   * 初始化不同的分支，即构建组合。
   *
   */
  void
  init_branches();

  /**
   * 遍历当前由enter_subsection()/leave_subsection()设置的部分，看看哪些条目是变体或数组条目。然后用这些信息填充multiple_choices变量。
   *
   */
  void
  init_branches_current_section();

  /**
   * 将一次运行的条目值转移到条目树上。
   *
   */
  void
  fill_entry_values(const unsigned int run_no);
};


// ---------------------- inline and template functions --------------------
template <class Archive>
inline void
ParameterHandler::save(Archive &ar, const unsigned int) const
{
  // Forward to serialization
  // function in the base class.
  ar &static_cast<const Subscriptor &>(*this);

  ar &*entries.get();

  std::vector<std::string> descriptions;

  for (const auto &pattern : patterns)
    descriptions.push_back(pattern->description());

  ar &descriptions;
}


template <class Archive>
inline void
ParameterHandler::load(Archive &ar, const unsigned int)
{
  // Forward to serialization
  // function in the base class.
  ar &static_cast<Subscriptor &>(*this);

  ar &*entries.get();

  std::vector<std::string> descriptions;
  ar &                     descriptions;

  patterns.clear();
  for (const auto &description : descriptions)
    patterns.push_back(Patterns::pattern_factory(description));
}


template <class ParameterType>
void
ParameterHandler::add_parameter(const std::string &          entry,
                                ParameterType &              parameter,
                                const std::string &          documentation,
                                const Patterns::PatternBase &pattern,
                                const bool                   has_to_be_set)
{
  static_assert(std::is_const<ParameterType>::value == false,
                "You tried to add a parameter using a type "
                "that is const. Use a non-const type.");

  declare_entry(entry,
                Patterns::Tools::Convert<ParameterType>::to_string(parameter,
                                                                   pattern),
                pattern,
                documentation,
                has_to_be_set);

  std::string        path = get_current_full_path(entry);
  const unsigned int pattern_index =
    entries->get<unsigned int>(path + path_separator + "pattern");

  auto action = [&, pattern_index](const std::string &val) {
    parameter = Patterns::Tools::Convert<ParameterType>::to_value(
      val, *patterns[pattern_index]);
  };
  add_action(entry, action);
}

DEAL_II_NAMESPACE_CLOSE

#endif


