//include/deal.II-translator/base/parameter_acceptor_0.txt
//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2021 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

#ifndef dealii_base_parameter_acceptor_h
#define dealii_base_parameter_acceptor_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/smartpointer.h>

#include <boost/signals2/signal.hpp>

#include <typeinfo>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个参数接受器基类。这个类用于定义一个公共接口，供那些想使用一个全局ParameterHandler来处理参数的类使用。这个类声明了一个静态的ParameterHandler，以及两个静态函数（declaration_all_parameters()和parse_all_parameters()）来管理所有的派生类。
 * 基本接口提供了两种订阅机制：一个*全局订阅机制**和一个*本地订阅机制**。
 * 全局订阅机制是这样的：每当一个从ParameterAcceptor派生的类的对象被创建，那么这个派生类型的对象的指针就会被注册，同时还有参数文件中的路径。
 * 这种注册表在调用单一函数
 * ParameterAcceptor::initialize("file.prm")
 * 时被遍历，该函数依次为每个注册的类调用方法
 * ParameterAcceptor::declare_parameters()
 * ，读取文件`file.prm`，随后为每个注册的类再次调用方法
 * ParameterAcceptor::parse_parameters(),
 * 。方法log_info()可以用来提取从ParameterAcceptor派生的类的信息，这些信息将在调用
 * ParameterAcceptor::initialize(). 时被解析。
 * ParameterAcceptor可以通过三种不同的方式使用：通过重载
 * ParameterAcceptor::declare_parameters() 和
 * ParameterAcceptor::parse_parameters()
 * 方法，为我们想要的每个参数调用其
 * ParameterAcceptor::add_parameter()
 * 方法，或者用你自己的类构造一个ParameterAcceptorProxy类，只要你的类实现了
 * @p declare_parameters 和 @p parse_parameters
 * 函数（这种情况下第一个可以是一个静态成员）。
 * 通过使用add_parameter方法，ParameterAcceptor确保给定的参数被注册在全局参数处理程序中（通过调用
 * ParameterHandler::add_parameter()), 的正确路径。如果你使用
 * ParameterAcceptor::add_parameter()
 * 方法定义所有的参数，那么你就不需要重载这个类的任何虚拟方法。
 * 如果需要对解析后的值进行一些后处理，用户可以给
 * ParameterAcceptor::declare_parameters_call_back 和
 * ParameterAcceptor::parse_parameters_call_back,
 * 附加一个信号，这些信号就在每个派生类的
 * declare_parameters() 和 parse_parameters() 函数之后调用。  step-69
 * 有一个这样做的例子。
 * 这个类的一个典型用法是这样的。
 *
 *
 * @code
 * // This is your own class, derived from ParameterAcceptor
 * class MyClass : public ParameterAcceptor
 * {
 * // The constructor of ParameterAcceptor requires a std::string,
 * // which defines the section name where the parameters of MyClass
 * // will be stored.
 * MyClass()
 *   : ParameterAcceptor("Some class name")
 * {
 *   add_parameter("A param", member_var);
 * }
 *
 * private:
 * std::vector<unsigned int> member_var;
 * ...
 * };
 *
 * int main()
 * {
 * // Make sure you create your object BEFORE calling
 * // ParameterAcceptor::initialize()
 * MyClass class;
 *
 * // With this call, all derived classes will have their
 * // parameters initialized
 * ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 *
 * 一个使用用户定义的声明和解析函数的实现是由下面的例子给出的。
 *
 *
 * @code
 * // Again your own class, derived from ParameterAcceptor
 * //
 * // If you don't pass anything to the constructor of
 * // ParameterAcceptor, then the class name is used, "MyClass"
 * // in this case
 * class MyClass : public ParameterAcceptor
 * {
 * virtual void declare_parameters(ParameterHandler &prm)
 * {
 *   ...
 * }
 *
 * virtual void parse_parameters(ParameterHandler &prm)
 * {
 *   ...
 * }
 * };
 *
 * int main()
 * {
 * // Make sure you create your object BEFORE calling
 * // ParameterAcceptor::initialize()
 * MyClass class;
 * ParameterAcceptor::initialize("file.prm");
 * class.run();
 * }
 * @endcode
 *
 *
 * 参数文件可以被组织成节/分节/子节。要做到这一点，在派生类的构造函数中传递给ParameterAcceptor的
 * std::string 需要包含分隔符"/"。事实上，"first/second/third/My
 * Class "将组织参数如下
 *
 *
 * @code
 * subsection first
 * subsection second
 *   subsection third
 *     subsection My Class
 *      ... # all the parameters
 *     end
 *   end
 * end
 * end
 * @endcode
 *
 * 在下面的例子中，我们提出了一些复杂度越来越高的用例。
 * MyClass是从ParameterAcceptor派生出来的，并且有一个成员对象本身是从ParameterAcceptor派生的。
 *
 * @code
 * class MyClass : public ParameterAcceptor
 * {
 * MyClass (std::string name);
 * virtual void declare_parameters(ParameterHandler &prm);
 * private:
 * SomeParsedClass<dim> my_subclass;
 * ...
 * };
 *
 * MyClass::MyClass(std::string name)
 * : ParameterAcceptor(name)
 * , my_subclass("Forcing term")
 * {}
 *
 * void MyClass::declare_parmeters(ParameterHandler &prm)
 * {
 * // many add_parameter(...);
 * }
 *
 * ...
 *
 * int main()
 * {
 * MyClass mc("My Class");
 * ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 *
 * 在这种情况下，参数的结构将是
 *
 * @code
 * subsection Forcing term
 * ... #parameters of SomeParsedClass
 * end
 * subsection My class
 * ... #all the parameters of MyClass defined in declare_parameters
 * end
 * @endcode
 *
 * 现在假设在主文件中我们需要两个或多个MyClass的对象
 *
 * @code
 * int main()
 * {
 * MyClass ca("Class A");
 * MyClass cb("Class B");
 * ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 *
 * 我们将在参数文件中读取的内容看起来像
 * @code
 * subsection Class A
 * ...
 * end
 * subsection Class B
 * ...
 * end
 * subsection Forcing term
 * ...
 * end
 * @endcode
 * 注意，只有一个 "强制术语
 * "部分，这是因为两个对象都为其SomeParsedClass的部分定义了相同的名称。有两种策略可以改变这种行为。第一个策略（不推荐）是改变SomeParsedClass部分的名称，使其也包含传递给MyClass构造函数的字符串。
 *
 * @code
 * MyClass::MyClass(std::string name)
 * : ParameterAcceptor(name)
 * , my_subclass(name+"
 *
 * --- forcing term")
 * {}
 * @endcode
 *
 * 另一种方法（推荐）是在主类中利用/section/subsection方法**。
 *
 * @code
 * int main()
 * {
 * MyClass ca("/Class A/Class");
 * MyClass cb("/Class B/Class");
 * ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 * 现在，在参数文件中我们可以找到
 *
 * @code
 * subsection Class A
 * subsection Class
 * ...
 * end
 * subsection Forcing term
 * ...
 * end
 * end
 * subsection Class B
 * subsection Class
 * ...
 * end
 * subsection Forcing term
 * ...
 * end
 * end
 * @endcode
 *
 * 注意字符串名称开头的"/"。这被ParameterAcceptor解释为像Unix系统中的根文件夹。"A类
 * "和 "B类
 * "部分不会嵌套在任何部分之下。另一方面，如果字符串不以"/"开头，就像前面的情况一样，将在当前路径下创建**节，这取决于先前定义的节/子节/分节。事实上，"强制条款
 * "部分被嵌套在 "A类 "或 "B类
 * "之下。为了使事情更清楚，让我们考虑以下两个例子
 *
 *
 * @code
 * int main()
 * {
 * MyClass ca("/Class A/Class");
 * MyClass cb("Class B/Class");
 * ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 * 参数文件将有以下结构
 * @code
 * subsection Class A
 * subsection Class
 * ...
 * end
 * subsection Forcing term
 * ...
 * end
 * subsection Class B
 *   subsection Class
 *   ...
 *   end
 *   subsection Forcing term
 *   ...
 *   end
 * end
 * end
 * @endcode
 *
 * 如果代替其中一个路径以"/"结尾，而不只是一个类的名称，那么后续的类将把它解释为一个完整的路径，把类的名称解释为一个目录名称。
 *
 * @code
 * int main()
 * {
 * MyClass ca("/Class A/Class/");
 * MyClass cb("Class B/Class");
 * ParameterAcceptor::initialize("file.prm");
 * }
 * @endcode
 * 参数文件将有以下结构
 *
 * @code
 * subsection Class A
 * subsection Class
 *    ...
 *    subsection Forcing term
 *    ...
 *    end
 *    subsection Class B
 *        subsection Class
 *        ...
 *        end
 *        subsection Forcing term
 *        ...
 *        end
 *    end
 * end
 * end
 * @endcode
 *
 * 作为最后的评论，为了允许对所有的节/子节进行适当的管理，对象的实例化和对
 * ParameterAcceptor::initialize()
 * 的调用不能在多个同时运行的线程中进行。
 * 如果你传递一个空的名字， boost::core::demangle()
 * 函数会被用来用类名本身的人类可读版本来填充节名。
 * 关于如何使用这个类，请看教程程序 step-60 中的例子。
 *
 *
 */
class ParameterAcceptor : public Subscriptor
{
public:
  /**
   * 构造函数将派生类添加到接受器的列表中。如果指定了一个部分的名称，那么这将用于在给定的部分范围内的参数，否则将使用派生类的一个漂亮的打印版本。
   *
   */
  ParameterAcceptor(const std::string &section_name = "");

  /**
   * 解构器。
   *
   */
  virtual ~ParameterAcceptor() override;

  /**
   * 调用
   * declare_all_parameters()，从`filename`中读取参数（只有当`filename`是一个非空的字符串时），然后调用parse_all_parameters()。
   * 如果参数`filename`是空字符串，那么就不会尝试读取参数文件。如果你可以使用默认值，并且不想读取外部文件来使用从ParameterAcceptor派生的类，这可能很有用。
   * 如果 @p output_filename
   * 不是空字符串，那么我们就把读到的内容写到 @p
   * output_filename 文件中，使用 @p output_style_for_output_filename.
   * 中指定的样式
   * 输入和输出文件的格式都是用文件本身的扩展名来选择。对于
   * @p filename, 可以是`prm`、`xml`或`json`，对于 @p output_filename.
   * 可以是任何支持的格式。 如果输入文件不存在，将按照
   * @p output_style_for_filename,
   * 中指定的样式为你创建一个同名的默认文件，并抛出一个异常。
   * 默认情况下，用于写入文件的文件格式是由文件名的扩展名推断出来的。如果相应的
   * ParameterHandler::OutputStyle
   * 指定了一个格式规范，这必须与文件扩展名兼容，否则将抛出一个异常。
   * 如果扩展名不被识别，并且你没有在相应的
   * ParameterHandler::OutputStyle,
   * 中指定一个格式，则会抛出一个断言。      @param  filename
   * 输入文件名  @param  output_filename 输出文件名  @param
   * output_style_for_output_filename 如何写入输出文件  @param  prm
   * 要使用的ParameterHandler  @param  output_style_for_filename
   * 如果默认的输入文件不存在如何写入
   *
   */
  static void
  initialize(const std::string &filename        = "",
             const std::string &output_filename = "",
             const ParameterHandler::OutputStyle
                                                 output_style_for_output_filename = ParameterHandler::Short,
             ParameterHandler &                  prm = ParameterAcceptor::prm,
             const ParameterHandler::OutputStyle output_style_for_filename =
               ParameterHandler::DefaultStyle);

  /**
   * 调用
   * declare_all_parameters()，从`input_stream`中读取`prm`格式的参数，然后调用parse_all_parameters()。
   * 如果 "input_stream "无效，则抛出一个异常。      @param
   * input_stream 输入流  @param  prm 要使用的ParameterHandler。
   *
   */
  static void
  initialize(std::istream &    input_stream,
             ParameterHandler &prm = ParameterAcceptor::prm);


  /**
   * 清除类列表和全局参数文件。
   *
   */
  static void
  clear();

  /**
   * 派生类可以使用这个方法来声明他们的参数。
   * ParameterAcceptor::initialize()
   * 为每个派生类调用它。默认实现是空的。
   *
   */
  virtual void
  declare_parameters(ParameterHandler &prm);

  /**
   * 申报参数的回调。这个信号在declare_parameters()被调用后立即被触发，以允许用户在参数被声明后立即准备他们的变量。默认实现是空的。
   *
   */
  boost::signals2::signal<void()> declare_parameters_call_back;

  /**
   * 派生类可以使用这个方法来解析他们的参数。
   * ParameterAcceptor::initialize()
   * 为每个派生类调用它。默认的实现是空的。
   *
   */
  virtual void
  parse_parameters(ParameterHandler &prm);

  /**
   * 解析参数的回调。这个函数在parse_parameters()的末尾被调用，以允许用户在参数被解析后立即处理他们的参数。默认实现是空的。
   * 你可以使用这个函数，例如，在你从参数文件中读取你想使用的正交点的数量后，创建一个正交规则。
   *
   */
  boost::signals2::signal<void()> parse_parameters_call_back;

  /**
   * 解析给定的ParameterHandler。这个函数为每个派生类进入get_section_name()返回的小节，并解析所有使用add_parameter()添加的参数。
   *
   */
  static void
  parse_all_parameters(ParameterHandler &prm = ParameterAcceptor::prm);

  /**
   * 用所有派生类的参数初始化全局ParameterHandler.这个函数进入每个派生类的get_section_name()返回的小节，并声明使用add_parameter()添加的所有参数。
   *
   */
  static void
  declare_all_parameters(ParameterHandler &prm = ParameterAcceptor::prm);

  /**
   * 返回这个类的分节名称。如果在构造时提供了一个名称，那么就返回这个名称，否则就返回这个类的拆分名称。
   *
   */
  std::string
  get_section_name() const;

  /**
   * 遍历所有注册的类，并找出我们需要进入的子段。
   *
   */
  std::vector<std::string>
  get_section_path() const;

  /**
   * 在正确的路径上添加一个参数。这个方法在输入正确的章节路径后，将所有参数转发给prm.add_parameter()方法。
   * 默认情况下，它使用 ParameterAcceptor::prm
   * 变量作为ParameterHandler。    更多信息请参见
   * ParameterHandler::add_parameter() 的文档。
   *
   */
  template <class ParameterType>
  void
  add_parameter(const std::string &          entry,
                ParameterType &              parameter,
                const std::string &          documentation = "",
                ParameterHandler &           prm_          = prm,
                const Patterns::PatternBase &pattern =
                  *Patterns::Tools::Convert<ParameterType>::to_pattern());

  /**
   * 全局参数处理程序。
   *
   */
  static ParameterHandler prm;

  /**
   * 将给定的 @p subsection
   * 添加到存储在这个类中的全局路径中。
   * 这个函数改变了enter_my_subsection()的行为，将一个新的小节追加到存储在这个类中的路径上。
   * 这个方法可以用来把这个类的参数分成若干小节，同时仍然保持这个类的一般行为。
   * 下面的片段给出了一个用法示例。
   * @code
   * class MyClass : public ParameterAcceptor
   * {
   * MyClass()
   *   : ParameterAcceptor("Main section")
   * {
   *   add_parameter("A param", member_var);
   *   enter_subsection("New section");
   *   add_parameter("Another param", another_member_var);
   *   leave_subsection();
   * }
   *
   * private:
   * std::vector<unsigned int> member_var = {1,2};
   * std::map<types::boundary_id, std::string> another_member_var;
   * ...
   * };
   *
   * int main()
   * {
   * // ParameterAcceptor::initialize()
   * MyClass class;
   *
   * // With this call, all derived classes will have their
   * // parameters initialized
   * ParameterAcceptor::initialize("file.prm");
   * }
   * @endcode
   * 这将产生一个参数文件，组织为
   * @code
   * subsection Main section
   * set A param = 1, 2
   * subsection New section
   *   set Another param =
   * end
   * end
   * @endcode
   *
   *
   */
  void
  enter_subsection(const std::string &subsection);

  /**
   * 留下通过调用enter_subsection()函数输入的分节。
   *
   */
  void
  leave_subsection();

  /**
   * 确保我们输入的是给定参数的正确分节。
   *
   */
  void
  enter_my_subsection(ParameterHandler &prm);

  /**
   * 这个函数撤销了enter_my_subsection()函数所做的事情。只有在enter_my_subsection()在这个函数之前被调用到`prm`时才有意义。
   *
   */
  void
  leave_my_subsection(ParameterHandler &prm);

private:
  /**
   * 一个包含所有ParameterAcceptor类型的构造类的列表。
   *
   */
  static std::vector<SmartPointer<ParameterAcceptor>> class_list;

   /** The index of this specific class within the class list. */ 
  const unsigned int acceptor_id;

  /**
   * 部分之间的分隔符。
   *
   */
  static const char sep = '/';

protected:
   /** The subsection name for this class. */ 
  const std::string section_name;

   /** The subsubsections that are currently active. */ 
  std::vector<std::string> subsections;
};



/**
 * 一个代理ParameterAcceptor的包装器，用于具有静态成员函数
 * @p declare_parameters, 和非虚拟 @p parse_parameters 方法的类。
 * 如果你不能或不想从ParameterAcceptor派生出你的 "参数接受
 * "类，例如，根据设计，你必须有一个静态成员函数 @p
 * declare_parameters 和一个成员 @p
 * parse_parameters，或者有人已经为你实现了这样一个类，并且只提供给你一个你无法修改的API，那么你仍然可以使用ParameterAcceptor设施，将你的类包装到ParameterAcceptorProxy。
 * 这个类实现了ParameterAcceptor的公共接口，同时它派生自模板类
 * @p SourceClass, ，允许你将你现有的 @p SourceClass
 * 注册为ParameterAcceptor类，而不需要你明确地从ParameterAcceptor派生你的
 * @p SourceClass 。
 * 下面的代码片段给出了一个使用方法的例子，使用
 * Functions::ParsedFunction 作为一个示例源类。
 *
 *
 * @code
 * ParameterAcceptorProxy<Functions::ParsedFunction<2> > fun("Some function");
 * ParameterAcceptor::initialize("test.prm");
 * @endcode
 *
 * 上面的代码片段将用 "某个函数 "部分初始化
 * ParameterAcceptor::prm
 * ，并将正确地解析和分配从文件`test.prm`中解析出的表达式给对象`fun`。如果不存在，程序将退出，并为你生成（这里你可以看到用上述片段生成的参数文件的结果的简短文本）。
 *
 *
 * @code
 * # Parameter file generated with
 * # DEAL_II_PACKAGE_VERSION = 9.0.0-pre
 * subsection Some function
 * set Function constants  =
 * set Function expression = 0
 * set Variable names      = x,y,t
 * end
 * @endcode
 *
 * 产生的`fun`对象，既是一个ParsedFunction对象，也是一个ParameterAcceptor对象，允许你用它来替代ParsedFunction类，自动声明和解析参数文件。
 * 关于如何使用这个类的例子，见教程程序 step-60 。
 *
 *
 */
template <class SourceClass>
class ParameterAcceptorProxy : public SourceClass, public ParameterAcceptor
{
public:
  /**
   * 默认构造函数。参数`section_name`被转发给ParameterAcceptor类的构造函数，而所有其他参数则被传递给SourceClass的构造函数。
   *
   */
  template <typename... Args>
  ParameterAcceptorProxy(const std::string &section_name, Args... args);

  /**
   * 重载 ParameterAcceptor::declare_parameters 函数，通过调用 @p
   * SourceClass::declare_parameters ，以 @p prm 作为参数。
   *
   */
  virtual void
  declare_parameters(ParameterHandler &prm) override;

  /**
   * 通过调用 @p SourceClass::parse_parameters 和 @p prm
   * 作为参数，重载了 ParameterAcceptor::parse_parameters 函数。
   *
   */
  virtual void
  parse_parameters(ParameterHandler &prm) override;
};



// Inline and template functions
template <class ParameterType>
void
ParameterAcceptor::add_parameter(const std::string &          entry,
                                 ParameterType &              parameter,
                                 const std::string &          documentation,
                                 ParameterHandler &           prm,
                                 const Patterns::PatternBase &pattern)
{
  enter_my_subsection(prm);
  prm.add_parameter(entry, parameter, documentation, pattern);
  leave_my_subsection(prm);
}



template <class SourceClass>
template <typename... Args>
ParameterAcceptorProxy<SourceClass>::ParameterAcceptorProxy(
  const std::string &section_name,
  Args... args)
  : SourceClass(args...)
  , ParameterAcceptor(section_name)
{}



template <class SourceClass>
void
ParameterAcceptorProxy<SourceClass>::declare_parameters(ParameterHandler &prm)
{
  SourceClass::declare_parameters(prm);
}



template <class SourceClass>
void
ParameterAcceptorProxy<SourceClass>::parse_parameters(ParameterHandler &prm)
{
  SourceClass::parse_parameters(prm);
}

DEAL_II_NAMESPACE_CLOSE

#endif


