//include/deal.II-translator/base/function_parser_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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

#ifndef dealii_function_parser_h
#define dealii_function_parser_h


#include <deal.II/base/config.h>

#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/thread_local_storage.h>

#include <map>
#include <memory>
#include <vector>

namespace mu
{
  class Parser;
}

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <typename>
class Vector;
#endif

/**
 * 这个类实现了一个函数对象，它通过解析描述这个函数的字符串来获得其值。它是muparser库的一个封装类（见https://beltoforion.de/en/muparser/）。这个类可以让你对
 * "x "和 "y "的给定值评估字符串，如
 * "sqrt（1-x^2+y^2）"。更多信息请参考muparser文档。 这个类在
 * step-33 和 step-36 的教程程序中使用（后者更容易理解）。
 * 除了muparser的内置函数外，即
 *
 * @code
 * sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh,
 * atan2, log2, log10, log, ln, exp, sqrt, sign, rint, abs, min, max, sum, avg
 * @endcode
 * 这个类还支持以下操作。
 *
 *
 *
 * -  <code>if(condition, then-value, else-value)</code>
 *
 *
 * -  <code>|</code> and <code>&</code>  (逻辑或和和)
 *
 *
 *
 * -  <code>int(x)</code>, <code>ceil(x)</code>, <code>floor(x)</code>  (四舍五入)
 *
 *
 *
 * -  <code>cot(x)</code>, <code>csc(x)</code>, <code>sec(x)</code>
 *
 *
 * -  <code>pow(x,n)</code>, <code>log(x)</code>
 *
 * -  <code>erfc(x)</code>
 *
 * -  <code>rand()</code>, <code>rand_seed(seed)</code>
 *
 * @note
 * 这个类通过扩展muparser，将刚才提到的函数列表实现为用户定义的函数。这特别意味着，`if(condition,
 * then-value,
 * else-value)`语法在确定条件是否为真之前，会评估所有三个参数，然后抛弃
 * "then "或 "else
 * "表达式。在几乎所有的情况下，这都不是问题，除非其中一个表达式的求值抛出了一个浮点异常，而这个异常随后会被丢弃。假设浮点异常被打开了，在大多数系统上，deal.II在调试模式下是默认的）。一个例子是表达式`if(x>0,
 * sqrt(x),
 * 0)`，它在数学上有很好的定义，但是在启用了这个功能的系统上，当评估到负的`x'时，程序会以浮点异常中止。这是因为`x`的平方根在`if`语句的条件被考虑之前就被计算了，以确定结果应该是第二个还是第三个参数。如果这种行为是一个问题，你可以求助于muparser内置的语法`(condition
 * ? then-value :
 * else-value)`，使用C++程序员熟悉的三元组语法。如果使用这种语法，muparser就会使用懒惰评估，其中只有一个分支被评估，取决于`condition`是真还是假。
 * 下面的例子展示了如何使用这个类。
 *
 * @code
 * // set up problem:
 * std::string variables = "x,y";
 * std::string expression = "cos(x) + sqrt(y)";
 * std::map<std::string, double> constants;
 *
 * // FunctionParser with 2 variables and 1 component:
 * FunctionParser<2> fp(1);
 * fp.initialize(variables,
 *             expression,
 *             constants);
 *
 * // Point at which we want to evaluate the function
 * Point<2> point(0.0, 4.0);
 *
 * // evaluate the expression at 'point':
 * double result = fp.value(point);
 *
 * deallog << "Function '" << expression << "'"
 *       << " @ " << point
 *       << " is " << result << std::endl;
 * @endcode
 * 第二个例子要复杂一些。
 *
 * @code
 * // Define some constants that will be used by the function parser
 * std::map<std::string, double> constants;
 * constants["pi"] = numbers::PI;
 *
 * // Define the variables that will be used inside the expressions
 * std::string variables = "x,y,z";
 *
 * // Define the expressions of the individual components of a
 * // vector valued function with two components:
 * std::vector<std::string> expressions(2);
 * expressions[0] = "sin(2*pi*x)+sinh(pi*z)";
 * expressions[1] = "sin(2*pi*y)*exp(x^2)";
 *
 * // function parser with 3 variables and 2 components
 * FunctionParser<3> vector_function(2);
 *
 * // And populate it with the newly created objects.
 * vector_function.initialize(variables,
 *                          expressions,
 *                          constants);
 *
 * // Point at which we want to evaluate the function
 * Point<3> point(0.0, 1.0, 1.0);
 *
 * // This Vector will store the result
 * Vector<double> result(2);
 *
 * // Fill 'result' by evaluating the function
 * vector_function.vector_value(point, result);
 *
 * // We can also only evaluate the 2nd component:
 * const double c = vector_function.value(point, 1);
 *
 * // Output the evaluated function
 * deallog << "Function '" << expressions[0] << "," << expressions[1] << "'"
 *       << " at " << point
 *       << " is " << result << std::endl;
 * @endcode
 *
 * 这个类重载了Function基类的虚拟方法value()和vector_value()，并使用了给初始化()方法的表达式的字节编译版本。请注意，除非你首先调用initialize()方法，该方法接受函数的文本描述作为参数（除其他外），否则该类将无法工作。
 * 描述一个函数的语法遵循通常的编程惯例，在底层muparser库的主页上有详细的解释：https://beltoforion.de/en/muparser/。
 * 关于支持ParameterHandler的FunctionParser类的包装器，见
 * Functions::ParsedFunction.  。
 * 矢量值函数可以使用字符串来声明，其中的函数成分用分号隔开，或者使用字符串的矢量，每个定义一个矢量成分。
 * 一个与时间有关的标量函数的例子如下。
 * @code
 *  // Empty constants object
 *  std::map<std::string,double> constants;
 *
 *  // Variables that will be used inside the expressions
 *  std::string variables = "x,y,t";
 *
 *  // Define the expression of the scalar time dependent function.
 *  std::string expression = "exp(y*x)*exp(-t)";
 *
 *  // Generate an empty scalar function
 *  FunctionParser<2> function;
 *
 *  // And populate it with the newly created objects.
 *  function.initialize(variables,
 *                      expression,
 *                      constants,
 * // Treat the last variable ("t") as time.
 *                      true);
 * @endcode
 *
 * 下面是另一个例子，说明如何通过使用单个字符串来实例化一个矢量值函数。
 *
 * @code
 *  // Empty constants object
 *  std::map<std::string,double> constants;
 *
 *  // Variables that will be used inside the expressions
 *  std::string variables = "x,y";
 *
 *  // Define the expression of the vector valued  function.
 *  std::string expression = "cos(2*pi*x)*y^2; sin(2*pi*x)*exp(y)";
 *
 *  // Generate an empty vector valued function
 *  FunctionParser<2> function(2);
 *
 *  // And populate it with the newly created objects.
 *  function.initialize(variables,
 *                      expression,
 *                      constants);
 * @endcode
 *
 *
 * @note
 * 这个类和SymbolicFunction类的区别是SymbolicFunction类允许计算一阶和二阶导数（以符号的方式），而这个类只计算一阶导数，使用有限差分。对于复杂的表达式，这个类通常比SymbolicFunction快。
 *
 *
 * @ingroup functions
 *
 *
 */
template <int dim>
class FunctionParser : public AutoDerivativeFunction<dim>
{
public:
  /**
   * 构造函数。它的参数与基类Function相同，附加参数 @p h,
   * 用于使用有限差分进行梯度计算。在你使用这个对象之前，需要用initialize()方法对其进行初始化。如果在调用initialize()方法之前试图使用这个函数，那么就会产生一个异常。
   *
   */
  FunctionParser(const unsigned int n_components = 1,
                 const double       initial_time = 0.0,
                 const double       h            = 1e-8);

  /**
   * 解析函数的构造器。直接接受一个分号分隔的表达式列表（函数的每个组成部分都有一个），一个可选的逗号分隔的常数列表，变量名称和通过有限差分计算一阶导数的步骤大小。
   *
   */
  FunctionParser(const std::string &expression,
                 const std::string &constants      = "",
                 const std::string &variable_names = default_variable_names() +
                                                     ",t",
                 const double h = 1e-8);

  /**
   * 复制构造函数。这种类型的对象不能被复制，因此这个构造函数被删除。
   *
   */
  FunctionParser(const FunctionParser &) = delete;

  /**
   * 移动构造函数。该类型的对象不能被移动，因此该构造函数被删除。
   *
   */
  FunctionParser(FunctionParser &&) = delete;

  /**
   * 解构器。
   *
   */
  virtual ~FunctionParser() override;

  /**
   * 复制操作符。该类型的对象不能被复制，因此该操作符被删除。
   *
   */
  FunctionParser &
  operator=(const FunctionParser &) = delete;

  /**
   * 移动操作符。此类型的对象不能被移动，因此该操作符被删除。
   *
   */
  FunctionParser &
  operator=(FunctionParser &&) = delete;

  /**
   * 常量图的类型。由initialize()方法使用。
   *
   */
  using ConstMap = std::map<std::string, double>;

  /**
   * 通过设置实际解析的函数来初始化对象。      @param[in]
   * vars
   * 一个字符串，包含将被评估的表达式所使用的变量。注意，这些变量可以有任何名字（当然与上面定义的函数名不同！），但顺序很重要。第一个变量将对应于函数被评估的点的第一个分量，第二个变量对应于第二个分量，以此类推。如果这个函数也与时间有关，那么有必要通过将
   * <code>time_dependent</code> 参数设置为真来指定它。
   * 如果这里指定的变量数与dim（如果这个函数不依赖于时间）或与dim+1（如果它依赖于时间）不同，就会产生一个异常。
   * @param[in]  expressions
   * 一个包含表达式的字符串列表，这些表达式将被内部分析器（muParser）进行字节编译。请注意，这个向量的大小必须与构造函数中声明的
   * FunctionParser
   * 的组件数量完全匹配。如果不是这样，就会抛出一个异常。
   * @param[in]  常数
   * 一个常数的映射，用来传递我们想在表达式中指定的任何必要的常数（在上面的例子中是数字pi）。当且仅当一个表达式只包含定义的变量和定义的常量（除了上面指定的函数）时，它才是有效的。如果给定的常量的名称是无效的（例如：
   * <code>constants["sin"] = 1.5;</code> ），就会抛出一个异常。
   * @param[in]  time_dependent
   * 如果这是一个依赖时间的函数，那么在 @p vars
   * 中声明的最后一个变量被假定为时间变量，并且
   * FunctionTime::get_time()
   * 被用来在评估函数时初始化它。当然，在这种情况下，initialize()解析的变量数量是dim+1。这个参数的值默认为false，也就是说，不考虑时间。
   *
   */
  void
  initialize(const std::string &             vars,
             const std::vector<std::string> &expressions,
             const ConstMap &                constants,
             const bool                      time_dependent = false);

  /**
   * 初始化函数。和上面一样，但接受一个字符串，而不是一个字符串的向量。如果这是一个向量值的函数，它的组成部分应该用分号来分隔。如果这个方法被调用，而成功解析的组件数与基函数的组件数不一致，就会产生一个异常。
   *
   */
  void
  initialize(const std::string &vars,
             const std::string &expression,
             const ConstMap &   constants,
             const bool         time_dependent = false);

  /**
   * 一个返回变量默认名称的函数，用于initialize()函数的第一个参数：它在1d中返回
   * "x"，在2d中返回 "x,y"，在3d中返回 "x,y,z"。
   *
   */
  static std::string
  default_variable_names();

  /**
   * 返回函数在给定点的值。除非只有一个分量（即函数是标量的），否则你应该说明你想评估的分量；它默认为零，即第一个分量。
   *
   */
  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;

  /**
   * 返回一个向量值函数在给定点的所有分量  @p  p.
   * <code>values</code>  应事先有正确的大小，即#n_components。
   *
   */
  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override;

  /**
   * 返回一个函数表达式的数组（每个组件一个），用于初始化这个函数。
   *
   */
  const std::vector<std::string> &
  get_expressions() const;

  /**
   * @addtogroup  Exceptions 
     * @{ 
   *
   */
  DeclException2(ExcParseError,
                 int,
                 std::string,
                 << "Parsing Error at Column " << arg1
                 << ". The parser said: " << arg2);

  DeclException2(ExcInvalidExpressionSize,
                 int,
                 int,
                 << "The number of components (" << arg1
                 << ") is not equal to the number of expressions (" << arg2
                 << ").");

  //@}

private:
#ifdef DEAL_II_WITH_MUPARSER
  /**
   * 每个线程的变量的位置
   *
   */
  mutable Threads::ThreadLocalStorage<std::vector<double>> vars;

  /**
   * 每个线程的muParser对象（每个组件也有一个）。我们存储的是unique_ptr，这样我们就不需要在这个头中包含
   * mu::Parser 的定义。
   *
   */
  mutable Threads::ThreadLocalStorage<std::vector<std::unique_ptr<mu::Parser>>>
    fp;

  /**
   * 一个数组来记录所有的常量，需要在每个线程中初始化fp。
   *
   */
  std::map<std::string, double> constants;

  /**
   * 一个用于记录变量名称的数组，需要在每个线程中初始化fp。
   *
   */
  std::vector<std::string> var_names;

  /**
   * 在当前线程中初始化fp和vars。这个函数在每个线程中只能被调用一次。一个线程可以通过测试'fp.get().size()==0'（未初始化）或>0（已初始化）来测试该函数是否已经被调用。
   *
   */
  void
  init_muparser() const;
#endif

  /**
   * 一个函数表达式的数组（每个组件一个），需要在每个线程中初始化fp。
   *
   */
  std::vector<std::string> expressions;

  /**
   * 可用性的状态。这个变量在每次函数被调用评估时都会被检查。它在initialize()方法中被设置为true。
   *
   */
  bool initialized;

  /**
   * 变量的数量。如果这也是一个时间的函数，那么变量的数量是dim+1，否则就是dim。如果这是一个依赖于时间的函数，那么时间应该是最后一个变量。如果#n_vars与initialize()方法解析出的变量数量不一致，那么就会产生一个异常。
   *
   */
  unsigned int n_vars;
};


template <int dim>
std::string
FunctionParser<dim>::default_variable_names()
{
  switch (dim)
    {
      case 1:
        return "x";
      case 2:
        return "x,y";
      case 3:
        return "x,y,z";
      default:
        Assert(false, ExcNotImplemented());
    }
  return "";
}



DEAL_II_NAMESPACE_CLOSE

#endif


