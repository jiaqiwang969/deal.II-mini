//include/deal.II-translator/base/tensor_function_parser_0.txt
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

#ifndef dealii_tensor_function_parser_h
#define dealii_tensor_function_parser_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_function.h>
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
 * 这个类实现了一个张量函数对象，它通过解析描述这个函数的字符串来获得其值。它是muparser库的一个封装类（见http://muparser.beltoforion.de/）。这个类本质上是FunctionParser类的一个扩展，用来读入一个TensorFunction。该类读入一个长度为dim<sup>rank</sup>的表达式（用分号隔开），其中张量函数的分量按照C++惯例填充（最快的索引是最右边的那个）。
 *
 *
 * @note
 * 与FunctionParser类相比，TensorFunctionParser类不支持自动区分。
 * 该类使用的一个最小的例子是。
 *
 * @code
 * // set up time dependent tensor function:
 * const std::string variables = "x,y,t";
 * const std::string expression =
 *     "exp(-t)*cos(x+y);-sin(pi*x*y-t);sin(pi*x*y-t);exp(t)*cos(x+y)";
 * std::map<std::string,double> constants;
 * constants["pi"] = numbers::PI;
 *
 * // TensorFunctionParser with 2+1 variables (space + time) in 2D of rank 2.
 * // It is necessary to tell the parser that there is an additional variable
 * // to be taken into account (t).
 * TensorFunctionParser<2,2> tfp;
 * tfp.initialize(variables,
 *             expression,
 *             constants,
 *             true); // flag for time dependence
 *
 * // Point at which we want to evaluate the function
 * Point<2> point(0.0, 1.0);
 *
 * // evaluate the expression at 'point':
 * double result = tfp.value(point);
 *
 * deallog << "Function '" << expression << "'"
 *       << " @ " << point
 *       << " is: "
 *       << std::endl
 *       << result[0][0] << " " << result[0][1] << std::endl
 *       << result[1][0] << " " << result[1][1]
 *       << std::endl;
 * @endcode
 *
 * 也请看FunctionParser类的文档。
 * 这个类重载了TensorFunction基类的虚拟方法value()和value_list()，其表达式的字节编译版本被赋予initialize()方法。请注意，除非你首先调用initialize()方法，该方法接受函数的文本描述作为参数（除其他外），否则该类将无法工作。
 * 描述一个函数的语法遵循通常的编程惯例，在底层muparser库的主页上有详细的解释：http://muparser.beltoforion.de/
 * 。
 *
 * 矢量值函数可以使用字符串来声明，其中函数成分用分号隔开，或者使用字符串的矢量，每个定义一个矢量成分。
 *
 *
 *
 * @ingroup functions
 *
 */
template <int rank, int dim, typename Number = double>
class TensorFunctionParser : public TensorFunction<rank, dim, Number>
{
public:
  /**
   * 标准构造函数。只设置初始时间。这个对象在使用前需要用initialize()方法进行初始化。如果在调用initialize()方法之前试图使用这个函数，那么就会产生一个异常。
   *
   */
  TensorFunctionParser(const double initial_time = 0.0);

  /**
   * 解析函数的构造器。这个对象在使用前需要用initialize()方法进行初始化。如果在调用initialize()方法之前试图使用这个函数，那么就会抛出一个异常。
   * 接收一个分号分隔的表达式列表（张量函数的每个分量一个），一个可选的逗号分隔的常量列表。
   *
   */
  TensorFunctionParser(
    const std::string &expression,
    const std::string &constants      = "",
    const std::string &variable_names = default_variable_names() + ",t");

  /**
   * 复制构造函数。这种类型的对象不能被复制，因此这个构造函数被删除。
   *
   */
  TensorFunctionParser(const TensorFunctionParser &) = delete;

  /**
   * 移动构造函数。该类型的对象不能被移动，因此该构造函数被删除。
   *
   */
  TensorFunctionParser(TensorFunctionParser &&) = delete;

  /**
   * 解构器。
   *
   */
  virtual ~TensorFunctionParser() override;

  /**
   * 复制操作符。该类型的对象不能被复制，因此该操作符被删除。
   *
   */
  TensorFunctionParser &
  operator=(const TensorFunctionParser &) = delete;

  /**
   * 移动操作符。此类型的对象不能被移动，因此该操作符被删除。
   *
   */
  TensorFunctionParser &
  operator=(TensorFunctionParser &&) = delete;

  /**
   * 常量图的类型。由initialize()方法使用。
   *
   */
  using ConstMap = std::map<std::string, double>;


  /**
   * 初始化张量函数。 这个方法接受以下参数。
   * @param[in]  vars
   * 一个包含变量的字符串，这些变量将被待评估的表达式使用。注意，变量可以有任何名字（当然与上面定义的函数名不同！），但顺序很重要。第一个变量将对应于函数被评估的点的第一个分量，第二个变量对应于第二个分量，以此类推。如果这个函数也是依赖于时间的，那么有必要通过将
   * <code>time_dependent</code> 参数设置为真来指定它。
   * 如果这里指定的变量数量与dim（如果这个函数不依赖于时间）或与dim+1（如果它依赖于时间）不同，就会出现异常。
   * @param[in]  expressions
   * 一个字符串向量，包含将被内部解析器（TensorFunctionParser）字节编译的表达式。请注意，这个向量的大小必须与构造函数中声明的TensorFunctionParser的组件数量完全匹配。如果不是这样，就会抛出一个异常。
   * @param[in]  常量
   * 一个常量图，用于传递我们想在表达式中指定的任何必要常量（在上面的例子中是数字pi）。当且仅当一个表达式只包含定义的变量和定义的常量（除了上面指定的函数之外）时，它才是有效的。如果给定的常量的名称是无效的（例如：
   * <code>constants["sin"] = 1.5;</code> ），就会抛出一个异常。
   * @param[in]  time_dependent
   * 如果这是一个依赖于时间的函数，那么在<b>vars</b>中声明的最后一个变量被假定为时间变量，并且在评估该函数时使用this->get_time()来初始化它。当然，在这种情况下，initialize()方法解析的变量数量是dim+1。这个参数的值默认为false，也就是说，不考虑时间。
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
   * 返回张量函数在给定点的值。
   *
   */
  virtual Tensor<rank, dim, Number>
  value(const Point<dim> &p) const override;

  /**
   * 返回张量函数在给定点的值。
   *
   */
  virtual void
  value_list(const std::vector<Point<dim>> &         p,
             std::vector<Tensor<rank, dim, Number>> &values) const override;

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
    tfp;

  /**
   * 一个数组来记录所有的常量，需要在每个线程中初始化tfp。
   *
   */
  std::map<std::string, double> constants;

  /**
   * 一个用于记录变量名称的数组，需要在每个线程中初始化tfp。
   *
   */
  std::vector<std::string> var_names;

  /**
   * 在当前线程中初始化tfp和vars。这个函数在每个线程中只能被调用一次。一个线程可以通过测试'tfp.get().size()==0'（未初始化）或>0（已初始化）来测试该函数是否已经被调用。
   *
   */
  void
  init_muparser() const;
#endif

  /**
   * 一个函数表达式的数组（每个组件一个），需要在每个线程中初始化tfp。
   *
   */
  std::vector<std::string> expressions;

  /**
   * 可用性的状态。这个变量在每次函数被调用评估时被检查。它在initialize()方法中被设置为true。
   *
   */
  bool initialized;

  /**
   * 变量的数量。如果这也是一个时间的函数，那么变量的数量是dim+1，否则就是dim。如果这是一个依赖于时间的函数，那么时间应该是最后一个变量。如果#n_vars与initialize()方法解析的变量数量不一致，则会抛出一个异常。
   *
   */
  unsigned int n_vars;

  /**
   * 组件的数量等于dim<sup>rank</sup>。
   *
   */
  unsigned int n_components;
};


template <int rank, int dim, typename Number>
std::string
TensorFunctionParser<rank, dim, Number>::default_variable_names()
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


