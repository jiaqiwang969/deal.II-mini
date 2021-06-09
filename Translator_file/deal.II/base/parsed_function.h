//include/deal.II-translator/base/parsed_function_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2021 by the deal.II authors
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


#ifndef dealii_parsed_function_h
#define dealii_parsed_function_h

#include <deal.II/base/config.h>

#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/parameter_handler.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  /**
   * FunctionParser类的友好接口。这个类是作为FunctionParser类的一个封装器。它被用于
   * step-34 的教程程序中。
   * 它提供了两个方法来声明和解析一个ParameterHandler对象，并创建在参数文件中声明的Function对象。这个类派生于AutoDerivativeFunction类，所以你不需要指定导数。这个类的一个使用例子如下。
   * @code
   * // A parameter handler
   * ParameterHandler prm;
   *
   * // Declare a section for the function we need
   * prm.enter_subsection("My vector function");
   * ParsedFunction<dim>::declare_parameters(prm, dim);
   * prm.leave_subsection();
   *
   * // Create a ParsedFunction
   * ParsedFunction<dim> my_vector_function(dim);
   *
   * // Parse an input file.
   * prm.parse_input(some_input_file);
   *
   * // Initialize the ParsedFunction object with the given file
   * prm.enter_subsection("My vector function");
   * my_vector_function.parse_parameters(prm);
   * prm.leave_subsection();
   *
   * @endcode
   * 下面是一个输入参数的例子（关于函数定义的语法的详细描述，请看FunctionParser类的文档）。
   * @code
   *
   * # A test two dimensional vector function, depending on time
   * subsection My vector function
   * set Function constants  = kappa=.1, lambda=2.
   * set Function expression = if(y>.5, kappa*x*(1-x),0); t^2*cos(lambda*pi*x)
   * set Variable names      = x,y,t
   * end
   *
   * @endcode
   *
   * @ingroup functions
   *
   */
  template <int dim>
  class ParsedFunction : public AutoDerivativeFunction<dim>
  {
  public:
    /**
     * 构建一个向量函数。生成的向量函数有 @p n_components
     * 个分量（默认为1）。参数 @p h
     * 用于初始化AutoDerivativeFunction类，该类来源于此。
     *
     */
    ParsedFunction(const unsigned int n_components = 1, const double h = 1e-8);

    /**
     * 声明本类所需的参数。额外的参数 @p
     * n_components用于根据将解析这个ParameterHandler的函数的组件数量来生成正确的代码。如果被解析的组件数量与此对象的组件数量不一致，就会抛出一个断言，并且程序被中止。
     * 这个类的默认行为是声明以下条目。
     * @code
     *
     * set Function constants  =
     * set Function expression = 0
     * set Variable names      = x,y,t
     *
     * @endcode
     *
     *
     */
    static void
    declare_parameters(ParameterHandler & prm,
                       const unsigned int n_components = 1);

    /**
     * 该类所需的解析参数。
     * 如果被解析的组件数量与此对象的组件数量不一致，就会抛出一个断言，并中止程序。
     * 为了使该类能够正常运行，我们遵循FunctionParser类中声明的相同约定（关于函数声明的语法，请看那里的详细说明）。
     * 可以从参数文件中解析的变量有以下三种。
     * @code
     *
     * set Function constants  =
     * set Function expression =
     * set Variable names      =
     *
     * @endcode
     * %函数常量是以name=value为形式的一对集合，用逗号分隔，例如。
     * @code
     *
     * set Function constants = lambda=1., alpha=2., gamma=3.
     *
     * @endcode
     * 这些常量可以在函数表达式的声明中使用，它遵循FunctionParser类的惯例。
     * 为了指定向量函数，必须用分号来分隔不同的组成部分，例如
     * @code
     *
     * set Function expression = cos(pi*x); cos(pi*y)
     *
     * @endcode
     * 变量名称条目可以用来自定义函数中使用的变量名称。它的默认值是
     * @code
     *
     * set Variable names      = x,t
     *
     * @endcode
     * 对于一维的问题。
     * @code
     *
     * set Variable names      = x,y,t
     *
     * @endcode
     * 用于二维问题和
     * @code
     *
     * set Variable names      = x,y,z,t
     *
     * @endcode
     * 适用于三维问题。
     * 时间变量可以根据FunctionTime基类的规格来设置。
     *
     */
    void
    parse_parameters(ParameterHandler &prm);

    /**
     * 返回一个向量值函数在给定点的所有分量  @p  p。
     *
     */
    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override;

    /**
     * 返回该函数在给定点的值。除非只有一个分量（即函数是标量的），否则你应该说明你想评估的分量；它默认为零，即第一个分量。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 对于时间相关的函数，将时间设置为一个特定的值。
     * 我们需要覆盖这一点，以便在访问器FunctionParser<dim>中也设置时间。
     *
     */
    virtual void
    set_time(const double newtime) override;

  private:
    /**
     * 我们用来做计算的对象。
     *
     */
    FunctionParser<dim> function_object;
  };
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif


