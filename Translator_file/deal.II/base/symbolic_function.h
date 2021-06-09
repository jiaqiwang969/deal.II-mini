//include/deal.II-translator/base/symbolic_function_0.txt
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

#ifndef dealii_symbolic_function_h
#define dealii_symbolic_function_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/sd.h>

#include <functional>
#include <iostream>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename number>
class Vector;
namespace Functions
{
  template <int dim, typename RangeNumberType>
  class SymbolicFunction;
}
#endif

namespace Functions
{
#ifdef DEAL_II_WITH_SYMENGINE
  /**
   * 一个函数类，利用符号微分来计算梯度、拉普拉斯、Hessians和时间导数。
   * 这个类可用于使用 Differentiation::SD
   * 命名空间提供的方法来定义函数。特别是，可以定义一个符号化的评价点（函数的参数），以及一个符号化的表达式。
   * 符号梯度和符号Hessians在构造时计算，当用户使用update_user_substitution_map()方法要求在符号函数中进行替换时，也会计算。
   * 每当一个评估方法被调用时，就会尝试进行替换，用评估点替换坐标符号参数，用get_time()方法返回的当前时间替换符号时间。用户必须确保在求值时，参数替换提供了一个完全求值的表达式（即除了数值之外，函数表达式中不包含其他符号），否则将抛出一个异常。额外的符号可以通过存储在用户提供的替换图中进行部分评估或替换，该替换图可以通过调用update_user_substitution_map（）或set_additional_function_arguments（）方法进行更新。
   * 这个类的最简单的使用情况在下面的例子中给出。
   * @code
   * SymbolicFunction<2> fun("x^2+y; t*x*y");
   * fun.set_time(3.0);
   * Point<2> p(1.0, 2.0);
   *
   * auto a = fun.value(p, / component / 0); // a = 3.0
   * auto b = fun.value(p, / component / 1); // b = 6.0
   *
   * auto df_dt = fun.time_derivative();
   *
   * auto c = df_dt.value(p, / component / 0); // c = 0.0
   * auto d = df_dt.value(p, / component / 1); // d = 2.0
   * @endcode
   * 在这个例子中，一个有两个组件的Function是用一个字符串定义的，这个字符串包含了用分号分隔的表达式。
   * 一个更复杂的例子，明确地使用了
   * Differentiation::SD::Expression 对象，由以下例子给出
   * @code
   * using namespace Differentiation::SD;
   * // Create a position Tensor<1,2,Differentiation::SD::Expression>
   * // with symbols "x" and "y", and the symbol "t"
   * const auto x = SymbolicFunction<2>::get_default_coordinate_symbols();
   * const auto t = make_symbol("t");
   *
   * // Use directly x[0] (the symbol "x"), x[1] (the symbol "y"), and t
   * // (the symbol "t").
   * Expression f = std::sin(x[0])*std::cos(x[1])*std::sin(t);
   * // Alternatively, you can achieve the same result parsing a string:
   * // Expression f("sin(x)*cos(y)*sin(t)", true);
   * SymbolicFunction<2> function({f}, x);
   *
   * // Evaluate the function, its gradient, and its Laplacian
   * Point<2> p(1.0, 2.0);
   * auto fp = function.value(p);
   * auto gradfp = function.gradient(p);
   * auto lapfp = function.laplacian(p);
   *
   * // Evaluate the time derivative of the function, its gradient, and its
   * // Laplacian
   * auto time_derivative = function.time_derivative();
   * auto dt_fp = time_derivative.value(p);
   * auto dt_gradfp = time_derivative.gradient(p);
   * auto dt_lapfp = time_derivative.laplacian(p);
   * @endcode
   * 部分替换是可能的（也就是说，你可以使用额外的符号来定义函数）。然而，一旦你评估该函数，你必须确保所有无关的符号（即那些没有提到空间
   * @p coordinate_symbols 或 @p time_symbol
   * 变量的符号）已经被数值，或空间或时间参数的表达式所替代，方法是调用update_user_substitution_map()或set_additional_function_arguments()方法。
   * 如果你的函数需要评估额外的参数，你可以通过调用set_additional_function_arguments()方法指定它们。
   * 如果你用相同的参数调用update_user_substitution_map()和set_additional_function_arguments()，对函数评估的影响将是相同的，然而，内部行为和函数导数将不同。update_user_substitution_map()方法执行一次替换（第一次需要替换），然后在内部存储一个结果表达式的副本，以及它的导数（如果需要）。然后在随后的所有评估中使用这些内容。调用set_additional_function_arguments()将在评估时间内，在*所有导数被计算之后，对传递的替换图进行实时评估。
   * @note
   * 这个类和FunctionParser类的区别在于，这个类允许计算一阶和二阶导数（以符号的方式），而FunctionParser类只计算一阶导数，使用有限差分。对于复杂的表达式，这个类可能比FunctionParser类要慢。
   * @ingroup functions
   *
   */
  template <int dim, typename RangeNumberType = double>
  class SymbolicFunction : public Function<dim, RangeNumberType>
  {
  public:
    /**
     * 构造函数。
     * 产生的Function对象将有与符号表达式向量 @p function.
     * 中的条目一样多的组件。向量 @p function
     * 应包含一个涉及坐标符号参数 @p coordinate_symbols
     * 和可能的符号时间参数 @p time_symbol.
     * 的符号表达式的列表，可以用其他符号定义它，只要可选参数
     * @p user_substitution_map  ] 替换除  @p coordinate_symbols  和  @p
     * time_symbol.
     * 以外的所有符号，这在以下情况下很有用：例如，你想用你想用符号命名的材料参数来表达公式，而不是在定义公式时通过它们的数值，或者你想用极坐标而不是笛卡尔坐标来表达你的公式，并且你希望符号引擎能为你计算导数。
     * 你以后可以通过调用update_user_substitution_map()来更新 @p
     * user_substitution_map 中包含的符号图。          @param  函数
     * 类型为 Differentiation::SD::Expression,
     * 的符号表达式的向量，代表此Function的组成部分。
     * @param  coordinate_symbols 代表坐标的符号张量，作为 @p
     * function 矢量中包含的符号表达式的输入参数。默认的
     * @p coordinate_symbols 是一个
     * Tensor<1,dim,Differentiation::SD::Expression> ，包含符号 "x
     * "代表`dim`等于1，"x"、"y "代表`dim`等于2，而 "x"、"y"、"z
     * "代表`dim`等于3。          @param  time_symbol
     * 一个代表时间的符号变量。它默认为一个名为 "t
     * "的符号变量。          @param  user_substitution_map
     * 任何可能包含在符号函数中的其他符号都需要在这个映射中指定。该地图可以是空的，函数仍然可以包含未评估的符号，只要你调用update_user_substitution_map()，并在任何评估发生之前提供除
     * @p coordinate_symbols 和 @p time_symbol 外的所有符号的替换。
     *
     */
    SymbolicFunction(
      const std::vector<Differentiation::SD::Expression> &function,
      const Tensor<1, dim, Differentiation::SD::Expression>
        &coordinate_symbols = get_default_coordinate_symbols(),
      const Differentiation::SD::Expression &time_symbol =
        Differentiation::SD::make_symbol("t"),
      const Differentiation::SD::types::substitution_map
        &user_substitution_map = {});

    /**
     * 构造函数，接收一个描述函数表达式的单一字符串，作为一个分号分隔的表达式列表。
     * 符号表达式可以使用默认的参数和默认的符号时间变量，再加上你可能需要的任何额外的符号，只要你在尝试评估函数或其导数之前，通过调用update_user_substitution_map()更新替换所有符号的用户替换图，并且使用set_additional_function_arguments()方法提供你函数的所有额外函数参数。
     *
     */
    SymbolicFunction(const std::string &expressions);

    /**
     * 存储并应用替换映射 @p substitutions
     * 到这个Function对象的每个符号组件。
     * 注意，这个方法将触发对每个分量的梯度、Hessians和Laplacians的重新计算。
     *
     */
    void
    update_user_substitution_map(
      const Differentiation::SD::types::substitution_map &substitutions);

    /**
     * 设置额外的 @p arguments ，在下一个评估步骤中被替换。
     * 注意， @p arguments 是在*评估 @p
     * permanent_user_substitution_map,
     * 和计算完所有导数之后被替换的。如果你传递的额外参数仍然依赖于坐标或时间符号，那么导数的评估将导致部分导数评估。
     * 这种方法提供了一种方法来评估依赖于更多参数而不仅仅是坐标和时间的函数。如果你想对复杂的符号表达式计算总导数，你应该调用update_user_substitution_map()代替。
     *
     */
    void
    set_additional_function_arguments(
      const Differentiation::SD::types::substitution_map &arguments);

    /**
     * 返回一个坐标符号的张量，可以用来定义这个符号函数对象的表达式。
     * 默认参数是一个 Tensor<1,dim,Differentiation::SD::Expression>
     * ，包含符号 "x "代表`dim`等于1，"x"，"y
     * "代表`dim`等于2，"x"，"y"，"z "代表`dim`等于3。
     *
     */
    static Tensor<1, dim, Differentiation::SD::Expression>
    get_default_coordinate_symbols();

    /**
     * 获取符号函数中用于坐标的实际参数。这个对象不包括任何用户定义的参数。
     *
     */
    const Tensor<1, dim, Differentiation::SD::Expression> &
    get_coordinate_symbols() const;

    /**
     * 获取该符号函数中实际使用的符号时间。
     *
     */
    const Differentiation::SD::Expression &
    get_time_symbol() const;

    /**
     * 获取该符号化函数中实际使用的符号化表达式。
     *
     */
    const std::vector<Differentiation::SD::Expression> &
    get_symbolic_function_expressions() const;

    /**
     * 获取当前存储的 @p user_substitution_map. 。
     *
     */
    const Differentiation::SD::types::substitution_map &
    get_user_substitution_map() const;

    /**
     * 返回一个SymbolicFunction对象，表示这个函数的时间导数。空间参数、符号时间和当前存储的用户替换图被转发给新函数。
     *
     */
    SymbolicFunction<dim, RangeNumberType>
    time_derivative() const;

    // documentation inherited from the base class
    virtual RangeNumberType
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    // documentation inherited from the base class
    virtual Tensor<1, dim, RangeNumberType>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    // documentation inherited from the base class
    virtual RangeNumberType
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    // documentation inherited from the base class
    virtual SymmetricTensor<2, dim, RangeNumberType>
    hessian(const Point<dim> & p,
            const unsigned int component = 0) const override;

    /**
     * 打印存储的参数和函数表达式，因为它将在调用方法value()时被评估。
     *
     */
    template <typename StreamType>
    StreamType &
    print(StreamType &out) const;

  private:
    /**
     * 返回一个替换图，用 @p point,
     * 的值替换参数，用this->get_time()的值替换符号时间，用
     * @p additional_function_arguments.
     * 给出的替换图替换任何其他参数。
     *
     */
    Differentiation::SD::types::substitution_map
    create_evaluation_substitution_map(const Point<dim> &point) const;

    /**
     * 重新计算函数的符号值，应用用户的置换图。这可能是一个昂贵的计算，只有在必要时才会调用。
     *
     */
    void
    update_values() const;

    /**
     * 重新计算函数的符号梯度，应用用户置换图。这可能是一个昂贵的计算，只有在必要时才会调用。
     *
     */
    void
    update_first_derivatives() const;

    /**
     * 重新计算函数的符号Hessian和符号Lapalacian。这可能是一个昂贵的计算，只有在必要时才会调用。
     *
     */
    void
    update_second_derivatives() const;

    /**
     * 这个符号函数的组件，在发生任何子结构之前。这是不可改变的，并且在构造时生成。
     * 在任何评估发生之前， @p user_substitution_map
     * 被应用于这个对象，其结果被存储在内部变量函数中。
     * 在评估过程中， @p symbolic_coordinate, 、 @p symbolic_time,
     * 和任何剩余的符号被替换为输入评估点、当前时间和
     * @p additional_function_arguments. 的内容。
     *
     */
    const std::vector<Differentiation::SD::Expression> user_function;

    /**
     * 存储用于表达式置换的用户置换图。这可以通过调用update_user_substitution_map()来更新。请注意，该函数仍然可以有未解决的符号，只要通过调用set_additional_function_arguments()来解决这些符号。
     *
     */
    Differentiation::SD::types::substitution_map user_substitution_map;

    /**
     * 存储一个用于额外参数替换的用户替换图。这将通过调用set_additional_function_arguments()来更新。
     *
     */
    Differentiation::SD::types::substitution_map additional_function_arguments;

    /**
     * 这个符号函数的实际成分。这是从应用 @p
     * user_substitution_map. 后的 @p user_function, 中得到的。
     *
     */
    mutable std::vector<Differentiation::SD::Expression> function;

    /**
     * 这个符号函数的每个分量的梯度。这是通过计算对象
     * @p function, 的符号梯度得到的，即应用 @p
     * user_substitution_map 到 @p user_function. 之后的梯度。
     *
     */
    mutable std::vector<Tensor<1, dim, Differentiation::SD::Expression>>
      function_gradient;

    /**
     * 这个符号函数的每个分量的Hessians。这是通过计算对象
     * @p function, 的符号Hessian得到的，也就是说，在应用 @p
     * user_substitution_map 到 @p user_function. 之后
     *
     */
    mutable std::vector<Tensor<2, dim, Differentiation::SD::Expression>>
      function_hessian;

    /**
     * 这个符号函数的每个分量的拉普拉斯系数。这是通过计算对象
     * @p function, 的符号拉普拉斯，即应用 @p user_substitution_map
     * 到 @p user_function. 之后得到的。
     *
     */
    mutable std::vector<Differentiation::SD::Expression> function_laplacian;

    /**
     * 函数的坐标符号参数。
     *
     */
    Tensor<1, dim, Differentiation::SD::Expression> coordinate_symbols;

    /**
     * 该函数的符号时间参数。
     *
     */
    mutable Differentiation::SD::Expression time_symbol;
  };

  /**
   * 允许使用位数左移运算符输出。
   *
   */
  template <int dim, typename RangeNumberType>
  inline std::ostream &
  operator<<(std::ostream &out, const SymbolicFunction<dim, RangeNumberType> &f)
  {
    return f.print(out);
  }



  // Inline and template functions
  template <int dim, typename RangeNumberType>
  template <typename StreamType>
  StreamType &
  SymbolicFunction<dim, RangeNumberType>::print(StreamType &out) const
  {
    for (unsigned int i = 0; i < dim; ++i)
      out << coordinate_symbols[i] << ", ";
    for (const auto &argument_pair : additional_function_arguments)
      out << argument_pair.first << ", ";
    out << time_symbol << " -> " << user_function[0];
    for (unsigned int i = 1; i < user_function.size(); ++i)
      out << "; " << user_function[i];
    if (!user_substitution_map.empty())
      {
        out << " # ( ";
        std::string sep = "";
        for (const auto &substitution : user_substitution_map)
          {
            out << sep << substitution.first << " = " << substitution.second;
            sep = ", ";
          }
        out << " )";
      }
    return out;
  }
#else
  template <int dim, typename RangeNumberType = double>
  class SymbolicFunction
  {
  public:
    SymbolicFunction()
    {
      AssertThrow(
        false,
        ExcMessage(
          "This class is not available if you did not enable SymEngine "
          "when compiling deal.II."));
    }
  };
#endif
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif


