//include/deal.II-translator/base/mutable_bind_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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

#ifndef dealii_base_mutable_bind_h
#define dealii_base_mutable_bind_h

#include <deal.II/base/config.h>

#include <deal.II/base/patterns.h>
#include <deal.II/base/std_cxx17/tuple.h>

#include <tuple>
#include <utility>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  /**
   * std::bind,
   * 的一个可变版本，它将一个函数指针的所有参数绑定到一个存储的元组上，并允许你在调用之间更新元组。
   * 这个类的一个使用实例是通过辅助函数mutable_bind()，该函数根据其参数在运行中创建一个MutableBind对象。
   * @code
   * void my_function(const int &a, const double &b);
   *
   * auto bound = mutable_bind(my_function, 1, 2.0);
   *
   * bound(); // will execute my_function(1, 2.0);
   *
   * bound.set_arguments(2, 3.0);
   * bound(); // will execute my_function(2, 3.0);
   *
   * bound.parse_arguments("3: 4.0");
   * bound(); // will execute my_function(3, 4.0);
   * @endcode
   * 参数被复制到元组中，并去除它们的引用和const属性。只有可复制的构造对象才被允许作为函数参数。如果你需要保留一些引用，你可以把你的函数包装成一个lambda函数。
   * @code
   * void
   * example_function(const Point<2> &p,
   *                 const double &d,
   *                 const unsigned int i = 3) {
   * ...
   * };
   *
   * const Point<2> p(1, 2);
   *
   * Utilities::MutableBind<void, double, unsigned int> exp = {
   *  [&p](const double &d,
   *       const unsigned int i)
   *  {
   *    example_function(p, d, i);
   *  },
   *  {}};
   *
   * exp.parse_arguments("3.0 : 4");
   * exp(); // calls example_function(p, 3.0, 4);
   * @endcode
   *
   *
   */
  template <typename ReturnType, class... FunctionArgs>
  class MutableBind
  {
  public:
    /**
     * 对存储的 std::tuple
     * 类型的别名。只有可复制构造的对象才允许作为元组成员。
     *
     */
    using TupleType = std::tuple<typename std::remove_cv<
      typename std::remove_reference<FunctionArgs>::type>::type...>;

    /**
     * 构建一个MutableBind对象，指定函数，以及每个参数分别。
     *
     */
    template <class FunctionType>
    MutableBind(FunctionType function, FunctionArgs &&... arguments);

    /**
     * 构造一个MutableBind对象，指定函数，以及作为一个元组的参数。
     *
     */
    template <class FunctionType>
    MutableBind(FunctionType function, TupleType &&arguments);

    /**
     * 构造一个只指定函数的MutableBind对象。默认情况下，参数被保留为其默认的构造函数值。
     *
     */
    template <class FunctionType>
    MutableBind(FunctionType function);

    /**
     * 调用原始函数，将绑定参数的元组作为参数传递。
     *
     */
    ReturnType
    operator()() const;

    /**
     * 设置参数，以便在下次调用operator()()时使用移动语义，在
     * @p function, 中使用。
     *
     */
    void
    set_arguments(TupleType &&arguments);

    /**
     * 设置在 @p function,
     * 中使用的参数，以便在下次调用operator()()时，使用移动语义。
     *
     */
    void
    set_arguments(FunctionArgs &&... arguments);

    /**
     * 将 @p function
     * 中的参数从字符串中解析出来，以便下次调用operator()()时使用。
     * 该转换是使用用户提供的 Patterns::PatternBase
     * 对象进行的。默认情况下，
     * Patterns::Tools::Convert<TupleType>::to_pattern()
     * 被用来确定如何从 @p value_string 转换为TupleType对象。
     * @param  value_string 要转换的字符串  @param  pattern
     * 执行转换时使用的模式的唯一指针
     *
     */
    void
    parse_arguments(const std::string &          value_string,
                    const Patterns::PatternBase &pattern =
                      *Patterns::Tools::Convert<TupleType>::to_pattern());

  private:
    /**
     * 一个存储原始函数的 std::function 。
     *
     */
    const std::function<ReturnType(FunctionArgs...)> function;

    /**
     * 目前存储的参数。当调用operator()()时，这些参数会被转发给上面的函数对象。
     *
     */
    TupleType arguments;
  };



  /**
   * 从一个函数指针和一个参数列表创建一个MutableBind对象。
   * 一个用法的例子是这样给出的。
   * @code
   * void my_function(const int &a, const double &b);
   *
   * auto bound = mutable_bind(my_function, 1, 2.0);
   *
   * bound(); // will execute my_function(1, 2.0);
   *
   * bound.set_arguments(2, 3.0);
   * bound(); // will execute my_function(2, 3.0);
   *
   * bound.parse_arguments("3: 4.0");
   * bound(); // will execute my_function(3, 4.0);
   * @endcode
   *
   *
   */
  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
  mutable_bind(ReturnType (*function)(FunctionArgs...),
               typename identity<FunctionArgs>::type &&... arguments);

  /**
   * 和上面一样，使用 std::function 对象。
   *
   */
  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
  mutable_bind(std::function<ReturnType(FunctionArgs...)>,
               typename identity<FunctionArgs>::type &&... arguments);

  /**
   * 从一个函数指针创建一个MutableBind对象，带有未初始化的参数。
   * 注意，如果你不调用 MutableBind::set_arguments()
   * 方法中的一个，或者在返回的对象上调用
   * MutableBind::parse_arguments()
   * 函数，那么传递给函数对象的参数将被初始化，其值来自每个参数的默认构造函数。
   *
   */
  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
    mutable_bind(ReturnType (*function)(FunctionArgs...));

  /**
   * 与上述相同，使用 std::function 对象。
   *
   */
  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
    mutable_bind(std::function<ReturnType(FunctionArgs...)>);



#ifndef DOXYGEN
  template <typename ReturnType, class... FunctionArgs>
  template <class FunctionType>
  MutableBind<ReturnType, FunctionArgs...>::MutableBind(
    FunctionType function,
    FunctionArgs &&... arguments)
    : function(function)
    , arguments(std::make_tuple(std::move(arguments)...))
  {}



  template <typename ReturnType, class... FunctionArgs>
  template <class FunctionType>
  MutableBind<ReturnType, FunctionArgs...>::MutableBind(FunctionType function,
                                                        TupleType && arguments)
    : function(function)
    , arguments(std::move(arguments))
  {}



  template <typename ReturnType, class... FunctionArgs>
  template <class FunctionType>
  MutableBind<ReturnType, FunctionArgs...>::MutableBind(FunctionType function)
    : function(function)
  {}



  template <typename ReturnType, class... FunctionArgs>
  ReturnType
  MutableBind<ReturnType, FunctionArgs...>::operator()() const
  {
    return std_cxx17::apply(function, arguments);
  }



  template <typename ReturnType, class... FunctionArgs>
  void
  MutableBind<ReturnType, FunctionArgs...>::set_arguments(
    FunctionArgs &&... args)
  {
    arguments = std::make_tuple(std::move(args)...);
  }



  template <typename ReturnType, class... FunctionArgs>
  void
  MutableBind<ReturnType, FunctionArgs...>::set_arguments(TupleType &&args)
  {
    arguments = std::move(args);
  }



  template <typename ReturnType, class... FunctionArgs>
  void
  MutableBind<ReturnType, FunctionArgs...>::parse_arguments(
    const std::string &          value_string,
    const Patterns::PatternBase &pattern)
  {
    arguments =
      Patterns::Tools::Convert<TupleType>::to_value(value_string, pattern);
  }



  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
  mutable_bind(ReturnType (*function)(FunctionArgs...),
               typename identity<FunctionArgs>::type &&... arguments)
  {
    return MutableBind<ReturnType, FunctionArgs...>(function,
                                                    std::move(arguments)...);
  }



  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
    mutable_bind(ReturnType (*function)(FunctionArgs...))
  {
    return MutableBind<ReturnType, FunctionArgs...>(function);
  }



  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
  mutable_bind(std::function<ReturnType(FunctionArgs...)> function,
               typename identity<FunctionArgs>::type &&... arguments)
  {
    return MutableBind<ReturnType, FunctionArgs...>(function,
                                                    std::move(arguments)...);
  }



  template <typename ReturnType, class... FunctionArgs>
  MutableBind<ReturnType, FunctionArgs...>
  mutable_bind(std::function<ReturnType(FunctionArgs...)> function)
  {
    return MutableBind<ReturnType, FunctionArgs...>(function);
  }
#endif
} // namespace Utilities

DEAL_II_NAMESPACE_CLOSE

#endif


