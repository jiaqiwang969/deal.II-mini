//include/deal.II-translator/base/patterns_0.txt
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

#ifndef dealii_patterns_h
#define dealii_patterns_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/component_mask.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/archive/basic_archive.hpp>
#include <boost/core/demangle.hpp>
#include <boost/property_tree/ptree_fwd.hpp>
#include <boost/property_tree/ptree_serialization.hpp>
#include <boost/serialization/split_member.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <algorithm>
#include <array>
#include <deque>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// forward declarations for interfaces and friendship
#ifndef DOXYGEN
class LogStream;
class MultipleParameterLoop;
template <int dim>
class FunctionParser;
#endif

/**
 * 一些作为ParameterHandler类的模式的类的命名空间。这些类实现了一个接口，用于检查输入文件中的参数是否符合某种模式，如
 * "是布尔值"、"整数值 "等。
 *
 *
 * @ingroup input
 *
 *
 */
namespace Patterns
{
  /**
   * 用于声明通用接口的基类。这个类的目的主要是为了定义模式的接口，并强制派生类有一个<tt>clone</tt>函数。因此，用《设计模式》一书（Gamma
   * et al.）的语言来说，它是一个 "原型"。
   *
   */
  class PatternBase
  {
  public:
    /**
     * 使这个类和所有派生类的析构器成为虚拟的。
     *
     */
    virtual ~PatternBase() = default;

    /**
     * 如果给定的字符串与模式匹配，返回<tt>true</tt>。
     *
     */
    virtual bool
    match(const std::string &test_string) const = 0;

    /**
     * 可能的描述输出格式的列表。        大写字母的选择与
     * ParameterHandler::OutputStyle. 相似。
     *
     */
    enum OutputStyle
    {
      /**
       * 适合在所有内置继承类的静态公共成员函数中进行机器解析的简单文本。
       * 最好是人类可读的，但机器解析更为关键。
       *
       */
      Machine,
      /**
       * 易于人类阅读的纯文本格式，适合纯文本文档。
       *
       */
      Text,
      /**
       * 易于人类阅读的LaTeX格式，适合在手册中打印。
       *
       */
      LaTeX
    };

    /**
     * 返回一个描述该模式的字符串。
     *
     */
    virtual std::string
    description(const OutputStyle style = Machine) const = 0;

    /**
     * 返回一个指向该对象的精确拷贝的指针。这是必要的，因为我们想把这种类型的对象存储在容器中，我们需要复制对象而不知道它们的实际数据类型（我们只有指向基类的指针）。
     * 由这个函数返回的对象的所有权被传递给这个函数的调用者。
     *
     */
    virtual std::unique_ptr<PatternBase>
    clone() const = 0;

    /**
     * 确定这个对象的内存消耗（以字节为单位）的估计值。为了避免不必要的开销，我们不强迫派生类提供这个函数作为一个虚拟的重载函数，而是尝试把现在的对象投给已知的派生类之一，如果失败了，那么就取这个基类的大小来代替，并加上32字节（这个值是任意的，它应该考虑到虚拟函数表，和一些可能的数据元素）。由于周围通常没有几千个这种类型的对象，并且由于memory_consumption机制是用来找出许多兆字节范围内的内存的，这似乎是一个合理的近似值。
     * 另一方面，如果你知道你的类大大偏离了这个假设，你仍然可以重载这个函数。
     *
     */
    virtual std::size_t
    memory_consumption() const;
  };

  /**
   * 根据描述返回指向正确派生类的指针。
   *
   */
  std::unique_ptr<PatternBase>
  pattern_factory(const std::string &description);

  namespace internal
  {
    /**
     * 对指定的 @p style 的字符串 @p input
     * 进行转义，使字符按预期出现。例如，像_这样的字符在LateX中不能按原样书写，必须转义为\_。
     *
     */
    std::string
    escape(const std::string &input, const PatternBase::OutputStyle style);

  } // namespace internal

  /**
   * 测试字符串是否是一个整数。如果构造函数给出了边界，那么给出的整数也需要在这些边界所指定的区间内。请注意，与C++标准库中的惯例不同，这个区间的两个边界都是包容的；原因是在大多数情况下，我们需要封闭的区间，但对于非整数值来说，这些只能用包容的边界来实现。因此，我们总是使用封闭区间来保持一致性。
   * 如果给构造函数的上界小于下界，那么每个整数都是允许的。
   * 如果一个值只能是正数并且小于一个合理的上限（例如要执行的细化步骤的数量），或者在许多其他情况下，给出边界可能是有用的。
   *
   */
  class Integer : public PatternBase
  {
  public:
    /**
     * 最小的整数值。如果numeric_limits类可用，则使用此信息来获得极值，否则设置它，使该类理解为允许所有值。
     *
     */
    static const int min_int_value;

    /**
     * 最大的整数值。如果numeric_limits类可用，则使用此信息来获得极值，否则设置它，使该类理解为允许所有的值。
     *
     */
    static const int max_int_value;

    /**
     * 构造函数。可以指定一个有效参数必须在其范围内的界限。如果上界小于下界，那么就意味着整个整数集。默认值的选择是不对参数强制执行界限。
     * 请注意，当前类型的对象所隐含的范围包括两个界限值，即
     * @p upper_bound
     * 是一个允许的值，而不是像其他情况下经常做的那样表示一个半开的值。
     *
     */
    Integer(const int lower_bound = min_int_value,
            const int upper_bound = max_int_value);

    /**
     * 如果字符串是一个整数且其值在指定范围内，则返回<tt>true</tt>。
     *
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * 返回对有效字符串预期匹配的模式的描述。如果向构造函数指定了边界，那么就把它们包括在这个描述中。
     *
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * 返回一个当前对象的副本，该对象是在堆上新分配的。该对象的所有权被转移到该函数的调用者。
     *
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * 如果描述的开始部分与description_init匹配，则创建一个新的对象。
     * 该对象的所有权将被转移给该函数的调用者。
     *
     */
    static std::unique_ptr<Integer>
    create(const std::string &description);

  private:
    /**
     * 下限的值。满足该类的 @ref
     * match
     * 操作的数字必须等于这个值或者更大，如果区间的边界为有效范围。
     *
     */
    const int lower_bound;

    /**
     * 上限的值。满足本类的 @ref
     * match
     * 操作的数字必须等于这个值或更小，如果区间的边界为有效范围。
     *
     */
    const int upper_bound;

    /**
     * 描述的初始部分
     *
     */
    static const char *description_init;
  };

  /**
   * 测试字符串是否为<tt>double</tt>。如果给构造函数设定了边界，那么给定的整数也需要在这些边界所指定的区间内。请注意，与C++标准库中的惯例不同，这个区间的两个边界都是包容的；原因是在大多数情况下，我们需要封闭的区间，但对于非整数值来说，这些只能用包容的边界来实现。因此，我们总是使用封闭区间来保持一致性。
   * 如果给构造函数的上界小于下界，那么每个双精度数字都是允许的。
   * 如果一个值只能是正数并且小于一个合理的上限（例如阻尼参数经常只有在0和1之间才是合理的），或者在许多其他情况下，给出边界可能是有用的。
   *
   */
  class Double : public PatternBase
  {
  public:
    /**
     * 用作默认值的最小双值，取自 <tt>std::numeric_limits</tt>.
     * 。
     *
     */
    static const double min_double_value;

    /**
     * 用作默认值的最大双倍值，取自
     * <tt>std::numeric_limits</tt>. 。
     *
     */
    static const double max_double_value;

    /**
     * 构造函数。可以指定一个有效参数必须在的范围内。如果上界小于下界，那么就意味着整个双精度的数字集。默认值的选择使参数不受任何约束。
     *
     */
    Double(const double lower_bound = min_double_value,
           const double upper_bound = max_double_value);

    /**
     * 如果字符串是一个数字且其值在指定范围内，则返回<tt>true</tt>。
     *
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * 返回对有效字符串预期匹配的模式的描述。如果向构造函数指定了边界，那么就把它们包括在这个描述中。
     *
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * 返回一个当前对象的副本，该对象是在堆上新分配的。该对象的所有权被转移到该函数的调用者。
     *
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * 如果给定的 @p description
     * 是一个有效的格式（例如通过调用现有对象的description()创建），则使用
     * @p new 在堆上创建一个新的对象，否则使用 @p nullptr
     * 。返回的对象的所有权被转移给这个函数的调用者，应该使用
     * @p delete. 释放它。
     *
     */
    static std::unique_ptr<Double>
    create(const std::string &description);

  private:
    /**
     * 下限的值。满足该类的 @ref
     * match
     * 操作的数字必须等于这个值，或者更大，如果区间的边界形成一个有效的范围。
     *
     */
    const double lower_bound;

    /**
     * 上限的值。如果区间的边界形成一个有效的范围，满足该类的 @ref
     * match 操作的数字必须等于这个值或小于这个值。
     *
     */
    const double upper_bound;

    /**
     * 描述的初始部分
     *
     */
    static const char *description_init;
  };

  /**
   *
   */
  class Selection : public PatternBase
  {
  public:
    /**
     * 构造函数。以给定的参数作为有效字符串的规范。
     *
     */
    Selection(const std::string &seq);

    /**
     * 如果该字符串是传递给构造函数的描述列表中的一个元素，则返回<tt>true</tt>。
     *
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * 返回对有效字符串预期匹配的模式的描述。这里，这是传递给构造函数的有效字符串的列表。
     *
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * 返回当前对象的副本，它是在堆上新分配的。该对象的所有权被转移到该函数的调用者。
     *
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * 确定这个对象的内存消耗（以字节为单位）的估计值。
     *
     */
    std::size_t
    memory_consumption() const override;

    /**
     * 如果description的开头与description_init匹配，则创建一个新对象。
     * 该对象的所有权将转移给该函数的调用者。
     *
     */
    static std::unique_ptr<Selection>
    create(const std::string &description);

  private:
    /**
     * 传递给构造函数的有效字符串的列表。我们不使这个字符串成为常量，因为我们在构造函数中对它进行了一定的处理。
     *
     */
    std::string sequence;

    /**
     * 描述的初始部分
     *
     */
    static const char *description_init;
  };


  /**
   * 这个模式匹配一个由逗号（或另一个字符串）隔开的值的列表，每个值都必须匹配给构造器的模式。
   * 通过两个额外的参数，可以指定这个列表必须有的元素数量。如果没有指定，该列表可以有零个或更多的条目。
   *
   */
  class List : public PatternBase
  {
  public:
    /**
     * 最大的整数值。如果numeric_limits类可用，则使用此信息获得极值，否则设置它，使此类理解为允许所有值。
     *
     */
    static const unsigned int max_int_value;

    /**
     * 构造函数。以给定的参数作为列表的有效元素的规范。
     * 其他三个参数可以用来表示列表的最小和最大允许长度，以及作为列表元素之间分隔符的字符串。
     *
     */
    List(const PatternBase &base_pattern,
         const unsigned int min_elements = 0,
         const unsigned int max_elements = max_int_value,
         const std::string &separator    = ",");


    /**
     * 返回内部存储的分隔符。
     *
     */
    const std::string &
    get_separator() const;

    /**
     * 返回内部存储的基本模式。
     *
     */
    const PatternBase &
    get_base_pattern() const;

    /**
     * 复制构造函数。
     *
     */
    List(const List &other);

    /**
     * 如果字符串是一个以逗号分隔的字符串列表，其中每一个都与给构造函数的模式匹配，则返回<tt>true</tt>。
     *
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * 返回有效字符串预期匹配的模式的描述。
     *
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * 返回一个当前对象的副本，该对象是在堆上新分配的。该对象的所有权被转移给该函数的调用者。
     *
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * 如果描述的开始部分与description_init匹配，则创建一个新对象。
     * 该对象的所有权将被转移给该函数的调用者。
     *
     */
    static std::unique_ptr<List>
    create(const std::string &description);

    /**
     * 确定这个对象的内存消耗（以字节为单位）的估计值。
     *
     */
    std::size_t
    memory_consumption() const override;

    /**
     * @addtogroup  Exceptions 
     * @{ 
     *
     */

    /**
     * 异常情况。
     *
     */
    DeclException2(ExcInvalidRange,
                   int,
                   int,
                   << "The values " << arg1 << " and " << arg2
                   << " do not form a valid range.");
    //@}
  private:
    /**
     * 列表中每个元素必须满足的模式的拷贝。
     *
     */
    std::unique_ptr<PatternBase> pattern;

    /**
     * 列表必须有的最小元素数。
     *
     */
    const unsigned int min_elements;

    /**
     * 列表必须有的最大元素数。
     *
     */
    const unsigned int max_elements;

    /**
     * 列表中各元素之间的分隔符。
     *
     */
    const std::string separator;

    /**
     * 描述的初始部分
     *
     */
    static const char *description_init;
  };


  /**
   * 这个模式匹配一个用逗号分隔的列表，每个列表表示一对键和值。key和value都必须与给构造函数的模式匹配。对于地图的每个条目，必须以
   * <code>key: value</code>
   * 的形式输入参数。换句话说，一个地图的描述形式是<code>key1:
   * value1, key2: value2, key3: value3,
   * ...</code>。两个构造函数参数允许选择除逗号以外的对之间的分隔符，以及除冒号以外的键和值之间的分隔符。
   * 通过两个额外的参数，可以指定这个列表必须有的元素数量。如果没有指定，地图可以有零个或多个条目。
   *
   */
  class Map : public PatternBase
  {
  public:
    /**
     * 最大的整数值。如果numeric_limits类可用，则使用此信息获得极值，否则设置它，使此类理解为允许所有值。
     *
     */
    static const unsigned int max_int_value;

    /**
     * 构造函数。以给定的参数作为列表的有效元素的规范。
     * 其他四个参数可以用来表示列表的最小和最大的允许长度，以及用来划分map的对的分隔符和用来分隔key和value的符号。
     *
     */
    Map(const PatternBase &key_pattern,
        const PatternBase &value_pattern,
        const unsigned int min_elements        = 0,
        const unsigned int max_elements        = max_int_value,
        const std::string &separator           = ",",
        const std::string &key_value_separator = ":");

    /**
     * 复制构造函数。
     *
     */
    Map(const Map &other);

    /**
     * 如果字符串是一个逗号分隔的字符串列表，其中每个字符串都符合构造函数所给的模式，则返回<tt>true</tt>。
     *
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * 返回有效字符串预期匹配的模式的描述。
     *
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * 返回一个当前对象的副本，该对象是在堆上新分配的。该对象的所有权被转移给该函数的调用者。
     *
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * 如果描述的开始部分与description_init匹配，则创建一个新对象。
     * 该对象的所有权将被转移给该函数的调用者。
     *
     */
    static std::unique_ptr<Map>
    create(const std::string &description);

    /**
     * 确定这个对象的内存消耗（以字节为单位）的估计值。
     *
     */
    std::size_t
    memory_consumption() const override;

    /**
     * 返回一个对关键模式的引用。
     *
     */
    const PatternBase &
    get_key_pattern() const;

    /**
     * 返回一个对值模式的引用。
     *
     */
    const PatternBase &
    get_value_pattern() const;

    /**
     * 返回地图条目的分隔符。
     *
     */
    const std::string &
    get_separator() const;

    /**
     * 返回键值分离器。
     *
     */
    const std::string &
    get_key_value_separator() const;

    /**
     * @addtogroup  Exceptions  
     * @{ 
     *
     */

    /**
     * 异常情况。
     *
     */
    DeclException2(ExcInvalidRange,
                   int,
                   int,
                   << "The values " << arg1 << " and " << arg2
                   << " do not form a valid range.");
    //@}
  private:
    /**
     * 复制地图的每个键和每个值必须满足的模式。
     *
     */
    std::unique_ptr<PatternBase> key_pattern;
    std::unique_ptr<PatternBase> value_pattern;

    /**
     * 列表必须有的最小元素数。
     *
     */
    const unsigned int min_elements;

    /**
     * 列表必须有的最大元素数。
     *
     */
    const unsigned int max_elements;

    /**
     * 列表中各元素之间的分隔符。
     *
     */
    const std::string separator;


    /**
     * 键和值之间的分隔符。
     *
     */
    const std::string key_value_separator;

    /**
     * 描述的初始部分
     *
     */
    static const char *description_init;
  };



  /**
   * 这个模式匹配任意类型的以冒号分隔的值。每个类型都必须与构造函数所给的模式相匹配。
   * 下面是一个用法示例。
   * @code
   * std::vector< std::unique_ptr<Patterns::PatternBase> > ps;
   *
   * ps.push_back(std::unique_ptr<Patterns::Integer>());
   * ps.push_back(std::unique_ptr<Patterns::Double>());
   * ps.push_back(std::unique_ptr<Patterns::Anything>());
   *
   * Patterns::Tuple pattern(ps, ":");
   *
   * bool check = ps.match("5 : 3.14 : Ciao"); // check = true
   * @endcode
   * 或者，如果你想利用 ParameterHandler::add_parameter(): 。
   * @code
   * using T = std::tuple<std::string, Point<3>, unsigned int>;
   *
   * T a = Patterns::Tools::Convert<T>::to_value("Ciao : 1.0, 2.0, 3.0 : 33");
   *
   * ParameterHandler prm;
   * prm.add_parameter("A tuple", a);
   *
   * prm.log_parameters(deallog);
   * // DEAL:parameters::A tuple: Ciao : 1.000000, 2.000000, 3.000000 : 33
   *
   * prm.set("A tuple", "Mondo : 2.0, 3.0, 4.0 : 34");
   * prm.log_parameters(deallog);
   * // DEAL:parameters::A tuple: Mondo : 2.0, 3.0, 4.0 : 34
   *
   * deallog << Patterns::Tools::Convert<T>::to_string(a) << std::endl;
   * // DEAL::Mondo : 2.000000, 3.000000, 4.000000 : 34
   * @endcode
   * 构造函数希望得到一个模式的向量，以及一个指定从字符串解析Tuple时要使用的分隔符的字符串。
   * 默认的分隔符是冒号，这是因为一对元素实际上是一个有两个元素的元组。
   *
   */
  class Tuple : public PatternBase
  {
  public:
    /**
     * 构造函数。使用一个指向模式的唯一指针的向量来构造元组。
     * @param  patterns 元组的每个对象应该匹配的模式  @param
     * separator 用于限定每个元素的可选字符串 构建器。
     *
     */
    Tuple(const std::vector<std::unique_ptr<PatternBase>> &patterns,
          const std::string &                              separator = ":");

    /**
     * 构造函数。和上面一样，专门用于 const
     * char。这是必要的，以避免编译器因下面提供的变量构造器而产生错误。
     *
     */
    Tuple(const std::vector<std::unique_ptr<PatternBase>> &patterns,
          const char *                                     separator);


    /**
     * 构造函数。从一个以上的PatternBase派生的类中创建一个Tuple。
     * @param  分隔符 要使用什么分隔符。      @param  patterns
     * 要使用的模式列表。
     *
     */
    template <class... PatternTypes>
    Tuple(const std::string &separator, const PatternTypes &... patterns);

    /**
     * 构造器。这是需要的，以允许用户直接指定分离器而不使用
     * std::string(";").
     * 因为我们支持纯变体模板版本，如果没有这个特殊化，编译器将以隐性错误失败。
     *
     */
    template <class... PatternTypes>
    Tuple(const char *separator, const PatternTypes &... patterns);

    /**
     * 构造函数。和上面一样，使用默认的分离器。
     * @param  patterns 要使用的模式列表
     *
     */
    template <typename... Patterns>
    Tuple(const Patterns &... patterns);

    /**
     * 复制构造函数。
     *
     */
    Tuple(const Tuple &other);

    /**
     * 如果字符串是一个字符串的列表，其中每一个都与给构造函数的模式相匹配，则返回<tt>true</tt>。
     *
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * 返回对有效字符串预期匹配的模式的描述。
     *
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * 返回一个当前对象的副本，该对象是在堆上新分配的。该对象的所有权被转移给该函数的调用者。
     *
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * 如果描述的开始部分与description_init匹配，则创建一个新对象。
     * 该对象的所有权将被转移给该函数的调用者。
     *
     */
    static std::unique_ptr<Tuple>
    create(const std::string &description);

    /**
     * 确定这个对象的内存消耗（以字节为单位）的估计值。
     *
     */
    std::size_t
    memory_consumption() const override;

    /**
     * 返回一个对元组中第i个模式的引用。
     *
     */
    const PatternBase &
    get_pattern(const unsigned int i) const;

    /**
     * 返回元组条目的分隔符。
     *
     */
    const std::string &
    get_separator() const;

  private:
    /**
     * 复制存储在Tuple中的模式。
     *
     */
    std::vector<std::unique_ptr<PatternBase>> patterns;

    /**
     * 列表中元素之间的分隔符。
     *
     */
    const std::string separator;

    /**
     * 描述的初始部分。
     *
     */
    static const char *description_init;
  };


  /**
   * 这个类很像选择类，但它允许输入是一个逗号分隔的列表，每个值都必须在构造函数参数中给出。输入可以是空的，也可以包含多个值，而且逗号周围可以有任意数量的空格。当然，给构造函数的值中不允许有逗号。
   * 例如，如果给构造函数的字符串是<tt>"ucd|gmv|eps"</tt>，那么以下是合法输入。"eps",
   * "gmv, eps", 或者""。
   *
   */
  class MultipleSelection : public PatternBase
  {
  public:
    /**
     * 构建器。  @p seq 是一个由"|"分隔的有效选项的列表。
     *
     */
    MultipleSelection(const std::string &seq);

    /**
     * 如果字符串是传递给构造函数的描述列表中的一个元素，则返回<tt>true</tt>。
     *
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * 返回对有效字符串预期匹配的模式的描述。这里，这是传递给构造函数的有效字符串的列表。
     *
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * 返回当前对象的副本，它是在堆上新分配的。该对象的所有权被转移到该函数的调用者。
     *
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * 如果描述的开始部分与description_init匹配，则创建一个新对象。
     * 该对象的所有权将被转移给该函数的调用者。
     *
     */
    static std::unique_ptr<MultipleSelection>
    create(const std::string &description);

    /**
     * 确定这个对象的内存消耗（以字节为单位）的估计值。
     *
     */
    std::size_t
    memory_consumption() const override;

    /**
     * @addtogroup  Exceptions  @{
     *
     */

    /**
     * 异常情况。
     *
     */
    DeclException1(
      ExcCommasNotAllowed,
      int,
      << "A comma was found at position " << arg1
      << " of your input string, but commas are not allowed here.");
    //@}
  private:
    /**
     * 传递给构造函数的有效字符串的列表。我们不使这个字符串成为常量，因为我们在构造函数中会对它进行一些处理。
     *
     */
    std::string sequence;

    /**
     * 描述的初始部分
     *
     */
    static const char *description_init;
  };

  /**
   * 测试该字符串是 "真 "还是 "假"。这被映射到选择类中。
   *
   */
  class Bool : public Selection
  {
  public:
    /**
     * 构造函数。
     *
     */
    Bool();

    /**
     * 返回对有效字符串预期匹配的模式的描述。
     *
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * 返回一个当前对象的副本，该对象是在堆上新分配的。该对象的所有权被转移给该函数的调用者。
     *
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * 如果描述的开始部分与description_init匹配，则创建一个新的对象。
     * 该对象的所有权将被转移给该函数的调用者。
     *
     */
    static std::unique_ptr<Bool>
    create(const std::string &description);

  private:
    /**
     * 描述的初始部分
     *
     */
    static const char *description_init;
  };

  /**
   * 当测试一个字符串时总是返回<tt>true</tt>。
   *
   */
  class Anything : public PatternBase
  {
  public:
    /**
     * 构造函数。(允许在这个类中至少有一个非虚拟函数，因为否则有时没有虚拟表被发出)。
     *
     */
    Anything() = default;

    /**
     * 如果字符串符合其约束条件，则返回<tt>true</tt>，即总是如此。
     *
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * 返回对有效字符串预期匹配的模式的描述。这里，这是一个字符串<tt>"[Anything]"</tt>。
     *
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * 返回当前对象的副本，该对象是在堆上新分配的。该对象的所有权被转移到该函数的调用者。
     *
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * 如果描述的开始部分与description_init匹配，则创建一个新对象。
     * 该对象的所有权将被转移给该函数的调用者。
     *
     */
    static std::unique_ptr<Anything>
    create(const std::string &description);

  private:
    /**
     * 描述的初始部分
     *
     */
    static const char *description_init;
  };


  /**
   * 一个模式，可以用来指示一个参数何时打算成为一个文件的名称。就其本身而言，这个类并不检查在参数文件中给出的字符串是否真的对应于一个现有的文件（例如，它可能是一个你想写入输出的文件的名称）。因此，该类在功能上等同于Anything类。然而，它允许指定一个参数的<i>intent</i>。给予构造函数的标志也允许指定该文件应该是一个输入或输出文件。
   * 这个类存在的原因是为了支持编辑参数文件的图形用户界面。如果文件名应该代表一个输入文件，这些文件可以打开一个文件选择对话框。
   *
   */
  class FileName : public PatternBase
  {
  public:
    /**
     * 文件可以被用于输入或输出。这可以在构造函数中通过选择标志<tt>类型</tt>来指定。
     *
     */
    enum FileType
    {
      /**
       * 打开用于输入。
       *
       */
      input = 0,
      /**
       * 为输出而开放。
       *
       */
      output = 1
    };

    /**
     * 构造器。 文件的类型可以通过选择标志来指定。
     *
     */
    FileName(const FileType type = input);

    /**
     * 如果字符串符合其约束条件，则返回<tt>true</tt>，即总是如此。
     *
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * 返回对有效字符串预期匹配的模式的描述。这里，这是一个字符串<tt>"[Filename]"</tt>。
     *
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * 返回当前对象的副本，它是在堆上新分配的。该对象的所有权被转移给该函数的调用者。
     *
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * 文件类型标志
     *
     */
    FileType file_type;

    /**
     * 如果描述的开始部分与description_init匹配，则创建一个新对象。
     * 该对象的所有权将转移给该函数的调用者。
     *
     */
    static std::unique_ptr<FileName>
    create(const std::string &description);

  private:
    /**
     * 描述的初始部分
     *
     */
    static const char *description_init;
  };


  /**
   * 一个模式，可用于指示参数何时打算成为一个目录的名称。就其本身而言，这个类并不检查在参数文件中给出的字符串是否真的对应于一个现有的目录。因此在功能上，该类等同于Anything类。然而，它允许指定一个参数的<i>intent</i>。
   * 这个类存在的原因是为了支持编辑参数文件的图形用户界面。这些可以打开一个文件选择对话框来选择或创建一个目录。
   *
   */
  class DirectoryName : public PatternBase
  {
  public:
    /**
     * 构造函数。
     *
     */
    DirectoryName() = default;

    /**
     * 如果字符串符合其约束条件，则返回<tt>true</tt>，即总是如此。
     *
     */
    virtual bool
    match(const std::string &test_string) const override;

    /**
     * 返回对有效字符串预期匹配的模式的描述。这里，这是一个字符串<tt>"[Filename]"</tt>。
     *
     */
    virtual std::string
    description(const OutputStyle style = Machine) const override;

    /**
     * 返回当前对象的副本，它是在堆上新分配的。该对象的所有权被转移给该函数的调用者。
     *
     */
    virtual std::unique_ptr<PatternBase>
    clone() const override;

    /**
     * 如果描述的开始部分与description_init匹配，则创建一个新对象。
     * 该对象的所有权将被转移给该函数的调用者。
     *
     */
    static std::unique_ptr<DirectoryName>
    create(const std::string &description);

  private:
    /**
     * 描述的初始部分
     *
     */
    static const char *description_init;
  };


  /**
   * 一些类和函数的命名空间，它们作用于值和模式，并允许从非基本类型转换为字符串，反之亦然。
   * 这些工具的一个典型用法是在下面的例子中。
   * @code
   * using T = std::vector<unsigned int>;
   *
   * T vec(3);
   * vec[0] = 1;
   * vec[1] = 3;
   * vec[2] = 5;
   *
   * auto pattern = Patterns::Tools::Convert<T>::to_pattern();
   *
   * std::cout << pattern->description() << std::endl;
   * // [List of <[Integer]> of length 0...4294967295 (inclusive)]
   *
   * auto s = Patterns::Tools::Convert<T>::to_string(vec);
   * std::cout << s << std::endl;
   * // 1, 2, 3
   *
   * auto vec = Patterns::Tools::Convert<T>::to_value("2,3,4,5");
   * // now vec has size 4, and contains the elements 2,3,4,5
   *
   * std::cout << internal::RankInfo<T>::list_rank << std::endl; // Outputs 1
   * std::cout << internal::RankInfo<T>::map_rank  << std::endl; // Outputs 0
   * @endcode
   * Convert<T>被这个命名空间中的函数
   * Patterns::Tools::add_parameter() 使用。在内部，它使用
   * internal::RankInfo<T>
   * 类来决定需要多少个不同的分隔符来将给定的类型转换为字符串。
   * 例如，要写向量的向量，默认是使用",
   * "作为第一个（内部）分隔符，而";
   * "作为第二个（外部）分隔符，即
   * @code
   * std::vector<std::vector<unsigned int>> vec;
   * vec = Convert<decltype(vec)>::to_value("1,2,3 ; 4,5,6");
   *
   * s = convert<decltype(vec[0])>::to_string(vec[0]);
   * // s now contains the string "1,2,3"
   * @endcode
   * Patterns::List 和 Patterns::Map
   * 兼容类型的分隔符根据列表和地图对象的等级来选择，使用数组
   * Patterns::Tools::internal::default_list_separator 和
   * Patterns::Tools::internal::default_map_separator.
   * 它们目前被设置为。
   * @code
   * default_list_separator{{","  ,  ";"  ,  "|"  ,   "%"}};
   * default_map_separator {{":"  ,  "="  ,  "@"  ,   "#"}};
   * @endcode
   * 当人们需要 Patterns::List 和 Patterns::Map
   * 类型的混合物时，它们的RankInfo是通过取Key和Value类型的vector_rank的最大值来计算的，因此，例如，可以有以下情况
   * @code
   * ... // Build compare class
   * std::map<std::vector<unsigned int>, std::vector<double>, compare> map;
   *
   * map = convert<decltype(map)>::to_value(
   * "1,2,3 : 5.0,6.0,7.0  ; 8,9,10 : 11.0,12.0,13.0");
   *
   * @endcode
   * 支持一些非基本类型，比如Point()，或者
   * std::complex<double>.
   * 如果你希望支持更多的类型，你必须对Convert结构以及RankInfo结构进行专业化。
   * @ingroup input
   *
   */
  namespace Tools
  {
    /**
     * 转换器类。这个类用于生成与给定类型相关的字符串和模式，以及从字符串转换到给定类型，反之亦然。
     * 第二个模板参数在内部使用，以允许先进的SFINAE（替换失败不是错误）技巧，用于为任意的STL容器和映射专门化这个类。
     *
     */
    template <class T, class Enable = void>
    struct Convert
    {
      /**
       * 返回一个 std::unique_ptr
       * 给Pattern，可以用来解释一个字符串为模板参数的类型，反之亦然。
       * 虽然当前的函数（在一般的Convert模板中）被删除了，但它在Convert类模板的特殊化中被实现并可用于特定种类的模板参数
       * @p T.  。
       *
       */
      static std::unique_ptr<Patterns::PatternBase>
      to_pattern() = delete;

      /**
       * 返回一个包含变量s的文本版本的字符串。使用传递的模式来执行转换，或者创建并使用一个默认的模式。
       * 虽然当前的函数（在一般的Convert模板中）被删除了，但它在Convert类模板的特殊化中被实现并可用于特定种类的模板参数
       * @p T.  。
       *
       */
      static std::string
      to_string(const T &                    s,
                const Patterns::PatternBase &p = *Convert<T>::to_pattern()) =
        delete;

      /**
       * 使用给定的模式，将一个字符串转换为一个值。使用传递的模式来执行转换，或者创建并使用一个默认的模式。
       * 虽然当前的函数（在一般的Convert模板中）被删除了，但它在Convert类模板的特殊类型的模板参数中被实现并可用
       * @p T.  。
       *
       */
      static T
      to_value(const std::string &          s,
               const Patterns::PatternBase &p = *Convert<T>::to_pattern()) =
        delete;
    };

    /**
     * 一个实用的函数，简化了对任意复杂类型的字符串的转换。
     * 这个函数以默认模式调用方法 Convert<T>::to_string()
     * 。下面是一个用法示例。
     * @code
     * auto t = std::make_tuple(1.0, std::make_pair(1, "ciao"));
     * auto s = Patterns::Tools::to_string(t);
     *
     * std::cout << s; // will print "1 % 1 : ciao""
     * @endcode
     * 参见类 Patterns::Tools::Convert, 和辅助类
     * Patterns::Tools::RankInfo
     * 的文档，以了解输出STL容器类型时选择分隔符的方式的细节。
     *
     */
    template <typename T>
    std::string
    to_string(const T &t);

    /**
     * 一个实用的函数，简化了从字符串到任意类型的转换。
     * 这个函数以默认模式调用方法 Convert<T>::to_value()
     * 。下面是一个用法示例。
     * @code
     * auto t = std::make_tuple(1.0, std::make_pair(1, "ciao"));
     * // replace the value of 't' by the parsed content of the string argument:
     * Patterns::Tools::to_value("2 % 3 : mondo", t);
     *
     * auto s = Patterns::Tools::to_string(t);
     * std::cout << s; // will print "2 % 3 : mondo""
     * @endcode
     * 参见类 Patterns::Tools::Convert, 和辅助类
     * Patterns::Tools::RankInfo
     * 的文档，以了解从字符串转换到容器类型时应该在字符串模式中使用的分隔符。
     * 注意，变量 @p t
     * 的当前内容被忽略了。它的类型是用来推断如何解释字符串的。如果字符串被成功解析，那么
     * @p t 将被设置为 @p s. 的解析内容。
     *
     */
    template <typename T>
    void
    to_value(const std::string &s, T &t);

    /**
     * @addtogroup  Exceptions  
     * @{ 
     *
     */

    /**
     * 异常情况。
     *
     */
    DeclException2(ExcNoMatch,
                   std::string,
                   std::string,
                   << "The string " << arg1 << " does not match the pattern \""
                   << arg2 << "\"");
    //@}
  } // namespace Tools
} // namespace Patterns


// ---------------------- inline and template functions --------------------
namespace Patterns
{
  template <class... PatternTypes>
  Tuple::Tuple(const char *separator, const PatternTypes &... ps)
    : // forward to the version with std::string argument
    Tuple(std::string(separator), ps...)
  {}



  template <class... PatternTypes>
  Tuple::Tuple(const std::string &separator, const PatternTypes &... ps)
    : separator(separator)
  {
    static_assert(is_base_of_all<PatternBase, PatternTypes...>::value,
                  "Not all of the input arguments of this function "
                  "are derived from PatternBase");
    static_assert(sizeof...(ps) > 0,
                  "The number of PatternTypes must be greater than zero!");
    const auto pattern_pointers = {(static_cast<const PatternBase *>(&ps))...};
    for (const auto p : pattern_pointers)
      patterns.push_back(p->clone());
  }



  template <class... PatternTypes>
  Tuple::Tuple(const PatternTypes &... ps)
    : // forward to the version with the separator argument
    Tuple(std::string(":"), ps...)
  {}



  namespace Tools
  {
    namespace internal
    {
      /**
       * 存储关于给定类的等级类型的信息。
       * 一个类的等级等于在一个字符串中唯一识别其元素所需的不同分隔符的数量。
       * 这个类用于检测类T是否与 Patterns::List 模式或与
       * Patterns::Map 模式兼容。            像Point()或
       * std::complex<double>
       * 这样的对象是向量类，并且有vector_rank
       * 1。基本类型，如 "int"、"unsigned int"、"double
       * "等，具有vector_rank 0。  `std::vector`,   `std::list`
       * 一般来说，容器的等级等于1+所含类型的vector_rank。对于地图类型也是如此。
       * 一个 list_rank::value
       * =0的类要么是基本类要么是地图。一个 map_rank::value
       * =0的类要么是一个List兼容类，要么是一个基本类型。
       * 基本类型与 Patterns::List,
       * 不兼容，但非基本类型，如Point()，或
       * std::complex<double>,
       * 则与List类型兼容。添加更多的兼容类型是为给定的类型添加此结构的特殊化问题。
       *
       */
      template <class T, class Enable = void>
      struct RankInfo
      {
        static constexpr int list_rank = 0;
        static constexpr int map_rank  = 0;
      };
    } // namespace internal

    // Arithmetic types
    template <class T>
    struct Convert<T,
                   typename std::enable_if<std::is_arithmetic<T>::value>::type>
    {
      template <typename Dummy = T>
      static
        typename std::enable_if<std::is_same<Dummy, T>::value &&
                                  std::is_same<T, bool>::value,
                                std::unique_ptr<Patterns::PatternBase>>::type
        to_pattern()
      {
        return std::make_unique<Patterns::Bool>();
      }

      template <typename Dummy = T>
      static
        typename std::enable_if<std::is_same<Dummy, T>::value &&
                                  !std::is_same<T, bool>::value &&
                                  std::is_integral<T>::value,
                                std::unique_ptr<Patterns::PatternBase>>::type
        to_pattern()
      {
        return std::make_unique<Patterns::Integer>(
          std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max());
      }

      template <typename Dummy = T>
      static
        typename std::enable_if<std::is_same<Dummy, T>::value &&
                                  !std::is_same<T, bool>::value &&
                                  std::is_floating_point<T>::value,
                                std::unique_ptr<Patterns::PatternBase>>::type
        to_pattern()
      {
        return std::make_unique<Patterns::Double>(
          std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max());
      }

      static std::string
      to_string(const T &                    value,
                const Patterns::PatternBase &p = *Convert<T>::to_pattern())
      {
        std::stringstream str;
        if (std::is_same<T, unsigned char>::value ||
            std::is_same<T, signed char>::value || std::is_same<T, char>::value)
          str << static_cast<int>(value);
        else if (std::is_same<T, bool>::value)
          str << (static_cast<bool>(value) ? "true" : "false");
        else
          str << value;
        AssertThrow(p.match(str.str()), ExcNoMatch(str.str(), p.description()));
        return str.str();
      }

      static T
      to_value(const std::string &          s,
               const Patterns::PatternBase &p = *Convert<T>::to_pattern())
      {
        AssertThrow(p.match(s), ExcNoMatch(s, p.description()));
        T value;
        if (std::is_same<T, bool>::value)
          value = (s == "true");
        else
          {
            std::istringstream is(s);
            if (std::is_same<T, unsigned char>::value ||
                std::is_same<T, signed char>::value ||
                std::is_same<T, char>::value)
              {
                int i;
                is >> i;
                value = i;
              }
            else
              is >> value;

            // If someone passes "123 abc" to the function, the method yields an
            // integer 123 alright, but the space terminates the read from the
            // string although there is more to come. This case, however, is
            // checked for in the call p->match(s) at the beginning of this
            // function, and would throw earlier. Here it is safe to assume that
            // if we didn't fail the conversion with the operator >>, then we
            // are good to go.
            AssertThrow(
              !is.fail(),
              ExcMessage("Failed to convert from \"" + s + "\" to the type \"" +
                         boost::core::demangle(typeid(T).name()) + "\""));
          }
        return value;
      }
    };

    namespace internal
    {
      constexpr std::array<const char *, 4> default_list_separator{
        {",", ";", "|", "%"}};
      constexpr std::array<const char *, 4> default_map_separator{
        {":", "=", "@", "#"}};

      // specialize a type for all of the STL containers and maps
      template <typename T>
      struct is_list_compatible : std::false_type
      {};
      template <typename T, std::size_t N>
      struct is_list_compatible<std::array<T, N>> : std::true_type
      {};
      template <typename... Args>
      struct is_list_compatible<std::vector<Args...>> : std::true_type
      {};
      template <typename... Args>
      struct is_list_compatible<std::deque<Args...>> : std::true_type
      {};
      template <typename... Args>
      struct is_list_compatible<std::list<Args...>> : std::true_type
      {};
      template <typename... Args>
      struct is_list_compatible<std::set<Args...>> : std::true_type
      {};
      template <typename... Args>
      struct is_list_compatible<std::multiset<Args...>> : std::true_type
      {};
      template <typename... Args>
      struct is_list_compatible<std::unordered_set<Args...>> : std::true_type
      {};
      template <typename... Args>
      struct is_list_compatible<std::unordered_multiset<Args...>>
        : std::true_type
      {};

      template <typename T>
      struct is_map_compatible : std::false_type
      {};
      template <class Key, class T, class Compare, class Allocator>
      struct is_map_compatible<std::map<Key, T, Compare, Allocator>>
        : std::true_type
      {};
      template <class Key, class T, class Compare, class Allocator>
      struct is_map_compatible<std::multimap<Key, T, Compare, Allocator>>
        : std::true_type
      {};
      template <class Key, class T, class Hash, class KeyEqual, class Allocator>
      struct is_map_compatible<
        std::unordered_map<Key, T, Hash, KeyEqual, Allocator>> : std::true_type
      {};
      template <class Key, class T, class Hash, class KeyEqual, class Allocator>
      struct is_map_compatible<
        std::unordered_multimap<Key, T, Hash, KeyEqual, Allocator>>
        : std::true_type
      {};
    } // namespace internal

    // type trait to use the implementation type traits as well as decay the
    // type
    template <typename T>
    struct is_list_compatible
    {
      static constexpr bool const value =
        internal::is_list_compatible<typename std::decay<T>::type>::value;
    };

    template <typename T>
    struct is_map_compatible
    {
      static constexpr bool const value =
        internal::is_map_compatible<typename std::decay<T>::type>::value;
    };

    namespace internal
    {
      // Helper function for list_rank
      template <class T>
      constexpr int
      max_list_rank()
      {
        return RankInfo<T>::list_rank;
      }

      template <class T1, class T2, class... Types>
      constexpr int
      max_list_rank()
      {
        return std::max(RankInfo<T1>::list_rank, max_list_rank<T2, Types...>());
      }

      // Helper function for map_rank
      template <class T>
      constexpr int
      max_map_rank()
      {
        return RankInfo<T>::map_rank;
      }

      template <class T1, class T2, class... Types>
      constexpr int
      max_map_rank()
      {
        return std::max(RankInfo<T1>::map_rank, max_map_rank<T2, Types...>());
      }

      // Rank of vector types
      template <class T>
      struct RankInfo<
        T,
        typename std::enable_if<is_list_compatible<T>::value>::type>
      {
        static constexpr int list_rank =
          RankInfo<typename T::value_type>::list_rank + 1;
        static constexpr int map_rank =
          RankInfo<typename T::value_type>::map_rank;
      };

      // Rank of map types
      template <class T>
      struct RankInfo<
        T,
        typename std::enable_if<is_map_compatible<T>::value>::type>
      {
        static constexpr int list_rank =
          max_list_rank<typename T::key_type, typename T::mapped_type>() + 1;
        static constexpr int map_rank =
          max_map_rank<typename T::key_type, typename T::mapped_type>() + 1;
      };

      // Rank of Tensor types
      template <int rank, int dim, class Number>
      struct RankInfo<Tensor<rank, dim, Number>>
      {
        static constexpr int list_rank = rank + RankInfo<Number>::list_rank;
        static constexpr int map_rank  = RankInfo<Number>::map_rank;
      };

      template <int dim, class Number>
      struct RankInfo<Point<dim, Number>> : RankInfo<Tensor<1, dim, Number>>
      {};

      // Rank of complex types
      template <class Number>
      struct RankInfo<std::complex<Number>>
      {
        static constexpr int list_rank = RankInfo<Number>::list_rank + 1;
        static constexpr int map_rank  = RankInfo<Number>::map_rank;
      };

      // Rank of FunctionParser
      template <int dim>
      struct RankInfo<std::unique_ptr<FunctionParser<dim>>>
      {
        static constexpr int list_rank = 1;
        static constexpr int map_rank  = 0;
      };

      // Rank of ComponentMask
      template <>
      struct RankInfo<ComponentMask>
      {
        static constexpr int list_rank = 1;
        static constexpr int map_rank  = 0;
      };

      // Rank of std::pair
      template <class Key, class Value>
      struct RankInfo<std::pair<Key, Value>>
      {
        static constexpr int list_rank =
          std::max(RankInfo<Key>::list_rank, RankInfo<Value>::list_rank);
        static constexpr int map_rank =
          std::max(RankInfo<Key>::map_rank, RankInfo<Value>::map_rank) + 1;
      };


      template <class... Types>
      struct RankInfo<std::tuple<Types...>>
      {
        static constexpr int list_rank = max_list_rank<Types...>();
        static constexpr int map_rank  = max_map_rank<Types...>() + 1;
      };
    } // namespace internal

    // stl containers
    template <class T>
    struct Convert<T,
                   typename std::enable_if<is_list_compatible<T>::value>::type>
    {
      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        static_assert(internal::RankInfo<T>::list_rank > 0,
                      "Cannot use this class for non List-compatible types.");
        return std::make_unique<Patterns::List>(
          *Convert<typename T::value_type>::to_pattern(),
          0,
          std::numeric_limits<unsigned int>::max(),
          internal::default_list_separator[internal::RankInfo<T>::list_rank -
                                           1]);
      }

      static std::string
      to_string(
        const T &                    t,
        const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        auto p = dynamic_cast<const Patterns::List *>(&pattern);
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a "
                               "string to a List type."));
        auto                     base_p = p->get_base_pattern().clone();
        std::vector<std::string> vec(t.size());

        std::transform(
          t.cbegin(), t.cend(), vec.begin(), [&base_p](const auto &entry) {
            return Convert<typename T::value_type>::to_string(entry, *base_p);
          });

        std::string s;
        if (vec.size() > 0)
          s = vec[0];
        for (unsigned int i = 1; i < vec.size(); ++i)
          s += p->get_separator() + " " + vec[i];

        AssertThrow(pattern.match(s), ExcNoMatch(s, p->description()));
        return s;
      }

      static T
      to_value(const std::string &          s,
               const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        AssertThrow(pattern.match(s), ExcNoMatch(s, pattern.description()));

        auto p = dynamic_cast<const Patterns::List *>(&pattern);
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a string "
                               "to a List type."));

        auto base_p = p->get_base_pattern().clone();
        T    t;

        auto v = Utilities::split_string_list(s, p->get_separator());
        for (const auto &str : v)
          t.insert(t.end(),
                   Convert<typename T::value_type>::to_value(str, *base_p));

        return t;
      }
    };

    // stl maps
    template <class T>
    struct Convert<T,
                   typename std::enable_if<is_map_compatible<T>::value>::type>
    {
      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        static_assert(internal::RankInfo<T>::list_rank > 0,
                      "Cannot use this class for non List-compatible types.");
        static_assert(internal::RankInfo<T>::map_rank > 0,
                      "Cannot use this class for non Map-compatible types.");
        return std::make_unique<Patterns::Map>(
          *Convert<typename T::key_type>::to_pattern(),
          *Convert<typename T::mapped_type>::to_pattern(),
          0,
          std::numeric_limits<unsigned int>::max(),
          internal::default_list_separator[internal::RankInfo<T>::list_rank -
                                           1],
          internal::default_map_separator[internal::RankInfo<T>::map_rank - 1]);
      }

      static std::string
      to_string(
        const T &                    t,
        const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        auto p = dynamic_cast<const Patterns::Map *>(&pattern);
        AssertThrow(p,
                    ExcMessage("I need a Map pattern to convert a string to "
                               "a Map compatible type."));
        auto                     key_p = p->get_key_pattern().clone();
        auto                     val_p = p->get_value_pattern().clone();
        std::vector<std::string> vec(t.size());

        unsigned int i = 0;
        for (const auto &ti : t)
          vec[i++] =
            Convert<typename T::key_type>::to_string(ti.first, *key_p) +
            p->get_key_value_separator() +
            Convert<typename T::mapped_type>::to_string(ti.second, *val_p);

        std::string s;
        if (vec.size() > 0)
          s = vec[0];
        for (unsigned int i = 1; i < vec.size(); ++i)
          s += p->get_separator() + " " + vec[i];

        AssertThrow(p->match(s), ExcNoMatch(s, p->description()));
        return s;
      }

      static T
      to_value(const std::string &          s,
               const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        AssertThrow(pattern.match(s), ExcNoMatch(s, pattern.description()));

        auto p = dynamic_cast<const Patterns::Map *>(&pattern);
        AssertThrow(p,
                    ExcMessage("I need a Map pattern to convert a "
                               "string to a Map compatible type."));

        auto key_p = p->get_key_pattern().clone();
        auto val_p = p->get_value_pattern().clone();
        T    t;

        auto v = Utilities::split_string_list(s, p->get_separator());
        for (const auto &str : v)
          {
            auto key_val =
              Utilities::split_string_list(str, p->get_key_value_separator());
            AssertDimension(key_val.size(), 2);
            t.insert(std::make_pair(
              Convert<typename T::key_type>::to_value(key_val[0], *key_p),
              Convert<typename T::mapped_type>::to_value(key_val[1])));
          }

        return t;
      }
    };

    // Tensors
    template <int rank, int dim, class Number>
    struct Convert<Tensor<rank, dim, Number>>
    {
      using T = Tensor<rank, dim, Number>;
      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        static_assert(internal::RankInfo<T>::list_rank > 0,
                      "Cannot use this class for non List-compatible types.");
        return std::make_unique<Patterns::List>(
          *Convert<typename T::value_type>::to_pattern(),
          dim,
          dim,
          internal::default_list_separator[internal::RankInfo<T>::list_rank -
                                           1]);
      }

      static std::string
      to_string(
        const T &                    t,
        const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        auto p = dynamic_cast<const Patterns::List *>(&pattern);
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a string "
                               "to a List compatible type."));
        auto                     base_p = p->get_base_pattern().clone();
        std::vector<std::string> vec(dim);

        for (unsigned int i = 0; i < dim; ++i)
          vec[i] = Convert<typename T::value_type>::to_string(t[i], *base_p);

        std::string s;
        if (vec.size() > 0)
          s = vec[0];
        for (unsigned int i = 1; i < vec.size(); ++i)
          s += p->get_separator() + " " + vec[i];

        AssertThrow(p->match(s), ExcNoMatch(s, p->description()));
        return s;
      }

      static T
      to_value(const std::string &          s,
               const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        AssertThrow(pattern.match(s), ExcNoMatch(s, pattern.description()));

        auto p = dynamic_cast<const Patterns::List *>(&pattern);
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a string "
                               "to a List compatible type."));

        auto base_p = p->get_base_pattern().clone();
        T    t;

        auto         v = Utilities::split_string_list(s, p->get_separator());
        unsigned int i = 0;
        for (const auto &str : v)
          t[i++] = Convert<typename T::value_type>::to_value(str, *base_p);

        return t;
      }
    };

    // Points
    template <int dim, class Number>
    struct Convert<Point<dim, Number>>
    {
      using T = Point<dim, Number>;

      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        return Convert<Tensor<1, dim, Number>>::to_pattern();
      }

      static std::string
      to_string(
        const T &                    t,
        const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        return Convert<Tensor<1, dim, Number>>::to_string(
          Tensor<1, dim, Number>(t), pattern);
      }

      static T
      to_value(const std::string &          s,
               const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        return T(Convert<Tensor<1, dim, Number>>::to_value(s, pattern));
      }
    };

    // Functions::FunctionParser
    template <int dim>
    struct Convert<std::unique_ptr<FunctionParser<dim>>>
    {
      using T = std::unique_ptr<FunctionParser<dim>>;

      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        static_assert(internal::RankInfo<T>::list_rank > 0,
                      "Cannot use this class for non List-compatible types.");

        return std::make_unique<Patterns::List>(
          Patterns::Anything(),
          1,
          Patterns::List::max_int_value,
          internal::default_list_separator[internal::RankInfo<T>::list_rank -
                                           1]);
      }

      static std::string
      to_string(
        const T &                    t,
        const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        auto p = dynamic_cast<const Patterns::List *>(&pattern);
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a string "
                               "to a List compatible type."));

        const auto &expressions = t->get_expressions();
        if (expressions.size() == 0)
          return std::string();

        std::string s = expressions[0];
        for (unsigned int i = 1; i < expressions.size(); ++i)
          s = s + p->get_separator() + expressions[i];

        AssertThrow(pattern.match(s), ExcNoMatch(s, p->description()));
        return s;
      }

      static T
      to_value(const std::string &          s,
               const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        AssertThrow(pattern.match(s), ExcNoMatch(s, pattern.description()));

        auto p = dynamic_cast<const Patterns::List *>(&pattern);
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a string "
                               "to a List compatible type."));

        const auto expressions =
          Utilities::split_string_list(s, p->get_separator());

        T t = std::make_unique<FunctionParser<dim>>(expressions.size());
        const std::string var =
          FunctionParser<dim>::default_variable_names() + ",t";
        const typename FunctionParser<dim>::ConstMap constants;
        t->initialize(var, expressions, constants, true);
        return t;
      }
    };

    // ComponentMask
    template <>
    struct Convert<ComponentMask>
    {
      using T = ComponentMask;

      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        return Convert<std::vector<bool>>::to_pattern();
      }

      static std::string
      to_string(
        const T &                    t,
        const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        std::vector<bool> mask(t.size());
        for (unsigned int i = 0; i < t.size(); ++i)
          mask[i] = t[i];

        return Convert<std::vector<bool>>::to_string(mask, pattern);
      }

      static T
      to_value(const std::string &          s,
               const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        const auto mask = Convert<std::vector<bool>>::to_value(s, pattern);
        return ComponentMask(mask);
      }
    };

    // Complex numbers
    template <class Number>
    struct Convert<std::complex<Number>>
    {
      using T = std::complex<Number>;

      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        static_assert(internal::RankInfo<T>::list_rank > 0,
                      "Cannot use this class for non List-compatible types.");
        return std::make_unique<Patterns::List>(
          *Convert<typename T::value_type>::to_pattern(),
          2,
          2,
          internal::default_list_separator[internal::RankInfo<T>::list_rank -
                                           1]);
      }

      static std::string
      to_string(
        const T &                    t,
        const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        auto p = dynamic_cast<const Patterns::List *>(&pattern);
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a string "
                               "to a List compatible type."));

        auto        base_p = p->get_base_pattern().clone();
        std::string s =
          Convert<typename T::value_type>::to_string(t.real(), *base_p) +
          p->get_separator() + " " +
          Convert<typename T::value_type>::to_string(t.imag(), *base_p);

        AssertThrow(pattern.match(s), ExcNoMatch(s, p->description()));
        return s;
      }

      /**
       * 使用给定的模式或默认模式，将字符串转换为数值。
       *
       */
      static T
      to_value(const std::string &          s,
               const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        AssertThrow(pattern.match(s), ExcNoMatch(s, pattern.description()));

        auto p = dynamic_cast<const Patterns::List *>(&pattern);
        AssertThrow(p,
                    ExcMessage("I need a List pattern to convert a string "
                               "to a List compatible type."));

        auto base_p = p->get_base_pattern().clone();

        auto v = Utilities::split_string_list(s, p->get_separator());
        AssertDimension(v.size(), 2);
        T t(Convert<typename T::value_type>::to_value(v[0], *base_p),
            Convert<typename T::value_type>::to_value(v[1], *base_p));
        return t;
      }
    };

    // Strings
    template <>
    struct Convert<std::string>
    {
      using T = std::string;

      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        return std::make_unique<Patterns::Anything>();
      }

      static std::string
      to_string(
        const T &                    t,
        const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        AssertThrow(pattern.match(t), ExcNoMatch(t, pattern.description()));
        return t;
      }

      static T
      to_value(const std::string &          s,
               const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        AssertThrow(pattern.match(s), ExcNoMatch(s, pattern.description()));
        return s;
      }
    };

    // Pairs
    template <class Key, class Value>
    struct Convert<std::pair<Key, Value>>
    {
      using T = std::pair<Key, Value>;

      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        static_assert(internal::RankInfo<T>::map_rank > 0,
                      "Cannot use this class for non Map-compatible types.");
        return std::make_unique<Patterns::Tuple>(
          internal::default_map_separator[internal::RankInfo<T>::map_rank - 1],
          *Convert<Key>::to_pattern(),
          *Convert<Value>::to_pattern());
      }

      static std::string
      to_string(
        const T &                    t,
        const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        std::tuple<Key, Value> m(t);
        std::string            s = Convert<decltype(m)>::to_string(m, pattern);
        AssertThrow(pattern.match(s), ExcNoMatch(s, pattern.description()));
        return s;
      }

      static T
      to_value(const std::string &          s,
               const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        std::tuple<Key, Value> m;
        m = Convert<decltype(m)>::to_value(s, pattern);
        return std::make_pair(std::get<0>(m), std::get<1>(m));
      }
    };

    // Tuples
    template <class... Args>
    struct Convert<std::tuple<Args...>>
    {
      using T = std::tuple<Args...>;

      static std::unique_ptr<Patterns::PatternBase>
      to_pattern()
      {
        static_assert(internal::RankInfo<T>::map_rank > 0,
                      "Cannot use this class for non tuple-compatible types.");
        return std::make_unique<Patterns::Tuple>(
          internal::default_map_separator[internal::RankInfo<T>::map_rank - 1],
          *Convert<Args>::to_pattern()...);
      }

      static std::string
      to_string(
        const T &                    t,
        const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        auto p = dynamic_cast<const Patterns::Tuple *>(&pattern);
        AssertThrow(p,
                    ExcMessage("I need a Tuple pattern to convert a tuple "
                               "to a string."));

        const auto  string_array = Convert<T>::to_string_internal_2(t, *p);
        std::string str;
        for (unsigned int i = 0; i < string_array.size(); ++i)
          str += (i ? " " + p->get_separator() + " " : "") + string_array[i];
        AssertThrow(p->match(str), ExcNoMatch(str, p->description()));
        return str;
      }

      static T
      to_value(const std::string &          s,
               const Patterns::PatternBase &pattern = *Convert<T>::to_pattern())
      {
        AssertThrow(pattern.match(s), ExcNoMatch(s, pattern.description()));

        auto p = dynamic_cast<const Patterns::Tuple *>(&pattern);
        AssertThrow(p,
                    ExcMessage("I need a Tuple pattern to convert a string "
                               "to a tuple type."));

        auto v = Utilities::split_string_list(s, p->get_separator());

        return Convert<T>::to_value_internal_2(v, *p);
      }

    private:
      template <std::size_t... U>
      static std::array<std::string, std::tuple_size<T>::value>
      to_string_internal_1(const T &              t,
                           const Patterns::Tuple &pattern,
                           std::index_sequence<U...>)
      {
        std::array<std::string, std::tuple_size<T>::value> a = {
          {Convert<typename std::tuple_element<U, T>::type>::to_string(
            std::get<U>(t), pattern.get_pattern(U))...}};
        return a;
      }

      static std::array<std::string, std::tuple_size<T>::value>
      to_string_internal_2(const T &t, const Patterns::Tuple &pattern)
      {
        return Convert<T>::to_string_internal_1(
          t, pattern, std::make_index_sequence<std::tuple_size<T>::value>{});
      }

      template <std::size_t... U>
      static T
      to_value_internal_1(const std::vector<std::string> &s,
                          const Patterns::Tuple &         pattern,
                          std::index_sequence<U...>)
      {
        return std::make_tuple(
          Convert<typename std::tuple_element<U, T>::type>::to_value(
            s[U], pattern.get_pattern(U))...);
      }

      static T
      to_value_internal_2(const std::vector<std::string> &s,
                          const Patterns::Tuple &         pattern)
      {
        return Convert<T>::to_value_internal_1(
          s, pattern, std::make_index_sequence<std::tuple_size<T>::value>{});
      }
    };

    // Utility function with default Pattern
    template <typename T>
    std::string
    to_string(const T &t)
    {
      return Convert<T>::to_string(t);
    }

    // Utility function with default Pattern
    template <typename T>
    void
    to_value(const std::string &s, T &t)
    {
      t = Convert<T>::to_value(s);
    }
  } // namespace Tools
} // namespace Patterns


DEAL_II_NAMESPACE_CLOSE

#endif


