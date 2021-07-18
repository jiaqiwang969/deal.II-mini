//include/deal.II-translator/algorithms/general_data_storage_0.txt
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

#ifndef dealii_algorithms_general_data_storage_h
#define dealii_algorithms_general_data_storage_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>

#include <boost/any.hpp>
#include <boost/core/demangle.hpp>

#include <algorithm>
#include <map>
#include <string>
#include <typeinfo>

DEAL_II_NAMESPACE_OPEN


/**
 * 这个类方便了任何一般数据的存储。它提供了一个机制来存储任何数量的数据，任何类型的数据，然后通过一个标识符字符串来访问。
 * 当使用这个类时，请引用
 *
 *
 * @code{.bib}
 * @article{SartoriGiulianiBardelloni-2018-a,
 * Author = {Sartori, Alberto and Giuliani, Nicola and
 *          Bardelloni, Mauro and Heltai, Luca},
 * Journal = {SoftwareX},
 * Pages = {318--327},
 * Title = {{deal2lkit: A toolkit library for high performance
 *          programming in deal.II}},
 * Volume = {7},
 * Year = {2018}}
 * @endcode
 *
 *
 *
 */
class GeneralDataStorage : public Subscriptor
{
public:
  /**
   * 默认构造函数。
   *
   */
  GeneralDataStorage() = default;

  /**
   * 复制构造函数。
   *
   */
  GeneralDataStorage(const GeneralDataStorage &) = default;

  /**
   * 移动构造函数。
   *
   */
  GeneralDataStorage(GeneralDataStorage &&) = default;

  /**
   * 这个类实例所存储的对象的数量。
   *
   */
  std::size_t
  size() const;

  /**
   * 将 @p other_data 的内容与此对象合并。
   *
   */
  void
  merge(const GeneralDataStorage &other_data);

  /**
   * 打印内部缓存的内容到 @p stream.   @p any_data
   * 映射中的每个键和值对都打印在单独的一行，
   * <tt>std::string</tt>
   * 键列在前面，后面是相关映射类型的拆分<tt>type_id</tt>。
   *
   */
  template <class Stream>
  void
  print_info(Stream &stream);

  /**
   * 清除存储在这个类实例中的所有数据。
   * 当你调用这个函数时，它会销毁所有你要求以副本形式存储的对象，并且忘记你要求以引用形式存储的数据的引用。因此，你现在可以在任何时候自由地销毁你所存储的引用的对象。
   *
   * - 在当前的`GeneralDataStorage`对象被销毁之前或之后。    为了澄清这一点，请考虑下面这个小例子。
   * @code
   * GeneralDataStorage data;
   *
   * {
   *   const double some_number = ...;
   *   data.add_unique_reference("value", some_number);
   *
   *   // Adding either of these next two lines could fix the
   *   // issue, by removing the association of some_number with data:
   *   // data.remove_object_with_name("value");
   *   // data.reset();
   * } // some_number goes out of scope here
   *
   * const double some_other_number
   *   = data.get_object_with_name<double>("value"); // Invalid call
   * @endcode
   * 在上面的代码中， @p data
   * 对象的范围比<tt>some_number</tt>长。当我们从 @p data
   * 中获取<tt>"value"</tt>时，对 @p some_number
   * 的引用已不再有效。
   * 同样地，对于复制到GeneralDataStorage对象中的数据，我们应该考虑它在哪个范围内保持有效。
   * @code
   * double* ptr_to_some_number = null_ptr;
   *
   * {
   *   GeneralDataStorage data;
   *   const double some_number = ...;
   *   data.add_unique_copy("value", some_number);
   *
   *   ptr_to_some_number = &(data.get_object_with_name<double>("value"));
   * } // The copy to some_number goes out of scope here
   *
   * const double some_other_number
   *   =ptr_to_some_number; // Invalid call
   * @endcode
   * 与第一个例子类似，我们必须意识到这样一个事实：由
   * @p data 存储的任何 @p Type
   * 的副本只有在存储它的GeneralDataStorage实例活着时才有效。
   * 此外，正如上一个例子所阐明的，当调用reset()函数时，或通过调用remove_object_with_name()删除时，被指向的类实例（由GeneralDataStorage拥有）的副本就不再是活的。
   * @code
   * GeneralDataStorage data;
   * double* ptr_to_some_number = null_ptr;
   *
   * {
   *   const double some_number = ...;
   *   data.add_unique_copy("value", some_number);
   *
   *   ptr_to_some_number = &(data.get_object_with_name<double>("value"));
   *
   *   // The copy to some_number would go out of scope when either of
   *   // following two calls are made:
   *   data.remove_object_with_name("value");
   *   data.reset();
   * }
   *
   * const double some_other_number
   *   =ptr_to_some_number; // Invalid call
   * @endcode
   *
   *
   */
  void
  reset();

  /**
   * @name  数据存储和访问
   *
   */
  //@{

  /**
   * 在内部存储一个给定对象的副本。被复制的对象为这个类所拥有，并且可以通过get_object_with_name()方法通过引用访问。
   * 这个函数确保没有 @p entry 与给定的 @p name
   * 已经被这个类实例所存储。
   *
   */
  template <typename Type>
  void
  add_unique_copy(const std::string &name, const Type &entry);

  /**
   * 在内部存储一个给定对象的副本。被复制的对象为这个类所拥有，并且可以通过引用get_object_with_name()方法来访问。
   * 这个函数不进行任何检查以确保具有给定 @p name 的 @p
   * entry 已经被这个类实例存储。如果 @p name
   * 确实指向现有的数据，那么这将被覆盖。
   *
   */
  template <typename Type>
  void
  add_or_overwrite_copy(const std::string &name, const Type &entry);

  /**
   * 为一个已经存在的对象添加一个引用。该对象不属于这个类，用户必须保证被引用的对象比这个类实例的寿命长。存储的引用可以通过get_object_with_name()方法访问。
   * 这个函数确保没有任何具有给定 @p name 的 @p entry
   * 已经被这个类实例所存储。
   *
   */
  template <typename Type>
  void
  add_unique_reference(const std::string &name, Type &entry);

  /**
   * 为一个已经存在的对象添加一个引用。该对象不属于这个类，而且用户必须保证被引用的对象比这个类实例的寿命长。存储的引用可以通过get_object_with_name()方法访问。
   * 这个函数不执行任何检查以确保具有给定 @p name 的 @p
   * entry 已经被这个类实例所存储。如果 @p name
   * 确实指向现有的数据，那么这将被覆盖。
   *
   */
  template <typename Type>
  void
  add_or_overwrite_reference(const std::string &name, Type &entry);

  /**
   * 返回一个给定名称的对象的引用。如果该对象不存在，那么输入的
   * @p arguments 将被用来构造一个给定的 @p Type
   * 的对象，然后返回这个新对象的引用。    一个 @p Type
   * 类型的对象的副本，是由这个类实例拥有的，通过调用它的构造函数和给定的参数集来生成。对于这个函数，
   * @p arguments 作为<tt>lvalue</tt>引用被传递。
   *
   */
  template <typename Type, typename Arg, typename... Args>
  Type &
  get_or_add_object_with_name(const std::string &name,
                              Arg &              argument,
                              Args &... arguments);

  /**
   * 返回一个给定名称的对象的引用。如果该对象不存在，那么输入的
   * @p arguments 将被用来构造一个给定的 @p Type
   * 的对象，然后返回这个新对象的引用。
   * 与上述单参数的情况相同。
   *
   */
  template <typename Type, typename Arg>
  Type &
  get_or_add_object_with_name(const std::string &name, Arg &argument);

  /**
   * 返回一个对给定名称的对象的引用。如果该对象不存在，那么输入的
   * @p arguments 将被用来构造一个给定的 @p Type
   * 的对象，然后返回这个新对象的引用。    一个 @p Type
   * 类型的对象的副本，是由这个类实例所拥有的，通过调用它的构造函数和给定的参数集来生成。与之前的同名函数不同的是，对于这个函数，
   * @p arguments 是作为<tt>rvalue</tt>引用传递的。
   *
   */
  template <typename Type, typename Arg, typename... Args>
  Type &
  get_or_add_object_with_name(const std::string &name,
                              Arg &&             argument,
                              Args &&... arguments);

  /**
   * 返回一个给定名称的对象的引用。如果该对象不存在，那么输入的
   * @p arguments 将被用来构造一个给定的 @p Type
   * 的对象，然后返回这个新对象的引用。
   * 与上述单参数的情况相同。
   *
   */
  template <typename Type, typename Arg>
  Type &
  get_or_add_object_with_name(const std::string &name, Arg &&argument);

  /**
   * 同上面的默认构造函数。
   *
   */
  template <typename Type>
  Type &
  get_or_add_object_with_name(const std::string &name);

  /**
   * 返回一个给定名称的对象的引用。
   * 如果具有给定名称的对象没有存储在这个类中，或者具有给定名称的对象既不是精确指定的
   * @p Type 也不是指向 @p Type.
   * 的指针，这个函数会抛出一个异常。
   *
   */
  template <typename Type>
  Type &
  get_object_with_name(const std::string &name);

  /**
   * 返回一个给定名称的对象的常数引用。
   * 如果具有给定名称的对象没有存储在这个类中，或者如果具有给定名称的对象既不是精确指定的
   * @p Type 也不是指向 @p Type.
   * 的指针，这个函数会抛出一个异常。
   *
   */
  template <typename Type>
  const Type &
  get_object_with_name(const std::string &name) const;

  /**
   * 找出我们是否存储了一个具有给定名称的对象。
   *
   */
  bool
  stores_object_with_name(const std::string &name) const;

  /**
   * 删除具有给定名称的对象。
   *
   */
  void
  remove_object_with_name(const std::string &name);

  //@}

  /**
   * 在内部 boost::any 地图中不存在一个具有此名称的条目。
   *
   */
  DeclException1(ExcNameNotFound,
                 std::string,
                 << "No entry with the name " << arg1 << " exists.");

  /**
   * 在内部 boost::any 地图中不存在具有此名称的条目。
   *
   */
  DeclException1(ExcNameHasBeenFound,
                 std::string,
                 << "An entry with the name " << arg1 << " already exists.");

  /**
   * 请求的类型和存储的类型不同。
   *
   */
  DeclException3(ExcTypeMismatch,
                 std::string,
                 const char *,
                 const char *,
                 << "The stored type for entry with name \"" << arg1 << "\" is "
                 << arg2 << " but you requested type " << arg3 << ".");

private:
  /**
   * 任意的用户数据，由一个字符串标识。
   *
   */
  std::map<std::string, boost::any> any_data;
};


 /*----------------- Inline and template methods -----------------*/ 


#ifndef DOXYGEN


template <class Stream>
void
GeneralDataStorage::print_info(Stream &os)
{
  for (const auto &it : any_data)
    {
      os << it.first << '\t' << '\t'
         << boost::core::demangle(it.second.type().name()) << std::endl;
    }
}


template <typename Type>
void
GeneralDataStorage::add_unique_copy(const std::string &name, const Type &entry)
{
  AssertThrow(!stores_object_with_name(name), ExcNameHasBeenFound(name));
  add_or_overwrite_copy(name, entry);
}


template <typename Type>
void
GeneralDataStorage::add_or_overwrite_copy(const std::string &name,
                                          const Type &       entry)
{
  any_data[name] = entry;
}


template <typename Type>
void
GeneralDataStorage::add_unique_reference(const std::string &name, Type &entry)
{
  AssertThrow(!stores_object_with_name(name), ExcNameHasBeenFound(name));
  add_or_overwrite_reference(name, entry);
}


template <typename Type>
void
GeneralDataStorage::add_or_overwrite_reference(const std::string &name,
                                               Type &             entry)
{
  Type *ptr      = &entry;
  any_data[name] = ptr;
}


template <typename Type>
Type &
GeneralDataStorage::get_object_with_name(const std::string &name)
{
  AssertThrow(stores_object_with_name(name), ExcNameNotFound(name));

  Type *p = nullptr;

  if (any_data[name].type() == typeid(Type *))
    {
      p = boost::any_cast<Type *>(any_data[name]);
    }
  else if (any_data[name].type() == typeid(Type))
    {
      p = boost::any_cast<Type>(&any_data[name]);
    }
  else
    {
      AssertThrow(false,
                  ExcTypeMismatch(name,
                                  any_data[name].type().name(),
                                  typeid(Type).name()));
    }

  return *p;
}


template <typename Type>
const Type &
GeneralDataStorage::get_object_with_name(const std::string &name) const
{
  AssertThrow(stores_object_with_name(name), ExcNameNotFound(name));

  const auto it = any_data.find(name);

  if (it->second.type() == typeid(Type *))
    {
      const Type *p = boost::any_cast<Type *>(it->second);
      return *p;
    }
  else if (it->second.type() == typeid(Type))
    {
      const Type *p = boost::any_cast<Type>(&it->second);
      return *p;
    }
  else
    {
      AssertThrow(false,
                  ExcTypeMismatch(name,
                                  it->second.type().name(),
                                  typeid(Type).name()));
      const Type *p = nullptr;
      return *p;
    }
}



template <typename Type, typename Arg>
Type &
GeneralDataStorage::get_or_add_object_with_name(const std::string &name,
                                                Arg &              argument)
{
  if (!stores_object_with_name(name))
    add_unique_copy(name, Type(argument));

  return get_object_with_name<Type>(name);
}



template <typename Type, typename Arg, typename... Args>
Type &
GeneralDataStorage::get_or_add_object_with_name(const std::string &name,
                                                Arg &              argument,
                                                Args &... arguments)
{
  if (!stores_object_with_name(name))
    add_unique_copy(name, Type(argument, arguments...));

  return get_object_with_name<Type>(name);
}



template <typename Type, typename Arg>
Type &
GeneralDataStorage::get_or_add_object_with_name(const std::string &name,
                                                Arg &&             argument)
{
  if (!stores_object_with_name(name))
    add_unique_copy(name, Type(std::forward<Arg>(argument)));

  return get_object_with_name<Type>(name);
}



template <typename Type, typename Arg, typename... Args>
Type &
GeneralDataStorage::get_or_add_object_with_name(const std::string &name,
                                                Arg &&             argument,
                                                Args &&... arguments)
{
  if (!stores_object_with_name(name))
    add_unique_copy(name,
                    Type(std::forward<Arg>(argument),
                         std::forward<Args>(arguments)...));

  return get_object_with_name<Type>(name);
}


template <typename Type>
Type &
GeneralDataStorage::get_or_add_object_with_name(const std::string &name)
{
  if (!stores_object_with_name(name))
    add_unique_copy(name, Type());

  return get_object_with_name<Type>(name);
}


#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_algorithms_general_data_storage_h


