//include/deal.II-translator/base/path_search_0.txt
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

#ifndef dealii_path_search_h
#define dealii_path_search_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>


DEAL_II_NAMESPACE_OPEN

/**
 * 支持在路径列表和后缀列表中搜索文件。
 * 为支持的每个文件类别维护一个搜索路径列表。一个文件类别是由一个唯一的字符串定义的。提供的类有：
 * <dl>  <dt> MESH  <dd>  各种格式的网格输入文件（见GridIn）
 * <dt> PARAMETER  <dd>  参数文件（<tt>.prm</tt>）  </dl>
 * 使用add_class()可以很容易地添加额外的文件类别。
 * 使用方法。首先，你为某个文件类构造一个PathSearch对象，比如说网格。然后，你使用find()方法获得一个完整的路径名称，你就可以打开文件了。
 *
 * @code
 * #include <deal.II/base/path_search.h>
 *
 * dealii::PathSearch search("MESH");
 * std::string full_name = search.find("grid");
 * std::ifstream in(full_name.c_str());
 * ...
 * @endcode
 *
 * 这段代码将首先遍历为文件类<tt>MESH</tt>设置的列表中的所有路径。如果它成功地打开了一个文件，它将返回<tt>istream</tt>对象。如果没有，它将尝试追加后缀列表中的第一个后缀，并进行同样的操作。以此类推。如果最后没有找到文件，会抛出一个异常。
 * 如果你想把你的搜索限制在某个网格格式，例如<tt>.inp</tt>，那么要么在上面的代码中使用<tt>"grid.inp"</tt>，要么使用替代的find(const
 * std::string&,const   std::string&,const  char*) 函数
 *
 * @code
 * std::string full_name = search.find("grid", ".inp");
 * @endcode
 *
 * 路径列表默认以当前目录（<tt>"./"<tt>）开始，后面可以选择以标准目录的deal.II开始。使用show()来找出一个给定类的路径列表。路径和后缀可以分别用函数add_path()和add_suffix()来添加。
 *
 *
 * @note
 * 路径列表中的目录应该总是以尾部的<tt>"/"<tt>结束，而后缀应该总是以点开始。这些字符不会被自动添加（允许你做一些真正的文件名编辑）。
 * @todo  增加对环境变量的支持，就像在kpathsea中一样。
 *
 *
 * @ingroup input
 *
 *
 */
class PathSearch
{
public:
  /**
   * 向列表中添加新项目的位置。
   *
   */
  enum Position
  {
    /// Add new item at end of list
    back,
    /// Add new item at front of list
    front,
    /// Add in path list after empty element
    after_none
  };

  /**
   * 构造函数。第一个参数是一个字符串，确定要搜索的文件的类别。
   * debug参数决定了该类文件的粗略程度。
   *
   */
  PathSearch(const std::string &cls, const unsigned int debug = 0);

  /**
   * 在构造函数指定的类中寻找一个文件，并返回其完整的路径名称（包括可能的后缀）。
   * 文件搜索的工作方式是实际尝试打开该文件。如果 @p
   * fopen 与提供的 @p open_mode,
   * 成功，那么文件就被找到了，否则搜索继续进行。
   * @warning  小心使用 @p open_mode!
   * 特别是，要非常小心地使用<tt>"w"</tt>!
   * 如果文件不存在，就无法找到它。如果它确实存在， @p
   * fopen 函数将把它截断成零长度。      @param  filename
   * 要找到的文件的基本名称，不包括路径成分和后缀。
   * @param  open_mode 移交给 @p fopen 函数的模式。
   *
   */
  std::string
  find(const std::string &filename, const char *open_mode = "r");

  /**
   * 在构造函数指定的类中找到一个文件，并返回其完整的路径名称。不使用标准的后缀列表，而只尝试应用给定的后缀。
   * 文件搜索的工作方式是实际尝试打开该文件。如果 @p
   * fopen 与所提供的 @p open_mode,
   * 成功，那么文件就被找到了，否则搜索继续进行。
   * @warning  小心使用 @p open_mode!
   * 特别是，要非常小心地使用<tt>"w"</tt>!
   * 如果文件不存在，就无法找到它。如果它确实存在， @p
   * fopen 函数将把它截断成零长度。      @param  filename
   * 要找到的文件的基本名称，不包括路径组件和后缀。
   * @param  后缀 打开时要使用的后缀。    @param  open_mode
   * 移交给  @p fopen  函数的模式。
   *
   */
  std::string
  find(const std::string &filename,
       const std::string &suffix,
       const char *       open_mode = "r");

  /**
   * 显示此对象使用的路径和后缀。
   *
   */
  template <class StreamType>
  void
  show(StreamType &stream) const;

  /**
   * 添加一个新的类。
   *
   */
  static void
  add_class(const std::string &cls);

  /**
   * 为当前的类添加一个路径。参见 PathSearch::Position
   * 中可能的位置参数。
   *
   */
  void
  add_path(const std::string &path, Position pos = back);

  /**
   * 为当前的类添加一个路径。见 PathSearch::Position
   * 中可能的位置参数。
   *
   */
  void
  add_suffix(const std::string &suffix, Position pos = back);

  /**
   * 这个类没有在路径搜索机制中注册。
   * @ingroup Exceptions
   *
   */
  DeclException1(ExcNoClass,
                 std::string,
                 << "The class " << arg1
                 << " must be registered before referring it in PathSearch.");
  /**
   * PathSearch类在其路径列表中找不到有此名称的文件。
   * @ingroup Exceptions
   *
   */
  DeclException2(ExcFileNotFound,
                 std::string,
                 std::string,
                 << "The file \"" << arg1 << "\" was not found in the path for "
                 << "files of class " << arg2 << ".");

private:
  /**
   * 类映射中的值的类型。
   *
   */
  using map_type = std::map<std::string, std::vector<std::string>>::value_type;

  /**
   * 初始化静态列表对象以便进一步使用。
   *
   */
  static void
  initialize_classes();

  /**
   * 获取某个类的路径列表。用来在构造函数中设置#my_path_list。
   *
   */
  static std::vector<std::string> &
  get_path_list(const std::string &cls);

  /**
   * 获取某类的后缀列表。用于在构造函数中设置#my_suffix_list。
   *
   */
  static std::vector<std::string> &
  get_suffix_list(const std::string &cls);

  /**
   * 这个对象所处理的文件类。
   *
   */
  const std::string cls;

  /**
   * 所有类的所有路径列表，这样我们就可以只建立一次。
   *
   */
  static std::map<std::string, std::vector<std::string>> path_lists;

  /**
   * 每个类的后缀列表。
   *
   */
  static std::map<std::string, std::vector<std::string>> suffix_lists;

  /**
   * 这个对象所属的类的路径列表。
   *
   */
  std::vector<std::string> &my_path_list;

  /**
   * 此对象所属类的后缀列表。
   *
   */
  std::vector<std::string> &my_suffix_list;

  /**
   * 调试标志。如果为零则没有输出。
   *
   */
  const unsigned int debug;

  /**
   * 空字符串。
   *
   */
  static std::string empty;
};


 /* ----------------------------- inline functions ------------------------- */ 


template <class StreamType>
inline void
PathSearch::show(StreamType &out) const
{
  out << "DEAL_II_" << cls << "PATH=\"";
  bool first = true;
  for (const auto &p : my_path_list)
    {
      if (!first)
        out << ':';
      out << p;
      first = false;
    }
  out << '"' << std::endl << " Suffixes";
  for (const auto &s : my_suffix_list)
    out << " \"" << s << '"';
  out << std::endl;
}

DEAL_II_NAMESPACE_CLOSE

#endif


