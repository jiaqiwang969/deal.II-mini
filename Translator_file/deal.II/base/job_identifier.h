//include/deal.II-translator/base/job_identifier_0.txt
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

#ifndef dealii_job_identifier_h
#define dealii_job_identifier_h


#include <deal.II/base/config.h>

#include <string>

DEAL_II_NAMESPACE_OPEN
/**
 * 识别一个程序的运行。<tt>JobIdentifier</tt>确定一个程序运行的开始时间，并将其作为一个程序标识符存储。存在一个该类的库对象<tt>dealjobid</tt>。这个对象可以被所有的输出函数访问，为当前的工作提供一个ID。
 *
 *
 * @ingroup utilities
 *
 *
 */
class JobIdentifier
{
public:
  /**
   * 构造函数。设置程序标识符为<tt>program_id</tt>与当前时间相连接的值。
   *
   */
  JobIdentifier();

  /**
   * 这个函数返回运行中的程序的标识符。目前，该库提供了一个返回
   * "JobID "的函数。
   * 用户可以在他们的源代码中定义这个函数的替换，避免链接库的版本。不幸的是，这种机制对共享库不起作用。
   *
   */
  static const char *
  program_id();

  /**
   * 获得作为参数传递的文件名的基本名称。也就是说，如果文件是<tt>mypath/file.cc</tt>，只返回<tt>file</tt>。例如，这个函数可以从一个用户程序中调用，参数为<tt>__FILE__</tt>，为正在运行的程序创建一个标识符。
   *
   */
  static std::string
  base_name(const std::string &filename);

  /**
   * 返回<tt>id</tt>的值。
   *
   */
  const std::string
  operator()() const;

  /**
   * 识别当前运行程序的%函数。
   *
   */
  static const JobIdentifier &
  get_dealjobid();

private:
  /**
   * 持有目前正在运行的程序的标识符的字符串。
   *
   */
  std::string id;
};

DEAL_II_NAMESPACE_CLOSE

#endif


