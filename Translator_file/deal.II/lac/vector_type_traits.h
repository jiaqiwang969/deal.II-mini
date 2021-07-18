//include/deal.II-translator/lac/vector_type_traits_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

#ifndef dealii_vector_type_traits_h
#define dealii_vector_type_traits_h

#include <deal.II/base/config.h>

#include <type_traits>


DEAL_II_NAMESPACE_OPEN


/**
 * 串行向量的类型特质，即不支持在进程中分配存储的向量类。
 * 专门化
 *
 * @code
 * template <>
 * struct is_serial_vector<VectorType> : std::true_type
 * {};
 * @endcode
 * 分别为一个串行向量类型。
 *
 * @code
 * template <>
 * struct is_serial_vector<VectorType> : std::false_type
 * {};
 * @endcode
 * 对于支持分布式存储的向量类型，必须在向量声明的头文件中完成。
 *
 *
 */
template <typename T>
struct is_serial_vector;


DEAL_II_NAMESPACE_CLOSE

#endif


