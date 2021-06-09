//include/deal.II-translator/lac/vector_operation_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

#ifndef dealii_lac_vector_operation_h
#define dealii_lac_vector_operation_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

/*!   @addtogroup Vectors  
     * @{  ．

 
*
*/

/**
 * 这个枚举记录了像向量和矩阵这样的并行线性代数对象中的当前操作。
 * 它被用于各种compress()函数中。为了兼容，它们也存在于串行代码中，在那里是空的。
 * 参见 @ref GlossCompress  "压缩分布式对象 "
 * 以获得更多信息。
 *
 *
 */
struct VectorOperation
{
  enum values
  {
    /**
     * 当前的操作是未知的。
     *
     */
    unknown,
    /**
     * 当前操作是一个插入。
     *
     */
    insert,
    /**
     * 当前的操作是一个加法。
     *
     */
    add,
    /**
     * 当前的操作是一个最小化操作。
     *
     */
    min,
    /**
     * 当前的操作是一个最大化的操作。
     *
     */
    max
  };
};

 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


