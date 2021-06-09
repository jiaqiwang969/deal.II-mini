//include/deal.II-translator/lac/parallel_block_vector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2020 by the deal.II authors
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

#ifndef dealii_parallel_block_vector_h
#define dealii_parallel_block_vector_h

#include <deal.II/base/config.h>

#include <deal.II/lac/la_parallel_block_vector.h>

DEAL_II_WARNING(
  "This file is deprecated. Use <deal.II/lac/la_parallel_block_vector.h> and LinearAlgebra::distributed::BlockVector instead.")

#include <cstring>
#include <iomanip>


DEAL_II_NAMESPACE_OPEN


namespace parallel
{
  namespace distributed
  {
    /*!   @addtogroup Vectors  
     * @{     
*
*/

    /**
     * 一个基于分布式deal.II向量的块向量的实现。虽然基类提供了大部分的接口，但这个类处理了向量的实际分配，并提供了底层向量类型的特定函数。
     * @note  这个模板的实例提供给<tt>  @<float@>  和  @<double@></tt>;  其他可以在应用程序中生成（见手册中的 @ref Instantiations 部分）。          @see   @ref GlossBlockLA  "块（线性代数）"
     * @deprecated  用 LinearAlgebra::distributed::BlockVector 代替。
     *
     */
    template <typename Number>
    using BlockVector DEAL_II_DEPRECATED =
      LinearAlgebra::distributed::BlockVector<Number>;

     /*@}*/ 
  } // namespace distributed
} // namespace parallel

DEAL_II_NAMESPACE_CLOSE

#endif


