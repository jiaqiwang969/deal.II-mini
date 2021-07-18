//include/deal.II-translator/differentiation/sd_0.txt
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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii_differentiation_sd_h
#define dealii_differentiation_sd_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/differentiation/sd/symengine_math.h>
#  include <deal.II/differentiation/sd/symengine_number_traits.h>
#  include <deal.II/differentiation/sd/symengine_number_types.h>
#  include <deal.II/differentiation/sd/symengine_optimizer.h>
#  include <deal.II/differentiation/sd/symengine_product_types.h>
#  include <deal.II/differentiation/sd/symengine_scalar_operations.h>
#  include <deal.II/differentiation/sd/symengine_tensor_operations.h>
#  include <deal.II/differentiation/sd/symengine_types.h>
#  include <deal.II/differentiation/sd/symengine_utilities.h>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  /**
   * 符号微分库的封装器。目前对以下库的支持。
   *
   *
   *
   *
   * - SymEngine
   * @ingroup auto_symb_diff
   *
   */
  namespace SD
  {
    /**
     * 这个命名空间定义了有助于为符号数和操作提供结构化接口的类和函数。
     *
     */
    namespace internal
    {}
  } // namespace SD
} // namespace Differentiation

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif // dealii_differentiation_sd_h


