//include/deal.II-translator/differentiation/ad_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
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

#ifndef dealii_differentiation_ad_h
#define dealii_differentiation_ad_h

#include <deal.II/base/config.h>

#include <deal.II/differentiation/ad/ad_helpers.h>
#include <deal.II/differentiation/ad/ad_number_traits.h>
#include <deal.II/differentiation/ad/ad_number_types.h>
#include <deal.II/differentiation/ad/adolc_math.h>
#include <deal.II/differentiation/ad/adolc_number_types.h>
#include <deal.II/differentiation/ad/adolc_product_types.h>
#include <deal.II/differentiation/ad/sacado_math.h>
#include <deal.II/differentiation/ad/sacado_number_types.h>
#include <deal.II/differentiation/ad/sacado_product_types.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个命名空间，封装了与自动和符号微分有关的各种类和辅助函数。
 *
 *
 * @ingroup auto_symb_diff
 *
 */
namespace Differentiation
{
  /**
   * 自动分化库的封装器。目前支持以下库。
   *
   *
   *
   * - ADOL-C
   *
   *
   *
   *
   *
   * - 萨卡多（Trilinos的一个组成部分）。
   * @ingroup auto_symb_diff
   *
   */
  namespace AD
  {}
} // namespace Differentiation

DEAL_II_NAMESPACE_CLOSE

#endif


