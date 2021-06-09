//include/deal.II-translator/numerics/data_component_interpretation_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2020 by the deal.II authors
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

#ifndef dealii_data_component_interpretation_h
#define dealii_data_component_interpretation_h



#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个专门用于声明
 * DataComponentInterpretation::DataComponentInterpretation
 * 枚举的命名空间。
 *
 *
 */
namespace DataComponentInterpretation
{
  /**
   * 这个枚举的成员是用来描述矢量值数据集的各个组成部分的逻辑解释的。例如，如果一个人有一个2d的斯托克斯方程的有限元，代表组件
   * $(u,v,p)$ ，我们希望表明前两个， $u$ 和 $v$
   * ，代表一个逻辑矢量，这样以后当我们生成图形输出时，我们可以把它们交给一个可视化程序，它将自动知道把它们作为一个矢量场来渲染，而不是作为两个独立的标量场。
   * 通过向 DataOut_DoFData::add_data_vector
   * 函数传递一组当前种类的枚举，这可以实现。    参见
   * step-22
   * 教程程序中关于如何在实践中使用这一信息的例子。
   *
   */
  enum DataComponentInterpretation
  {
    /**
     * 表示一个数据集的某个分量对应于一个独立于其他分量的标量场。
     *
     */
    component_is_scalar,

    /**
     * 表示一个数据集的某个分量是一个矢量值量的一部分。
     *
     */
    component_is_part_of_vector,

    /**
     * 表示一个数据集的某个分量是张量值（二阶）数量的一部分。
     *
     */
    component_is_part_of_tensor
  };
} // namespace DataComponentInterpretation


DEAL_II_NAMESPACE_CLOSE

#endif


