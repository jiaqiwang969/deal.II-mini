//include/deal.II-translator/numerics/dof_output_operator_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2018 by the deal.II authors
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


#ifndef dealii_dof_output_operator_h
#define dealii_dof_output_operator_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/any_data.h>
#include <deal.II/algorithms/operator.h>

#include <deal.II/base/event.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  /**
   * 一个输出运算器，在每一步中写一个单独的文件，并将向量写成相对于给定DoFHandler的有限元函数。
   *
   */
  template <typename VectorType, int dim, int spacedim = dim>
  class DoFOutputOperator : public OutputOperator<VectorType>
  {
  public:
    /*构造器。<tt>filename</tt>是所有文件的共同基础名称，参数<tt>digits</tt>应该是序列中最高数字的位数。文件名默认为 "outputNNN "的形式，其中NNN是由上一步命令设定的数字。位数较少的数字从左边开始用零填充。   
*
*/
    DoFOutputOperator(const std::string &filename_base = std::string("output"),
                      const unsigned int digits        = 3);

    void
    parse_parameters(ParameterHandler &param);
    void
    initialize(const DoFHandler<dim, spacedim> &dof_handler);

    virtual OutputOperator<VectorType> &
    operator<<(const AnyData &vectors) override;

  private:
    SmartPointer<const DoFHandler<dim, spacedim>,
                 DoFOutputOperator<VectorType, dim, spacedim>>
      dof;

    const std::string  filename_base;
    const unsigned int digits;

    DataOut<dim> out;
  };

  template <typename VectorType, int dim, int spacedim>
  inline void
  DoFOutputOperator<VectorType, dim, spacedim>::initialize(
    const DoFHandler<dim, spacedim> &dof_handler)
  {
    dof = &dof_handler;
  }
} // namespace Algorithms


DEAL_II_NAMESPACE_CLOSE

#endif


