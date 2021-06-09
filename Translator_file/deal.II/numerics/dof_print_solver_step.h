//include/deal.II-translator/numerics/dof_print_solver_step_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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

#ifndef dealii_dof_print_solver_step_h
#define dealii_dof_print_solver_step_h

#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iomanip>
#include <sstream>

DEAL_II_NAMESPACE_OPEN


/**
 * 打印求解器中的中间解。
 * 这是从作为模板参数提供的求解器类中派生出来的。
 * 它使用DoFHandler实现了求解器的 @p print_vector
 * 功能。这样一来，中间向量可以被看作是有限元函数。这个类首先可以用来理解求解器是如何工作的（例如，可视化各种求解器的平滑特性，例如在多网格背景下），其次可以研究为什么和如何求解器不能解决某些类别的问题。
 * 这个类的对象通过模板参数提供一个求解器类，并提供一个文件名（作为一个字符串），在每次迭代中用这个文件构建一个新的文件（命名为<tt>basename.[step].[senix]</tt>），并使用DataOut类将解决方案作为一个有限元场写入其中。请注意，这个类可能会产生巨大的数据量!
 *
 *
 * @ingroup output
 *
 *
 */
template <int dim, typename SolverType, class VectorType = Vector<double>>
class DoFPrintSolverStep : public SolverType
{
public:
  /**
   * 构造函数。 首先，我们接受求解器所需的参数。  @p
   * data_out是作为有限元函数做输出的对象。
   * 每个迭代步骤将产生一个名称为<tt>basename.[step].[senix]</tt>的输出文件。
   *
   */
  DoFPrintSolverStep(SolverControl &           control,
                     VectorMemory<VectorType> &mem,
                     DataOut<dim> &            data_out,
                     const std::string &       basename);

  /**
   * 迭代方法的回调函数。
   *
   */
  virtual void
  print_vectors(const unsigned int step,
                const VectorType & x,
                const VectorType & r,
                const VectorType & d) const;

private:
  /**
   * 输出对象。
   *
   */
  DataOut<dim> &out;

  /**
   * 文件名的基础。
   *
   */
  const std::string basename;
};


 /* ----------------------- template functions --------------- */ 

template <int dim, typename SolverType, class VectorType>
DoFPrintSolverStep<dim, SolverType, VectorType>::DoFPrintSolverStep(
  SolverControl &           control,
  VectorMemory<VectorType> &mem,
  DataOut<dim> &            data_out,
  const std::string &       basename)
  : SolverType(control, mem)
  , out(data_out)
  , basename(basename)
{}


template <int dim, typename SolverType, class VectorType>
void
DoFPrintSolverStep<dim, SolverType, VectorType>::print_vectors(
  const unsigned int step,
  const VectorType & x,
  const VectorType & r,
  const VectorType & d) const
{
  out.clear_data_vectors();
  out.add_data_vector(x, "solution");
  out.add_data_vector(r, "residual");
  out.add_data_vector(d, "update");

  std::ostringstream filename;
  filename << basename << std::setw(3) << std::setfill('0') << step
           << out.default_suffix();

  const std::string fname = filename.str();

  deallog << "Writing file:" << fname << std::endl;

  out.build_patches();
  std::ofstream of(fname.c_str());
  out.write(of);
}

DEAL_II_NAMESPACE_CLOSE

#endif


