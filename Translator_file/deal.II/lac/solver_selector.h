//include/deal.II-translator/lac/solver_selector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#ifndef dealii_solver_selector_h
#define dealii_solver_selector_h


#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/vector.h>

DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup Solvers */ 
 /*@{*/ 

/**
 * 通过改变一个参数来选择一个求解器。 通过调用这个 @p
 * SolverSelector, 的 @p solve
 * 函数，它选择在这个类的构造函数中指定的那个 @p Solver
 * 的 @p solve 函数。 <h3>Usage</h3>
 * 这个类的最简单的使用方法是如下。
 *
 * @code
 * // generate a @p SolverSelector that calls the @p SolverCG
 * SolverControl control;
 * SolverSelector<Vector<double> > solver_selector ("cg", control);
 *
 * // generate e.g. a @p PreconditionRelaxation
 * PreconditionRelaxation<SparseMatrix<double>, Vector<double> >
 * preconditioning(
 *   A, &SparseMatrix<double>::template precondition_SSOR<double>,0.8);
 *
 * // call the @p solve function with this preconditioning as last argument
 * solver_selector.solve(A,x,b,preconditioning);
 * @endcode
 * 但 @p SolverSelector
 * 类的全部用处要到下面的例子介绍时才清楚，该例子假设用户使用
 * @p ParameterHandler 类，并且已经声明了一个 "解算器
 * "条目，例如，用
 *
 * @code
 * Parameter_Handler prm;
 * prm.declare_entry ("solver", "none",
 *                  Patterns::Selection(
 *                    SolverSelector<>::get_solver_names()));
 * ...
 * @endcode
 * 假设在用户的参数文件中，存在这样一行
 *
 * @code
 * set solver = cg
 * @endcode
 * 那么上面的例子中的构造器调用可以写成
 *
 * @code
 * SolverSelector<SparseMatrix<double>, Vector<double> >
 * solver_selector(prm.get("solver"), control);
 * @endcode
 *
 *  如果在某个时候存在一个新的求解器
 * "xyz"，那么用户不需要改变他们的程序。只有在 @p
 * SolverSelector
 * 的实现中，才需要增加对这个求解器的调用，每个拥有上述程序行的用户只需要在他们的参数文件中'set
 * solver = xyz'就可以获得对这个新求解器的访问。
 *
 *
 */
template <typename VectorType = Vector<double>>
class SolverSelector : public Subscriptor
{
public:
  /**
   * 底层矢量类型的别名
   *
   */
  using vector_type = VectorType;

  /**
   * 构造函数，填入默认值
   *
   */
  SolverSelector() = default;

  /**
   * 构造函数，选择解算器 @p name 和解算器控制对象 @p
   * control 了。
   *
   */
  SolverSelector(const std::string &name, SolverControl &control);

  /**
   * 破坏器
   *
   */
  virtual ~SolverSelector() override = default;

  /**
   * 解算器程序。调用 @p solver 的 @p solve 函数，其 @p
   * SolverName是在构造函数中指定的。
   *
   */
  template <class Matrix, class Preconditioner>
  void
  solve(const Matrix &        A,
        VectorType &          x,
        const VectorType &    b,
        const Preconditioner &precond) const;

  /**
   * 选择一个新的求解器。请注意，这个类中使用的所有求解器名称都是小写的。
   *
   */
  void
  select(const std::string &name);

  /**
   * 设置一个新的SolverControl。这需要在解算前设置。
   *
   */
  void
  set_control(SolverControl &ctrl);

  /**
   * 设置附加数据。更多信息见 @p Solver 类。
   *
   */
  void
  set_data(const typename SolverRichardson<VectorType>::AdditionalData &data);

  /**
   * 设置附加数据。更多信息见 @p Solver 类。
   *
   */
  void
  set_data(const typename SolverCG<VectorType>::AdditionalData &data);

  /**
   * 设置附加数据。更多信息见 @p Solver 类。
   *
   */
  void
  set_data(const typename SolverMinRes<VectorType>::AdditionalData &data);

  /**
   * 设置附加数据。更多信息见 @p Solver 类。
   *
   */
  void
  set_data(const typename SolverBicgstab<VectorType>::AdditionalData &data);

  /**
   * 设置附加数据。更多信息见 @p Solver 类。
   *
   */
  void
  set_data(const typename SolverGMRES<VectorType>::AdditionalData &data);

  /**
   * 设置附加数据。更多信息见 @p Solver 类。
   *
   */
  void
  set_data(const typename SolverFGMRES<VectorType>::AdditionalData &data);

  /**
   * 获取所有实现的求解器的名称。可能的选项列表包括。    <ul>   <li>  "Richardson"  </li>   <li>  "cg"  </li>   <li>  "bicgstab"  </li>   <li>  "gmres"  </li>   <li>  "fgres"  </li>   <li>  "minres"  </li>   </ul>
   *
   */
  static std::string
  get_solver_names();

  /**
   * 异常情况。
   *
   */
  DeclException1(ExcSolverDoesNotExist,
                 std::string,
                 << "Solver " << arg1 << " does not exist. Use one of "
                 << std::endl
                 << get_solver_names());



protected:
  /**
   * 存储每个 @p 解算器类的构造函数中需要的 @p SolverControl
   * 。这可以用 @p set_control(). 来改变。
   *
   */
  SmartPointer<SolverControl, SolverSelector<VectorType>> control;

  /**
   * 存储解算器的名称。
   *
   */
  std::string solver_name;

private:
  /**
   * 存储额外的数据。
   *
   */
  typename SolverRichardson<VectorType>::AdditionalData richardson_data;

  /**
   * 存储额外的数据。
   *
   */
  typename SolverCG<VectorType>::AdditionalData cg_data;

  /**
   * 存储额外的数据。
   *
   */
  typename SolverMinRes<VectorType>::AdditionalData minres_data;

  /**
   * 存储额外的数据。
   *
   */
  typename SolverBicgstab<VectorType>::AdditionalData bicgstab_data;

  /**
   * 存储额外的数据。
   *
   */
  typename SolverGMRES<VectorType>::AdditionalData gmres_data;

  /**
   * 存储额外的数据。
   *
   */
  typename SolverFGMRES<VectorType>::AdditionalData fgmres_data;
};

 /*@}*/ 
 /* --------------------- Inline and template functions ------------------- */ 


template <typename VectorType>
SolverSelector<VectorType>::SolverSelector(const std::string &name,
                                           SolverControl &    solver_control)
  : solver_name(name)
  , control(&solver_control)
{}



template <typename VectorType>
void
SolverSelector<VectorType>::select(const std::string &name)
{
  solver_name = name;
}



template <typename VectorType>
template <class Matrix, class Preconditioner>
void
SolverSelector<VectorType>::solve(const Matrix &        A,
                                  VectorType &          x,
                                  const VectorType &    b,
                                  const Preconditioner &precond) const
{
  if (solver_name == "richardson")
    {
      SolverRichardson<VectorType> solver(*control, richardson_data);
      solver.solve(A, x, b, precond);
    }
  else if (solver_name == "cg")
    {
      SolverCG<VectorType> solver(*control, cg_data);
      solver.solve(A, x, b, precond);
    }
  else if (solver_name == "minres")
    {
      SolverMinRes<VectorType> solver(*control, minres_data);
      solver.solve(A, x, b, precond);
    }
  else if (solver_name == "bicgstab")
    {
      SolverBicgstab<VectorType> solver(*control, bicgstab_data);
      solver.solve(A, x, b, precond);
    }
  else if (solver_name == "gmres")
    {
      SolverGMRES<VectorType> solver(*control, gmres_data);
      solver.solve(A, x, b, precond);
    }
  else if (solver_name == "fgmres")
    {
      SolverFGMRES<VectorType> solver(*control, fgmres_data);
      solver.solve(A, x, b, precond);
    }
  else
    Assert(false, ExcSolverDoesNotExist(solver_name));
}



template <typename VectorType>
void
SolverSelector<VectorType>::set_control(SolverControl &ctrl)
{
  control = &ctrl;
}



template <typename VectorType>
std::string
SolverSelector<VectorType>::get_solver_names()
{
  return "richardson|cg|bicgstab|gmres|fgmres|minres";
}



template <typename VectorType>
void
SolverSelector<VectorType>::set_data(
  const typename SolverGMRES<VectorType>::AdditionalData &data)
{
  gmres_data = data;
}



template <typename VectorType>
void
SolverSelector<VectorType>::set_data(
  const typename SolverFGMRES<VectorType>::AdditionalData &data)
{
  fgmres_data = data;
}



template <typename VectorType>
void
SolverSelector<VectorType>::set_data(
  const typename SolverRichardson<VectorType>::AdditionalData &data)
{
  richardson_data = data;
}



template <typename VectorType>
void
SolverSelector<VectorType>::set_data(
  const typename SolverCG<VectorType>::AdditionalData &data)
{
  cg_data = data;
}



template <typename VectorType>
void
SolverSelector<VectorType>::set_data(
  const typename SolverMinRes<VectorType>::AdditionalData &data)
{
  minres_data = data;
}



template <typename VectorType>
void
SolverSelector<VectorType>::set_data(
  const typename SolverBicgstab<VectorType>::AdditionalData &data)
{
  bicgstab_data = data;
}

DEAL_II_NAMESPACE_CLOSE

#endif


