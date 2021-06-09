//include/deal.II-translator/base/function_cspline_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

#ifndef dealii_function_cspline_h
#define dealii_function_cspline_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_GSL
#  include <deal.II/base/function.h>
#  include <deal.II/base/point.h>
#  include <deal.II/base/thread_management.h>

#  include <gsl/gsl_spline.h>

#  include <memory>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  DeclException1(ExcCSplineEmpty,
                 int,
                 << "Interpolation points vector size can not be <" << arg1
                 << ">.");

  DeclException2(ExcCSplineSizeMismatch,
                 int,
                 int,
                 << "The size of interpolation points <" << arg1
                 << "> is different from the size of interpolation values <"
                 << arg2 << ">.");


  DeclException3(ExcCSplineOrder,
                 int,
                 double,
                 double,
                 << "The input interpolation points are not strictly ordered : "
                 << std::endl
                 << "x[" << arg1 << "] = " << arg2 << " >= x[" << (arg1 + 1)
                 << "] = " << arg3 << ".");

  DeclException3(
    ExcCSplineRange,
    double,
    double,
    double,
    << "Spline function can not be evaluated outside of the interpolation range: "
    << std::endl
    << arg1 << " is not in [" << arg2 << ";" << arg3 << "].");

  /**
   * 使用GNU科学库的立方样条曲线函数。
   * 产生的曲线在每个区间上都是片状立体的，在提供的数据点上有匹配的第一和第二导数。第二导数在第一点和最后一点被选择为零。
   * @note  这个函数只在dim==1时实现。
   * @ingroup functions
   *
   */
  template <int dim>
  class CSpline : public Function<dim>
  {
  public:
    /**
     * 构造函数，应提供一组要进行插值的点  @p
     * interpolation_points  和一组函数值  @p interpolation_values  。
     *
     */
    CSpline(const std::vector<double> &interpolation_points,
            const std::vector<double> &interpolation_values);

    virtual double
    value(const Point<dim> & point,
          const unsigned int component = 0) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    virtual SymmetricTensor<2, dim>
    hessian(const Point<dim> & p,
            const unsigned int component = 0) const override;

    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    /**
     * 返回这个对象的内存消耗估计值，以字节为单位。
     *
     */
    virtual std::size_t
    memory_consumption() const override;

  private:
    /**
     * 提供插值的点
     *
     */
    const std::vector<double> interpolation_points;

    /**
     * 函数在插值点的值
     *
     */
    const std::vector<double> interpolation_values;

    /**
     * 用于花键插值的GSL加速器
     *
     */
    std::unique_ptr<gsl_interp_accel, void (*)(gsl_interp_accel *)> acc;

    /**
     * GSL三次样条插值器
     *
     */
    std::unique_ptr<gsl_spline, void (*)(gsl_spline *)> cspline;

    /**
     * 一个用于加速器对象的突变器。
     *
     */
    mutable Threads::Mutex acc_mutex;
  };
} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif

#endif


