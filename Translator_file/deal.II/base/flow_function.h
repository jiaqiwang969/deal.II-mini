//include/deal.II-translator/base/flow_function_0.txt
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

#ifndef dealii_flow_function_h
#define dealii_flow_function_h


#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/thread_management.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{
  /**
   * 用于分析解决不可压缩的流动问题的基类。
   * 除了Function接口之外，这个函数还提供了压力的偏移：如果计算出的解决方案的压力有一个不同于零的积分平均值，这个值可以给pressure_adjustment()，以便计算出正确的压力误差。
   * @note  派生类应该始终实现压力的积分平均值为零。
   * @note
   * 线程安全。一些函数利用内部数据来计算数值。因此，每个线程都应该获得自己的派生类对象。
   * @ingroup functions
   *
   */
  template <int dim>
  class FlowFunction : public Function<dim>
  {
  public:
    /**
     * 构造函数，设置一些内部数据结构。
     *
     */
    FlowFunction();

    /**
     * 虚拟解构器。
     *
     */
    virtual ~FlowFunction() override = default;

    /**
     * 存储压力函数的调整值，使其平均值为<tt>p</tt>。
     *
     */
    void
    pressure_adjustment(double p);

    /**
     * 在一个更适合于矢量值函数的结构中的值。外层向量以解分量为索引，内层以正交点为索引。
     *
     */
    virtual void
    vector_values(const std::vector<Point<dim>> &   points,
                  std::vector<std::vector<double>> &values) const override = 0;
    /**
     * 在一个更适合于矢量值函数的结构中的梯度。外侧矢量以解分量为索引，内部以正交点为索引。
     *
     */
    virtual void
    vector_gradients(
      const std::vector<Point<dim>> &           points,
      std::vector<std::vector<Tensor<1, dim>>> &gradients) const override = 0;
    /**
     * 更适合于矢量值函数的结构中的力项。
     * 外侧向量以解分量为索引，内部以正交点为索引。
     * @warning
     * 这不是真正的拉普拉斯，而是斯托克斯方程中作为右手边的力项。
     *
     */
    virtual void
    vector_laplacians(const std::vector<Point<dim>> &   points,
                      std::vector<std::vector<double>> &values) const = 0;

    virtual void
    vector_value(const Point<dim> &points,
                 Vector<double> &  value) const override;
    virtual double
    value(const Point<dim> & points,
          const unsigned int component) const override;
    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;
    virtual void
    vector_gradient_list(
      const std::vector<Point<dim>> &           points,
      std::vector<std::vector<Tensor<1, dim>>> &gradients) const override;
    /**
     * 动量方程中的力项。
     *
     */
    virtual void
    vector_laplacian_list(const std::vector<Point<dim>> &points,
                          std::vector<Vector<double>> &  values) const override;

    /**
     * 返回这个对象的内存消耗估计值，以字节为单位。
     *
     */
    virtual std::size_t
    memory_consumption() const override;

  protected:
    /**
     * 派生类要增加的压力的平均值。
     *
     */
    double mean_pressure;

  private:
    /**
     * 一个突扰器，可以守护以下的抓取数组。
     *
     */
    mutable Threads::Mutex mutex;

    /**
     * 通常Function接口的辅助值。
     *
     */
    mutable std::vector<std::vector<double>> aux_values;

    /**
     * 通常Function接口的辅助值。
     *
     */
    mutable std::vector<std::vector<Tensor<1, dim>>> aux_gradients;
  };

  /**
   * 二维和三维的层流管道流动。通道沿<i>x</i>轴延伸，半径为
   * @p radius. 。 @p Reynolds
   * 数字用于适当缩放纳维-斯托克斯问题的压力。
   * @ingroup functions
   *
   */
  template <int dim>
  class PoisseuilleFlow : public FlowFunction<dim>
  {
  public:
    /**
     * 为给定的通道半径<tt>r</tt>和雷诺数<tt>Re</tt>构建一个对象。
     *
     */
    PoisseuilleFlow<dim>(const double r, const double Re);

    virtual ~PoisseuilleFlow() override = default;

    virtual void
    vector_values(const std::vector<Point<dim>> &   points,
                  std::vector<std::vector<double>> &values) const override;
    virtual void
    vector_gradients(
      const std::vector<Point<dim>> &           points,
      std::vector<std::vector<Tensor<1, dim>>> &gradients) const override;
    virtual void
    vector_laplacians(const std::vector<Point<dim>> &   points,
                      std::vector<std::vector<double>> &values) const override;

  private:
    const double radius;
    const double Reynolds;
  };


  /**
   * 在立方体[-1,1]<sup>dim</sup>上具有同质边界条件的人工无发散函数。    该函数在二维中是@f[
   * \left(\begin{array}{c}u\\v\\p\end{array}\right)
   * \left(\begin{array}{c}\cos^2x \sin y\cos y\\-\sin x\cos x\cos^2y\\
   * \sin x\cos x\sin y\cos y\end{array}\right)
   * @f]。
   * @ingroup functions
   *
   */
  template <int dim>
  class StokesCosine : public FlowFunction<dim>
  {
  public:
    /**
     * 构造函数设置压力计算所需的雷诺数和右手边的缩放比例。
     *
     */
    StokesCosine(const double viscosity = 1., const double reaction = 0.);
    /**
     * 改变粘度和反应参数。
     *
     */
    void
    set_parameters(const double viscosity, const double reaction);

    virtual ~StokesCosine() override = default;

    virtual void
    vector_values(const std::vector<Point<dim>> &   points,
                  std::vector<std::vector<double>> &values) const override;
    virtual void
    vector_gradients(
      const std::vector<Point<dim>> &           points,
      std::vector<std::vector<Tensor<1, dim>>> &gradients) const override;
    virtual void
    vector_laplacians(const std::vector<Point<dim>> &   points,
                      std::vector<std::vector<double>> &values) const override;

  private:
    /// The viscosity
    double viscosity;
    /// The reaction parameter
    double reaction;
  };


  /**
   * 斯托克斯方程在2d L形域上的奇异解。    该函数满足
   * $-\triangle \mathbf{u} + \nabla p = 0$
   * ，代表L形域的重入角周围的典型奇异解，可以用
   * GridGenerator::hyper_L().
   * 创建，速度在重入角的两个面上消失， $\nabla\mathbf{u}$
   * ]和 $p$
   * 在原点是奇异的，而在域的其余部分是光滑的，因为它们可以写成一个光滑函数与项
   * $r^{\lambda-1}$ 的乘积，其中 $r$ 是半径， $\lambda \approx
   * 0.54448$ 是一个固定参数。    摘自Houston, Sch&ouml;tzau,
   * Wihler, 着手ENUMATH 2003。
   * @ingroup functions
   *
   */
  class StokesLSingularity : public FlowFunction<2>
  {
  public:
    /// Constructor setting up some data.
    StokesLSingularity();

    virtual void
    vector_values(const std::vector<Point<2>> &     points,
                  std::vector<std::vector<double>> &values) const override;
    virtual void
    vector_gradients(
      const std::vector<Point<2>> &           points,
      std::vector<std::vector<Tensor<1, 2>>> &gradients) const override;
    virtual void
    vector_laplacians(const std::vector<Point<2>> &     points,
                      std::vector<std::vector<double>> &values) const override;

  private:
    /// The auxiliary function Psi.
    double
    Psi(double phi) const;
    /// The derivative of Psi()
    double
    Psi_1(double phi) const;
    /// The 2nd derivative of Psi()
    double
    Psi_2(double phi) const;
    /// The 3rd derivative of Psi()
    double
    Psi_3(double phi) const;
    /// The 4th derivative of Psi()
    double
    Psi_4(double phi) const;
    /// The angle of the reentrant corner, set to 3*pi/2
    const double omega;
    /// The exponent of the radius, computed as the solution to
    /// $\sin(\lambda\omega)+\lambda \sin(\omega)=0$
    static const double lambda;
    /// Cosine of lambda times omega
    const double coslo;
    /// Auxiliary variable 1+lambda
    const double lp;
    /// Auxiliary variable 1-lambda
    const double lm;
  };

  /**
   * Kovasznay(1947)的二维流解。
   * 该函数在直线<i>x=1/2</i>右侧的半平面上有效。
   * @ingroup functions
   *
   */
  class Kovasznay : public FlowFunction<2>
  {
  public:
    /**
     * 构建一个给定雷诺数<tt>Re</tt>的对象。如果参数<tt>Stokes</tt>为真，由vector_laplacians()返回的动量方程的右侧包含非线性，这样就可以得到Kovasznay解作为Stokes问题的解。
     *
     */
    Kovasznay(const double Re, bool Stokes = false);

    virtual ~Kovasznay() override = default;

    virtual void
    vector_values(const std::vector<Point<2>> &     points,
                  std::vector<std::vector<double>> &values) const override;
    virtual void
    vector_gradients(
      const std::vector<Point<2>> &           points,
      std::vector<std::vector<Tensor<1, 2>>> &gradients) const override;
    virtual void
    vector_laplacians(const std::vector<Point<2>> &     points,
                      std::vector<std::vector<double>> &values) const override;

    /// The value of lambda.
    double
    lambda() const;

  private:
    const double Reynolds;
    double       lbda;
    double       p_average;
    const bool   stokes;
  };

} // namespace Functions

DEAL_II_NAMESPACE_CLOSE

#endif


