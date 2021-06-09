//include/deal.II-translator/base/function_lib_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2021 by the deal.II authors
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

#ifndef dealii_function_lib_h
#define dealii_function_lib_h


#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/table.h>

#include <array>

DEAL_II_NAMESPACE_OPEN

/**
 * 实现一些从Function类派生的具体类的命名空间，这些类描述实际的函数。这是我们自己的程序曾经需要的类的集合，并且认为在某些时候可能对其他人也是有用的。
 *
 *
 * @ingroup functions
 *
 *
 */
namespace Functions
{
  /**
   * 到原点的距离的平方。
   * 这个函数返回一个点的半径向量的平方准则。
   * 连同该函数，它的导数和拉普拉斯也被定义。
   * @ingroup functions
   *
   */
  template <int dim>
  class SquareFunction : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;
    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override;
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;
    virtual void
    vector_gradient(const Point<dim> &           p,
                    std::vector<Tensor<1, dim>> &gradient) const override;
    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component = 0) const override;
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;
    virtual void
    laplacian_list(const std::vector<Point<dim>> &points,
                   std::vector<double> &          values,
                   const unsigned int             component = 0) const override;
  };



  /**
   * 2d和3d中的函数<tt>xy</tt>，在1d中没有实现。这个函数可以作为一个消失的拉普拉斯的例子。
   * @ingroup functions
   *
   */
  template <int dim>
  class Q1WedgeFunction : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component = 0) const override;

    virtual void
    vector_gradient_list(
      const std::vector<Point<dim>> &,
      std::vector<std::vector<Tensor<1, dim>>> &) const override;

    /**
     * 该函数在一个点上的拉普拉斯系数。
     *
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    /**
     * 多点上的拉普拉斯函数。
     *
     */
    virtual void
    laplacian_list(const std::vector<Point<dim>> &points,
                   std::vector<double> &          values,
                   const unsigned int             component = 0) const override;
  };



  /**
   * 单位超立方体上的d-quadratic枕头。
   * 这是一个用于测试实现的函数。它在域 $(-1,1)^d$
   * 上具有零迪里切特边界值。在内部，它是 $1-x_i^2$
   * 在所有空间维度上的积。
   * 向构造函数提供一个非零参数，整个函数可以被一个常数抵消。
   * 与该函数一起，它的导数和拉普拉斯也被定义。
   * @ingroup functions
   *
   */
  template <int dim>
  class PillowFunction : public Function<dim>
  {
  public:
    /**
     * 构造函数。提供一个常数，它将被添加到每个函数值中。
     *
     */
    PillowFunction(const double offset = 0.);

    /**
     * 在一个点上的值。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 在多个点上的值。
     *
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * 单点上的梯度。
     *
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * 多点上的梯度。
     *
     */
    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component = 0) const override;

    /**
     * 单点的拉普拉斯系数。
     *
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    /**
     * 多点上的拉普拉斯。
     *
     */
    virtual void
    laplacian_list(const std::vector<Point<dim>> &points,
                   std::vector<double> &          values,
                   const unsigned int             component = 0) const override;

  private:
    const double offset;
  };



  /**
   * 余弦形的枕头函数。这是另一个在 $[-1,1]^d$
   * 上具有零边界值的函数。在内部，它是 $\cos(\pi/2 x_i)$
   * 的积。
   * @ingroup functions
   *
   */
  template <int dim>
  class CosineFunction : public Function<dim>
  {
  public:
    /**
     * 构造函数允许选择性地生成一个矢量值的余弦函数，每个分量的值都相同。
     *
     */
    CosineFunction(const unsigned int n_components = 1);

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component = 0) const override;

    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    virtual void
    laplacian_list(const std::vector<Point<dim>> &points,
                   std::vector<double> &          values,
                   const unsigned int             component = 0) const override;

    /**
     * 单点的二阶导数。
     *
     */
    virtual SymmetricTensor<2, dim>
    hessian(const Point<dim> & p,
            const unsigned int component = 0) const override;

    /**
     * 多点上的二阶导数。
     *
     */
    virtual void
    hessian_list(const std::vector<Point<dim>> &       points,
                 std::vector<SymmetricTensor<2, dim>> &hessians,
                 const unsigned int component = 0) const override;
  };



  /**
   * 余弦形枕头函数的梯度。    这是一个具有 @p dim
   * 分量的向量值函数，是余弦函数的梯度。在正方形[-1,1]上，它的切向边界条件为零。因此，它可以用来测试麦克斯韦算子的实现，而不必理会边界项。
   * @ingroup functions
   *
   */
  template <int dim>
  class CosineGradFunction : public Function<dim>
  {
  public:
    /**
     * 构造函数，创建一个具有 @p dim 成分的函数。
     *
     */
    CosineGradFunction();

    virtual double
    value(const Point<dim> &p, const unsigned int component) const override;
    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override;
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component) const override;

    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim> &p, const unsigned int component) const override;

    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component) const override;

    virtual void
    vector_gradient_list(
      const std::vector<Point<dim>> &           points,
      std::vector<std::vector<Tensor<1, dim>>> &gradients) const override;

    virtual double
    laplacian(const Point<dim> &p, const unsigned int component) const override;
  };



  /**
   * 每个坐标方向上的指数函数的乘积。
   * @ingroup functions
   *
   */
  template <int dim>
  class ExpFunction : public Function<dim>
  {
  public:
    /**
     * 在一个点上的值。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 在多点的数值。
     *
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * 单点上的梯度。
     *
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * 多点上的梯度。
     *
     */
    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component = 0) const override;

    /**
     * 单点的拉普拉斯系数。
     *
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    /**
     * 多点上的拉普拉契亚。
     *
     */
    virtual void
    laplacian_list(const std::vector<Point<dim>> &points,
                   std::vector<double> &          values,
                   const unsigned int             component = 0) const override;
  };



  /**
   * 一个能解决拉普拉斯方程的函数（有特定的边界值，但右手边为零），该函数在二维的L形域的中心有一个奇点（即在这个非凸域的重心角的位置）。
   * 该函数在极坐标中由 $r^{\frac{2}{3}}
   * \sin(\frac{2}{3} \phi)$ 给出，其奇点在原点，应与 GridGenerator::hyper_L(). 一起使用。这里， $\phi$ 被定义为对正 $x$ 轴的顺时针*角。    这个函数经常被用来说明拉普拉斯方程@f[
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   *
   * -\Delta u = 0
   * @f]的解可以是奇异的，即使边界值是平滑的。（这里，如果域是L形域 $(-1,1)^2 \backslash [0,1]^2$ ， $u$ 的边界值在邻近原点的两条线段上是零，而在边界的其余部分等于 $r^{\frac{2}{3}} \sin(\frac{2}{3} \phi)$ ）。该函数本身在域上仍然是有界的，但它的梯度在原点附近是 $r^{-1/3}$ 的形式，因此在接近原点时发散了。
   * @ingroup functions
   *
   */
  class LSingularityFunction : public Function<2>
  {
  public:
    virtual double
    value(const Point<2> &p, const unsigned int component = 0) const override;

    virtual void
    value_list(const std::vector<Point<2>> &points,
               std::vector<double> &        values,
               const unsigned int           component = 0) const override;

    virtual void
    vector_value_list(const std::vector<Point<2>> &points,
                      std::vector<Vector<double>> &values) const override;

    virtual Tensor<1, 2>
    gradient(const Point<2> &   p,
             const unsigned int component = 0) const override;

    virtual void
    gradient_list(const std::vector<Point<2>> &points,
                  std::vector<Tensor<1, 2>> &  gradients,
                  const unsigned int           component = 0) const override;

    virtual void
    vector_gradient_list(
      const std::vector<Point<2>> &,
      std::vector<std::vector<Tensor<1, 2>>> &) const override;

    virtual double
    laplacian(const Point<2> &   p,
              const unsigned int component = 0) const override;

    virtual void
    laplacian_list(const std::vector<Point<2>> &points,
                   std::vector<double> &        values,
                   const unsigned int           component = 0) const override;
  };



  /**
   * 二维L型域上的谐波奇点的梯度。
   * LSingularityFunction的梯度，它是一个卷曲和发散都消失的矢量值函数。
   * @ingroup functions
   *
   */
  class LSingularityGradFunction : public Function<2>
  {
  public:
    /**
     * 默认构造函数将维度设置为2。
     *
     */
    LSingularityGradFunction();
    virtual double
    value(const Point<2> &p, const unsigned int component) const override;

    virtual void
    value_list(const std::vector<Point<2>> &points,
               std::vector<double> &        values,
               const unsigned int           component) const override;

    virtual void
    vector_value_list(const std::vector<Point<2>> &points,
                      std::vector<Vector<double>> &values) const override;

    virtual Tensor<1, 2>
    gradient(const Point<2> &p, const unsigned int component) const override;

    virtual void
    gradient_list(const std::vector<Point<2>> &points,
                  std::vector<Tensor<1, 2>> &  gradients,
                  const unsigned int           component) const override;

    virtual void
    vector_gradient_list(
      const std::vector<Point<2>> &,
      std::vector<std::vector<Tensor<1, 2>>> &) const override;

    virtual double
    laplacian(const Point<2> &p, const unsigned int component) const override;

    virtual void
    laplacian_list(const std::vector<Point<2>> &points,
                   std::vector<double> &        values,
                   const unsigned int           component) const override;
  };



  /**
   * 二维和三维狭缝域上的奇异性。
   * @ingroup functions
   *
   */
  template <int dim>
  class SlitSingularityFunction : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component = 0) const override;

    virtual void
    vector_gradient_list(
      const std::vector<Point<dim>> &,
      std::vector<std::vector<Tensor<1, dim>>> &) const override;

    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    virtual void
    laplacian_list(const std::vector<Point<dim>> &points,
                   std::vector<double> &          values,
                   const unsigned int             component = 0) const override;
  };


  /**
   * 二维中具有一个诺伊曼边界的狭缝域上的奇异性。
   * @ingroup functions
   *
   */
  class SlitHyperSingularityFunction : public Function<2>
  {
  public:
    virtual double
    value(const Point<2> &p, const unsigned int component = 0) const override;

    virtual void
    value_list(const std::vector<Point<2>> &points,
               std::vector<double> &        values,
               const unsigned int           component = 0) const override;

    virtual void
    vector_value_list(const std::vector<Point<2>> &points,
                      std::vector<Vector<double>> &values) const override;

    virtual Tensor<1, 2>
    gradient(const Point<2> &   p,
             const unsigned int component = 0) const override;

    virtual void
    gradient_list(const std::vector<Point<2>> &points,
                  std::vector<Tensor<1, 2>> &  gradients,
                  const unsigned int           component = 0) const override;

    virtual void
    vector_gradient_list(
      const std::vector<Point<2>> &,
      std::vector<std::vector<Tensor<1, 2>>> &) const override;

    virtual double
    laplacian(const Point<2> &   p,
              const unsigned int component = 0) const override;

    virtual void
    laplacian_list(const std::vector<Point<2>> &points,
                   std::vector<double> &        values,
                   const unsigned int           component = 0) const override;
  };



  /**
   * 在x方向上的跳跃被输送到某个方向。
   * 如果平移平行于y轴，函数为<tt>-atan(sx)</tt>，其中<tt>s</tt>是构造函数中提供的陡度参数。
   * 对于不同的平流方向，这个函数将在参数空间中转动。
   * 与该函数一起，其导数和拉普拉斯也被定义。
   * @ingroup functions
   *
   */
  template <int dim>
  class JumpFunction : public Function<dim>
  {
  public:
    /**
     * 构造函数。在这里提供平流方向和斜率的陡度。
     *
     */
    JumpFunction(const Point<dim> &direction, const double steepness);

    /**
     * 在一个点上的函数值。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 多点上的函数值。
     *
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * 在一个点上的梯度。
     *
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * 多点上的梯度。
     *
     */
    virtual void
    gradient_list(const std::vector<Point<dim>> &points,
                  std::vector<Tensor<1, dim>> &  gradients,
                  const unsigned int             component = 0) const override;

    /**
     * 一个点上的函数的拉普拉斯系数。
     *
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    /**
     * 多点上的拉普拉斯函数。
     *
     */
    virtual void
    laplacian_list(const std::vector<Point<dim>> &points,
                   std::vector<double> &          values,
                   const unsigned int             component = 0) const override;

    /**
     * 返回这个对象的内存消耗估计值，以字节为单位。这不是精确的（但通常会很接近），因为计算树的内存使用量（例如，
     * <tt>std::map</tt>) 是很困难的。
     *
     */
    virtual std::size_t
    memory_consumption() const override;

  protected:
    /**
     * 对流矢量。
     *
     */
    const Point<dim> direction;

    /**
     * 斜坡的陡度（最大导数）。
     *
     */
    const double steepness;

    /**
     * 平流角。
     *
     */
    double angle;

    /**
     * <tt>角的正弦</tt>。
     *
     */
    double sine;

    /**
     * <tt>角的余弦</tt>。
     *
     */
    double cosine;
  };



  /**
   * 给出一个文数向量，生成一个余弦函数。文数系数以傅里叶空间中的
   * $d$  -维点 $k$ 的形式给出，然后函数被恢复为 $f(x) =
   * \cos(\sum_i k_i x_i) = Re(\exp(i k.x))$  。
   * 该类的名称来自于它类似于傅里叶余弦分解的一个组成部分。
   * @ingroup functions
   *
   */
  template <int dim>
  class FourierCosineFunction : public Function<dim>
  {
  public:
    /**
     * 构造函数。以每个空间方向上的傅里叶系数作为参数。
     *
     */
    FourierCosineFunction(const Tensor<1, dim> &fourier_coefficients);

    /**
     * 返回该函数在给定点的值。除非只有一个分量（即函数是标量的），否则你应该说明你想要评估的分量；它默认为零，即第一个分量。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 返回函数的指定分量在给定点的梯度。
     *
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * 计算给定分量在<tt>p</tt>点的拉普拉斯。
     *
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

  private:
    /**
     * 存储的傅里叶系数。
     *
     */
    const Tensor<1, dim> fourier_coefficients;
  };



  /**
   * 给定一个文数向量生成一个正弦函数。文数系数以傅里叶空间中的
   * $d$  -维点 $k$ 的形式给出，然后函数被恢复为 $f(x) =
   * \sin(\sum_i k_i x_i) = Im(\exp(i k.x))$  。
   * 该类的名称来自于它类似于傅里叶正弦分解的一个组成部分。
   * @ingroup functions
   *
   */
  template <int dim>
  class FourierSineFunction : public Function<dim>
  {
  public:
    /**
     * 构造函数。以每个空间方向上的傅里叶系数作为参数。
     *
     */
    FourierSineFunction(const Tensor<1, dim> &fourier_coefficients);

    /**
     * 返回该函数在给定点的值。除非只有一个分量（即函数是标量的），否则你应该说明你想要评估的分量；它默认为零，即第一个分量。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 返回函数的指定分量在给定点的梯度。
     *
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * 计算给定分量在<tt>p</tt>点的拉普拉斯。
     *
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

  private:
    /**
     * 存储的傅里叶系数。
     *
     */
    const Tensor<1, dim> fourier_coefficients;
  };


  /**
   * 给出一串文数向量和权重，生成正弦函数之和。每个文数系数在傅里叶空间中作为
   * $d$ -维点 $k$ 给出，然后整个函数被恢复为 $f(x) = \sum_j
   * w_j sin(\sum_i k_i x_i) = Im(\sum_j w_j \exp(i k.x))$  。
   * @ingroup functions
   *
   */
  template <int dim>
  class FourierSineSum : public Function<dim>
  {
  public:
    /**
     * 构造器。以每个空间方向上的傅里叶系数为参数。
     *
     */
    FourierSineSum(const std::vector<Point<dim>> &fourier_coefficients,
                   const std::vector<double> &    weights);

    /**
     * 返回该函数在给定点的值。除非只有一个分量（即函数是标量的），否则你应该说明你想要评估的分量；它默认为零，即第一个分量。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 返回函数的指定分量在给定点的梯度。
     *
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * 计算给定分量在<tt>p</tt>点的拉普拉斯。
     *
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

  private:
    /**
     * 存储的傅里叶系数和权重。
     *
     */
    const std::vector<Point<dim>> fourier_coefficients;
    const std::vector<double>     weights;
  };



  /**
   * 给出一串文数向量和权重，生成余弦函数之和。每个文数系数在傅里叶空间中以
   * $d$  -维点 $k$ 的形式给出，然后整个函数恢复为 $f(x) =
   * \sum_j w_j cos(\sum_i k_i x_i) = Re(\sum_j w_j \exp(i k.x))$  。
   * @ingroup functions
   *
   */
  template <int dim>
  class FourierCosineSum : public Function<dim>
  {
  public:
    /**
     * 构造器。以每个空间方向上的傅里叶系数为参数。
     *
     */
    FourierCosineSum(const std::vector<Point<dim>> &fourier_coefficients,
                     const std::vector<double> &    weights);

    /**
     * 返回该函数在给定点的值。除非只有一个分量（即函数是标量的），否则你应该说明你想要评估的分量；它默认为零，即第一个分量。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 返回函数的指定分量在给定点的梯度。
     *
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * 计算给定分量在<tt>p</tt>点的拉普拉斯。
     *
     */
    virtual double
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

  private:
    /**
     * 存储的傅里叶系数和权重。
     *
     */
    const std::vector<Point<dim>> fourier_coefficients;
    const std::vector<double>     weights;
  };


  /**
   * 切断函数的基础函数。该类存储了一个截断函数的中心和支持球的半径。如果该函数是矢量值的，它还存储了非零分量的数量。
   * 这个类也可以用于近似的狄拉克三角函数。这些是特殊的截断函数，其积分总是等于1，与支持球的半径无关。
   * @ingroup functions
   *
   */
  template <int dim>
  class CutOffFunctionBase : public Function<dim>
  {
  public:
    /**
     * 在这个和派生类的构造函数中使用的值，表示没有选择任何组件。
     *
     */
    static const unsigned int no_component = numbers::invalid_unsigned_int;

    /**
     * 构造函数。          @param[in]  radius 球的半径  @param[in]
     * center 球的中心  @param[in]  n_components
     * 这个函数对象的组件数量  @param[in]  select 如果这与
     * CutOffFunctionBase<dim>::no_component,
     * 不同，那么函数将只对这个组件非零  @param[in]  ]
     * integrate_to_one
     * 每当设置一个新的半径时，都要重新调整函数的值，以保证积分等于1
     * @param[in]  unitary_integral_value
     * 当半径等于1.0时的积分值。派生类将需要提供这个值，以保证正确地执行重新缩放。
     *
     */
    CutOffFunctionBase(
      const double       radius       = 1.,
      const Point<dim>   center       = Point<dim>(),
      const unsigned int n_components = 1,
      const unsigned int select       = CutOffFunctionBase<dim>::no_component,
      const bool         integrate_to_one       = false,
      const double       unitary_integral_value = 1.0);

    /**
     * 虚拟解构器。
     *
     */
    virtual ~CutOffFunctionBase() = default;

    /**
     * 将球的中心设置为 @p p. 点。
     *
     */
    virtual void
    set_center(const Point<dim> &p);

    /**
     * 设置球的半径为 @p r 。
     *
     */
    virtual void
    set_radius(const double r);

    /**
     * 返回存储在这个对象中的中心。
     *
     */
    const Point<dim> &
    get_center() const;

    /**
     * 返回存储在此对象中的半径。
     *
     */
    double
    get_radius() const;

    /**
     * 返回一个布尔值，表示此函数是否积分到一。
     *
     */
    bool
    integrates_to_one() const;

  protected:
    /**
     * 积分球的中心。
     *
     */
    Point<dim> center;

    /**
     * 该球的半径。
     *
     */
    double radius;

    /**
     * 选择的组件。如果<tt>no_component</tt>，函数在所有组件中是相同的。
     *
     */
    const unsigned int selected;

    /**
     * 控制半径变化时我们是否重新调整数值的标志。
     *
     */
    bool integrate_to_one;

    /**
     * 参考积分值。当 @p radius
     * =1.0时，派生类应指定其积分是什么。
     *
     */
    const double unitary_integral_value;

    /**
     * 当前的重定比例，以应用截止函数。
     *
     */
    double rescaling;
  };


  /**
   * CutOffFunctionBase对象的张量积。
   * 该类不是使用距离来计算截断函数，而是在每个坐标方向上对同一CutOffFunctionBase对象进行张量乘积。
   * @ingroup functions
   *
   */
  template <int dim>
  class CutOffFunctionTensorProduct : public CutOffFunctionBase<dim>
  {
  public:
    /**
     * 构建一个空的CutOffFunctionTensorProduct对象。
     * 在你使用这个类之前，你必须用一个从CutOffFunctionBase对象派生的类来调用set_base()方法。
     * 如果你在调用set_base()方法之前试图使用这个类，将会触发异常。
     *
     */
    CutOffFunctionTensorProduct(
      double             radius       = 1.0,
      const Point<dim> & center       = Point<dim>(),
      const unsigned int n_components = 1,
      const unsigned int select       = CutOffFunctionBase<dim>::no_component,
      const bool         integrate_to_one = false);

    /**
     * 用一个类型为 @tparam
     * CutOffFunctionBaseType<1>的对象初始化这个类。
     *
     */
    template <template <int> class CutOffFunctionBaseType>
    void
    set_base();

    /**
     * 设置新的中心。
     *
     */
    virtual void
    set_center(const Point<dim> &center) override;

    /**
     * 设置新的半径。
     *
     */
    virtual void
    set_radius(const double radius) override;

    /**
     * 在一个点上的函数值。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 在一个点上的函数梯度。
     *
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

  private:
    std::array<std::unique_ptr<CutOffFunctionBase<1>>, dim> base;

    bool initialized;
  };



  /**
   * 任意球的L-无穷大中的截止函数。
   * 这个函数是围绕<tt>中心</tt>的球的特征函数，具有指定的<tt>半径</tt>，也就是说，\f[
   * f = \chi(B_r(c)). \f]
   * 如果是矢量值，它可以被限制为一个分量。
   * @ingroup functions
   *
   */
  template <int dim>
  class CutOffFunctionLinfty : public CutOffFunctionBase<dim>
  {
  public:
    /**
     * 构造函数。参数是球的中心和它的半径。
     * 如果一个参数<tt>select</tt>被给定，并且不是
     *
     * - ，则截断函数将只对这个分量非零。
     *
     */
    CutOffFunctionLinfty(
      const double radius             = 1.,
      const Point<dim>                = Point<dim>(),
      const unsigned int n_components = 1,
      const unsigned int select       = CutOffFunctionBase<dim>::no_component,
      const bool         integrate_to_one = false);

    /**
     * 在一个点上的函数值。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 多点上的函数值。
     *
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * 多点上的函数值。
     *
     */
    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;
  };


  /**
   * 任意球的截断函数。这个函数是一个圆锥体，支持在一个围绕<tt>中心</tt>的特定</tt>球中。最大值是1。如果是矢量值，它可以被限制为一个分量。
   * @ingroup functions
   *
   */
  template <int dim>
  class CutOffFunctionW1 : public CutOffFunctionBase<dim>
  {
  public:
    /**
     * 构造函数。参数是球的中心和它的半径。
     * 如果给了一个参数<tt>select</tt>，截断函数将只对这个分量非零。
     *
     */
    CutOffFunctionW1(
      const double radius             = 1.,
      const Point<dim>                = Point<dim>(),
      const unsigned int n_components = 1,
      const unsigned int select       = CutOffFunctionBase<dim>::no_component,
      const bool         integrate_to_one = false);

    /**
     * 在一个点上的函数值。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 多点上的函数值。
     *
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * 多点上的函数值。
     *
     */
    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;
  };


  /**
   * 一个任意大小的球的截断函数，它在空间 $C^1$
   * （即连续可微）。这是一个在沉浸边界法的文献中经常使用的截止函数。
   * 该函数在径向坐标中的表达式为 $f(r)=1/2(cos(\pi r/s)+1)$
   * ，其中 $r<s$ 是到中心的距离， $s$
   * 是球体的半径。如果是矢量值，它可以被限制为一个单一分量。
   * @ingroup functions
   *
   */
  template <int dim>
  class CutOffFunctionC1 : public CutOffFunctionBase<dim>
  {
  public:
    /**
     * 构造函数。
     *
     */
    CutOffFunctionC1(
      const double radius             = 1.,
      const Point<dim>                = Point<dim>(),
      const unsigned int n_components = 1,
      const unsigned int select       = CutOffFunctionBase<dim>::no_component,
      bool               integrate_to_one = false);

    /**
     * 在一个点上的函数值。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 多点上的函数值。
     *
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * 多点上的函数值。
     *
     */
    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;

    /**
     * 在一个点上的函数梯度。
     *
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;
  };


  /**
   * 任意球的截断函数。这是传统的C-无穷大的截止函数，用于围绕<tt>中心</tt>的一定半径的球，
   * $f(r)=exp(1-1/(1-r**2/s**2))$  ，其中 $r$ 是到中心的距离， $s$
   * 是球体的半径。如果是矢量值，它可以被限制为一个单一成分。
   * @ingroup functions
   *
   */
  template <int dim>
  class CutOffFunctionCinfty : public CutOffFunctionBase<dim>
  {
  public:
    /**
     * 构造函数。参数是球的中心和它的半径。
     * 如果给了一个参数<tt>select</tt>，截断函数将只对这个分量非零。
     *
     */
    CutOffFunctionCinfty(
      const double radius             = 1.,
      const Point<dim>                = Point<dim>(),
      const unsigned int n_components = 1,
      const unsigned int select       = CutOffFunctionBase<dim>::no_component,
      bool               integrate_to_one = false);

    /**
     * 在一个点上的函数值。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 多点上的函数值。
     *
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * 多点上的函数值。
     *
     */
    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<double>> &  values) const override;

    /**
     * 在一个点上的函数梯度。
     *
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;
  };



  /**
   * 一个代表单项式的函数对象的类。单项式是只有一个单项的多项式，即在1-d中它们的形式为
   * $x^\alpha$ ，在2-d中的形式为 $x_1^{\alpha_1}x_2^{\alpha_2}$
   * ，而在3-d中的形式为 $x_1^{\alpha_1}x_2^{\alpha_2}x_3^{\alpha_3}$
   * 。因此，单项式是由指数的 $dim$
   * 元组来描述的。因此，该类的构造函数需要一个Tensor<1,dim>来描述指数的集合。大多数时候，这些指数当然是整数，但实数指数当然也同样有效。当基数是负数时，指数不可能是实数。
   * @ingroup functions
   *
   */
  template <int dim, typename Number = double>
  class Monomial : public Function<dim, Number>
  {
  public:
    /**
     * 构造函数。第一个参数在该类的一般描述中解释。第二个参数表示这个对象应代表的向量分量的数量。所有的向量分量将有相同的值。
     *
     */
    Monomial(const Tensor<1, dim, Number> &exponents,
             const unsigned int            n_components = 1);

    /**
     * 在一个点上的函数值。
     *
     */
    virtual Number
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 返回一个矢量值函数在某一点的所有分量。
     * <tt>值</tt>应事先有正确的大小，即#n_components。
     *
     */
    virtual void
    vector_value(const Point<dim> &p, Vector<Number> &values) const override;

    /**
     * 在多个点上的函数值。
     *
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<Number> &          values,
               const unsigned int             component = 0) const override;

    /**
     * 在一个点上的函数梯度。
     *
     */
    virtual Tensor<1, dim, Number>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

  private:
    /**
     * 指数的集合。
     *
     */
    const Tensor<1, dim, Number> exponents;
  };



  /**
   * 一个标量函数，通过（双，三）线性插值从一组点数据中计算其值，这些数据排列在一个可能是不均匀的张量积网中。换句话说，考虑到三维的情况，假设有点
   * $x_0,\ldots, x_{K-1}$ ， $y_0,\ldots,y_{L-1}$ ， $z_1,\ldots,z_{M-1}$
   * ，以及定义在点 $(x_k,y_l,z_m)^T$ 的数据 $d_{klm}$
   * ，那么在点 $\mathbf x=(x,y,z)$
   * 上评估函数将找到盒子，以便 $x_k\le x\le x_{k+1}, y_l\le y\le
   * y_{l+1}, z_m\le z\le z_{m+1}$
   * ，并对该单元的数据做三线性内插。类似的操作在更低的维度上也会进行。
   * 这一类最常被用于评估在域内一些点上通过实验提供的系数或右手边，或者用于比较有限元网格上的解的输出与先前在网格上定义的数据。
   * @note  如果点 $x_i$ 实际上是在一个区间 $[x_0,x_1]$
   * 上等距分布的，并且在更高的维度上其他数据点也是如此，你应该使用InterpolatedUniformGridData类来代替。
   * 如果一个点被要求在坐标阵列的端点所定义的框外，那么该函数被假定为在每个坐标方向的最后一个数据点之外简单地扩展常量值。
   * (如果一个点位于框外，该类不会抛出一个错误，因为经常发生的情况是，一个点正好位于框外，其数量与数字四舍五入的数量相同。)
   * @note  相关类InterpolatedUniformGridData的使用在  step-53
   * 中讨论。      <h3>Dealing with large data sets</h3>
   * 这个类经常被用来插值由相当大的数据表提供的数据，这些数据表从磁盘上读取的成本很高，而且在并行（MPI）程序的每个进程上复制时都会占用大量的内存。
   * 表类可以通过使用共享内存只在必要时存储数据来帮助摊薄这一成本
   *
   * - 参见TableBase类的文档。一旦我们获得了这样一个表对象，它使用共享内存来存储数据，并且只在必要的时候使用，我们就必须避免当前的类复制*。
   * 表到它自己的成员变量中。相反，有必要使用这个类的move*构造函数来接管表和其共享内存空间的所有权。这可以通过TableBase类的文档中显示的代码片段的以下扩展来实现。
   * @code
   *  const unsigned int N=..., M=...;     // table sizes, assumed known
   *  Table<2,double>    data_table;
   *  const unsigned int root_rank = 0;
   *
   *  if (Utilities::MPI::this_mpi_process(mpi_communicator) == root_rank)
   *  {
   *    data_table.resize (N,M);
   *
   *    std::ifstream input_file ("data_file.dat");
   *    ...;                               // read the data from the file
   *  }
   *
   *  // Now distribute to all processes
   *  data_table.replicate_across_communicator (mpi_communicator, root_rank);
   *
   *  // Set up the x- and y-coordinates of the points stored in the
   *  // data table
   *  std::array<std::vector<double>, dim> coordinate_values;
   *  ...;                                 // do what needs to be done
   *
   *  // And finally set up the interpolation object. The calls
   *  // to std::move() make sure that the tables are moved into
   *  // the memory space of the InterpolateTensorProductGridData
   *  // object:
   *  InterpolatedTensorProductGridData<2>
   *        interpolation_function (std::move(coordinate_values),
   *                                std::move(data_table));
   * @endcode
   * @ingroup functions
   *
   */
  template <int dim>
  class InterpolatedTensorProductGridData : public Function<dim>
  {
  public:
    /**
     * 构造函数，用 @p
     * data_values中给出的数据初始化这个类实例。
     * @param  coordinate_values
     * 一个由dim数组组成的数组。每个内部数组包含坐标值
     * $x_0,\ldots, x_{K-1}$
     * ，类似地，其他坐标方向也是如此。这些数组不需要有相同的大小。很明显，对于一个dim-维的函数对象，我们需要dim这样的数组。这个数组内的坐标值被假定为严格的升序，以便于有效的查找。
     * @param  data_values
     * 一个由上述坐标数组定义的每个网格点的二维数据表格。传入的数据会被复制到内部数据结构中。请注意，Table类有一些转换构造函数，可以将其他数据类型转换为你指定这个参数的表格。
     *
     */
    InterpolatedTensorProductGridData(
      const std::array<std::vector<double>, dim> &coordinate_values,
      const Table<dim, double> &                  data_values);

    /**
     * 像前面的构造函数一样，但是把参数作为r值引用，并移动*，而不是复制*数据。在这些表中存储的数据很大，并且不再需要单独初始化当前对象的信息的情况下，这通常是有用的。换句话说，没有必要保留原始对象，这个对象可以从中复制其信息，但它不妨接管（"移动"）数据。
     *
     */
    InterpolatedTensorProductGridData(
      std::array<std::vector<double>, dim> &&coordinate_values,
      Table<dim, double> &&                  data_values);

    /**
     * 通过对给定的数据集进行双线性插值，计算出函数集的值。
     * @param  p 要对函数进行评估的点。      @param  component
     * 矢量分量。因为这个函数是标量的，所以只有0是这里的有效参数。
     * @return
     * 该点的内插值。如果该点位于坐标集之外，则该函数将以一个常数扩展。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 计算由给定数据集的双线性内插定义的函数的梯度。
     * @param  p 要对函数梯度进行评估的点。      @param
     * component
     * 矢量分量。因为这个函数是标量的，所以只有0是这里的有效参数。
     * @return
     * 内插函数的梯度值，在此点。如果该点位于坐标集之外，则该函数被扩展为常数，因此其梯度被扩展为0。
     *
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * 返回这个对象的内存消耗估计值，单位是字节。
     *
     */
    virtual std::size_t
    memory_consumption() const override;

    /**
     * 返回一个对内部存储数据的引用。
     *
     */
    const Table<dim, double> &
    get_data() const;

  protected:
    /**
     * 查找包含输入点的矩形表格中的索引
     *
     */
    TableIndices<dim>
    table_index_of_point(const Point<dim> &p) const;

    /**
     * 每个坐标方向上的坐标值集合。
     *
     */
    const std::array<std::vector<double>, dim> coordinate_values;

    /**
     * 要内插的数据。
     *
     */
    const Table<dim, double> data_values;
  };


  /**
   * 一个标量函数，通过（双、三）线性内插计算其数值，这些数据是在一个均匀间隔的张量积网格上排列的一组点数据。换句话说，考虑到三维的情况，让点
   * $x_0,\ldots, x_{K-1}$ 从区间 $[x_0,x_{K-1}]$ 统一细分为大小为
   * $\Delta x = (x_{K-1}-x_0)/(K-1)$ 的子区间，以及类似的
   * $y_0,\ldots,y_{L-1}$ ， $z_1,\ldots,z_{M-1}$ 。同样考虑数据
   * $d_{klm}$  定义在点  $(x_k,y_l,z_m)^T$  上，然后在点  $\mathbf
   * x=(x,y,z)$  上评估函数将找到盒子，从而  $x_k\le x\le
   * x_{k+1}, y_l\le y\le y_{l+1}, z_m\le z\le z_{m+1}$
   * ，并对该单元的数据做三线插值。类似的操作在更低的维度上也能完成。
   * 该类最常被用于评估在域内一些点上通过实验提供的系数或右手边，或者用于比较有限元网格上的解的输出与先前在网格上定义的数据。
   * @note  如果你有一个问题，其中的点 $x_i$
   * 不是等距的（例如，它们是在一个分级网格上计算的结果，该网格更接近一个边界的密度），那么使用InterpolatedTensorProductGridData类来代替。
   * 如果一个点被要求在坐标阵列的端点所定义的范围之外，那么该函数被假定为简单地通过常量值扩展到每个坐标方向的最后一个数据点之外。
   * (如果一个点位于框外，该类不会抛出一个错误，因为经常发生的情况是，一个点正好位于框外，其数量与数字四舍五入的数量相同。)
   * @note  这个类的使用在  step-53  中讨论。      <h3>Dealing with
   * large data sets</h3>
   * 该类支持与InterpolatedTensorProductGridData类一样的处理大型数据集的设施。更多信息和示例代码见那里。
   * @ingroup functions
   *
   */
  template <int dim>
  class InterpolatedUniformGridData : public Function<dim>
  {
  public:
    /**
     * 构造函数  @param  interval_endpoints
     * 每个坐标方向上的（均匀细分的）区间的左右端点。
     * @param  n_subintervals
     * 每个坐标方向上的子间隔的数量。一个坐标的值为1，意味着该区间被认为是由整个范围组成的一个子区间。值为2意味着有两个子区间，每个子区间有一半的范围，等等。
     * @param  data_values
     * 一个由上述坐标阵列定义的每个网格点的数据的二维表。请注意，Table类有一些转换构造函数，可以将其他数据类型转换为你指定这个参数的表格。
     *
     */
    InterpolatedUniformGridData(
      const std::array<std::pair<double, double>, dim> &interval_endpoints,
      const std::array<unsigned int, dim> &             n_subintervals,
      const Table<dim, double> &                        data_values);

    /**
     * 像前面的构造函数一样，但是把参数作为r值引用，并移动*，而不是复制*数据。在这些表中存储的数据很大，并且不再需要单独初始化当前对象的信息的情况下，这通常是有用的。换句话说，没有必要保留原始对象，这个对象可以从中复制其信息，但它不妨接管（"移动"）数据。
     *
     */
    InterpolatedUniformGridData(
      std::array<std::pair<double, double>, dim> &&interval_endpoints,
      std::array<unsigned int, dim> &&             n_subintervals,
      Table<dim, double> &&                        data_values);

    /**
     * 通过对给定的数据集进行双线性插值，计算出函数集的值。
     * @param  p 要对函数进行评估的点。      @param  component
     * 矢量分量。因为这个函数是标量的，所以只有零是这里的有效参数。
     * @return
     * 该点的内插值。如果该点位于坐标集之外，则该函数将以一个常数扩展。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 通过对给定的数据集进行双线性内插，计算函数集的梯度。
     * @param  p 要对函数进行评估的点。      @param  component
     * 矢量分量。因为这个函数是标量的，所以只有0是这里的有效参数。
     * @return
     * 内插函数在该点的梯度。如果该点位于坐标集之外，则该函数被扩展为一个常数，其梯度当然为零。
     *
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * 返回这个对象的内存消耗估计值，单位是字节。
     *
     */
    virtual std::size_t
    memory_consumption() const override;

    /**
     * 返回一个对内部存储数据的引用。
     *
     */
    const Table<dim, double> &
    get_data() const;

  private:
    /**
     * 每个坐标方向上的区间端点的集合。
     *
     */
    const std::array<std::pair<double, double>, dim> interval_endpoints;

    /**
     * 每个坐标方向上的子区间数量。
     *
     */
    const std::array<unsigned int, dim> n_subintervals;

    /**
     * 要进行内插的数据。
     *
     */
    const Table<dim, double> data_values;
  };


  /**
   * 一个代表多项式的函数对象的类。多项式是由多个单项式的相加组成的。如果多项式有n个单项式并且维度等于dim，那么多项式可以写成
   * $\sum_{i=1}^{n} a_{i}(\prod_{d=1}^{dim} x_{d}^{\alpha_{i,d}})$ ，其中
   * $a_{i}$ 是单项式的系数， $\alpha_{i,d}$
   * 是其指数。该类的构造函数需要一个Table<2,double>来描述指数集，一个Vector<double>来描述系数集。
   * @ingroup functions
   *
   */
  template <int dim>
  class Polynomial : public Function<dim>
  {
  public:
    /**
     * 构造函数。多项式的系数和指数被作为参数传递。表<2,
     * double>的指数的行数等于多项式的单项式的数目，列数等于dim。指数表的第i行包含第i个单项式的
     * ${\alpha_{i,d}}$ 指数 $a_{i}\prod_{d=1}^{dim} x_{d}^{\alpha_{i,d}}$
     * 。系数向量的第i个元素包含第i个单项式的系数 $a_{i}$
     * 。
     *
     */
    Polynomial(const Table<2, double> &   exponents,
               const std::vector<double> &coefficients);

    /**
     * 在一个点上的函数值。
     *
     */
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override;


    /**
     * 多点的函数值。
     *
     */
    virtual void
    value_list(const std::vector<Point<dim>> &points,
               std::vector<double> &          values,
               const unsigned int             component = 0) const override;

    /**
     * 在一个点上的函数梯度。
     *
     */
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * 返回这个对象的内存消耗估计值，以字节为单位。
     *
     */
    virtual std::size_t
    memory_consumption() const override;

  private:
    /**
     * 指数的集合。
     *
     */
    const Table<2, double> exponents;

    /**
     * 系数的集合。
     *
     */
    const std::vector<double> coefficients;
  };

#ifndef DOXYGEN



  // Template definitions
  template <int dim>
  template <template <int> class CutOffFunctionBaseType>
  void
  CutOffFunctionTensorProduct<dim>::set_base()
  {
    initialized = true;
    static_assert(
      std::is_base_of<CutOffFunctionBase<1>, CutOffFunctionBaseType<1>>::value,
      "You can only construct a CutOffFunctionTensorProduct function from "
      "a class derived from CutOffFunctionBase.");

    for (unsigned int i = 0; i < dim; ++i)
      base[i].reset(new CutOffFunctionBaseType<1>(this->radius,
                                                  Point<1>(this->center[i]),
                                                  this->n_components,
                                                  this->selected,
                                                  this->integrate_to_one));
  }



#endif

} // namespace Functions
DEAL_II_NAMESPACE_CLOSE

#endif


