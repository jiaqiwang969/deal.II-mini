CCTest_file/step-67.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2020 - 2021 by the deal.II authors 
 * 
 * This file is part of the deal.II library. 
 * 
 * The deal.II library is free software; you can use it, redistribute 
 * it, and/or modify it under the terms of the GNU Lesser General 
 * Public License as published by the Free Software Foundation; either 
 * version 2.1 of the License, or (at your option) any later version. 
 * The full text of the license can be found in the file LICENSE.md at 
 * the top level directory of deal.II. 
 * 
 * --------------------------------------------------------------------- 
 * 
 * Author: Martin Kronbichler, 2020 
 */ 



// 包含文件与之前的无矩阵教程程序 step-37 、 step-48 和 step-59 相似。
#include <deal.II/base/conditional_ostream.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/base/time_stepping.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/base/vectorization.h> 

#include <deal.II/distributed/tria.h> 

#include <deal.II/dofs/dof_handler.h> 

#include <deal.II/fe/fe_dgq.h> 
#include <deal.II/fe/fe_system.h> 

#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/tria.h> 

#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/la_parallel_vector.h> 

#include <deal.II/matrix_free/fe_evaluation.h> 
#include <deal.II/matrix_free/matrix_free.h> 

#include <deal.II/numerics/data_out.h> 

#include <fstream> 
#include <iomanip> 
#include <iostream> 

// 下面的文件包括CellwiseInverseMassMatrix数据结构，我们将在质量矩阵反演中使用它，这是本教程程序中唯一的新包含文件。

#include <deal.II/matrix_free/operators.h> 

namespace Euler_DG 
{ 
  using namespace dealii; 

// 与其他无矩阵教程程序类似，我们在文件的顶部收集所有控制程序执行的参数。除了我们想要运行的维度和多项式程度，我们还指定了我们想要用于欧拉方程中非线性项的高斯正交公式的点数。此外，我们指定了随时间变化的问题的时间间隔，并实现了两个不同的测试案例。第一个是二维的分析解，而第二个是介绍中描述的围绕圆柱体的通道流。根据测试案例，我们还改变了运行模拟的最终时间，以及一个变量`output_tick`，它指定了我们要在哪个时间间隔内写入输出（假设tick大于时间步长）。

  constexpr unsigned int testcase             = 0; 
  constexpr unsigned int dimension            = 2; 
  constexpr unsigned int n_global_refinements = 3; 
  constexpr unsigned int fe_degree            = 5; 
  constexpr unsigned int n_q_points_1d        = fe_degree + 2; 

  using Number = double; 

  constexpr double gamma       = 1.4; 
  constexpr double final_time  = testcase == 0 ? 10 : 2.0; 
  constexpr double output_tick = testcase == 0 ? 1 : 0.05; 

// 接下来是时间积分器的一些细节，即用公式 $\Delta t =
//  \text{Cr} n_\text{stages} \frac{h}{(p+1)^{1.5} (\|\mathbf{u} +
//  c)_\text{max}}$ 来衡量时间步长的库朗数，以及选择一些低存储量的Runge--Kutta方法。我们指定Runge--Kutta方案每级的Courant数，因为这对不同级数的方案给出了一个更实际的数值成本表达。

  const double courant_number = 0.15 / std::pow(fe_degree, 1.5); 
  enum LowStorageRungeKuttaScheme 
  { 
    stage_3_order_3, /* Kennedy, Carpenter, Lewis, 2000 */ 


    stage_5_order_4, /* Kennedy, Carpenter, Lewis, 2000 */ 


    stage_7_order_4, /* Tselios, Simos, 2007 */ 


    stage_9_order_5, /* Kennedy, Carpenter, Lewis, 2000 */ 


  }; 
  constexpr LowStorageRungeKuttaScheme lsrk_scheme = stage_5_order_4; 

// 最终，我们选择了空间离散化的一个细节，即单元间面的数值通量（黎曼求解器）。在这个程序中，我们实现了Lax--Friedrichs通量和Harten--Lax--van Leer(HLL)通量的一个改进版本。

  enum EulerNumericalFlux 
  { 
    lax_friedrichs_modified, 
    harten_lax_vanleer, 
  }; 
  constexpr EulerNumericalFlux numerical_flux_type = lax_friedrichs_modified; 

//  @sect3{Equation data}  

// 我们现在定义了一个带有测试情况0的精确解的类和一个带有测试情况1的通道背景流场的类。鉴于欧拉方程是一个在 $d$ 维度上有 $d+2$ 个方程的问题，我们需要告诉函数基类正确的分量数量。

  template <int dim> 
  class ExactSolution : public Function<dim> 
  { 
  public: 
    ExactSolution(const double time) 
      : Function<dim>(dim + 2, time) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 
  }; 

// 就实际实现的函数而言，分析性测试案例是一个等熵涡旋案例（例如参见Hesthaven和Warburton的书，第209页第6.6节中的例6.1），它满足欧拉方程，右侧的力项为零。考虑到这个定义，我们返回密度、动量或能量，这取决于所要求的成分。请注意，密度的原始定义涉及一些表达式的 $\frac{1}{\gamma -1}$ -次方。由于 `std::pow()` 在某些系统上的实现相当慢，我们用对数和指数（以2为底）来代替它，这在数学上是等价的，但通常优化得更好。与 `std::pow()`, 相比，对于非常小的数字，这个公式可能会在最后一位数字上失去准确性，但我们还是很高兴，因为小数字映射为接近1的数据。

// 对于通道测试案例，我们简单地选择密度为1， $x$ 方向的速度为0.4，其他方向的速度为0，以及对应于背景速度场测量的1.3声速的能量，根据关系 $E = \frac{c^2}{\gamma (\gamma -1)} + \frac 12 \rho \|u\|^2$ 计算得出。

  template <int dim> 
  double ExactSolution<dim>::value(const Point<dim> & x, 
                                   const unsigned int component) const 
  { 
    const double t = this->get_time(); 

    switch (testcase) 
      { 
        case 0: 
          { 
            Assert(dim == 2, ExcNotImplemented()); 
            const double beta = 5; 

            Point<dim> x0; 
            x0[0] = 5.; 
            const double radius_sqr = 
              (x - x0).norm_square() - 2. * (x[0] - x0[0]) * t + t * t; 
            const double factor = 
              beta / (numbers::PI * 2) * std::exp(1. - radius_sqr); 
            const double density_log = std::log2( 
              std::abs(1. - (gamma - 1.) / gamma * 0.25 * factor * factor)); 
            const double density = std::exp2(density_log * (1. / (gamma - 1.))); 
            const double u       = 1. - factor * (x[1] - x0[1]); 
            const double v       = factor * (x[0] - t - x0[0]); 

            if (component == 0) 
              return density; 
            else if (component == 1) 
              return density * u; 
            else if (component == 2) 
              return density * v; 
            else 
              { 
                const double pressure = 
                  std::exp2(density_log * (gamma / (gamma - 1.))); 
                return pressure / (gamma - 1.) + 
                       0.5 * (density * u * u + density * v * v); 
              } 
          } 

        case 1: 
          { 
            if (component == 0) 
              return 1.; 
            else if (component == 1) 
              return 0.4; 
            else if (component == dim + 1) 
              return 3.097857142857143; 
            else 
              return 0.; 
          } 

        default: 
          Assert(false, ExcNotImplemented()); 
          return 0.; 
      } 
  } 

//  @sect3{Low-storage explicit Runge--Kutta time integrators}  

// 接下来的几行实现了一些低存储量的Runge--Kutta方法的变体。这些方法有特定的布彻表，系数为 $b_i$ 和 $a_i$ ，如介绍中所示。如同Runge--Kutta方法的惯例，我们可以从这些系数中推导出时间步骤 $c_i = \sum_{j=1}^{i-2} b_i + a_{i-1}$ 。这种方案的主要优点是每个阶段只需要两个向量，即解的累积部分 $\mathbf{w}$ （在最后一个阶段后的新时间 $t^{n+1}$ 保持解 $\mathbf{w}^{n+1}$ ），在各阶段被评估的更新向量 $\mathbf{r}_i$ ，加上一个向量 $\mathbf{k}_i$ 来保持算子评估。这样的Runge--Kutta设置减少了内存存储和内存访问。由于内存带宽通常是现代硬件上的性能限制因素，当微分算子的评估得到很好的优化时，性能可以比标准的时间积分器得到改善。考虑到传统的Runge--Kutta方案可能允许稍大的时间步长，因为更多的自由参数可以获得更好的稳定性，这一点也是真实的。

// 在本教程中，我们集中讨论Kennedy, Carpenter和Lewis(2000)文章中定义的低存储方案的几个变体，以及Tselios和Simos(2007)描述的一个变体。还有一大系列的其他方案，可以通过额外的系数集或稍微不同的更新公式来解决。

// 我们为这四种积分器定义了一个单一的类，用上述的枚举来区分。对每个方案，我们再将 $b_i$ 和 $a_i$ 的向量填充到类中的给定变量。

  class LowStorageRungeKuttaIntegrator 
  { 
  public: 
    LowStorageRungeKuttaIntegrator(const LowStorageRungeKuttaScheme scheme) 
    { 
      TimeStepping::runge_kutta_method lsrk; 

// 首先是Kennedy等人（2000）提出的三阶方案。虽然它的稳定区域比其他方案小得多，但它只涉及三个阶段，所以在每个阶段的工作方面很有竞争力。

      switch (scheme) 
        { 
          case stage_3_order_3: 
            { 
              lsrk = TimeStepping::LOW_STORAGE_RK_STAGE3_ORDER3; 
              break; 
            } 

// 下一个方案是四阶的五级方案，同样在Kennedy等人（2000）的论文中定义。

          case stage_5_order_4: 
            { 
              lsrk = TimeStepping::LOW_STORAGE_RK_STAGE5_ORDER4; 
              break; 
            } 

// 下面这个七级和四阶的方案已经明确地推导出用于声学问题。它在四阶方案中兼顾了虚特征值的精度，并结合了一个大的稳定区域。由于DG方案在最高频率之间是耗散的，这不一定转化为每级可能的最高时间步长。在本教程方案的背景下，数值通量在耗散中起着至关重要的作用，因此也是最大的稳定时间步长。对于修改后的Lax--Friedrichs通量，如果只考虑稳定性，该方案在每级步长方面与`stage_5_order_4`方案相似，但对于HLL通量来说，效率稍低。

          case stage_7_order_4: 
            { 
              lsrk = TimeStepping::LOW_STORAGE_RK_STAGE7_ORDER4; 
              break; 
            } 

// 这里包括的最后一个方案是Kennedy等人（2000）的五阶九级方案。它是这里使用的方案中最精确的，但是较高的精度牺牲了一些稳定性，所以每级的归一化步长比四阶方案要小。

          case stage_9_order_5: 
            { 
              lsrk = TimeStepping::LOW_STORAGE_RK_STAGE9_ORDER5; 
              break; 
            } 

          default: 
            AssertThrow(false, ExcNotImplemented()); 
        } 
      TimeStepping::LowStorageRungeKutta< 
        LinearAlgebra::distributed::Vector<Number>> 
        rk_integrator(lsrk); 
      rk_integrator.get_coefficients(ai, bi, ci); 
    } 

    unsigned int n_stages() const 
    { 
      return bi.size(); 
    } 

// 时间积分器的主要功能是通过阶段，评估算子，为下一次评估准备  $\mathbf{r}_i$  矢量，并更新解决方案矢量  $\mathbf{w}$  。我们把工作交给所涉及的`pde_operator`，以便能够把Runge--Kutta设置的矢量操作与微分算子的评估合并起来，以获得更好的性能，所以我们在这里所做的就是委托矢量和系数。

// 我们单独调用第一阶段的算子，因为我们需要稍微修改一下那里的参数。我们从旧的解决方案 $\mathbf{w}^n$ 而不是 $\mathbf r_i$ 向量中评估解决方案，所以第一个参数是`solution`。我们在这里让阶段向量 $\mathbf{r}_i$ 也持有评估的临时结果，因为它在其他情况下不会被使用。对于所有后续阶段，我们使用向量`vec_ki`作为第二个向量参数来存储运算符的求值结果。最后，当我们到了最后一个阶段，我们必须跳过对向量 $\mathbf{r}_{s+1}$ 的计算，因为没有系数 $a_s$ 可用（也不会用到）。

    template <typename VectorType, typename Operator> 
    void perform_time_step(const Operator &pde_operator, 
                           const double    current_time, 
                           const double    time_step, 
                           VectorType &    solution, 
                           VectorType &    vec_ri, 
                           VectorType &    vec_ki) const 
    { 
      AssertDimension(ai.size() + 1, bi.size()); 

      pde_operator.perform_stage(current_time, 
                                 bi[0] * time_step, 
                                 ai[0] * time_step, 
                                 solution, 
                                 vec_ri, 
                                 solution, 
                                 vec_ri); 

      for (unsigned int stage = 1; stage < bi.size(); ++stage) 
        { 
          const double c_i = ci[stage]; 
          pde_operator.perform_stage(current_time + c_i * time_step, 
                                     bi[stage] * time_step, 
                                     (stage == bi.size() - 1 ? 
                                        0 : 
                                        ai[stage] * time_step), 
                                     vec_ri, 
                                     vec_ki, 
                                     solution, 
                                     vec_ri); 
        } 
    } 

  private: 
    std::vector<double> bi; 
    std::vector<double> ai; 
    std::vector<double> ci; 
  }; 

//  @sect3{Implementation of point-wise operations of the Euler equations}  

// 在下面的函数中，我们实现了与欧拉方程有关的各种特定问题的运算。每个函数都作用于我们在解向量中持有的守恒变量向量 $[\rho, \rho\mathbf{u}, E]$ ，并计算各种派生量。

// 首先是速度的计算，我们从动量变量 $\rho \mathbf{u}$ 除以 $\rho$ 得出。这里需要注意的是，我们用关键字`DEAL_II_ALWAYS_INLINE`来装饰所有这些函数。这是一个特殊的宏，映射到一个编译器专用的关键字，告诉编译器永远不要为这些函数创建一个函数调用，而是将实现<a href="https:en.wikipedia.org/wiki/Inline_function">inline</a>移到它们被调用的地方。这对性能至关重要，因为我们对其中一些函数的调用达到了几百万甚至几十亿次。例如，我们既使用速度来计算通量，也使用速度来计算压力，而这两个地方都要在每个单元的每个正交点进行评估。确保这些函数是内联的，不仅可以确保处理器不必执行跳转指令进入函数（以及相应的返回跳转），而且编译器可以在调用函数的地方之后的代码中重新使用一个函数的上下文的中间信息。(我们注意到，编译器通常很善于自己找出哪些函数要内联。这里有一个地方，编译器可能是自己想出来的，也可能不是，但我们可以肯定的是，内联是一种胜利。)

// 我们应用的另一个技巧是为反密度设置一个单独的变量  $\frac{1}{\rho}$  。这使得编译器只对通量进行一次除法，尽管除法在多个地方使用。由于除法的费用大约是乘法或加法的10到20倍，避免多余的除法对性能至关重要。我们注意到，由于四舍五入的影响，在浮点运算中，先取反数，后与之相乘并不等同于除法，所以编译器不允许用标准的优化标志来交换一种方式。然而，以正确的方式编写代码也不是特别困难。

// 总而言之，所选择的总是内联和仔细定义昂贵的算术运算的策略使我们能够写出紧凑的代码，而不需要将所有的中间结果传递出去，尽管要确保代码映射到优秀的机器码。

  template <int dim, typename Number> 
  inline DEAL_II_ALWAYS_INLINE // 
    Tensor<1, dim, Number> 
    euler_velocity(const Tensor<1, dim + 2, Number> &conserved_variables) 
  { 
    const Number inverse_density = Number(1.) / conserved_variables[0]; 

    Tensor<1, dim, Number> velocity; 
    for (unsigned int d = 0; d < dim; ++d) 
      velocity[d] = conserved_variables[1 + d] * inverse_density; 

    return velocity; 
  } 

// 下一个函数从保守变量的矢量中计算压力，使用公式  $p = (\gamma - 1) \left(E - \frac 12 \rho \mathbf{u}\cdot \mathbf{u}\right)$  。如上所述，我们使用来自`euler_velocity()`函数的速度。注意，我们需要在这里指定第一个模板参数`dim`，因为编译器无法从张量的参数中推导出它，而第二个参数（数字类型）可以自动推导出来。

  template <int dim, typename Number> 
  inline DEAL_II_ALWAYS_INLINE // 
    Number 
    euler_pressure(const Tensor<1, dim + 2, Number> &conserved_variables) 
  { 
    const Tensor<1, dim, Number> velocity = 
      euler_velocity<dim>(conserved_variables); 

    Number rho_u_dot_u = conserved_variables[1] * velocity[0]; 
    for (unsigned int d = 1; d < dim; ++d) 
      rho_u_dot_u += conserved_variables[1 + d] * velocity[d]; 

    return (gamma - 1.) * (conserved_variables[dim + 1] - 0.5 * rho_u_dot_u); 
  } 

// 这里是欧拉通量函数的定义，也就是实际方程的定义。考虑到速度和压力（编译器的优化将确保只做一次），考虑到介绍中所说的方程，这是直截了当的。

  template <int dim, typename Number> 
  inline DEAL_II_ALWAYS_INLINE // 
    Tensor<1, dim + 2, Tensor<1, dim, Number>> 
    euler_flux(const Tensor<1, dim + 2, Number> &conserved_variables) 
  { 
    const Tensor<1, dim, Number> velocity = 
      euler_velocity<dim>(conserved_variables); 
    const Number pressure = euler_pressure<dim>(conserved_variables); 

    Tensor<1, dim + 2, Tensor<1, dim, Number>> flux; 
    for (unsigned int d = 0; d < dim; ++d) 
      { 
        flux[0][d] = conserved_variables[1 + d]; 
        for (unsigned int e = 0; e < dim; ++e) 
          flux[e + 1][d] = conserved_variables[e + 1] * velocity[d]; 
        flux[d + 1][d] += pressure; 
        flux[dim + 1][d] = 
          velocity[d] * (conserved_variables[dim + 1] + pressure); 
      } 

    return flux; 
  } 

// 接下来的这个函数是一个简化数值通量实现的助手，它实现了一个张量的张量（具有大小为`dim + 2`的非标准外维，所以deal.II的张量类提供的标准重载在此不适用）与另一个相同内维的张量的作用，即一个矩阵-向量积。

  template <int n_components, int dim, typename Number> 
  inline DEAL_II_ALWAYS_INLINE // 
    Tensor<1, n_components, Number> 
    operator*(const Tensor<1, n_components, Tensor<1, dim, Number>> &matrix, 
              const Tensor<1, dim, Number> &                         vector) 
  { 
    Tensor<1, n_components, Number> result; 
    for (unsigned int d = 0; d < n_components; ++d) 
      result[d] = matrix[d] * vector; 
    return result; 
  } 

// 这个函数实现了数值通量（黎曼求解器）。它从一个界面的两边获得状态，并获得法向量，从解的一边  $\mathbf{w}^-$  向解  $\mathbf{w}^+$  的方向。在依赖片断恒定数据的有限体积方法中，数值通量是核心成分，因为它是唯一输入物理信息的地方。在DG方法中，由于元素内部的多项式和那里使用的物理通量，数值通量就不那么核心了。由于在连续解的极限中，两边的数值一致的高阶插值，数值通量可以被看作是对两边解的跳跃的控制，以弱化连续性。必须认识到，在存在冲击的情况下，仅靠数值通量是无法稳定高阶DG方法的，因此任何DG方法都必须与进一步的冲击捕捉技术相结合，以处理这些情况。在本教程中，我们将重点讨论欧拉方程在没有强不连续的亚声速体系中的波状解，我们的基本方案已经足够了。

// 尽管如此，数值通量对整个方案的数值耗散起着决定性作用，并影响到显式Runge-Kutta方法的可接受的时间步长。我们考虑两种选择，一种是改良的Lax-Friedrichs方案，另一种是广泛使用的Harten-Lax-van Leer（HLL）通量。对于这两种方案，我们首先需要得到界面两边的速度和压力，并评估物理欧拉通量。

// 对于局部Lax--Friedrichs通量，其定义是 $\hat{\mathbf{F}}
//  =\frac{\mathbf{F}(\mathbf{w}^-)+\mathbf{F}(\mathbf{w}^+)}{2} +
//  \frac{\lambda}{2}\left[\mathbf{w}^--\mathbf{w}^+\right]\otimes
//  \mathbf{n^-}$  ，其中因子 $\lambda =
//  \max\left(\|\mathbf{u}^-\|+c^-, \|\mathbf{u}^+\|+c^+\right)$ 给出了最大波速， $c = \sqrt{\gamma p / \rho}$ 是音速。在这里，考虑到通量对解的影响很小，为了计算效率的原因，我们选择了该表达式的两个修改。对于上述因子 $\lambda$ 的定义，我们需要取四个平方根，两个用于两个速度规范，两个用于两侧的声速。因此，第一个修改是宁可使用 $\sqrt{\|\mathbf{u}\|^2+c^2}$ 作为最大速度的估计（如介绍中所示，它与实际最大速度最多相差2倍）。这使我们能够从最大速度中提取平方根，并且只需进行一次平方根计算就可以了。第二个修改是进一步放宽参数 $\lambda$ --它越小，耗散系数就越小（与 $\mathbf{w}$ 的跳跃相乘，最终可能导致耗散变小或变大）。这使得我们可以用更大的时间步长将频谱纳入显式Runge--Kutta积分器的稳定区域。然而，我们不能使耗散太小，因为否则假想的特征值会越来越大。最后，目前的保守公式在 $\lambda\to 0$ 的极限中不是能量稳定的，因为它不是偏斜对称的，在这种情况下需要额外的措施，如分裂形式的DG方案。

// 对于HLL通量，我们遵循文献中的公式，通过一个参数 $s$ 引入Lax--Friedrichs的两个状态的额外加权。它是由欧拉方程的物理传输方向得出的，以当前的速度方向和声速为准。对于速度，我们在此选择一个简单的算术平均数，这对危险情况和材料参数的适度跳跃是足够的。

// 由于数值通量在弱形式下是与法向量相乘的，因此我们对方程中的所有项都用法向量来乘以结果。在这些乘法中，上面定义的 "操作符*"可以实现类似于数学定义的紧凑符号。

// 在这个函数和下面的函数中，我们使用变量后缀`_m`和`_p`来表示从 $\mathbf{w}^-$ 和 $\mathbf{w}^+$ 得出的量，即在观察相邻单元时相对于当前单元的 "这里 "和 "那里 "的数值。

  template <int dim, typename Number> 
  inline DEAL_II_ALWAYS_INLINE // 
    Tensor<1, dim + 2, Number> 
    euler_numerical_flux(const Tensor<1, dim + 2, Number> &u_m, 
                         const Tensor<1, dim + 2, Number> &u_p, 
                         const Tensor<1, dim, Number> &    normal) 
  { 
    const auto velocity_m = euler_velocity<dim>(u_m); 
    const auto velocity_p = euler_velocity<dim>(u_p); 

    const auto pressure_m = euler_pressure<dim>(u_m); 
    const auto pressure_p = euler_pressure<dim>(u_p); 

    const auto flux_m = euler_flux<dim>(u_m); 
    const auto flux_p = euler_flux<dim>(u_p); 

    switch (numerical_flux_type) 
      { 
        case lax_friedrichs_modified: 
          { 
            const auto lambda = 
              0.5 * std::sqrt(std::max(velocity_p.norm_square() + 
                                         gamma * pressure_p * (1. / u_p[0]), 
                                       velocity_m.norm_square() + 
                                         gamma * pressure_m * (1. / u_m[0]))); 

            return 0.5 * (flux_m * normal + flux_p * normal) + 
                   0.5 * lambda * (u_m - u_p); 
          } 

        case harten_lax_vanleer: 
          { 
            const auto avg_velocity_normal = 
              0.5 * ((velocity_m + velocity_p) * normal); 
            const auto   avg_c = std::sqrt(std::abs( 
              0.5 * gamma * 
              (pressure_p * (1. / u_p[0]) + pressure_m * (1. / u_m[0])))); 
            const Number s_pos = 
              std::max(Number(), avg_velocity_normal + avg_c); 
            const Number s_neg = 
              std::min(Number(), avg_velocity_normal - avg_c); 
            const Number inverse_s = Number(1.) / (s_pos - s_neg); 

            return inverse_s * 
                   ((s_pos * (flux_m * normal) - s_neg * (flux_p * normal)) - 
                    s_pos * s_neg * (u_m - u_p)); 
          } 

        default: 
          { 
            Assert(false, ExcNotImplemented()); 
            return {}; 
          } 
      } 
  } 

// 这个函数和下一个函数是辅助函数，提供紧凑的评估调用，因为多个点通过VectorizedArray参数被分批放在一起（详见 step-37 教程）。这个函数用于亚音速外流边界条件，我们需要将能量分量设置为一个规定值。下一个函数请求所有分量上的解，用于流入边界，其中解的所有分量都被设置。

  template <int dim, typename Number> 
  VectorizedArray<Number> 
  evaluate_function(const Function<dim> &                      function, 
                    const Point<dim, VectorizedArray<Number>> &p_vectorized, 
                    const unsigned int                         component) 
  { 
    VectorizedArray<Number> result; 
    for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v) 
      { 
        Point<dim> p; 
        for (unsigned int d = 0; d < dim; ++d) 
          p[d] = p_vectorized[d][v]; 
        result[v] = function.value(p, component); 
      } 
    return result; 
  } 

  template <int dim, typename Number, int n_components = dim + 2> 
  Tensor<1, n_components, VectorizedArray<Number>> 
  evaluate_function(const Function<dim> &                      function, 
                    const Point<dim, VectorizedArray<Number>> &p_vectorized) 
  { 
    AssertDimension(function.n_components, n_components); 
    Tensor<1, n_components, VectorizedArray<Number>> result; 
    for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v) 
      { 
        Point<dim> p; 
        for (unsigned int d = 0; d < dim; ++d) 
          p[d] = p_vectorized[d][v]; 
        for (unsigned int d = 0; d < n_components; ++d) 
          result[d][v] = function.value(p, d); 
      } 
    return result; 
  } 

//  @sect3{The EulerOperation class}  

// 这个类实现了欧拉问题的评估器，类似于  step-37  或  step-59  的 `LaplaceOperator` 类。由于本算子是非线性的，不需要矩阵接口（交给预处理程序），我们跳过了无矩阵算子中的各种`vmult`函数，只实现了`apply`函数以及`apply`与上述低存储Runge-Kutta时间积分器所需的矢量更新的组合（称为`perform_stage`）。此外，我们还增加了三个涉及无矩阵例程的额外函数，即一个是根据元素中的速度和声速计算时间步长的估计值（与实际时间步长的Courant数相结合），一个是解的投影（专门针对DG情况的 VectorTools::project() ），还有一个是计算与可能的分析解或与某些背景状态的规范的误差。

// 该课的其余部分与其他无矩阵教程相似。正如介绍中所讨论的，我们提供了几个函数，允许用户在由 types::boundary_id 变量标记的领域边界的不同部分传递各种形式的边界条件，以及可能的体力。

  template <int dim, int degree, int n_points_1d> 
  class EulerOperator 
  { 
  public: 
    static constexpr unsigned int n_quadrature_points_1d = n_points_1d; 

    EulerOperator(TimerOutput &timer_output); 

    void reinit(const Mapping<dim> &   mapping, 
                const DoFHandler<dim> &dof_handler); 

    void set_inflow_boundary(const types::boundary_id       boundary_id, 
                             std::unique_ptr<Function<dim>> inflow_function); 

    void set_subsonic_outflow_boundary( 
      const types::boundary_id       boundary_id, 
      std::unique_ptr<Function<dim>> outflow_energy); 

    void set_wall_boundary(const types::boundary_id boundary_id); 

    void set_body_force(std::unique_ptr<Function<dim>> body_force); 

    void apply(const double                                      current_time, 
               const LinearAlgebra::distributed::Vector<Number> &src, 
               LinearAlgebra::distributed::Vector<Number> &      dst) const; 

    void 
    perform_stage(const Number cur_time, 
                  const Number factor_solution, 
                  const Number factor_ai, 
                  const LinearAlgebra::distributed::Vector<Number> &current_ri, 
                  LinearAlgebra::distributed::Vector<Number> &      vec_ki, 
                  LinearAlgebra::distributed::Vector<Number> &      solution, 
                  LinearAlgebra::distributed::Vector<Number> &next_ri) const; 

    void project(const Function<dim> &                       function, 
                 LinearAlgebra::distributed::Vector<Number> &solution) const; 

    std::array<double, 3> compute_errors( 
      const Function<dim> &                             function, 
      const LinearAlgebra::distributed::Vector<Number> &solution) const; 

    double compute_cell_transport_speed( 
      const LinearAlgebra::distributed::Vector<Number> &solution) const; 

    void 
    initialize_vector(LinearAlgebra::distributed::Vector<Number> &vector) const; 

  private: 
    MatrixFree<dim, Number> data; 

    TimerOutput &timer; 

    std::map<types::boundary_id, std::unique_ptr<Function<dim>>> 
      inflow_boundaries; 
    std::map<types::boundary_id, std::unique_ptr<Function<dim>>> 
                                   subsonic_outflow_boundaries; 
    std::set<types::boundary_id>   wall_boundaries; 
    std::unique_ptr<Function<dim>> body_force; 

    void local_apply_inverse_mass_matrix( 
      const MatrixFree<dim, Number> &                   data, 
      LinearAlgebra::distributed::Vector<Number> &      dst, 
      const LinearAlgebra::distributed::Vector<Number> &src, 
      const std::pair<unsigned int, unsigned int> &     cell_range) const; 

    void local_apply_cell( 
      const MatrixFree<dim, Number> &                   data, 
      LinearAlgebra::distributed::Vector<Number> &      dst, 
      const LinearAlgebra::distributed::Vector<Number> &src, 
      const std::pair<unsigned int, unsigned int> &     cell_range) const; 

    void local_apply_face( 
      const MatrixFree<dim, Number> &                   data, 
      LinearAlgebra::distributed::Vector<Number> &      dst, 
      const LinearAlgebra::distributed::Vector<Number> &src, 
      const std::pair<unsigned int, unsigned int> &     face_range) const; 

    void local_apply_boundary_face( 
      const MatrixFree<dim, Number> &                   data, 
      LinearAlgebra::distributed::Vector<Number> &      dst, 
      const LinearAlgebra::distributed::Vector<Number> &src, 
      const std::pair<unsigned int, unsigned int> &     face_range) const; 
  }; 

  template <int dim, int degree, int n_points_1d> 
  EulerOperator<dim, degree, n_points_1d>::EulerOperator(TimerOutput &timer) 
    : timer(timer) 
  {} 

// 对于欧拉算子的初始化，我们设置了类中包含的MatrixFree变量。这可以通过给定一个描述可能的弯曲边界的映射以及一个描述自由度的DoFHandler对象来完成。由于我们在这个教程程序中使用的是不连续的Galerkin离散化，没有对解场施加强烈的约束，所以我们不需要传入AffineConstraints对象，而是使用一个假的来构造。关于正交，我们要选择两种不同的方式来计算基础积分。第一种是灵活的，基于模板参数`n_points_1d`（将被分配到本文件顶部指定的`n_q_points_1d`值）。更精确的积分是必要的，以避免由于欧拉算子中的可变系数而产生的混叠问题。第二个不太精确的正交公式是一个基于`fe_degree+1`的严密公式，需要用于反质量矩阵。虽然该公式只在仿生元素形状上提供了精确的反，而在变形元素上则没有，但它可以通过张量积技术快速反转质量矩阵，这对于确保整体的最佳计算效率是必要的。

  template <int dim, int degree, int n_points_1d> 
  void EulerOperator<dim, degree, n_points_1d>::reinit( 
    const Mapping<dim> &   mapping, 
    const DoFHandler<dim> &dof_handler) 
  { 
    const std::vector<const DoFHandler<dim> *> dof_handlers = {&dof_handler}; 
    const AffineConstraints<double>            dummy; 
    const std::vector<const AffineConstraints<double> *> constraints = {&dummy}; 
    const std::vector<Quadrature<1>> quadratures = {QGauss<1>(n_q_points_1d), 
                                                    QGauss<1>(fe_degree + 1)}; 

    typename MatrixFree<dim, Number>::AdditionalData additional_data; 
    additional_data.mapping_update_flags = 
      (update_gradients | update_JxW_values | update_quadrature_points | 
       update_values); 
    additional_data.mapping_update_flags_inner_faces = 
      (update_JxW_values | update_quadrature_points | update_normal_vectors | 
       update_values); 
    additional_data.mapping_update_flags_boundary_faces = 
      (update_JxW_values | update_quadrature_points | update_normal_vectors | 
       update_values); 
    additional_data.tasks_parallel_scheme = 
      MatrixFree<dim, Number>::AdditionalData::none; 

    data.reinit( 
      mapping, dof_handlers, constraints, quadratures, additional_data); 
  } 

  template <int dim, int degree, int n_points_1d> 
  void EulerOperator<dim, degree, n_points_1d>::initialize_vector( 
    LinearAlgebra::distributed::Vector<Number> &vector) const 
  { 
    data.initialize_dof_vector(vector); 
  } 

// 随后的四个成员函数是必须从外部调用的，以指定各种类型的边界。对于一个流入的边界，我们必须以密度  $\rho$  、动量  $\rho \mathbf{u}$  和能量  $E$  来指定所有成分。考虑到这些信息，我们将函数与各自的边界ID一起存储在这个类的地图成员变量中。同样，我们对亚音速外流边界（我们也要求一个函数，用来检索能量）和壁面（无穿透）边界进行处理，在壁面上我们施加零法线速度（不需要函数，所以我们只要求边界ID）。对于目前的DG代码来说，边界条件只作为弱形式的一部分被应用（在时间积分期间），设置边界条件的调用可以出现在对这个类的`reinit()`调用之前或之后。这与连续有限元代码不同，在连续有限元代码中，边界条件决定了被送入MatrixFree初始化的AffineConstraints对象的内容，因此需要在无矩阵数据结构的初始化之前设置。

// 在四个函数中的每一个中添加的检查是用来确保边界条件在边界的各个部分是相互排斥的，也就是说，用户不会意外地将一个边界既指定为流入边界，又指定为亚声速流出边界。

  template <int dim, int degree, int n_points_1d> 
  void EulerOperator<dim, degree, n_points_1d>::set_inflow_boundary( 
    const types::boundary_id       boundary_id, 
    std::unique_ptr<Function<dim>> inflow_function) 
  { 
    AssertThrow(subsonic_outflow_boundaries.find(boundary_id) == 
                    subsonic_outflow_boundaries.end() && 
                  wall_boundaries.find(boundary_id) == wall_boundaries.end(), 
 
 
 
 
 
                ExcMessage("Expected function with dim+2 components")); 

    inflow_boundaries[boundary_id] = std::move(inflow_function); 
  } 

  template <int dim, int degree, int n_points_1d> 
  void EulerOperator<dim, degree, n_points_1d>::set_subsonic_outflow_boundary( 
    const types::boundary_id       boundary_id, 
    std::unique_ptr<Function<dim>> outflow_function) 
  { 
    AssertThrow(inflow_boundaries.find(boundary_id) == 
                    inflow_boundaries.end() && 
                  wall_boundaries.find(boundary_id) == wall_boundaries.end(), 
                ExcMessage("You already set the boundary with id " + 
                           std::to_string(static_cast<int>(boundary_id)) + 
                           " to another type of boundary before now setting " + 
                           "it as subsonic outflow")); 
    AssertThrow(outflow_function->n_components == dim + 2, 
                ExcMessage("Expected function with dim+2 components")); 

    subsonic_outflow_boundaries[boundary_id] = std::move(outflow_function); 
  } 

  template <int dim, int degree, int n_points_1d> 
  void EulerOperator<dim, degree, n_points_1d>::set_wall_boundary( 
    const types::boundary_id boundary_id) 
  { 
    AssertThrow(inflow_boundaries.find(boundary_id) == 
                    inflow_boundaries.end() && 
                  subsonic_outflow_boundaries.find(boundary_id) == 
                    subsonic_outflow_boundaries.end(), 
                ExcMessage("You already set the boundary with id " + 
                           std::to_string(static_cast<int>(boundary_id)) + 
                           " to another type of boundary before now setting " + 
                           "it as wall boundary")); 

    wall_boundaries.insert(boundary_id); 
  } 

  template <int dim, int degree, int n_points_1d> 
  void EulerOperator<dim, degree, n_points_1d>::set_body_force( 
    std::unique_ptr<Function<dim>> body_force) 
  { 
    AssertDimension(body_force->n_components, dim); 

    this->body_force = std::move(body_force); 
  } 

//  @sect4{Local evaluators}  

// 现在我们开始研究欧拉问题的局部评估器。评估器相对简单，遵循  step-37  、  step-48  或  step-59  中提出的内容。第一个显著的区别是，我们使用的是具有非标准正交点数量的FEE评估。以前我们总是将正交点的数量设置为等于多项式度数加1（确保在仿生元素形状上的精确积分），现在我们将正交点的数量设置为一个单独的变量（例如多项式度数加多项式度数的二分之一或三分之一），以更准确地处理非线性项。由于评估器通过模板参数输入了适当的循环长度，并在变量 FEEvaluation::n_q_points, 中保留了整个单元格的正交点数量，所以我们现在自动操作更精确的公式，而无需进一步修改。

// 第二个区别是由于我们现在评估的是一个多分量系统，而不是之前考虑的标量系统。无矩阵框架提供了几种方法来处理多成分的情况。这里显示的变体是利用一个嵌入了多个分量的FEEvaluation对象，由第四个模板参数`dim + 2`指定欧拉系统中的分量。因此， FEEvaluation::get_value() 的返回类型不再是一个标量（这将返回一个VectorizedArray类型，收集几个元素的数据），而是一个`dim+2`组件的张量。该功能与标量的情况类似；它由一个基类的模板专业化处理，称为FEEvaluationAccess。另一个变体是使用几个FEEvaluation对象，一个标量对象用于密度，一个带`dim`分量的矢量值对象用于动量，另一个标量评价器用于能量。为了确保这些分量指向解决方案的正确部分，FEEvaluation的构造函数在所需的MatrixFree字段之后需要三个可选的整数参数，即多DoFHandler系统的DoFHandler编号（默认取第一个），如果有多个Quadrature对象，则取正交点的编号（见下文），以及作为第三个参数的矢量系统中的分量。由于我们有一个单一的矢量来表示所有的分量，我们将使用第三个参数，并将其设置为`0`表示密度，`1`表示矢量值的动量，`dim+1`表示能量槽。然后FEEvaluation在 FEEvaluationBase::read_dof_values() 和 FEEvaluation::distributed_local_to_global() 或更紧凑的 FEEvaluation::gather_evaluate() 和 FEEvaluation::integrate_scatter() 调用中挑选适当的解矢量子范围。

// 当涉及到身体力向量的评估时，为了效率，我们区分了两种情况。如果我们有一个常数函数（源自 Functions::ConstantFunction), ），我们可以在正交点的循环外预先计算出数值，并简单地在所有地方使用该数值。对于一个更通用的函数，我们反而需要调用我们上面提供的`evaluate_function()`方法；这个路径更昂贵，因为我们需要访问与正交点数据有关的内存。

// 其余部分沿用其他教程的程序。由于我们已经在单独的`euler_flux()`函数中实现了欧拉方程的所有物理学，我们在这里所要做的就是给定在正交点评估的当前解，由`phi.get_value(q)`返回，并告诉FEEvaluation对象，通过形状函数的梯度（这是一个外部`dim+2`分量的张量，每个张量持有一个`dim`分量的 $x,y,z$  ] 欧拉通量的分量）。) 最后值得一提的是，在我们得到一个外部函数的情况下，我们通过测试函数`phi.submit_value()`的值来排队测试数据的顺序。我们必须在调用`phi.get_value(q)'之后进行，因为`get_value()'（读取解决方案）和`submit_value()'（排队等待测试函数的乘法和正交点的求和）访问同一个底层数据域。这里很容易实现没有临时变量`w_q`，因为值和梯度之间没有混合。对于更复杂的设置，必须首先复制出例如正交点的值和梯度，然后通过 FEEvaluationBase::submit_value() 和 FEEvaluationBase::submit_gradient(). 再次排列结果。

// 作为最后的说明，我们提到我们没有使用这个函数的第一个MatrixFree参数，这是一个来自 MatrixFree::loop(). 的回调，接口规定了现在的参数列表，但是由于我们在一个成员函数中，MatrixFree对象已经可以作为`data`变量，我们坚持使用，以避免混淆。

  template <int dim, int degree, int n_points_1d> 
  void EulerOperator<dim, degree, n_points_1d>::local_apply_cell( 
    const MatrixFree<dim, Number> &, 
    LinearAlgebra::distributed::Vector<Number> &      dst, 
    const LinearAlgebra::distributed::Vector<Number> &src, 
    const std::pair<unsigned int, unsigned int> &     cell_range) const 
  { 
    FEEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data); 

    Tensor<1, dim, VectorizedArray<Number>> constant_body_force; 
    const Functions::ConstantFunction<dim> *constant_function = 
      dynamic_cast<Functions::ConstantFunction<dim> *>(body_force.get()); 

    if (constant_function) 
      constant_body_force = evaluate_function<dim, Number, dim>( 
        *constant_function, Point<dim, VectorizedArray<Number>>()); 

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) 
      { 
        phi.reinit(cell); 
        phi.gather_evaluate(src, EvaluationFlags::values); 

        for (unsigned int q = 0; q < phi.n_q_points; ++q) 
          { 
            const auto w_q = phi.get_value(q); 
            phi.submit_gradient(euler_flux<dim>(w_q), q); 
            if (body_force.get() != nullptr) 
              { 
                const Tensor<1, dim, VectorizedArray<Number>> force = 
                  constant_function ? constant_body_force : 
                                      evaluate_function<dim, Number, dim>( 
                                        *body_force, phi.quadrature_point(q)); 

                Tensor<1, dim + 2, VectorizedArray<Number>> forcing; 
                for (unsigned int d = 0; d < dim; ++d) 
                  forcing[d + 1] = w_q[0] * force[d]; 
                for (unsigned int d = 0; d < dim; ++d) 
                  forcing[dim + 1] += force[d] * w_q[d + 1]; 

                phi.submit_value(forcing, q); 
              } 
          } 

        phi.integrate_scatter(((body_force.get() != nullptr) ? 
                                 EvaluationFlags::values : 
                                 EvaluationFlags::nothing) | 
                                EvaluationFlags::gradients, 
                              dst); 
      } 
  } 

// 下一个函数涉及到内部面的积分计算，在这里我们需要与面相邻的两个单元的评估器。我们将变量`phi_m`与解分量 $\mathbf{w}^-$ 相关联，将变量`phi_p`与解分量 $\mathbf{w}^+$ 相关联。我们在FEFaceEvaluation的构造函数中通过第二个参数来区分两边，`true`表示内侧，`false`表示外侧，内侧和外侧表示相对于法向量的方向。

// 注意调用 FEFaceEvaluation::gather_evaluate() 和 FEFaceEvaluation::integrate_scatter() 结合了对向量的访问和因式分解部分。这种合并操作不仅节省了一行代码，而且还包含了一个重要的优化。鉴于我们在Gauss-Lobatto正交公式的点上使用拉格朗日多项式的节点基础，在每个面上只有 $(p+1)^{d-1}$ 的基础函数评估为非零。因此，评估器只访问了向量中的必要数据，而跳过了乘以零的部分。如果我们首先读取向量，我们就需要从向量中加载所有的数据，因为孤立的调用不知道后续操作中需要哪些数据。如果随后的 FEFaceEvaluation::evaluate() 调用要求数值和导数，确实需要每个分量的所有 $(p+1)^d$ 向量条目，因为所有基函数的法向导数都是非零的。

// 评价器的参数以及程序与单元评价相似。由于非线性项的存在，我们再次使用更精确的（过度）积分方案，指定为列表中第三个模板参数。在正交点上，我们再去找我们的自由函数来计算数值通量。它从两边（即 $\mathbf{w}^-$ 和 $\mathbf{w}^+$ ）接收在正交点评估的解决方案，以及到减去一边的法向量。正如上面所解释的，数值通量已经乘以来自减法侧的法向量了。我们需要转换符号，因为在引言中得出的弱形式中，边界项带有一个减号。然后，通量被排队在减号和加号上进行测试，由于加号上的法向量与减号上的法向量正好相反，所以要调换符号。

  template <int dim, int degree, int n_points_1d> 
  void EulerOperator<dim, degree, n_points_1d>::local_apply_face( 
    const MatrixFree<dim, Number> &, 
    LinearAlgebra::distributed::Vector<Number> &      dst, 
    const LinearAlgebra::distributed::Vector<Number> &src, 
    const std::pair<unsigned int, unsigned int> &     face_range) const 
  { 
    FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi_m(data, 
                                                                      true); 
    FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi_p(data, 
                                                                      false); 

    for (unsigned int face = face_range.first; face < face_range.second; ++face) 
      { 
        phi_p.reinit(face); 
        phi_p.gather_evaluate(src, EvaluationFlags::values); 

        phi_m.reinit(face); 
        phi_m.gather_evaluate(src, EvaluationFlags::values); 

        for (unsigned int q = 0; q < phi_m.n_q_points; ++q) 
          { 
            const auto numerical_flux = 
              euler_numerical_flux<dim>(phi_m.get_value(q), 
                                        phi_p.get_value(q), 
                                        phi_m.get_normal_vector(q)); 
            phi_m.submit_value(-numerical_flux, q); 
            phi_p.submit_value(numerical_flux, q); 
          } 

        phi_p.integrate_scatter(EvaluationFlags::values, dst); 
        phi_m.integrate_scatter(EvaluationFlags::values, dst); 
      } 
  } 

// 对于位于边界的面，我们需要施加适当的边界条件。在这个教程程序中，我们实现了上述的四种情况。第五种情况，即超音速流出条件，将在下面的 "结果 "部分讨论）。不连续的Galerkin方法对边界条件的施加不是作为约束条件，而只是弱化。因此，各种条件是通过找到一个适当的<i>exterior</i>量 $\mathbf{w}^+$ 来施加的，然后将其交给也用于内部面的数值通量函数。实质上，我们在域外 "假装 "一个状态，如果那是现实，PDE的解将满足我们想要的边界条件。

// 对于墙的边界，我们需要对动量变量施加一个无正态通量的条件，而对于密度和能量，我们使用的是诺伊曼条件  $\rho^+ = \rho^-$  和  $E^+ = E^-$  。为了实现无正态通量条件，我们将外部数值设定为内部数值，并减去墙面法线方向，即法线矢量方向上的速度的2倍。

// 对于流入边界，我们简单地将给定的Dirichlet数据 $\mathbf{w}_\mathrm{D}$ 作为边界值。另一种方法是使用 $\mathbf{w}^+ = -\mathbf{w}^- + 2 \mathbf{w}_\mathrm{D}$  ，即所谓的镜像原理。

// 强加外流本质上是一个诺伊曼条件，即设定  $\mathbf{w}^+ = \mathbf{w}^-$  。对于亚声速流出的情况，我们仍然需要强加一个能量值，我们从各自的函数中得出这个值。对于<i>backflow</i>的情况，即在Neumann部分有动量通入域的情况，需要一个特殊的步骤。根据文献（这一事实可以通过适当的能量论证得出），我们必须切换到流入部分的通量的另一个变体，见Gravemeier, Comerford, Yoshihara, Ismail, Wall, "A novel formulation for Neumann inflow conditions in biomechanics", Int. J. Numer. Meth. 生物医学。Eng., vol. 28 (2012). 这里，动量项需要再次添加，这相当于去除动量变量上的通量贡献。我们在后处理步骤中这样做，而且只适用于我们都处于外流边界且法向量与动量（或等同于速度）之间的点积为负的情况。由于我们在SIMD矢量化中一次处理多个正交点的数据，这里需要明确地在SIMD数组的条目上循环。

// 在下面的实现中，我们在正交点的层面上检查各种类型的边界。当然，我们也可以将决定权移出正交点循环，将整个面孔视为同类，这就避免了在正交点的内循环中进行一些地图/集合的查找。然而，效率的损失并不明显，所以我们在这里选择了更简单的代码。还要注意的是，最后的 "else "子句会捕捉到这样的情况，即边界的某些部分没有通过 `EulerOperator::set_..._boundary(...)`. 分配任何边界条件。
  template <int dim, int degree, int n_points_1d> 
  void EulerOperator<dim, degree, n_points_1d>::local_apply_boundary_face( 
    const MatrixFree<dim, Number> &, 
    LinearAlgebra::distributed::Vector<Number> &      dst, 
    const LinearAlgebra::distributed::Vector<Number> &src, 
    const std::pair<unsigned int, unsigned int> &     face_range) const 
  { 
    FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data, true); 

    for (unsigned int face = face_range.first; face < face_range.second; ++face) 
      { 
        phi.reinit(face); 
        phi.gather_evaluate(src, EvaluationFlags::values); 

        for (unsigned int q = 0; q < phi.n_q_points; ++q) 
          { 
            const auto w_m    = phi.get_value(q); 
            const auto normal = phi.get_normal_vector(q); 

            auto rho_u_dot_n = w_m[1] * normal[0]; 
            for (unsigned int d = 1; d < dim; ++d) 
              rho_u_dot_n += w_m[1 + d] * normal[d]; 

            bool at_outflow = false; 

            Tensor<1, dim + 2, VectorizedArray<Number>> w_p; 
            const auto boundary_id = data.get_boundary_id(face); 
            if (wall_boundaries.find(boundary_id) != wall_boundaries.end()) 
              { 
                w_p[0] = w_m[0]; 
                for (unsigned int d = 0; d < dim; ++d) 
                  w_p[d + 1] = w_m[d + 1] - 2. * rho_u_dot_n * normal[d]; 
                w_p[dim + 1] = w_m[dim + 1]; 
              } 
            else if (inflow_boundaries.find(boundary_id) != 
                     inflow_boundaries.end()) 
              w_p = 
                evaluate_function(*inflow_boundaries.find(boundary_id)->second, 
                                  phi.quadrature_point(q)); 
            else if (subsonic_outflow_boundaries.find(boundary_id) != 
                     subsonic_outflow_boundaries.end()) 
              { 
                w_p          = w_m; 
                w_p[dim + 1] = evaluate_function( 
                  *subsonic_outflow_boundaries.find(boundary_id)->second, 
                  phi.quadrature_point(q), 
                  dim + 1); 
                at_outflow = true; 
              } 
            else 
              AssertThrow(false, 
                          ExcMessage("Unknown boundary id, did " 
                                     "you set a boundary condition for " 
                                     "this part of the domain boundary?")); 

            auto flux = euler_numerical_flux<dim>(w_m, w_p, normal); 

            if (at_outflow) 
              for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v) 
                { 
                  if (rho_u_dot_n[v] < -1e-12) 
                    for (unsigned int d = 0; d < dim; ++d) 
                      flux[d + 1][v] = 0.; 
                } 

            phi.submit_value(-flux, q); 
          } 

        phi.integrate_scatter(EvaluationFlags::values, dst); 
      } 
  } 

// 下一个函数实现了质量矩阵的逆运算。在介绍中已经广泛讨论了算法和原理，所以我们在这里只讨论 MatrixFreeOperators::CellwiseInverseMassMatrix 类的技术问题。它所做的操作与质量矩阵的正向评估类似，只是使用了不同的插值矩阵，代表逆 $S^{-1}$ 因子。这些代表了从指定的基础（在这种情况下，高斯--洛巴托正交公式点中的拉格朗日基础）到高斯正交公式点中的拉格朗日基础的改变。在后者的基础上，我们可以应用点的逆向`JxW`因子，即正交权重乘以从参考坐标到实坐标的映射的雅各布系数。一旦完成了这一操作，基数将再次变回节点高斯-洛巴托基数。所有这些操作都由下面的 "apply() "函数完成。我们需要提供的是要操作的局部场（我们通过一个FEEvaluation对象从全局向量中提取），并将结果写回质量矩阵操作的目标向量。

// 需要注意的一点是，我们在FEEvaluation的构造函数中添加了两个整数参数（可选），第一个是0（在多DoFHandler系统中选择DoFHandler；在这里，我们只有一个），第二个是1，用于进行正交公式选择。由于我们将正交公式0用于非线性项的过度积分，我们使用公式1与默认的 $p+1$ （或变量名称中的`fe_degree+1`）点用于质量矩阵。这导致了对质量矩阵的平方贡献，并确保了精确的积分，正如介绍中所解释的。

  template <int dim, int degree, int n_points_1d> 
  void EulerOperator<dim, degree, n_points_1d>::local_apply_inverse_mass_matrix( 
    const MatrixFree<dim, Number> &, 
    LinearAlgebra::distributed::Vector<Number> &      dst, 
    const LinearAlgebra::distributed::Vector<Number> &src, 
    const std::pair<unsigned int, unsigned int> &     cell_range) const 
  { 
    FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1); 
    MatrixFreeOperators::CellwiseInverseMassMatrix<dim, degree, dim + 2, Number> 
      inverse(phi); 

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) 
      { 
        phi.reinit(cell); 
        phi.read_dof_values(src); 

        inverse.apply(phi.begin_dof_values(), phi.begin_dof_values()); 

        phi.set_dof_values(dst); 
      } 
  } 

//  @sect4{The apply() and related functions}  

// 我们现在来到实现欧拉算子整体评估的函数，即 $\mathcal M^{-1} \mathcal L(t, \mathbf{w})$  ，调用上面介绍的局部评估器。这些步骤在前面的代码中应该是清楚的。需要注意的一点是，我们需要调整与边界各部分相关的函数中的时间，以便在边界数据与时间相关的情况下与方程一致。然后，我们调用 MatrixFree::loop() 来执行单元和面的积分，包括在`src`向量中进行必要的ghost数据交换。该函数的第七个参数，"true"，指定我们要在开始向其累积积分之前，将 "dst "向量作为循环的一部分归零。这个变体比在循环之前明确调用`dst = 0.;`要好，因为归零操作是在矢量的子范围内完成的，其部分是由附近的积分写入的。这加强了数据的定位，并允许缓存，节省了向量数据到主内存的一次往返，提高了性能。循环的最后两个参数决定了哪些数据被交换：由于我们只访问一个面的形状函数的值，这是典型的一阶双曲问题，并且由于我们有一个节点基础，节点位于参考元素表面，我们只需要交换这些部分。这又节省了宝贵的内存带宽。

// 一旦应用了空间算子 $\mathcal L$ ，我们需要进行第二轮操作，应用反质量矩阵。这里，我们调用 MatrixFree::cell_loop() ，因为只有单元格积分出现。单元循环比全循环更便宜，因为只访问与本地拥有的单元相关的自由度，这只是DG离散化的本地拥有的自由度。因此，这里不需要鬼魂交换。

// 在所有这些函数的周围，我们设置了定时器范围来记录计算时间，以统计各部分的贡献。

  template <int dim, int degree, int n_points_1d> 
  void EulerOperator<dim, degree, n_points_1d>::apply( 
    const double                                      current_time, 
    const LinearAlgebra::distributed::Vector<Number> &src, 
    LinearAlgebra::distributed::Vector<Number> &      dst) const 
  { 
    { 
      TimerOutput::Scope t(timer, "apply - integrals"); 

      for (auto &i : inflow_boundaries) 
        i.second->set_time(current_time); 
      for (auto &i : subsonic_outflow_boundaries) 
        i.second->set_time(current_time); 

      data.loop(&EulerOperator::local_apply_cell, 
                &EulerOperator::local_apply_face, 
                &EulerOperator::local_apply_boundary_face, 
                this, 
                dst, 
                src, 
                true, 
                MatrixFree<dim, Number>::DataAccessOnFaces::values, 
                MatrixFree<dim, Number>::DataAccessOnFaces::values); 
    } 

    { 
      TimerOutput::Scope t(timer, "apply - inverse mass"); 

      data.cell_loop(&EulerOperator::local_apply_inverse_mass_matrix, 
                     this, 
                     dst, 
                     dst); 
    } 
  } 

// 让我们转到做Runge--Kutta更新的整个阶段的函数。它调用 EulerOperator::apply() ，然后对向量进行一些更新，即`next_ri = solution + factor_ai * k_i`和`solution += factor_solution * k_i`。与其通过向量接口执行这些步骤，我们在这里提出了一个替代策略，在基于缓存的架构上速度更快。由于向量所消耗的内存往往比缓存所能容纳的要大得多，因此数据必须有效地来自缓慢的RAM内存。这种情况可以通过循环融合来改善，即在一次扫描中对`next_ki`和`solution`进行更新。在这种情况下，我们将读取两个向量`rhs`和`solution`并写入`next_ki`和`solution`，而在基线情况下，至少有4次读取和两次写入。在这里，我们更进一步，当质量矩阵反转在向量的某一部分完成后，立即执行循环。  MatrixFree::cell_loop() 提供了一种机制，在单元格的循环第一次接触到一个向量条目之前，附加一个 `std::function` （我们在这里没有使用，但用于例如向量的归零），以及在循环最后接触到一个条目之后，调用第二个 `std::function` 。回调的形式是给定向量上的一个范围（就MPI宇宙中的本地索引编号而言），可以由`local_element()`函数来处理。

// 对于这个第二个回调，我们创建一个lambda，在一个范围内工作，并在这个范围内写入相应的更新。理想情况下，我们会在本地循环之前添加`DEAL_II_OPENMP_SIMD_PRAGMA`，以建议编译器对这个循环进行SIMD并行化（这意味着在实践中我们要确保在循环内部使用的指针的索引范围之间没有重叠，也称为别名）。事实证明，在写这篇文章的时候，GCC 7.2无法编译lambda函数中的OpenMP pragma，所以我们在下面注释了这个pragma。如果你的编译器比较新，你应该可以再次取消注释这些行。

// 注意，当我们不需要更新`next_ri`向量时，我们为最后的Runge--Kutta阶段选择不同的代码路径。这个策略带来了相当大的速度提升。在40核机器上，默认矢量更新时，逆质量矩阵和矢量更新需要60%以上的计算时间，而在更优化的变体中，这一比例约为35%。换句话说，这是一个大约三分之一的速度提升。

  template <int dim, int degree, int n_points_1d> 
  void EulerOperator<dim, degree, n_points_1d>::perform_stage( 
    const Number                                      current_time, 
    const Number                                      factor_solution, 
    const Number                                      factor_ai, 
    const LinearAlgebra::distributed::Vector<Number> &current_ri, 
    LinearAlgebra::distributed::Vector<Number> &      vec_ki, 
    LinearAlgebra::distributed::Vector<Number> &      solution, 
    LinearAlgebra::distributed::Vector<Number> &      next_ri) const 
  { 
    { 
      TimerOutput::Scope t(timer, "rk_stage - integrals L_h"); 

      for (auto &i : inflow_boundaries) 
        i.second->set_time(current_time); 
      for (auto &i : subsonic_outflow_boundaries) 
        i.second->set_time(current_time); 

      data.loop(&EulerOperator::local_apply_cell, 
                &EulerOperator::local_apply_face, 
                &EulerOperator::local_apply_boundary_face, 
                this, 
                vec_ki, 
                current_ri, 
                true, 
                MatrixFree<dim, Number>::DataAccessOnFaces::values, 
                MatrixFree<dim, Number>::DataAccessOnFaces::values); 
    } 

    { 
      TimerOutput::Scope t(timer, "rk_stage - inv mass + vec upd"); 
      data.cell_loop( 
        &EulerOperator::local_apply_inverse_mass_matrix, 
        this, 
        next_ri, 
        vec_ki, 
        std::function<void(const unsigned int, const unsigned int)>(), 
        [&](const unsigned int start_range, const unsigned int end_range) { 
          const Number ai = factor_ai; 
          const Number bi = factor_solution; 
          if (ai == Number()) 
            { 

          /* DEAL_II_OPENMP_SIMD_PRAGMA  */ 
              for (unsigned int i = start_range; i < end_range; ++i) 
                { 
                  const Number k_i          = next_ri.local_element(i); 
                  const Number sol_i        = solution.local_element(i); 
                  solution.local_element(i) = sol_i + bi * k_i; 
                } 
            } 
          else 
            { 

              /* DEAL_II_OPENMP_SIMD_PRAGMA  */ 
              for (unsigned int i = start_range; i < end_range; ++i) 
                { 
                  const Number k_i          = next_ri.local_element(i); 
                  const Number sol_i        = solution.local_element(i); 
                  solution.local_element(i) = sol_i + bi * k_i; 
                  next_ri.local_element(i)  = sol_i + ai * k_i; 
                } 
            } 
        }); 
    } 
  } 

// 在讨论了将解提前一个时间步长的函数的实现后，现在让我们来看看实现其他辅助性操作的函数。具体来说，这些是计算投影、评估误差和计算单元上信息传输速度的函数。

// 这些函数中的第一个基本上等同于 VectorTools::project(), ，只是速度快得多，因为它是专门针对DG元素的，不需要设置和解决线性系统，因为每个元素都有独立的基函数。我们在这里展示代码的原因，除了这个非关键操作的小幅提速之外，还因为它显示了 MatrixFreeOperators::CellwiseInverseMassMatrix. 提供的额外功能。

// 投影操作的工作原理如下。如果我们用 $S$ 表示在正交点评估的形状函数矩阵，那么在单元格 $K$ 上的投影是一个形式为 $\underbrace{S J^K S^\mathrm T}_{\mathcal M^K} \mathbf{w}^K = S J^K \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$ 的操作，其中 $J^K$ 是包含雅各布系数乘以正交权重（JxW）的对角矩阵， $\mathcal M^K$ 是单元格的质量矩阵， $\tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$ 是要投影到正交点的领域评估。实际上，矩阵 $S$ 通过张量积有额外的结构，如介绍中所解释的）。这个系统现在可以等效地写成 $\mathbf{w}^K = \left(S J^K S^\mathrm T\right)^{-1} S J^K \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q} = S^{-\mathrm T} \left(J^K\right)^{-1} S^{-1} S J^K \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$  。现在，项 $S^{-1} S$ 和 $\left(J^K\right)^{-1} J^K$ 相抵消，导致最后的表达式 $\mathbf{w}^K = S^{-\mathrm T} \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$  。这个操作由 MatrixFreeOperators::CellwiseInverseMassMatrix::transform_from_q_points_to_basis(). 实现。这个名字来自于这个投影只是乘以 $S^{-\mathrm T}$ ，一个从高斯正交点的节点基到给定的有限元基的基数变化。请注意，我们调用 FEEvaluation::set_dof_values() 将结果写入矢量，覆盖之前的内容，而不是像典型的积分任务那样累积结果--我们可以这样做，因为对于不连续的Galerkin离散，每个矢量条目都只有一个单元的贡献。

  template <int dim, int degree, int n_points_1d> 
  void EulerOperator<dim, degree, n_points_1d>::project( 
    const Function<dim> &                       function, 
    LinearAlgebra::distributed::Vector<Number> &solution) const 
  { 
    FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1); 
    MatrixFreeOperators::CellwiseInverseMassMatrix<dim, degree, dim + 2, Number> 
      inverse(phi); 
    solution.zero_out_ghost_values(); 
    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell) 
      { 
        phi.reinit(cell); 
        for (unsigned int q = 0; q < phi.n_q_points; ++q) 
          phi.submit_dof_value(evaluate_function(function, 
                                                 phi.quadrature_point(q)), 
                               q); 
        inverse.transform_from_q_points_to_basis(dim + 2, 
                                                 phi.begin_dof_values(), 
                                                 phi.begin_dof_values()); 
        phi.set_dof_values(solution); 
      } 
  } 

// 下一个函数再次重复了同样由deal.II库提供的功能，即 VectorTools::integrate_difference(). 我们在这里展示了明确的代码，以强调跨几个单元的矢量化是如何工作的，以及如何通过该接口累积结果。回顾一下，每个<i>lane</i>的矢量化数组持有来自不同单元的数据。通过对当前MPI进程所拥有的所有单元批的循环，我们就可以填充一个结果的VectorizedArray；为了得到一个全局的总和，我们需要进一步去对SIMD阵列中的条目进行求和。然而，这样的程序并不稳定，因为SIMD数组事实上可能并不持有其所有通道的有效数据。当本地拥有的单元的数量不是SIMD宽度的倍数时，就会发生这种情况。为了避免无效数据，我们必须在访问数据时明确地跳过那些无效的通道。虽然人们可以想象，我们可以通过简单地将空车道设置为零（从而不对总和做出贡献）来使其工作，但情况比这更复杂。如果我们要从动量中计算出一个速度呢？那么，我们就需要除以密度，而密度是零--结果就会是NaN，并污染结果。当我们在单元格批次中循环时，使用函数 MatrixFree::n_active_entries_per_cell_batch() 给我们提供有效数据的通道数，累积有效SIMD范围内的结果，就可以避免这种陷阱。它在大多数单元上等于 VectorizedArray::size() ，但如果单元数与SIMD宽度相比有余数，则在最后一个单元批上可能会更少。

  template <int dim, int degree, int n_points_1d> 
  std::array<double, 3> EulerOperator<dim, degree, n_points_1d>::compute_errors( 
    const Function<dim> &                             function, 
    const LinearAlgebra::distributed::Vector<Number> &solution) const 
  { 
    TimerOutput::Scope t(timer, "compute errors"); 
    double             errors_squared[3] = {}; 
    FEEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data, 0, 0); 

    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell) 
      { 
        phi.reinit(cell); 
        phi.gather_evaluate(solution, EvaluationFlags::values); 
        VectorizedArray<Number> local_errors_squared[3] = {}; 
        for (unsigned int q = 0; q < phi.n_q_points; ++q) 
          { 
            const auto error = 
              evaluate_function(function, phi.quadrature_point(q)) - 
              phi.get_value(q); 
            const auto JxW = phi.JxW(q); 

            local_errors_squared[0] += error[0] * error[0] * JxW; 
            for (unsigned int d = 0; d < dim; ++d) 
              local_errors_squared[1] += (error[d + 1] * error[d + 1]) * JxW; 
            local_errors_squared[2] += (error[dim + 1] * error[dim + 1]) * JxW; 
          } 
        for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell); 
             ++v) 
          for (unsigned int d = 0; d < 3; ++d) 
            errors_squared[d] += local_errors_squared[d][v]; 
      } 

    Utilities::MPI::sum(errors_squared, MPI_COMM_WORLD, errors_squared); 

    std::array<double, 3> errors; 
    for (unsigned int d = 0; d < 3; ++d) 
      errors[d] = std::sqrt(errors_squared[d]); 

    return errors; 
  } 

// EulerOperator类的最后一个函数是用来估计传输速度的，由网格大小缩放，这与设置显式时间积分器的时间步长有关。在欧拉方程中，有两种传输速度，即对流速度 $\mathbf{u}$ 和相对于以速度 $\mathbf u$ 运动的介质而言，声波的传播速度 $c = \sqrt{\gamma p/\rho}$  。

// 在时间步长的公式中，我们感兴趣的不是这些绝对速度，而是信息穿过一个单元所需的时间量。对于与介质一起传输的信息， $\mathbf u$ 是由网格大小缩放的，所以最大速度的估计可以通过计算 $\|J^{-\mathrm T} \mathbf{u}\|_\infty$  得到，其中 $J$ 是实域到参考域的转换的雅各布。请注意， FEEvaluationBase::inverse_jacobian() 返回的是反转和转置的雅各布，代表从实数到参考坐标的度量项，所以我们不需要再次转置。我们在下面的代码中把这个极限存储在变量`convective_limit`中。

// 声音的传播是各向同性的，所以我们需要考虑到任何方向的网格尺寸。然后，适当的网格大小比例由 $J$ 的最小奇异值给出，或者，等同于 $J^{-1}$ 的最大奇异值。请注意，当忽略弯曲的单元时，可以用单元顶点之间的最小距离来近似这个量。为了得到Jacobian的最大奇异值，一般的策略是使用一些LAPACK函数。由于我们在这里需要的只是一个估计值，所以我们可以避免将一个向量数组的张量分解成几个矩阵的麻烦，并在没有向量的情况下进入一个（昂贵的）特征值函数，而是使用应用于 $J^{-1}J^{-\mathrm T}$ 的幂方法进行几次迭代（在下面的代码中为五次）。这种方法的收敛速度取决于最大特征值与次大特征值的比率以及初始猜测，即所有1的矢量。这可能表明，我们在接近立方体形状的单元上得到缓慢的收敛，在这种情况下，所有的长度几乎都是一样的。然而，这种缓慢的收敛意味着结果将位于两个最大的奇异值之间，而这两个奇异值无论如何都是接近最大值的。在所有其他情况下，收敛将是快速的。因此，我们可以只在这里硬编码5次迭代，并确信结果是好的。

  template <int dim, int degree, int n_points_1d> 
  double EulerOperator<dim, degree, n_points_1d>::compute_cell_transport_speed( 
    const LinearAlgebra::distributed::Vector<Number> &solution) const 
  { 
    TimerOutput::Scope t(timer, "compute transport speed"); 
    Number             max_transport = 0; 
    FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1); 

    for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell) 
      { 
        phi.reinit(cell); 
        phi.gather_evaluate(solution, EvaluationFlags::values); 
        VectorizedArray<Number> local_max = 0.; 
        for (unsigned int q = 0; q < phi.n_q_points; ++q) 
          { 
            const auto solution = phi.get_value(q); 
            const auto velocity = euler_velocity<dim>(solution); 
            const auto pressure = euler_pressure<dim>(solution); 

            const auto inverse_jacobian = phi.inverse_jacobian(q); 
            const auto convective_speed = inverse_jacobian * velocity; 
            VectorizedArray<Number> convective_limit = 0.; 
            for (unsigned int d = 0; d < dim; ++d) 
              convective_limit = 
                std::max(convective_limit, std::abs(convective_speed[d])); 

            const auto speed_of_sound = 
              std::sqrt(gamma * pressure * (1. / solution[0])); 

            Tensor<1, dim, VectorizedArray<Number>> eigenvector; 
            for (unsigned int d = 0; d < dim; ++d) 
              eigenvector[d] = 1.; 
            for (unsigned int i = 0; i < 5; ++i) 
              { 
                eigenvector = transpose(inverse_jacobian) * 
                              (inverse_jacobian * eigenvector); 
                VectorizedArray<Number> eigenvector_norm = 0.; 
                for (unsigned int d = 0; d < dim; ++d) 
                  eigenvector_norm = 
                    std::max(eigenvector_norm, std::abs(eigenvector[d])); 
                eigenvector /= eigenvector_norm; 
              } 
            const auto jac_times_ev   = inverse_jacobian * eigenvector; 
            const auto max_eigenvalue = std::sqrt( 
              (jac_times_ev * jac_times_ev) / (eigenvector * eigenvector)); 
            local_max = 
              std::max(local_max, 
                       max_eigenvalue * speed_of_sound + convective_limit); 
          } 

// 与前面的函数类似，我们必须确保只在一个单元格批次的有效单元格上积累速度。

 
             ++v) 
          for (unsigned int d = 0; d < 3; ++d) 
            max_transport = std::max(max_transport, local_max[v]); 
      } 

    max_transport = Utilities::MPI::max(max_transport, MPI_COMM_WORLD); 

    return max_transport; 
  } 

//  @sect3{The EulerProblem class}  

// 该类将EulerOperator类与时间积分器和通常的全局数据结构（如FiniteElement和DoFHandler）相结合，以实际运行Euler问题的模拟。

// 成员变量是一个三角形、一个有限元、一个映射（用于创建高阶曲面，见 step-10 ），以及一个描述自由度的DoFHandler。此外，我们还保留了上面描述的EulerOperator的实例，它将完成所有积分方面的繁重工作，以及一些时间积分的参数，如当前时间或时间步长。

// 此外，我们使用一个PostProcessor实例来向输出文件写入一些额外的信息，这与  step-33  中的做法类似。DataPostprocessor类的接口很直观，要求我们提供关于需要评估的信息（通常只有解决方案的值，除了Schlieren图，我们只在二维中启用它是有意义的），以及被评估的东西的名称。请注意，也可以通过可视化程序（如ParaView）中的计算器工具来提取大部分信息，但在写输出时就已经做了，这要方便得多。

  template <int dim> 
  class EulerProblem 
  { 
  public: 
    EulerProblem(); 

    void run(); 

  private: 
    void make_grid_and_dofs(); 

    void output_results(const unsigned int result_number); 

    LinearAlgebra::distributed::Vector<Number> solution; 

    ConditionalOStream pcout; 

#ifdef DEAL_II_WITH_P4EST 
    parallel::distributed::Triangulation<dim> triangulation; 
#else 
    Triangulation<dim> triangulation; 
#endif 

    FESystem<dim>        fe; 
    MappingQGeneric<dim> mapping; 
    DoFHandler<dim>      dof_handler; 

    TimerOutput timer; 

    EulerOperator<dim, fe_degree, n_q_points_1d> euler_operator; 

    double time, time_step; 

    class Postprocessor : public DataPostprocessor<dim> 
    { 
    public: 
      Postprocessor(); 

      virtual void evaluate_vector_field( 
        const DataPostprocessorInputs::Vector<dim> &inputs, 
        std::vector<Vector<double>> &computed_quantities) const override; 

      virtual std::vector<std::string> get_names() const override; 

      virtual std::vector< 
        DataComponentInterpretation::DataComponentInterpretation> 
      get_data_component_interpretation() const override; 

      virtual UpdateFlags get_needed_update_flags() const override; 

    private: 
      const bool do_schlieren_plot; 
    }; 
  }; 

  template <int dim> 
  EulerProblem<dim>::Postprocessor::Postprocessor() 
    : do_schlieren_plot(dim == 2) 
  {} 

// 对于字段变量的主要评估，我们首先检查数组的长度是否等于预期值（长度`2*dim+4`或`2*dim+5`来自我们在下面get_names()函数中指定的名字的大小）。然后我们在所有的评估点上循环，填充相应的信息。首先，我们填写密度 $\rho$ 、动量 $\rho \mathbf{u}$ 和能量 $E$ 的原始解变量，然后我们计算得出速度 $\mathbf u$ 、压力 $p$ 、声速 $c=\sqrt{\gamma p / \rho}$ ，以及显示 $s = |\nabla \rho|^2$ 的Schlieren图，如果它被启用。参见 step-69 中另一个创建Schlieren图的例子）。

  template <int dim> 
  void EulerProblem<dim>::Postprocessor::evaluate_vector_field( 
    const DataPostprocessorInputs::Vector<dim> &inputs, 
    std::vector<Vector<double>> &               computed_quantities) const 
  { 
    const unsigned int n_evaluation_points = inputs.solution_values.size(); 

    if (do_schlieren_plot == true) 
      Assert(inputs.solution_gradients.size() == n_evaluation_points, 
             ExcInternalError()); 

    Assert(computed_quantities.size() == n_evaluation_points, 
           ExcInternalError()); 
    Assert(inputs.solution_values[0].size() == dim + 2, ExcInternalError()); 
    Assert(computed_quantities[0].size() == 
             dim + 2 + (do_schlieren_plot == true ? 1 : 0), 
           ExcInternalError()); 

    for (unsigned int q = 0; q < n_evaluation_points; ++q) 
      { 
        Tensor<1, dim + 2> solution; 
        for (unsigned int d = 0; d < dim + 2; ++d) 
          solution[d] = inputs.solution_values[q](d); 

        const double         density  = solution[0]; 
        const Tensor<1, dim> velocity = euler_velocity<dim>(solution); 
        const double         pressure = euler_pressure<dim>(solution); 

        for (unsigned int d = 0; d < dim; ++d) 
          computed_quantities[q](d) = velocity[d]; 
        computed_quantities[q](dim)     = pressure; 
        computed_quantities[q](dim + 1) = std::sqrt(gamma * pressure / density); 

        if (do_schlieren_plot == true) 
          computed_quantities[q](dim + 2) = 
            inputs.solution_gradients[q][0] * inputs.solution_gradients[q][0]; 
      } 
  } 

  template <int dim> 
  std::vector<std::string> EulerProblem<dim>::Postprocessor::get_names() const 
  { 
    std::vector<std::string> names; 
    for (unsigned int d = 0; d < dim; ++d) 
      names.emplace_back("velocity"); 
    names.emplace_back("pressure"); 
    names.emplace_back("speed_of_sound"); 

    if (do_schlieren_plot == true) 
      names.emplace_back("schlieren_plot"); 

    return names; 
  } 

// 对于量的解释，我们有标量密度、能量、压力、声速和Schlieren图，以及动量和速度的向量。

  template <int dim> 
  std::vector<DataComponentInterpretation::DataComponentInterpretation> 
  EulerProblem<dim>::Postprocessor::get_data_component_interpretation() const 
  { 
    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      interpretation; 
    for (unsigned int d = 0; d < dim; ++d) 
      interpretation.push_back( 
        DataComponentInterpretation::component_is_part_of_vector); 
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); 
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); 

    if (do_schlieren_plot == true) 
      interpretation.push_back( 
        DataComponentInterpretation::component_is_scalar); 

    return interpretation; 
  } 

// 关于必要的更新标志，我们只需要所有数量的值，但Schlieren图除外，它是基于密度梯度的。

  template <int dim> 
  UpdateFlags EulerProblem<dim>::Postprocessor::get_needed_update_flags() const 
  { 
    if (do_schlieren_plot == true) 
      return update_values | update_gradients; 
    else 
      return update_values; 
  } 

// 这个类的构造函数并不令人惊讶。我们设置了一个基于 "MPI_COMM_WORLD "通信器的平行三角形，一个具有 "dim+2 "分量的密度、动量和能量的矢量有限元，一个与底层有限元相同程度的高阶映射，并将时间和时间步长初始化为零。

  template <int dim> 
  EulerProblem<dim>::EulerProblem() 
    : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
#ifdef DEAL_II_WITH_P4EST 
    , triangulation(MPI_COMM_WORLD) 
#endif 
    , fe(FE_DGQ<dim>(fe_degree), dim + 2) 
    , mapping(fe_degree) 
    , dof_handler(triangulation) 
    , timer(pcout, TimerOutput::never, TimerOutput::wall_times) 
    , euler_operator(timer) 
    , time(0) 
    , time_step(0) 
  {} 

// 作为一个网格，本教程程序实现了两种选择，取决于全局变量`testcase`。对于分析型变量（`testcase==0`），域是 $(0, 10) \times (-5, 5)$ ，域的四周都有迪里希特边界条件（流入）。对于 "testcase==1"，我们将域设置为矩形箱中的圆柱体，源自Sch&auml;fer和Turek（1996）对不可压缩的粘性流动的圆柱体的流动测试案例。在这里，我们有更多种类的边界。通道左侧的流入部分是给定的流入类型，为此我们选择了一个恒定的流入轮廓，而我们在右侧设置了一个亚声速的流出。对于圆柱体周围的边界（边界id等于2）以及通道壁（边界id等于3），我们使用壁的边界类型，即无正态流。此外，对于三维圆柱体，我们还在垂直方向上增加了一个重力。有了基础网格（包括由 GridGenerator::channel_with_cylinder()), 设置的流形），我们就可以执行指定数量的全局细化，从DoFHandler创建未知的编号，并将DoFHandler和Mapping对象交给EulerOperator的初始化。

  template <int dim> 
  void EulerProblem<dim>::make_grid_and_dofs() 
  { 
    switch (testcase) 
      { 
        case 0: 
          { 
            Point<dim> lower_left; 
            for (unsigned int d = 1; d < dim; ++d) 
              lower_left[d] = -5; 

            Point<dim> upper_right; 
            upper_right[0] = 10; 
            for (unsigned int d = 1; d < dim; ++d) 
              upper_right[d] = 5; 

            GridGenerator::hyper_rectangle(triangulation, 
                                           lower_left, 
                                           upper_right); 
            triangulation.refine_global(2); 

            euler_operator.set_inflow_boundary( 
              0, std::make_unique<ExactSolution<dim>>(0)); 

            break; 
          } 

        case 1: 
          { 
            GridGenerator::channel_with_cylinder( 
              triangulation, 0.03, 1, 0, true); 

            euler_operator.set_inflow_boundary( 
              0, std::make_unique<ExactSolution<dim>>(0)); 
            euler_operator.set_subsonic_outflow_boundary( 
              1, std::make_unique<ExactSolution<dim>>(0)); 

            euler_operator.set_wall_boundary(2); 
            euler_operator.set_wall_boundary(3); 

            if (dim == 3) 
              euler_operator.set_body_force( 
                std::make_unique<Functions::ConstantFunction<dim>>( 
                  std::vector<double>({0., 0., -0.2}))); 

            break; 
          } 

        default: 
          Assert(false, ExcNotImplemented()); 
      } 

    triangulation.refine_global(n_global_refinements); 

    dof_handler.distribute_dofs(fe); 

    euler_operator.reinit(mapping, dof_handler); 
    euler_operator.initialize_vector(solution); 

// 在下文中，我们输出一些关于问题的统计数据。因为我们经常会出现相当多的单元格或自由度，所以我们希望用逗号来分隔每一组的三位数来打印它们。这可以通过 "locales "来实现，尽管这种工作方式不是特别直观。  step-32 对此有稍微详细的解释。

    std::locale s = pcout.get_stream().getloc(); 
    pcout.get_stream().imbue(std::locale("")); 
    pcout << "Number of degrees of freedom: " << dof_handler.n_dofs() 
          << " ( = " << (dim + 2) << " [vars] x " 
          << triangulation.n_global_active_cells() << " [cells] x " 
          << Utilities::pow(fe_degree + 1, dim) << " [dofs/cell/var] )" 
          << std::endl; 
    pcout.get_stream().imbue(s); 
  } 

// 对于输出，我们首先让欧拉算子计算出数值结果的误差。更确切地说，对于分析解的情况，我们计算与分析结果的误差，而对于第二个测试情况，我们计算与密度和能量恒定的背景场以及 $x$ 方向的恒定速度的偏差。

// 下一步是创建输出。这与 step-33 中的做法类似：我们让上面定义的后处理器控制大部分的输出，除了我们直接写的原始场。对于分析解的测试案例，我们还对分析解进行了另一次投影，并打印出该场和数值解之间的差异。一旦我们定义了所有要写的量，我们就建立输出的补丁。与 step-65 类似，我们通过设置适当的标志来创建一个高阶VTK输出，这使我们能够可视化高多项式度的场。最后，我们调用 `DataOutInterface::write_vtu_in_parallel()` 函数，将结果写入给定的文件名。这个函数使用了特殊的MPI并行写设施，与其他大多数教程程序中使用的标准库的 `std::ofstream` 变体相比，它通常对并行文件系统更加优化。`write_vtu_in_parallel()`函数的一个特别好的特点是，它可以将所有MPI行列的输出合并到一个文件中，使得没有必要有一个所有此类文件的中央记录（即 "pvtu "文件）。

// 对于并行程序来说，看一下单元在处理器之间的划分往往是有启发的。为此，我们可以向 DataOut::add_data_vector() 传递一个数字向量，其中包含与当前处理器拥有的活动单元一样多的条目；然后这些数字应该是拥有这些单元的处理器的等级。例如，这样一个向量可以从 GridTools::get_subdomain_association(). 中获得。另一方面，在每个MPI进程中，DataOut将只读取那些对应于本地拥有的单元的条目，这些条目当然都有相同的值：即当前进程的等级。矢量的其余条目中的内容实际上并不重要，因此我们可以用一个廉价的技巧逃脱。我们只是把我们给 DataOut::add_data_vector() 的向量的所有*值都填上当前MPI进程的等级。关键是在每个进程中，只有对应于本地拥有的单元格的条目会被读取，而忽略其他条目中的（错误）值。事实上，每个进程提交的向量中的条目子集是正确的，这就足够了。

  template <int dim> 
  void EulerProblem<dim>::output_results(const unsigned int result_number) 
  { 
    const std::array<double, 3> errors = 
      euler_operator.compute_errors(ExactSolution<dim>(time), solution); 
    const std::string quantity_name = testcase == 0 ? "error" : "norm"; 

    pcout << "Time:" << std::setw(8) << std::setprecision(3) << time 
          << ", dt: " << std::setw(8) << std::setprecision(2) << time_step 
          << ", " << quantity_name << " rho: " << std::setprecision(4) 
          << std::setw(10) << errors[0] << ", rho * u: " << std::setprecision(4) 
          << std::setw(10) << errors[1] << ", energy:" << std::setprecision(4) 
          << std::setw(10) << errors[2] << std::endl; 

    { 
      TimerOutput::Scope t(timer, "output"); 

      Postprocessor postprocessor; 
      DataOut<dim>  data_out; 

      DataOutBase::VtkFlags flags; 
      flags.write_higher_order_cells = true; 
      data_out.set_flags(flags); 

      data_out.attach_dof_handler(dof_handler); 
      { 
        std::vector<std::string> names; 
        names.emplace_back("density"); 
        for (unsigned int d = 0; d < dim; ++d) 
          names.emplace_back("momentum"); 
        names.emplace_back("energy"); 

        std::vector<DataComponentInterpretation::DataComponentInterpretation> 
          interpretation; 
        interpretation.push_back( 
          DataComponentInterpretation::component_is_scalar); 
        for (unsigned int d = 0; d < dim; ++d) 
          interpretation.push_back( 
            DataComponentInterpretation::component_is_part_of_vector); 
        interpretation.push_back( 
          DataComponentInterpretation::component_is_scalar); 

        data_out.add_data_vector(dof_handler, solution, names, interpretation); 
      } 
      data_out.add_data_vector(solution, postprocessor); 

      LinearAlgebra::distributed::Vector<Number> reference; 
      if (testcase == 0 && dim == 2) 
        { 
          reference.reinit(solution); 
          euler_operator.project(ExactSolution<dim>(time), reference); 
          reference.sadd(-1., 1, solution); 
          std::vector<std::string> names; 
          names.emplace_back("error_density"); 
          for (unsigned int d = 0; d < dim; ++d) 
            names.emplace_back("error_momentum"); 
          names.emplace_back("error_energy"); 

          std::vector<DataComponentInterpretation::DataComponentInterpretation> 
            interpretation; 
          interpretation.push_back( 
            DataComponentInterpretation::component_is_scalar); 
          for (unsigned int d = 0; d < dim; ++d) 
            interpretation.push_back( 
              DataComponentInterpretation::component_is_part_of_vector); 
          interpretation.push_back( 
            DataComponentInterpretation::component_is_scalar); 

          data_out.add_data_vector(dof_handler, 
                                   reference, 
                                   names, 
                                   interpretation); 
        } 

      Vector<double> mpi_owner(triangulation.n_active_cells()); 
      mpi_owner = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD); 
      data_out.add_data_vector(mpi_owner, "owner"); 

      data_out.build_patches(mapping, 
                             fe.degree, 
                             DataOut<dim>::curved_inner_cells); 

      const std::string filename = 
        "solution_" + Utilities::int_to_string(result_number, 3) + ".vtu"; 
      data_out.write_vtu_in_parallel(filename, MPI_COMM_WORLD); 
    } 
  } 

//  EulerProblem::run() 函数将所有的部分组合起来。它首先调用创建网格和设置数据结构的函数，然后初始化时间积分器和低存储积分器的两个临时向量。我们称这些向量为`rk_register_1`和`rk_register_2`，并使用第一个向量表示 $\mathbf{r}_i$ ，第二个向量表示 $\mathbf{k}_i$ ，在介绍中概述的Runge--Kutta方案的公式。在我们开始时间循环之前，我们通过 `EulerOperator::compute_cell_transport_speed()` 函数计算时间步长。为了便于比较，我们将那里得到的结果与最小网格尺寸进行比较，并将它们打印到屏幕上。对于像本教程程序中接近于统一的声速和速度，预测的有效网格尺寸将是接近的，但如果缩放比例不同，它们可能会有变化。

  template <int dim> 
  void EulerProblem<dim>::run() 
  { 
    { 
      const unsigned int n_vect_number = VectorizedArray<Number>::size(); 
      const unsigned int n_vect_bits   = 8 * sizeof(Number) * n_vect_number; 

      pcout << "Running with " 
            << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) 
            << " MPI processes" << std::endl; 
      pcout << "Vectorization over " << n_vect_number << " " 
            << (std::is_same<Number, double>::value ? "doubles" : "floats") 
            << " = " << n_vect_bits << " bits (" 
            << Utilities::System::get_current_vectorization_level() << ")" 
            << std::endl; 
    } 

    make_grid_and_dofs(); 

    const LowStorageRungeKuttaIntegrator integrator(lsrk_scheme); 

    LinearAlgebra::distributed::Vector<Number> rk_register_1; 
    LinearAlgebra::distributed::Vector<Number> rk_register_2; 
    rk_register_1.reinit(solution); 
    rk_register_2.reinit(solution); 

    euler_operator.project(ExactSolution<dim>(time), solution); 

    double min_vertex_distance = std::numeric_limits<double>::max(); 
    for (const auto &cell : triangulation.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        min_vertex_distance = 
          std::min(min_vertex_distance, cell->minimum_vertex_distance()); 
    min_vertex_distance = 
      Utilities::MPI::min(min_vertex_distance, MPI_COMM_WORLD); 

    time_step = courant_number * integrator.n_stages() / 
                euler_operator.compute_cell_transport_speed(solution); 
    pcout << "Time step size: " << time_step 
          << ", minimal h: " << min_vertex_distance 
          << ", initial transport scaling: " 
          << 1. / euler_operator.compute_cell_transport_speed(solution) 
          << std::endl 
          << std::endl; 

    output_results(0); 

// 现在我们准备开始时间循环，我们一直运行到时间达到预期的结束时间。每隔5个时间步长，我们就计算一个新的时间步长估计值--由于解决方案是非线性的，在模拟过程中调整这个值是最有效的。如果Courant数选择得过于激进，模拟通常会在时间步数为NaN时爆炸，所以在这里很容易发现。有一点需要注意的是，由于不同的时间步长选择的相互作用，四舍五入的误差可能会传播到前几位数，从而导致略有不同的解决方案。为了降低这种敏感性，通常的做法是将时间步长四舍五入或截断到几位数，例如在这种情况下是3。如果当前时间接近规定的输出 "刻度 "值（如0.02），我们也会写出输出。在时间循环结束后，我们通过打印一些统计数据来总结计算，这主要由 TimerOutput::print_wall_time_statistics() 函数完成。

    unsigned int timestep_number = 0; 

    while (time < final_time - 1e-12) 
      { 
        ++timestep_number; 
        if (timestep_number % 5 == 0) 
          time_step = 
            courant_number * integrator.n_stages() / 
            Utilities::truncate_to_n_digits( 
              euler_operator.compute_cell_transport_speed(solution), 3); 

        { 
          TimerOutput::Scope t(timer, "rk time stepping total"); 
          integrator.perform_time_step(euler_operator, 
                                       time, 
                                       time_step, 
                                       solution, 
                                       rk_register_1, 
                                       rk_register_2); 
        } 

        time += time_step; 

        if (static_cast<int>(time / output_tick) != 
              static_cast<int>((time - time_step) / output_tick) || 
            time >= final_time - 1e-12) 
          output_results( 
            static_cast<unsigned int>(std::round(time / output_tick))); 
      } 

    timer.print_wall_time_statistics(MPI_COMM_WORLD); 
    pcout << std::endl; 
  } 

} // namespace Euler_DG 

// main()函数并不令人惊讶，它遵循了以前所有MPI程序中的做法。当我们运行一个MPI程序时，我们需要调用`MPI_Init()`和`MPI_Finalize()`，我们通过 Utilities::MPI::MPI_InitFinalize 数据结构来完成。请注意，我们只用MPI来运行程序，并将线程数设置为1。

int main(int argc, char **argv) 
{ 
  using namespace Euler_DG; 
  using namespace dealii; 

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 

  try 
    { 
      deallog.depth_console(0); 

      EulerProblem<dimension> euler_problem; 
      euler_problem.run(); 
    } 
  catch (std::exception &exc) 
    { 
      std::cerr << std::endl 
                << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      std::cerr << "Exception on processing: " << std::endl 
                << exc.what() << std::endl 
                << "Aborting!" << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 

      return 1; 
    } 
  catch (...) 
    { 
      std::cerr << std::endl 
                << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      std::cerr << "Unknown exception!" << std::endl 
                << "Aborting!" << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      return 1; 
    } 

  return 0; 
} 


