

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2012 - 2021 by the deal.II authors 
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
 * Authors: Joerg Frohne, Texas A&M University and 
 *                        University of Siegen, 2012, 2013 
 *          Wolfgang Bangerth, Texas A&M University, 2012, 2013 
 *          Timo Heister, Texas A&M University, 2013 
 */ 


// @sect3{Include files}  这组包含文件在这个时候已经没有什么惊喜了。

#include <deal.II/base/conditional_ostream.h> 
#include <deal.II/base/parameter_handler.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/base/index_set.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/timer.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparsity_tools.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/block_sparsity_pattern.h> 
#include <deal.II/lac/solver_bicgstab.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/trilinos_sparse_matrix.h> 
#include <deal.II/lac/trilinos_block_sparse_matrix.h> 
#include <deal.II/lac/trilinos_vector.h> 
#include <deal.II/lac/trilinos_parallel_block_vector.h> 
#include <deal.II/lac/trilinos_precondition.h> 
#include <deal.II/lac/trilinos_solver.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_tools.h> 
#include <deal.II/grid/manifold_lib.h> 

#include <deal.II/distributed/tria.h> 
#include <deal.II/distributed/grid_refinement.h> 
#include <deal.II/distributed/solution_transfer.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 
#include <deal.II/numerics/fe_field_function.h> 

#include <fstream> 
#include <iostream> 

// 最后，我们包括两个系统头文件，让我们为输出文件创建一个目录。第一个头文件提供了 <code>mkdir</code> 的功能，第二个头文件让我们确定在 <code>mkdir</code> 失败时发生了什么。

#include <sys/stat.h> 
#include <cerrno> 

namespace Step42 
{ 
  using namespace dealii; 
// @sect3{The <code>ConstitutiveLaw</code> class template}  

// 该类提供了一个构成法的接口，即应变  $\varepsilon(\mathbf u)$  和应力  $\sigma$  之间的关系。在这个例子中，我们使用的是具有线性、各向同性硬化的弹塑性材料行为。这种材料的特点是杨氏模量  $E$  ，泊松比  $\nu$  ，初始屈服应力  $\sigma_0$  和各向同性硬化参数  $\gamma$  。 对于 $\gamma = 0$ ，我们得到完美的弹塑性行为。

// 正如描述这个程序的论文所解释的那样，第一个牛顿步骤是用一个完全弹性材料模型来解决的，以避免同时处理两种非线性（塑性和接触）。为此，这个类有一个函数 <code>set_sigma_0()</code> ，我们在后面使用这个函数，简单地将 $\sigma_0$ 设置为一个非常大的值--基本上保证了实际应力不会超过它，从而产生一个弹性材料。当我们准备使用塑性模型时，我们使用相同的函数将 $\sigma_0$ 设置回其适当的值。 由于这种方法，我们需要将 <code>sigma_0</code> 作为这个类的唯一非静态成员变量。

  template <int dim> 
  class ConstitutiveLaw 
  { 
  public: 
    ConstitutiveLaw(const double E, 
                    const double nu, 
                    const double sigma_0, 
                    const double gamma); 

    void set_sigma_0(double sigma_zero); 

    bool get_stress_strain_tensor( 
      const SymmetricTensor<2, dim> &strain_tensor, 
      SymmetricTensor<4, dim> &      stress_strain_tensor) const; 

    void get_linearized_stress_strain_tensors( 
      const SymmetricTensor<2, dim> &strain_tensor, 
      SymmetricTensor<4, dim> &      stress_strain_tensor_linearized, 
      SymmetricTensor<4, dim> &      stress_strain_tensor) const; 

  private: 
    const double kappa; 
    const double mu; 
    double       sigma_0; 
    const double gamma; 

    const SymmetricTensor<4, dim> stress_strain_tensor_kappa; 
    const SymmetricTensor<4, dim> stress_strain_tensor_mu; 
  }; 

// ConstitutiveLaw类的构造函数为我们的可变形体设置所需的材料参数。弹性各向同性介质的材料参数可以用多种方式定义，如一对 $E, \nu$ （弹性模量和泊松数），使用Lam&eacute;参数 $\lambda,mu$ 或其他几种常用的约定。在这里，构造器采用 $E,\nu$ 形式的材料参数描述，但由于这证明这些不是出现在塑性投影仪方程中的系数，我们立即将它们转换为更合适的体模和剪模集合 $\kappa,\mu$ 。 此外，构造器以 $\sigma_0$ （无任何塑性应变的屈服应力）和 $\gamma$ （硬化参数）作为参数。在这个构造函数中，我们还计算了应力-应变关系的两个主成分及其线性化。

  template <int dim> 
  ConstitutiveLaw<dim>::ConstitutiveLaw(double E, 
                                        double nu, 
                                        double sigma_0, 
                                        double gamma) 
    : kappa(E / (3 * (1 - 2 * nu))) 
    , mu(E / (2 * (1 + nu))) 
    , sigma_0(sigma_0) 
    , gamma(gamma) 
    , stress_strain_tensor_kappa(kappa * 
                                 outer_product(unit_symmetric_tensor<dim>(), 
                                               unit_symmetric_tensor<dim>())) 
    , stress_strain_tensor_mu( 
        2 * mu * 
        (identity_tensor<dim>() - outer_product(unit_symmetric_tensor<dim>(), 
                                                unit_symmetric_tensor<dim>()) / 
                                    3.0)) 
  {} 

  template <int dim> 
  void ConstitutiveLaw<dim>::set_sigma_0(double sigma_zero) 
  { 
    sigma_0 = sigma_zero; 
  } 
// @sect4{ConstitutiveLaw::get_stress_strain_tensor}  

// 这是构成法则的主成分。它计算的是四阶对称张量，根据上面给出的投影，当在一个特定的应变点上评估时，该张量将应变与应力联系起来。我们需要这个函数来计算 <code>PlasticityContactProblem::residual_nl_system()</code> 中的非线性残差，我们将这个张量与正交点的应变相乘。计算遵循介绍中列出的公式。在比较那里的公式和下面的实现时，记得 $C_\mu : \varepsilon = \tau_D$ 和 $C_\kappa : \varepsilon = \kappa \text{trace}(\varepsilon) I = \frac 13 \text{trace}(\tau) I$  。

// 该函数返回正交点是否是塑性的，以便在下游对有多少正交点是塑性的，有多少是弹性的进行一些统计。

  template <int dim> 
  bool ConstitutiveLaw<dim>::get_stress_strain_tensor( 
    const SymmetricTensor<2, dim> &strain_tensor, 
    SymmetricTensor<4, dim> &      stress_strain_tensor) const 
  { 
    Assert(dim == 3, ExcNotImplemented()); 

    SymmetricTensor<2, dim> stress_tensor; 
    stress_tensor = 
      (stress_strain_tensor_kappa + stress_strain_tensor_mu) * strain_tensor; 

    const SymmetricTensor<2, dim> deviator_stress_tensor = 
      deviator(stress_tensor); 
    const double deviator_stress_tensor_norm = deviator_stress_tensor.norm(); 

    stress_strain_tensor = stress_strain_tensor_mu; 
    if (deviator_stress_tensor_norm > sigma_0) 
      { 
        const double beta = sigma_0 / deviator_stress_tensor_norm; 
        stress_strain_tensor *= (gamma + (1 - gamma) * beta); 
      } 

    stress_strain_tensor += stress_strain_tensor_kappa; 

    return (deviator_stress_tensor_norm > sigma_0); 
  } 
// @sect4{ConstitutiveLaw::get_linearized_stress_strain_tensors}  

// 该函数返回线性化的应力应变张量，围绕前一个牛顿步骤 $u^{i-1}$ 的解进行线性化  $i-1$  。 参数 <code>strain_tensor</code> （通常表示为 $\varepsilon(u^{i-1})$ ）必须作为参数传递，并作为线性化点。该函数在变量stress_strain_tensor中返回非线性构成法的导数，在stress_strain_tensor_linearized中返回线性化问题的应力-应变张量。 参见 PlasticityContactProblem::assemble_nl_system ，其中使用了这个函数。

  template <int dim> 
  void ConstitutiveLaw<dim>::get_linearized_stress_strain_tensors( 
    const SymmetricTensor<2, dim> &strain_tensor, 
    SymmetricTensor<4, dim> &      stress_strain_tensor_linearized, 
    SymmetricTensor<4, dim> &      stress_strain_tensor) const 
  { 
    Assert(dim == 3, ExcNotImplemented()); 

    SymmetricTensor<2, dim> stress_tensor; 
    stress_tensor = 
      (stress_strain_tensor_kappa + stress_strain_tensor_mu) * strain_tensor; 

    stress_strain_tensor            = stress_strain_tensor_mu; 
    stress_strain_tensor_linearized = stress_strain_tensor_mu; 

    SymmetricTensor<2, dim> deviator_stress_tensor = deviator(stress_tensor); 
    const double deviator_stress_tensor_norm = deviator_stress_tensor.norm(); 

    if (deviator_stress_tensor_norm > sigma_0) 
      { 
        const double beta = sigma_0 / deviator_stress_tensor_norm; 
        stress_strain_tensor *= (gamma + (1 - gamma) * beta); 
        stress_strain_tensor_linearized *= (gamma + (1 - gamma) * beta); 
        deviator_stress_tensor /= deviator_stress_tensor_norm; 
        stress_strain_tensor_linearized -= 
          (1 - gamma) * beta * 2 * mu * 
          outer_product(deviator_stress_tensor, deviator_stress_tensor); 
      } 

    stress_strain_tensor += stress_strain_tensor_kappa; 
    stress_strain_tensor_linearized += stress_strain_tensor_kappa; 
  } 
//<h3>Equation data: boundary forces, boundary values, obstacles</h3>

// 下面的内容应该是比较标准的。我们需要边界强迫项（我们在此选择为零）和不属于接触面的边界部分的边界值（在此也选择为零）的类。

  namespace EquationData 
  { 
    template <int dim> 
    class BoundaryForce : public Function<dim> 
    { 
    public: 
      BoundaryForce(); 

      virtual double value(const Point<dim> & p, 
                           const unsigned int component = 0) const override; 

      virtual void vector_value(const Point<dim> &p, 
                                Vector<double> &  values) const override; 
    }; 

    template <int dim> 
    BoundaryForce<dim>::BoundaryForce() 
      : Function<dim>(dim) 
    {} 

    template <int dim> 
    double BoundaryForce<dim>::value(const Point<dim> &, 
                                     const unsigned int) const 
    { 
      return 0.; 
    } 

    template <int dim> 
    void BoundaryForce<dim>::vector_value(const Point<dim> &p, 
                                          Vector<double> &  values) const 
    { 
      for (unsigned int c = 0; c < this->n_components; ++c) 
        values(c) = BoundaryForce<dim>::value(p, c); 
    } 

    template <int dim> 
    class BoundaryValues : public Function<dim> 
    { 
    public: 
      BoundaryValues(); 

      virtual double value(const Point<dim> & p, 
                           const unsigned int component = 0) const override; 
    }; 

    template <int dim> 
    BoundaryValues<dim>::BoundaryValues() 
      : Function<dim>(dim) 
    {} 

    template <int dim> 
    double BoundaryValues<dim>::value(const Point<dim> &, 
                                      const unsigned int) const 
    { 
      return 0.; 
    } 

//  @sect4{The <code>SphereObstacle</code> class}  

// 下面这个类是可以从输入文件中选择的两个障碍物中的第一个。它描述了一个以位置 $x=y=0.5, z=z_{\text{surface}}+0.59$ 和半径 $r=0.6$ 为中心的球体，其中 $z_{\text{surface}}$ 是可变形体的（平）表面的垂直位置。该函数的 <code>value</code> 返回给定 $x,y$ 值的障碍物位置，如果该点实际位于球体下方，则返回一个不可能干扰变形的大正值，如果它位于球体的 "阴影 "之外。

    template <int dim> 
    class SphereObstacle : public Function<dim> 
    { 
    public: 
      SphereObstacle(const double z_surface); 

      virtual double value(const Point<dim> & p, 
                           const unsigned int component = 0) const override; 

      virtual void vector_value(const Point<dim> &p, 
                                Vector<double> &  values) const override; 

    private: 
      const double z_surface; 
    }; 

    template <int dim> 
    SphereObstacle<dim>::SphereObstacle(const double z_surface) 
      : Function<dim>(dim) 
      , z_surface(z_surface) 
    {} 

    template <int dim> 
    double SphereObstacle<dim>::value(const Point<dim> & p, 
                                      const unsigned int component) const 
    { 
      if (component == 0) 
        return p(0); 
      else if (component == 1) 
        return p(1); 
      else if (component == 2) 
        { 
          if ((p(0) - 0.5) * (p(0) - 0.5) + (p(1) - 0.5) * (p(1) - 0.5) < 0.36) 
            return (-std::sqrt(0.36 - (p(0) - 0.5) * (p(0) - 0.5) - 
                               (p(1) - 0.5) * (p(1) - 0.5)) + 
                    z_surface + 0.59); 
          else 
            return 1000; 
        } 

      Assert(false, ExcNotImplemented()); 
      return 1e9; // an unreasonable value; ignored in debug mode because of the 

// 前面的断言

    } 

    template <int dim> 
    void SphereObstacle<dim>::vector_value(const Point<dim> &p, 
                                           Vector<double> &  values) const 
    { 
      for (unsigned int c = 0; c < this->n_components; ++c) 
        values(c) = SphereObstacle<dim>::value(p, c); 
    } 
// @sect4{The <code>BitmapFile</code> and <code>ChineseObstacle</code> classes}  

// 下面两个类描述了介绍中概述的障碍物，即汉字。两个中的第一个， <code>BitmapFile</code> 负责从一个以pbm ascii格式存储的图片文件中读入数据。这个数据将被双线性插值，从而提供一个描述障碍物的函数。(下面的代码显示了如何通过在给定的数据点之间进行内插来构造一个函数。人们可以使用在这个教程程序写完后引入的 Functions::InterpolatedUniformGridData, ，它正是我们在这里想要的，但看看如何手工操作是有启发的）。)

// 我们从文件中读取的数据将被存储在一个名为 obstacle_data 的双 std::vector 中。 这个向量构成了计算单片双线性函数的基础，作为一个多项式插值。我们将从文件中读取的数据由零（白色）和一（黑色）组成。

//  <code>hx,hy</code> 变量表示 $x$ 和 $y$ 方向的像素之间的间距。  <code>nx,ny</code> 是这些方向上的像素的数量。   <code>get_value()</code> 返回图像在给定位置的值，由相邻像素值插值而成。

    template <int dim> 
    class BitmapFile 
    { 
    public: 
      BitmapFile(const std::string &name); 

      double get_value(const double x, const double y) const; 

    private: 
      std::vector<double> obstacle_data; 
      double              hx, hy; 
      int                 nx, ny; 

      double get_pixel_value(const int i, const int j) const; 
    }; 

// 该类的构造函数从给定的文件名中读入描述障碍物的数据。

    template <int dim> 
    BitmapFile<dim>::BitmapFile(const std::string &name) 
      : obstacle_data(0) 
      , hx(0) 
      , hy(0) 
      , nx(0) 
      , ny(0) 
    { 
      std::ifstream f(name); 
      AssertThrow(f, 
                  ExcMessage(std::string("Can't read from file <") + name + 
                             ">!")); 

      std::string temp; 
      f >> temp >> nx >> ny; 

      AssertThrow(nx > 0 && ny > 0, ExcMessage("Invalid file format.")); 

      for (int k = 0; k < nx * ny; ++k) 
        { 
          double val; 
          f >> val; 
          obstacle_data.push_back(val); 
        } 

      hx = 1.0 / (nx - 1); 
      hy = 1.0 / (ny - 1); 

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
        std::cout << "Read obstacle from file <" << name << ">" << std::endl 
                  << "Resolution of the scanned obstacle picture: " << nx 
                  << " x " << ny << std::endl; 
    } 

// 下面两个函数返回坐标为 $i,j$ 的给定像素的值，我们将其与定义在位置 <code>i*hx, j*hy</code> 的函数值和任意坐标 $x,y$ 的函数值相识别，在这里我们对两个函数中第一个函数返回的点值进行双线性内插。在第二个函数中，对于每个 $x,y$ ，我们首先计算离 $x,y$ 左下方最近的像素坐标的（整数）位置，然后计算这个像素内的坐标 $\xi,\eta$ 。我们从下方和上方截断这两种变量，以避免在评估函数时超出其定义的范围而可能发生的舍入误差问题。

    template <int dim> 
    double BitmapFile<dim>::get_pixel_value(const int i, const int j) const 
    { 
      assert(i >= 0 && i < nx); 
      assert(j >= 0 && j < ny); 
      return obstacle_data[nx * (ny - 1 - j) + i]; 
    } 

    template <int dim> 
    double BitmapFile<dim>::get_value(const double x, const double y) const 
    { 
      const int ix = std::min(std::max(static_cast<int>(x / hx), 0), nx - 2); 
      const int iy = std::min(std::max(static_cast<int>(y / hy), 0), ny - 2); 

      const double xi  = std::min(std::max((x - ix * hx) / hx, 1.), 0.); 
      const double eta = std::min(std::max((y - iy * hy) / hy, 1.), 0.); 

      return ((1 - xi) * (1 - eta) * get_pixel_value(ix, iy) + 
              xi * (1 - eta) * get_pixel_value(ix + 1, iy) + 
              (1 - xi) * eta * get_pixel_value(ix, iy + 1) + 
              xi * eta * get_pixel_value(ix + 1, iy + 1)); 
    } 

// 最后，这是一个实际使用上面的类的类。它有一个BitmapFile对象作为成员，描述障碍物的高度。如上所述，BitmapFile类将为我们提供一个掩码，即要么是0，要么是1的值（如果你要求的是像素之间的位置，则是在0和1之间插值的值）。这个类将其转化为高度，即低于可变形体表面的0.001（如果BitmapFile类在此位置报告为1）或高于障碍物的0.999（如果BitmapFile类报告为0）。那么下面的函数应该是不言自明的。

    template <int dim> 
    class ChineseObstacle : public Function<dim> 
    { 
    public: 
      ChineseObstacle(const std::string &filename, const double z_surface); 

      virtual double value(const Point<dim> & p, 
                           const unsigned int component = 0) const override; 

      virtual void vector_value(const Point<dim> &p, 
                                Vector<double> &  values) const override; 

    private: 
      const BitmapFile<dim> input_obstacle; 
      double                z_surface; 
    }; 

    template <int dim> 
    ChineseObstacle<dim>::ChineseObstacle(const std::string &filename, 
                                          const double       z_surface) 
      : Function<dim>(dim) 
      , input_obstacle(filename) 
      , z_surface(z_surface) 
    {} 

    template <int dim> 
    double ChineseObstacle<dim>::value(const Point<dim> & p, 
                                       const unsigned int component) const 
    { 
      if (component == 0) 
        return p(0); 
      if (component == 1) 
        return p(1); 
      else if (component == 2) 
        { 
          if (p(0) >= 0.0 && p(0) <= 1.0 && p(1) >= 0.0 && p(1) <= 1.0) 
            return z_surface + 0.999 - input_obstacle.get_value(p(0), p(1)); 
        } 

      Assert(false, ExcNotImplemented()); 
      return 1e9; // an unreasonable value; ignored in debug mode because of the 

// 前面的断言

    } 

    template <int dim> 
    void ChineseObstacle<dim>::vector_value(const Point<dim> &p, 
                                            Vector<double> &  values) const 
    { 
      for (unsigned int c = 0; c < this->n_components; ++c) 
        values(c) = ChineseObstacle<dim>::value(p, c); 
    } 
  } // namespace EquationData 
// @sect3{The <code>PlasticityContactProblem</code> class template}  

// 这是本程序的主类，提供了描述非线性接触问题所需的所有函数和变量。它接近于 step-41 ，但有一些额外的功能，如处理悬挂节点，牛顿方法，使用Trilinos和p4est进行并行分布式计算。处理悬空节点使生活变得有点复杂，因为我们现在需要另一个AffineConstraints对象。我们为接触情况下的主动集合方法创建一个牛顿方法，并处理构成法的非线性算子。

// 这个类的总体布局与其他大多数教程程序非常相似。为了使我们的生活更容易一些，这个类从输入文件中读取一组输入参数。这些参数，使用ParameterHandler类，在 <code>declare_parameters</code> 函数中声明（该函数是静态的，因此它可以在我们创建当前类型的对象之前被调用），然后一个已经用于读取输入文件的ParameterHandler对象将被传递给该类的构造函数。

// 其余的成员函数大体上与我们在其他几个教程程序中看到的一样，虽然为当前的非线性系统增加了一些内容。我们将在下文中对它们的用途进行评论。

  template <int dim> 
  class PlasticityContactProblem 
  { 
  public: 
    PlasticityContactProblem(const ParameterHandler &prm); 

    void run(); 

    static void declare_parameters(ParameterHandler &prm); 

  private: 
    void make_grid(); 
    void setup_system(); 
    void compute_dirichlet_constraints(); 
    void update_solution_and_constraints(); 
    void 
         assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &mass_matrix); 
    void assemble_newton_system( 
      const TrilinosWrappers::MPI::Vector &linearization_point); 
    void compute_nonlinear_residual( 
      const TrilinosWrappers::MPI::Vector &linearization_point); 
    void solve_newton_system(); 
    void solve_newton(); 
    void refine_grid(); 
    void move_mesh(const TrilinosWrappers::MPI::Vector &displacement) const; 
    void output_results(const unsigned int current_refinement_cycle); 

    void output_contact_force() const; 

// 就成员变量而言，我们先用一个变量来表示这个程序运行的MPI宇宙，一个我们用来让确切的一个处理器产生输出到控制台的流（见 step-17  ）和一个用来为程序的各个部分计时的变量。

    MPI_Comm           mpi_communicator; 
    ConditionalOStream pcout; 
    TimerOutput        computing_timer; 

// 下一组描述网格和有限元空间。特别是，对于这个并行程序，有限元空间有与之相关的变量，表明哪些自由度存在于当前的处理器上（索引集，也见 step-40 和 @ref distributed 文档模块），以及各种约束：那些由悬挂节点，由Dirichlet边界条件，以及由接触节点的活动集施加的约束。在这里定义的三个AffineConstraints变量中，第一个变量只包含悬挂节点的约束，第二个变量也包含与Dirichlet边界条件相关的约束，第三个变量包含这些约束和接触约束。

// 变量 <code>active_set</code> 包括那些由接触约束的自由度，我们用 <code>fraction_of_plastic_q_points_per_cell</code> 来跟踪每个单元上应力等于屈服应力的正交点的分数。后者仅用于创建显示塑性区的图形输出，但不用于任何进一步的计算；该变量是该类的成员变量，因为该信息是作为计算残差的副产品计算的，但仅在很晚的时候使用。(注意，该向量是一个长度等于<i>local mesh</i>上活动单元数量的向量；它从未被用来在处理器之间交换信息，因此可以是一个普通的deal.II向量)。

    const unsigned int                        n_initial_global_refinements; 
    parallel::distributed::Triangulation<dim> triangulation; 

    const unsigned int fe_degree; 
    FESystem<dim>      fe; 
    DoFHandler<dim>    dof_handler; 

    IndexSet locally_owned_dofs; 
    IndexSet locally_relevant_dofs; 

    AffineConstraints<double> constraints_hanging_nodes; 
    AffineConstraints<double> constraints_dirichlet_and_hanging_nodes; 
    AffineConstraints<double> all_constraints; 

    IndexSet      active_set; 
    Vector<float> fraction_of_plastic_q_points_per_cell; 

// 下一个变量块对应的是解决方案和我们需要形成的线性系统。特别是，这包括牛顿矩阵和右手边；与残差（即牛顿右手边）相对应的向量，但我们没有消除其中的各种约束，该向量用于确定在下一次迭代中需要约束哪些自由度；以及一个与介绍中简要提到的 $B$ 矩阵的对角线相对应的向量，并在随文中讨论。

    TrilinosWrappers::SparseMatrix newton_matrix; 

    TrilinosWrappers::MPI::Vector solution; 
    TrilinosWrappers::MPI::Vector newton_rhs; 
    TrilinosWrappers::MPI::Vector newton_rhs_uncondensed; 
    TrilinosWrappers::MPI::Vector diag_mass_matrix_vector; 

// 下一个块包含描述材料响应的变量。

    const double         e_modulus, nu, gamma, sigma_0; 
    ConstitutiveLaw<dim> constitutive_law; 

// 然后是各种各样的其他变量，用于识别参数文件所选择的要求我们建立的网格，被推入可变形体的障碍物，网格细化策略，是否将解决方案从一个网格转移到下一个网格，以及要执行多少个网格细化循环。在可能的情况下，我们将这些类型的变量标记为 <code>const</code> ，以帮助读者识别哪些变量以后可能会被修改，哪些可能不会被修改（输出目录是一个例外--它在构造函数之外从不被修改，但在构造函数中冒号后面的成员初始化列表中初始化是很尴尬的，因为在那里我们只有一次机会设置它；网格细化准则也是如此）。

    const std::string                          base_mesh; 
    const std::shared_ptr<const Function<dim>> obstacle; 

    struct RefinementStrategy 
    { 
      enum value 
      { 
        refine_global, 
        refine_percentage, 
        refine_fix_dofs 
      }; 
    }; 
    typename RefinementStrategy::value refinement_strategy; 

    const bool         transfer_solution; 
    std::string        output_dir; 
    const unsigned int n_refinement_cycles; 
    unsigned int       current_refinement_cycle; 
  }; 
// @sect3{Implementation of the <code>PlasticityContactProblem</code> class}  
// @sect4{PlasticityContactProblem::declare_parameters}  

// 让我们从声明可在输入文件中选择的运行时参数开始。这些值将在本类的构造函数中读回，以初始化本类的成员变量。

  template <int dim> 
  void PlasticityContactProblem<dim>::declare_parameters(ParameterHandler &prm) 
  { 
    prm.declare_entry( 
      "polynomial degree", 
      "1", 
      Patterns::Integer(), 
      "Polynomial degree of the FE_Q finite element space, typically 1 or 2."); 
    prm.declare_entry("number of initial refinements", 
                      "2", 
                      Patterns::Integer(), 
                      "Number of initial global mesh refinement steps before " 
                      "the first computation."); 
    prm.declare_entry( 
      "refinement strategy", 
      "percentage", 
      Patterns::Selection("global|percentage"), 
      "Mesh refinement strategy:\n" 
      " global: one global refinement\n" 
      " percentage: a fixed percentage of cells gets refined using the Kelly estimator."); 
    prm.declare_entry("number of cycles", 
                      "5", 
                      Patterns::Integer(), 
                      "Number of adaptive mesh refinement cycles to run."); 
    prm.declare_entry( 
      "obstacle", 
      "sphere", 
      Patterns::Selection("sphere|read from file"), 
      "The name of the obstacle to use. This may either be 'sphere' if we should " 
      "use a spherical obstacle, or 'read from file' in which case the obstacle " 
      "will be read from a file named 'obstacle.pbm' that is supposed to be in " 
      "ASCII PBM format."); 
    prm.declare_entry( 
      "output directory", 
      "", 
      Patterns::Anything(), 
      "Directory for output files (graphical output and benchmark " 
      "statistics). If empty, use the current directory."); 
    prm.declare_entry( 
      "transfer solution", 
      "false", 
      Patterns::Bool(), 
      "Whether the solution should be used as a starting guess " 
      "for the next finer mesh. If false, then the iteration starts at " 
      "zero on every mesh."); 
    prm.declare_entry("base mesh", 
                      "box", 
                      Patterns::Selection("box|half sphere"), 
                      "Select the shape of the domain: 'box' or 'half sphere'"); 
  } 
// @sect4{The <code>PlasticityContactProblem</code> constructor}  

// 鉴于成员变量的声明以及从输入文件中读取的运行时参数的声明，在这个构造函数中没有任何令人惊讶的地方。在正文中，我们初始化了网格细化策略和输出目录，必要时创建这样一个目录。

  template <int dim> 
  PlasticityContactProblem<dim>::PlasticityContactProblem( 
    const ParameterHandler &prm) 
    : mpi_communicator(MPI_COMM_WORLD) 
    , pcout(std::cout, 
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)) 
    , computing_timer(MPI_COMM_WORLD, 
                      pcout, 
                      TimerOutput::never, 
                      TimerOutput::wall_times) 

    , n_initial_global_refinements( 
        prm.get_integer("number of initial refinements")) 
    , triangulation(mpi_communicator) 
    , fe_degree(prm.get_integer("polynomial degree")) 
    , fe(FE_Q<dim>(QGaussLobatto<1>(fe_degree + 1)), dim) 
    , dof_handler(triangulation) 

    , e_modulus(200000) 
    , nu(0.3) 
    , gamma(0.01) 
    , sigma_0(400.0) 
    , constitutive_law(e_modulus, nu, sigma_0, gamma) 

    , base_mesh(prm.get("base mesh")) 
    , obstacle(prm.get("obstacle") == "read from file" ? 
                 static_cast<const Function<dim> *>( 
                   new EquationData::ChineseObstacle<dim>( 
                     "obstacle.pbm", 
                     (base_mesh == "box" ? 1.0 : 0.5))) : 
                 static_cast<const Function<dim> *>( 
                   new EquationData::SphereObstacle<dim>( 
                     base_mesh == "box" ? 1.0 : 0.5))) 

    , transfer_solution(prm.get_bool("transfer solution")) 
    , n_refinement_cycles(prm.get_integer("number of cycles")) 
    , current_refinement_cycle(0) 

  { 
    std::string strat = prm.get("refinement strategy"); 
    if (strat == "global") 
      refinement_strategy = RefinementStrategy::refine_global; 
    else if (strat == "percentage") 
      refinement_strategy = RefinementStrategy::refine_percentage; 
    else 
      AssertThrow(false, ExcNotImplemented()); 

    output_dir = prm.get("output directory"); 
    if (output_dir != "" && *(output_dir.rbegin()) != '/') 
      output_dir += "/"; 

// 如果有必要，为输出创建一个新的目录。

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0) 
      { 
        const int ierr = mkdir(output_dir.c_str(), 0777); 
        AssertThrow(ierr == 0 || errno == EEXIST, ExcIO()); 
      } 

    pcout << "    Using output directory '" << output_dir << "'" << std::endl; 
    pcout << "    FE degree " << fe_degree << std::endl; 
    pcout << "    transfer solution " << (transfer_solution ? "true" : "false") 
          << std::endl; 
  } 

//  @sect4{PlasticityContactProblem::make_grid}  

// 下一个区块是关于构建起始网格的。我们将使用下面的辅助函数和 <code>make_grid()</code> 的第一个块来构造一个对应于半球形的网格。deal.II有一个函数可以创建这样的网格，但是它的位置和方向都是错误的，所以我们需要在使用它之前对它进行一些位移和旋转。

// 供以后参考，如 GridGenerator::half_hyper_ball(), 文件中所述，半球体的平坦表面的边界指标为零，而其余部分的边界指标为一。

  Point<3> rotate_half_sphere(const Point<3> &in) 
  { 
    return {in(2), in(1), -in(0)}; 
  } 

  template <int dim> 
  void PlasticityContactProblem<dim>::make_grid() 
  { 
    if (base_mesh == "half sphere") 
      { 
        const Point<dim> center(0, 0, 0); 
        const double     radius = 0.8; 
        GridGenerator::half_hyper_ball(triangulation, center, radius); 

// 由于我们将在下面附加一个不同的流形，我们立即清除默认的流形描述。

        triangulation.reset_all_manifolds(); 

        GridTools::transform(&rotate_half_sphere, triangulation); 
        GridTools::shift(Point<dim>(0.5, 0.5, 0.5), triangulation); 

        SphericalManifold<dim> manifold_description(Point<dim>(0.5, 0.5, 0.5)); 
        GridTools::copy_boundary_to_manifold_id(triangulation); 
        triangulation.set_manifold(0, manifold_description); 
      } 

// 或者，创建一个超立方体网格。创建后，按如下方式分配边界指标。
// @code
//  >     _______
//  >    /  1    /|
//  >   /______ / |
//  >  |       | 8|
//  >  |   8   | /
//  >  |_______|/
//  >      6
//  @endcode
  // 换句话说，立方体的边的边界指标是8。底部的边界指标是6，顶部的指标是1。我们通过循环所有面的所有单元并查看单元中心的坐标值来设置这些指标，并在以后评估哪个边界将携带迪里希特边界条件或将受到潜在接触时使用这些指标。(在目前的情况下，网格只包含一个单元，它的所有面都在边界上，所以严格来说，所有单元的循环和查询一个面是否在边界上都是不必要的；我们保留它们只是出于习惯：这种代码可以在许多程序中找到，基本上都是这种形式。)

    else 
      { 
        const Point<dim> p1(0, 0, 0); 
        const Point<dim> p2(1.0, 1.0, 1.0); 

        GridGenerator::hyper_rectangle(triangulation, p1, p2); 

        for (const auto &cell : triangulation.active_cell_iterators()) 
          for (const auto &face : cell->face_iterators()) 
            if (face->at_boundary()) 
              { 
                if (std::fabs(face->center()[2] - p2[2]) < 1e-12) 
                  face->set_boundary_id(1); 
                if (std::fabs(face->center()[0] - p1[0]) < 1e-12 || 
                    std::fabs(face->center()[0] - p2[0]) < 1e-12 || 
                    std::fabs(face->center()[1] - p1[1]) < 1e-12 || 
                    std::fabs(face->center()[1] - p2[1]) < 1e-12) 
                  face->set_boundary_id(8); 
                if (std::fabs(face->center()[2] - p1[2]) < 1e-12) 
                  face->set_boundary_id(6); 
              } 
      } 

    triangulation.refine_global(n_initial_global_refinements); 
  } 

//  @sect4{PlasticityContactProblem::setup_system}  

// 谜题的下一块是设置DoFHandler，调整向量大小，并处理其他各种状态变量，如索引集和约束矩阵。

// 在下面的内容中，每一组操作都被放入一个大括号封闭的块中，该块的顶部声明的变量正在进行计时（ TimerOutput::Scope 变量的构造器开始计时部分，在块的末端调用的析构器再次停止计时）。

  template <int dim> 
  void PlasticityContactProblem<dim>::setup_system() 
  { 

/* 设置dofs，并为本地拥有的相关dofs获取索引集  */ 

    { 
      TimerOutput::Scope t(computing_timer, "Setup: distribute DoFs"); 
      dof_handler.distribute_dofs(fe); 

      locally_owned_dofs = dof_handler.locally_owned_dofs(); 
      locally_relevant_dofs.clear(); 
      DoFTools::extract_locally_relevant_dofs(dof_handler, 
                                              locally_relevant_dofs); 
    } 

/*设置悬挂节点和Dirichlet约束 */ 

 
    { 
      TimerOutput::Scope t(computing_timer, "Setup: constraints"); 
      constraints_hanging_nodes.reinit(locally_relevant_dofs); 
      DoFTools::make_hanging_node_constraints(dof_handler, 
                                              constraints_hanging_nodes); 
      constraints_hanging_nodes.close(); 

      pcout << "   Number of active cells: " 
            << triangulation.n_global_active_cells() << std::endl 
            << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
            << std::endl; 

      compute_dirichlet_constraints(); 
    } 

/* 初始化向量和活动集  */ 

    { 
      TimerOutput::Scope t(computing_timer, "Setup: vectors"); 
      solution.reinit(locally_relevant_dofs, mpi_communicator); 
      newton_rhs.reinit(locally_owned_dofs, mpi_communicator); 
      newton_rhs_uncondensed.reinit(locally_owned_dofs, mpi_communicator); 
      diag_mass_matrix_vector.reinit(locally_owned_dofs, mpi_communicator); 
      fraction_of_plastic_q_points_per_cell.reinit( 
        triangulation.n_active_cells()); 

      active_set.clear(); 
      active_set.set_size(dof_handler.n_dofs()); 
    } 

// 最后，我们设置了稀疏模式和矩阵。我们暂时（ab）用系统矩阵来同时建立（对角线）矩阵，用于消除与障碍物接触的自由度，但我们随后立即将牛顿矩阵设回零。

    { 
      TimerOutput::Scope                t(computing_timer, "Setup: matrix"); 
      TrilinosWrappers::SparsityPattern sp(locally_owned_dofs, 
                                           mpi_communicator); 

      DoFTools::make_sparsity_pattern(dof_handler, 
                                      sp, 
                                      constraints_dirichlet_and_hanging_nodes, 
                                      false, 
                                      Utilities::MPI::this_mpi_process( 
                                        mpi_communicator)); 
      sp.compress(); 
      newton_matrix.reinit(sp); 

      TrilinosWrappers::SparseMatrix &mass_matrix = newton_matrix; 

      assemble_mass_matrix_diagonal(mass_matrix); 

      const unsigned int start = (newton_rhs.local_range().first), 
                         end   = (newton_rhs.local_range().second); 
      for (unsigned int j = start; j < end; ++j) 
        diag_mass_matrix_vector(j) = mass_matrix.diag_element(j); 
      diag_mass_matrix_vector.compress(VectorOperation::insert); 

      mass_matrix = 0; 
    } 
  } 
// @sect4{PlasticityContactProblem::compute_dirichlet_constraints}  

// 这个函数从前面的函数中分离出来，计算与迪里切特型边界条件相关的约束，并通过与来自悬挂节点的约束合并，将其放入 <code>constraints_dirichlet_and_hanging_nodes</code> 变量。

// 正如在介绍中所阐述的，我们需要区分两种情况。

// - 如果域是一个盒子，我们将底部的位移设置为零，并允许沿侧面的Z方向的垂直运动。如 <code>make_grid()</code> 函数所示，前者对应于边界指标6，后者对应于8。

// - 如果域是一个半球形，那么我们沿边界的弯曲部分施加零位移，与边界指标0相关。

  template <int dim> 
  void PlasticityContactProblem<dim>::compute_dirichlet_constraints() 
  { 
    constraints_dirichlet_and_hanging_nodes.reinit(locally_relevant_dofs); 
    constraints_dirichlet_and_hanging_nodes.merge(constraints_hanging_nodes); 

    if (base_mesh == "box") 
      { 

//插值解决方案的所有组成部分

        VectorTools::interpolate_boundary_values( 
          dof_handler, 
          6, 
          EquationData::BoundaryValues<dim>(), 
          constraints_dirichlet_and_hanging_nodes, 
          ComponentMask()); 

//对解决方案的X和Y分量进行插值（这是一个位掩码，所以应用运算器|）。

        const FEValuesExtractors::Scalar x_displacement(0); 
        const FEValuesExtractors::Scalar y_displacement(1); 
        VectorTools::interpolate_boundary_values( 
          dof_handler, 
          8, 
          EquationData::BoundaryValues<dim>(), 
          constraints_dirichlet_and_hanging_nodes, 
          (fe.component_mask(x_displacement) | 
           fe.component_mask(y_displacement))); 
      } 
    else 
      VectorTools::interpolate_boundary_values( 
        dof_handler, 
        0, 
        EquationData::BoundaryValues<dim>(), 
        constraints_dirichlet_and_hanging_nodes, 
        ComponentMask()); 

    constraints_dirichlet_and_hanging_nodes.close(); 
  } 

//  @sect4{PlasticityContactProblem::assemble_mass_matrix_diagonal}  

// 下一个辅助函数计算（对角线）质量矩阵，用于确定我们在接触算法中使用的主动集合方法的主动集合。这个矩阵是质量矩阵类型的，但与标准质量矩阵不同，我们可以通过使用正交公式使其成为对角线（即使在高阶元素的情况下），该公式的正交点与有限元插值点的位置完全相同。我们通过使用QGaussLobatto正交公式来实现这一点，同时用一组从同一正交公式得出的插值点初始化有限元。该函数的其余部分相对简单：我们将得到的矩阵放入给定的参数中；因为我们知道矩阵是对角线的，所以只需在 $i$ 而不是 $j$ 上有一个循环即可。严格来说，我们甚至可以避免在正交点 <code>q_point</code> 处将形状函数的值与自身相乘，因为我们知道形状值是一个恰好有一个的向量，当与自身相点时产生1。由于这个函数不是时间关键，为了清楚起见，我们添加了这个术语。

  template <int dim> 
  void PlasticityContactProblem<dim>::assemble_mass_matrix_diagonal( 
    TrilinosWrappers::SparseMatrix &mass_matrix) 
  { 
    QGaussLobatto<dim - 1> face_quadrature_formula(fe.degree + 1); 

    FEFaceValues<dim> fe_values_face(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_JxW_values); 

    const unsigned int dofs_per_cell   = fe.n_dofs_per_cell(); 
    const unsigned int n_face_q_points = face_quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    const FEValuesExtractors::Vector displacement(0); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        for (const auto &face : cell->face_iterators()) 
          if (face->at_boundary() && face->boundary_id() == 1) 
            { 
              fe_values_face.reinit(cell, face); 
              cell_matrix = 0; 

              for (unsigned int q_point = 0; q_point < n_face_q_points; 
                   ++q_point) 
                for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                  cell_matrix(i, i) += 
                    (fe_values_face[displacement].value(i, q_point) * 
                     fe_values_face[displacement].value(i, q_point) * 
                     fe_values_face.JxW(q_point)); 

              cell->get_dof_indices(local_dof_indices); 

              for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                mass_matrix.add(local_dof_indices[i], 
                                local_dof_indices[i], 
                                cell_matrix(i, i)); 
            } 
    mass_matrix.compress(VectorOperation::add); 
  } 
// @sect4{PlasticityContactProblem::update_solution_and_constraints}  

// 下面的函数是我们在 <code>solve_newton()</code> 函数中每次牛顿迭代时调用的第一个函数。它的作用是将解决方案投射到可行集上，并更新接触或穿透障碍物的自由度的活动集。

// 为了实现这个功能，我们首先需要做一些记账工作。我们需要写入解决方案向量（我们只能用没有鬼魂元素的完全分布的向量来做），我们需要从各自的向量中读取拉格朗日乘数和对角线质量矩阵的元素（我们只能用有鬼魂元素的向量来做），所以我们创建各自的向量。然后我们还要初始化约束对象，该对象将包含来自接触和所有其他来源的约束，以及一个包含所有属于接触的本地自由度的索引集的对象。

  template <int dim> 
  void PlasticityContactProblem<dim>::update_solution_and_constraints() 
  { 
    std::vector<bool> dof_touched(dof_handler.n_dofs(), false); 

    TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, 
                                                       mpi_communicator); 
    distributed_solution = solution; 

    TrilinosWrappers::MPI::Vector lambda(locally_relevant_dofs, 
                                         mpi_communicator); 
    lambda = newton_rhs_uncondensed; 

    TrilinosWrappers::MPI::Vector diag_mass_matrix_vector_relevant( 
      locally_relevant_dofs, mpi_communicator); 
    diag_mass_matrix_vector_relevant = diag_mass_matrix_vector; 

    all_constraints.reinit(locally_relevant_dofs); 
    active_set.clear(); 

// 第二部分是在所有单元格上的循环，在这个循环中，我们看每一个自由度被定义的点的活动集条件是否为真，我们需要把这个自由度加入到接触节点的活动集中。正如我们一直所做的，如果我们想在单个点上评估函数，我们用一个FEValues对象（或者，这里是FEFaceValues对象，因为我们需要检查表面的接触）和一个适当选择的正交对象来做。我们通过选择定义在单元格面上的形状函数的 "支持点 "来创建这个面的正交对象（关于支持点的更多信息，请参见这个 @ref GlossSupport "词汇表条目"）。因此，我们有多少个正交点，就有多少个面的形状函数，在正交点上循环就相当于在面的形状函数上循环。有了这个，代码看起来如下。

    Quadrature<dim - 1> face_quadrature(fe.get_unit_face_support_points()); 
    FEFaceValues<dim>   fe_values_face(fe, 
                                     face_quadrature, 
                                     update_quadrature_points); 

    const unsigned int dofs_per_face   = fe.n_dofs_per_face(); 
    const unsigned int n_face_q_points = face_quadrature.size(); 

    std::vector<types::global_dof_index> dof_indices(dofs_per_face); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (!cell->is_artificial()) 
        for (const auto &face : cell->face_iterators()) 
          if (face->at_boundary() && face->boundary_id() == 1) 
            { 
              fe_values_face.reinit(cell, face); 
              face->get_dof_indices(dof_indices); 

              for (unsigned int q_point = 0; q_point < n_face_q_points; 
                   ++q_point) 
                { 

// 在每个正交点（即位于接触边界上的自由度的每个支持点），我们再询问它是否是z-位移自由度的一部分，如果我们还没有遇到这个自由度（对于那些位于面之间的边缘的自由度可能发生），我们需要评估变形物体与障碍物之间的间隙。如果活动集条件为真，那么我们在AffineConstraints对象中添加一个约束，下一次牛顿更新需要满足这个约束，将求解向量的相应元素设置为正确的值，并将索引添加到IndexSet对象中，该索引存储哪个自由度是接触的一部分。

                  const unsigned int component = 
                    fe.face_system_to_component_index(q_point).first; 

                  const unsigned int index_z = dof_indices[q_point]; 

                  if ((component == 2) && (dof_touched[index_z] == false)) 
                    { 
                      dof_touched[index_z] = true; 

                      const Point<dim> this_support_point = 
                        fe_values_face.quadrature_point(q_point); 

                      const double obstacle_value = 
                        obstacle->value(this_support_point, 2); 
                      const double solution_here = solution(index_z); 
                      const double undeformed_gap = 
                        obstacle_value - this_support_point(2); 

                      const double c = 100.0 * e_modulus; 
                      if ((lambda(index_z) / 
                               diag_mass_matrix_vector_relevant(index_z) + 
                             c * (solution_here - undeformed_gap) > 
                           0) && 
                          !constraints_hanging_nodes.is_constrained(index_z)) 
                        { 
                          all_constraints.add_line(index_z); 
                          all_constraints.set_inhomogeneity(index_z, 
                                                            undeformed_gap); 
                          distributed_solution(index_z) = undeformed_gap; 

                          active_set.add_index(index_z); 
                        } 
                    } 
                } 
            } 

// 在这个函数的最后，我们在处理器之间交换数据，更新 <code>solution</code> 变量中那些已经被其他处理器写入的幽灵元素。然后我们将Dirichlet约束和那些来自悬挂节点的约束合并到已经包含活动集的AffineConstraints对象中。我们通过输出主动约束自由度的总数来结束这个函数，对于这个自由度，我们对每个处理器拥有的主动约束自由度的数量进行加总。这个本地拥有的受限自由度的数量当然是活动集和本地拥有的自由度集的交集的元素数量，我们可以通过在两个IndexSets上使用 <code>operator&</code> 得到。

    distributed_solution.compress(VectorOperation::insert); 
    solution = distributed_solution; 

    all_constraints.close(); 
    all_constraints.merge(constraints_dirichlet_and_hanging_nodes); 

 
          << Utilities::MPI::sum((active_set & locally_owned_dofs).n_elements(), 
                                 mpi_communicator) 
          << std::endl; 
  } 
// @sect4{PlasticityContactProblem::assemble_newton_system}  

// 鉴于问题的复杂性，可能会让人感到惊讶的是，在每次牛顿迭代中组装我们要解决的线性系统实际上是相当简单的。下面的函数建立了牛顿的右手边和牛顿矩阵。它看起来相当简单，因为繁重的工作发生在对 <code>ConstitutiveLaw::get_linearized_stress_strain_tensors()</code> 的调用中，特别是在 AffineConstraints::distribute_local_to_global(), 中使用我们之前计算的约束。

  template <int dim> 
  void PlasticityContactProblem<dim>::assemble_newton_system( 
    const TrilinosWrappers::MPI::Vector &linearization_point) 
  { 
    TimerOutput::Scope t(computing_timer, "Assembling"); 

    QGauss<dim>     quadrature_formula(fe.degree + 1); 
    QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_JxW_values); 

    FEFaceValues<dim> fe_values_face(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_quadrature_points | 
                                       update_JxW_values); 

    const unsigned int dofs_per_cell   = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points      = quadrature_formula.size(); 
    const unsigned int n_face_q_points = face_quadrature_formula.size(); 

    const EquationData::BoundaryForce<dim> boundary_force; 
    std::vector<Vector<double>> boundary_force_values(n_face_q_points, 
                                                      Vector<double>(dim)); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    const FEValuesExtractors::Vector displacement(0); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          fe_values.reinit(cell); 
          cell_matrix = 0; 
          cell_rhs    = 0; 

          std::vector<SymmetricTensor<2, dim>> strain_tensor(n_q_points); 
          fe_values[displacement].get_function_symmetric_gradients( 
            linearization_point, strain_tensor); 

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
            { 
              SymmetricTensor<4, dim> stress_strain_tensor_linearized; 
              SymmetricTensor<4, dim> stress_strain_tensor; 
              constitutive_law.get_linearized_stress_strain_tensors( 
                strain_tensor[q_point], 
                stress_strain_tensor_linearized, 
                stress_strain_tensor); 

              for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                { 

// 在计算了应力-应变张量及其线性化之后，我们现在可以把矩阵和右手边的部分放在一起。在这两部分中，我们需要线性化的应力-应变张量乘以 $\varphi_i$ 的对称梯度，即 $I_\Pi\varepsilon(\varphi_i)$ 项，因此我们引入这个项的缩写。回顾一下，该矩阵对应于随附出版物的符号中的双线性形式 $A_{ij}=(I_\Pi\varepsilon(\varphi_i),\varepsilon(\varphi_j))$ ，而右手边是 $F_i=([I_\Pi-P_\Pi C]\varepsilon(\varphi_i),\varepsilon(\mathbf u))$ ，其中 $u$ 是当前的线性化点（通常是最后的解）。这可能表明，如果材料是完全弹性的（其中 $I_\Pi=P_\Pi$ ），右手边将为零，但这忽略了一个事实，即右手边还将包含由于接触而产生的非均质约束的贡献。                
//接下来的代码块增加了由于边界力的贡献，如果有的话。

                  const SymmetricTensor<2, dim> stress_phi_i = 
                    stress_strain_tensor_linearized * 
                    fe_values[displacement].symmetric_gradient(i, q_point); 

                  for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                    cell_matrix(i, j) += 
                      (stress_phi_i * 
                       fe_values[displacement].symmetric_gradient(j, q_point) * 
                       fe_values.JxW(q_point)); 

                  cell_rhs(i) += 
                    ((stress_phi_i - 
                      stress_strain_tensor * 
                        fe_values[displacement].symmetric_gradient(i, 
                                                                   q_point)) * 
                     strain_tensor[q_point] * fe_values.JxW(q_point)); 
                } 
            } 

          for (const auto &face : cell->face_iterators()) 
            if (face->at_boundary() && face->boundary_id() == 1) 
              { 
                fe_values_face.reinit(cell, face); 

                boundary_force.vector_value_list( 
                  fe_values_face.get_quadrature_points(), 
                  boundary_force_values); 

                for (unsigned int q_point = 0; q_point < n_face_q_points; 
                     ++q_point) 
                  { 
                    Tensor<1, dim> rhs_values; 
                    rhs_values[2] = boundary_force_values[q_point][2]; 
                    for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                      cell_rhs(i) += 
                        (fe_values_face[displacement].value(i, q_point) * 
                         rhs_values * fe_values_face.JxW(q_point)); 
                  } 
              } 

          cell->get_dof_indices(local_dof_indices); 
          all_constraints.distribute_local_to_global(cell_matrix, 
                                                     cell_rhs, 
                                                     local_dof_indices, 
                                                     newton_matrix, 
                                                     newton_rhs, 
                                                     true); 
        } 

    newton_matrix.compress(VectorOperation::add); 
    newton_rhs.compress(VectorOperation::add); 
  } 

//  @sect4{PlasticityContactProblem::compute_nonlinear_residual}  

// 下面的函数计算给定当前解（或任何其他线性化点）的方程的非线性残差。这在线性搜索算法中是需要的，我们需要尝试之前和当前（试验）解的各种线性组合来计算当前牛顿步骤的（真实的、全局化的）解。

// 说到这里，在稍微滥用函数名称的情况下，它实际上做了很多事情。例如，它还计算出与牛顿残差相对应的矢量，但没有消除受限自由度。我们需要这个向量来计算接触力，并最终计算出下一个活动集。同样，通过跟踪我们在每个单元上遇到的显示塑性屈服的正交点的数量，我们也可以计算出 <code>fraction_of_plastic_q_points_per_cell</code> 矢量，随后我们可以输出这个矢量来可视化塑性区。在这两种情况下，作为线条搜索的一部分，这些结果是不必要的，因此我们可能会浪费少量的时间来计算它们。同时，无论如何，这些信息是我们在这里需要做的事情的自然副产品，而且我们想在每个牛顿步骤结束时收集一次，所以我们不妨在这里做。

// 这个函数的实际实现应该是相当明显的。

  template <int dim> 
  void PlasticityContactProblem<dim>::compute_nonlinear_residual( 
    const TrilinosWrappers::MPI::Vector &linearization_point) 
  { 
    QGauss<dim>     quadrature_formula(fe.degree + 1); 
    QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_JxW_values); 

    FEFaceValues<dim> fe_values_face(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_quadrature_points | 
                                       update_JxW_values); 

    const unsigned int dofs_per_cell   = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points      = quadrature_formula.size(); 
    const unsigned int n_face_q_points = face_quadrature_formula.size(); 

    const EquationData::BoundaryForce<dim> boundary_force; 
    std::vector<Vector<double>> boundary_force_values(n_face_q_points, 
                                                      Vector<double>(dim)); 

    Vector<double> cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    const FEValuesExtractors::Vector displacement(0); 

    newton_rhs             = 0; 
    newton_rhs_uncondensed = 0; 

    fraction_of_plastic_q_points_per_cell = 0; 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          fe_values.reinit(cell); 
          cell_rhs = 0; 

          std::vector<SymmetricTensor<2, dim>> strain_tensors(n_q_points); 
          fe_values[displacement].get_function_symmetric_gradients( 
            linearization_point, strain_tensors); 

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
            { 
              SymmetricTensor<4, dim> stress_strain_tensor; 
              const bool              q_point_is_plastic = 
                constitutive_law.get_stress_strain_tensor( 
                  strain_tensors[q_point], stress_strain_tensor); 
              if (q_point_is_plastic) 
                ++fraction_of_plastic_q_points_per_cell( 
                  cell->active_cell_index()); 

              for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                { 
                  cell_rhs(i) -= 
                    (strain_tensors[q_point] * stress_strain_tensor * 
                     fe_values[displacement].symmetric_gradient(i, q_point) * 
                     fe_values.JxW(q_point)); 

                  Tensor<1, dim> rhs_values; 
                  rhs_values = 0; 
                  cell_rhs(i) += (fe_values[displacement].value(i, q_point) * 
                                  rhs_values * fe_values.JxW(q_point)); 
                } 
            } 

          for (const auto &face : cell->face_iterators()) 
            if (face->at_boundary() && face->boundary_id() == 1) 
              { 
                fe_values_face.reinit(cell, face); 

                boundary_force.vector_value_list( 
                  fe_values_face.get_quadrature_points(), 
                  boundary_force_values); 

                for (unsigned int q_point = 0; q_point < n_face_q_points; 
                     ++q_point) 
                  { 
                    Tensor<1, dim> rhs_values; 
                    rhs_values[2] = boundary_force_values[q_point][2]; 
                    for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                      cell_rhs(i) += 
                        (fe_values_face[displacement].value(i, q_point) * 
                         rhs_values * fe_values_face.JxW(q_point)); 
                  } 
              } 

          cell->get_dof_indices(local_dof_indices); 
          constraints_dirichlet_and_hanging_nodes.distribute_local_to_global( 
            cell_rhs, local_dof_indices, newton_rhs); 

          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            newton_rhs_uncondensed(local_dof_indices[i]) += cell_rhs(i); 
        } 

    fraction_of_plastic_q_points_per_cell /= quadrature_formula.size(); 
    newton_rhs.compress(VectorOperation::add); 
    newton_rhs_uncondensed.compress(VectorOperation::add); 
  } 

//  @sect4{PlasticityContactProblem::solve_newton_system}  

// 在我们讨论单个网格上的实际牛顿迭代之前的最后一块是线性系统的求解器。有几个复杂的问题使代码略显模糊，但大多数情况下，它只是设置然后求解。在这些复杂的问题中，包括。



// 对于悬空节点，我们必须将 AffineConstraints::set_zero 函数应用于newton_rhs。  如果一个求解值为 $x_0$ 的悬空节点有一个与障碍物接触的数值为 $x_1$ 的邻居和一个没有接触的邻居 $x_2$ ，这就有必要。因为前者的更新将是规定的，所以悬挂的节点约束将有一个不均匀性，看起来像  $x_0 = x_1/2 +   \text{gap}/2$  。所以右侧的相应条目是无意义的非零值。这些值我们必须设置为零。

// - 就像在  step-40  中一样，在求解或使用解决方案时，我们需要在有和没有鬼魂元素的向量之间进行洗牌。

// 该函数的其余部分与 step-40 和 step-41 类似，只是我们使用BiCGStab求解器而不是CG。这是由于对于非常小的硬化参数 $\gamma$ ，线性系统变得几乎是半无限的，尽管仍然是对称的。BiCGStab似乎更容易处理这种线性系统。

  template <int dim> 
  void PlasticityContactProblem<dim>::solve_newton_system() 
  { 
    TimerOutput::Scope t(computing_timer, "Solve"); 

    TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, 
                                                       mpi_communicator); 
    distributed_solution = solution; 

    constraints_hanging_nodes.set_zero(distributed_solution); 
    constraints_hanging_nodes.set_zero(newton_rhs); 

 
    { 
      TimerOutput::Scope t(computing_timer, "Solve: setup preconditioner"); 

      std::vector<std::vector<bool>> constant_modes; 
      DoFTools::extract_constant_modes(dof_handler, 
                                       ComponentMask(), 
                                       constant_modes); 

      TrilinosWrappers::PreconditionAMG::AdditionalData additional_data; 
      additional_data.constant_modes        = constant_modes; 
      additional_data.elliptic              = true; 
      additional_data.n_cycles              = 1; 
      additional_data.w_cycle               = false; 
      additional_data.output_details        = false; 
      additional_data.smoother_sweeps       = 2; 
      additional_data.aggregation_threshold = 1e-2; 

      preconditioner.initialize(newton_matrix, additional_data); 
    } 

    { 
      TimerOutput::Scope t(computing_timer, "Solve: iterate"); 

      TrilinosWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator); 

      const double relative_accuracy = 1e-8; 
      const double solver_tolerance = 
        relative_accuracy * 
        newton_matrix.residual(tmp, distributed_solution, newton_rhs); 

      SolverControl solver_control(newton_matrix.m(), solver_tolerance); 
      SolverBicgstab<TrilinosWrappers::MPI::Vector> solver(solver_control); 
      solver.solve(newton_matrix, 
                   distributed_solution, 
                   newton_rhs, 
                   preconditioner); 

      pcout << "         Error: " << solver_control.initial_value() << " -> " 
            << solver_control.last_value() << " in " 
            << solver_control.last_step() << " Bicgstab iterations." 
            << std::endl; 
    } 

    all_constraints.distribute(distributed_solution); 

    solution = distributed_solution; 
  } 
// @sect4{PlasticityContactProblem::solve_newton}  

// 最后，这是在当前网格上实现阻尼牛顿方法的函数。这里有两个嵌套的循环：外循环用于牛顿迭代，内循环用于直线搜索，只有在必要时才会使用。为了获得一个好的和合理的起始值，我们在每个网格上的第一个牛顿步骤中解决一个弹性问题（如果我们在网格之间转移解决方案，则只在第一个网格上解决）。我们通过在这些迭代中将屈服应力设置为一个不合理的大值，然后在随后的迭代中将其设置为正确值。

// 除此以外，这个函数的顶部部分应该是相当明显的。我们将变量 <code>previous_residual_norm</code> 初始化为可以用双精度数字表示的最大负值，以便在第一步中比较当前残差是否小于前一步的残差时总是失败。

  template <int dim> 
  void PlasticityContactProblem<dim>::solve_newton() 
  { 
    TrilinosWrappers::MPI::Vector old_solution(locally_owned_dofs, 
                                               mpi_communicator); 
    TrilinosWrappers::MPI::Vector residual(locally_owned_dofs, 
                                           mpi_communicator); 
    TrilinosWrappers::MPI::Vector tmp_vector(locally_owned_dofs, 
                                             mpi_communicator); 
    TrilinosWrappers::MPI::Vector locally_relevant_tmp_vector( 
      locally_relevant_dofs, mpi_communicator); 
    TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, 
                                                       mpi_communicator); 

    double residual_norm; 
    double previous_residual_norm = -std::numeric_limits<double>::max(); 

    const double correct_sigma = sigma_0; 

    IndexSet old_active_set(active_set); 

    for (unsigned int newton_step = 1; newton_step <= 100; ++newton_step) 
      { 
        if (newton_step == 1 && 
            ((transfer_solution && current_refinement_cycle == 0) || 
             !transfer_solution)) 
          constitutive_law.set_sigma_0(1e+10); 
        else if (newton_step == 2 || current_refinement_cycle > 0 || 
                 !transfer_solution) 
          constitutive_law.set_sigma_0(correct_sigma); 

        pcout << " " << std::endl; 
        pcout << "   Newton iteration " << newton_step << std::endl; 
        pcout << "      Updating active set..." << std::endl; 

        { 
          TimerOutput::Scope t(computing_timer, "update active set"); 
          update_solution_and_constraints(); 
        } 

        pcout << "      Assembling system... " << std::endl; 
        newton_matrix = 0; 
        newton_rhs    = 0; 
        assemble_newton_system(solution); 

        pcout << "      Solving system... " << std::endl; 
        solve_newton_system(); 

// 在我们计算了当前牛顿步骤的试解 $\tilde{\mathbf u}$ 之后，情况就变得有点棘手了。我们处理的是一个高度非线性的问题，所以我们必须用直线搜索的方式来抑制牛顿方法。为了理解我们如何做到这一点，请回顾一下，在我们的表述中，我们在每一个牛顿步骤中计算一个试解，而不是在新旧解之间进行更新。由于解集是一个凸集，我们将使用直线搜索，尝试以前的解和试验解的线性组合，以保证阻尼解再次出现在我们的解集中。我们最多应用5个阻尼步骤。

// 在我们使用直线搜索的时候有一些例外情况。首先，如果这是任何网格上的第一个牛顿步骤，那么我们就没有任何点来比较残差，所以我们总是接受一个完整的步骤。同样地，如果这是第一个网格上的第二个牛顿步骤（如果我们不在网格之间转移解决方案，则是任何网格上的第二个牛顿步骤），则我们只用弹性模型计算了其中的第一个步骤（见上文我们如何将屈服应力σ设置为一个不合理的大值）。在这种情况下，第一个牛顿解是一个纯粹的弹性解，第二个牛顿解是一个塑性解，任何线性组合都不一定会位于可行的集合中--所以我们只是接受我们刚刚得到的解。

// 在这两种情况下，我们绕过直线搜索，只是在必要时更新残差和其他向量。

        if ((newton_step == 1) || 
            (transfer_solution && newton_step == 2 && 
             current_refinement_cycle == 0) || 
            (!transfer_solution && newton_step == 2)) 
          { 
            compute_nonlinear_residual(solution); 
            old_solution = solution; 

            residual                     = newton_rhs; 
            const unsigned int start_res = (residual.local_range().first), 
                               end_res   = (residual.local_range().second); 
            for (unsigned int n = start_res; n < end_res; ++n) 
              if (all_constraints.is_inhomogeneously_constrained(n)) 
                residual(n) = 0; 

            residual.compress(VectorOperation::insert); 

            residual_norm = residual.l2_norm(); 

            pcout << "      Accepting Newton solution with residual: " 
                  << residual_norm << std::endl; 
          } 
        else 
          { 
            for (unsigned int i = 0; i < 5; ++i) 
              { 
                distributed_solution = solution; 

                const double alpha = std::pow(0.5, static_cast<double>(i)); 
                tmp_vector         = old_solution; 
                tmp_vector.sadd(1 - alpha, alpha, distributed_solution); 

                TimerOutput::Scope t(computing_timer, "Residual and lambda"); 

                locally_relevant_tmp_vector = tmp_vector; 
                compute_nonlinear_residual(locally_relevant_tmp_vector); 
                residual = newton_rhs; 

                const unsigned int start_res = (residual.local_range().first), 
                                   end_res   = (residual.local_range().second); 
                for (unsigned int n = start_res; n < end_res; ++n) 
                  if (all_constraints.is_inhomogeneously_constrained(n)) 
                    residual(n) = 0; 

                residual.compress(VectorOperation::insert); 

                residual_norm = residual.l2_norm(); 

 
                  << "      Residual of the non-contact part of the system: " 
                  << residual_norm << std::endl 
                  << "         with a damping parameter alpha = " << alpha 
                  << std::endl; 

                if (residual_norm < previous_residual_norm) 
                  break; 
              } 

            solution     = tmp_vector; 
            old_solution = solution; 
          } 

        previous_residual_norm = residual_norm; 

// 最后一步是检查收敛情况。如果活动集在所有处理器中都没有变化，并且残差小于阈值 $10^{-10}$  ，那么我们就终止对当前网格的迭代。

        if (Utilities::MPI::sum((active_set == old_active_set) ? 0 : 1, 
                                mpi_communicator) == 0) 
          { 
            pcout << "      Active set did not change!" << std::endl; 
            if (residual_norm < 1e-10) 
              break; 
          } 

        old_active_set = active_set; 
      } 
  } 
// @sect4{PlasticityContactProblem::refine_grid}  

// 如果你已经在deal.II教程中做到了这一点，下面这个细化网格的函数应该不会再对你构成任何挑战。它对网格进行细化，可以是全局的，也可以是使用Kelly误差估计器的，如果这样要求的话，还可以将上一个网格的解转移到下一个网格。在后一种情况下，我们还需要再次计算活动集和其他数量，为此我们需要由  <code>compute_nonlinear_residual()</code>  计算的信息。

  template <int dim> 
  void PlasticityContactProblem<dim>::refine_grid() 
  { 
    if (refinement_strategy == RefinementStrategy::refine_global) 
      { 
        for (typename Triangulation<dim>::active_cell_iterator cell = 
               triangulation.begin_active(); 
             cell != triangulation.end(); 
             ++cell) 
          if (cell->is_locally_owned()) 
            cell->set_refine_flag(); 
      } 
    else 
      { 
        Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
        KellyErrorEstimator<dim>::estimate( 
          dof_handler, 
          QGauss<dim - 1>(fe.degree + 2), 
          std::map<types::boundary_id, const Function<dim> *>(), 
          solution, 
          estimated_error_per_cell); 

        parallel::distributed::GridRefinement ::refine_and_coarsen_fixed_number( 
          triangulation, estimated_error_per_cell, 0.3, 0.03); 
      } 

    triangulation.prepare_coarsening_and_refinement(); 

    parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> 
      solution_transfer(dof_handler); 
    if (transfer_solution) 
      solution_transfer.prepare_for_coarsening_and_refinement(solution); 

    triangulation.execute_coarsening_and_refinement(); 

    setup_system(); 

    if (transfer_solution) 
      { 
        TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, 
                                                           mpi_communicator); 
        solution_transfer.interpolate(distributed_solution); 

// 强制执行约束条件，使插值后的解决方案在新的网格上符合要求。

        constraints_hanging_nodes.distribute(distributed_solution); 

        solution = distributed_solution; 
        compute_nonlinear_residual(solution); 
      } 
  } 
// @sect4{PlasticityContactProblem::move_mesh}  

// 在我们到达 <code>run()</code> 之前的其余三个函数都与生成输出有关。下面一个是尝试显示变形体的变形构造。为此，这个函数接收一个位移矢量场，通过先前计算的位移来移动网格（局部）的每个顶点。在生成图形输出之前，我们将以当前的位移场调用该函数，在生成图形输出之后，我们将以负的位移场再次调用该函数，以撤销对网格所做的修改。

// 这个函数本身是非常简单的。我们所要做的就是跟踪我们已经接触过的顶点，因为我们在单元格上循环时多次遇到相同的顶点。

  template <int dim> 
  void PlasticityContactProblem<dim>::move_mesh( 
    const TrilinosWrappers::MPI::Vector &displacement) const 
  { 
    std::vector<bool> vertex_touched(triangulation.n_vertices(), false); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        for (const auto v : cell->vertex_indices()) 
          if (vertex_touched[cell->vertex_index(v)] == false) 
            { 
              vertex_touched[cell->vertex_index(v)] = true; 

              Point<dim> vertex_displacement; 
              for (unsigned int d = 0; d < dim; ++d) 
                vertex_displacement[d] = 
                  displacement(cell->vertex_dof_index(v, d)); 

              cell->vertex(v) += vertex_displacement; 
            } 
  } 

//  @sect4{PlasticityContactProblem::output_results}  

// 接下来是我们用来实际生成图形输出的函数。这个函数有点繁琐，但实际上并不特别复杂。它在顶部移动网格（最后再把它移回来），然后计算沿接触面的接触力。我们可以通过取未处理的残差向量，并通过询问它们是否有与之相关的不均匀约束来确定哪些自由度对应于有接触的自由度（如随文所示）。一如既往，我们需要注意的是，我们只能写进完全分布的向量（即没有鬼魂元素的向量），但当我们想产生输出时，我们需要的向量确实对所有局部相关的自由度都有鬼魂项。

  template <int dim> 
  void PlasticityContactProblem<dim>::output_results( 
    const unsigned int current_refinement_cycle) 
  { 
    TimerOutput::Scope t(computing_timer, "Graphical output"); 

    pcout << "      Writing graphical output... " << std::flush; 

    move_mesh(solution); 

// 接触力的计算

    TrilinosWrappers::MPI::Vector distributed_lambda(locally_owned_dofs, 
                                                     mpi_communicator); 
    const unsigned int start_res = (newton_rhs_uncondensed.local_range().first), 
                       end_res = (newton_rhs_uncondensed.local_range().second); 
    for (unsigned int n = start_res; n < end_res; ++n) 
      if (all_constraints.is_inhomogeneously_constrained(n)) 
        distributed_lambda(n) = 
          newton_rhs_uncondensed(n) / diag_mass_matrix_vector(n); 
    distributed_lambda.compress(VectorOperation::insert); 
    constraints_hanging_nodes.distribute(distributed_lambda); 

    TrilinosWrappers::MPI::Vector lambda(locally_relevant_dofs, 
                                         mpi_communicator); 
    lambda = distributed_lambda; 

    TrilinosWrappers::MPI::Vector distributed_active_set_vector( 
      locally_owned_dofs, mpi_communicator); 
    distributed_active_set_vector = 0.; 
    for (const auto index : active_set) 
      distributed_active_set_vector[index] = 1.; 
    distributed_lambda.compress(VectorOperation::insert); 

    TrilinosWrappers::MPI::Vector active_set_vector(locally_relevant_dofs, 
                                                    mpi_communicator); 
    active_set_vector = distributed_active_set_vector; 

    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 

    const std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      data_component_interpretation( 
        dim, DataComponentInterpretation::component_is_part_of_vector); 
    data_out.add_data_vector(solution, 
                             std::vector<std::string>(dim, "displacement"), 
                             DataOut<dim>::type_dof_data, 
                             data_component_interpretation); 
    data_out.add_data_vector(lambda, 
                             std::vector<std::string>(dim, "contact_force"), 
                             DataOut<dim>::type_dof_data, 
                             data_component_interpretation); 
    data_out.add_data_vector(active_set_vector, 
                             std::vector<std::string>(dim, "active_set"), 
                             DataOut<dim>::type_dof_data, 
                             data_component_interpretation); 

    Vector<float> subdomain(triangulation.n_active_cells()); 
    for (unsigned int i = 0; i < subdomain.size(); ++i) 
      subdomain(i) = triangulation.locally_owned_subdomain(); 
    data_out.add_data_vector(subdomain, "subdomain"); 

    data_out.add_data_vector(fraction_of_plastic_q_points_per_cell, 
                             "fraction_of_plastic_q_points"); 

    data_out.build_patches(); 

// 在函数的其余部分，我们在每个处理器上生成一个VTU文件，以这个处理器的子域ID为索引。在第一个处理器上，我们随后还创建了一个 <code>.pvtu</code> 文件，对VTU文件的<i>all</i>进行索引，这样就可以一次性读取整个输出文件集。这些 <code>.pvtu</code> 被Paraview用来描述整个并行计算的输出文件。然后我们再为Paraview的竞争者--VisIt可视化程序做同样的事情，创建一个匹配的 <code>.visit</code> 文件。

    const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record( 
      output_dir, "solution", current_refinement_cycle, mpi_communicator, 2); 
    pcout << pvtu_filename << std::endl; 

    TrilinosWrappers::MPI::Vector tmp(solution); 
    tmp *= -1; 
    move_mesh(tmp); 
  } 
// @sect4{PlasticityContactProblem::output_contact_force}  

// 这最后一个辅助函数通过计算接触面积上Z方向的接触压力的积分来计算接触力。为此，我们将所有非活动因子的接触压力lambda设置为0（一个自由度是否是接触的一部分，就像我们在前一个函数中做的那样）。对于所有活动的自由度，lambda包含非线性残差（newton_rhs_uncondensed）和质量矩阵（diag_mass_matrix_vector）的相应对角线条目的商数。因为悬空节点出现在接触区的可能性不小，所以对分布式_lambda向量应用constraints_hanging_nodes.distribution是很重要的。

  template <int dim> 
  void PlasticityContactProblem<dim>::output_contact_force() const 
  { 
    TrilinosWrappers::MPI::Vector distributed_lambda(locally_owned_dofs, 
                                                     mpi_communicator); 
    const unsigned int start_res = (newton_rhs_uncondensed.local_range().first), 
                       end_res = (newton_rhs_uncondensed.local_range().second); 
    for (unsigned int n = start_res; n < end_res; ++n) 
      if (all_constraints.is_inhomogeneously_constrained(n)) 
        distributed_lambda(n) = 
          newton_rhs_uncondensed(n) / diag_mass_matrix_vector(n); 
      else 
        distributed_lambda(n) = 0; 
    distributed_lambda.compress(VectorOperation::insert); 
    constraints_hanging_nodes.distribute(distributed_lambda); 

    TrilinosWrappers::MPI::Vector lambda(locally_relevant_dofs, 
                                         mpi_communicator); 
    lambda = distributed_lambda; 

    double contact_force = 0.0; 

    QGauss<dim - 1>   face_quadrature_formula(fe.degree + 1); 
    FEFaceValues<dim> fe_values_face(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_JxW_values); 

    const unsigned int n_face_q_points = face_quadrature_formula.size(); 

    const FEValuesExtractors::Vector displacement(0); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        for (const auto &face : cell->face_iterators()) 
          if (face->at_boundary() && face->boundary_id() == 1) 
            { 
              fe_values_face.reinit(cell, face); 

              std::vector<Tensor<1, dim>> lambda_values(n_face_q_points); 
              fe_values_face[displacement].get_function_values(lambda, 
                                                               lambda_values); 

              for (unsigned int q_point = 0; q_point < n_face_q_points; 
                   ++q_point) 
                contact_force += 
                  lambda_values[q_point][2] * fe_values_face.JxW(q_point); 
            } 
    contact_force = Utilities::MPI::sum(contact_force, MPI_COMM_WORLD); 

    pcout << "Contact force = " << contact_force << std::endl; 
  } 
// @sect4{PlasticityContactProblem::run}  

// 和其他所有的教程程序一样， <code>run()</code> 函数包含了整体逻辑。这里没有太多的内容：本质上，它在所有的网格细化循环中执行循环，并在每个循环中，将事情交给 <code>solve_newton()</code> 中的牛顿求解器，并调用函数来创建如此计算的解决方案的图形输出。然后输出一些关于运行时间和内存消耗的统计数据，这些数据是在这个网格的计算过程中收集的。

  template <int dim> 
  void PlasticityContactProblem<dim>::run() 
  { 
    computing_timer.reset(); 
    for (; current_refinement_cycle < n_refinement_cycles; 
         ++current_refinement_cycle) 
      { 
        { 
          TimerOutput::Scope t(computing_timer, "Setup"); 

          pcout << std::endl; 
          pcout << "Cycle " << current_refinement_cycle << ':' << std::endl; 

          if (current_refinement_cycle == 0) 
            { 
              make_grid(); 
              setup_system(); 
            } 
          else 
            { 
              TimerOutput::Scope t(computing_timer, "Setup: refine mesh"); 
              refine_grid(); 
            } 
        } 

        solve_newton(); 

        output_results(current_refinement_cycle); 

        computing_timer.print_summary(); 
        computing_timer.reset(); 

        Utilities::System::MemoryStats stats; 
        Utilities::System::get_memory_stats(stats); 
        pcout << "Peak virtual memory used, resident in kB: " << stats.VmSize 
              << " " << stats.VmRSS << std::endl; 

        if (base_mesh == "box") 
          output_contact_force(); 
      } 
  } 
} // namespace Step42 
// @sect3{The <code>main</code> function}  

//  <code>main()</code> 函数真的没有什么内容。看起来他们总是这样做。

int main(int argc, char *argv[]) 
{ 
  using namespace dealii; 
  using namespace Step42; 

  try 
    { 
      ParameterHandler prm; 
      PlasticityContactProblem<3>::declare_parameters(prm); 
      if (argc != 2) 
        { 
          std::cerr << "*** Call this program as <./step-42 input.prm>" 
                    << std::endl; 
          return 1; 
        } 

      prm.parse_input(argv[1]); 
      Utilities::MPI::MPI_InitFinalize mpi_initialization( 
        argc, argv, numbers::invalid_unsigned_int); 
      { 
        PlasticityContactProblem<3> problem(prm); 
        problem.run(); 
      } 
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


