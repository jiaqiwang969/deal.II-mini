

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2018 - 2021 by the deal.II authors 
 * 
 * This file is part of the deal.II library. 
 * 
 * The deal.II library is free software; you can use it, redistribute 
 * it, and/or modify it under the terms of the GNU Lesser General 
 * Public License as published by the Free Software Foundation; either 
 * version 2.1 of the License, or (at your option) any later version. 
 * The full text of the license can be found in the file LICENSE at 
 * the top level of the deal.II distribution. 
 * 
 * --------------------------------------------------------------------- 
 * 
 * Author: Daniel Garcia-Sanchez, CNRS, 2019 
 */ 


// @sect3{Include files}  

// 我们在这个程序中需要的大部分包含文件已经在以前的程序中讨论过了，特别是在  step-40  .

#include <deal.II/base/conditional_ostream.h> 
#include <deal.II/base/function.h> 

#include <deal.II/base/index_set.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/generic_linear_algebra.h> 
#include <deal.II/lac/petsc_solver.h> 
#include <deal.II/lac/vector.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

#include <fstream> 
#include <iostream> 

// 下面的标头提供了我们用来表示材料属性的张量类。

#include <deal.II/base/tensor.h> 

// 下面的标头对于deal.II的HDF5接口是必要的。

#include <deal.II/base/hdf5.h> 

// 这个头是我们用来评估模拟结果的函数 VectorTools::point_value 所需要的。

#include <deal.II/numerics/vector_tools.h> 

// 我们在函数 GridTools::find_active_cell_around_point 中使用的函数 `ElasticWave::store_frequency_step_data()` 需要这些头文件。
#include <deal.II/grid/grid_tools.h> 
#include <deal.II/grid/grid_tools_cache.h> 

namespace step62 
{ 
  using namespace dealii; 
// @sect3{Auxiliary classes and functions}  下列类用于存储模拟的参数。

//  @sect4{The `RightHandSide` class}  该类用于定义结构左侧的力脉冲。

  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    RightHandSide(HDF5::Group &data); 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component) const override; 

  private: 

// 变量`data`是 HDF5::Group ，所有的模拟结果都将被储存在其中。请注意， `RightHandSide::data`, 变量 
// `PML::data`,  
// `Rho::data` 和 `Parameters::data` 指向HDF5文件的同一个组。当 HDF5::Group 被复制时，它将指向HDF5文件的同一组。

    HDF5::Group data; 

// 仿真参数作为HDF5属性存储在`data`中。以下属性在jupyter笔记本中定义，作为HDF5属性存储在`data`中，然后由构造函数读取。

    const double     max_force_amplitude; 
    const double     force_sigma_x; 
    const double     force_sigma_y; 
    const double     max_force_width_x; 
    const double     max_force_width_y; 
    const Point<dim> force_center; 

  public: 

// 在这个特定的模拟中，力只有一个 $x$ 分量， $F_y=0$  。

    const unsigned int force_component = 0; 
  }; 
// @sect4{The `PML` class}  这个类是用来定义完美匹配层（PML）的形状，以吸收向边界传播的波。

  template <int dim> 
  class PML : public Function<dim, std::complex<double>> 
  { 
  public: 
    PML(HDF5::Group &data); 

    virtual std::complex<double> 
    value(const Point<dim> &p, const unsigned int component) const override; 

  private: 
// HDF5::Group ，所有的模拟结果将被存储在其中。

    HDF5::Group data; 

// 和以前一样，以下属性在jupyter笔记本中定义，作为HDF5属性存储在`data`中，然后由构造函数读取。

    const double pml_coeff; 
    const int    pml_coeff_degree; 
    const double dimension_x; 
    const double dimension_y; 
    const bool   pml_x; 
    const bool   pml_y; 
    const double pml_width_x; 
    const double pml_width_y; 
    const double a_coeff_x; 
    const double a_coeff_y; 
  }; 

//  @sect4{The `Rho` class}  这个类是用来定义质量密度的。

  template <int dim> 
  class Rho : public Function<dim> 
  { 
  public: 
    Rho(HDF5::Group &data); 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

  private: 
// HDF5::Group ，所有的模拟结果将被存储在其中。

    HDF5::Group data; 

// 和以前一样，以下属性在jupyter笔记本中定义，作为HDF5属性存储在`data`中，然后由构造函数读取。

    const double       lambda; 
    const double       mu; 
    const double       material_a_rho; 
    const double       material_b_rho; 
    const double       cavity_resonance_frequency; 
    const unsigned int nb_mirror_pairs; 
    const double       dimension_y; 
    const unsigned int grid_level; 
    double             average_rho_width; 
  }; 

//  @sect4{The `Parameters` class}  该类包含所有将在模拟中使用的参数。

  template <int dim> 
  class Parameters 
  { 
  public: 
    Parameters(HDF5::Group &data); 
// HDF5::Group ，所有的模拟结果将被存储在其中。

    HDF5::Group data; 

// 和以前一样，以下属性在jupyter笔记本中定义，作为HDF5属性存储在`data`中，然后由构造函数读取。

    const std::string        simulation_name; 
    const bool               save_vtu_files; 
    const double             start_frequency; 
    const double             stop_frequency; 
    const unsigned int       nb_frequency_points; 
    const double             lambda; 
    const double             mu; 
    const double             dimension_x; 
    const double             dimension_y; 
    const unsigned int       nb_probe_points; 
    const unsigned int       grid_level; 
    const Point<dim>         probe_start_point; 
    const Point<dim>         probe_stop_point; 
    const RightHandSide<dim> right_hand_side; 
    const PML<dim>           pml; 
    const Rho<dim>           rho; 

  private: 
    const double comparison_float_constant = 1e-12; 
  }; 

//  @sect4{The `QuadratureCache` class}  质量和刚度矩阵的计算是非常昂贵的。这些矩阵对所有的频率步骤都是一样的。右手边的向量对所有的频率步长也是一样的。我们用这个类来存储这些对象，并在每个频率步骤中重新使用它们。请注意，这里我们不存储集合的质量和刚度矩阵以及右手边，而是存储单个单元的数据。QuadratureCache "类与在  step-18  中使用过的 "PointHistory "类非常相似。

  template <int dim> 
  class QuadratureCache 
  { 
  public: 
    QuadratureCache(const unsigned int dofs_per_cell); 

  private: 
    unsigned int dofs_per_cell; 

  public: 

// 我们在变量mass_coefficient和stiffness_coefficient中存储质量和刚度矩阵。我们还存储了右手边和JxW值，这些值对所有的频率步骤都是一样的。

    FullMatrix<std::complex<double>>  mass_coefficient; 
    FullMatrix<std::complex<double>>  stiffness_coefficient; 
    std::vector<std::complex<double>> right_hand_side; 
    double                            JxW; 
  }; 

//  @sect4{The `get_stiffness_tensor()` function}  

// 该函数返回材料的刚度张量。为了简单起见，我们认为刚度是各向同性和同质的；只有密度  $\rho$  取决于位置。正如我们之前在  step-8  中所表明的，如果刚度是各向同性和均质的，那么刚度系数  $c_{ijkl}$  可以表示为两个系数  $\lambda$  和  $\mu$  的函数。系数张量简化为 
// @f[
//    c_{ijkl}
//    =
//    \lambda \delta_{ij} \delta_{kl} +
//    \mu (\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}).
//  @f] 。

  template <int dim> 
  SymmetricTensor<4, dim> get_stiffness_tensor(const double lambda, 
                                               const double mu) 
  { 
    SymmetricTensor<4, dim> stiffness_tensor; 
    for (unsigned int i = 0; i < dim; ++i) 
      for (unsigned int j = 0; j < dim; ++j) 
        for (unsigned int k = 0; k < dim; ++k) 
          for (unsigned int l = 0; l < dim; ++l) 
            stiffness_tensor[i][j][k][l] = 
              (((i == k) && (j == l) ? mu : 0.0) + 
               ((i == l) && (j == k) ? mu : 0.0) + 
               ((i == j) && (k == l) ? lambda : 0.0)); 
    return stiffness_tensor; 
  } 

//  @sect3{The `ElasticWave` class}  

// 接下来让我们声明这个程序的主类。它的结构与 step-40 的教程程序非常相似。主要的区别是。

// - 扫过的频率值。

// - 我们将刚度和质量矩阵保存在`quadrature_cache`中，并在每个频率步骤中使用它们。

// - 我们在HDF5文件中存储每个频率步骤的探头测量的能量。

  template <int dim> 
  class ElasticWave 
  { 
  public: 
    ElasticWave(const Parameters<dim> &parameters); 
    void run(); 

  private: 
    void setup_system(); 
    void assemble_system(const double omega, 
                         const bool   calculate_quadrature_data); 
    void solve(); 
    void initialize_probe_positions_vector(); 
    void store_frequency_step_data(const unsigned int frequency_idx); 
    void output_results(); 

// 在每个频率步骤之前都会调用这个，以便为缓存变量设置一个原始状态。

    void setup_quadrature_cache(); 

// 这个函数在频率向量上循环，并在每个频率步骤上运行模拟。

    void frequency_sweep(); 

// 参数存储在这个变量中。

    Parameters<dim> parameters; 

    MPI_Comm mpi_communicator; 

    parallel::distributed::Triangulation<dim> triangulation; 

    QGauss<dim> quadrature_formula; 

// 我们将每个单元的质量和刚度矩阵存储在这个向量中。

    std::vector<QuadratureCache<dim>> quadrature_cache; 

    FESystem<dim>   fe; 
    DoFHandler<dim> dof_handler; 

    IndexSet locally_owned_dofs; 
    IndexSet locally_relevant_dofs; 

    AffineConstraints<std::complex<double>> constraints; 

    LinearAlgebraPETSc::MPI::SparseMatrix system_matrix; 
    LinearAlgebraPETSc::MPI::Vector       locally_relevant_solution; 
    LinearAlgebraPETSc::MPI::Vector       system_rhs; 

// 这个向量包含我们要模拟的频率范围。

    std::vector<double> frequency; 

// 这个向量包含了测量探头各点的坐标 $(x,y)$ 。

    FullMatrix<double> probe_positions; 

// HDF5数据集来存储频率和`探头位置`向量。

    HDF5::DataSet frequency_dataset; 
    HDF5::DataSet probe_positions_dataset; 

// HDF5数据集，存储探头测量的能量值。

    HDF5::DataSet displacement; 

    ConditionalOStream pcout; 
    TimerOutput        computing_timer; 
  }; 

//  @sect3{Implementation of the auxiliary classes}  
// @sect4{The `RightHandSide` class implementation}  

// 构造函数使用 HDF5::Group  `data`函数从 HDF5::Group::get_attribute()  读取所有参数。

  template <int dim> 
  RightHandSide<dim>::RightHandSide(HDF5::Group &data) 
    : Function<dim>(dim) 
    , data(data) 
    , max_force_amplitude(data.get_attribute<double>("max_force_amplitude")) 
    , force_sigma_x(data.get_attribute<double>("force_sigma_x")) 
    , force_sigma_y(data.get_attribute<double>("force_sigma_y")) 
    , max_force_width_x(data.get_attribute<double>("max_force_width_x")) 
    , max_force_width_y(data.get_attribute<double>("max_force_width_y")) 
    , force_center(Point<dim>(data.get_attribute<double>("force_x_pos"), 
                              data.get_attribute<double>("force_y_pos"))) 
  {} 

//这个函数定义了力矢量脉冲的空间形状，它采取高斯函数
// @f{align*}
//  F_x &=
//  \left\{
//  \begin{array}{ll}
//    a \exp(- (\frac{(x-b_x)^2 }{ 2 \sigma_x^2}+\frac{(y-b_y)^2 }{ 2
//    \sigma_y^2}))
//  & \text{if}\, x_\textrm{min} <x<x_\textrm{max}\, \text{and}\,
//  y_\textrm{min} <y<y_\textrm{max}  \\ 0 & \text{otherwise},
//  \end{array}
//  \right.\\ F_y &= 0
//  @f}
//  的形式，其中 $a$ 是取力的最大振幅， $\sigma_x$ 和 $\sigma_y$ 是 $x$ 和 $y$ 分量的标准偏差。请注意，脉冲已被裁剪为 $x_\textrm{min}<x<x_\textrm{max}$ 和 $y_\textrm{min} <y<y_\textrm{max}$  。

  template <int dim> 
  double RightHandSide<dim>::value(const Point<dim> & p, 
                                   const unsigned int component) const 
  { 
    if (component == force_component) 
      { 
        if (std::abs(p[0] - force_center[0]) < max_force_width_x / 2 && 
            std::abs(p[1] - force_center[1]) < max_force_width_y / 2) 
          { 
            return max_force_amplitude * 
                   std::exp(-(std::pow(p[0] - force_center[0], 2) / 
                                (2 * std::pow(force_sigma_x, 2)) + 
                              std::pow(p[1] - force_center[1], 2) / 
                                (2 * std::pow(force_sigma_y, 2)))); 
          } 
        else 
          { 
            return 0; 
          } 
      } 
    else 
      { 
        return 0; 
      } 
  } 

//  @sect4{The `PML` class implementation}  

// 和以前一样，构造函数使用 HDF5::Group 函数从 HDF5::Group::get_attribute() `data`中读取所有参数。正如我们所讨论的，在jupyter笔记本中已经定义了PML的二次开机。通过改变参数`pml_coeff_degree`，可以使用线性、立方或其他幂度。参数`pml_x`和`pml_y`可以用来开启和关闭`x`和`y`PML。

  template <int dim> 
  PML<dim>::PML(HDF5::Group &data) 
    : Function<dim, std::complex<double>>(dim) 
    , data(data) 
    , pml_coeff(data.get_attribute<double>("pml_coeff")) 
    , pml_coeff_degree(data.get_attribute<int>("pml_coeff_degree")) 
    , dimension_x(data.get_attribute<double>("dimension_x")) 
    , dimension_y(data.get_attribute<double>("dimension_y")) 
    , pml_x(data.get_attribute<bool>("pml_x")) 
    , pml_y(data.get_attribute<bool>("pml_y")) 
    , pml_width_x(data.get_attribute<double>("pml_width_x")) 
    , pml_width_y(data.get_attribute<double>("pml_width_y")) 
    , a_coeff_x(pml_coeff / std::pow(pml_width_x, pml_coeff_degree)) 
    , a_coeff_y(pml_coeff / std::pow(pml_width_y, pml_coeff_degree)) 
  {} 

// `x`部分的PML系数的形式为  $s'_x = a_x x^{\textrm{degree}}$  。
  template <int dim> 
  std::complex<double> PML<dim>::value(const Point<dim> & p, 
                                       const unsigned int component) const 
  { 
    double calculated_pml_x_coeff = 0; 
    double calculated_pml_y_coeff = 0; 

    if ((component == 0) && pml_x) 
      { 
        const double pml_x_start_position = dimension_x / 2 - pml_width_x; 
        if (std::abs(p[0]) > pml_x_start_position) 
          { 
            const double x_prime = std::abs(p[0]) - pml_x_start_position; 
            calculated_pml_x_coeff = 
              a_coeff_x * std::pow(x_prime, pml_coeff_degree); 
          } 
      } 

    if ((component == 1) && pml_y) 
      { 
        const double pml_y_start_position = dimension_y / 2 - pml_width_y; 
        if (std::abs(p[1]) > pml_y_start_position) 
          { 
            const double y_prime = std::abs(p[1]) - pml_y_start_position; 
            calculated_pml_y_coeff = 
              a_coeff_y * std::pow(y_prime, pml_coeff_degree); 
          } 
      } 

    return 1. + std::max(calculated_pml_x_coeff, calculated_pml_y_coeff) * 
                  std::complex<double>(0., 1.); 
  } 

//  @sect4{The `Rho` class implementation}  

// 这个类是用来定义质量密度的。正如我们之前所解释的，一个声学超晶格空腔是由两个[分布式反射器](https:en.wikipedia.org/wiki/Band_gap)、镜子和一个 $\lambda/2$ 空腔组成的，其中 $\lambda$ 是声波长。声学DBRs是一种周期性结构，其中一组具有对比性物理特性（声速指数）的双层堆栈被重复 $N$ 次。波速的变化是由具有不同密度的层交替产生的。

  template <int dim> 
  Rho<dim>::Rho(HDF5::Group &data) 
    : Function<dim>(1) 
    , data(data) 
    , lambda(data.get_attribute<double>("lambda")) 
    , mu(data.get_attribute<double>("mu")) 
    , material_a_rho(data.get_attribute<double>("material_a_rho")) 
    , material_b_rho(data.get_attribute<double>("material_b_rho")) 
    , cavity_resonance_frequency( 
        data.get_attribute<double>("cavity_resonance_frequency")) 
    , nb_mirror_pairs(data.get_attribute<int>("nb_mirror_pairs")) 
    , dimension_y(data.get_attribute<double>("dimension_y")) 
    , grid_level(data.get_attribute<int>("grid_level")) 
  { 

// 为了提高精度，我们使用[subpixel smoothing]（https:meep.readthedocs.io/en/latest/Subpixel_Smoothing/）。

    average_rho_width = dimension_y / (std::pow(2.0, grid_level)); 
    data.set_attribute("average_rho_width", average_rho_width); 
  } 

  template <int dim> 
  double Rho<dim>::value(const Point<dim> &p, 
                         const unsigned int /*component*/) const 
  { 

// 声速由
// @f[
//   c = \frac{K_e}{\rho}
//  @f]
//  定义，其中 $K_e$ 是有效弹性常数， $\rho$ 是密度。这里我们考虑的是波导宽度远小于波长的情况。在这种情况下，可以证明对于二维的情况
//  @f[
//   K_e = 4\mu\frac{\lambda +\mu}{\lambda+2\mu}
//  @f]
//  和三维的情况 $K_e$ 等于杨氏模量。
//  @f[
//   K_e = \mu\frac{3\lambda +2\mu}{\lambda+\mu}
//  @f]

    double elastic_constant; 
    if (dim == 2) 
      { 
        elastic_constant = 4 * mu * (lambda + mu) / (lambda + 2 * mu); 
      } 
    else if (dim == 3) 
      { 
        elastic_constant = mu * (3 * lambda + 2 * mu) / (lambda + mu); 
      } 
    else 
      { 
        Assert(false, ExcInternalError()); 
      } 
    const double material_a_speed_of_sound = 
      std::sqrt(elastic_constant / material_a_rho); 
    const double material_a_wavelength = 
      material_a_speed_of_sound / cavity_resonance_frequency; 
    const double material_b_speed_of_sound = 
      std::sqrt(elastic_constant / material_b_rho); 
    const double material_b_wavelength = 
      material_b_speed_of_sound / cavity_resonance_frequency; 

//密度 $\rho$ 采取以下形式 <img alt="声学超晶格空腔" src="https:www.dealii.org/images/steps/developer/  step-62  .04.svg" height="200" //其中棕色代表材料_a，绿色代表材料_b。

    for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++) 
      { 
        const double layer_transition_center = 
          material_a_wavelength / 2 + 
          idx * (material_b_wavelength / 4 + material_a_wavelength / 4); 
        if (std::abs(p[0]) >= 
              (layer_transition_center - average_rho_width / 2) && 
            std::abs(p[0]) <= (layer_transition_center + average_rho_width / 2)) 
          { 
            const double coefficient = 
              (std::abs(p[0]) - 
               (layer_transition_center - average_rho_width / 2)) / 
              average_rho_width; 
            return (1 - coefficient) * material_a_rho + 
                   coefficient * material_b_rho; 
          } 
      } 

// 这里我们定义了[subpixel smoothing](https:meep.readthedocs.io/en/latest/Subpixel_Smoothing/)，它可以提高模拟的精度。

    for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++) 
      { 
        const double layer_transition_center = 
          material_a_wavelength / 2 + 
          idx * (material_b_wavelength / 4 + material_a_wavelength / 4) + 
          material_b_wavelength / 4; 
        if (std::abs(p[0]) >= 
              (layer_transition_center - average_rho_width / 2) && 
            std::abs(p[0]) <= (layer_transition_center + average_rho_width / 2)) 
          { 
            const double coefficient = 
              (std::abs(p[0]) - 
               (layer_transition_center - average_rho_width / 2)) / 
              average_rho_width; 
            return (1 - coefficient) * material_b_rho + 
                   coefficient * material_a_rho; 
          } 
      } 

// 然后是腔体

    if (std::abs(p[0]) <= material_a_wavelength / 2) 
      { 
        return material_a_rho; 
      } 

// 材料层_a

    for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++) 
      { 
        const double layer_center = 
          material_a_wavelength / 2 + 
          idx * (material_b_wavelength / 4 + material_a_wavelength / 4) + 
          material_b_wavelength / 4 + material_a_wavelength / 8; 
        const double layer_width = material_a_wavelength / 4; 
        if (std::abs(p[0]) >= (layer_center - layer_width / 2) && 
            std::abs(p[0]) <= (layer_center + layer_width / 2)) 
          { 
            return material_a_rho; 
          } 
      } 

// material_b层

    for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++) 
      { 
        const double layer_center = 
          material_a_wavelength / 2 + 
          idx * (material_b_wavelength / 4 + material_a_wavelength / 4) + 
          material_b_wavelength / 8; 
        const double layer_width = material_b_wavelength / 4; 
        if (std::abs(p[0]) >= (layer_center - layer_width / 2) && 
            std::abs(p[0]) <= (layer_center + layer_width / 2)) 
          { 
            return material_b_rho; 
          } 
      } 

// 最后，默认的是 material_a。

    return material_a_rho; 
  } 

//  @sect4{The `Parameters` class implementation}  

// 构造函数使用 HDF5::Group 函数从 HDF5::Group::get_attribute() `data`中读取所有参数。

  template <int dim> 
  Parameters<dim>::Parameters(HDF5::Group &data) 
    : data(data) 
    , simulation_name(data.get_attribute<std::string>("simulation_name")) 
    , save_vtu_files(data.get_attribute<bool>("save_vtu_files")) 
    , start_frequency(data.get_attribute<double>("start_frequency")) 
    , stop_frequency(data.get_attribute<double>("stop_frequency")) 
    , nb_frequency_points(data.get_attribute<int>("nb_frequency_points")) 
    , lambda(data.get_attribute<double>("lambda")) 
    , mu(data.get_attribute<double>("mu")) 
    , dimension_x(data.get_attribute<double>("dimension_x")) 
    , dimension_y(data.get_attribute<double>("dimension_y")) 
    , nb_probe_points(data.get_attribute<int>("nb_probe_points")) 
    , grid_level(data.get_attribute<int>("grid_level")) 
    , probe_start_point(data.get_attribute<double>("probe_pos_x"), 
                        data.get_attribute<double>("probe_pos_y") - 
                          data.get_attribute<double>("probe_width_y") / 2) 
    , probe_stop_point(data.get_attribute<double>("probe_pos_x"), 
                       data.get_attribute<double>("probe_pos_y") + 
                         data.get_attribute<double>("probe_width_y") / 2) 
    , right_hand_side(data) 
    , pml(data) 
    , rho(data) 
  {} 

//  @sect4{The `QuadratureCache` class implementation}  

// 我们需要为质量和刚度矩阵以及右手边的矢量保留足够的空间。

  template <int dim> 
  QuadratureCache<dim>::QuadratureCache(const unsigned int dofs_per_cell) 
    : dofs_per_cell(dofs_per_cell) 
    , mass_coefficient(dofs_per_cell, dofs_per_cell) 
    , stiffness_coefficient(dofs_per_cell, dofs_per_cell) 
    , right_hand_side(dofs_per_cell) 
  {} 

//  @sect3{Implementation of the `ElasticWave` class}  
// @sect4{Constructor}  

// 这与  step-40  的构造函数非常相似。此外，我们还创建了HDF5数据集`frequency_dataset`，`position_dataset`和`displacement`。注意在创建HDF5数据集时使用了 "模板 "关键字。这是C++的要求，使用`template`关键字是为了将`create_dataset`作为一个依赖的模板名称。

  template <int dim> 
  ElasticWave<dim>::ElasticWave(const Parameters<dim> &parameters) 
    : parameters(parameters) 
    , mpi_communicator(MPI_COMM_WORLD) 
    , triangulation(mpi_communicator, 
                    typename Triangulation<dim>::MeshSmoothing( 
                      Triangulation<dim>::smoothing_on_refinement | 
                      Triangulation<dim>::smoothing_on_coarsening)) 
    , quadrature_formula(2) 
    , fe(FE_Q<dim>(1), dim) 
    , dof_handler(triangulation) 
    , frequency(parameters.nb_frequency_points) 
    , probe_positions(parameters.nb_probe_points, dim) 
    , frequency_dataset(parameters.data.template create_dataset<double>( 
        "frequency", 
        std::vector<hsize_t>{parameters.nb_frequency_points})) 
    , probe_positions_dataset(parameters.data.template create_dataset<double>( 
        "position", 
        std::vector<hsize_t>{parameters.nb_probe_points, dim})) 
    , displacement( 
        parameters.data.template create_dataset<std::complex<double>>( 
          "displacement", 
          std::vector<hsize_t>{parameters.nb_probe_points, 
                               parameters.nb_frequency_points})) 
    , pcout(std::cout, 
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)) 
    , computing_timer(mpi_communicator, 
                      pcout, 
                      TimerOutput::summary, 
                      TimerOutput::wall_times) 
  {} 

//  @sect4{ElasticWave::setup_system}  

// 这个函数没有什么新内容，与 step-40 的唯一区别是，我们不需要应用边界条件，因为我们使用PML来截断域。

  template <int dim> 
  void ElasticWave<dim>::setup_system() 
  { 
    TimerOutput::Scope t(computing_timer, "setup"); 

    dof_handler.distribute_dofs(fe); 

    locally_owned_dofs = dof_handler.locally_owned_dofs(); 
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs); 

    locally_relevant_solution.reinit(locally_owned_dofs, 
                                     locally_relevant_dofs, 
                                     mpi_communicator); 

    system_rhs.reinit(locally_owned_dofs, mpi_communicator); 

    constraints.clear(); 
    constraints.reinit(locally_relevant_dofs); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 

    constraints.close(); 

    DynamicSparsityPattern dsp(locally_relevant_dofs); 

    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false); 
    SparsityTools::distribute_sparsity_pattern(dsp, 
                                               locally_owned_dofs, 
                                               mpi_communicator, 
                                               locally_relevant_dofs); 

    system_matrix.reinit(locally_owned_dofs, 
                         locally_owned_dofs, 
                         dsp, 
                         mpi_communicator); 
  } 

//  @sect4{ElasticWave::assemble_system}  

// 这个函数也与 step-40 非常相似，尽管有明显的区别。我们为每个频率/欧米茄步骤组装系统。在第一步中，我们设置`calculate_quadrature_data = True`，然后我们计算质量和刚度矩阵以及右手边的矢量。在随后的步骤中，我们将使用这些数据来加速计算。

  template <int dim> 
  void ElasticWave<dim>::assemble_system(const double omega, 
                                         const bool   calculate_quadrature_data) 
  { 
    TimerOutput::Scope t(computing_timer, "assembly"); 

    FEValues<dim>      fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<std::complex<double>> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<std::complex<double>>     cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

// 这里我们存储右手边的值，rho和PML的值。

    std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim)); 
    std::vector<double>         rho_values(n_q_points); 
    std::vector<Vector<std::complex<double>>> pml_values( 
      n_q_points, Vector<std::complex<double>>(dim)); 

// 我们计算已经在jupyter笔记本中定义的 $\lambda$ 和 $\mu$ 的刚度张量。请注意，与 $\rho$ 相反，刚度在整个领域中是恒定的。

    const SymmetricTensor<4, dim> stiffness_tensor = 
      get_stiffness_tensor<dim>(parameters.lambda, parameters.mu); 

// 我们使用与 step-20 相同的方法处理矢量值问题。

    const FEValuesExtractors::Vector displacement(0); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          cell_matrix = 0; 
          cell_rhs    = 0; 

// 只有当我们要计算质量和刚度矩阵时，我们才必须计算右手边的rho和PML的值。否则我们可以跳过这个计算，这样可以大大减少总的计算时间。

          if (calculate_quadrature_data) 
            { 
              fe_values.reinit(cell); 

              parameters.right_hand_side.vector_value_list( 
                fe_values.get_quadrature_points(), rhs_values); 
              parameters.rho.value_list(fe_values.get_quadrature_points(), 
                                        rho_values); 
              parameters.pml.vector_value_list( 
                fe_values.get_quadrature_points(), pml_values); 
            } 

// 我们已经在  step-18  中做了这个工作。获得一个指向当前单元本地正交缓存数据的指针，作为防御措施，确保这个指针在全局数组的范围内。

          QuadratureCache<dim> *local_quadrature_points_data = 
            reinterpret_cast<QuadratureCache<dim> *>(cell->user_pointer()); 
          Assert(local_quadrature_points_data >= &quadrature_cache.front(), 
                 ExcInternalError()); 
          Assert(local_quadrature_points_data <= &quadrature_cache.back(), 
                 ExcInternalError()); 
          for (unsigned int q = 0; q < n_q_points; ++q) 
            { 

// quadrature_data变量用于存储质量和刚度矩阵、右手边向量和`JxW`的值。

              QuadratureCache<dim> &quadrature_data = 
                local_quadrature_points_data[q]; 

// 下面我们声明力向量和PML的参数  $s$  和  $\xi$  。

              Tensor<1, dim>                       force; 
              Tensor<1, dim, std::complex<double>> s; 
              std::complex<double>                 xi(1, 0); 

// 下面的块只在第一个频率步骤中计算。

              if (calculate_quadrature_data) 
                { 

// 存储`JxW`的值。

                  quadrature_data.JxW = fe_values.JxW(q); 

                  for (unsigned int component = 0; component < dim; ++component) 
                    { 

// 将向量转换为张量，并计算出xi

                      force[component] = rhs_values[q][component]; 
                      s[component]     = pml_values[q][component]; 
                      xi *= s[component]; 
                    } 

// 这里我们计算 $\alpha_{mnkl}$ 和 $\beta_{mnkl}$ 张量。

                  Tensor<4, dim, std::complex<double>> alpha; 
                  Tensor<4, dim, std::complex<double>> beta; 
                  for (unsigned int m = 0; m < dim; ++m) 
                    for (unsigned int n = 0; n < dim; ++n) 
                      for (unsigned int k = 0; k < dim; ++k) 
                        for (unsigned int l = 0; l < dim; ++l) 
                          { 
                            alpha[m][n][k][l] = xi * 
                                                stiffness_tensor[m][n][k][l] / 
                                                (2.0 * s[n] * s[k]); 
                            beta[m][n][k][l] = xi * 
                                               stiffness_tensor[m][n][k][l] / 
                                               (2.0 * s[n] * s[l]); 
                          } 

                  for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                    { 
                      const Tensor<1, dim> phi_i = 
                        fe_values[displacement].value(i, q); 
                      const Tensor<2, dim> grad_phi_i = 
                        fe_values[displacement].gradient(i, q); 

                      for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                        { 
                          const Tensor<1, dim> phi_j = 
                            fe_values[displacement].value(j, q); 
                          const Tensor<2, dim> grad_phi_j = 
                            fe_values[displacement].gradient(j, q); 

// 计算质量矩阵的值。

                          quadrature_data.mass_coefficient[i][j] = 
                            rho_values[q] * xi * phi_i * phi_j; 

//在刚度张量的 $mnkl$ 指数上循环。

                          std::complex<double> stiffness_coefficient = 0; 
                          for (unsigned int m = 0; m < dim; ++m) 
                            for (unsigned int n = 0; n < dim; ++n) 
                              for (unsigned int k = 0; k < dim; ++k) 
                                for (unsigned int l = 0; l < dim; ++l) 
                                  { 

// 这里我们计算刚度矩阵。                          
//注意，由于PML的存在，刚度矩阵不是对称的。我们使用梯度函数（见[文档](https:www.dealii.org/current/doxygen/deal.II/group__vector__valued.html)），它是一个  <code>Tensor@<2,dim@></code>  。                          
// 矩阵 $G_{ij}$ 由条目
                          // @f[
                          //  G_{ij}=
                          //  \frac{\partial\phi_i}{\partial x_j}
                          //  =\partial_j \phi_i
                          // @f]
                          // 组成 注意指数 $i$ 和 $j$ 的位置以及我们在本教程中使用的符号。  $\partial_j\phi_i$  . 由于刚度张量不是对称的，所以很容易出错。

                                    stiffness_coefficient += 
                                      grad_phi_i[m][n] * 
                                      (alpha[m][n][k][l] * grad_phi_j[l][k] + 
                                       beta[m][n][k][l] * grad_phi_j[k][l]); 
                                  } 

// 我们将刚度矩阵的值保存在quadrature_data中。

                          quadrature_data.stiffness_coefficient[i][j] = 
                            stiffness_coefficient; 
                        } 

// 和正交数据中的右手边的值。

 
                        phi_i * force * fe_values.JxW(q); 
                    } 
                } 

// 我们再次循环单元的自由度来计算系统矩阵。这些循环非常快，因为我们已经计算了刚度和质量矩阵，只有 $\omega$ 的值发生了变化。

              for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                { 
                  for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                    { 
                      std::complex<double> matrix_sum = 0; 
                      matrix_sum += -std::pow(omega, 2) * 
                                    quadrature_data.mass_coefficient[i][j]; 
                      matrix_sum += quadrature_data.stiffness_coefficient[i][j]; 
                      cell_matrix(i, j) += matrix_sum * quadrature_data.JxW; 
                    } 
                  cell_rhs(i) += quadrature_data.right_hand_side[i]; 
                } 
            } 
          cell->get_dof_indices(local_dof_indices); 
          constraints.distribute_local_to_global(cell_matrix, 
                                                 cell_rhs, 
                                                 local_dof_indices, 
                                                 system_matrix, 
                                                 system_rhs); 
        } 

    system_matrix.compress(VectorOperation::add); 
    system_rhs.compress(VectorOperation::add); 
  } 
// @sect4{ElasticWave::solve}  

// 这比  step-40  更加简单。我们使用并行的直接求解器MUMPS，它比迭代求解器需要更少的选项。缺点是它不能很好地扩展。用迭代求解器来解决Helmholtz方程并不简单。移位拉普拉斯多网格法是一种众所周知的预处理该系统的方法，但这超出了本教程的范围。

  template <int dim> 
  void ElasticWave<dim>::solve() 
  { 
    TimerOutput::Scope              t(computing_timer, "solve"); 
    LinearAlgebraPETSc::MPI::Vector completely_distributed_solution( 
      locally_owned_dofs, mpi_communicator); 

    SolverControl                    solver_control; 
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator); 
    solver.solve(system_matrix, completely_distributed_solution, system_rhs); 

    pcout << "   Solved in " << solver_control.last_step() << " iterations." 
          << std::endl; 
    constraints.distribute(completely_distributed_solution); 
    locally_relevant_solution = completely_distributed_solution; 
  } 
// @sect4{ElasticWave::initialize_position_vector}  

// 我们用这个函数来计算位置向量的值。

  template <int dim> 
  void ElasticWave<dim>::initialize_probe_positions_vector() 
  { 
    for (unsigned int position_idx = 0; 
         position_idx < parameters.nb_probe_points; 
         ++position_idx) 
      { 

// 由于运算符+和

// -被重载来减去两个点，所以必须做如下操作。`Point_b<dim> + (-Point_a<dim>)`。

        const Point<dim> p = 
          (position_idx / ((double)(parameters.nb_probe_points - 1))) * 
            (parameters.probe_stop_point + (-parameters.probe_start_point)) + 
          parameters.probe_start_point; 
        probe_positions[position_idx][0] = p[0]; 
        probe_positions[position_idx][1] = p[1]; 
        if (dim == 3) 
          { 
            probe_positions[position_idx][2] = p[2]; 
          } 
      } 
  } 
// @sect4{ElasticWave::store_frequency_step_data}  

// 该函数在HDF5文件中存储探头测量的能量。

  template <int dim> 
  void 
  ElasticWave<dim>::store_frequency_step_data(const unsigned int frequency_idx) 
  { 
    TimerOutput::Scope t(computing_timer, "store_frequency_step_data"); 

// 我们存储 $x$ 方向的位移； $y$ 方向的位移可以忽略不计。

    const unsigned int probe_displacement_component = 0; 

// 向量坐标包含HDF5文件中位于本地所有单元中的探测点的坐标。向量displacement_data包含这些点的位移值。

    std::vector<hsize_t>              coordinates; 
    std::vector<std::complex<double>> displacement_data; 

    const auto &mapping = get_default_linear_mapping(triangulation); 
    GridTools::Cache<dim, dim> cache(triangulation, mapping); 
    typename Triangulation<dim, dim>::active_cell_iterator cell_hint{}; 
    std::vector<bool>                                      marked_vertices = {}; 
    const double                                           tolerance = 1.e-10; 

    for (unsigned int position_idx = 0; 
         position_idx < parameters.nb_probe_points; 
         ++position_idx) 
      { 
        Point<dim> point; 
        for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx) 
          { 
            point[dim_idx] = probe_positions[position_idx][dim_idx]; 
          } 
        bool point_in_locally_owned_cell = false; 
        { 
          auto cell_and_ref_point = GridTools::find_active_cell_around_point( 
            cache, point, cell_hint, marked_vertices, tolerance); 
          if (cell_and_ref_point.first.state() == IteratorState::valid) 
            { 
              cell_hint = cell_and_ref_point.first; 
              point_in_locally_owned_cell = 
                cell_and_ref_point.first->is_locally_owned(); 
            } 
        } 
        if (point_in_locally_owned_cell) 
          { 

// 然后，我们可以在`displacement_data`中存储探头各点的位移值。

            Vector<std::complex<double>> tmp_vector(dim); 
            VectorTools::point_value(dof_handler, 
                                     locally_relevant_solution, 
                                     point, 
                                     tmp_vector); 
            coordinates.emplace_back(position_idx); 
            coordinates.emplace_back(frequency_idx); 
            displacement_data.emplace_back( 
              tmp_vector(probe_displacement_component)); 
          } 
      } 

// 我们在HDF5文件中写入位移数据。调用 HDF5::DataSet::write_selection() 是MPI集体的，这意味着所有进程都要参与。

    if (coordinates.size() > 0) 
      { 
        displacement.write_selection(displacement_data, coordinates); 
      } 

// 因此，即使进程没有数据可写，它也必须参与集体调用。为此我们可以使用  HDF5::DataSet::write_none().  注意，我们必须指定数据类型，在这种情况下  `std::complex<double>`.  。
    else 
      { 
        displacement.write_none<std::complex<double>>(); 
      } 

// 如果输入文件中的变量`save_vtu_files`等于`True`，那么所有数据将被保存为vtu。写入`vtu'文件的过程已经在  step-40  中描述。

    if (parameters.save_vtu_files) 
      { 
        std::vector<std::string> solution_names(dim, "displacement"); 
        std::vector<DataComponentInterpretation::DataComponentInterpretation> 
          interpretation( 
            dim, DataComponentInterpretation::component_is_part_of_vector); 

        DataOut<dim> data_out; 
        data_out.add_data_vector(dof_handler, 
                                 locally_relevant_solution, 
                                 solution_names, 
                                 interpretation); 
        Vector<float> subdomain(triangulation.n_active_cells()); 
        for (unsigned int i = 0; i < subdomain.size(); ++i) 
          subdomain(i) = triangulation.locally_owned_subdomain(); 
        data_out.add_data_vector(subdomain, "subdomain"); 

        std::vector<Vector<double>> force( 
          dim, Vector<double>(triangulation.n_active_cells())); 
        std::vector<Vector<double>> pml( 
          dim, Vector<double>(triangulation.n_active_cells())); 
        Vector<double> rho(triangulation.n_active_cells()); 

        for (auto &cell : triangulation.active_cell_iterators()) 
          { 
            if (cell->is_locally_owned()) 
              { 
                for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx) 
                  { 
                    force[dim_idx](cell->active_cell_index()) = 
                      parameters.right_hand_side.value(cell->center(), dim_idx); 
                    pml[dim_idx](cell->active_cell_index()) = 
                      parameters.pml.value(cell->center(), dim_idx).imag(); 
                  } 
                rho(cell->active_cell_index()) = 
                  parameters.rho.value(cell->center()); 
              } 

// 在我们不感兴趣的单元格上，将各自的值设置为一个假值，以确保如果我们的假设有什么错误，我们会通过查看图形输出发现。

            else 
              { 
                for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx) 
                  { 
                    force[dim_idx](cell->active_cell_index()) = -1e+20; 
                    pml[dim_idx](cell->active_cell_index())   = -1e+20; 
                  } 
                rho(cell->active_cell_index()) = -1e+20; 
              } 
          } 

        for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx) 
          { 
            data_out.add_data_vector(force[dim_idx], 
                                     "force_" + std::to_string(dim_idx)); 
            data_out.add_data_vector(pml[dim_idx], 
                                     "pml_" + std::to_string(dim_idx)); 
          } 
        data_out.add_data_vector(rho, "rho"); 

        data_out.build_patches(); 

        std::stringstream  frequency_idx_stream; 
        const unsigned int nb_number_positions = 
          ((unsigned int)std::log10(parameters.nb_frequency_points)) + 1; 
        frequency_idx_stream << std::setw(nb_number_positions) 
                             << std::setfill('0') << frequency_idx; 
        std::string filename = (parameters.simulation_name + "_" + 
                                frequency_idx_stream.str() + ".vtu"); 
        data_out.write_vtu_in_parallel(filename.c_str(), mpi_communicator); 
      } 
  } 

//  @sect4{ElasticWave::output_results}  

// 该函数写入尚未写入的数据集。

  template <int dim> 
  void ElasticWave<dim>::output_results() 
  { 

// 向量`频率`和`位置`对所有进程都是一样的。因此任何一个进程都可以写入相应的`数据集'。因为调用 HDF5::DataSet::write 是MPI集体的，其余进程将不得不调用 HDF5::DataSet::write_none.  。
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0) 
      { 
        frequency_dataset.write(frequency); 
        probe_positions_dataset.write(probe_positions); 
      } 
    else 
      { 
        frequency_dataset.write_none<double>(); 
        probe_positions_dataset.write_none<double>(); 
      } 
  } 

//  @sect4{ElasticWave::setup_quadrature_cache}  

// 我们在计算开始时使用这个函数来设置缓存变量的初始值。这个函数在  step-18  中已经描述过。与  step-18  的函数没有区别。

  template <int dim> 
  void ElasticWave<dim>::setup_quadrature_cache() 
  { 
    triangulation.clear_user_data(); 

    { 
      std::vector<QuadratureCache<dim>> tmp; 
      quadrature_cache.swap(tmp); 
    } 

    quadrature_cache.resize(triangulation.n_locally_owned_active_cells() * 
                              quadrature_formula.size(), 
                            QuadratureCache<dim>(fe.n_dofs_per_cell())); 
    unsigned int cache_index = 0; 
    for (const auto &cell : triangulation.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          cell->set_user_pointer(&quadrature_cache[cache_index]); 
          cache_index += quadrature_formula.size(); 
        } 
    Assert(cache_index == quadrature_cache.size(), ExcInternalError()); 
  } 

//  @sect4{ElasticWave::frequency_sweep}  

// 为了清楚起见，我们将 step-40 的函数`run`分为函数`run`和`frequency_sweep`。在函数`frequency_sweep`中，我们把迭代放在频率向量上。

  template <int dim> 
  void ElasticWave<dim>::frequency_sweep() 
  { 
    for (unsigned int frequency_idx = 0; 
         frequency_idx < parameters.nb_frequency_points; 
         ++frequency_idx) 
      { 
        pcout << parameters.simulation_name + " frequency idx: " 
              << frequency_idx << '/' << parameters.nb_frequency_points - 1 
              << std::endl; 

        setup_system(); 
        if (frequency_idx == 0) 
          { 
            pcout << "   Number of active cells :       " 
                  << triangulation.n_active_cells() << std::endl; 
            pcout << "   Number of degrees of freedom : " 
                  << dof_handler.n_dofs() << std::endl; 
          } 

        if (frequency_idx == 0) 
          { 

// 只写一次模拟参数

            parameters.data.set_attribute("active_cells", 
                                          triangulation.n_active_cells()); 
            parameters.data.set_attribute("degrees_of_freedom", 
                                          dof_handler.n_dofs()); 
          } 

// 我们计算出这个特定步骤的频率和欧米茄值。

        const double current_loop_frequency = 
          (parameters.start_frequency + 
           frequency_idx * 
             (parameters.stop_frequency - parameters.start_frequency) / 
             (parameters.nb_frequency_points - 1)); 
        const double current_loop_omega = 
          2 * numbers::PI * current_loop_frequency; 

// 在第一个频率步骤中，我们计算出质量和刚度矩阵以及右手边的数据。在随后的频率步骤中，我们将使用这些值。这大大改善了计算时间。

        assemble_system(current_loop_omega, 
                        (frequency_idx == 0) ? true : false); 
        solve(); 

        frequency[frequency_idx] = current_loop_frequency; 
        store_frequency_step_data(frequency_idx); 

        computing_timer.print_summary(); 
        computing_timer.reset(); 
        pcout << std::endl; 
      } 
  } 

//  @sect4{ElasticWave::run}  

// 这个函数与  step-40  中的函数非常相似。

  template <int dim> 
  void ElasticWave<dim>::run() 
  { 
#ifdef DEBUG 
    pcout << "Debug mode" << std::endl; 
#else 
    pcout << "Release mode" << std::endl; 
#endif 

    { 
      Point<dim> p1; 
      p1(0) = -parameters.dimension_x / 2; 
      p1(1) = -parameters.dimension_y / 2; 
      if (dim == 3) 
        { 
          p1(2) = -parameters.dimension_y / 2; 
        } 
      Point<dim> p2; 
      p2(0) = parameters.dimension_x / 2; 
      p2(1) = parameters.dimension_y / 2; 
      if (dim == 3) 
        { 
          p2(2) = parameters.dimension_y / 2; 
        } 
      std::vector<unsigned int> divisions(dim); 
      divisions[0] = int(parameters.dimension_x / parameters.dimension_y); 
      divisions[1] = 1; 
      if (dim == 3) 
        { 
          divisions[2] = 1; 
        } 
      GridGenerator::subdivided_hyper_rectangle(triangulation, 
                                                divisions, 
                                                p1, 
                                                p2); 
    } 

    triangulation.refine_global(parameters.grid_level); 

    setup_quadrature_cache(); 

    initialize_probe_positions_vector(); 

    frequency_sweep(); 

    output_results(); 
  } 
} // namespace step62 

//  @sect4{The main function}  

// 主函数与  step-40  中的函数非常相似。

int main(int argc, char *argv[]) 
{ 
  try 
    { 
      using namespace dealii; 
      const unsigned int dim = 2; 

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 

      HDF5::File data_file("results.h5", 
                           HDF5::File::FileAccessMode::create, 
                           MPI_COMM_WORLD); 
      auto       data = data_file.create_group("data"); 

// 每个模拟（位移和校准）都存储在一个单独的HDF5组中。

      const std::vector<std::string> group_names = {"displacement", 
                                                    "calibration"}; 
      for (auto group_name : group_names) 
        { 

// 对于这两个组名中的每一个，我们现在创建组并将属性放入这些组。具体来说，这些是。

// - 波导的尺寸（在 $x$ 和 $y$ 方向）。

// - 探头的位置（在 $x$ 和 $y$ 方向）。

// - 探针中的点的数量

// - 全局细化水平

// - 腔体谐振频率

// - 镜像对的数量 

// - 镜子的数量 

// - 材料特性

// - 力的参数

// - PML参数

// - 频率参数

          auto group = data.create_group(group_name); 

          group.set_attribute<double>("dimension_x", 2e-5); 
          group.set_attribute<double>("dimension_y", 2e-8); 
          group.set_attribute<double>("probe_pos_x", 8e-6); 
          group.set_attribute<double>("probe_pos_y", 0); 
          group.set_attribute<double>("probe_width_y", 2e-08); 
          group.set_attribute<unsigned int>("nb_probe_points", 5); 
          group.set_attribute<unsigned int>("grid_level", 1); 
          group.set_attribute<double>("cavity_resonance_frequency", 20e9); 
          group.set_attribute<unsigned int>("nb_mirror_pairs", 15); 

          group.set_attribute<double>("poissons_ratio", 0.27); 
          group.set_attribute<double>("youngs_modulus", 270000000000.0); 
          group.set_attribute<double>("material_a_rho", 3200); 

          if (group_name == std::string("displacement")) 
            group.set_attribute<double>("material_b_rho", 2000); 
          else 
            group.set_attribute<double>("material_b_rho", 3200); 

          group.set_attribute( 
            "lambda", 
            group.get_attribute<double>("youngs_modulus") * 
              group.get_attribute<double>("poissons_ratio") / 
              ((1 + group.get_attribute<double>("poissons_ratio")) * 
               (1 - 2 * group.get_attribute<double>("poissons_ratio")))); 
          group.set_attribute("mu", 
                              group.get_attribute<double>("youngs_modulus") / 
                                (2 * (1 + group.get_attribute<double>( 
                                            "poissons_ratio")))); 

          group.set_attribute<double>("max_force_amplitude", 1e26); 
          group.set_attribute<double>("force_sigma_x", 1e-7); 
          group.set_attribute<double>("force_sigma_y", 1); 
          group.set_attribute<double>("max_force_width_x", 3e-7); 
          group.set_attribute<double>("max_force_width_y", 2e-8); 
          group.set_attribute<double>("force_x_pos", -8e-6); 
          group.set_attribute<double>("force_y_pos", 0); 

          group.set_attribute<bool>("pml_x", true); 
          group.set_attribute<bool>("pml_y", false); 
          group.set_attribute<double>("pml_width_x", 1.8e-6); 
          group.set_attribute<double>("pml_width_y", 5e-7); 
          group.set_attribute<double>("pml_coeff", 1.6); 
          group.set_attribute<unsigned int>("pml_coeff_degree", 2); 

          group.set_attribute<double>("center_frequency", 20e9); 
          group.set_attribute<double>("frequency_range", 0.5e9); 
          group.set_attribute<double>( 
            "start_frequency", 
            group.get_attribute<double>("center_frequency") - 
              group.get_attribute<double>("frequency_range") / 2); 
          group.set_attribute<double>( 
            "stop_frequency", 
            group.get_attribute<double>("center_frequency") + 
              group.get_attribute<double>("frequency_range") / 2); 
          group.set_attribute<unsigned int>("nb_frequency_points", 400); 

          if (group_name == std::string("displacement")) 
            group.set_attribute<std::string>( 
              "simulation_name", std::string("phononic_cavity_displacement")); 
          else 
            group.set_attribute<std::string>( 
              "simulation_name", std::string("phononic_cavity_calibration")); 

          group.set_attribute<bool>("save_vtu_files", false); 
        } 

      { 

// 位移模拟。参数从位移HDF5组中读取，结果保存在同一HDF5组中。

        auto                    displacement = data.open_group("displacement"); 
        step62::Parameters<dim> parameters(displacement); 

        step62::ElasticWave<dim> elastic_problem(parameters); 
        elastic_problem.run(); 
      } 

      { 

// 校准模拟。参数从校准HDF5组中读取，结果保存在同一HDF5组中。

        auto                    calibration = data.open_group("calibration"); 
        step62::Parameters<dim> parameters(calibration); 

        step62::ElasticWave<dim> elastic_problem(parameters); 
        elastic_problem.run(); 
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

