CCTest_file/step-28.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2009 - 2021 by the deal.II authors 
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
 * Author: Yaqi Wang, Texas A&M University, 2009, 2010 
 */ 


// @sect3{Include files}  

// 我们从一堆包含文件开始，这些文件已经在以前的教程程序中解释过了。一个新的文件是  <code>timer.h</code>  : 这是第一个使用Timer类的例子程序。Timer同时记录了经过的挂钟时间（即安装在墙上的时钟所测量的时间）和CPU时钟时间（当前进程在CPU上使用的时间）。我们将在下面使用一个Timer来测量每个网格细化周期所需的CPU时间。

#include <deal.II/base/timer.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/parameter_handler.h> 
#include <deal.II/base/thread_management.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparsity_pattern.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/grid_generator.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

#include <fstream> 
#include <iostream> 

// 我们使用下一个include文件来访问块向量，它为我们提供了一种方便的方式来管理所有能量组的解和右手向量。

#include <deal.II/lac/block_vector.h> 

// 这个文件是用来将解从一个网格转移到另一个不同的网格。我们在每次网格迭代后初始化解法时使用它。

#include <deal.II/numerics/solution_transfer.h> 

// 当在一个网格上定义的函数与另一个网格上定义的形状函数进行整合时，我们需要一个函数  @p get_finest_common_cells  (在介绍中已经讨论过)，它定义在以下头文件中。

#include <deal.II/grid/grid_tools.h> 

// 我们使用一个来自boost的小工具类来保存输出流的状态（见下面的 <code>run</code> 函数）。

#include <boost/io/ios_state.hpp> 

// 这里还有两个C++标准头，我们用它们来定义列表数据类型，以及微调我们生成的输出。

#include <list> 
#include <iomanip> 

// 最后一步和以前所有的程序一样。

namespace Step28 
{ 
  using namespace dealii; 
// @sect3{Material data}  

// 首先，我们需要定义一个类，为主类提供材料数据（包括扩散系数、清除截面、散射截面、裂变截面和裂变光谱）。

// 构造函数的参数决定了我们为多少个能量组设置了相关的表格。目前，这个程序只包括2个能量组的数据，但是一个更复杂的程序可能也能为更多的能量组初始化数据结构，这取决于在参数文件中选择了多少个能量组。

// 对于每个不同的系数类型，都有一个函数来返回该系数在特定能量组（或能量组的组合，如分布截面 $\chi_g\nu\Sigma_{f,g'}$ 或散射截面 $\Sigma_{s,g'\to g}$ ）的值。除了能量组之外，这些系数还取决于燃料或控制棒的类型，正如介绍中所解释的那样。因此，这些函数需要一个额外的参数， @p  material_id，以确定特定种类的棒。在这个程序中，我们使用 <code>n_materials=8</code> 不同种类的棒子。

// 除了散射截面，每个系数都可以表示为一个二维浮点数组中的一个条目，该数组由能量组编号以及材料ID索引。表类模板是存储此类数据的理想方式。最后，散射系数取决于两个能量组的索引，因此需要存储在一个三维数组中，为此我们再次使用表类，这时第一个模板参数（表示数组的维度）当然需要是三。

  class MaterialData 
  { 
  public: 
    MaterialData(const unsigned int n_groups); 

    double get_diffusion_coefficient(const unsigned int group, 
                                     const unsigned int material_id) const; 
    double get_removal_XS(const unsigned int group, 
                          const unsigned int material_id) const; 
    double get_fission_XS(const unsigned int group, 
                          const unsigned int material_id) const; 
    double get_fission_dist_XS(const unsigned int group_1, 
                               const unsigned int group_2, 
                               const unsigned int material_id) const; 
    double get_scattering_XS(const unsigned int group_1, 
                             const unsigned int group_2, 
                             const unsigned int material_id) const; 
    double get_fission_spectrum(const unsigned int group, 
                                const unsigned int material_id) const; 

  private: 
    const unsigned int n_groups; 
    const unsigned int n_materials; 

    Table<2, double> diffusion; 
    Table<2, double> sigma_r; 
    Table<2, double> nu_sigma_f; 
    Table<3, double> sigma_s; 
    Table<2, double> chi; 
  }; 

// 该类的构造函数用于初始化所有材料数据数组。它需要能量组的数量作为参数（如果该值不等于2，就会抛出一个错误，因为目前只实现了两个能量组的数据；但是，使用这个参数，该函数仍然是灵活的，可以扩展到未来）。在开始的成员初始化部分，它也将数组的大小调整为正确的大小。

// 目前，材料数据被存储为8种不同类型的材料。这一点在将来也可以很容易地被扩展。

  MaterialData::MaterialData(const unsigned int n_groups) 
    : n_groups(n_groups) 
    , n_materials(8) 
    , diffusion(n_materials, n_groups) 
    , sigma_r(n_materials, n_groups) 
    , nu_sigma_f(n_materials, n_groups) 
    , sigma_s(n_materials, n_groups, n_groups) 
    , chi(n_materials, n_groups) 
  { 
    switch (this->n_groups) 
      { 
        case 2: 
          { 
            for (unsigned int m = 0; m < n_materials; ++m) 
              { 
                diffusion[m][0] = 1.2; 
                diffusion[m][1] = 0.4; 
                chi[m][0]       = 1.0; 
                chi[m][1]       = 0.0; 
                sigma_r[m][0]   = 0.03; 
                for (unsigned int group_1 = 0; group_1 < n_groups; ++group_1) 
                  for (unsigned int group_2 = 0; group_2 < n_groups; ++group_2) 
                    sigma_s[m][group_1][group_2] = 0.0; 
              } 

            diffusion[5][1] = 0.2; 

            sigma_r[4][0] = 0.026; 
            sigma_r[5][0] = 0.051; 
            sigma_r[6][0] = 0.026; 
            sigma_r[7][0] = 0.050; 

            sigma_r[0][1] = 0.100; 
            sigma_r[1][1] = 0.200; 
            sigma_r[2][1] = 0.250; 
            sigma_r[3][1] = 0.300; 
            sigma_r[4][1] = 0.020; 
            sigma_r[5][1] = 0.040; 
            sigma_r[6][1] = 0.020; 
            sigma_r[7][1] = 0.800; 

            nu_sigma_f[0][0] = 0.0050; 
            nu_sigma_f[1][0] = 0.0075; 
            nu_sigma_f[2][0] = 0.0075; 
            nu_sigma_f[3][0] = 0.0075; 
            nu_sigma_f[4][0] = 0.000; 
            nu_sigma_f[5][0] = 0.000; 
            nu_sigma_f[6][0] = 1e-7; 
            nu_sigma_f[7][0] = 0.00; 

            nu_sigma_f[0][1] = 0.125; 
            nu_sigma_f[1][1] = 0.300; 
            nu_sigma_f[2][1] = 0.375; 
            nu_sigma_f[3][1] = 0.450; 
            nu_sigma_f[4][1] = 0.000; 
            nu_sigma_f[5][1] = 0.000; 
            nu_sigma_f[6][1] = 3e-6; 
            nu_sigma_f[7][1] = 0.00; 

            sigma_s[0][0][1] = 0.020; 
            sigma_s[1][0][1] = 0.015; 
            sigma_s[2][0][1] = 0.015; 
            sigma_s[3][0][1] = 0.015; 
            sigma_s[4][0][1] = 0.025; 
            sigma_s[5][0][1] = 0.050; 
            sigma_s[6][0][1] = 0.025; 
            sigma_s[7][0][1] = 0.010; 

            break; 
          } 

        default: 
          Assert(false, 
                 ExcMessage( 
                   "Presently, only data for 2 groups is implemented")); 
      } 
  } 

// 接下来是返回给定材料和能量组的系数值的函数。它们所做的就是确保给定的参数在允许的范围内，然后在相应的表格中查找相应的值。

  double 
  MaterialData::get_diffusion_coefficient(const unsigned int group, 
                                          const unsigned int material_id) const 
  { 
    Assert(group < n_groups, ExcIndexRange(group, 0, n_groups)); 
    Assert(material_id < n_materials, 
           ExcIndexRange(material_id, 0, n_materials)); 

    return diffusion[material_id][group]; 
  } 

  double MaterialData::get_removal_XS(const unsigned int group, 
                                      const unsigned int material_id) const 
  { 
    Assert(group < n_groups, ExcIndexRange(group, 0, n_groups)); 
    Assert(material_id < n_materials, 
           ExcIndexRange(material_id, 0, n_materials)); 

    return sigma_r[material_id][group]; 
  } 

  double MaterialData::get_fission_XS(const unsigned int group, 
                                      const unsigned int material_id) const 
  { 
    Assert(group < n_groups, ExcIndexRange(group, 0, n_groups)); 
    Assert(material_id < n_materials, 
           ExcIndexRange(material_id, 0, n_materials)); 

    return nu_sigma_f[material_id][group]; 
  } 

  double MaterialData::get_scattering_XS(const unsigned int group_1, 
                                         const unsigned int group_2, 
                                         const unsigned int material_id) const 
  { 
    Assert(group_1 < n_groups, ExcIndexRange(group_1, 0, n_groups)); 
    Assert(group_2 < n_groups, ExcIndexRange(group_2, 0, n_groups)); 
    Assert(material_id < n_materials, 
           ExcIndexRange(material_id, 0, n_materials)); 

    return sigma_s[material_id][group_1][group_2]; 
  } 

  double 
  MaterialData::get_fission_spectrum(const unsigned int group, 
                                     const unsigned int material_id) const 
  { 
    Assert(group < n_groups, ExcIndexRange(group, 0, n_groups)); 
    Assert(material_id < n_materials, 
           ExcIndexRange(material_id, 0, n_materials)); 

    return chi[material_id][group]; 
  } 

// 计算裂变分布截面的函数略有不同，因为它将其值计算为另外两个系数的乘积。我们不需要在这里检查参数，因为这在我们调用其他两个相关函数时已经发生了，尽管这样做可能也无妨。

  double MaterialData::get_fission_dist_XS(const unsigned int group_1, 
                                           const unsigned int group_2, 
                                           const unsigned int material_id) const 
  { 
    return (get_fission_spectrum(group_1, material_id) * 
            get_fission_XS(group_2, material_id)); 
  } 

//  @sect3{The <code>EnergyGroup</code> class}  

// 第一个有趣的类是包含所有特定于单个能量组的东西。为了将那些属于同一个对象的东西分组，我们声明了一个结构，该结构包含了用于单个能量组的网格的Triangulation和DoFHandler对象，以及一些其他对象和成员函数，我们将在下面的章节中讨论。

// 这个类的主要原因如下：对于正向问题（有指定的右手边）和特征值问题，人们通常要解决一连串的问题，而不是完全耦合的问题，每个能量组。一旦意识到单一能量组的系统矩阵是对称和正定的（它只是一个扩散算子），而完全耦合问题的矩阵通常是非对称和非定值的，这就可以理解了。如果涉及几个以上的能量组，它也是非常大和相当完整的。

// 让我们先看看在有外部右手的情况下要解决的方程（对于时间无关的情况）。
// @f{eqnarray*} -\nabla
//  \cdot(D_g(x) \nabla \phi_g(x)) + \Sigma_{r,g}(x)\phi_g(x) =
//  \chi_g\sum_{g'=1}^G\nu\Sigma_{f,g'}(x)\phi_{g'}(x) + \sum_{g'\ne
//  g}\Sigma_{s,g'\to g}(x)\phi_{g'}(x) + s_{\mathrm{ext},g}(x) 
//  @f}

// 我们通常会通过将右手边的所有项与  $g'=g$  移到左手边来解决这个方程，并求出  $\phi_g$  。当然，我们还不知道 $\phi_{g'}$ ，因为这些变量的方程包括涉及 $\phi_g$ 的右侧项。在这种情况下，通常的做法是进行迭代：计算 
// @f{eqnarray*} -\nabla \cdot(D_g(x) \nabla \phi^{(n)}_g(x)) &+&
//  \Sigma_{r,g}(x)\phi^{(n)}_g(x) \\ &=&
//  \chi_g\sum_{g'=1}^{g-1}\nu\Sigma_{f,g'}(x)\phi^{(n)}_{g'}(x) +
//  \chi_g\sum_{g'=g}^G\nu\Sigma_{f,g'}(x)\phi^{(n-1)}_{g'}(x) + \sum_{g'\ne
//  g, g'<g}\Sigma_{s,g'\to g}(x)\phi^{(n)}_{g'}(x) + \sum_{g'\ne g,
//  g'>g}\Sigma_{s,g'\to g}(x)\phi^{(n-1)}_{g'}(x) + s_{\mathrm{ext},g}(x)
//  @f} 。

// 换句话说，我们一个一个地解方程，如果 $g'\ge g$ ，就用上一次迭代的 $\phi_{g'}$ 的值，如果 $g'<g$ ，就用本次迭代已经计算的 $\phi_{g'}$ 的值。

// 在计算特征值时，我们做了一个非常类似的迭代，只是我们没有外部的右手边，而且每次迭代后的解都会被缩放，正如在介绍中所解释的。

// 在任何一种情况下，如果我们所做的只是让下面这一类人具备这些能力，那么这两种情况就可以共同处理。(i) 形成左手边的矩阵，(ii) 形成组内右手边的贡献，即涉及不相干的来源，(iii) 形成源于组  $g'$  的对右手边的贡献。这个类正是做这些工作（以及一些簿记工作，如网格细化、设置矩阵和向量等）。另一方面，这个类本身并不知道有多少个能量组，特别是它们之间的相互作用，也就是说，外部迭代的样子（以及因此我们是解决一个特征值还是一个直接问题）的决定是留给本程序下面的NeutronDiffusionProblem类。

// 所以让我们来看看这个类和它的接口。

  template <int dim> 
  class EnergyGroup 
  { 
  public: 
// @sect5{<code>EnergyGroup</code> public member functions}  

// 该类有相当数量的公共成员函数，因为其操作方式是由外部控制的，因此所有做重要事情的函数都需要从另一个类中调用。让我们从记账开始：该类显然需要知道它所代表的能量组，使用哪些材料数据，以及从哪个粗略的网格开始。构造函数接收这些信息，并通过这些信息初始化相关的成员变量（见下文）。

// 然后，我们还需要设置线性系统的函数，即在给定的有限元对象的情况下，正确地确定矩阵的大小和它的稀疏模式等。 <code>setup_linear_system</code> 函数就是这样做的。最后，对于这个初始块，有两个函数可以返回这个对象中使用的活动单元和自由度的数量--利用这一点，我们可以使三角形和DoF处理成员变量成为私有的，不必授予外部使用，增强了封装性。

    EnergyGroup(const unsigned int        group, 
                const MaterialData &      material_data, 
                const Triangulation<dim> &coarse_grid, 
                const FiniteElement<dim> &fe); 

    void setup_linear_system(); 

    unsigned int n_active_cells() const; 
    unsigned int n_dofs() const; 

// 然后是为每个迭代和当前能量组组装线性系统的函数。请注意，该矩阵与迭代次数无关，因此在每个细化周期只需计算一次。对于必须在每次逆功率迭代中更新的右手边来说，情况就有点复杂了，而且由于计算它可能涉及到几个不同的网格，正如介绍中所解释的那样，这就更复杂了。为了使事情在解决正向或特征值问题方面更加灵活，我们将右手边的计算分成一个函数，将无关的源和组内贡献（我们将其称为零函数，作为特征值问题的源项）和一个计算来自另一个能量组的右手边的贡献。

    void assemble_system_matrix(); 
    void assemble_ingroup_rhs(const Function<dim> &extraneous_source); 
    void assemble_cross_group_rhs(const EnergyGroup<dim> &g_prime); 

// 接下来我们需要一组函数来实际计算线性系统的解，并对其进行处理（比如计算介绍中提到的裂变源贡献，将图形信息写入输出文件，计算误差指标，或者根据这些标准和阈值实际细化和粗化网格）。所有这些函数以后都可以从驱动类 <code>NeutronDiffusionProblem</code> 中调用，或者你想实现的任何其他类来解决涉及中子通量方程的问题。

    void solve(); 

    double get_fission_source() const; 

    void output_results(const unsigned int cycle) const; 

    void estimate_errors(Vector<float> &error_indicators) const; 

    void refine_grid(const Vector<float> &error_indicators, 
                     const double         refine_threshold, 
                     const double         coarsen_threshold); 
// @sect5{<code>EnergyGroup</code> public data members}  

// 作为面向对象编程的良好实践，我们通过使它们成为私有的来隐藏大多数数据成员。然而，我们必须允许驱动进程的类访问解向量以及上一次迭代的解，因为在幂迭代中，解向量在每次迭代中都被我们正在寻找的特征值的当前猜测所缩放。

  public: 
    Vector<double> solution; 
    Vector<double> solution_old; 
// @sect5{<code>EnergyGroup</code> private data members}  

// 其余的数据成员是私有的。与之前所有的教程程序相比，唯一的新数据成员是一个存储此对象所代表的能量组的整数，以及此对象的构造函数从驱动类中得到的材料数据对象的引用。同样地，构造函数得到了我们要使用的有限元对象的引用。

// 最后，我们必须在每次迭代中对线性系统应用边界值，即相当频繁。我们不是每次都插值，而是在每个新的网格上插值一次，然后和这个类的所有其他数据一起存储。

  private: 
    const unsigned int  group; 
    const MaterialData &material_data; 

    Triangulation<dim>        triangulation; 
    const FiniteElement<dim> &fe; 
    DoFHandler<dim>           dof_handler; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> system_rhs; 

    std::map<types::global_dof_index, double> boundary_values; 
    AffineConstraints<double>                 hanging_node_constraints; 
// @sect5{<code>EnergyGroup</code> private member functions}  

// 在这个类中有一个私有成员函数。它递归地走过两个网格的单元，以计算跨组的右手边项。这个算法在本程序的介绍中已经解释过了。这个函数的参数是对一个对象的引用，该对象代表了我们要整合的右手项的能量组，一个指向用于当前能量组的网格单元的迭代器，一个指向另一个网格上相应单元的迭代器，以及将自由度从两个单元中较粗的单元插补到较细的单元的矩阵。

  private: 
    void assemble_cross_group_rhs_recursive( 
      const EnergyGroup<dim> &                       g_prime, 
      const typename DoFHandler<dim>::cell_iterator &cell_g, 
      const typename DoFHandler<dim>::cell_iterator &cell_g_prime, 
      const FullMatrix<double> &                     prolongation_matrix); 
  }; 
// @sect4{Implementation of the <code>EnergyGroup</code> class}  

// 这个类的前几个函数大部分是不言自明的。构造函数只设置了几个数据成员，并创建了一个给定三角形的副本，作为该能量组使用的三角形的基础。接下来的两个函数只是从私有数据成员中返回数据，从而使我们能够使这些数据成员私有化。

  template <int dim> 
  EnergyGroup<dim>::EnergyGroup(const unsigned int        group, 
                                const MaterialData &      material_data, 
                                const Triangulation<dim> &coarse_grid, 
                                const FiniteElement<dim> &fe) 
    : group(group) 
    , material_data(material_data) 
    , fe(fe) 
    , dof_handler(triangulation) 
  { 
    triangulation.copy_triangulation(coarse_grid); 
    dof_handler.distribute_dofs(fe); 
  } 

  template <int dim> 
  unsigned int EnergyGroup<dim>::n_active_cells() const 
  { 
    return triangulation.n_active_cells(); 
  } 

  template <int dim> 
  unsigned int EnergyGroup<dim>::n_dofs() const 
  { 
    return dof_handler.n_dofs(); 
  } 

//  @sect5{<code>EnergyGroup::setup_linear_system</code>}  

// 第一个 "真正的 "函数是在新的网格上或网格细化后设置网格、矩阵等。我们用这个函数来初始化稀疏系统矩阵，以及右手边的向量。如果求解向量之前从未被设置过（如用零大小表示），我们也会初始化它并将其设置为默认值。如果它已经有一个非零的大小，我们就不这么做了（也就是说，这个函数是在网格细化之后调用的），因为在这种情况下，我们希望在不同的网格细化中保留解决方案（这一点我们在 <code>EnergyGroup::refine_grid</code> 函数中做到了）。

  template <int dim> 
  void EnergyGroup<dim>::setup_linear_system() 
  { 
    const unsigned int n_dofs = dof_handler.n_dofs(); 

    hanging_node_constraints.clear(); 
    DoFTools::make_hanging_node_constraints(dof_handler, 
                                            hanging_node_constraints); 
    hanging_node_constraints.close(); 

    system_matrix.clear(); 

    DynamicSparsityPattern dsp(n_dofs, n_dofs); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 
    hanging_node_constraints.condense(dsp); 
    sparsity_pattern.copy_from(dsp); 

    system_matrix.reinit(sparsity_pattern); 

    system_rhs.reinit(n_dofs); 

    if (solution.size() == 0) 
      { 
        solution.reinit(n_dofs); 
        solution_old.reinit(n_dofs); 
        solution_old = 1.0; 
        solution     = solution_old; 
      } 

// 在这个函数的最后，我们更新边界节点列表和它们的数值，首先清除这个列表和重新插值的边界数值（记住，这个函数是在第一次设置网格后，以及每次网格细化后调用）。

// 为了理解这段代码，有必要认识到我们使用 <code>GridGenerator::subdivided_hyper_rectangle</code> 函数来创建网格(在 <code>NeutronDiffusionProblem::initialize_problem</code> )，其中我们将最后一个参数设置为 <code>true</code>  。这意味着域的边界被 "着色"，也就是说，域的四个（或六个，在3D）边被分配了不同的边界指标。结果是，底部边界得到指标0，顶部一个边界指标1，左右边界分别得到指标2和3。

// 在这个程序中，我们只模拟一个，即右上角的反应器的四分之一。也就是说，我们只想在顶部和右侧边界插值边界条件，而在底部和左侧边界不做任何事情（即施加自然的、无流量的诺伊曼边界条件）。这很容易被推广到任意维度，即我们想在指标为1、3、......的边界上插值，我们在下面的循环中这样做（注意，对 <code>VectorTools::interpolate_boundary_values</code> 的调用是加法的，即它们不首先清除边界值图）。

    boundary_values.clear(); 

    for (unsigned int i = 0; i < dim; ++i) 
      VectorTools::interpolate_boundary_values(dof_handler, 
                                               2 * i + 1, 
                                               Functions::ZeroFunction<dim>(), 
                                               boundary_values); 
  } 

//  @sect5{<code>EnergyGroup::assemble_system_matrix</code>}  

// 接下来我们需要函数来组装系统矩阵和右手边。考虑到介绍中列出的方程以及我们在以前的例子程序中看到的内容，组装矩阵是很简单的。注意使用 <code>cell->material_id()</code> 来获取一个单元的材料种类。还请注意我们是如何设置正交公式的顺序的，以便它总是适合于正在使用的有限元。

// 最后，请注意，由于我们在这里只组装了系统矩阵，所以我们还不能消除边界值（我们需要右边的向量来实现）。我们将其推迟到 <code>EnergyGroup::solve</code> 函数中，这时所有的信息都可以得到。

  template <int dim> 
  void EnergyGroup<dim>::assemble_system_matrix() 
  { 
    const QGauss<dim> quadrature_formula(fe.degree + 1); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_matrix = 0; 

        fe_values.reinit(cell); 

        const double diffusion_coefficient = 
          material_data.get_diffusion_coefficient(group, cell->material_id()); 
        const double removal_XS = 
          material_data.get_removal_XS(group, cell->material_id()); 

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            for (unsigned int j = 0; j < dofs_per_cell; ++j) 
              cell_matrix(i, j) += 
                ((diffusion_coefficient * fe_values.shape_grad(i, q_point) * 
                    fe_values.shape_grad(j, q_point) + 
                  removal_XS * fe_values.shape_value(i, q_point) * 
                    fe_values.shape_value(j, q_point)) * 
                 fe_values.JxW(q_point)); 

        cell->get_dof_indices(local_dof_indices); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
            system_matrix.add(local_dof_indices[i], 
                              local_dof_indices[j], 
                              cell_matrix(i, j)); 
      } 

    hanging_node_constraints.condense(system_matrix); 
  } 

//  @sect5{<code>EnergyGroup::assemble_ingroup_rhs</code>}  

// 正如 <code>EnergyGroup</code> 类的文档中所解释的，我们将组装右手边分成两部分：组内耦合和跨组耦合。首先，我们需要一个函数来组装这里的一个特定组的右手边，即包括一个无关的源（我们将在特征值问题上设置为零）以及组内裂变贡献。 组内散射已经在清除截面的定义中得到了考虑）。这个函数的工作原理就组装右手边而言是非常标准的，因此不需要更多的评论，只是我们要提到在函数的开头将右手边的向量设置为零--这一点我们不打算为跨组项做，这些跨组项只是添加到右手边的向量。

  template <int dim> 
  void 
  EnergyGroup<dim>::assemble_ingroup_rhs(const Function<dim> &extraneous_source) 
  { 
    system_rhs.reinit(dof_handler.n_dofs()); 

    const QGauss<dim> quadrature_formula(fe.degree + 1); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_quadrature_points | 
                              update_JxW_values); 

    Vector<double>      cell_rhs(dofs_per_cell); 
    std::vector<double> extraneous_source_values(n_q_points); 
    std::vector<double> solution_old_values(n_q_points); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_rhs = 0; 

        fe_values.reinit(cell); 

        const double fission_dist_XS = 
          material_data.get_fission_dist_XS(group, group, cell->material_id()); 

        extraneous_source.value_list(fe_values.get_quadrature_points(), 
                                     extraneous_source_values); 

        fe_values.get_function_values(solution_old, solution_old_values); 

        cell->get_dof_indices(local_dof_indices); 

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            cell_rhs(i) += 
              ((extraneous_source_values[q_point] + 
                fission_dist_XS * solution_old_values[q_point]) * 
               fe_values.shape_value(i, q_point) * fe_values.JxW(q_point)); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          system_rhs(local_dof_indices[i]) += cell_rhs(i); 
      } 
  } 

//  @sect5{<code>EnergyGroup::assemble_cross_group_rhs</code>}  

// 对于组装单一能量组方程的右手向量来说，更有趣的函数是将能量组  $g$  和  $g'$  耦合起来。正如介绍中所解释的，我们首先要找到两个能量组的网格所共有的单元集。首先我们调用 <code>get_finest_common_cells</code> 来获得这一对来自两个网格的共同单元的列表。一对单元格中的两个单元格可能都不活跃，但至少有一个是活跃的。然后我们将这些单元对中的每一个交给一个函数，该函数将递归地计算右手边的项。

// 注意，组内耦合在之前已经处理过了，所以我们提前退出这个函数，如果  $g=g'$  。

  template <int dim> 
  void 
  EnergyGroup<dim>::assemble_cross_group_rhs(const EnergyGroup<dim> &g_prime) 
  { 
    if (group == g_prime.group) 
      return; 

    const std::list<std::pair<typename DoFHandler<dim>::cell_iterator, 
                              typename DoFHandler<dim>::cell_iterator>> 
      cell_list = 
        GridTools::get_finest_common_cells(dof_handler, g_prime.dof_handler); 

    for (const auto &cell_pair : cell_list) 
      { 
        FullMatrix<double> unit_matrix(fe.n_dofs_per_cell()); 
        for (unsigned int i = 0; i < unit_matrix.m(); ++i) 
          unit_matrix(i, i) = 1; 
        assemble_cross_group_rhs_recursive(g_prime, 
                                           cell_pair.first, 
                                           cell_pair.second, 
                                           unit_matrix); 
      } 
  } 

//  @sect5{<code>EnergyGroup::assemble_cross_group_rhs_recursive</code>}  

// 这是最后一个处理在潜在的不同网格上递归组装右手边条款的函数，使用介绍中描述的算法。该函数接收一个代表能量组 $g'$ 的对象的引用，以及能量组 $g$ 和 $g'$ 的网格中相应单元的迭代器。起初，即从上面的函数中调用这个函数时，这两个单元将是两个网格上的匹配单元；然而，这两个单元中的一个可能被进一步细化，我们将递归地调用这个函数，两个迭代器中的一个被原始单元的一个子单元所取代。

// 最后一个参数是介绍中的矩阵乘积矩阵 $B_{c^{(k)}}^T \cdots B_{c'}^T B_c^T$ ，它从两个单元中较粗的单元插值到较细的单元。如果这两个单元格匹配，那么这就是身份矩阵--正是我们最初传递给这个函数的。

// 该函数必须考虑两种情况：两种单元格都没有进一步细化，即没有子代，在这种情况下，我们可以最终组装这对单元格的右侧贡献；以及两种单元格中的一种被进一步细化，在这种情况下，我们必须通过循环未被激活的单元格的子代来不断地进行循环。下面将讨论这两种情况。

  template <int dim> 
  void EnergyGroup<dim>::assemble_cross_group_rhs_recursive( 
    const EnergyGroup<dim> &                       g_prime, 
    const typename DoFHandler<dim>::cell_iterator &cell_g, 
    const typename DoFHandler<dim>::cell_iterator &cell_g_prime, 
    const FullMatrix<double> &                     prolongation_matrix) 
  { 

// 第一种情况是，两个单元格都没有进一步的细化。在这种情况下，我们可以组装相关条款（见介绍）。这涉及到在两个单元中较细的单元上组装质量矩阵（事实上有两个具有不同系数的质量矩阵，一个用于裂变分布截面 $\chi_g\nu\Sigma_{f,g'}$ ，一个用于散射截面 $\Sigma_{s,g'\to g}$ ）。这是直截了当的，但请注意我们如何通过查看两个单元的细化程度来确定哪个是更细的单元。

    if (!cell_g->has_children() && !cell_g_prime->has_children()) 
      { 
        const QGauss<dim>  quadrature_formula(fe.degree + 1); 
        const unsigned int n_q_points = quadrature_formula.size(); 

        FEValues<dim> fe_values(fe, 
                                quadrature_formula, 
                                update_values | update_JxW_values); 

        if (cell_g->level() > cell_g_prime->level()) 
          fe_values.reinit(cell_g); 
        else 
          fe_values.reinit(cell_g_prime); 

        const double fission_dist_XS = 
          material_data.get_fission_dist_XS(group, 
                                            g_prime.group, 
                                            cell_g_prime->material_id()); 

        const double scattering_XS = 
          material_data.get_scattering_XS(g_prime.group, 
                                          group, 
                                          cell_g_prime->material_id()); 

        FullMatrix<double> local_mass_matrix_f(fe.n_dofs_per_cell(), 
                                               fe.n_dofs_per_cell()); 
        FullMatrix<double> local_mass_matrix_g(fe.n_dofs_per_cell(), 
                                               fe.n_dofs_per_cell()); 

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i) 
            for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j) 
              { 
                local_mass_matrix_f(i, j) += 
                  (fission_dist_XS * fe_values.shape_value(i, q_point) * 
                   fe_values.shape_value(j, q_point) * fe_values.JxW(q_point)); 
                local_mass_matrix_g(i, j) += 
                  (scattering_XS * fe_values.shape_value(i, q_point) * 
                   fe_values.shape_value(j, q_point) * fe_values.JxW(q_point)); 
              } 

// 现在我们有了所有的插值（延长）矩阵以及局部质量矩阵，所以我们只需要根据两个单元中哪一个更细，形成积
//  @f[
//  F_i|_{K_{cc'\cdots c^{(k)}}} = [B_c B_{c'} \cdots B_{c^{(k)}}
//  M_{K_{cc'\cdots c^{(k)}}}]^{ij} \phi_{g'}^j, 
//  @f]

//  或

//  @f[
//  F_i|_{K_{cc'\cdots c^{(k)}}} = [(B_c B_{c'} \cdots B_{c^{(k)}}
//  M_{K_{cc'\cdots c^{(k)}}})^T]^{ij} \phi_{g'}^j, 
//  @f]
//  。我们使用 <code>vmult</code> 函数提供的矩阵-向量乘积，或者使用 <code>Tvmult</code> 与转置矩阵的乘积来完成。这样做之后，我们将结果转移到能量组的全局右侧向量中  $g$  。


        Vector<double> g_prime_new_values(fe.n_dofs_per_cell()); 
        Vector<double> g_prime_old_values(fe.n_dofs_per_cell()); 
        cell_g_prime->get_dof_values(g_prime.solution_old, g_prime_old_values); 
        cell_g_prime->get_dof_values(g_prime.solution, g_prime_new_values); 

        Vector<double> cell_rhs(fe.n_dofs_per_cell()); 
        Vector<double> tmp(fe.n_dofs_per_cell()); 

        if (cell_g->level() > cell_g_prime->level()) 
          { 
            prolongation_matrix.vmult(tmp, g_prime_old_values); 
            local_mass_matrix_f.vmult(cell_rhs, tmp); 

            prolongation_matrix.vmult(tmp, g_prime_new_values); 
            local_mass_matrix_g.vmult_add(cell_rhs, tmp); 
          } 
        else 
          { 
            local_mass_matrix_f.vmult(tmp, g_prime_old_values); 
            prolongation_matrix.Tvmult(cell_rhs, tmp); 

            local_mass_matrix_g.vmult(tmp, g_prime_new_values); 
            prolongation_matrix.Tvmult_add(cell_rhs, tmp); 
          } 

        std::vector<types::global_dof_index> local_dof_indices( 
          fe.n_dofs_per_cell()); 
        cell_g->get_dof_indices(local_dof_indices); 

        for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i) 
          system_rhs(local_dof_indices[i]) += cell_rhs(i); 
      } 

// 另一种情况是，两个单元中的一个被进一步细化。在这种情况下，我们必须在所有的子单元上循环，将现有的矩阵的插值（延长）乘以从现在的单元到其子单元的插值（使用矩阵-矩阵乘法函数 <code>mmult</code> ），然后将结果再次交给这个非常相同的函数，但将有子单元替换为其子单元之一。

    else 
      for (unsigned int child = 0; 
           child < GeometryInfo<dim>::max_children_per_cell; 
           ++child) 
        { 
          FullMatrix<double> new_matrix(fe.n_dofs_per_cell(), 
                                        fe.n_dofs_per_cell()); 
          fe.get_prolongation_matrix(child).mmult(new_matrix, 
                                                  prolongation_matrix); 

          if (cell_g->has_children()) 
            assemble_cross_group_rhs_recursive(g_prime, 
                                               cell_g->child(child), 
                                               cell_g_prime, 
                                               new_matrix); 
          else 
            assemble_cross_group_rhs_recursive(g_prime, 
                                               cell_g, 
                                               cell_g_prime->child(child), 
                                               new_matrix); 
        } 
  } 
// @sect5{<code>EnergyGroup::get_fission_source</code>}  

// 在（反）功率迭代中，我们使用综合裂变源来更新 $k$  -特征值。鉴于其定义，以下函数基本上是不言自明的。

  template <int dim> 
  double EnergyGroup<dim>::get_fission_source() const 
  { 
    const QGauss<dim>  quadrature_formula(fe.degree + 1); 
    const unsigned int n_q_points = quadrature_formula.size(); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_JxW_values); 

    std::vector<double> solution_values(n_q_points); 

    double fission_source = 0; 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_values.reinit(cell); 

        const double fission_XS = 
          material_data.get_fission_XS(group, cell->material_id()); 

        fe_values.get_function_values(solution, solution_values); 

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          fission_source += 
            (fission_XS * solution_values[q_point] * fe_values.JxW(q_point)); 
      } 

    return fission_source; 
  } 
// @sect5{<code>EnergyGroup::solve</code>}  

// 接下来是一个解决之前组装的线性系统的函数。事情基本是标准的，只是我们把应用边界值的时间推迟到了这里，因为在之前的所有函数中，我们还是在为右边的向量做加法。

  template <int dim> 
  void EnergyGroup<dim>::solve() 
  { 
    hanging_node_constraints.condense(system_rhs); 
    MatrixTools::apply_boundary_values(boundary_values, 
                                       system_matrix, 
                                       solution, 
                                       system_rhs); 

    SolverControl            solver_control(system_matrix.m(), 
                                 1e-12 * system_rhs.l2_norm()); 
    SolverCG<Vector<double>> cg(solver_control); 

    PreconditionSSOR<SparseMatrix<double>> preconditioner; 
    preconditioner.initialize(system_matrix, 1.2); 

    cg.solve(system_matrix, solution, system_rhs, preconditioner); 

    hanging_node_constraints.distribute(solution); 
  } 

//  @sect5{<code>EnergyGroup::estimate_errors</code>}  

// 网格细化被分成两个函数。第一个函数估计每个单元的误差，用解的大小对其进行归一化处理，并将其返回到作为参数的向量中。调用函数收集所有能量组的所有误差指标，并计算出细化和粗化单元的阈值。

  template <int dim> 
  void EnergyGroup<dim>::estimate_errors(Vector<float> &error_indicators) const 
  { 
    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      QGauss<dim - 1>(fe.degree + 1), 
      std::map<types::boundary_id, const Function<dim> *>(), 
      solution, 
      error_indicators); 
    error_indicators /= solution.linfty_norm(); 
  } 

//  @sect5{<code>EnergyGroup::refine_grid</code>}  

// 第二部分是细化网格，给定前一个函数中计算的误差指标和误差阈值，超过这个阈值的单元应被细化，低于这个阈值的单元应被粗化。注意，我们在这里没有使用 <code>GridRefinement</code> 中的任何函数，而是自己设置细化标志。

// 在设置完这些标志后，我们使用SolutionTransfer类将求解向量从旧网格转移到新网格。这里使用的程序在该类的文档中已有详细描述。

  template <int dim> 
  void EnergyGroup<dim>::refine_grid(const Vector<float> &error_indicators, 
                                     const double         refine_threshold, 
                                     const double         coarsen_threshold) 
  { 
    for (const auto &cell : triangulation.active_cell_iterators()) 
      if (error_indicators(cell->active_cell_index()) > refine_threshold) 
        cell->set_refine_flag(); 
      else if (error_indicators(cell->active_cell_index()) < coarsen_threshold) 
        cell->set_coarsen_flag(); 

    SolutionTransfer<dim> soltrans(dof_handler); 

    triangulation.prepare_coarsening_and_refinement(); 
    soltrans.prepare_for_coarsening_and_refinement(solution); 

    triangulation.execute_coarsening_and_refinement(); 
    dof_handler.distribute_dofs(fe); 
    setup_linear_system(); 

    solution.reinit(dof_handler.n_dofs()); 
    soltrans.interpolate(solution_old, solution); 

// 强制执行约束条件，使插值后的解决方案在新的网格上符合要求。

    hanging_node_constraints.distribute(solution); 

    solution_old.reinit(dof_handler.n_dofs()); 
    solution_old = solution; 
  } 
// @sect5{<code>EnergyGroup::output_results</code>}  

// 该类的最后一个函数在每次网格迭代后输出网格和解。这在以前已经展示过很多次了。唯一值得指出的是使用 <code>Utilities::int_to_string</code> 函数将一个整数转换为其字符串表示。该函数的第二个参数表示我们应使用多少个数字 -- 如果这个值大于1，那么数字将被填充前导零。

  template <int dim> 
  void EnergyGroup<dim>::output_results(const unsigned int cycle) const 
  { 
    const std::string filename = std::string("solution-") + 
                                 Utilities::int_to_string(group, 2) + "." + 
                                 Utilities::int_to_string(cycle, 2) + ".vtu"; 

    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "solution"); 
    data_out.build_patches(); 

    std::ofstream output(filename); 
    data_out.write_vtu(output); 
  } 

//  @sect3{The <code>NeutronDiffusionProblem</code> class template}  

// 这是程序的主类，并不是因为它实现了所有的功能（事实上，大部分的功能都在 <code>EnergyGroup</code> 类中实现了），而是因为它包含了决定什么时候计算的驱动算法。它主要是像其他许多教程程序中显示的那样，它有一个公共的 <code>run</code> 函数和私有函数来做其他的事情。在一些地方，我们必须为所有能源组做一些事情，在这种情况下，如果deal.II被配置为多线程，我们将为每个组启动任务，让这些事情并行运行。 关于并行化的策略，请看一下 @ref threads 模块。

// 与以前的例子程序最大的不同是，我们还声明了一个嵌套类，该类有成员变量，用于所有可在输入文件中传递给程序的运行时参数。现在，这些参数是能量组的数量、细化周期的数量、要使用的有限元的多项式程度以及用于确定反幂迭代何时收敛的公差。此外，我们有一个该类的构造函数，将所有这些值设置为默认值，还有一个函数 <code>declare_parameters</code> 向ParameterHandler类描述输入文件中接受哪些参数，还有一个函数 <code>get_parameters</code> 可以从ParameterHandler对象中提取这些参数的值。参见 step-29 ，了解另一个使用ParameterHandler的例子。

  template <int dim> 
  class NeutronDiffusionProblem 
  { 
  public: 
    class Parameters 
    { 
    public: 
      Parameters(); 

      static void declare_parameters(ParameterHandler &prm); 
      void        get_parameters(ParameterHandler &prm); 

      unsigned int n_groups; 
      unsigned int n_refinement_cycles; 

      unsigned int fe_degree; 

      double convergence_tolerance; 
    }; 

    NeutronDiffusionProblem(const Parameters &parameters); 

    void run(); 

  private: 
// @sect5{<code>NeutronDiffusionProblem</code> private member functions}  

// 这个类中没有那么多的成员函数，因为大部分的功能已经被移到了 <code>EnergyGroup</code> 类中，只是从这个类的 <code>run()</code> 成员函数中调用。留下来的成员函数有不言自明的名字。

    void initialize_problem(); 

    void refine_grid(); 

    double get_total_fission_source() const; 
// @sect5{<code>NeutronDiffusionProblem</code> private member variables}  

// 接下来，我们有几个成员变量。特别是，这些是（i）对参数对象的引用（由本程序的主函数拥有，并传递给本类的构造函数），（ii）描述输入文件中要求的能量组数量的材料参数的对象，以及（iii）所有能量组将使用的有限元。

    const Parameters & parameters; 
    const MaterialData material_data; 
    FE_Q<dim>          fe; 

// 此外，我们有(iv)目前迭代时计算的特征值的值。事实上，这是在所有能量组之间共享的解决方案的唯一部分--解决方案的所有其他部分，如中子通量是特定于一个或另一个能量组的，因此被存储在描述单一能量组的对象中。

    double k_eff; 

// 最后一个计算对象（v）是一个指向能量组对象的数组。当然，这个数组的长度等于参数文件中指定的能量组的数量。

    std::vector<std::unique_ptr<EnergyGroup<dim>>> energy_groups; 

// 最后(vi)我们有一个文件流，我们将把总结的输出保存到这个文件中。

    std::ofstream convergence_table_stream; 
  }; 
// @sect4{Implementation of the <code>Parameters</code> class}  

// 在继续实现外层类之前，我们必须实现参数结构的功能。这是很直接的，事实上，对于所有使用ParameterHandler功能的这类参数类来说，看起来都是一样的。因此，我们将不再对此进行评论。

  template <int dim> 
  NeutronDiffusionProblem<dim>::Parameters::Parameters() 
    : n_groups(2) 
    , n_refinement_cycles(5) 
    , fe_degree(2) 
    , convergence_tolerance(1e-12) 
  {} 

  template <int dim> 
  void NeutronDiffusionProblem<dim>::Parameters::declare_parameters( 
    ParameterHandler &prm) 
  { 
    prm.declare_entry("Number of energy groups", 
                      "2", 
                      Patterns::Integer(), 
                      "The number of energy different groups considered"); 
    prm.declare_entry("Refinement cycles", 
                      "5", 
                      Patterns::Integer(), 
                      "Number of refinement cycles to be performed"); 
    prm.declare_entry("Finite element degree", 
                      "2", 
                      Patterns::Integer(), 
                      "Polynomial degree of the finite element to be used"); 
    prm.declare_entry( 
      "Power iteration tolerance", 
      "1e-12", 
      Patterns::Double(), 
      "Inner power iterations are stopped when the change in k_eff falls " 
      "below this tolerance"); 
  } 

  template <int dim> 
  void NeutronDiffusionProblem<dim>::Parameters::get_parameters( 
    ParameterHandler &prm) 
  { 
    n_groups              = prm.get_integer("Number of energy groups"); 
    n_refinement_cycles   = prm.get_integer("Refinement cycles"); 
    fe_degree             = prm.get_integer("Finite element degree"); 
    convergence_tolerance = prm.get_double("Power iteration tolerance"); 
  } 

//  @sect4{Implementation of the <code>NeutronDiffusionProblem</code> class}  

// 现在是 <code>NeutronDiffusionProblem</code> 类。构造函数和析构函数没有什么值得注意的地方。

  template <int dim> 
  NeutronDiffusionProblem<dim>::NeutronDiffusionProblem( 
    const Parameters &parameters) 
    : parameters(parameters) 
    , material_data(parameters.n_groups) 
    , fe(parameters.fe_degree) 
    , k_eff(std::numeric_limits<double>::quiet_NaN()) 
  {} 

//  @sect5{<code>NeutronDiffusionProblem::initialize_problem</code>}  

// 第一个感兴趣的函数是设置反应堆核心的几何形状的函数。这在介绍中会有更详细的描述。

// 该函数的第一部分定义了几何数据，然后创建了一个粗略的网格，其单元数与我们模拟的那部分反应堆堆芯中的燃料棒（或针状单元）的数量相当。正如上面插值边界值时提到的， <code>GridGenerator::subdivided_hyper_rectangle</code> 函数的最后一个参数指定域的两侧应具有唯一的边界指标，这将使我们能够以简单的方式确定哪些边界具有诺伊曼条件，哪些边界具有迪里希特条件。

  template <int dim> 
  void NeutronDiffusionProblem<dim>::initialize_problem() 
  { 
    const unsigned int rods_per_assembly_x = 17, rods_per_assembly_y = 17; 
    const double       pin_pitch_x = 1.26, pin_pitch_y = 1.26; 
    const double       assembly_height = 200; 

    const unsigned int assemblies_x = 2, assemblies_y = 2, assemblies_z = 1; 

    const Point<dim> bottom_left = Point<dim>(); 
    const Point<dim> upper_right = 
      (dim == 2 ? Point<dim>(assemblies_x * rods_per_assembly_x * pin_pitch_x, 
                             assemblies_y * rods_per_assembly_y * pin_pitch_y) : 
                  Point<dim>(assemblies_x * rods_per_assembly_x * pin_pitch_x, 
                             assemblies_y * rods_per_assembly_y * pin_pitch_y, 
                             assemblies_z * assembly_height)); 

    std::vector<unsigned int> n_subdivisions; 
    n_subdivisions.push_back(assemblies_x * rods_per_assembly_x); 
    if (dim >= 2) 
      n_subdivisions.push_back(assemblies_y * rods_per_assembly_y); 
    if (dim >= 3) 
      n_subdivisions.push_back(assemblies_z); 

    Triangulation<dim> coarse_grid; 
    GridGenerator::subdivided_hyper_rectangle( 
      coarse_grid, n_subdivisions, bottom_left, upper_right, true); 

// 该函数的第二部分涉及每种类型的组件的引脚单元的材料数量。在这里，我们定义了四种不同类型的组件，对于这些组件，我们在以下表格中描述了燃料棒的排列。

// 这里描述的装配体来自于介绍中提到的基准，它们是（按照这个顺序）。  <ol>  
// <li>  'UX'组件。二氧化铀燃料组件，带有24个导向管和一个中央可移动裂变室  <li>  'UA' 组件。带有24个AIC的二氧化铀燃料组件和一个中央可移动裂变室  <li>  'PX'组件。MOX燃料组件，带有24个导向管和一个中央可移动裂变室  <li>  'R'组件：一个反射器。   </ol>  

// 注意这里列出的数字和从基准描述中提取的数字，以良好的老Fortran方式，是基于一的。我们以后在给各个单元分配材料时将从每个数字中减去1，以便将事情转换为C语言风格的零基索引。

    const unsigned int n_assemblies = 4; 
    const unsigned int assembly_materials 
      [n_assemblies][rods_per_assembly_x][rods_per_assembly_y] = { 
        {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 5, 1, 1, 5, 1, 1, 7, 1, 1, 5, 1, 1, 5, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}}, 
        {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 8, 1, 1, 8, 1, 1, 7, 1, 1, 8, 1, 1, 8, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
         {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}}, 
        {{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}, 
         {2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2}, 
         {2, 3, 3, 3, 3, 5, 3, 3, 5, 3, 3, 5, 3, 3, 3, 3, 2}, 
         {2, 3, 3, 5, 3, 4, 4, 4, 4, 4, 4, 4, 3, 5, 3, 3, 2}, 
         {2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2}, 
         {2, 3, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 3, 2}, 
         {2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2}, 
         {2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2}, 
         {2, 3, 5, 4, 4, 5, 4, 4, 7, 4, 4, 5, 4, 4, 5, 3, 2}, 
         {2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2}, 
         {2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2}, 
         {2, 3, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 3, 2}, 
         {2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2}, 
         {2, 3, 3, 5, 3, 4, 4, 4, 4, 4, 4, 4, 3, 5, 3, 3, 2}, 
         {2, 3, 3, 3, 3, 5, 3, 3, 5, 3, 3, 5, 3, 3, 3, 3, 2}, 
         {2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2}, 
         {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}}, 
        {{6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
         {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}}}; 

// 在描述了组成装配体的材料之后，我们必须指定装配体在核心中的排列。我们使用一个对称的模式，实际上只使用'UX'和'PX'装配体。

    const unsigned int core[assemblies_x][assemblies_y][assemblies_z] = { 
      {{0}, {2}}, {{2}, {0}}}; 

// 我们现在可以为每个单元实际设置材料ID。为此，我们在所有单元中循环，查看单元中心的位置，并确定它将在哪个组件和燃料棒中。我们增加了一些检查，以确保我们计算的位置在我们必须查找材料的数组的范围内）。在循环结束时，我们相应地设置材料标识符。

    for (auto &cell : coarse_grid.active_cell_iterators()) 
      { 
        const Point<dim> cell_center = cell->center(); 

        const unsigned int tmp_x = int(cell_center[0] / pin_pitch_x); 
        const unsigned int ax    = tmp_x / rods_per_assembly_x; 
        const unsigned int cx    = tmp_x - ax * rods_per_assembly_x; 

        const unsigned     tmp_y = int(cell_center[1] / pin_pitch_y); 
        const unsigned int ay    = tmp_y / rods_per_assembly_y; 
        const unsigned int cy    = tmp_y - ay * rods_per_assembly_y; 

        const unsigned int az = 
          (dim == 2 ? 0 : int(cell_center[dim - 1] / assembly_height)); 

        Assert(ax < assemblies_x, ExcInternalError()); 
        Assert(ay < assemblies_y, ExcInternalError()); 
        Assert(az < assemblies_z, ExcInternalError()); 

        Assert(core[ax][ay][az] < n_assemblies, ExcInternalError()); 

        Assert(cx < rods_per_assembly_x, ExcInternalError()); 
        Assert(cy < rods_per_assembly_y, ExcInternalError()); 

        cell->set_material_id(assembly_materials[core[ax][ay][az]][cx][cy] - 1); 
      } 

// 有了这样初始化的粗网格，我们创建适当数量的能量组对象，让它们用上面生成的粗网格初始化各自的网格。

    for (unsigned int group = 0; group < parameters.n_groups; ++group) 
      energy_groups.emplace_back(std::make_unique<EnergyGroup<dim>>( 
        group, material_data, coarse_grid, fe)); 
    convergence_table_stream.open("convergence_table"); 
    convergence_table_stream.precision(12); 
  } 
// @sect5{<code>NeutronDiffusionProblem::get_total_fission_source</code>}  

// 在特征值计算中，我们需要在每次功率迭代后计算裂变中子源总量。然后用总功率来更新k-effective。

// 由于总裂变源是所有能量组的总和，而且每个总和都可以独立计算，所以我们实际上是以并行方式进行的。其中一个问题是， <code>EnergyGroup</code> 类中计算裂变源的函数会返回一个值。我们想在循环本身中把这些值加在一起：理想的情况是，每个任务计算它的值，然后立即把它加到总数中。以这种方式对数值进行加总需要两个功能。  <ol>  
// <li>  我们需要一种存储数值的方式，使多个线程能够以防止数据竞赛的方式并发地读写（即线程安全的读写）。 </li>  
// <li>  我们需要一种方法来增加这样一个值，而且是线程安全的。 </li>  
// </ol>  

// 第一个特性可以通过模板类实现  <code>std::atomic</code>  。然而，第二个特性，由 <code>std::atomic<double>::fetch_add()</code> 实现，只在C++20及以后的版本中可用：由于deal.II支持旧版本的C++语言标准，我们还不能使用这个特性。因此，取而代之的是，我们简单地将每个组的值写成向量中的一个条目，并在函数的最后将这些值相加。

  template <int dim> 
  double NeutronDiffusionProblem<dim>::get_total_fission_source() const 
  { 
    std::vector<double>  fission_sources(parameters.n_groups); 
    Threads::TaskGroup<> tasks; 
    for (unsigned int group = 0; group < parameters.n_groups; ++group) 
      tasks += Threads::new_task<>([&, group]() { 
        fission_sources[group] = energy_groups[group]->get_fission_source(); 
      }); 
    tasks.join_all(); 

    return std::accumulate(fission_sources.begin(), fission_sources.end(), 0.0); 
  } 

//  @sect5{<code>NeutronDiffusionProblem::refine_grid</code>}  

// 下一个函数让各个能量组对象细化其网格。这其中的大部分，也是可以独立并行完成的任务：首先，让所有的能量组对象并行计算它们的误差指标，然后计算所有能量组的最大误差指标，并确定细化和粗化单元的阈值，然后要求所有的能量组相应地细化它们的网格，也是并行的。

  template <int dim> 
  void NeutronDiffusionProblem<dim>::refine_grid() 
  { 
    std::vector<types::global_dof_index> n_cells(parameters.n_groups); 
    for (unsigned int group = 0; group < parameters.n_groups; ++group) 
      n_cells[group] = energy_groups[group]->n_active_cells(); 

    BlockVector<float> group_error_indicators(n_cells); 

    { 
      Threads::TaskGroup<> tasks; 
      for (unsigned int group = 0; group < parameters.n_groups; ++group) 
        tasks += Threads::new_task([&, group]() { 
          energy_groups[group]->estimate_errors( 
            group_error_indicators.block(group)); 
        }); 
    } 

//  Threads::TaskGroup 的析构器连接所有线程，所以我们知道在我们退出范围时，计算已经完成。

    const float max_error         = group_error_indicators.linfty_norm(); 
    const float refine_threshold  = 0.3 * max_error; 
    const float coarsen_threshold = 0.01 * max_error; 

    { 
      Threads::TaskGroup<void> tasks; 
      for (unsigned int group = 0; group < parameters.n_groups; ++group) 
        tasks += Threads::new_task([&, group]() { 
          energy_groups[group]->refine_grid(group_error_indicators.block(group), 
                                            refine_threshold, 
                                            coarsen_threshold); 
        }); 
    } 
  } 
// @sect5{<code>NeutronDiffusionProblem::run</code>}  

// 最后，这就是肉的函数：在一连串的网格上进行迭代，并对每一个网格进行幂级迭代，以计算特征值。

// 鉴于介绍中对算法的描述，实际上没有什么可评论的。

  template <int dim> 
  void NeutronDiffusionProblem<dim>::run() 
  { 

// 我们希望只为这个函数改变输出精度，并在这个函数返回时恢复到 <code>std::cout</code> 的状态。因此，我们需要一种方法来撤销输出格式的改变。Boost提供了一种方便的方法来保存输出流的状态，并在当前块结束时（当调用 <code>restore_flags</code> 的析构器时）用 <code>ios_flags_saver</code> 类来恢复它，我们在这里使用了这种方法。

    boost::io::ios_flags_saver restore_flags(std::cout); 
    std::cout << std::setprecision(12) << std::fixed; 

// 我们通过k_eff的变化来计算下面的误差（即k_eff_old的区别。

    double k_eff_old = 0.0; 

    for (unsigned int cycle = 0; cycle < parameters.n_refinement_cycles; 
         ++cycle) 
      { 

// 我们将在下面测量每个周期所需的CPU时间。计时器的构造函数调用 Timer::start(), ，所以一旦我们创建了一个计时器，就可以查询它的信息。由于这个循环的许多部分是用任务并行化的，所以我们测量的CPU时间（如果我们用一个以上的线程运行）将大于墙的时间。

        Timer timer; 

        std::cout << "Cycle " << cycle << ':' << std::endl; 

        if (cycle == 0) 
          { 
            initialize_problem(); 
            for (unsigned int group = 0; group < parameters.n_groups; ++group) 
              energy_groups[group]->setup_linear_system(); 
          } 

        else 
          { 
            refine_grid(); 
            for (unsigned int group = 0; group < parameters.n_groups; ++group) 
              energy_groups[group]->solution *= k_eff; 
          } 

        std::cout << "   Numbers of active cells:       "; 
        for (unsigned int group = 0; group < parameters.n_groups; ++group) 
          std::cout << energy_groups[group]->n_active_cells() << ' '; 
        std::cout << std::endl; 
        std::cout << "   Numbers of degrees of freedom: "; 
        for (unsigned int group = 0; group < parameters.n_groups; ++group) 
          std::cout << energy_groups[group]->n_dofs() << ' '; 
        std::cout << std::endl << std::endl; 

        Threads::TaskGroup<> tasks; 
        for (unsigned int group = 0; group < parameters.n_groups; ++group) 
          tasks += Threads::new_task( 
            [&, group]() { energy_groups[group]->assemble_system_matrix(); }); 
        tasks.join_all(); 

        double       error; 
        unsigned int iteration = 1; 
        do 
          { 
            for (unsigned int group = 0; group < parameters.n_groups; ++group) 
              { 
                energy_groups[group]->assemble_ingroup_rhs( 
                  Functions::ZeroFunction<dim>()); 

                for (unsigned int bgroup = 0; bgroup < parameters.n_groups; 
                     ++bgroup) 
                  energy_groups[group]->assemble_cross_group_rhs( 
                    *energy_groups[bgroup]); 

                energy_groups[group]->solve(); 
              } 

            k_eff = get_total_fission_source(); 
            error = std::abs(k_eff - k_eff_old) / std::abs(k_eff); 
            const double flux_ratio = energy_groups[0]->solution.linfty_norm() / 
                                      energy_groups[1]->solution.linfty_norm(); 
            const double max_thermal = energy_groups[1]->solution.linfty_norm(); 
            std::cout << "Iter number:" << std::setw(2) << std::right 
                      << iteration << " k_eff=" << k_eff 
                      << " flux_ratio=" << flux_ratio 
                      << " max_thermal=" << max_thermal << std::endl; 
            k_eff_old = k_eff; 

            for (unsigned int group = 0; group < parameters.n_groups; ++group) 
              { 
                energy_groups[group]->solution_old = 
                  energy_groups[group]->solution; 
                energy_groups[group]->solution_old /= k_eff; 
              } 

            ++iteration; 
          } 
        while ((error > parameters.convergence_tolerance) && (iteration < 500)); 
        convergence_table_stream << cycle << " " << energy_groups[0]->n_dofs() 
                                 << " " << energy_groups[1]->n_dofs() << " " 
                                 << k_eff << " " 
                                 << energy_groups[0]->solution.linfty_norm() / 
                                      energy_groups[1]->solution.linfty_norm() 
                                 << '\n'; 

        for (unsigned int group = 0; group < parameters.n_groups; ++group) 
          energy_groups[group]->output_results(cycle); 

// 打印出关于模拟的信息以及耗费的CPU时间。我们可以不先调用 Timer::cpu_time() ，而直接调用 Timer::stop() ，以获得调用该函数时的已用CPU时间。

        std::cout << std::endl; 
        std::cout << "   Cycle=" << cycle << ", n_dofs=" 
                  << energy_groups[0]->n_dofs() + energy_groups[1]->n_dofs() 
                  << ",  k_eff=" << k_eff << ", time=" << timer.cpu_time() 
                  << std::endl; 

        std::cout << std::endl << std::endl; 
      } 
  } 
} // namespace Step28 

//  @sect3{The <code>main()</code> function}  

// 程序中的最后一件事在 <code>main()</code> 函数中。其结构与其他大多数教程程序一样，唯一的例外是我们在这里处理一个参数文件。 为此，我们首先看一下传递给这个函数的命令行参数：如果在命令行上没有指定输入文件，那么就使用 "project.prm"，否则就取命令行上作为第一个参数给出的文件名。

// 有了这个，我们创建一个ParameterHandler对象，让 <code>NeutronDiffusionProblem::Parameters</code> 类声明它想在输入文件中看到的所有参数（或者，采取默认值，如果参数文件中没有列出任何参数），然后读取输入文件，要求参数对象提取数值，最后把所有东西交给 <code>NeutronDiffusionProblem</code> 类型的对象来计算特征值。

int main(int argc, char **argv) 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step28; 

      std::string filename; 
      if (argc < 2) 
        filename = "project.prm"; 
      else 
        filename = argv[1]; 

      const unsigned int dim = 2; 

      ParameterHandler parameter_handler; 

      NeutronDiffusionProblem<dim>::Parameters parameters; 
      parameters.declare_parameters(parameter_handler); 

      parameter_handler.parse_input(filename); 

      parameters.get_parameters(parameter_handler); 

      NeutronDiffusionProblem<dim> neutron_diffusion_problem(parameters); 
      neutron_diffusion_problem.run(); 
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

