//include/deal.II-translator/numerics/matrix_tools_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_matrix_tools_h
#define dealii_matrix_tools_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>

#include <deal.II/lac/affine_constraints.h>

#include <map>

#ifdef DEAL_II_WITH_PETSC
#  include <petscsys.h>
#endif

DEAL_II_NAMESPACE_OPEN


// forward declarations
#ifndef DOXYGEN
template <int dim>
class Quadrature;


template <typename number>
class Vector;
template <typename number>
class FullMatrix;
template <typename number>
class SparseMatrix;

template <typename number>
class BlockSparseMatrix;
template <typename Number>
class BlockVector;

template <int dim, int spacedim>
class Mapping;
template <int dim, int spacedim>
class DoFHandler;

namespace hp
{
  template <int>
  class QCollection;
  template <int, int>
  class MappingCollection;
} // namespace hp


#  ifdef DEAL_II_WITH_PETSC
namespace PETScWrappers
{
  class MatrixBase;
  class VectorBase;
  namespace MPI
  {
    class BlockSparseMatrix;
    class BlockVector;
  } // namespace MPI
} // namespace PETScWrappers
#  endif

#  ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  class SparseMatrix;
  class BlockSparseMatrix;
  namespace MPI
  {
    class Vector;
    class BlockVector;
  } // namespace MPI
} // namespace TrilinosWrappers
#  endif
#endif


/**
 * 这个命名空间提供了一些函数，这些函数为一个给定的三角形，使用一个给定的有限元，一个给定的映射和一个正交公式，组装某些标准矩阵。
 *
 *  <h3>Conventions for all functions</h3>
 * 几乎所有的函数都有两个版本，一个需要一个明确的Mapping参数，另一个不需要。第二个版本一般会调用第一个版本的隐式
 * $Q_1$
 * 参数（即用一个MappingQGeneric(1)类型的参数）。如果你的代码打算使用不同于（双/三）线性的映射，那么你需要调用应该使用<b>with</b>映射参数的函数。
 * 所有的函数都需要一个稀疏矩阵对象来保存要创建的矩阵。这些函数假定矩阵被初始化为与给定自由度处理程序相对应的稀疏模式（SparsityPattern），也就是说，稀疏结构已经是需要的。你可以通过调用
 * DoFTools::make_sparsity_pattern() 函数做到这一点。
 * 此外，假定矩阵中没有相关数据。如果矩阵之前不是空的，有些条目会被覆盖，有些条目会包含无效的数据。因此，你可能想在组装前清除矩阵。
 * 默认情况下，所有创建的矩阵都是 "原始
 * "的：它们没有被压缩，也就是说，悬挂的节点没有被消除。原因是你可能想添加几个矩阵，然后可以在之后只压缩一次，而不是对每个矩阵都压缩。要实际对这些矩阵进行计算，你必须使用
 * AffineConstraints::condense
 * 函数对矩阵进行凝结；你还必须对右手边进行相应的凝结，并在事后分发解。另外，你可以给一个可选的参数AffineConstraints，它将单元格矩阵（和向量）条目用distribution_local_to_global写入全局矩阵和向量中。这样一来，从不同的来源添加几个矩阵就比较复杂了，你应该确保不混合应用约束的不同方式。当给定的AffineConstraints对象包含不均匀的约束时，需要特别小心谨慎。在这种情况下，这样组装的矩阵必须是唯一的矩阵（或者你需要为你生成的<b>every</b>矩阵组装<b>same</b>的右手边并加在一起）。
 * 如果你想在可能的AffineConstraints对象中使用边界条件，除了这个命名空间的函数生成的矩阵之外，你必须使用一个类似<tt>apply_boundary_values</tt>的函数，其中包括矩阵、解和右手边。
 *
 *  <h3>Supported matrices</h3>
 * 目前，有一些函数可以创建以下矩阵。  <ul>   <li>   @p create_mass_matrix:  通过数字正交创建带有条目 $m_{ij} =
 * \int_\Omega \phi_i(x) \phi_j(x) dx$ 的矩阵。这里， $\phi_i$
 * 是给出的有限元空间的基函数。
 * 可以给出一个系数来代替评估 $m_{ij} = \int_\Omega a(x)
 * \phi_i(x) \phi_j(x) dx$ 。
 * <li>   @p create_laplace_matrix: 通过数字正交法创建具有条目
 * $a_{ij} = \int_\Omega \nabla\phi_i(x) \nabla\phi_j(x) dx$ 的矩阵。
 * 同样，可以给一个系数来代替评估 $a_{ij} = \int_\Omega a(x)
 * \nabla\phi_i(x) \nabla\phi_j(x) dx$ 。  </ul>
 * 确保给这些函数的正交公式的阶数足够高，以便以要求的精度计算矩阵。对于这个正交规则的选择，你需要考虑到FiniteElement基础函数的多项式程度，系数
 * @p a, 的粗糙度以及给定 @p Mapping 的程度（如果有）。
 * 注意，对于矢量值元素，质量矩阵和拉普拉斯矩阵的实现方式是每个组件只与自己耦合，也就是说，属于不同组件的形状函数没有耦合。如果自由度已经根据其矢量分量进行了排序（例如，使用
 * DoFRenumbering::component_wise()),
 * ，那么产生的矩阵将是块状对角线。
 * 如果要建立质量矩阵或拉普拉斯矩阵的有限元有一个以上的分量，这些函数接受单一系数以及矢量值的系数函数。对于后一种情况，分量的数量必须与系统有限元的分量数量相吻合。
 *
 *  <h3>Matrices on the boundary</h3>
 * create_boundary_mass_matrix()创建条目为 $m_{ij} = \int_{\Gamma} \phi_i
 * \phi_j dx$ 的矩阵，其中 $\Gamma$
 * 是边界部分的联合，其指标包含在一个
 * std::map<types::boundary_id,  ] const
 * Function<spacedim,number>*>传递给函数（例如，如果你想为具有指标0和2的边界部分设置质量矩阵，你传递给函数一个键类型为
 * types::boundary_id 的地图，作为包含键0和2的参数 @p
 * boundary_functions
 * ）。矩阵的大小等于在边界上有支持的自由度的数量，也就是说，它
 * <em> 不是 </em>
 * 所有自由度的矩阵，而只是一个子集。公式中的 $\phi_i$
 * 是至少部分支持在 $\Gamma$
 * 上的基函数的子集）。为了确定要考虑哪些形状函数，以及为了确定哪个顺序，该函数需要一个
 * @p dof_to_boundary_mapping;
 * 这个对象将全局DoF数映射为位于边界上的自由度的编号，可以用函数
 * DoFTools::map_dof_to_boundary_indices(). 得到。
 * 为了工作，该函数需要一个正确大小的矩阵，建立在相应的稀疏模式之上。由于我们只在自由度的一个子集上工作，我们不能使用为整个自由度集创建的矩阵和稀疏模式。相反，你应该使用
 * DoFHandler::make_boundary_sparsity_pattern()
 * 函数来创建正确的稀疏模式，并在其上建立一个矩阵。
 * 请注意，目前没有任何函数可以计算 <em> 所有 </em>
 * 形状函数的质量矩阵，尽管这样一个函数的实现是微不足道的。
 *
 *  <h3>Right hand sides</h3>
 * 在许多情况下，你不仅要建立矩阵，还要建立一个右手边，这将给出一个具有
 * $f_i = \int_\Omega f(x) \phi_i(x) dx$
 * 的向量。为此，每个函数都有两个版本，一个只建立矩阵，一个也建立右手边的向量。如果你想创建一个右手边的向量而不创建矩阵，你可以使用
 * VectorTools::create_right_hand_side()
 * 函数。如果你想创建许多右手向量，使用后者可能很有用。
 *
 *
 * @ingroup numerics
 *
 *
 */
namespace MatrixCreator
{
  /**
   * 组装质量矩阵。如果没有给出系数（即，如果函数对象的指针是零，因为它是默认的），那么系数就被认为是常数，等于1。
   * 如果你想指定  @p constraints
   * 并使用默认的系数参数，你必须指定（未使用的）系数参数为
   * <code>(const Function<spacedim,number>const)nullptr</code>  。
   * 如果库被配置为使用多线程，这个函数就会并行工作。
   * 可选的参数 @p constraints
   * 允许直接在结果矩阵上应用约束。然而，请注意，当你有不均匀的约束，并且后来想添加几个这样的矩阵时，例如在时间相关的设置中，例如
   * step-26  的主循环，这就变得很困难。
   * 更多信息请参见该命名空间的一般文档。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(
    const Mapping<dim, spacedim> &          mapping,
    const DoFHandler<dim, spacedim> &       dof,
    const Quadrature<dim> &                 q,
    SparseMatrix<number> &                  matrix,
    const Function<spacedim, number> *const a    = nullptr,
    const AffineConstraints<number> &constraints = AffineConstraints<number>());

  /**
   * 调用create_mass_matrix()函数，见上文，使用<tt>mapping=MappingQGeneric
   * @<dim@>(1)</tt>.  。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(
    const DoFHandler<dim, spacedim> &       dof,
    const Quadrature<dim> &                 q,
    SparseMatrix<number> &                  matrix,
    const Function<spacedim, number> *const a    = nullptr,
    const AffineConstraints<number> &constraints = AffineConstraints<number>());

  /**
   * 组建质量矩阵和一个右手边的向量。如果没有给出系数（也就是说，如果函数对象的指针是零，因为它是默认的），那么系数被认为是常数，等于1。
   * 如果你想指定  @p constraints
   * 并使用默认的系数参数，你必须指定（未使用的）系数参数为
   * <code>(const Function <spacedim,number>const)nullptr</code>  。
   * 如果库被配置为使用多线程，这个函数就会并行工作。
   * 可选的参数 @p constraints
   * 允许直接在结果矩阵上应用约束。然而，请注意，当你有不均匀的约束，并且后来想添加几个这样的矩阵时，例如在时间相关的设置中，例如
   * step-26  的主循环，这就变得很困难。
   * 更多信息请参见该命名空间的一般文档。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(
    const Mapping<dim, spacedim> &          mapping,
    const DoFHandler<dim, spacedim> &       dof,
    const Quadrature<dim> &                 q,
    SparseMatrix<number> &                  matrix,
    const Function<spacedim, number> &      rhs,
    Vector<number> &                        rhs_vector,
    const Function<spacedim, number> *const a    = nullptr,
    const AffineConstraints<number> &constraints = AffineConstraints<number>());

  /**
   * 调用create_mass_matrix()函数，见上文，使用<tt>mapping=MappingQGeneric
   * @<dim@>(1)</tt>.  。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(
    const DoFHandler<dim, spacedim> &       dof,
    const Quadrature<dim> &                 q,
    SparseMatrix<number> &                  matrix,
    const Function<spacedim, number> &      rhs,
    Vector<number> &                        rhs_vector,
    const Function<spacedim, number> *const a    = nullptr,
    const AffineConstraints<number> &constraints = AffineConstraints<number>());

  /**
   * 与上面的函数相同，但用于hp-objects。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &           dof,
    const hp::QCollection<dim> &                q,
    SparseMatrix<number> &                      matrix,
    const Function<spacedim, number> *const     a = nullptr,
    const AffineConstraints<number> &constraints = AffineConstraints<number>());

  /**
   * 与上面的功能相同，但对hp-objects而言。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(
    const DoFHandler<dim, spacedim> &       dof,
    const hp::QCollection<dim> &            q,
    SparseMatrix<number> &                  matrix,
    const Function<spacedim, number> *const a    = nullptr,
    const AffineConstraints<number> &constraints = AffineConstraints<number>());

  /**
   * 与上述功能相同，但对hp-objects而言。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &           dof,
    const hp::QCollection<dim> &                q,
    SparseMatrix<number> &                      matrix,
    const Function<spacedim, number> &          rhs,
    Vector<number> &                            rhs_vector,
    const Function<spacedim, number> *const     a = nullptr,
    const AffineConstraints<number> &constraints = AffineConstraints<number>());

  /**
   * 与上面的功能相同，但对hp-objects而言。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  create_mass_matrix(
    const DoFHandler<dim, spacedim> &       dof,
    const hp::QCollection<dim> &            q,
    SparseMatrix<number> &                  matrix,
    const Function<spacedim, number> &      rhs,
    Vector<number> &                        rhs_vector,
    const Function<spacedim, number> *const a    = nullptr,
    const AffineConstraints<number> &constraints = AffineConstraints<number>());


  /**
   * 沿着边界组装质量矩阵和一个右手边的矢量。
   * 假设该矩阵已经被初始化为合适的稀疏模式（DoFHandler提供了一个合适的函数）。
   * 如果库被配置为使用多线程，这个函数就可以并行工作。
   * @arg   @p weight:
   * 用于计算质量矩阵的可选权重。如果没有给出权重，它将被设置为1。
   * 如果你想指定  @p component_mapping
   * 并使用默认的系数参数，你必须将（未使用的）系数参数指定为
   * <code>(const Function <spacedim,number>const)nullptr</code>  。
   * @arg   @p component_mapping:  如果 @p boundary_functions 和 @p dof
   * 中的成分不重合，这个向量允许它们被重新映射。如果这个向量不是空的，它必须有一个条目代表
   * @p 中的每个分量。这个条目是 @p boundary_functions
   * 中的分量编号，应该用于 @p dof. 中的这个分量。
   * 默认情况下，不应用重映射。      @todo
   * 这个函数对具有单元依赖形状函数的有限元不起作用。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  create_boundary_mass_matrix(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const Quadrature<dim - 1> &      q,
    SparseMatrix<number> &           matrix,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                     boundary_functions,
    Vector<number> &                        rhs_vector,
    std::vector<types::global_dof_index> &  dof_to_boundary_mapping,
    const Function<spacedim, number> *const weight            = 0,
    std::vector<unsigned int>               component_mapping = {});


  /**
   * 调用create_boundary_mass_matrix()函数，见上文，使用<tt>mapping=MappingQGeneric
   * @<dim@>(1)</tt>.  。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  create_boundary_mass_matrix(
    const DoFHandler<dim, spacedim> &dof,
    const Quadrature<dim - 1> &      q,
    SparseMatrix<number> &           matrix,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                     boundary_functions,
    Vector<number> &                        rhs_vector,
    std::vector<types::global_dof_index> &  dof_to_boundary_mapping,
    const Function<spacedim, number> *const a                 = nullptr,
    std::vector<unsigned int>               component_mapping = {});

  /**
   * 与上面的函数相同，但用于hp-objects。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  create_boundary_mass_matrix(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &           dof,
    const hp::QCollection<dim - 1> &            q,
    SparseMatrix<number> &                      matrix,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                     boundary_functions,
    Vector<number> &                        rhs_vector,
    std::vector<types::global_dof_index> &  dof_to_boundary_mapping,
    const Function<spacedim, number> *const a                 = nullptr,
    std::vector<unsigned int>               component_mapping = {});

  /**
   * 与上面的功能相同，但对hp-objects而言。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  create_boundary_mass_matrix(
    const DoFHandler<dim, spacedim> &dof,
    const hp::QCollection<dim - 1> & q,
    SparseMatrix<number> &           matrix,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                     boundary_functions,
    Vector<number> &                        rhs_vector,
    std::vector<types::global_dof_index> &  dof_to_boundary_mapping,
    const Function<spacedim, number> *const a                 = nullptr,
    std::vector<unsigned int>               component_mapping = {});

  /**
   * 组装拉普拉斯矩阵。如果没有给出系数（也就是说，如果函数对象的指针是零，因为它是默认的），那么系数就被认为是常数，等于1。
   * 如果你想指定  @p constraints
   * 并使用默认的系数参数，你必须指定（未使用的）系数参数为
   * <code>(const Function<spacedim>const)nullptr</code>  。
   * 如果库被配置为使用多线程，这个函数就会并行工作。
   * 可选的参数 @p constraints
   * 允许直接在结果矩阵上应用约束。然而，请注意，当你有不均匀的约束，并且后来想添加几个这样的矩阵时，例如在时间相关的设置中，例如
   * step-26  的主循环，这就变得困难了。
   * 更多信息请参见该命名空间的一般文档。
   *
   */
  template <int dim, int spacedim>
  void
  create_laplace_matrix(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const Quadrature<dim> &          q,
    SparseMatrix<double> &           matrix,
    const Function<spacedim> *const  a           = nullptr,
    const AffineConstraints<double> &constraints = AffineConstraints<double>());

  /**
   * 调用create_laplace_matrix()函数，见上文，使用<tt>mapping=MappingQGeneric
   * @<dim@>(1)</tt>.  。
   *
   */
  template <int dim, int spacedim>
  void
  create_laplace_matrix(
    const DoFHandler<dim, spacedim> &dof,
    const Quadrature<dim> &          q,
    SparseMatrix<double> &           matrix,
    const Function<spacedim> *const  a           = nullptr,
    const AffineConstraints<double> &constraints = AffineConstraints<double>());

  /**
   * 组建拉普拉斯矩阵和一个右手边的向量。如果没有给出系数，则假定其为常数一。
   * 如果你想指定  @p constraints
   * 并使用默认的系数参数，你必须指定（未使用的）系数参数为
   * <code>(const Function<spacedim>const)nullptr</code>  。
   * 如果库被配置为使用多线程，这个函数就会并行工作。
   * 可选的参数 @p constraints
   * 允许直接在结果矩阵上应用约束。然而，请注意，当你有不均匀的约束，并且后来想添加几个这样的矩阵时，例如在时间相关的设置中，例如
   * step-26  的主循环，这就变得很困难。
   * 更多信息请参见该命名空间的一般文档。
   *
   */
  template <int dim, int spacedim>
  void
  create_laplace_matrix(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const Quadrature<dim> &          q,
    SparseMatrix<double> &           matrix,
    const Function<spacedim> &       rhs,
    Vector<double> &                 rhs_vector,
    const Function<spacedim> *const  a           = nullptr,
    const AffineConstraints<double> &constraints = AffineConstraints<double>());

  /**
   * 调用create_laplace_matrix()函数，见上文，使用<tt>mapping=MappingQGeneric
   * @<dim@>(1)</tt>.  。
   *
   */
  template <int dim, int spacedim>
  void
  create_laplace_matrix(
    const DoFHandler<dim, spacedim> &dof,
    const Quadrature<dim> &          q,
    SparseMatrix<double> &           matrix,
    const Function<spacedim> &       rhs,
    Vector<double> &                 rhs_vector,
    const Function<spacedim> *const  a           = nullptr,
    const AffineConstraints<double> &constraints = AffineConstraints<double>());

  /**
   * 像上面的函数一样，但对hp-objects而言。
   *
   */
  template <int dim, int spacedim>
  void
  create_laplace_matrix(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &           dof,
    const hp::QCollection<dim> &                q,
    SparseMatrix<double> &                      matrix,
    const Function<spacedim> *const             a = nullptr,
    const AffineConstraints<double> &constraints = AffineConstraints<double>());

  /**
   * 像上面的函数一样，但对hp-objects而言。
   *
   */
  template <int dim, int spacedim>
  void
  create_laplace_matrix(
    const DoFHandler<dim, spacedim> &dof,
    const hp::QCollection<dim> &     q,
    SparseMatrix<double> &           matrix,
    const Function<spacedim> *const  a           = nullptr,
    const AffineConstraints<double> &constraints = AffineConstraints<double>());

  /**
   * 像上面的函数一样，但对hp-objects而言。
   *
   */
  template <int dim, int spacedim>
  void
  create_laplace_matrix(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim> &           dof,
    const hp::QCollection<dim> &                q,
    SparseMatrix<double> &                      matrix,
    const Function<spacedim> &                  rhs,
    Vector<double> &                            rhs_vector,
    const Function<spacedim> *const             a = nullptr,
    const AffineConstraints<double> &constraints = AffineConstraints<double>());

  /**
   * 像上面的函数一样，但对hp-objects而言。
   *
   */
  template <int dim, int spacedim>
  void
  create_laplace_matrix(
    const DoFHandler<dim, spacedim> &dof,
    const hp::QCollection<dim> &     q,
    SparseMatrix<double> &           matrix,
    const Function<spacedim> &       rhs,
    Vector<double> &                 rhs_vector,
    const Function<spacedim> *const  a           = nullptr,
    const AffineConstraints<double> &constraints = AffineConstraints<double>());

  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcComponentMismatch,
                   "You are providing either a right hand side function or a "
                   "coefficient with a number of vector components that is "
                   "inconsistent with the rest of the arguments. If you do "
                   "provide a coefficient or right hand side function, then "
                   "it either needs to have as many components as the finite "
                   "element in use, or only a single vector component. In "
                   "the latter case, the same value will be taken for "
                   "each vector component of the finite element.");
} // namespace MatrixCreator



/**
 * 提供一个对矩阵进行操作的函数集合。这些包括对线性方程组的边界条件的应用和其他。
 *
 *  <h3>Boundary conditions</h3>
 * apply_boundary_values()函数修改了一个线性系统，以纳入由Dirichlet型边界条件（或者，更具体地说："强
 * "边界条件）所产生的约束。为了真正做到这一点，当前命名空间中的这个名字的函数需要一个自由度指数的列表，以及这些自由度应该具有的值。要了解如何获得这样一个列表，请看
 * VectorTools::interpolate_boundary_values()
 * 函数的讨论，作为一个例子。
 * 有两种方法可以将固定自由度（如边界节点）纳入线性系统，如下文所述。这两种方法都是在局部对全局线性系统的贡献层面上操作，或者是全局系统本身。第三种方法，使用
 * AffineConstraints::copy_local_to_global(),
 * 执行相同的过程，作为将一个单元的局部贡献加入全局线性系统的一部分（"装配
 * "步骤），是目前教程程序中最主要的方法。
 * @dealiiVideoLecture{21.6,21.65}
 *
 * <h3>Global elimination</h3>
 * 在第一种方法中，我们首先在不尊重固定自由度的情况下装配全局线性系统，并在第二步中再次从线性系统中消除它们。纳入组装过程如下：当矩阵和向量被设置好后，要列出受Dirichlet边界条件约束的节点，并对矩阵和向量进行相应修改。这是通过删除矩阵中该自由度行的所有条目，将主对角线条目设置为一个合适的正值，并将右手边的元素设置为一个值，这样线性系统的解将在该节点有边界值。为了使剩余的线性方程组解耦并使系统再次对称（至少以前是这样的），要对这条线进行一个高斯消除步骤，将这条线（现在几乎是空的）加入到所有与给定自由度相耦合的其他线中，从而消除这个自由度与其他自由度之间的所有耦合。现在，除了主要的对角线条目之外，各自的列也只由零组成。另外，这个命名空间的函数有一个布尔参数，允许省略这最后一步，如果所产生的线性系统的对称性不需要的话。请注意，通常情况下，即使是CG也能应对具有这种特殊结构的非对称线性系统。
 * 寻找哪些行包含我们目前正在执行的高斯消除步骤的列中的条目是困难的，或者是非常简单的，取决于情况。如果稀疏模式是对称的（矩阵是否对称在这里无关紧要），那么我们可以通过查看当前行中哪些列是不空的来推断出在当前列中有一个非零条目的行。在这种情况下，我们只需要查看固定数量的行，不需要搜索所有的行。另一方面，如果稀疏模式是非对称的，那么我们需要使用一个迭代求解器，它在任何情况下都可以处理非对称矩阵，所以可能无论如何都不需要做高斯消除。事实上，这个函数就是这样工作的：它接受一个参数（
 * @p eliminate_columns)
 * ，指定稀疏模式是否是对称的；如果是，那么列就会被消除，右手边也会相应修改。如果不是，则只删除该行，完全不涉及该列，而且除了与当前行对应的值之外，所有右手边的值保持不变。
 * 如果你的矩阵的稀疏模式是非对称的，你必须把这个参数的值设置为
 * @p false
 * ，因为这样我们就不能在不搜索所有行的情况下消除这一列，这样就太昂贵了（如果
 * @p N  ]为行数， @p m
 * 为每行的非零元素数，那么消除一列就是<tt>O(N*log(m))</tt>操作，因为在每行搜索需要<tt>log(m)</tt>操作）。)
 * 如果你的稀疏性模式是对称的，但你的矩阵不是，那么你也可以指定
 * @p false
 * 。如果你的稀疏模式和矩阵都是对称的，你可能想指定 @p
 * true
 * （那么消除一行的复杂性是<tt>O(m*log(m))</tt>，因为我们只需要搜索
 * @p m 行的列的相应元素）。鉴于 @p m
 * 大致是恒定的，不管离散化如何，边界节点的数量在2d中是<tt>sqrt(N)</tt>，对称稀疏模式的算法是<tt>O(sqrt(N)*m*log(m))</tt>，而对于一般情况，它将是<tt>O(N*sqrt(N)*log(m))</tt>；后者太昂贵，无法执行。
 * 似乎我们必须明确，在做高斯消除步骤时不要覆盖其他边界节点的线。然而，由于我们在通过这样的节点时重置了右手边，所以改变其他尚未处理的边界节点的右手边的值并不是问题。改变已经处理过的节点的那些条目将是一个问题，但是由于已经处理过的节点的行上的本列矩阵条目是零，高斯步骤不会改变右手边。因此，我们不需要对其他边界节点进行特别的处理。
 * 为了使解题速度更快，我们用正确的边界值预设了解向量（至于为什么要这样做，请看下面关于局部消除的描述）。不清楚在使用迭代求解器时，删除边界自由度和其他自由度之间的耦合是否真的能迫使解向量中的相应条目具有正确的值，因为它们的搜索方向可能包含边界节点方向上的分量。出于这个原因，我们进行了一个非常简单的线平衡，不是将主对角线条目设置为一，而是设置为删除此线之前的值，如果主对角线条目由于某种原因为零，则设置为第一个非零主对角线条目。当然，我们必须适当改变右手边的内容。这不是一个很好的策略，但它至少应该给主对角线条目一个正确的维度顺序的值，这使得求解过程更稳定一些。一个精炼的算法会将该条目设置为其他对角线条目的平均值，但这似乎太昂贵了。
 * 在某些情况下，对同一矩阵进行多次求解，但对不同的右手边或边界值进行求解，可能是有趣的。一个典型的例子是解决一个随时间变化的问题，其中边界值或右手边发生变化，但矩阵本身没有变化。这时，人们可能会想只组装一次矩阵，然后在同一个矩阵对象上重复调用
 * MatrixTools::apply_boundary_values()
 * 函数，在每个时间步长中新形成一个右手边向量。然而，由于对右手向量的边界值的修改取决于原始矩阵，如果不把原始矩阵存储在某个地方，并在每个时间步长中用存储在其他地方的未修改的矩阵初始化系统矩阵，这是不可能的。
 * step-26
 * 通过存储组成系统矩阵的构件，对这一过程进行了变通，但一般原理是相同的。另外，我们可以使用constrained_linear_operator()函数。在它的文档中，你还可以找到一个正式的（数学）描述，即修改矩阵和右手边向量的边界值的过程。
 *
 *  <h3>Local elimination</h3>
 * 处理边界值的第二种方式是在将本地矩阵和向量的贡献转移到全局稀疏矩阵和向量之前，适当地修改它们。这就是local_apply_boundary_values()的作用。这样做的好处是，我们省去了对apply_boundary_values函数的调用（这个函数很昂贵，因为它必须在稀疏数据结构上工作）。另一方面，local_apply_boundary_values()函数被多次调用，即使我们只有非常少的固定边界节点，主要的缺点是，如果有悬挂的节点也需要处理，这个函数就不能像预期的那样工作。这个函数不起作用的原因是，它是要在分布到全局矩阵之前运行的，也就是在悬空节点分布之前；由于悬空节点可以被约束到一个边界节点，对悬空节点的处理会在对应于边界值的行和列上再次增加条目，而这些条目我们已经在局部消除步骤中腾出了。更糟糕的是，在3D中受约束的节点甚至可以位于边界上。因此，当务之急是边界节点的消除发生在悬空节点消除之后
 * @em
 * ，但这无法通过边界节点的局部消除来实现，除非根本就没有悬空节点约束。
 * 局部消除有一个额外的缺点：我们无法获得解向量，只能获得对矩阵和右侧的局部贡献。这方面的问题很微妙，但会导致非常难以发现的困难：当我们消除一个自由度时，我们删除这个未知数的行和列，并将对角线条目设置为某个正值。为了使问题或多或少具有良好的条件，我们将这个对角线条目设置为其先验值的绝对值（如果该值非零），或者设置为所有其他非零对角线元素的平均大小。然后，我们设置右边的值，使得到的解条目具有边界值所给的正确值。由于我们将这些贡献加到所有的局部贡献上，对角线条目和右手边的各自数值也相应地加起来，所以线性系统的解中的条目仍然有效。
 * 然而，如果这样选择的对角线条目不适合线性系统，就会出现问题。例如，考虑一个带有矩阵<tt>[[A
 * B][C^T
 * 0]]</tt>的混合拉普拉斯问题，我们只为解的第二部分指定边界值。在混合公式中，应力应变张量只出现在矩阵
 * @p B 或 @p C,
 * 中，所以其中一个可能明显大于或小于另一个。现在，如果我们消除边界值，就会删除一些行和列，但也会在右下角块的对角线上引入一些条目，这样我们就会得到系统<tt>[[A'
 * B'][C'^T X]]</tt>。矩阵 @p X 中的对角线项将与 @p A.
 * 中的对角线项具有相同的数量级。
 * 现在，如果我们用Schur补码公式解决这个系统，我们必须反转矩阵<tt>X-C'^TA'^{-1}B'</tt>。删除上面的行和列可以确保边界节点在Schur补数中确实也有空行和空列，除了
 * @p X. 中的条目。然而， @p X
 * 中的条目可能与<tt>C'^TA'^{-1}B'</tt>中的条目的数量级明显不同!
 * 如果是这样的话，我们可能会遇到迭代求解器的麻烦。例如，假设我们从解向量的零条目开始，而
 * @p X
 * 中的条目太小了几个数量级；在这种情况下，迭代求解器将在每一步计算残差向量并形成修正向量，但由于
 * @p X
 * 中的条目太小，边界节点的残差贡献确实很小，尽管边界节点的值仍然接近于零，不符合规定的边界值。由于残差如此之小，迭代求解器计算的修正值也非常小，最后求解器会显示收敛到一个小的总残差，而边界值仍然有明显的错误。
 * 在上述的全局消除过程中，我们通过给边界节点的正确值
 * "打底
 * "来避免这个问题。但是，在局部消除过程中，我们不能这样做。因此，如果你遇到类似上述的问题，你需要将
 * @p X
 * 中的对角线项增加到与舒尔补数的其他部分相匹配的大小，或者更简单，在启动求解器之前给解向量打底。
 * 总之，边界节点的局部消除只有在没有悬空节点的情况下才有效，即使如此也不一定能完全令人满意地工作。
 *
 *
 * @ingroup numerics
 *
 *
 */
namespace MatrixTools
{
  /**
   * 导入命名空间MatrixCreator，以便向后兼容旧版本的deal.II，其中这些命名空间是类，类MatrixTools是公开派生自类MatrixCreator。
   *
   */
  using namespace MatrixCreator;

  /**
   * 对系统矩阵和向量应用Dirichlet边界条件，如该命名空间的一般文档中所述。
   *
   */
  template <typename number>
  void
  apply_boundary_values(
    const std::map<types::global_dof_index, number> &boundary_values,
    SparseMatrix<number> &                           matrix,
    Vector<number> &                                 solution,
    Vector<number> &                                 right_hand_side,
    const bool                                       eliminate_columns = true);

  /**
   * 对系统矩阵和向量应用迪里切特边界条件，如本命名空间的一般文档中所述。这个函数适用于块状稀疏矩阵和块状向量。
   *
   */
  template <typename number>
  void
  apply_boundary_values(
    const std::map<types::global_dof_index, number> &boundary_values,
    BlockSparseMatrix<number> &                      matrix,
    BlockVector<number> &                            solution,
    BlockVector<number> &                            right_hand_side,
    const bool                                       eliminate_columns = true);

#ifdef DEAL_II_WITH_PETSC
  /**
   * 对系统矩阵和向量应用Dirichlet边界条件，如本命名空间的一般文档中所述。这个函数对用于包裹PETSc对象的类起作用。
   * <b>Important:</b>
   * 这个函数的效率不高：它需要交替地读和写到矩阵中，这种情况PETSc处理得不好。此外，我们只删除了与边界节点相对应的行，但删除相应的列的情况（即如果
   * @p eliminate_columns 是 @p true)
   * ，目前还没有实现，而且可能永远不会实现，因为如果不直接访问PETSc数据结构，成本太高。这就导致了最后一个参数的默认值所表示的动作实际上没有实现；该参数的默认值是
   * <code>true</code>
   * ，以保持与该命名空间中其他同名函数的一致性）。
   * 这个函数在  step-17  和  step-18  中使用。
   * @note
   * 如果矩阵是用MPI在多个处理器之间并行存储的，这个函数只触及本地存储的行，而简单地忽略所有其他行。换句话说，每个处理器负责自己的行，
   * @p boundary_values
   * 参数需要包含你想处理的矩阵中所有本地拥有的行。(但是它也可以包含不属于本地的自由度的条目；这些条目将被简单地忽略)。此外，在并行计算的背景下，如果你处理了某一行，而其他处理器还在等待对同一行的写入或添加，你就会陷入困境。换句话说，如果另一个处理器仍然想向某行的某个元素添加东西，而你调用这个函数将该行清零，那么你下次调用compress()时可能会将远程值添加到你刚刚创建的零值中。因此，你要在对矩阵进行最后一次修改后，在开始清空行之前调用compress()。
   *
   */
  void
  apply_boundary_values(
    const std::map<types::global_dof_index, PetscScalar> &boundary_values,
    PETScWrappers::MatrixBase &                           matrix,
    PETScWrappers::VectorBase &                           solution,
    PETScWrappers::VectorBase &                           right_hand_side,
    const bool eliminate_columns = true);

  /**
   * 和上面一样，但对于并行的BlockSparseMatrix。
   *
   */
  void
  apply_boundary_values(
    const std::map<types::global_dof_index, PetscScalar> &boundary_values,
    PETScWrappers::MPI::BlockSparseMatrix &               matrix,
    PETScWrappers::MPI::BlockVector &                     solution,
    PETScWrappers::MPI::BlockVector &                     right_hand_side,
    const bool eliminate_columns = true);

#endif

#ifdef DEAL_II_WITH_TRILINOS
  /**
   * 对系统矩阵和向量应用Dirichlet边界条件，如本命名空间的一般文档中所述。这个函数对用于包裹特里诺斯对象的类起作用。
   * <b>Important:</b>
   * 这个函数的效率不高：它需要交替地读和写到矩阵中，这种情况Trilinos处理得不好。此外，我们只删除了与边界节点相对应的行，但删除相应的列的情况（即如果
   * @p eliminate_columns 是 @p true)
   * ，目前没有实现，可能永远不会实现，因为如果不直接访问Trilinos数据结构，成本太高。这导致最后一个参数的默认值所表示的动作实际上没有实现；该参数的默认值为
   * <code>true</code>
   * ，以保持与本命名空间中其他同名函数的一致性）。
   * @note
   * 如果矩阵是用MPI在多个处理器之间并行存储的，这个函数只触及本地存储的行，而简单地忽略所有其他行。换句话说，每个处理器负责自己的行，
   * @p boundary_values
   * 参数需要包含你想要处理的矩阵的所有本地拥有的行。(但是它也可以包含不属于本地的自由度的条目；这些条目将被简单地忽略)。此外，在并行计算的背景下，如果你处理了某一行，而其他处理器还在等待对同一行的写入或添加，你就会陷入困境。换句话说，如果另一个处理器仍然想向某行的某个元素添加东西，而你调用这个函数将该行清零，那么你下次调用compress()时可能会将远程值添加到你刚刚创建的零值中。因此，你要在对矩阵进行最后一次修改后，在开始清空行之前调用compress()。
   *
   */
  void
  apply_boundary_values(
    const std::map<types::global_dof_index, TrilinosScalar> &boundary_values,
    TrilinosWrappers::SparseMatrix &                         matrix,
    TrilinosWrappers::MPI::Vector &                          solution,
    TrilinosWrappers::MPI::Vector &                          right_hand_side,
    const bool eliminate_columns = true);

  /**
   * 这个函数的作用与上面的函数相同，只是现在对块结构进行处理。
   *
   */
  void
  apply_boundary_values(
    const std::map<types::global_dof_index, TrilinosScalar> &boundary_values,
    TrilinosWrappers::BlockSparseMatrix &                    matrix,
    TrilinosWrappers::MPI::BlockVector &                     solution,
    TrilinosWrappers::MPI::BlockVector &                     right_hand_side,
    const bool eliminate_columns = true);
#endif

  /**
   * 这个函数不是在创建全局矩阵后将边界值应用于全局矩阵和向量，而是在装配过程中通过修改局部矩阵和向量的贡献来实现。如果你在所有的局部贡献上调用这个函数，得到的矩阵将有相同的条目，最后在全局系统上调用apply_boundary_values()就没有必要。
   * 由于这个函数不需要在稀疏矩阵的复杂数据结构上工作，所以它相对便宜。因此，如果你有很多固定的自由度（比如边界节点），或者对稀疏矩阵的访问很昂贵（比如对于块状稀疏矩阵，或者对于PETSc或Trilinos矩阵），它可能是一种胜利。然而，如果还有悬空节点需要考虑，它就不能如期工作了。更多的注意事项列在这个命名空间的一般文档中。
   * @dealiiVideoLecture{21.6,21.65}
   *
   */
  template <typename number>
  void
  local_apply_boundary_values(
    const std::map<types::global_dof_index, number> &boundary_values,
    const std::vector<types::global_dof_index> &     local_dof_indices,
    FullMatrix<number> &                             local_matrix,
    Vector<number> &                                 local_rhs,
    const bool                                       eliminate_columns);

  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcBlocksDontMatch,
                   "You are providing a matrix whose subdivision into "
                   "blocks in either row or column direction does not use "
                   "the same blocks sizes as the solution vector or "
                   "right hand side vectors, respectively.");
} // namespace MatrixTools



DEAL_II_NAMESPACE_CLOSE

#endif


