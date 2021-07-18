//include/deal.II-translator/fe/fe_tools_0.txt
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

#ifndef dealii_fe_tools_H
#define dealii_fe_tools_H



#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/component_mask.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <memory>
#include <string>
#include <vector>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename number>
class FullMatrix;
template <int dim>
class Quadrature;
template <int dim, int spacedim>
class FiniteElement;
template <int dim, int spacedim>
class DoFHandler;
template <int dim>
class FiniteElementData;
template <typename number>
class AffineConstraints;
#endif


 /*!@addtogroup feall */ 
 /*@{*/ 


/**
 * 这个命名空间提供了一个 @p FiniteElement  @p fe1
 * 的离散函数向另一个 @p FiniteElement  @p fe2的内插和外推。
 * 它还提供了在每个单元上进行插值的局部插值矩阵。此外，它还提供了差分矩阵
 * $id-I_h$ ，这是评估 $(id-I_h)z$ 所需要的，例如，对偶解 $z$
 * 。
 * 关于<tt>spacedim</tt>模板参数的更多信息，请查阅FiniteElement或Triangulation的文档。
 *
 *
 */
namespace FETools
{
  /**
   * 一个用于创建给定程度的有限元素的工厂对象的基类。每当人们想用一种透明的方式来创建有限元对象时，就会调用派生类。
   * 该类在 FETools::get_fe_by_name() 和 FETools::add_fe_name()
   * 函数中使用。
   *
   */
  template <int dim, int spacedim = dim>
  class FEFactoryBase : public Subscriptor
  {
  public:
    /**
     * 创建一个FiniteElement，并返回一个指向它的指针。
     *
     */
    virtual std::unique_ptr<FiniteElement<dim, spacedim>>
    get(const unsigned int degree) const = 0;

    /**
     * 从正交公式创建一个FiniteElement（目前只对FE_Q实现），并返回一个指针。
     *
     */

    virtual std::unique_ptr<FiniteElement<dim, spacedim>>
    get(const Quadrature<1> &quad) const = 0;

    /**
     * 虚拟析构器，除了让编译器高兴，什么都不做。
     *
     */
    virtual ~FEFactoryBase() override = default;
  };

  /**
   * 一个具体的类，用于创建特定程度的有限元的工厂对象。
   * 该类的get()函数生成一个作为模板参数的类型的有限元对象，并将度数（无论有限元类希望如何解释这个数字）作为get()的参数给出。
   *
   */
  template <class FE>
  class FEFactory : public FEFactoryBase<FE::dimension, FE::space_dimension>
  {
  public:
    /**
     * 创建一个FiniteElement并返回一个指针。
     *
     */
    virtual std::unique_ptr<FiniteElement<FE::dimension, FE::space_dimension>>
    get(const unsigned int degree) const override;

    /**
     * 从正交公式中创建一个FiniteElement（目前只对FE_Q实现）并返回一个指针。
     *
     */
    virtual std::unique_ptr<FiniteElement<FE::dimension, FE::space_dimension>>
    get(const Quadrature<1> &quad) const override;
  };

  /**
   * @warning  在大多数情况下，你可能想使用
   * compute_base_renumbering()。
   * 计算出一个单元格按分量重新编号所需的向量。
   * 此外，计算存储本地块向量中每个组件的起始索引的向量。
   * 第二个向量是这样组织的：每个基元都有一个向量，包含这个基元所服务的每个分量的起始索引。
   * 当第一个向量被检查到有正确的大小时，第二个向量被重新初始化以方便使用。
   *
   */
  template <int dim, int spacedim>
  void
  compute_component_wise(const FiniteElement<dim, spacedim> &    fe,
                         std::vector<unsigned int> &             renumbering,
                         std::vector<std::vector<unsigned int>> &start_indices);

  /**
   * 计算按区块对单元格的道夫重新编号所需的向量。
   * 此外，计算存储每个本地块向量的起始索引或大小的向量。
   * 如果 @p bool 参数为真， @p block_data
   * 将被填充为每个局部块的起始指数。如果它是假的，那么就返回块的大小。
   * 向量<tt>renumbering</tt>将以本地自由度的标准编号为索引，即第一个顶点，然后是第二个顶点，之后是顶点线、四边形和六边形。对于每个索引，条目表示这个自由度在一个编号方案中得到的索引，其中第一块的编号完全在第二块之前。
   *
   */
  template <int dim, int spacedim>
  void
  compute_block_renumbering(const FiniteElement<dim, spacedim> &  fe,
                            std::vector<types::global_dof_index> &renumbering,
                            std::vector<types::global_dof_index> &block_data,
                            bool return_start_indices = true);

  /**
   * @name 生成局部矩阵  
     * @{ 
   *
   */
  /**
   * 计算内插矩阵，将每个单元上的 @p fe1-function 内插到 @p
   * fe2-function 。插值矩阵的大小为<tt>(fe2.n_dofs_per_cell(),
   * fe1.n_dofs_per_cell())/tt>。    注意，如果有限元空间 @p fe1
   * 是有限元空间 @p fe2 的一个子集，那么 @p
   * interpolation_matrix 就是一个嵌入矩阵。
   *
   */
  template <int dim, typename number, int spacedim>
  void
  get_interpolation_matrix(const FiniteElement<dim, spacedim> &fe1,
                           const FiniteElement<dim, spacedim> &fe2,
                           FullMatrix<number> &interpolation_matrix);

  /**
   * 计算内插矩阵，将一个 @p fe1-function 内插到一个 @p
   * fe2-function, ，并将此内插到每个单元上的第二个 @p
   * fe1-function 。插值矩阵的大小为<tt>(fe1.n_dofs_per_cell(),
   * fe1.n_dofs_per_cell())/tt>。    注意，只有当 @p fe1
   * 的有限元空间不是 @p fe2,
   * 的有限元空间的子集时，这个函数才有意义，因为如果它是一个子集，那么
   * @p interpolation_matrix 就只有单位矩阵了。
   *
   */
  template <int dim, typename number, int spacedim>
  void
  get_back_interpolation_matrix(const FiniteElement<dim, spacedim> &fe1,
                                const FiniteElement<dim, spacedim> &fe2,
                                FullMatrix<number> &interpolation_matrix);

  /**
   * 计算身份矩阵减去回补矩阵。  在这个函数之后， @p
   * difference_matrix 的大小将是<tt>(fe1.n_dofs_per_cell(),
   * fe1.n_dofs_per_cell())/tt>。之前的参数内容将被覆盖。
   * 这个函数计算将 @p fe1 函数 $z$ 转换到 $z-I_hz$
   * 的矩阵，其中 $I_h$ 表示从 @p fe1 空间到 @p fe2
   * 空间的插值运算符。因此，这个矩阵对于评估误差呈现是有用的，其中
   * $z$ 表示对偶解。
   *
   */
  template <int dim, typename number, int spacedim>
  void
  get_interpolation_difference_matrix(const FiniteElement<dim, spacedim> &fe1,
                                      const FiniteElement<dim, spacedim> &fe2,
                                      FullMatrix<number> &difference_matrix);

  /**
   * 计算从fe1到fe2的局部 $L^2$ -投影矩阵。
   *
   */
  template <int dim, typename number, int spacedim>
  void
  get_projection_matrix(const FiniteElement<dim, spacedim> &fe1,
                        const FiniteElement<dim, spacedim> &fe2,
                        FullMatrix<number> &                matrix);

  /**
   * 这是一个相当专业的函数，在构建有限元对象时使用。它用于为一个元素建立形状函数的基础，给定一组多项式和插值点。该函数只针对正好有
   * @p dim
   * 个矢量分量的有限元实现。特别是，这适用于从FE_PolyTensor类派生的类。
   * 具体来说，这个函数的目的是这样的。FE_PolyTensor从其派生类中接收一个描述多项式空间的参数。这个空间可以用单项式或其他方式进行参数化，但一般来说，不是我们用于有限元的形式，在有限元中，我们通常要使用从某种节点函数（例如，特定点的内插）派生的基。
   * 具体来说，假设多项式空间使用的基是
   * $\{\tilde\varphi_j(\mathbf x)\}_{j=1}^N$
   * ，而有限元的节点函数是 $\{\Psi_i\}_{i=1}^N$
   * 。然后我们想为有限元空间计算一个基  $\{\varphi_j(\mathbf
   * x)\}_{j=1}^N$  ，以便  $\Psi_i[\varphi_j] = \delta_{ij}$
   * 。为了做到这一点，我们可以设置  $\varphi_j(\mathbf x) =
   * \sum_{k=1}^N c_{jk} \tilde\varphi_k(\mathbf x)$
   * ，我们需要确定扩展系数  $c_{jk}$  。我们通过将 $\Psi_i$
   * 应用于方程的两边来做到这一点，得到
   * @f{align*}{
   * \Psi_i [\varphi_j] = \sum_{k=1}^N c_{jk} \Psi_i[\tilde\varphi_k],
   * @f}
   * 我们知道左手边等于  $\delta_{ij}$  。
   * 如果你认为这是一个 $N\times N$
   * 的方程系统，左边和右边的矩阵的元素，那么这可以写为
   * @f{align*}{
   * I = C X^T
   * @f}
   * 其中 $C$ 是系数 $c_{jk}$ 和 $X_{ik} = \Psi_i[\tilde\varphi_k]$
   * 的矩阵。因此，为了计算膨胀系数 $C=X^{-T}$
   * ，我们需要将节点函数应用于多项式空间的 "原始
   * "基础的所有函数。    在有限元收到这个矩阵 $X$
   * 回来之前，它描述其形状函数（例如，在
   * FiniteElement::shape_value()) 中的形式 $\tilde\varphi_j$
   * 。在它调用这个函数后，它有了膨胀系数，可以描述它的形状函数为
   * $\varphi_j$  。    因此，这个函数计算这个矩阵  $X$
   * ，具体情形如下。
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
   * - 有限元素正好有 @p dim 个矢量分量。
   *
   *
   *
   *
   *
   * - 该函数 $f_j$ 是由任何元素通过 FiniteElement::convert_generalized_support_point_values_to_dof_values() 函数实现的。      @param  fe 要对其进行上述操作的有限元。    @return  如上所述的矩阵 $X$ 。
   *
   */
  template <int dim, int spacedim>
  FullMatrix<double>
  compute_node_matrix(const FiniteElement<dim, spacedim> &fe);

  /**
   * 对于所有可能的（各向同性和各向异性）细化情况，计算从一个粗单元到子单元的嵌入矩阵。所得矩阵的每一列都包含了粗网格基函数在细网格基上的表示；矩阵被分割，使得每个子单元都有一个矩阵。
   * 这个函数在足够多的正交点中计算粗网格函数，并使用最小二乘法近似拟合细网格函数。因此，这个函数的使用只限于有限元空间实际上是嵌套的情况。
   * 注意， <code>matrices[refinement_case-1][child]</code>
   * 包括了RefinementCase  <code>refinement_case</code> 的子
   * <code>child</code> 的嵌入（或延长）矩阵。这里，我们使用
   * <code>refinement_case-1</code> instead of <code>refinement_case</code>
   * ，因为对于 RefinementCase::no_refinement(=0)
   * ，没有可用的延长矩阵。
   * 通常，这个函数被FiniteElement类的各种实现所调用，以填充各自的
   * FiniteElement::prolongation 矩阵。      @param  fe
   * 我们为其计算嵌入矩阵的有限元类。      @param  矩阵
   * 对FullMatrix对象的 RefinementCase<dim>::isotropic_refinement
   * 向量的引用。每个向量对应一个RefinementCase  @p
   * refinement_case ，向量大小为
   * GeometryInfo<dim>::n_children(refinement_case).
   * 这是在FiniteElement中使用的格式，我们主要想在这里使用这个函数。
   * @param  isotropic_only 设置为 <code>true</code>
   * 如果你只想计算各向同性的细化矩阵。      @param
   * threshold是计算嵌入的最小二乘法中允许的差距。
   *
   */
  template <int dim, typename number, int spacedim>
  void
  compute_embedding_matrices(
    const FiniteElement<dim, spacedim> &          fe,
    std::vector<std::vector<FullMatrix<number>>> &matrices,
    const bool                                    isotropic_only = false,
    const double                                  threshold      = 1.e-12);

  /**
   * 计算约束矩阵所需的面的嵌入矩阵。      @param  fe
   * 要计算这些矩阵的有限元。      @param  matrices
   * 一个<i>GeometryInfo<dim>::subfaces_per_face =
   * 2<sup>dim-1</sup></i>
   * FullMatrix对象的数组，保存每个子面的嵌入矩阵。
   * @param  face_coarse 计算该面的粗面的面的编号。      @param
   * face_fine 脸部精细化的脸部数量，这是为其计算的。
   * @param
   * threshold是计算嵌入的最小二乘法算法中允许的差距。
   * @warning
   * 这个函数将被用于计算约束矩阵。它还没有得到充分的测试。
   *
   */
  template <int dim, typename number, int spacedim>
  void
  compute_face_embedding_matrices(
    const FiniteElement<dim, spacedim> &fe,
    FullMatrix<number> (&matrices)[GeometryInfo<dim>::max_children_per_face],
    const unsigned int face_coarse,
    const unsigned int face_fine,
    const double       threshold = 1.e-12);

  /**
   * 对于所有可能的（各向同性和各向异性）细化情况，计算从子单元到粗略单元的<i>L<sup>2</sup></i>投影矩阵。
   * 注意， <code>matrices[refinement_case-1][child]</code>
   * 包括RefinementCase  <code>refinement_case</code> 的子集
   * <code>child</code> 的投影（或限制）矩阵。这里，我们使用
   * <code>refinement_case-1</code> instead of <code>refinement_case</code>
   * ，因为对于 RefinementCase::no_refinement(=0)
   * ，没有可用的投影矩阵。
   * 通常，这个函数被FiniteElement类的各种实现所调用，以填充各自的
   * FiniteElement::restriction 矩阵。      @arg[in]  fe
   * 我们为其计算投影矩阵的有限元类。      @arg[out]  矩阵
   * 对一组 <tt>RefinementCase<dim>::isotropic_refinement</tt>
   * FullMatrix对象向量的引用。每个向量对应一个RefinementCase
   * @p refinement_case ，向量大小为
   * <tt>GeometryInfo<dim>::n_children(refinement_case)</tt>.
   * 这是在FiniteElement中使用的格式，我们主要想在这里使用这个函数。
   * @arg[in]  isotropic_only 如果设置为  <code>true</code>
   * ，那么这个函数只计算各向同性的细化情况的数据。输出向量的其他元素将不被触及（但仍然存在）。
   *
   */
  template <int dim, typename number, int spacedim>
  void
  compute_projection_matrices(
    const FiniteElement<dim, spacedim> &          fe,
    std::vector<std::vector<FullMatrix<number>>> &matrices,
    const bool                                    isotropic_only = false);

  /**
   * 将定义在正交点的标量数据投射到单个单元上的有限元空间。
   * 这个函数的作用如下：假设有标量数据<tt>u<sub>q</sub>, 0
   * <= q <
   * Q:=quadrature.size()</tt>定义在一个单元的正交点，这些点由给定的<tt>rhs_quadrature</tt>对象定义。然后，我们可能想问在给定的FE对象所定义的有限元空间中的有限元函数（在单个单元上）<tt>v<sub>h</sub></tt>，它是<tt>u</tt>在下列意义上的投影。
   * 通常，投影<tt>v<sub>h</sub></tt>是满足<tt>(v<sub>h</sub>,w)=(u,w)</tt>对所有离散测试函数<tt>w</tt>的函数。在目前的情况下，我们无法评估右侧，因为<tt>u</tt>只在<tt>rhs_quadrature</tt>给出的正交点中定义，所以我们用正交近似值代替它。同样，左手边也是用<tt>lhs_quadrature</tt>对象来近似；如果这个正交对象选择得当，那么左手边的积分就可以准确完成，不需要任何近似。如果右手边的正交对象的正交点太少，则有必要使用不同的正交对象
   *
   * --例如，如果数据<tt>q</tt>只定义在单元格中心，那么相应的一点正交公式显然不足以用定式来近似左手边的标量积。
   * 经过这些正交近似，我们最终得到一个<tt>V<sub>h</sub></tt>的<tt>v<sub>h</sub></tt>的节点表示，满足以下线性方程组。<tt>M
   * V<sub>h</sub> = Q
   * U</tt>，其中<tt>M<sub>ij</sub>=(phi_i,phi_j)</tt>是由<tt>lhs_quadrature</tt>近似的质量矩阵。和<tt>Q</tt>是矩阵<tt>Q<sub>iq</sub>=phi<sub>i</sub>(x<sub>q</sub>)
   * w<sub>q</sub></tt>
   * 其中<tt>w<sub>q</sub></tt>是正交权重。<tt>U</tt>是正交点数据的矢量
   * <tt>u<sub>q</sub></tt>。
   * 为了得到<tt>V<sub>h</sub></tt>的投影的节点表示，我们计算<tt>V<sub>h</sub>
   * = X U, X=M<sup>-1</sup>
   * Q</tt>。这个函数的目的是计算矩阵<tt>X</tt>并通过这个函数的最后一个参数返回。
   * 请注意，这个函数目前只支持标量数据。质量矩阵的扩展当然是微不足道的，但如果在所有正交点都包含矢量值的数据，就必须定义矢量<tt>U</tt>中的数据顺序。
   * 在 step-18 示例程序的介绍中描述了这个函数的用途。
   * 这个函数的反面，将有限元函数插值到正交点上，基本上就是
   * <tt>FEValues::get_function_values</tt>
   * 函数的作用；为了使事情变得简单一点，
   * <tt>FETools::compute_interpolation_to_quadrature_points_matrix</tt>
   * 提供了这个函数的矩阵形式。
   * 请注意，这个函数是在单个单元上工作的，而不是在整个三角形上。因此，实际上，如果你使用连续或不连续版本的有限元，这并不重要。
   * 值得注意的是，这个函数有几个令人困惑的情况。
   * 第一个是，真正有意义的是投射到每个单元的自由度最多与正交点一样多的有限元上；N个正交点数据投射到一个有M>N个未知数的空间是定义明确的，但往往会产生有趣和不直观的结果。其次，人们认为如果正交点数据在有限元的支持点中定义，即<tt>ths_quadrature</tt>的正交点等于<tt>fe.get_unit_support_points()</tt>，那么投影应该是同一的，即有限元的每个自由度等于相应形状函数的支持点中给定数据的值。然而，一般情况下不是这样的：虽然在这种情况下矩阵<tt>Q</tt>是身份矩阵，但质量矩阵<tt>M</tt>不等于身份矩阵，除了正交公式<tt>lhs_quadrature</tt>在有限元的支持点上也有其正交点的特殊情况。
   * 最后，这个函数只定义了一个单元的投影，而人们经常希望将其应用于三角结构中的所有单元。然而，如果它被应用于一个又一个的单元，如果自由度存在于单元之间的界面上，后面的单元的结果可能会覆盖之前的单元已经计算出来的节点值。因此，该函数对不连续的元素最有用。
   *
   */
  template <int dim, int spacedim>
  void
  compute_projection_from_quadrature_points_matrix(
    const FiniteElement<dim, spacedim> &fe,
    const Quadrature<dim> &             lhs_quadrature,
    const Quadrature<dim> &             rhs_quadrature,
    FullMatrix<double> &                X);

  /**
   * 给定一个（标量）局部有限元函数，计算矩阵，将节点值的矢量映射到该函数在正交点的值的矢量上，如第二个参数所给的。在某种意义上，这个函数与
   * FETools::compute_projection_from_quadrature_points_matrix
   * 函数的作用相反。
   *
   */
  template <int dim, int spacedim>
  void
  compute_interpolation_to_quadrature_points_matrix(
    const FiniteElement<dim, spacedim> &fe,
    const Quadrature<dim> &             quadrature,
    FullMatrix<double> &                I_q);

  /**
   * 计算存储在正交点 @p vector_of_tensors_at_qp
   * 的张量（一阶张量）数据对单元支持点的数据 @p
   * 张量的矢量_节点的投影。  @p vector_of_tensors_at_qp
   * 中的数据是按照正交点编号顺序排列的。  @p
   * vector_of_tensors_at_qp 的大小必须对应于 @p projection_matrix.
   * 的列数， @p vector_of_tensors_at_nodes 的大小必须对应于 @p
   * 向量的行数。 投影矩阵 @p projection_matrix
   * 描述了来自正交点的标量数据的投影，可以通过
   * FETools::compute_projection_from_quadrature_points_matrix 函数得到。
   *
   */
  template <int dim>
  void
  compute_projection_from_quadrature_points(
    const FullMatrix<double> &         projection_matrix,
    const std::vector<Tensor<1, dim>> &vector_of_tensors_at_qp,
    std::vector<Tensor<1, dim>> &      vector_of_tensors_at_nodes);



  /**
   * 与上一个函数相同，但对一个 @p SymmetricTensor 。
   *
   */
  template <int dim>
  void
  compute_projection_from_quadrature_points(
    const FullMatrix<double> &                  projection_matrix,
    const std::vector<SymmetricTensor<2, dim>> &vector_of_tensors_at_qp,
    std::vector<SymmetricTensor<2, dim>> &      vector_of_tensors_at_nodes);



  /**
   * 这个方法实现了
   * FETools::compute_projection_from_quadrature_points_matrix
   * 方法，用于网格的面。
   * 它返回的矩阵，X，是针对面的，其大小是fe.n_dofs_per_cell()乘以rhs_quadrature.size()。
   * 该类的尺寸，dim必须大于1，因为需要Quadrature<dim-1>对象。更多信息请参见正交类的文档。
   *
   */
  template <int dim, int spacedim>
  void
  compute_projection_from_face_quadrature_points_matrix(
    const FiniteElement<dim, spacedim> &fe,
    const Quadrature<dim - 1> &         lhs_quadrature,
    const Quadrature<dim - 1> &         rhs_quadrature,
    const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
    const unsigned int                                              face,
    FullMatrix<double> &                                            X);



  /**
   * 围绕
   * FiniteElement::convert_generalized_support_point_values_to_dof_values()
   * 的封装器，可用于任意的数字类型。      @param[in]
   * finite_element 用于计算dof值的FiniteElement。    @param[in]
   * support_point_values 一个大小为 @p dofs_per_cell
   * 的数组（相当于get_generalized_support_points()函数将返回的点的数量），其中每个元素是一个向量，其条目数量与该元素的向量分量相同。这个数组应该包含有限元的广义支持点上的函数值。
   * @param[out]  dof_values 一个大小为 @p dofs_per_cell
   * 的数组，包含应用于给定函数的元素的节点函数值。
   *
   */
  template <int dim, int spacedim, typename number>
  void
  convert_generalized_support_point_values_to_dof_values(
    const FiniteElement<dim, spacedim> &finite_element,
    const std::vector<Vector<number>> & support_point_values,
    std::vector<number> &               dof_values);



  //@}
  /**
   * @name  应该在DoFTools中出现的函数
   *
   */
  //@{
  /**
   * 计算一个 @p dof1-function  @p u1 到 @p Dof2-函数 @p u2.  @p dof1
   * 和 @p dof2 的内插，需要基于相同的三角形的DoFHandler。
   * 如果元素 @p fe1 和 @p fe2
   * 都是连续的或者都是不连续的，那么这个插值就是通常的点插值。
   * 如果 @p fe1 是一个连续的和 @p fe2
   * 是一个不连续的有限元也是如此。对于 @p fe1
   * 是不连续的， @p fe2
   * 是连续有限元的情况，在不连续处没有定义点插值。
   * 因此，平均值是在不连续点上的DoF值处取的。
   * 请注意，对于有悬挂节点的网格上的连续元素（即局部细化的网格），这个函数并不能给出预期的输出。
   * 事实上，产生的输出矢量在悬空节点处不一定尊重连续性要求：例如，如果你是将一个Q2场插值到Q1场，那么在悬空节点处，输出场将具有输入场的函数值，然而这通常不是两个相邻节点的平均值。因此，它不是整个三角形上Q1函数空间的一部分，尽管它当然是每个单元上的Q1。
   * 对于这种情况（带有悬挂节点的网格上的连续元素），请使用
   * @p interpolate()
   * 函数，并以额外的AffineConstraints对象作为参数，见下文，或者通过调用悬挂节点约束对象的
   * @p distribute 函数使场符合。
   *
   */
  template <int dim, int spacedim, class InVector, class OutVector>
  void
  interpolate(const DoFHandler<dim, spacedim> &dof1,
              const InVector &                 u1,
              const DoFHandler<dim, spacedim> &dof2,
              OutVector &                      u2);

  /**
   * 计算 @p dof1-function  @p u1 到 @p Dof2-函数的插值 @p u2.  @p
   * dof1 和 @p dof2 需要是基于同一三角的DoFHandlers。  @p
   * constraints 是对应于 @p dof2.
   * 的悬空节点约束对象，当插值到具有悬空节点的网格（局部细化网格）上的连续元素时，该对象特别重要。
   * 如果元素 @p fe1 和 @p fe2
   * 都是连续的或者都是不连续的，那么这个内插就是通常的点内插。
   * 如果 @p fe1 是一个连续的和 @p fe2
   * 是一个不连续的有限元也是如此。对于 @p fe1
   * 是不连续的和 @p fe2
   * 是连续有限元的情况，在不连续处没有定义点插值。
   * 因此，均值取自不连续处的DoF值。
   *
   */
  template <int dim, int spacedim, class InVector, class OutVector>
  void
  interpolate(
    const DoFHandler<dim, spacedim> &                        dof1,
    const InVector &                                         u1,
    const DoFHandler<dim, spacedim> &                        dof2,
    const AffineConstraints<typename OutVector::value_type> &constraints,
    OutVector &                                              u2);

  /**
   * 计算 @p fe1-function  @p u1 到 @p
   * fe2-函数的内插，并将其内插到第二个 @p fe1-function 名为
   * @p 的u1_interpolated。
   * 注意，这个函数对悬挂节点的连续元素不起作用。对于这种情况，请使用下面的
   * @p back_interpolate 函数，它需要一个额外的 @p
   * AffineConstraints 对象。    此外，请注意，对于 @p fe1
   * 对应的有限元空间是 @p fe2,
   * 对应的有限元空间的一个子集的特殊情况，这个函数只是一个身份映射。
   *
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  back_interpolate(const DoFHandler<dim, spacedim> &   dof1,
                   const InVector &                    u1,
                   const FiniteElement<dim, spacedim> &fe2,
                   OutVector &                         u1_interpolated);

  /**
   * 计算 @p dof1-function  @p u1 到 @p
   * dof2-函数的内插，并将其内插到第二个 @p dof1-function
   * ，命名为 @p u1_interpolated.   @p constraints1 和 @p constraints2
   * 分别是对应于 @p dof1 和 @p dof2, 的悬挂节点约束。
   * 当涉及到带有悬挂节点的网格（局部细化网格）上的连续元素时，这些对象特别重要。
   * 此外，请注意，对于 @p dof1 对应的有限元空间是 @p dof2,
   * 对应的有限元空间的一个子集的特定情况，这个函数只是一个身份映射。
   *
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  back_interpolate(
    const DoFHandler<dim, spacedim> &                        dof1,
    const AffineConstraints<typename OutVector::value_type> &constraints1,
    const InVector &                                         u1,
    const DoFHandler<dim, spacedim> &                        dof2,
    const AffineConstraints<typename OutVector::value_type> &constraints2,
    OutVector &                                              u1_interpolated);

  /**
   * 对于给定的 @p dof1-function  $z_1$ 计算 $(Id-I_h)z_1$ ，其中
   * $I_h$ 是从 @p fe1 到 @p fe2. 的插值，结果 $(Id-I_h)z_1$
   * 被写入 @p z1_difference.
   * 注意，这个函数对于悬挂节点的连续元素不起作用。对于这种情况，请使用下面的
   * @p interpolation_difference 函数，它需要一个额外的 @p
   * AffineConstraints 对象。
   *
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  interpolation_difference(const DoFHandler<dim, spacedim> &   dof1,
                           const InVector &                    z1,
                           const FiniteElement<dim, spacedim> &fe2,
                           OutVector &                         z1_difference);

  /**
   * 对于给定的 @p dof1-function  $z_1$ 计算 $(Id-I_h)z_1$ ，其中
   * $I_h$ 是从 @p fe1 到 @p fe2. 的插值，结果 $(Id-I_h)z_1$
   * 被写入 @p z1_difference.   @p constraints1 和 @p constraints2
   * 是悬挂节点约束，分别对应 @p dof1  和  @p dof2,
   * 。当涉及到带有悬挂节点的网格（局部细化网格）上的连续元素时，这些对象特别重要。
   * 对于并行计算，提供 @p z1 和 @p
   * z1_difference，不含鬼魂元素。
   *
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  interpolation_difference(
    const DoFHandler<dim, spacedim> &                        dof1,
    const AffineConstraints<typename OutVector::value_type> &constraints1,
    const InVector &                                         z1,
    const DoFHandler<dim, spacedim> &                        dof2,
    const AffineConstraints<typename OutVector::value_type> &constraints2,
    OutVector &                                              z1_difference);



  /**
   * $L^2$  不连续元素的投影。操作方向与插值相同。
   * 如果有限元空间是不连续的，全局投影可以通过局部矩阵来计算。对于连续元素，这是不可能的，因为全局质量矩阵必须被倒置。
   *
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  project_dg(const DoFHandler<dim, spacedim> &dof1,
             const InVector &                 u1,
             const DoFHandler<dim, spacedim> &dof2,
             OutVector &                      u2);

  /**
   * 计算 @p dof1 函数 @p z1 到 @p Dof2函数 @p z2. 的修补外推 @p
   * dof1 和 @p dof2
   * 需要是基于同一三角的DoFHandler对象。这个函数用于，例如，将逐片线性解外推到逐片二次解。
   * 这个函数的名字是历史性的，可能不是特别好的选择。该函数一个接一个地执行以下操作。
   *
   *
   *
   *
   *
   *
   * - 它直接从 @p dof1 的每个单元内插到`dof2`的相应单元，使用这些单元上使用的有限元空间的内插矩阵，并由相关的有限元对象提供。这一步是使用 FETools::interpolate() 函数完成的。
   *
   *
   *
   *
   *
   *
   * - 然后它在`dof2`的所有非活动单元上执行循环。  如果这样的非活动单元至少有一个活动的子单元，那么我们称这个单元的子单元为 "补丁"。然后我们使用与`dof2`相关的有限元空间，从这个补丁的子单元插补到补丁，并立即插补回子单元。从本质上讲，这个信息抛开了解向量中所有生活在比补丁单元更小范围内的信息。
   *
   *
   *
   *
   *
   * - 由于我们从最粗的层次到最细的层次遍历非活动单元，如果网格被自适应地细化，我们可能会发现对应于先前处理过的补丁的子单元的补丁（如果网格被全局细化，这种情况不会发生，因为那里的补丁的子单元都是活动的）。我们也对这些补丁进行上述操作，但很容易看出，在作为先前处理过的补丁的子单元的补丁上，该操作现在是身份操作（因为它从当前补丁的子单元插值一个先前从更粗的补丁插值到这些子单元的函数）。因此，这不会再改变解向量。    这个函数的名字源于这样一个事实：它可以用来在一个更粗的网格上构造一个更高多项式程度的函数的表示。例如，如果你想象你从一个全局细化的网格上的 $Q_1$ 函数开始，并且 @p dof2 与 $Q_2$ 元素相关联，那么这个函数计算出在一个网格尺寸为 $2h$ 的一次粗化的网格上将原始片状线性函数插值为二次函数的等效算子 $I_{2h}^{(2)}$ （但在原始网格上表示这个函数尺寸为 $h$  ）。  如果精确解足够平滑，那么 $u^\ast=I_{2h}^{(2)}u_h$ 通常比 $u_h$ 更能接近PDE的精确解 $u$ 。换句话说，这个函数提供了一个后处理步骤，以类似于人们经常通过外推一连串的解决方案获得的方式改进解决方案，解释了该函数名称的由来。
   * @note
   * 如果使用上述算法，产生的场不满足给定有限元的连续性要求。当你在有悬挂节点的网格上使用连续元素时，请使用带有额外AffineConstraints参数的
   * @p  外推函数，见下文。
   * @note
   * 由于这个函数在单元的补丁上操作，它要求底层网格对每个粗略的网格单元至少进行一次细化。
   * 如果不是这样，就会产生一个异常。
   *
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  extrapolate(const DoFHandler<dim, spacedim> &dof1,
              const InVector &                 z1,
              const DoFHandler<dim, spacedim> &dof2,
              OutVector &                      z2);

  /**
   * 计算 @p dof1 函数 @p z1 对 @p Dof2函数 @p z2. 的修补外推 @p
   * dof1 和 @p dof2 需要是基于同一三角的DoFHandler对象。   @p
   * constraints 是对应于 @p dof2.
   * 的悬空节点约束对象，当插值到有悬空节点的网格上的连续元素时，这个对象是必要的（局部细化网格）。
   * 否则，该函数的作用与上述另一个 @p extrapolate
   * 函数相同（对于该函数，文档中提供了大量的操作描述）。
   *
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  extrapolate(
    const DoFHandler<dim, spacedim> &                        dof1,
    const InVector &                                         z1,
    const DoFHandler<dim, spacedim> &                        dof2,
    const AffineConstraints<typename OutVector::value_type> &constraints,
    OutVector &                                              z2);

  //@}
  /**
   * 连续有限元中自由度的编号是分层次的，也就是说，我们首先按照三角形定义的顶点顺序对顶点自由度进行编号，然后按照顺序和尊重线的方向对线自由度进行编号，然后对四边形的自由度进行编号，等等。然而，我们也可以用词法对它们进行编号，即指数先按X方向，然后按Y方向，最后按Z方向。例如，FE_DGQ()类的不连续元素就是以这种方式编号的。
   * 这个函数返回一个向量，其中包含了分级编号中每个自由度对连续有限元的给定程度的词汇索引的信息。
   *
   */
  template <int dim>
  std::vector<unsigned int>
  hierarchic_to_lexicographic_numbering(unsigned int degree);

  /**
   * 这是与上述函数相反的函数，生成连续有限元给定多项式程度的词典式编号到层次式编号的映射。所有关于上述函数的评论在这里也是有效的。
   *
   */
  template <int dim>
  std::vector<unsigned int>
  lexicographic_to_hierarchic_numbering(unsigned int degree);

  /**
   * 一个命名空间，包含了在实现FiniteElement时帮助设置内部数据结构的函数，FiniteElement是由更简单的（"基"）元素构建的，例如FESystem。这些函数计算出来的东西通常作为被构造的派生有限元对象的FiniteElement基类的构造器参数。    一般来说，有两种方法可以构建更复杂的元素，这反映在这个命名空间的几个函数的参数被称为  <code>do_tensor_product</code>  ：  <ol>   <li>  张量积构造（  <code>do_tensor_product=true</code>  ）。  张量积构造，在最简单的情况下，从标量元素中建立一个矢量值元素（更多信息见 @ref vector_valued  "本文档模块 "和 @ref GlossComponent  "本词汇条"）。  举个例子，考虑创建一个有两个矢量分量的矢量值元素，其中第一个应该有线性形状函数，第二个有二次形状函数。在1d中，基础元素的形状函数（在参考单元上）则是
   * @f{align*}{
   * Q_1 &= \{ 1-x, x \},
   * \\  Q_2 &= \{ 2(\frac 12
   *
   * - x)(1-x), 2(x
   *
   * - \frac 12)x, 4x(1-x) \},
   * @f}
   * 其中形状函数以通常的方式排序（首先在第一个顶点上，然后在第二个顶点上，然后在单元的内部）。张量乘积结构将创建一个具有以下形状函数的元素。
   * @f{align*}{
   * Q_1 \times Q_2 &=
   * \left\{
   *   \begin{pmatrix} 1-x \\ 0 \end{pmatrix},
   *   \begin{pmatrix} 0 \\ 2(\frac 12
   *
   * - x)(1-x)  \end{pmatrix},
   *   \begin{pmatrix} x \\ 0 \end{pmatrix},
   *   \begin{pmatrix} 0 \\ 2(x
   *
   * - \frac 12)x \end{pmatrix},
   *   \begin{pmatrix} 0 \\ 4x(1-x) \end{pmatrix}
   * \right\}.
   * @f}
   * 这里的列表又是按标准顺序排列的。
   * 当然，如果基数元素本身已经是矢量值，这个过程也是有效的：在这种情况下，组成的元素只是具有与基数元素加起来一样多的矢量成分。
   * <li>  组合形状函数（ <code>do_tensor_product=false</code>
   * ）。与之前的策略相反，组合形状函数只是将<i>all</i>的形状函数合在一起。在上面的例子中，这将产生以下元素。
   * @f{align*}{
   * Q_1 + Q_2 &= \{ 1-x, 2(\frac 12
   *
   * - x)(1-x),
   *                 x, 2(x
   *
   * - \frac 12)x, 4x(1-x) \}.
   * @f}
   * 换句话说，如果基础元素是标量的，产生的元素也将是标量的。一般来说，基本元素都必须有相同数量的矢量成分。    当然，上面构建的元素不再有一组线性独立的形状函数。因此，通过以同样的方式处理组成元素的所有形状函数而创建的任何矩阵将是奇异的。因此，在实践中，这种策略通常用于明确确保某些形状函数被区别对待的情况（例如，通过与权重函数相乘），或者用于所组合的形状函数不是线性独立的情况。      </ol>
   *
   */
  namespace Compositing
  {
    /**
     * 取有限元和乘数的向量，乘出组成的元素在每个顶点、线条等方面有多少自由度。
     * 如果 @p do_tensor_product
     * 为真，在FiniteElementData对象中返回的构件数是每个有限元中的构件数乘以相应的倍数的乘积之和。
     * 否则，组件数取自第一个具有非零倍数的有限元素，所有其他具有非零倍数的元素需要具有相同数量的矢量组件。
     * 关于 FETools::Compositing
     * 参数的更多信息，请参见命名空间 @p do_tensor_product
     * 的文档。
     *
     */
    template <int dim, int spacedim>
    FiniteElementData<dim>
    multiply_dof_numbers(
      const std::vector<const FiniteElement<dim, spacedim> *> &fes,
      const std::vector<unsigned int> &                        multiplicities,
      const bool do_tensor_product = true);

    /**
     * 对于任意数量的
     * <code>std::pair<std::unique_ptr<FiniteElement<dim,  spacedim>>,
     * unsigned int></code>和 <code>do_tensor_product = true</code>
     * 类型的参数，与上述相同。
     *
     */
    template <int dim, int spacedim>
    FiniteElementData<dim>
    multiply_dof_numbers(
      const std::initializer_list<
        std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>>
        &fe_systems);

    /**
     * 与上述相同，但对特定数量的子元素。
     *
     */
    template <int dim, int spacedim>
    FiniteElementData<dim>
    multiply_dof_numbers(const FiniteElement<dim, spacedim> *fe1,
                         const unsigned int                  N1,
                         const FiniteElement<dim, spacedim> *fe2 = nullptr,
                         const unsigned int                  N2  = 0,
                         const FiniteElement<dim, spacedim> *fe3 = nullptr,
                         const unsigned int                  N3  = 0,
                         const FiniteElement<dim, spacedim> *fe4 = nullptr,
                         const unsigned int                  N4  = 0,
                         const FiniteElement<dim, spacedim> *fe5 = nullptr,
                         const unsigned int                  N5  = 0);

    /**
     * 计算 "限制是加法的
     * "标志（见FiniteElement类的文档），用于第二个参数中给出的乘数的有限元列表。
     * 限制是相加的
     * "标志是单个形状函数的属性，不取决于组成元素是否使用
     * FETools::Composition
     * 命名空间的文档中概述的张量乘积或组合策略。因此，这个函数没有
     * @p do_tensor_product 参数。
     *
     */
    template <int dim, int spacedim>
    std::vector<bool>
    compute_restriction_is_additive_flags(
      const std::vector<const FiniteElement<dim, spacedim> *> &fes,
      const std::vector<unsigned int> &                        multiplicities);

    /**
     * 与上述相同，用于任意数量的
     * <code>std::pair<std::unique_ptr<FiniteElement<dim,  spacedim>>,
     * unsigned int></code>类型的参数。
     *
     */
    template <int dim, int spacedim>
    std::vector<bool>
    compute_restriction_is_additive_flags(
      const std::initializer_list<
        std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>>
        &fe_systems);

    /**
     * 取一个 @p FiniteElement
     * 对象并返回一个布尔向量，描述由 @p N1,  @p N2,
     * ...子元素 @p fe1,  @p fe2,
     * 的副本组成的混合元素的每个形状函数的 @p
     * restriction_is_additive_flags （见FiniteElement类的文档）。
     * "限制是加法的
     * "标志是单个形状函数的属性，不取决于组成元素是否使用
     * FETools::Composition
     * 命名空间的文档中概述的张量积或组合策略。因此，这个函数没有
     * @p do_tensor_product 参数。
     *
     */
    template <int dim, int spacedim>
    std::vector<bool>
    compute_restriction_is_additive_flags(
      const FiniteElement<dim, spacedim> *fe1,
      const unsigned int                  N1,
      const FiniteElement<dim, spacedim> *fe2 = nullptr,
      const unsigned int                  N2  = 0,
      const FiniteElement<dim, spacedim> *fe3 = nullptr,
      const unsigned int                  N3  = 0,
      const FiniteElement<dim, spacedim> *fe4 = nullptr,
      const unsigned int                  N4  = 0,
      const FiniteElement<dim, spacedim> *fe5 = nullptr,
      const unsigned int                  N5  = 0);


    /**
     * 计算由有限元列表描述的组成有限元的每个形状函数的非零分量，其乘数在第二个参数中给出。
     * 如果 @p do_tensor_product
     * 为真，则分量的数量（以及ComponentMask对象的大小）是每个有限元中的分量数量乘以相应倍数的乘积之和。
     * 否则，组件数从第一个具有非零倍数的有限元中提取，所有其他具有非零倍数的元素都需要有相同数量的向量组件。
     * 关于 @p do_tensor_product
     * 参数的更多信息，请参见命名空间 FETools::Compositing
     * 的文档。
     *
     */
    template <int dim, int spacedim>
    std::vector<ComponentMask>
    compute_nonzero_components(
      const std::vector<const FiniteElement<dim, spacedim> *> &fes,
      const std::vector<unsigned int> &                        multiplicities,
      const bool do_tensor_product = true);

    /**
     * 对于任意数量的
     * <code>std::pair<std::unique_ptr<FiniteElement<dim,  spacedim>>,
     * unsigned int></code>和 <code>do_tensor_product = true</code>
     * 类型的参数，与上述相同。
     *
     */
    template <int dim, int spacedim>
    std::vector<ComponentMask>
    compute_nonzero_components(
      const std::initializer_list<
        std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>>
        &fe_systems);

    /**
     * 计算组成有限元的非零向量分量。这个函数与前一个函数类似，只是指针表示要组成的元素，而参数
     * @p N1,   @p N2,
     * ......是乘数。空指针表示一个参数将被跳过。
     * 如果 @p do_tensor_product
     * 为真，组件的数量（以及ComponentMask对象的大小）是每个有限元素中的组件数量乘以相应的倍数的乘积之和。
     * 否则，组件数从第一个具有非零倍数的有限元中提取，所有其他具有非零倍数的元素都需要有相同数量的向量组件。
     * 关于 FETools::Compositing
     * 参数的更多信息，请参见命名空间 @p do_tensor_product
     * 的文档。
     *
     */
    template <int dim, int spacedim>
    std::vector<ComponentMask>
    compute_nonzero_components(
      const FiniteElement<dim, spacedim> *fe1,
      const unsigned int                  N1,
      const FiniteElement<dim, spacedim> *fe2               = nullptr,
      const unsigned int                  N2                = 0,
      const FiniteElement<dim, spacedim> *fe3               = nullptr,
      const unsigned int                  N3                = 0,
      const FiniteElement<dim, spacedim> *fe4               = nullptr,
      const unsigned int                  N4                = 0,
      const FiniteElement<dim, spacedim> *fe5               = nullptr,
      const unsigned int                  N5                = 0,
      const bool                          do_tensor_product = true);

    /**
     * 对于一个给定的（复合） @p finite_element 构建 @p
     * 系统_to_component_table， @p system_to_base_table 和 @p
     * component_to_base_table。        如果 @p do_tensor_product
     * 为真，用于复合元素的组件数是每个有限元素中的组件数乘以相应的倍数的乘积之和。
     * 否则，组件数取自第一个具有非零倍数的有限元，所有其他具有非零倍数的元素需要有相同数量的矢量组件。
     * 关于 FETools::Compositing
     * 参数的更多信息，请参见命名空间 @p do_tensor_product
     * 的文档。
     *
     */
    template <int dim, int spacedim>
    void
    build_cell_tables(
      std::vector<std::pair<std::pair<unsigned int, unsigned int>,
                            unsigned int>> &system_to_base_table,
      std::vector<std::pair<unsigned int, unsigned int>>
        &                                   system_to_component_table,
      std::vector<std::pair<std::pair<unsigned int, unsigned int>,
                            unsigned int>> &component_to_base_table,
      const FiniteElement<dim, spacedim> &  finite_element,
      const bool                            do_tensor_product = true);

    /**
     * 对于给定的（复合） @p finite_element 建立 @p
     * face_system_to_base_table, 和 @p face_system_to_component_table.
     * 如果 @p do_tensor_product
     * 为真，用于复合元素的组件数是每个有限元素中的组件数乘以相应倍数的乘积之和。
     * 否则，组件数取自第一个具有非零倍数的有限元，所有其他具有非零倍数的元素需要有相同数量的矢量组件。
     * 关于 @p do_tensor_product
     * 参数的更多信息，请参见命名空间 FETools::Compositing
     * 的文档。
     *
     */
    template <int dim, int spacedim>
    void
    build_face_tables(
      std::vector<std::pair<std::pair<unsigned int, unsigned int>,
                            unsigned int>> &face_system_to_base_table,
      std::vector<std::pair<unsigned int, unsigned int>>
        &                                 face_system_to_component_table,
      const FiniteElement<dim, spacedim> &finite_element,
      const bool                          do_tensor_product = true,
      const unsigned int                  face_no           = 0  /*TODO*/ );

  } // namespace Compositing


  /**
   * 解析有限元的名称并相应地生成一个有限元对象。解析器忽略单词之间的空格字符（与正则表达式[A-Za-z0-9_]相匹配的东西）。
   * 名称必须是由 FiniteElement::get_name
   * 函数返回的形式，其中尺寸模板参数&lt;2&gt;等可以省略。另外，明确的数字可以用<tt>dim</tt>或<tt>d</tt>代替。如果给出一个数字，它<b>must</b>与该函数的模板参数相匹配。
   * FESystem元素的名称遵循
   * <code>FESystem[FE_Base1^p1-FE_Base2^p2]</code> The powers <code>p1</code>
   * 的模式，可以是数字，也可以用<tt>dim</tt>或<tt>d</tt>代替。
   * 如果不能从这个字符串中重建有限元，将抛出一个 @p
   * FETools::ExcInvalidFEName 类型的异常。    该函数返回一个
   * std::unique_ptr
   * 的新创建的有限元，意味着调用者获得了对返回对象的所有权。
   * 由于模板参数的值不能从给这个函数的（字符串）参数中推导出来，你必须在调用这个函数时明确指定它。
   * 这个函数知道库中定义的所有标准元素。然而，它默认不知道你可能在你的程序中定义的元素。要使这个函数知道你自己的元素，请使用add_fe_name()函数。
   * 如果想得到一个模数维为1的有限元，这个函数就不起作用。
   *
   */
  template <int dim, int spacedim = dim>
  std::unique_ptr<FiniteElement<dim, spacedim>>
  get_fe_by_name(const std::string &name);

  /**
   * 用给定的 @p name.
   * 扩展可由get_fe_by_name()生成的有限元列表，如果以后用这个名字调用get_fe_by_name()，它将使用给定的作为第二个参数的对象来创建一个有限元对象。
   * @p name
   * 参数的格式应包括有限元的名称。然而，单独使用类名或使用
   * FiniteElement::get_name
   * 的结果（包括空间维度以及多项式程度）都是安全的，因为第一个非名称字符之后的所有内容将被忽略。
   * FEFactory对象应该是一个用<tt>new</tt>新建的对象。
   * FETools将拥有这个对象的所有权，一旦不再使用它，就将其删除。
   * 在大多数情况下，如果你想在给get_fe_by_name提供名称
   * <code>my_fe</code> 时创建 <code>MyFE</code>
   * 类型的对象，你会希望这个函数的第二个参数是FEFactory
   * @<MyFE@>,
   * 类型，但你当然可以创建你的自定义有限元工厂类。
   * 这个函数接管了作为第二个参数的对象的所有权，也就是说，你不应该试图在以后销毁它。该对象将在程序的生命周期结束时被删除。
   * 如果该元素的名称已经被使用，则会抛出一个异常。
   * 因此，get_fe_by_name()的功能只能被添加，不能被改变。
   * @note  这个函数操作一个全局表（每个空间维度有一个表）。它是线程安全的，因为对这个表的每一次访问都有一个锁来保证。然而，由于每个名字只能被添加一次，用户代码必须确保只有一个线程添加一个新元素。    还要注意，这个表对每个空间维度都存在一次。如果你有一个在不同空间维度上处理有限元的程序（例如， @ref step_4  "  step-4  "
   * 做这样的事情），那么你应该为每个空间维度调用这个函数，你希望你的有限元被添加到地图中。
   *
   */
  template <int dim, int spacedim>
  void
  add_fe_name(const std::string &                 name,
              const FEFactoryBase<dim, spacedim> *factory);

  /**
   * 用于get_fe_by_name()的字符串不能被翻译成有限元。
   * 要么是字符串的格式不好，要么是你使用的是自定义元素，必须先用add_fe_name()添加。
   * @ingroup Exceptions
   *
   */
  DeclException1(ExcInvalidFEName,
                 std::string,
                 << "Can't re-generate a finite element from the string '"
                 << arg1 << "'.");

  /**
   * 用于get_fe_by_name()的字符串不能被翻译成有限元。
   * 应该避免在有限元名称中出现尺寸参数。如果有的话，尺寸应该是<tt>dim</tt>或者<tt>d</tt>。这里，你给了一个数字尺寸参数，这与有限元类的模板尺寸不匹配。
   * @ingroup Exceptions
   *
   */
  DeclException2(ExcInvalidFEDimension,
                 char,
                 int,
                 << "The dimension " << arg1
                 << " in the finite element string must match "
                 << "the space dimension " << arg2 << ".");

  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcInvalidFE);

  /**
   * 有限元必须是 @ref GlossPrimitive "原始的"
   * 。
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcFENotPrimitive);
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcTriangulationMismatch);

  /**
   * 在有悬挂节点的网格上使用连续元素，但缺少约束矩阵。
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(ExcHangingNodesNotAllowed,
                   "You are using continuous elements on a grid with "
                   "hanging nodes but without providing hanging node "
                   "constraints. Use the respective function with "
                   "additional AffineConstraints argument(s), instead.");
  /**
   * 你需要至少两个网格级别。
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcGridNotRefinedAtLeastOnce);
  /**
   * 所用矩阵的尺寸与预期的尺寸不一致。
   * @ingroup Exceptions
   *
   */
  DeclException4(ExcMatrixDimensionMismatch,
                 int,
                 int,
                 int,
                 int,
                 << "This is a " << arg1 << "x" << arg2 << " matrix, "
                 << "but should be a " << arg3 << "x" << arg4 << " matrix.");

  /**
   * 如果嵌入矩阵的计算不准确，则抛出异常。
   * @ingroup Exceptions
   *
   */
  DeclException1(ExcLeastSquaresError,
                 double,
                 << "Least squares fit leaves a gap of " << arg1);

  /**
   * 如果一个变量可能不大于另一个变量，则抛出异常。
   * @ingroup Exceptions
   *
   */
  DeclException2(ExcNotGreaterThan,
                 int,
                 int,
                 << arg1 << " must be greater than " << arg2);
} // namespace FETools


#ifndef DOXYGEN

namespace FETools
{
  template <class FE>
  std::unique_ptr<FiniteElement<FE::dimension, FE::space_dimension>>
  FEFactory<FE>::get(const unsigned int degree) const
  {
    return std::make_unique<FE>(degree);
  }

  namespace Compositing
  {
    template <int dim, int spacedim>
    std::vector<bool>
    compute_restriction_is_additive_flags(
      const std::initializer_list<
        std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>>
        &fe_systems)
    {
      std::vector<const FiniteElement<dim, spacedim> *> fes;
      std::vector<unsigned int>                         multiplicities;

      const auto extract =
        [&fes, &multiplicities](
          const std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
                          unsigned int> &fe_system) {
          fes.push_back(fe_system.first.get());
          multiplicities.push_back(fe_system.second);
        };

      for (const auto &p : fe_systems)
        extract(p);

      return compute_restriction_is_additive_flags(fes, multiplicities);
    }



    template <int dim, int spacedim>
    FiniteElementData<dim>
    multiply_dof_numbers(
      const std::initializer_list<
        std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>>
        &fe_systems)
    {
      std::vector<const FiniteElement<dim, spacedim> *> fes;
      std::vector<unsigned int>                         multiplicities;

      const auto extract =
        [&fes, &multiplicities](
          const std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
                          unsigned int> &fe_system) {
          fes.push_back(fe_system.first.get());
          multiplicities.push_back(fe_system.second);
        };

      for (const auto &p : fe_systems)
        extract(p);

      return multiply_dof_numbers(fes, multiplicities, true);
    }



    template <int dim, int spacedim>
    std::vector<ComponentMask>
    compute_nonzero_components(
      const std::initializer_list<
        std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>>
        &fe_systems)
    {
      std::vector<const FiniteElement<dim, spacedim> *> fes;
      std::vector<unsigned int>                         multiplicities;

      const auto extract =
        [&fes, &multiplicities](
          const std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
                          unsigned int> &fe_system) {
          fes.push_back(fe_system.first.get());
          multiplicities.push_back(fe_system.second);
        };

      for (const auto &p : fe_systems)
        extract(p);

      return compute_nonzero_components(fes, multiplicities, true);
    }
  } // namespace Compositing
} // namespace FETools

#endif

 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif  /* dealii_fe_tools_H */ 


