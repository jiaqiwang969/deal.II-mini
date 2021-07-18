//include/deal.II-translator/fe/fe_enriched_0.txt
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

#ifndef dealii_fe_enriched_h
#define dealii_fe_enriched_h

#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/symmetric_tensor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_update_flags.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>

#include <map>
#include <numeric>
#include <utility>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 实现了Babuska和Melenk的分区统一有限元方法（PUM），该方法将标准有限元与另一个（通常是线性）有限元相乘，丰富了标准有限元的内容。\f[
 * U(\mathbf x) = \sum_i N_i(\mathbf x) U_i + \sum_j N_j(\mathbf x) \sum_k
 * F_k(\mathbf x) U_{jk} \f]其中 $ N_i(\mathbf x) $ 和 $ N_j(\mathbf x) $
 * 是基础有限元（包括从等参元到实元的映射）； $
 * F_k(\mathbf x) $ 是实空间的标量丰富函数（例如 $ 1/r $  ， $
 * \exp(-r) $  ，等等）； $ U_i $ 和 $ U_{jk} $
 * 是标准和丰富的DoF。这允许在有限元空间中包括关于被解决的偏微分方程的先验知识，这反过来又改善了空间的局部近似特性。这对高度振荡解、有域角或无界域的问题或边界条件的突然变化很有用。PUM方法使用满足统一分区属性的有限元空间（如FE_Q）。在其他属性中，这使得所产生的空间能够准确地重现富集函数。
 * 这个类的最简单的构造函数需要两个有限元对象和一个要使用的富集函数。比如说
 *
 *
 * @code
 * FE_Enriched<dim> fe(FE_Q<dim>(2),
 *                   FE_Q<dim>(1),
 *                   function)
 * @endcode
 *
 * 在这种情况下，标准的DoFs由 <code>FE_Q<dim>(2)</code>
 * 分配，而丰富的DoFs来自一个单一的有限元
 * <code>FE_Q<dim>(1)</code> ，使用一个丰富的函数
 * <code>function</code>
 * 。在这种情况下，丰富元素上的总DoF数是来自
 * <code>FE_Q<dim>(2)</code> 和 <code>FE_Q<dim>(1)</code> 的DoF数之和。
 * 作为丰富函数的一个例子，考虑 $ \exp(-x) $
 * ，它导致了单位元素上的以下形状函数。  <table> <tr> <td
 * align="center">
       @image html fe_enriched_1d.png
 * </td> <td align="center">
       @image html fe_enriched_h-refinement.png
 * </td> </tr> <tr> <td align="center"> 1d element, base and enriched shape
 * functions. </td> <td align="center"> enriched shape function corresponding
 * to the central vertex. </td> </tr> </table>
 * 请注意，对富集形状函数或有限元场的梯度（Hessians）的评估需要对富集函数的梯度（梯度和Hessians）进行评估。
 *
 * @f{align*}{
 * U(\mathbf x)
 *   &= \sum_i N_i(\mathbf x) U_i
 *   + \sum_{j,k} N_j(\mathbf x) F_k(\mathbf x) U_{jk} \\
 * \mathbf \nabla U(\mathbf x)
 *   &= \sum_i \mathbf \nabla N_i(\mathbf x) U_i
 *   + \sum_{j,k} \left[\mathbf \nabla N_j(\mathbf x) F_k(\mathbf x) +
 *                      N_j(\mathbf x) \mathbf \nabla F_k(\mathbf x) \right]
 * U_{jk} \\ \mathbf \nabla \mathbf \nabla U(\mathbf x)
 *   &= \sum_i \mathbf \nabla \mathbf \nabla N_i(\mathbf x) U_i
 *   + \sum_{j,k} \left[\mathbf \nabla \mathbf \nabla N_j(\mathbf x)
 * F_k(\mathbf x) + \mathbf \nabla F_k(\mathbf x) \mathbf \nabla N_j(\mathbf x)
 * + \mathbf \nabla N_j(\mathbf x) \mathbf \nabla F_k(\mathbf x) + N_j(\mathbf
 * x) \mathbf \nabla \mathbf \nabla F_k(\mathbf x) \right] U_{jk}
 * @f}
 *
 * <h3>Using enriched and non-enriched FEs together</h3>
 * 在大多数应用中，只在领域的某些部分（如裂纹尖端周围）引入富集函数，而在其他地方使用标准FE（如FE_Q）是有益的。这可以通过使用deal.II中的hp-finite
 * element框架来实现，该框架允许在不同单元上使用不同的元素。为了使得到的空间
 * $C^0$ 连续，那么就需要DoFHandler类和
 * DoFTools::make_hanging_node_constraints()
 * 函数能够弄清楚在富集单元和非富集单元之间的界面上该做什么。具体来说，我们希望对应于富集形状函数的自由度在这些界面上为零。这些类和函数不能自动做到这一点，但是可以通过在没有富集的细胞上使用普通的FE_Q来实现这个效果，而是将FE_Q包装成FE_Enriched对象<i>without
 * actually enriching it</i>。这可以按以下方法进行。
 *
 * @code
 * FE_Enriched<dim> fe_non_enriched(FE_Q<dim>(1));
 * @endcode
 * 这个构造函数等同于调用
 *
 * @code
 * FE_Enriched<dim> fe_non_enriched(FE_Q<dim>(1),
 *                                  FE_Nothing<dim>(1,true),
 *                                  nullptr);
 * @endcode
 * 并将产生正确的约束，用于归属于两个区域之间界面上的支持点的丰富的DoF。
 * <h3>References</h3> 当使用这个类时，请引用  @cite davydov2017hp
 * 。PUM是在  @cite melenk1996  和  @cite babuska1997  中引入的。
 * <h3>Implementation</h3>
 * 该类的实现是基于FESystem的，它被聚合为一个私有成员。最简单的构造函数
 * <code> FE_Enriched<dim> fe(FE_Q<dim>(2),
 * FE_Q<dim>(1),function)</code>将内部初始化FESystem为
 *
 *
 * @code
 * FESystem<dim> fe_system(FE_Q<dim>(2),1,
 *                       FE_Q<dim>(1),1);
 * @endcode
 *
 * 请注意，让这个类派生自FESystem是不明智的，因为后者将给定的元素串联成一个向量元素的不同组成部分，而当前的类将给定的元素合并成相同的组成部分。例如，如果给定两个标量元素，当用FESystem做同样的事情时，产生的元素将是标量的，而不是有两个分量。
 * 形状函数、 @p interface_constrains, 、 @p prolongation
 * （嵌入）和 @p restriction 矩阵的排序取自FESystem类。
 *
 *
 * @ingroup fe
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_Enriched : public FiniteElement<dim, spacedim>
{
public:
  /**
   * 构造函数，接收基础有限元素 @p fe_base
   * 和丰富的有限元素 @p fe_enriched ，它将被乘以 @p
   * enrichment_function. ，如果 @p fe_enriched 不是FE_Nothing， @p
   * enrichment_function 的寿命必须至少和FE_Enriched对象一样长。
   *
   */
  FE_Enriched(const FiniteElement<dim, spacedim> &fe_base,
              const FiniteElement<dim, spacedim> &fe_enriched,
              const Function<spacedim> *          enrichment_function);

  /**
   * 构造函数，它只包装基础FE  @p fe_base.
   * 至于丰富的有限元空间，则使用FE_Nothing。
   * 当这个非富集元素与富集有限元素在具有hp-capabilities的DoFHandler中一起使用时，将自动生成连续性约束。
   * 关于如何在hp-finite
   * element方法的背景下使用这个元素，请参见类文档中的讨论。
   *
   */
  FE_Enriched(const FiniteElement<dim, spacedim> &fe_base);

  /**
   * 构造函数，接收指向基础FiniteElement  @p fe_base
   * 的指针和丰富的FiniteElement的向量  @p fe_enriched  。  @p
   * fe_enriched[i]  有限元将用 @p functions[i].
   * 中的函数来充实。这是最通用的公共构造函数，它也允许在域的不同不相干部分有不同的充实函数。
   * 为此，最后一个参数提供了一个单元格迭代器与一个函数的关联。这样做是为了简化这个类的使用，当具有不同功能的互不相连的域的数量超过几个时。
   * 否则，我们将不得不为每个互不相干的富集域使用这个类的不同实例。
   * 如果你不打算使用这个功能，你可以利用C++11的lambdas来定义假函数。下面是一个例子，在第一个元素上使用两个函数，在第二个元素上使用一个函数。
   * @code
   * FE_Enriched<dim> fe(
   * &fe_base,
   * {&fe_1, &fe_2},
   * {{[=] (const typename Triangulation<dim>::cell_iterator &)
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
   * -> const Function<dim>
   *   {
   *     return &fe_1_function1;
   *   },
   *   [=] (const typename Triangulation<dim>::cell_iterator &)
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
   * -> const Function<dim>
   *   {
   *     return &fe_1_function2;
   *   }},
   *  {[=] (const typename Triangulation<dim>::cell_iterator &)
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
   * -> const Function<dim>
   *   {
   *     return &fe_2_function;
   *   }
   *  }});
   * @endcode
   *
   * @note
   * 当使用同一个有限元与N个不同的函数进行充实时，建议第二个参数的大小为1，最后一个参数的大小为1×N，同样可以通过提供N个相等的充实元素而保持最后一个参数的大小为N×1来实现。
   * @note
   * 当在不相交的域上使用不同的充实函数时，该类没有检查这些域是否真的不相交。
   *
   */
  FE_Enriched(
    const FiniteElement<dim, spacedim> *                     fe_base,
    const std::vector<const FiniteElement<dim, spacedim> *> &fe_enriched,
    const std::vector<std::vector<std::function<const Function<spacedim> *(
      const typename Triangulation<dim, spacedim>::cell_iterator &)>>>
      &functions);

private:
  /**
   * 最一般的私有构造函数。前两个输入参数与FESystem中的参数一致。它在内部只与
   * <code>multiplicities[0]=1</code>
   * 一起使用，这是对这个有限元的逻辑要求。
   *
   */
  FE_Enriched(
    const std::vector<const FiniteElement<dim, spacedim> *> &fes,
    const std::vector<unsigned int> &                        multiplicities,
    const std::vector<std::vector<std::function<const Function<spacedim> *(
      const typename Triangulation<dim, spacedim>::cell_iterator &)>>>
      &functions);

public:
  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  /**
   * 返回一个标识有限元的字符串。
   *
   */
  virtual std::string
  get_name() const override;

  /**
   * 访问一个组成元素。索引需要小于基础元素的数量。在这个类的上下文中，基础元素的数量总是多于一个：一个非富集元素加上一个要富集的元素，这可能是FE_Nothing。
   *
   */
  virtual const FiniteElement<dim, spacedim> &
  base_element(const unsigned int index) const override;

  /**
   * 返回 @p ith 形状函数在 @p p.  @p p
   * 点的值，该点是参考元素上的一个点。
   * 这个函数只对非富集元素返回有意义的值，因为实空间的富集需要在实空间的点上对函数进行评估。
   *
   */
  virtual double
  shape_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * @name 转移矩阵 
     * @{ 
   *
   */

  /**
   * 从细网格空间投射到粗网格空间。
   * 这个功能只有在所有的子元素也使用与父元素相同的函数进行充实时才有意义。
   *
   */
  virtual const FullMatrix<double> &
  get_restriction_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * 网格之间的嵌入矩阵。
   * 这个函数只有在所有的子元素也使用与父元素相同的函数进行充实时才有意义。
   *
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  //@}

  /**
   * @name 支持hp的函数  
     * @{ 
   *
   */

  /**
   * 返回这个元素是否实现了hp-constraints。
   * 当且仅当其所有基础元素都返回 @p true
   * 时，此函数才会返回 @p true 。
   *
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * 返回从一个元素的面插值到邻近元素的面的矩阵。
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。
   * 这个元素的基础元素将不得不实现这个功能。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现给定元素的插值，那么他们必须抛出一个类型为
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented,
   * 的异常，这个异常将从这个元素传播出去。
   *
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double> &                matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * 返回从一个元素的面插值到邻近元素的子面的矩阵。
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。
   * 这个元素的基础元素将不得不实现这个功能。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现从一个给定元素的内插，那么他们必须抛出一个类型为
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented,
   * 的异常，这个异常将从这个元素传播出去。
   *
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source,
    const unsigned int                  subface,
    FullMatrix<double> &                matrix,
    const unsigned int                  face_no = 0) const override;

  /**
   * 如果在一个顶点上，有几个有限元被激活，hp-code首先为这些FEs的自由度分配不同的全局索引。然后调用这个函数来找出其中哪些应该得到相同的值，从而可以得到相同的全局自由度指数。
   * 因此，该函数返回当前有限元对象的自由度与 @p fe_other,
   * 的自由度之间的相同性列表，后者是对代表在该特定顶点上活动的其他有限元之一的有限元对象的引用。该函数计算两个有限元对象的哪些自由度是相等的，这两个自由度的编号都在零和两个有限元的n_dofs_per_vertex()的相应值之间。每一对的第一个索引表示本元素的一个顶点自由度，而第二个是另一个有限元素的相应索引。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理线上自由度。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理四边形上的自由度。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  /**
   * @copydoc   FiniteElement::compare_for_domination() .
   *
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

  //@}


  /**
   * 返回充实函数
   *
   */
  const std::vector<std::vector<std::function<const Function<spacedim> *(
    const typename Triangulation<dim, spacedim>::cell_iterator &)>>>
  get_enrichments() const;

  /**
   * 返回底层的FESystem对象。
   *
   */
  const FESystem<dim, spacedim> &
  get_fe_system() const;

protected:
  /**
   * 一个用于保存在正交点评估该FE所需的内部数据的类。
   *
   */
  class InternalData : public FiniteElement<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * 对于每个有限元(基数)和每个增强函数(base_index)，这个结构将包含增强函数的值、梯度和希斯系数。
     *
     */
    struct EnrichmentValues
    {
      std::vector<double>                       values;
      std::vector<Tensor<1, spacedim>>          gradients;
      std::vector<SymmetricTensor<2, spacedim>> hessians;
    };

    /**
     * 构造器。在setup_data中使用，用于包装FESystem的内部数据对象。前者是由FE_Enriched必须实现的get_data、get_subface_data和get_face_data所调用。
     * 由于 FESystem::get_data(),  FESystem::get_face_data() 和
     * FESystem::get_subface_data()
     * 只是创建一个对象并返回一个指针（也就是说，它们不保留所有权），我们将投射结果存储在
     * std::unique_ptr 中，以表明InternalData拥有该对象。
     *
     */
    InternalData(std::unique_ptr<typename FESystem<dim, spacedim>::InternalData>
                   fesystem_data);

    /**
     * 给予FESystem数据的 @p  <code>base_no</code>
     * 第1个基本元素的指针的读访问权。
     *
     */
    typename FiniteElement<dim, spacedim>::InternalDataBase &
    get_fe_data(const unsigned int base_no) const;

    /**
     * 当调用 FiniteElement::fill_fe_values() 和类似的函数时，对
     * <code>base_no</code>
     * 第1个基本元素将写入其输出的对象的指针给予读访问权。
     *
     */
    internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &
    get_fe_output_object(const unsigned int base_no) const;

    /**
     * 聚合FESystem的内部数据。每次我们调用FESystem的fill_fe_values()和类似的功能时都会用到它。
     *
     */
    std::unique_ptr<typename FESystem<dim, spacedim>::InternalData>
      fesystem_data;

    /**
     * 对于每个用于富集的FE（基数  <code>i</code>
     * ）和每个富集函数（基数倍数  <code>j</code>  ），
     * <code>enrichment_values[i][j]</code>
     * 将被用来存储可能要求的值、梯度和富集函数的
     * hessians  <code>j</code>  。
     * 由于提供给fill_fe_values的InternalData是常数，所以该变量是可变的。
     * @note
     * 我们不希望将这些信息存储在有限元对象本身，因为这意味着(i)一次只能有一个FEValues对象使用一个有限元对象，(ii)这些对象不能在多线程环境下使用。
     *
     */
    mutable std::vector<std::vector<EnrichmentValues>> enrichment;
  };

  /**
   * 对于在富集中使用的每个有限元 @p i
   * 和与之相关的每个富集函数 @p j （基本上是其倍数），
   * @p base_no_mult_local_enriched_dofs[i][j]
   * 包含FE_Enriched有限元上的相关局部DoF。
   *
   */
  std::vector<std::vector<std::vector<unsigned int>>>
    base_no_mult_local_enriched_dofs;

  /**
   * 丰富的函数。
   * 第一个向量的大小与使用富集的有限元空间的数量相同。而内向量的大小对应于与单个FiniteElement相关的富集函数的数量。
   *
   */
  const std::vector<std::vector<std::function<const Function<spacedim> *(
    const typename Triangulation<dim, spacedim>::cell_iterator &)>>>
    enrichments;

  /**
   * 辅助变量，用于区分我们做丰富化的情况和类简单包裹另一个FiniteElement的情况。
   * 这个变量在构造函数中被初始化，方法是在一个丰富元素的向量上循环，并检查它们是否都是FE_Nothing。如果是这样，那么值就被设置为
   * <code>false</code>  ，否则就是 <code>true</code>  。
   *
   */
  const bool is_enriched;

  /**
   * 从get_data、get_face_data和get_subface_data调用的辅助函数。它获取
   * @p fes_data 中FESystem对象的内部数据和正交规则 @p qudrature.
   * 。这个函数基本上是从FESystem类的一个实例中获取内部数据，并将其封装到我们自己的InternalData类中，该类中还有一些对象用于保存每个正交点的富集函数的值/梯度/希斯值，取决于
   * @p flags.  。
   *
   */
  template <int dim_1>
  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
  setup_data(
    std::unique_ptr<typename FESystem<dim, spacedim>::InternalData> fes_data,
    const UpdateFlags                                               flags,
    const Quadrature<dim_1> &quadrature) const;

  /**
   * 准备内部数据结构并填入独立于单元格的值。返回一个指向对象的指针，然后该函数的调用者（FEValues）必须承担该对象的所有权（这包括在不再需要时的销毁）。
   *
   */
  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_data(
    const UpdateFlags             flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim> &       quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  using FiniteElement<dim, spacedim>::get_face_data;

  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_face_data(
    const UpdateFlags               update_flags,
    const Mapping<dim, spacedim> &  mapping,
    const hp::QCollection<dim - 1> &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_subface_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim - 1> &   quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  using FiniteElement<dim, spacedim>::fill_fe_face_values;

  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const hp::QCollection<dim - 1> &                            quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          sub_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

private:
  /**
   * 这个函数设置了系统的索引表，以及 @p 限制和 @p
   * prolongation 矩阵。
   *
   */
  void
  initialize(const std::vector<const FiniteElement<dim, spacedim> *> &fes,
             const std::vector<unsigned int> &multiplicities);

  /**
   * 底层的FESystem对象。
   *
   */
  const std::unique_ptr<const FESystem<dim, spacedim>> fe_system;

  /**
   * 在调用fill_fe_(face/subface_)values后，这个函数实现了链式规则，将存储的形状值/梯度/黑度乘以在正交点评估的富集函数的值。
   *
   */
  template <int dim_1>
  void
  multiply_by_enrichment(
    const Quadrature<dim_1> &quadrature,
    const InternalData &     fe_data,
    const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &                                                         mapping_data,
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
      &output_data) const;
};



/**
 * 该命名空间包括一个创建FE_Enriched有限元 (hp::FECollection)
 * 集合所需的类，该集合可与DoFHandler在hp模式下用于具有多个可能重叠的子域的单独富集函数。
 * 为了创建 hp::FECollection, ，在创建 hp::FECollection.
 * 之前，使用了一种图形着色算法来为富集函数分配颜色，因此而得名。
 *
 *
 */
namespace ColorEnriched
{
  /**
   * 谓词函数的别名模板，它为一个
   * Triangulation<dim,spacedim>::cell_iterator 对象返回一个布尔值。
   * 这被辅助函数和 ColorEnriched::Helper 类的实现所使用。
   *
   */
  template <int dim, int spacedim = dim>
  using predicate_function = std::function<bool(
    const typename Triangulation<dim, spacedim>::cell_iterator &)>;

#ifndef DOXYGEN
  namespace internal
  {
    /**
     * 如果与 @p dof_handler
     * 相关的网格中的子域之间存在连接，即，如果子域至少共享一个顶点，则返回true。两个子域由
     * @p predicate_1 和 @p predicate_2.
     * 提供的谓词定义。谓词是一个函数（或带有operator()的类型的对象），它接收一个单元格迭代器并给出一个布尔值。如果它的返回值为
     * "真"，则表示它在一个单元格中是有效的。
     * 一个自定义谓词的例子是检查与一个固定点的距离。请注意，operator()接收了一个单元格迭代器。使用构造函数，可以选择固定点和距离。
     * @code
     * <int dim>
     * struct predicate
     * {
     *   predicate(const Point<dim> p, const int radius)
     *   :p(p),radius(radius){}
     *
     *   template <class Iterator>
     *   bool operator () (const Iterator &i)
     *   {
     *       return ( (i->center()
     *
     * - p).norm() < radius);
     *   }
     *
     * private:
     *   Point<dim> p;
     *   int radius;
     *
     * };
     * @endcode
     * 然后可以用下面的函数来寻找子域是否相连。
     * @code
     * find_connection_between_subdomains
     * (dof_handler,
     * predicate<dim>(Point<dim>(0,0), 1)
     * predicate<dim>(Point<dim>(2,2), 1));
     * @endcode
     * @param[in]  dof_handler DoFHandler对象  @param[in]  predicate_1
     * 一个定义子域1的函数（或带有operator()的类型的对象）。
     * 该函数接收一个单元格并返回一个布尔值。
     * @param[in]  predicate_2 与 @p predicate_1
     * 相同，但定义了子域2。  @return
     * 如果子域至少共享一个顶点，则是一个布尔值 "真"。
     *
     */
    template <int dim, int spacedim>
    bool
    find_connection_between_subdomains(
      const DoFHandler<dim, spacedim> &        dof_handler,
      const predicate_function<dim, spacedim> &predicate_1,
      const predicate_function<dim, spacedim> &predicate_2);

    /**
     * 使用图形着色算法为子域分配颜色，每个子域被视为一个图形节点。连接的子域，即至少共享一个顶点的子域有不同的颜色。每个子域都是用
     * @p predicates.   @param[in]  dof_handler 一个DoFHandler对象
     * @param[in]  predicates 定义子域的谓词  @param[out]
     * predicate_colors 与每个子域有关的颜色（无符号int）。
     *
     */
    template <int dim, int spacedim>
    unsigned int
    color_predicates(
      const DoFHandler<dim, spacedim> &                     dof_handler,
      const std::vector<predicate_function<dim, spacedim>> &predicates,
      std::vector<unsigned int> &                           predicate_colors);

    /**
     * 用于构造Helper类的数据成员 @p cellwise_color_predicate_map
     * 和 @p fe_sets
     * 。输入是DoFHandler对象，谓词的向量和与之相关的颜色。在调用这个函数之前，可以使用函数color_predicates将颜色分配给谓词（即子域）。
     * 每个活跃的FE索引都有一组与之相关的颜色。
     * 一个具有活动FE索引i的单元格有一组颜色，由
     * <code>fe_sets[i]</code>
     * 给出。一个具有{a,b}颜色的活动FE索引意味着该单元有两个颜色为a和b的活动谓词（即它们对该单元返回真）。}意味着具有活动FE索引0的单元格没有相关的谓词。
     * 索引为1的单元格有一个颜色为1的活动谓词。
     * 索引为2的单元格有一个颜色为2的活动谓词。
     * 索引为3的单元格有颜色为1和颜色为2的活动谓词。
     * 一个地图cellwise_color_predicate_map被用来将单元格中的predicate颜色与predicate
     * id联系起来。为此，每个单元格都被赋予一个唯一的id，这个id现在被存储在材料id中。
     * 当网格被细化时，材料ID会被继承到子代，所以将材料ID与颜色地图联系起来的地图仍将是相关的。
     * 现在我们可以用一个例子来解释颜色地图。如果材料ID为100的单元格有活动的谓词4（颜色=1）和5（颜色=2），地图将在键100处插入对（1，4）和（2，5）（即单元格的唯一ID被映射到将颜色与谓词ID关联的地图上）。
     * @param[in]  dof_handler DoFHandler对象  @param[in]
     * 定义子域的谓词向量。      <code>@p predicates[i]</code>
     * 如果一个单元格属于索引为i的子域，则返回true。
     * @param[in]  predicate_colors
     * 与每个子域相关的颜色向量（无符号int）。
     * @param[out]  cellwise_color_predicate_map
     * 用于将单元格中的谓词颜色与谓词id相关联的地图。
     * @param[out]  fe_sets 一个颜色列表的向量
     *
     */
    template <int dim, int spacedim>
    void
    set_cellwise_color_set_and_fe_index(
      DoFHandler<dim, spacedim> &                           dof_handler,
      const std::vector<predicate_function<dim, spacedim>> &predicates,
      const std::vector<unsigned int> &                     predicate_colors,
      std::map<unsigned int, std::map<unsigned int, unsigned int>>
        &                                  cellwise_color_predicate_map,
      std::vector<std::set<unsigned int>> &fe_sets);

    /**
     * 一个函数，用于返回与颜色相对应的富集函数的向量。该向量的大小等于与谓词（即子域）相关的不同颜色的总数。
     * 假设一个细胞有一个id为4（颜色=1）和5（颜色=2）的活动谓词。只要我们知道材料id，cellwise_color_predicate_map就有这个信息。
     * 构建的color_enrichments是这样的：color_enrichments[color=1](cell)将返回一个指向id=4的富集函数的指针，即enrichments[4]。
     * 换句话说，使用这个函数中先前收集的信息，我们将用户提供的丰富函数的向量翻译成适合FE_Enriched类的函数向量。
     * @param[in]  n_colors 谓词的颜色数量  @param[in]  enrichments
     * 丰富函数的向量  @param[in]  cellwise_color_predicate_map
     * 用于将单元格中的谓词颜色与谓词id联系起来的地图。
     * @param[out]  color_enrichments
     * 一个接受单元格并返回函数指针的函数向量。
     *
     */
    template <int dim, int spacedim>
    void
    make_colorwise_enrichment_functions(
      const unsigned int                                      n_colors,
      const std::vector<std::shared_ptr<Function<spacedim>>> &enrichments,
      const std::map<unsigned int, std::map<unsigned int, unsigned int>>
        &cellwise_color_predicate_map,
      std::vector<std::function<const Function<spacedim> *(
        const typename Triangulation<dim, spacedim>::cell_iterator &)>>
        &color_enrichments);


    /**
     * 创建一个使用FE_Enriched元素构建的 hp::FECollection
     * 对象，该对象本身是使用颜色丰富的函数构建的，其大小等于颜色的数量。
     * @param[in]  fe_sets 一个颜色列表的向量  @param[in]
     * color_enrichments
     * 一个函数的向量，接收单元格并返回一个函数指针。
     * @param[in]  fe_base 基础有限元素  @param[in]  fe_enriched
     * 丰富的有限元素  @param[in]  fe_nothing
     * 自由度为零的有限元素  @param[out]  fe_collection
     * 有限元素的集合
     *
     */
    template <int dim, int spacedim>
    void
    make_fe_collection_from_colored_enrichments(
      const unsigned int n_colors,
      const std::vector<std::set<unsigned int>>
        &fe_sets, // total list of color sets possible
      const std::vector<std::function<const Function<spacedim> *(
        const typename Triangulation<dim, spacedim>::cell_iterator &)>>
        &color_enrichments, // color wise enrichment functions
      const FiniteElement<dim, spacedim> &fe_base, // basic FE element
      const FiniteElement<dim, spacedim>
        &fe_enriched, // FE multiplied by enrichment function
      const FE_Nothing<dim, spacedim> &fe_nothing,
      hp::FECollection<dim, spacedim> &fe_collection);
  }    // namespace internal
#endif // DOXYGEN



  /**
   * ColorEnriched::Helper  类创建了一个FE_Enriched有限元的集合
   * (hp::FECollection)
   * ，与DoFHandler一起用于具有多个可能重叠的子域的单独富集函数的域。请注意，重叠的区域可能有多个与之相关的富集函数。这是用FE_Enriched对象的一般构造函数来实现的，它允许不同的充实功能。
   * 考虑到一个域有多个丰富的子域，这些子域是不相交的，也就是说，彼此之间没有联系。
   * 为了确保 $C^0$
   * 富集子域（由单一富集函数表征）和非富集域之间界面的连续性，我们可以在富集子域中使用FE_Enriched对象，在非富集域中使用标准有限元（例如：FE_Q）包装成FE_Enriched对象（其内部使用主导的FE_Nothing对象）。
   * 更多信息请参考FE_Enriched文件。需要注意的是，一个FE_Enriched对象是用一个基础FE（FiniteElement对象）和一个或多个丰富的FE构建的。FE_Nothing是一个假的丰富的FE。
   * 当两个充实子域共享一个接口时，情况会变得更加复杂。当子域的丰富函数的数量相同时，一个子域的FE_Enriched对象被构造成这样，每个丰富的FE与另一个子域的FE_Enriched对象中的FE_Nothing配对（形象化）。
   * 例如，让FEs
   * fe_enr1和fe_enr2对应两个子域，这两个子域将被用于充实函数。那么这两个子域的FE_Enriched对象就分别用[fe_base,
   * fe_enr1, fe_nothing]和[fe_base, fe_nothing, fe_enr2]来建立。
   * 注意，丰富的FE的向量的大小（在FE_Enriched构造函数中使用）等于2，与丰富函数的数量相同。当充实函数的数量不一样时，额外的充实FE会与FE_Nothing配对。这就保证了接口处的丰富的DOF被
   * DoFTools::make_hanging_node_constraints() 函数设置为零。
   * 使用这两种策略，我们使用一般构造函数构造适当的FE_Enriched。注意，这是在没有悬挂节点的网格上完成的。
   * 现在考虑一个有多个子域的域，这些子域之间可能共享一个接口。如前所述，每个子域的FE_Enriched对象中丰富的FE的数量需要等于子域的数量。这是因为我们没有使用域的连接方式的信息，任何子域都可能与任何其他子域共享接口（暂时不考虑重叠！）。然而，一般来说，一个给定的子域只与几个子域共享一个接口。这就需要使用图形着色算法来减少丰富的FE的向量的大小（在FE_Enriched构造函数中使用）。通过给没有共享接口的子域使用相同的颜色，可以构建一个为每个子域返回不同的充实函数的
   * 'std::function'
   * 。那么丰富的FE的向量的大小就等于用于谓词（或子域）的不同颜色的数量。
*  @note  图形着色函数， SparsityTools::color_sparsity_pattern, 用于为子域分配颜色需要MPI（使用 Utilities::MPI::MPI_InitFinalize 来初始化MPI和必要的Zoltan设置）。  基于Zoltan的着色函数是一种并行着色算法，但通过 SparsityTools::color_sparsity_pattern. 串行使用。 Helper类的构建需要一个基础FiniteElement @p fe_base, 一个丰富的FiniteElement @p fe_enriched （用于所有丰富函数），一个谓词函数向量（用于定义子域）以及相应的丰富函数。FECollection对象，是与DoFHandler对象一起使用的FE_Enriched对象的集合，可以使用成员函数build_fe_collection来检索，该函数还修改了DoFHandler对象的活动FE索引（作为参数提供给build_fe_collection函数）。    <h3>Simple example</h3> 考虑一个由谓词函数定义的具有三个子域的域。  不同的单元与FE指数相关联，如下图所示。可以看到三个大小相等的方形子域'a'、'b'和'c'。与这些子域相关的谓词也被标记为'a'、'b'和'c'。  子域'a'和'b'与标有FE指数3的单元相交。c'中的细胞被标记为FE指数1。可以看出，'a'和'b'、'b'和'c'之间存在连接，但'a'和'c'没有连接。    <style>div.image img[src="3source_fe_indices.png"]{width:25%;}</style>\endhtmlonly  @image html 3source_fe_indices.png "Active FE indices"  如前所述，谓词的颜色是用图着色算法分配的。每个谓词都是图中的一个节点，如果两个子域共享一个接口，那么相应的谓词应该被赋予不同的颜色。  谓词的颜色与图像中显示的不同。图片中的颜色是按FE指数计算的）。)   谓词'a'和'c'可以被赋予相同的颜色，因为它们没有联系，但是赋予'b'的颜色必须与'a'和'c'不同。    索引(i)为 @p fe_collection  (hp::FECollection) 的有限元的名称可以通过 <code>fe_collection[index].get_name()</code> 得到，如下表所示。请注意，所有的FE_Enriched元素都是相同的尺寸，并且像之前讨论的那样使用FE_Nothing<2>(dominating)。      <table>
   * <tr> <th>Active FE index</th> <th>FiniteElement name</th> </tr> <tr>
   * <td>0</td>
   * <td><code>FE_Enriched<2>[FE_Q<2>(2)-FE_Nothing<2>(dominating)-FE_Nothing<2>(dominating)]</code></td>
   * </tr> <tr> <td>1</td>
   * <td><code>FE_Enriched<2>[FE_Q<2>(2)-FE_Q<2>(1)-FE_Nothing<2>(dominating)]</code></td>
   * </tr> <tr> <td>2</td>
   * <td><code>FE_Enriched<2>[FE_Q<2>(2)-FE_Q<2>(1)-FE_Q<2>(1)]</code></td>
   * </tr> <tr> <td>3</td>
   * <td><code>FE_Enriched<2>[FE_Q<2>(2)-FE_Nothing<2>(dominating)-FE_Q<2>(1)]</code></td>
   * </tr>
   * </table>
   * 这个类所使用的内部数据成员在问题被解决时需要是可用的。这可以通过声明该对象为静态对象来保证，只有在程序终止时才会被取消分配。另一种方法是将其作为包含类的数据成员。由于在构造Helper时，谓词和丰富函数的向量可能不可用，所以可以使用一个
   * 'std::shared_ptr'
   * 到Helper对象，并在谓词和丰富函数可用时构造。
   * @warning
   * 目前的实现依赖于为每个单元分配一个材料ID，该ID在设置和h-adaptive细化后不应被修改。对于一个给定的单元，材料ID被用来定义颜色谓词图，它不会随着细化而改变。
   * <h3>Example usage:</h3>
   * @code
   * FE_Q<dim> fe_base(2);
   * FE_Q<dim> fe_enriched(1);
   * std::vector< predicate_function<dim> > predicates;
   * std::vector< std::shared_ptr<Function<dim>> > enrichments;
   *
   * Triangulation<dim> triangulation;
   * DoFHandler<dim>    dof_handler(triangulation);
   *
   * static ColorEnriched::Helper<dim> FE_helper(fe_base,
   *                                           fe_enriched,
   *                                           predicates,
   *                                           enrichments);
   * const hp::FECollection<dim>&
   * fe_collection(FE_helper.build_fe_collection(dof_handler));
   * @endcode
   *
   */
  template <int dim, int spacedim = dim>
  struct Helper
  {
    /**
     * 助手类的构造函数。          @param  fe_base
     * 一个基本的有限元素  @param  fe_enriched
     * 一个丰富的有限元素  @param  谓词  std::vector
     * 定义子域的谓词。      <code>@p predicates[i]</code>
     * 如果一个单元属于索引为(i)的子域，则该单元返回真。
     * @param 充实函数 std::vector 的充实函数
     *
     */
    Helper(const FiniteElement<dim, spacedim> &                    fe_base,
           const FiniteElement<dim, spacedim> &                    fe_enriched,
           const std::vector<predicate_function<dim, spacedim>> &  predicates,
           const std::vector<std::shared_ptr<Function<spacedim>>> &enrichments);

    /**
     * 准备一个DoFHandler对象。网格单元的活动FE指数被初始化，以便与
     * ColorEnriched::Helper<dim,spacedim>::fe_collection.   @param
     * dof_handler一个DoFHandler对象 @return   hp::FECollection,
     * 所需的有限元素集合 @p dof_handler.  一起工作。
     *
     */
    const hp::FECollection<dim, spacedim> &
    build_fe_collection(DoFHandler<dim, spacedim> &dof_handler);

  private:
    /**
     * 包含DoFHandler对象所需的有限元对象的集合。
     *
     */
    hp::FECollection<dim, spacedim> fe_collection;

    /**
     * 用于构建 ColorEnriched::Helper<dim,spacedim>::fe_collection.
     * 所需的FE_Enriched对象的基础FiniteElement。
     *
     */
    const FiniteElement<dim, spacedim> &fe_base;

    /**
     * 用于构建 ColorEnriched::Helper<dim,spacedim>::fe_collection.
     * 所要求的FE_Enriched对象的丰富的FiniteElement。
     *
     */
    const FiniteElement<dim, spacedim> &fe_enriched;

    /**
     * 用于构建 ColorEnriched::Helper<dim,spacedim>::fe_collection
     * 所要求的FE_Enriched对象的零自由度的有限元素。
     *
     */
    const FE_Nothing<dim, spacedim> fe_nothing;

    /**
     * std::vector 定义子域的谓词。      <code>@p predicates[i]</code>
     * 如果一个单元属于索引为(i)的子域，则返回真。
     *
     */
    const std::vector<predicate_function<dim, spacedim>> predicates;

    /**
     * std::vector 对应于谓词的富集函数。在构建
     * ColorEnriched::Helper<dim,spacedim>::fe_collection. 时需要这些。
     *
     */
    const std::vector<std::shared_ptr<Function<spacedim>>> enrichments;

    /**
     * 任何可调用目标的别名模板，如函数、lambda表达式、接受
     * Triangulation<dim,spacedim>::cell_iterator
     * 并返回指向Function<dim>指针的函数对象。这被用来定义
     * Helper<dim,spacedim>::color_enrichments
     * ，它为Triangulation<dim,spacedim>中的单元格返回一个富集函数。
     *
     */
    using cell_iterator_function = std::function<const Function<spacedim> *(
      const typename Triangulation<dim, spacedim>::cell_iterator &)>;

    /**
     * std::vector
     * 的函数，它接收一个单元格并返回一个函数指针。color_enrichments[i](cell_iterator)
     * 返回一个指向该单元格的正确富集函数（即其对应的谓词具有i的颜色）的指针。
     *
     */
    std::vector<cell_iterator_function> color_enrichments;

    /**
     * std::vector
     * 的颜色（无符号int）与每个子域相关。没有两个相连的子域（即共享一个顶点的子域）具有相同的颜色。
     *
     */
    std::vector<unsigned int> predicate_colors;

    /**
     * predicate_colors中不同颜色的总数
     *
     */
    unsigned int n_colors;

    /**
     * 用于关联单元格的地图，该地图反过来将单元格中的活动谓词的颜色与相应的谓词id关联起来。
     *
     */
    std::map<unsigned int, std::map<unsigned int, unsigned int>>
      cellwise_color_predicate_map;

    /**
     * 为一组给定的谓词和DoFHandler对象提供不同的可能颜色集的向量。
     *
     */
    std::vector<std::set<unsigned int>> fe_sets;
  };
} // namespace ColorEnriched

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_fe_enriched_h


