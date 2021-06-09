//include/deal.II-translator/numerics/fe_field_function_0.txt
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

#ifndef dealii_fe_function_h
#define dealii_fe_function_h

#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/std_cxx17/optional.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/thread_local_storage.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_tools_cache.h>

#include <deal.II/lac/vector.h>


DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
namespace VectorTools
{
  class ExcPointNotAvailableHere;
}
#endif

namespace Functions
{
  /**
   * 这是一个针对给定的dof处理程序和给定的解决方案向量的插值函数。这个函数可以被评估的点必须是在dof处理程序的域内，但除了这一点，没有给出其他要求。这个函数相当慢，因为它需要为你想评估有限元函数的点（或点组）构造一个正交对象。为了做到这一点，它需要找出这些点的位置。
   * 如果你事先知道你的点位于哪个单元，你可以在要求函数的值或梯度之前调用set_active_cell()来加快事情的进展。如果你不这样做，而你的点不在当前存储的单元格中，那么就会调用函数
   * GridTools::find_active_cell_around_point 来找出点的位置。
   * 你可以指定一个可选的映射，以便在网格中寻找点时使用。如果你不这样做，这个函数就会使用一个Q1的映射。
   * 一旦 FEFieldFunction
   * 知道了点的位置，它就会为这些点创建一个正交公式，然后用给定的正交点调用
   * FEValues::get_function_values  或  FEValues::get_function_gradients  。
   * 如果你只需要正交点，而不需要有限元函数的值（你可能想要这个用于邻接插值），你也可以单独使用函数compute_point_locations()。
   * 下面是使用这个函数的一个例子。
   * @code
   *
   * // Generate two triangulations
   * Triangulation<dim> tria_1;
   * Triangulation<dim> tria_2;
   *
   * // Read the triangulations from files, or build them up, or get them
   * // from some place. Assume that tria_2 isentirely* included in tria_1.
   *
   * // Associate a dof handler and a solution to the first triangulation
   * DoFHandler<dim> dh1 (tria_1);
   * Vector<double> solution_1;
   *
   * // On this first domain, set up the various data structures,
   * // assemble matrices, solve the linear system, and get a Nobel
   * // prize for the work we have done here:
   * [...]
   *
   * // Then create a DoFHandler and solution vector for the second domain:
   * DoFHandler<dim> dh2 (tria_2);
   * Vector<double> solution_2;
   *
   * // Finally, project the solution on the first domain onto the
   * // second domain, assuming that this does not require querying
   * // values from outside the first domain:
   * Functions::FEFieldFunction<dim> fe_function_1 (dh_1, solution_1);
   * VectorTools::project (dh_2, constraints_2, quad,
   *                     fe_function_1, solution_2);
   *
   * // Alternatively, we could have also interpolated it:
   * Vector<double> solution_3;
   * VectorTools::interpolate (dh_2, fe_function_1, solution_3);
   * @endcode
   * 上面的代码片断将工作，假设第二个三角形完全包含在第一个三角形中。    FEFieldFunction被设计成一种简单的方式来获取你在不同的、可能不匹配的网格中的计算结果。在这个类中没有假定点的位置的知识，这使得它的工作完全依赖于 GridTools::find_active_cell_around_point 工具。然而，该类可以通过使用 FEFieldFunction::set_active_cell 的方法，对将要计算的点的位置进行 "有根据的猜测"
   * ，所以如果你有一个聪明的方法来告诉你的点在哪里，你可以让这个类知道，从而节省大量的计算时间。
   * <h3>Using FEFieldFunction with
   * parallel::distributed::Triangulation</h3>当使用这个类与并行分布式三角测量对象并在一个特定的点上评估解决方案时，不是每个处理器都会拥有评估解决方案的单元。相反，可能发现这个点的单元实际上是一个幽灵或人工单元（见
   * @ref GlossArtificialCell 和 @ref GlossGhostCell ）。
   * 解决方案可以在幽灵单元上进行评估，但对于人工单元，我们无法访问那里的解决方案，在这样的点上评估解决方案的函数将触发一个
   * VectorTools::ExcPointNotAvailableHere. 类型的异常。
   * 为了处理这种情况，你将希望使用如下代码，例如，在原点评估解决方案时（这里使用一个平行的TrilinosWrappers向量来保存解决方案）。
   * @code
   * Functions::FEFieldFunction<dim,TrilinosWrappers::MPI::Vector>
   *   solution_function (dof_handler, solution);
   * Point<dim> origin = Point<dim>();
   *
   * double solution_at_origin;
   * bool   point_found = true;
   * try
   *   {
   *     solution_at_origin = solution_function.value (origin);
   *   }
   * catch (const VectorTools::ExcPointNotAvailableHere &)
   *   {
   *     point_found = false;
   *   }
   *
   * if (point_found == true)
   *   ...do something...;
   * @endcode
   *
   * @ingroup functions
   *
   */
  template <int dim, typename VectorType = Vector<double>, int spacedim = dim>
  class FEFieldFunction : public Function<dim, typename VectorType::value_type>
  {
  public:
    /**
     * 构建一个向量函数。一个智能指针被存储到dof处理程序中，所以你必须确保它在这个对象的整个生命周期内是有意义的。这个函数的组件数量等于有限元对象的组件数量。如果指定了一个映射，这就是用来找出点的位置。否则就使用标准的Q1映射。
     *
     */
    FEFieldFunction(
      const DoFHandler<dim, spacedim> &dh,
      const VectorType &               data_vector,
      const Mapping<dim> &             mapping = StaticMappingQ1<dim>::mapping);

    /**
     * 设置当前单元。如果你事先知道你的点在哪里，你可以通过调用这个函数告诉这个对象。这将使事情变得更快一些。
     *
     */
    void
    set_active_cell(
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &newcell);

    /**
     * 在给定的点上获得一个矢量值。使用单点的效率很低。如果你一次需要多于一个，请使用vector_value_list()函数。出于效率的考虑，如果所有的点都位于同一个单元格上，那就更好了。这并不是强制性的，但是它确实能加快事情的进展。
     * @note
     * 当在 parallel::distributed::Triangulation
     * 上使用这个函数时，当你试图在一个位于人工单元上的点上评估解决方案时，你可能会得到一个异常（见
     * @ref GlossLocallyOwnedCell  ）。
     * 更多信息请参见该类的一般文档中的章节。
     *
     */
    virtual void
    vector_value(
      const Point<dim> &                       p,
      Vector<typename VectorType::value_type> &values) const override;

    /**
     * 返回函数在给定点的值。除非只有一个分量（即该函数是标量的），否则你应该说明你想要评估的分量；它默认为零，即第一个分量。使用单点的效率很低。如果你一次需要多于一个，请使用vector_value_list()函数。出于效率的考虑，如果所有的点都位于同一个单元格上，那就更好了。这并不是强制性的，但是它确实能加快事情的进展。
     * @note
     * 当在一个 parallel::distributed::Triangulation
     * 上使用这个函数时，当你试图在一个位于人工单元上的点上评估解决方案时，你可能会得到一个异常（见
     * @ref GlossLocallyOwnedCell  ）。
     * 更多信息请参见该类的一般文档中的章节。
     *
     */
    virtual typename VectorType::value_type
    value(const Point<dim> &p, const unsigned int component = 0) const override;

    /**
     * 将 @p values 设置为 @p points. 处的函数指定分量的点值
     * 假设 @p values
     * 已经有了合适的大小，即与点阵的大小相同。如果所有的点都位于同一个单元格上，这就相当有效。如果不是这样的话，事情可能会变得有点慢。
     * @note
     * 当在一个 parallel::distributed::Triangulation
     * 上使用这个函数时，你可能会在试图评估位于一个人造单元上的点的解决方案时得到一个异常（见
     * @ref GlossLocallyOwnedCell  ）。
     * 更多信息请参见该类的一般文档中的章节。
     *
     */
    virtual void
    value_list(const std::vector<Point<dim>> &               points,
               std::vector<typename VectorType::value_type> &values,
               const unsigned int component = 0) const override;


    /**
     * 将 @p values 设置为 @p points. 处的函数的点值 假设 @p
     * values
     * 已经有了合适的大小，即与点阵列的大小相同。如果所有的点都位于同一个单元格上，这就相当有效。如果不是这样的话，事情可能会变得有点慢。
     * @note
     * 当在一个 parallel::distributed::Triangulation
     * 上使用这个函数时，当你试图在一个位于人工单元上的点上评估解决方案时，你可能会得到一个异常（见
     * @ref GlossLocallyOwnedCell  ）。
     * 更多信息请参见该类的一般文档中的章节。
     *
     */
    virtual void
    vector_value_list(const std::vector<Point<dim>> &points,
                      std::vector<Vector<typename VectorType::value_type>>
                        &values) const override;

    /**
     * 返回函数在给定点的所有分量的梯度。
     * 使用单点的效率很低。如果你一次需要多于一个，请使用vector_value_list()函数。出于效率的考虑，如果所有的点都位于同一个单元格上，那就更好了。这并不是强制性的，但是它确实能加快事情的进展。
     * @note
     * 当在一个 parallel::distributed::Triangulation
     * 上使用这个函数时，当你试图在一个位于人工单元上的点上评估解决方案时，你可能会得到一个异常（见
     * @ref GlossLocallyOwnedCell ）。
     * 更多信息请参见该类的一般文档中的章节。
     *
     */
    virtual void
    vector_gradient(const Point<dim> &p,
                    std::vector<Tensor<1, dim, typename VectorType::value_type>>
                      &gradients) const override;

    /**
     * 返回函数的指定分量在给定点的梯度。使用单点的效率很低。如果你一次需要多于一个，请使用vector_value_list()函数。出于效率的考虑，如果所有的点都位于同一个单元格上，那就更好了。这并不是强制性的，但是它确实能加快事情的进展。
     * @note
     * 当在一个 parallel::distributed::Triangulation
     * 上使用这个函数时，当你试图在一个位于人工单元上的点上评估解决方案时，你可能会得到一个异常（见
     * @ref GlossLocallyOwnedCell ）。
     * 更多信息请参见该类的一般文档中的章节。
     *
     */
    virtual Tensor<1, dim, typename VectorType::value_type>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;

    /**
     * 返回所有给定点的函数的所有分量的梯度。如果所有的点都在同一个单元格上，这就相当有效。如果不是这种情况，事情可能会变得有点慢。
     * @note
     * 当在 parallel::distributed::Triangulation
     * 上使用这个函数时，当你试图在一个位于人工单元上的点上求解时，你可能会得到一个异常（见
     * @ref GlossLocallyOwnedCell  ）。
     * 更多信息请参见该类的一般文档中的章节。
     *
     */
    virtual void
    vector_gradient_list(
      const std::vector<Point<dim>> &p,
      std::vector<std::vector<Tensor<1, dim, typename VectorType::value_type>>>
        &gradients) const override;

    /**
     * 返回函数的指定分量在所有给定点的梯度。
     * 如果所有的点都在同一个单元格上，这就相当有效。如果不是这种情况，事情可能会慢一点。
     * @note
     * 当在 parallel::distributed::Triangulation
     * 上使用这个函数时，当你试图在一个位于人工单元上的点上求解时可能会得到一个异常（见
     * @ref GlossLocallyOwnedCell ）。
     * 更多信息请参见该类的一般文档中的章节。
     *
     */
    virtual void
    gradient_list(
      const std::vector<Point<dim>> &                               p,
      std::vector<Tensor<1, dim, typename VectorType::value_type>> &gradients,
      const unsigned int component = 0) const override;


    /**
     * 计算点<tt>p</tt>处给定分量的拉普拉斯。
     * @note
     * 当在一个 parallel::distributed::Triangulation
     * 上使用这个函数时，当你试图在一个位于人工单元上的点上评估解决方案时，你可能会得到一个异常（见
     * @ref GlossLocallyOwnedCell  ）。
     * 更多信息请参见该类的一般文档中的章节。
     *
     */
    virtual typename VectorType::value_type
    laplacian(const Point<dim> & p,
              const unsigned int component = 0) const override;

    /**
     * 计算点<tt>p</tt>上所有分量的拉普拉斯，并将其存储在<tt>values</tt>中。
     * @note 当在一个
     * parallel::distributed::Triangulation
     * 上使用这个函数时，当你试图在一个位于人工单元上的点上评估解决方案时，你可能会得到一个异常（见
     * @ref GlossLocallyOwnedCell  ）。
     * 更多信息请参见该类的一般文档中的章节。
     *
     */
    virtual void
    vector_laplacian(
      const Point<dim> &                       p,
      Vector<typename VectorType::value_type> &values) const override;

    /**
     * 计算一组点上的一个分量的拉普拉斯（Laplacian）。
     * @note 当在一个
     * parallel::distributed::Triangulation
     * 上使用这个函数时，当你试图在一个位于人工单元上的点上求解时，可能会得到一个异常（见
     * @ref GlossLocallyOwnedCell ）。
     * 更多信息请参见该类的一般文档中的章节。
     *
     */
    virtual void
    laplacian_list(const std::vector<Point<dim>> &               points,
                   std::vector<typename VectorType::value_type> &values,
                   const unsigned int component = 0) const override;

    /**
     * 计算一组点上的所有分量的拉普拉斯系数。
     * @note 当在一个
     * parallel::distributed::Triangulation
     * 上使用这个函数时，当你试图在一个位于人工单元上的点上求解时，可能会得到一个异常（见
     * @ref GlossLocallyOwnedCell ）。
     * 更多信息请参见该类的一般文档中的章节。
     *
     */
    virtual void
    vector_laplacian_list(const std::vector<Point<dim>> &points,
                          std::vector<Vector<typename VectorType::value_type>>
                            &values) const override;

    /**
     * 给出一组位于域中的点（或者，在并行三角计算的情况下，位于域的局部拥有的部分或当前处理器的幽灵单元上），将这些点分类到至少有一个点所在的每个单元的桶中。
     * 这个函数填充了三个输出向量。  @p cells,   @p qpoints  和
     * @p maps.
     * 第一个是包含这些点的单元格的列表，第二个是与第一个列表的每个单元格相匹配的正交点的列表，第三个包含给定的正交点的索引，即
     * @p points[maps[3][4]]
     * 最后是第四个单元格的第五个正交点。          @return
     * 这个函数返回共同包含给定 @p
     * 点的单元格的数量。这也等同于输出数组的长度。
     * 这个函数简单地调用了 GridTools::compute_point_locations
     * ：使用原始函数避免了在每次函数调用时计算一个新的Cache。
     *
     */
    unsigned int
    compute_point_locations(
      const std::vector<Point<dim>> &points,
      std::vector<typename DoFHandler<dim, spacedim>::active_cell_iterator>
        &                                     cells,
      std::vector<std::vector<Point<dim>>> &  qpoints,
      std::vector<std::vector<unsigned int>> &maps) const;

  private:
    /**
     * 持有本地cell_hint的类型定义。
     *
     */
    using cell_hint_t = Threads::ThreadLocalStorage<
      typename DoFHandler<dim, spacedim>::active_cell_iterator>;

    /**
     * 指向dof处理程序的指针。
     *
     */
    SmartPointer<const DoFHandler<dim, spacedim>,
                 FEFieldFunction<dim, VectorType, spacedim>>
      dh;

    /**
     * 对实际数据向量的引用。
     *
     */
    const VectorType &data_vector;

    /**
     * 一个对正在使用的映射的引用。
     *
     */
    const Mapping<dim> &mapping;

    /**
     * 缓存对象
     *
     */
    GridTools::Cache<dim, spacedim> cache;

    /**
     * 最新的单元格提示。
     *
     */
    mutable cell_hint_t cell_hint;

    /**
     * 给定一个单元格，如果它确实位于该单元格内，则返回该单元格内给定点的参考坐标。否则返回一个未初始化的
     * std_cxx17::optional 对象。
     *
     */
    std_cxx17::optional<Point<dim>>
    get_reference_coordinates(
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
      const Point<dim> &point) const;
  };
} // namespace Functions

namespace Legacy
{
  namespace Functions
  {
    /**
     * @deprecated  使用没有DoFHandlerType模板的
     * dealii::Functions::FEFieldFunction 代替。
     *
     */
    template <int dim,
              typename DoFHandlerType = DoFHandler<dim>,
              typename VectorType     = Vector<double>>
    using FEFieldFunction DEAL_II_DEPRECATED = dealii::Functions::
      FEFieldFunction<dim, VectorType, DoFHandlerType::space_dimension>;
  } // namespace Functions
} // namespace Legacy


DEAL_II_NAMESPACE_CLOSE

#endif


