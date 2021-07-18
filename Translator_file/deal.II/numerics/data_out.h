//include/deal.II-translator/numerics/data_out_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#ifndef dealii_data_out_h
#define dealii_data_out_h



#include <deal.II/base/config.h>

#include <deal.II/grid/filtered_iterator.h>

#include <deal.II/numerics/data_out_dof_data.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace DataOutImplementation
  {
    /**
     * 一个派生类，用于DataOut类中。这是一个用于WorkStream上下文的文档中讨论的AdditionalData类的数据结构。
     *
     */
    template <int dim, int spacedim>
    struct ParallelData : public ParallelDataBase<dim, spacedim>
    {
      ParallelData(
        const unsigned int               n_datasets,
        const unsigned int               n_subdivisions,
        const std::vector<unsigned int> &n_postprocessor_outputs,
        const dealii::hp::MappingCollection<dim, spacedim> &mapping,
        const std::vector<
          std::shared_ptr<dealii::hp::FECollection<dim, spacedim>>>
          &                                           finite_elements,
        const UpdateFlags                             update_flags,
        const std::vector<std::vector<unsigned int>> &cell_to_patch_index_map);

      std::vector<Point<spacedim>> patch_evaluation_points;

      const std::vector<std::vector<unsigned int>> *cell_to_patch_index_map;
    };
  } // namespace DataOutImplementation
} // namespace internal



/**
 * 该类是提供由定义在单元格集合上的有限元场描述的数据输出的主要类。
 * 该类是DataOut_DoFData类所提出的功能的实际实现。它提供了一个函数build_patches()来生成要以某种图形格式写入的数据。这个基类的文档中描述了大部分的接口和一个使用实例。
 * 这个类唯一提供的是函数build_patches()，它循环处理由基类的attach_dof_handler()函数存储的三角形的所有单元（不包括不属于当前处理器的
 * parallel::distributed::Triangulation
 * 对象的单元），并将这些上的数据转换为实际的补丁，这些补丁是后来被基类的函数输出的对象。你可以给函数一个参数，决定在每个坐标方向上进行多少次细分，也就是说，每个补丁由多少个子单元组成。默认值是1，但你可能想为高阶元素选择一个更高的数字，例如，二阶元素为2，三阶元素为3，依此类推。(见
 * step-11
 * 中的例子。)这个参数的目的是由于大多数图形程序不允许在文件格式中指定高阶多项式函数。只有顶点的数据可以被绘制出来，然后在单元格内部显示为双线性插值。如果你有更高阶的有限元，这可能是不够的，获得更好的输出的唯一方法是将网格的每个单元细分为几个单元进行图形输出。当然，由于输出格式的限制，你所看到的仍然是输出中每个单元的双线性插值（其中这些单元不是使用中的三角形单元的细分），但至少是更细的网格上的高阶多项式的双线性插值。
 * 请注意，在调用了一次build_patches()之后，你可以调用DataOutInterface的一个或多个write()函数。因此，你可以以一种以上的格式输出相同的数据，而不需要重建补丁。
 *
 *  <h3>User interface information</h3>
 * 这个类的基类，DataOutBase、DataOutInterface和DataOut_DoFData提供了它们自己的几个接口。请参考DataOutBase类的文档，以了解目前支持的不同的输出格式，DataOutInterface是在运行时选择输出时使用哪种格式的方法，当新的格式可用时，不需要调整你的程序，以及决定输出方面的标志。DataOut_DoFData()类的文档中有一个使用节点数据来生成输出的例子。
 *
 *  <h3>Extensions</h3>
 * 默认情况下，这个类为所有活动单元产生补丁。有时候，这并不是你想要的，也许是因为数量太多（而且太小，无法单独看到），或者是因为你只想看到域的某个区域（例如在
 * step-46 中只看到域的流体部分），或者其他一些原因。
 * 为此，在内部build_patches()不会自己生成要转换为补丁的单元格序列，而是依赖于两个私有的
 * std::function
 * 对象first_cell_function()和next_cell_function()。默认情况下，它们分别返回第一个活动单元，和下一个活动单元。但这可以通过set_cell_selection()函数来改变，该函数允许你替换这一行为。set_cell_selection()想知道的是，你想如何挑选出应该产生输出的第一个单元，以及给定一个产生输出的单元，你想如何挑选下一个单元。
 * 例如，这可能只包括位于域的一部分的单元格（例如。如果你不关心其他地方的解决方案，例如在完美匹配层方法中衰减出射波的缓冲区），或者如果你不希望在自适应细化网格的所有层次上产生输出，因为这会产生太多的数据（在这种情况下，由你实现set_cell_selection()的`first_cell`和`next_cell`参数返回的单元集将包括非活动单元，而
 * DataOut::build_patches()  ]
 * 将简单地采取解决方案的内插值，而不是这些单元格上的精确值来输出）。)
 *
 *
 * @ingroup output
 *
 *
 */
template <int dim, int spacedim = dim>
class DataOut : public DataOut_DoFData<dim, dim, spacedim, spacedim>
{
public:
  /**
   * 对所考虑的dof处理程序类的迭代器类型的类型定义。
   *
   */
  using cell_iterator =
    typename DataOut_DoFData<dim, dim, spacedim, spacedim>::cell_iterator;

  /**
   * set_cell_selection()中使用的返回第一个单元格的函数对象的类型。
   *
   */
  using FirstCellFunctionType =
    typename std::function<cell_iterator(const Triangulation<dim, spacedim> &)>;

  /**
   * 在set_cell_selection()中使用的返回下一个单元格的函数对象的类型。
   *
   */
  using NextCellFunctionType =
    typename std::function<cell_iterator(const Triangulation<dim, spacedim> &,
                                         const cell_iterator &)>;

  /**
   * 枚举，描述域中应写入弯曲边界的单元的部分。在现实中，我们知道没有任何文件格式真正支持弯曲的边界，但是如果调用
   * DataOut::build_patches()
   * 时细分的数量大于1，就可以通过将边缘绘制成一连串的直线（在三维中将面绘制成双线性补丁的集合）来模拟这一点。
   *
   */
  enum CurvedCellRegion
  {
    /**
     * 几何体或边界描述将永远不会被查询到弯曲的几何体。这意味着，即使你每个单元有一个以上的细分（具体含义见
     * DataOut::build_patches()
     * ），即使几何体真的是弯曲的，每个单元仍然会被细分，好像它只是一个双线或三线单元。
     *
     */
    no_curved_cells,

    /**
     * 对于位于边界的单元，即至少有一个面位于边界的单元，几何体或边界描述将被查询到弯曲的几何体。如果您没有为单元格内部附加流形描述，而只为边界处的面附加流形描述，这就足够了。
     *
     */
    curved_boundary,

    /**
     * 将对所有单元和所有面进行几何体描述的查询，无论它们是否在边界上。如果你将流形对象附加到单元格（而不是只附加到边界面），这个选项是合适的。
     *
     */
    curved_inner_cells
  };

  /**
   * 构造函数。
   *
   */
  DataOut();

  /**
   * 这是该类的核心功能，因为它建立了基类的低级函数所要写入的补丁列表。从本质上讲，补丁是三角形和DoFHandler对象的每个单元上的数据的一些中间表示，然后可以用来以某种可视化程序可读的格式编写文件。
   * 你可以在这个类的一般文档中找到关于这个函数的使用概述。在这个类的基类DataOut_DoFData的文档中也提供了一个例子。
   * @param  n_subdivisions
   * 一个参数，决定了这个函数将在每个单元中建立多少个
   * "补丁"。如果你在调用时没有指定这个值，或者提供默认值为0，那么这个值将被解释为
   * DataOutInterface::default_subdivisions
   * ，在大多数情况下它将等于1（除非你把它设置为其他值）。这个参数的目的是将网格的每个单元细分为2d的
   * $2\times 2, 3\times 3, \ldots$ "补丁"，以及 $2\times 2\times 2,
   * 3\times 3\times 3, \ldots$
   * （如果传递值为2，3，等等），其中每个补丁代表单元的常规细分的数据，分成相等的部分。大多数时候，这是没有必要的，每个单元格输出一个补丁正是你想要绘制的解决方案。也就是说，我们写进文件用于可视化的数据只能代表每个细胞上的（双、三）线性数据，而大多数可视化程序事实上只能可视化这种数据。如果你用（双、三）线性有限元工作，这已经很好了，在这种情况下，你能看到的正是已经计算出来的数据。另一方面，如果你使用(双、三)二次元，那么写入输出文件的只是一个(双、三)线性插值到当前网格上，也就是说，只有顶点上的数值。如果这还不够好，你可以，例如，指定
   * @p n_subdivisions
   * 等于2，在一次精炼的网格上绘制解决方案，或者如果设置为3，在一个每个单元由3乘3补丁代表的网格上。在这些较小的斑块上，考虑到输出格式的限制，数据仍然是线性内插的，但是在更细的网格上对二次元数据进行线性内插，仍然比在原始网格上更好地表示实际的二次元曲面。
   * 换句话说，使用这个参数不能帮助你准确地绘制出解决方案，但如果你使用更高多项式程度的有限元，它可以让你更接近。
   * @note
   * 当使用高阶有限元时，指定`n_subdivisions>1`是有用的，但一般来说，它实际上不会导致可视化显示高阶多项式曲面
   *
   * 相反，你只是得到一个更细的网格上的高阶曲面的（双，三）线性插值。然而，当通过
   * DataOutInterface::write_vtk() 或 DataOutInterface::write_vtu()
   * （其中DataOutInterface是当前类的基类）以VTK和VTU文件格式输出解决方案时，正如我们在教程中经常做的，你可以通过
   * DataOutBase::VtkFlags 结构提供一组标志，其中包括
   * DataOutBase::VtkFlags::write_higher_order_cells
   * 标志。当设置时，这个函数产生的细分将被解释为高阶多项式的支持点，然后将实际被可视化。例如，这在
   * step-11
   * 中显示。然而，值得注意的是，这需要一个足够新的基于VTK的可视化程序的版本。
   *
   */
  virtual void
  build_patches(const unsigned int n_subdivisions = 0);

  /**
   * 与上述相同，只是额外的第一个参数定义了一个映射，该映射将用于输出的生成。如果<tt>n_subdivisions>1</tt>，则可以使用该映射来计算源于域的边界单元的细分斑块内部的点，即高阶映射导致通过使用更多的细分来表示一个弯曲的边界。一些像MappingQEulerian这样的映射会导致域内部的单元出现弯曲。  如果你将流形描述附加到三角化的单元上，情况也是如此（详见 @ref manifold  "流形"
   * ）。然而，没有简单的方法来查询映射或流形是否真的导致了弯曲的单元。
   * 因此，最后一个参数 @p curved_region
   * 取三个值之一，导致完全没有弯曲的单元，边界上有弯曲的单元（默认）或整个域中有弯曲的单元。关于这三个选项的更多信息，请参阅CurvedCellRegion枚举的描述。
   * 即使对于非曲线单元，映射参数也可以用于欧拉映射（见类MappingQ1Eulerian），其中映射不仅用于确定单元内部的点的位置，而且还用于确定顶点的位置。
   * 它提供了一个机会，可以在实际进行计算的变形三角形上观察解决方案，即使网格在内部是以未变形的配置存储的，而变形只是通过一个额外的向量来跟踪每个顶点的变形情况。
   *
   */
  virtual void
  build_patches(const Mapping<dim, spacedim> &mapping,
                const unsigned int            n_subdivisions = 0,
                const CurvedCellRegion        curved_region  = curved_boundary);

  /**
   * 与上述相同，但对于 hp::MappingCollection. 而言。
   *
   */
  virtual void
  build_patches(const hp::MappingCollection<dim, spacedim> &mapping,
                const unsigned int                          n_subdivisions = 0,
                const CurvedCellRegion curved_region = curved_boundary);

  /**
   * 一个允许选择应该为哪些单元格生成输出的函数。这个函数需要两个参数，都是
   * `std::function`
   * 对象，可用于生成输出的第一个单元格是什么，以及给定一个单元格的下一个函数是什么。通过这些函数对象，可以选择一个应该产生输出的单元格子集（例如，只选择那些属于域的一部分的单元格
   *
   * 例如，在诸如 step-46
   * 这样的代码中，只选择属于域的一部分的单元），或者完全改变输出的位置*（例如，在多网格层次结构的非活动单元上产生输出，或者如果网格的最细层次非常细，产生图形输出会导致数据量过大）。
   * @param[in]  first_cell
   * 一个函数对象，以该类工作的三角形为参数，应该返回应该产生输出的第一个单元。
   * @param[in]  next_cell
   * 一个函数对象，以三角形以及生成输出的最后一个单元为参数，并应返回应生成输出的下一个单元。如果没有下一个单元，即如果`next_cell`函数对象的输入参数是要生成输出的最后一个单元，那么`next_cell`必须返回`triangulation.end()`。
   * 这些函数对象并不难写，但也不是一目了然。因此，这个函数还有一个变体，它接受一个IteratorFilter参数并自己生成相应的函数。
   * @note
   * 如果你有单元数据（与节点或dof数据相反），如错误指示器，那么你必须确保`first_cell`和`next_cell`函数对象只在活动单元上行走，因为单元数据不能被内插到更粗的单元。如果你确实有单元格数据并使用这对函数，而它们返回的是一个非活动单元格，那么将抛出一个异常。
   *
   */
  void
  set_cell_selection(
    const std::function<cell_iterator(const Triangulation<dim, spacedim> &)>
      &                                                        first_cell,
    const std::function<cell_iterator(const Triangulation<dim, spacedim> &,
                                      const cell_iterator &)> &next_cell);

  /**
   * 前一个函数的变体，根据作为参数的FilteredIterator对象中编码的过滤器，选择所有单元格的一个子集进行输出。产生该参数的典型方法是通过make_filtered_iterator()函数。
   * 另外，由于FilteredIterator对象可以从一个谓词（即一个返回
   * "bool
   * "的函数对象）中创建，所以可以只用一个lambda函数来调用这个函数，然后它将自动转换为FilteredIterator对象。例如，下面这段代码就可以工作。
   * @code
   * DataOut<dim> data_out;
   * data_out.set_cell_selection(
   *        [](const typename Triangulation<dim>::cell_iterator &cell) {
   *            return (cell->is_active() && cell->subdomain_id() == 0);
   *        });
   * @endcode
   * 在这种情况下，lambda函数选择了所有那些 @ref GlossActive "活动 "
   * 的、子域id为0的单元格。然后，这些将是唯一产生输出的单元。
   *
   */
  void
  set_cell_selection(const FilteredIterator<cell_iterator> &filtered_iterator);

  /**
   * 返回由set_cell_selection()设置的用于确定第一个和下一个单元格的两个函数对象。
   *
   */
  const std::pair<FirstCellFunctionType, NextCellFunctionType>
  get_cell_selection() const;

private:
  /**
   * 一个函数对象，用于选择第一个单元格是什么，以生成图形输出。参见set_cell_selection()函数以了解更多信息。
   *
   */
  std::function<cell_iterator(const Triangulation<dim, spacedim> &)>
    first_cell_function;

  /**
   * 一个函数对象，用于选择下一个单元格，并在其上生成图形输出，给定一个前一个单元格。参见set_cell_selection()函数以了解更多信息。
   *
   */
  std::function<cell_iterator(const Triangulation<dim, spacedim> &,
                              const cell_iterator &)>
    next_cell_function;

  /**
   * 建立一个补丁。这个函数在WorkStream上下文中被调用。
   * 这里的第一个参数是迭代器，第二个是抓取数据对象。在调用
   * WorkStream::run().
   * 时，所有下面的参数都与特定的值联系在一起。该函数不接受CopyData对象，而是在自己的堆栈中分配一个，以提高内存访问效率。
   *
   */
  void
  build_one_patch(
    const std::pair<cell_iterator, unsigned int> *cell_and_index,
    internal::DataOutImplementation::ParallelData<dim, spacedim> &scratch_data,
    const unsigned int     n_subdivisions,
    const CurvedCellRegion curved_cell_region);
};

namespace Legacy
{
  /**
   * @deprecated  使用没有DoFHandlerType模板的 dealii::DataOut
   * 代替。
   *
   */
  template <int dim, typename DoFHandlerType = DoFHandler<dim>>
  using DataOut DEAL_II_DEPRECATED =
    dealii::DataOut<dim, DoFHandlerType::space_dimension>;
} // namespace Legacy


DEAL_II_NAMESPACE_CLOSE

#endif


