//include/deal.II-translator/distributed/cell_weights_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2021 by the deal.II authors
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

#ifndef dealii_distributed_cell_weights_h
#define dealii_distributed_cell_weights_h

#include <deal.II/base/config.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_handler.h>


DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  /**
   * 任何时候 parallel::TriangulationBase
   * 被重新分区，不管是根据要求还是通过细化/粗化，单元格将被分配到所有子域中，以达到平等的平衡工作量。如果每个单元的工作量不同，对于具有hp-capabilities的DoFHandler对象来说通常是这样的，我们可以通过为不同的单元引入单独的权重来考虑到这一点。
   * 这个类允许通过查询与DoFHandler的每个单元相关联的FiniteElement来计算这些用于负载平衡的权重。人们可以从该类提供的预定义权重算法中选择，也可以提供一个自定义的算法。
   * 如果相关的DoFHandler还没有被初始化，即它的
   * hp::FECollection 是空的，所有单元的权重将被评估为零。
   * 该类提供了两种不同的方式来连接所选择的加权函数和链接的
   * parallel::TriangulationBase.
   * 的相应信号。推荐的方式包括创建一个该类的对象，该对象将在创建时自动负责注册加权函数，并在销毁时取消其注册。为了达到令人满意的工作平衡结果，每个与我们工作的三角形相关的DoFHandler都需要存在这个类的一个对象。
   * 连接的加权函数可以使用 CellWeights::reinit()
   * 函数随时改变。下面的代码片断演示了如何实现每个单元按其当前自由度数加权。我们选择了一个`1000`的系数，对应于每个单元在创建时被分配到的初始权重。
   * @code
   * parallel::CellWeights<dim, spacedim> cell_weights(
   * hp_dof_handler,
   * parallel::CellWeights<dim, spacedim>::ndofs_weighting({1000, 1}));
   * @endcode
   * 另一方面，你也能够通过使用这个类的静态成员函数来手动处理信号连接。在这种情况下，一个类似的代码例子看起来如下。
   * @code
   * boost::signals2::connection connection =
   * hp_dof_handler.get_triangulation().signals.cell_weight.connect(
   *   parallel::CellWeights<dim, spacedim>::make_weighting_callback(
   *     hp_dof_handler,
   *     parallel::CellWeights<dim, spacedim>::ndofs_weighting(
   *       {1000, 1}));
   * @endcode
   * 这个类的使用在  step-75  中演示。
   * @note  参见 Triangulation::Signals::cell_weight
   * ，了解更多关于加权和负载平衡的信息。
   * @note
   * 注意这个类在这个类的构造函数中把权重函数连接到Triangulation。如果与DoFHandler相关的Triangulation在后者的生命周期内通过
   * DoFHandler::reinit(),
   * 发生变化，将在weight_callback()函数中触发断言。使用
   * CellWeights::reinit()
   * 取消对旧三角函数的注册，并将其连接到新的三角函数。
   * @ingroup distributed
   *
   */
  template <int dim, int spacedim = dim>
  class CellWeights
  {
  public:
    /**
     * 一个别名，它定义了一个函数的特性，可以在负载平衡期间用于加权单元。
     * 这种加权函数的参数是一个单元的迭代器和重新分区后将分配给它的未来有限元。
     * 它们返回一个无符号整数，被解释为单元的权重，或者说，与之相关的额外计算负荷。
     *
     */
    using WeightingFunction = std::function<
      unsigned int(const typename DoFHandler<dim, spacedim>::cell_iterator &,
                   const FiniteElement<dim, spacedim> &)>;

    /**
     * 构造函数。          @param[in]  dof_handler
     * 将被用于确定每个单元的有限元的DoFHandler。
     * @param[in]  weighting_function
     * 在负载平衡期间确定每个单元重量的函数。
     *
     */
    CellWeights(const dealii::DoFHandler<dim, spacedim> &dof_handler,
                const WeightingFunction &                weighting_function);

    /**
     * 解构器。        断开之前连接到权重信号的函数。
     *
     */
    ~CellWeights();

    /**
     * 将不同的 @p weighting_function 连接到与 @p dof_handler.
     * 相关的三角图上 断开先前连接到加权信号的函数。
     *
     */
    void
    reinit(const DoFHandler<dim, spacedim> &dof_handler,
           const WeightingFunction &        weighting_function);

    /**
     * 将 @p weighting_function
     * 转换为符合回调函数条件的不同类型，可以连接到Triangulation的加权信号。
     * 这个函数确实<b>not</b>将转换后的函数连接到与 @p
     * dof_handler. 相关的Triangulation。
     *
     */
    static std::function<unsigned int(
      const typename dealii::Triangulation<dim, spacedim>::cell_iterator &cell,
      const typename dealii::Triangulation<dim, spacedim>::CellStatus status)>
    make_weighting_callback(const DoFHandler<dim, spacedim> &dof_handler,
                            const WeightingFunction &weighting_function);

    /**
     * @name  选择加权函数 
     * @{ 
     *
     */

    /**
     * 在每个单元格上选择一个常数权重 @p factor 。
     *
     */
    static WeightingFunction
    constant_weighting(const unsigned int factor = 1000);

    /**
     * 通过 @p coefficients 提供的一对浮点数 $(a,b)$
     * 以下列方式确定每个具有 $n_K$ 自由度的单元 $K$
     * 的权重。\f[ w_K = a \, n_K^b
     * \f]由于要求细胞权重为整数，所以右边将被四舍五入为最接近的整数。
     *
     */
    static WeightingFunction
    ndofs_weighting(const std::pair<float, float> &coefficients);

    /**
     * 容器 @p coefficients 提供了成对的浮点数 $(a_i, b_i)$
     * ，它们以下列方式确定每个具有 $n_K$ 自由度的单元 $K$
     * 的权重 $w_K$ 。\f[ w_K = \sum_i a_i \, n_K^{b_i}
     * \f]由于要求单元格权重为整数，所以右手边将被四舍五入为最近的整数。
     *
     */
    static WeightingFunction
    ndofs_weighting(const std::vector<std::pair<float, float>> &coefficients);

    /**
     * @}
     *
     */

  private:
    /**
     * 与连接到DoFHandler的Triangulation的相应cell_weight信号的连接。
     *
     */
    boost::signals2::connection connection;

    /**
     * 一个回调函数，它将连接到 @p triangulation,
     * 的cell_weight信号， @p dof_handler
     * 被连接到该信号。最终返回每个单元的重量，由作为参数提供的
     * @p weighting_function 决定。如果 @p dof_handler
     * 还没有被初始化，则返回0。
     *
     */
    static unsigned int
    weighting_callback(
      const typename dealii::Triangulation<dim, spacedim>::cell_iterator &cell,
      const typename dealii::Triangulation<dim, spacedim>::CellStatus status,
      const DoFHandler<dim, spacedim> &                 dof_handler,
      const parallel::TriangulationBase<dim, spacedim> &triangulation,
      const WeightingFunction &                         weighting_function);
  };
} // namespace parallel


DEAL_II_NAMESPACE_CLOSE

#endif


