//include/deal.II-translator/meshworker/vector_selector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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

#ifndef dealii_mesh_worker_vector_selector_h
#define dealii_mesh_worker_vector_selector_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/any_data.h>
#include <deal.II/algorithms/named_selection.h>

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/tensor.h>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int, int>
class FEValuesBase;
#endif

namespace MeshWorker
{
  /**
   * 一个从命名向量列表中选择向量的类。
   * 由于AnyData对象中的向量数量可能随着应用程序或循环的嵌套而增加，因此能够选择那些实际用于计算残差等的向量非常重要。这个类组织了这种选择。
   * 例如，在IntegrationWorker中，它被用来确定哪些值、导数或二次导数被实际计算。
   * @ingroup MeshWorker
   *
   */
  class VectorSelector : public Subscriptor
  {
  public:
    /**
     * 在有限元函数的选择中添加一个矢量。参数是矢量的名称和指标，即要从矢量中提取哪些信息。名称指的是AnyData对象中的一个条目，它将被initialize()所识别。
     * 三个bool参数表明，是否要在每个单元或面计算有限元函数的值、梯度和Hessians。
     *
     */
    void
    add(const std::string &name,
        const bool         values    = true,
        const bool         gradients = false,
        const bool         hessians  = false);

    /**
     * 与上面的函数相同，但有可能选择全局矢量的一个块。
     *
     */
    //      void add(const std::string& name,
    //               const unsigned int selected_block,
    //             bool values = true,
    //             bool gradients = false,
    //             bool hessians = false);

    /**
     * 用一个数据向量初始化选择区域。add()只输入向量的名称，这些名称将在单元格和面的积分循环中使用，而这个函数将这些名称与AnyData对象中的实际向量联系起来。
     * @note
     * 该函数缓存了与名称相关的索引。因此，在AnyData对象改变后，每次都必须调用它。
     *
     */
    void
    initialize(const AnyData &);

    /**
     * 检查是否有任何矢量被选中。
     *
     */
    bool
    empty() const;

    /**
     * 如果任何向量的值被选中，则返回true。
     *
     */
    bool
    has_values() const;

    /**
     * 如果任何向量的梯度被选中，则返回true。
     *
     */
    bool
    has_gradients() const;

    /**
     * 如果为任何矢量选择了 hessians，则返回 true。
     *
     */
    bool
    has_hessians() const;

    /**
     * 数值的向量数
     *
     */
    unsigned int
    n_values() const;

    /**
     * 梯度的向量数
     *
     */
    unsigned int
    n_gradients() const;

    /**
     * Hessians的向量数
     *
     */
    unsigned int
    n_hessians() const;

    /**
     * 第1个值的向量索引
     *
     */
    unsigned int
    value_index(const unsigned int i) const;

    /**
     * 第1个梯度的向量索引
     *
     */
    unsigned int
    gradient_index(const unsigned int i) const;

    /**
     * 第一位Hessian的向量索引
     *
     */
    unsigned int
    hessian_index(const unsigned int i) const;

    /**
     * 打印选择的内容到流中。
     *
     */
    template <class StreamType, typename DATA>
    void
    print(StreamType &s, const AnyData &v) const;

    /**
     * 打印选择的数量到流中。
     *
     */
    template <class StreamType>
    void
    print(StreamType &s) const;

    /**
     * 这个对象所使用的内存。
     *
     */
    std::size_t
    memory_consumption() const;

  protected:
    /**
     * 用来计算数值的向量的选择。
     *
     */
    NamedSelection value_selection;

    /**
     * 选择用于计算梯度的向量。
     *
     */
    NamedSelection gradient_selection;

    /**
     * 选择用于计算斜率的向量。
     *
     */
    NamedSelection hessian_selection;
  };

  /**
   * 基于VectorSelector，这是IntegrationInfo用来计算正交点的源向量值的类。
   * @ingroup MeshWorker
   *
   */
  template <int dim, int spacedim = dim, typename Number = double>
  class VectorDataBase : public VectorSelector
  {
  public:
    /**
     * 构造函数
     *
     */
    VectorDataBase() = default;

    /**
     * 来自基类对象的构造函数
     *
     */
    VectorDataBase(const VectorSelector &);

    /**
     * 用一个AnyData对象初始化，并在VectorSelector基类中缓存索引。
     * @note
     * 在调用此函数之前，确保VectorSelector基类被填充了合理的数据。
     *
     */
    void
    initialize(const AnyData &);

    /**
     * 虚拟的，但空的析构器。
     *
     */
    virtual ~VectorDataBase() override = default;

    /**
     * 添加到VectorSelector的唯一函数是一个抽象的虚拟函数，在派生类模板中实现，由IntegrationInfo调用。
     * 根据我们基类中的选择，它用有限元函数的局部数据填充前三个参数。它通常是为整个FES系统，或者为每个基元单独调用。
     * @param  values 是填充有正交点的有限元函数值的向量。
     * @param  梯度是填充有限元函数在正交点的导数的向量。
     * @param
     * hessians是填充有限元函数在正交点的二次导数的向量。
     * @param
     * fe是FEValuesBase对象，用于计算函数值。它的UpdateFlags必须被适当地设置。
     * @param  index是本地索引向量。如果 @p fe
     * 指的是系统的基本元素，这个向量应该按块排序，下面的参数
     * @p start 和 @p size 指定使用 @p indices 的子集。
     * @param  分量是 @p values,  @p gradients 和 @p
     * 在此函数中输入的第一个索引。          @param
     * n_comp是要填充的组件的数量。          @param
     * start是这个块在 @p indices,
     * 中的第一个索引，如果没有使用基元，则为零。
     * @param
     * size是当前元素或基础元素的每个单元的道夫数量。
     *
     */
    virtual void
    fill(std::vector<std::vector<std::vector<Number>>> &values,
         std::vector<std::vector<std::vector<Tensor<1, spacedim, Number>>>>
           &gradients,
         std::vector<std::vector<std::vector<Tensor<2, spacedim, Number>>>>
           &                                         hessians,
         const FEValuesBase<dim, spacedim> &         fe,
         const std::vector<types::global_dof_index> &index,
         const unsigned int                          component,
         const unsigned int                          n_comp,
         const unsigned int                          start,
         const unsigned int                          size) const;

    /**
     * 从水平向量填充本地数据向量。执行与其他fill()相同的操作，但使用单元格层次来访问层次向量中的一个层次，而不是活动单元格上的全局数据向量。
     *
     */
    virtual void
    mg_fill(std::vector<std::vector<std::vector<Number>>> &values,
            std::vector<std::vector<std::vector<Tensor<1, spacedim, Number>>>>
              &gradients,
            std::vector<std::vector<std::vector<Tensor<2, spacedim, Number>>>>
              &                                         hessians,
            const FEValuesBase<dim, spacedim> &         fe,
            const unsigned int                          level,
            const std::vector<types::global_dof_index> &index,
            const unsigned int                          component,
            const unsigned int                          n_comp,
            const unsigned int                          start,
            const unsigned int                          size) const;

  protected:
    AnyData data;
  };


  /**
   * 基于VectorSelector，这是为某种类型的向量实现函数
   * VectorDataBase::fill() 的类，使用AnyData按名称识别向量。
   * @ingroup MeshWorker
   *
   */
  template <typename VectorType, int dim, int spacedim = dim>
  class VectorData
    : public VectorDataBase<dim, spacedim, typename VectorType::value_type>
  {
  public:
    /**
     * 构造函数。
     *
     */
    VectorData() = default;

    /**
     * 使用预填充VectorSelector的构造函数
     *
     */
    VectorData(const VectorSelector &);

    /**
     * 用一个命名向量的对象进行初始化。
     *
     */
    void
    initialize(const AnyData &);

    /**
     * 用一个单一的向量初始化，并在VectorSelector基类中缓存索引。
     * @note
     * 在调用此函数之前，确保VectorSelector基类被填充了合理的数据。
     *
     */
    void
    initialize(const VectorType *, const std::string &name);

    virtual void
    fill(std::vector<std::vector<std::vector<typename VectorType::value_type>>>
           &values,
         std::vector<std::vector<
           std::vector<Tensor<1, spacedim, typename VectorType::value_type>>>>
           &gradients,
         std::vector<std::vector<
           std::vector<Tensor<2, spacedim, typename VectorType::value_type>>>>
           &                                         hessians,
         const FEValuesBase<dim, spacedim> &         fe,
         const std::vector<types::global_dof_index> &index,
         const unsigned int                          component,
         const unsigned int                          n_comp,
         const unsigned int                          start,
         const unsigned int                          size) const override;

    virtual void
    mg_fill(
      std::vector<std::vector<std::vector<typename VectorType::value_type>>>
        &values,
      std::vector<std::vector<
        std::vector<Tensor<1, spacedim, typename VectorType::value_type>>>>
        &gradients,
      std::vector<std::vector<
        std::vector<Tensor<2, spacedim, typename VectorType::value_type>>>>
        &                                         hessians,
      const FEValuesBase<dim, spacedim> &         fe,
      const unsigned int                          level,
      const std::vector<types::global_dof_index> &index,
      const unsigned int                          component,
      const unsigned int                          n_comp,
      const unsigned int                          start,
      const unsigned int                          size) const override;

    /**
     * 这个对象所使用的内存。
     *
     */
    std::size_t
    memory_consumption() const;
  };


  /**
   * 基于VectorSelector，这是实现某类多级向量的函数
   * VectorDataBase::fill() 的类，使用AnyData按名称识别向量。
   * @ingroup MeshWorker
   *
   */
  template <typename VectorType, int dim, int spacedim = dim>
  class MGVectorData : public VectorData<VectorType, dim, spacedim>
  {
  public:
    /**
     * 构造函数。
     *
     */
    MGVectorData() = default;

    /**
     * 使用预填充VectorSelector的构造函数
     *
     */
    MGVectorData(const VectorSelector &);

    /**
     * 用一个命名向量的对象进行初始化
     *
     */
    void
    initialize(const AnyData &);

    /**
     * 用一个单一的向量初始化，并在VectorSelector基类中缓存索引。
     * @note
     * 在调用此函数之前，确保VectorSelector基类被填充了合理的数据。
     *
     */
    void
    initialize(const MGLevelObject<VectorType> *, const std::string &name);
  };


  //----------------------------------------------------------------------//

  inline void
  VectorSelector::add(const std::string &name,
                      const bool         values,
                      const bool         gradients,
                      const bool         hessians)
  {
    if (values)
      value_selection.add(name);
    if (gradients)
      gradient_selection.add(name);
    if (hessians)
      hessian_selection.add(name);
  }


  // inline void
  // VectorSelector::add(const std::string& name,
  //   const unsigned int block,
  //   bool values, bool gradients, bool hessians)
  //{
  //  if (values) value_selection.add(name, block);
  //  if (gradients) gradient_selection.add(name, block);
  //  if (hessians) hessian_selection.add(name, block);
  //}


  inline void
  VectorSelector::initialize(const AnyData &src)
  {
    value_selection.initialize(src);
    gradient_selection.initialize(src);
    hessian_selection.initialize(src);
  }

  inline bool
  VectorSelector::empty() const
  {
    return (value_selection.size() == 0 && gradient_selection.size() == 0 &&
            hessian_selection.size() == 0);
  }


  inline bool
  VectorSelector::has_values() const
  {
    return value_selection.size() != 0;
  }


  inline bool
  VectorSelector::has_gradients() const
  {
    return gradient_selection.size() != 0;
  }


  inline bool
  VectorSelector::has_hessians() const
  {
    return hessian_selection.size() != 0;
  }


  inline unsigned int
  VectorSelector::n_values() const
  {
    return value_selection.size();
  }


  inline unsigned int
  VectorSelector::n_gradients() const
  {
    return gradient_selection.size();
  }


  inline unsigned int
  VectorSelector::n_hessians() const
  {
    return hessian_selection.size();
  }


  inline unsigned int
  VectorSelector::value_index(const unsigned int i) const
  {
    return value_selection(i);
  }


  inline unsigned int
  VectorSelector::gradient_index(const unsigned int i) const
  {
    return gradient_selection(i);
  }


  inline unsigned int
  VectorSelector::hessian_index(const unsigned int i) const
  {
    return hessian_selection(i);
  }


  template <class StreamType>
  inline void
  VectorSelector::print(StreamType &s) const
  {
    s << "values: " << n_values() << " gradients: " << n_gradients()
      << " hessians: " << n_hessians() << std::endl;
  }


  template <class StreamType, typename DATA>
  inline void
  VectorSelector::print(StreamType &s, const AnyData &v) const
  {
    s << "values:   ";
    for (unsigned int i = 0; i < n_values(); ++i)
      s << " '" << v.name(value_selection(i)) << '\'';
    s << std::endl << "gradients:";
    for (unsigned int i = 0; i < n_gradients(); ++i)
      s << " '" << v.name(gradient_selection(i)) << '\'';
    s << std::endl << "hessians: ";
    for (unsigned int i = 0; i < n_hessians(); ++i)
      s << " '" << v.name(hessian_selection(i)) << '\'';
    s << std::endl;
  }


  inline std::size_t
  VectorSelector::memory_consumption() const
  {
    return sizeof(*this);
  }
} // namespace MeshWorker


DEAL_II_NAMESPACE_CLOSE

#endif


