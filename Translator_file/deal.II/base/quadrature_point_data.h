//include/deal.II-translator/base/quadrature_point_data_0.txt
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

#ifndef dealii_quadrature_point_data_h
#define dealii_quadrature_point_data_h

#include <deal.II/base/config.h>

#include <deal.II/base/quadrature.h>
#include <deal.II/base/std_cxx17/optional.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/vector.h>

#include <map>
#include <type_traits>
#include <vector>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup Quadrature */ 
 /*@{*/ 

/**
 * 一个用于在每个由类型为  @p CellIteratorType
 * 的迭代器代表的单元格存储数据向量的类  @p DataType
 * 。这个类的底层结构和初始化()方法是这样设计的：人们可以使用从基础DataType类派生出来的不同子类来存储给定单元的数据。这意味着指针的使用，在我们的例子中
 *
 * -  std::shared_ptr.
 * 类型 @p DataType
 * 是任意的，但是当使用从TransferableQuadraturePointData派生的类时，可以使用
 * parallel::distributed::ContinuousQuadratureDataTransfer. 的设施。
 *
 *
 * @note
 * 每个单元上存储的数据类型可以是不同的。然而，在单元格内，该类存储的是单一数据类型的对象向量。由于这个原因，当例如采用水平集方法来描述材料行为时，这个类可能不够灵活。
 *
 *
 */
template <typename CellIteratorType, typename DataType>
class CellDataStorage : public Subscriptor
{
public:
  /**
   * 默认构造函数。
   *
   */
  CellDataStorage() = default;

  /**
   * 默认的解构器。
   *
   */
  ~CellDataStorage() override = default;

  /**
   * 在  @p cell  上初始化数据，以存储  @p
   * number_of_data_points_per_cell  类型为  @p T  的对象。  类型名
   * @p T  可能是另一个类，它是由基  @p DataType
   * 类派生的。为了初始化对象的向量，我们必须假设 @p T
   * 类有一个默认的构造函数。
   * 这个函数必须在每个要存储数据的单元格中调用。
   * 数据被初始化后，可以使用get_data()函数对其进行修改。
   * @note  此函数的后续调用与 @p cell
   * 相同，不会改变与之相关的对象。为了删除存储的数据，请使用erase()函数。
   * @note  可以为不同的单元使用不同的类型 @p T
   * ，这可能反映出，例如，在领域的不同部分，连续体力学的不同构成模型。
   * @note
   * 该方法第一次被调用时，会存储一个SmartPointer到拥有该单元的Triangulation对象。以后调用此方法时，希望该单元来自于同一存储的三角剖面。
   * @pre  类型 @p T 需要等于 @p DataType, 或是从 @p DataType.
   * 派生的类  @p T 需要是可默认构造的。
   *
   */
  template <typename T = DataType>
  void
  initialize(const CellIteratorType &cell,
             const unsigned int      number_of_data_points_per_cell);

  /**
   * 同上，但对于所有本地拥有的单元格，即`cell->is_locally_owned()==true`的迭代器范围，从
   * @p cell_start 开始，直到但不包括 @p cell_end  。
   *
   */
  template <typename T = DataType>
  void
  initialize(const CellIteratorType &cell_start,
             const CellIteratorType &cell_end,
             const unsigned int      number_of_data_points_per_cell);

  /**
   * 移除存储在 @p cell. 的数据
   * 如果数据被移除，则返回true。  如果没有数据附在 @p
   * cell, 上，这个函数将不做任何事情，并返回false。
   * @note
   * 此函数还将检查是否有对存储在此单元格上的数据的未决引用。也就是说，对存储的数据的唯一引用是由这个类作出的。
   *
   */
  bool
  erase(const CellIteratorType &cell);

  /**
   * 清除存储在此对象中的所有数据。
   *
   */
  void
  clear();

  /**
   * 获取位于  @p cell  的数据的向量。
   * 一个可能的附加类型名 @p T
   * 是基类DataType可以被投到的类。由于 @p DataType
   * 是以共享指针的形式存储的，所以通过值而不是引用返回向量的开销最小。
   * 如果 @p T 类与 @p DataType
   * 在逐个单元的基础上不相同，这允许灵活性。      @pre
   * 类型 @p T 需要与提供给initialize()的类匹配。      @pre   @p
   * cell
   * 必须来自用于初始化()单元数据的同一个三角结构。
   *
   */
  template <typename T = DataType>
  std::vector<std::shared_ptr<T>>
  get_data(const CellIteratorType &cell);

  /**
   * 获取位于  @p cell  的数据的常数指针的向量。
   * 一个可能的附加类型名 @p T
   * 是基类DataType可以被投向的类。由于 @p DataType
   * 是以共享指针的形式存储的，所以通过值而不是引用来返回一个向量的开销最小。
   * 如果 @p T 类与 @p DataType
   * 在逐个单元的基础上不相同，这允许灵活性。      @pre
   * 类型 @p T 需要与提供给initialize()的类匹配。      @pre   @p
   * cell
   * 必须来自用于初始化()单元数据的同一个三角结构。
   *
   */
  template <typename T = DataType>
  std::vector<std::shared_ptr<const T>>
  get_data(const CellIteratorType &cell) const;

  /**
   * 返回一个 std_cxx17::optional ，表明 @p cell
   * 是否包含相关数据。如果数据是可用的，解除引用
   * std_cxx17::optional
   * 会显示一个指向正交点的基础数据的指针向量。
   * 一个可能的附加类型名 @p T
   * 是基类DataType可以被投射到的类。由于 @p DataType
   * 是以共享指针的形式存储的，所以通过值而不是引用来返回矢量的开销最小。
   * 如果 @p T 类与 @p DataType
   * 在逐个单元的基础上不相同，这允许灵活性。      @pre
   * 类型 @p T 需要与提供给initialize()的类匹配。    @pre   @p
   * cell
   * 必须来自用于初始化()单元数据的同一个三角结构。
   *
   */
  template <typename T = DataType>
  std_cxx17::optional<std::vector<std::shared_ptr<T>>>
  try_get_data(const CellIteratorType &cell);

  /**
   * 返回一个 std_cxx17::optional ，表明 @p cell
   * 是否包含相关数据。如果数据可用，解除引用
   * std_cxx17::optional
   * 会显示一个指向正交点的基础数据的常数指针向量。
   * 一个可能的附加类型名 @p T
   * 是基类DataType可以被投射到的类。由于 @p DataType
   * 是以共享指针的形式存储的，所以通过值而不是引用来返回矢量的开销最小。
   * 如果 @p T 类与 @p DataType
   * 在逐个单元的基础上不相同，这允许灵活性。      @pre
   * 类型 @p T 需要与提供给initialize()的类匹配。    @pre   @p
   * cell
   * 必须来自用于初始化()单元数据的同一个三角结构。
   *
   */
  template <typename T = DataType>
  std_cxx17::optional<std::vector<std::shared_ptr<const T>>>
  try_get_data(const CellIteratorType &cell) const;

private:
  /**
   * 尺寸的数量
   *
   */
  static constexpr unsigned int dimension =
    CellIteratorType::AccessorType::dimension;

  /**
   * 空间维度的数量
   *
   */
  static constexpr unsigned int space_dimension =
    CellIteratorType::AccessorType::space_dimension;

  /**
   * 为了确保CellDataStorage中的所有单元都来自同一个三角形，我们需要在该类中存储一个对该三角形的引用。
   *
   */
  SmartPointer<const Triangulation<dimension, space_dimension>,
               CellDataStorage<CellIteratorType, DataType>>
    tria;

  /**
   * 一个用于存储每个单元格的数据向量的地图。
   * 我们需要使用CellId作为键，因为它在自适应细化过程中保持唯一。
   *
   */
  std::map<CellId, std::vector<std::shared_ptr<DataType>>> map;

  /**
   * @addtogroup  Exceptions
   *
   */
  DeclExceptionMsg(
    ExcCellDataTypeMismatch,
    "Cell data is being retrieved with a type which is different than the type used to initialize it");

  /**
   * @addtogroup  Exceptions
   *
   */
  DeclExceptionMsg(
    ExcTriangulationMismatch,
    "The provided cell iterator does not belong to the triangulation that corresponds to the CellDataStorage object.");
};


/**
 * 一个抽象的类，它规定了在细化或重新划分过程中，单个正交点的数据可以转移的要求。
 * 这个类提供了一个框架，代表正交点数据的派生类可以声明他们存储了多少标量值，然后实现将这些标量打包和解压到数组的函数。这些数组用于将数据从一个单元的正交点转移到另一个单元的正交点，以及在网格细化和重新划分时在处理器之间转移。父单元和子单元之间的正交点数据传输需要某种投影和/或内插。一个可能的实现是通过
 * parallel::distributed::ContinuousQuadratureDataTransfer
 * 类中实现的二级投影和延长矩阵。
 * 要存储和访问从该类派生的类的实例，请参阅CellDataStorage类。
 *
 *
 */
class TransferableQuadraturePointData
{
public:
  /**
   * 默认构造函数。
   *
   */
  TransferableQuadraturePointData() = default;

  /**
   * 默认的虚拟解构器。
   *
   */
  virtual ~TransferableQuadraturePointData() = default;

  /**
   * 返回将从用户的DataType类中打包/解包的数值总数。因此，它也是pack_values()和unpack_values()中向量的大小。
   *
   */
  virtual unsigned int
  number_of_values() const = 0;

  /**
   * 一个必须在派生类中实现的虚拟函数，用于将派生类中存储的所有数据打包成一个向量
   * @p values  。
   * 这个向量将包含每个正交点的所有标量和/或量纲数据。
   * @note  该函数将以大小为number_of_values()的 @p values
   * 被调用。
   * 实现时可能仍有一个断言来检查是否确实如此。
   *
   */
  virtual void
  pack_values(std::vector<double> &values) const = 0;

  /**
   * 与上面的情况相反，即把一个向量 @p values
   * 解压为存储在这个类中的数据。
   * @note 该函数将以大小为number_of_values()的 @p values
   * 被调用。
   * 实现可能仍然有一个断言来检查是否确实如此。
   *
   */
  virtual void
  unpack_values(const std::vector<double> &values) = 0;
};


#ifdef DEAL_II_WITH_P4EST
namespace parallel
{
  namespace distributed
  {
    /**
     * 在执行  parallel::distributed::Triangulation  的h-adaptive
     * refinement时，用于传输存储在正交点的连续数据的类。
     * <h3>Implementation details</h3>
     * 该类实现了在使用L2投影进行自适应细化的情况下，单元格之间正交点数据的传输。这也包括不同处理器之间的信息自动运输。
     * 为此，该类的构造函数提供了三个主要对象：标量FiniteElement
     * @p projection_fe,   @p mass_quadrature  和  @p data_quadrature
     * 正交规则。    首先，位于每个单元的 @p data_quadrature
     * 的数据被L2投影到由单个FiniteElement  @p projection_fe
     * 定义的连续空间。    这是用
     * FETools::compute_projection_from_quadrature_points_matrix().
     * 实现的。在这样做时，需要这个元素的质量矩阵，这将用
     * @p mass_quadrature
     * 规则计算。如果该单元现在属于另一个处理器，那么数据就会被发送到这个处理器。该类利用了p4est（和
     * parallel::distributed::Triangulation)
     * 的一个特性，允许人们在网格细化和再平衡过程中给单元附加信息。在接收到目标单元的信息后，使用
     * FETools::compute_interpolation_to_quadrature_points_matrix()
     * 计算的矩阵将数据投射回正交点。
     * 在进行局部细化的情况下，该类首先将父元素的局部DoF值投影到每个子元素。
     * 这个类是由 @p DataType 类型模板化的，然而用户的 @p
     * DataType
     * 类必须派生自TransferableQuadraturePointData类。在实践中，这相当于为一个有2个标量的正交点数据实现以下三个函数。
     * @code
     * class MyQData : public TransferableQuadraturePointData
     * {
     * public:
     * double elasticity_parameter_lambda;
     * double elasticity_parameter_mu;
     * unsigned int number_of_values() const
     * {
     *    return 2;
     * }
     *
     * // a function to pack scalars into a vector
     * void pack_values(std::vector<double> &values) const
     * {
     *   Assert (values.size()==2, ExcInternalError());
     *   values[0] = elasticity_parameter_lambda;
     *   values[1] = elasticity_parameter_mu;
     * }
     *
     * void unpack_values(const std::vector<double> &values)
     * {
     *   Assert (values.size() ==2, ExcInternalError());
     *   elasticity_parameter_lambda = values[0];
     *   elasticity_parameter_mu     = values[1];
     * }
     * };
     * @endcode
     * 注意，打包和解包的顺序必须相同。
     * 然后，这个类可以通过以下方式与CellDataStorage一起使用。
     * @code
     * CellDataStorage<typename Triangulation<dim,dim>::cell_iterator,MyQData>
     * data_storage;
     * parallel::distributed::ContinuousQuadratureDataTransfer<dim,MyQData>
     * data_transfer(FE_Q<dim>(2),QGauss<dim>(3),QGauss<dim>(4));
     * //...populate data for all active cells in data_storage
     * //...mark cells for refinement...
     * data_transfer.prepare_for_coarsening_and_refinement(triangulation,data_storage);
     * triangulation.execute_coarsening_and_refinement();
     * //...initialize quadrature point data on new cells by calling
     * // CellDataStorage::initialize()
     * data_transfer.interpolate();
     * @endcode
     * 这种方法可以扩展到具有任意顺序的张量对象的正交点数据，尽管在MyQData类中对数据的打包和解包要多做一点工作。
     * @note  目前不支持粗化。
     * @note  这个类所提供的功能也可以用
     * parallel::distributed::SolutionTransfer.
     * 来实现，但是，这需要以下步骤。(i)
     * 创建一个具有（非连续Galerkin）FiniteElement的辅助DoFHandler，该FiniteElement有足够的分量来表示存储在正交点的所有数据；(ii)
     * 将数据投影到FiniteElement空间，从而将结果存储在全局向量中；(iii)
     * 使用 parallel::distributed::SolutionTransfer
     * 将FE向量投影到新网格上；(iv)
     * 最后通过FEValues类将数据投影回新网格的正交点。ContinuousQuadratureDataTransfer类的目的是简化整个过程，只要求正交点数据类派生自TransferableQuadraturePointData。其他一切都将自动完成。
     * @note
     * 这个类不太适合存储在正交点的值代表不连续场的样本的情况。这种情况的一个例子是，存储在正交点的数据代表了材料的弹性或塑性状态，即在固体中不连续变化的属性。在这种情况下，试图将数据从正交点转移到连续的有限元场（至少在一个单元内）可能会产生过冲和下冲，一旦在不同的正交点集（在子单元或父单元上）进行评估，就会产生没有多大意义的数值。
     *
     */
    template <int dim, typename DataType>
    class ContinuousQuadratureDataTransfer
    {
    public:
      static_assert(
        std::is_base_of<TransferableQuadraturePointData, DataType>::value,
        "User's DataType class should be derived from TransferableQuadraturePointData");

      /**
       * 一个单元格的别名。
       *
       */
      using CellIteratorType =
        typename parallel::distributed::Triangulation<dim>::cell_iterator;

      /**
       * 构造函数，它接收有限元素 @p projection_fe
       * ，用于积分其局部质量矩阵的正交规则 @p
       * mass_quadrature ，最后用于存储 @p DataType.   @pre   @p
       * projection_fe 的正交规则必须是标度值。
       * @note 由于该类在单元格基础上进行投影， @p
       * projection_fe 只需要在单元格内是连续的。
       *
       */
      ContinuousQuadratureDataTransfer(const FiniteElement<dim> &projection_fe,
                                       const Quadrature<dim> &mass_quadrature,
                                       const Quadrature<dim> &data_quadrature);

      /**
       * 为粗化和细化三角形做准备  @p tria  。        @p
       * data_storage
       * 代表应该转移的单元数据，它应该对每个本地拥有的活动单元进行初始化。
       * @note
       * 尽管CellDataStorage类允许在不同的单元上存储来自基类的不同对象，这里我们假设
       * @p data_storage
       * 包含相同类型的对象，更确切地说，它们打包/解包相同的数据。
       *
       */
      void
      prepare_for_coarsening_and_refinement(
        parallel::distributed::Triangulation<dim> &  tria,
        CellDataStorage<CellIteratorType, DataType> &data_storage);

      /**
       * 在网格被细化或粗化之前，将先前存储在该对象中的数据插值到当前活动单元组的正交点上。
       * @note 在调用此函数之前，用户应使用
       * CellDataStorage::initialize().
       * 在新的单元格中填充存储在提供给prepare_for_coarsening_and_refinement()的
       * @p data_storage
       * 对象中的数据，如果不是这样，在调试模式下将会抛出一个异常。
       *
       */
      void
      interpolate();

    private:
      /**
       * 一个回调函数，用于将当前网格上的数据打包成对象，以后可以在细化、粗化和重新划分后进行检索。
       *
       */
      std::vector<char>
      pack_function(
        const typename parallel::distributed::Triangulation<dim>::cell_iterator
          &cell,
        const typename parallel::distributed::Triangulation<dim>::CellStatus
          status);

      /**
       * 一个回调函数，用来解开当前网格上的数据，这些数据在细化、粗化和重新划分之前，已经被打包在网格上了。
       *
       */
      void
      unpack_function(
        const typename parallel::distributed::Triangulation<dim>::cell_iterator
          &cell,
        const typename parallel::distributed::Triangulation<dim>::CellStatus
          status,
        const boost::iterator_range<std::vector<char>::const_iterator>
          &data_range);

      /**
       * 用于投影正交点的数据的FiniteElement。
       *
       */
      const std::unique_ptr<const FiniteElement<dim>> projection_fe;

      /**
       * 将被发送的数据的大小，这取决于数据类型类。
       *
       */
      std::size_t data_size_in_bytes;

      /**
       * 存储DataType的正交点的数量。
       *
       */
      const unsigned int n_q_points;

      /**
       * 从正交点到单个标量的本地DoF的投影矩阵。
       *
       */
      FullMatrix<double> project_to_fe_matrix;

      /**
       * 单个标量从局部DoF到正交点的投影矩阵。
       *
       */
      FullMatrix<double> project_to_qp_matrix;

      /**
       * 辅助矩阵，表示存储在正交点（第二个索引）的每个内部值对
       * @p projection_fe （第一个索引）的本地DoF的投影。
       *
       */
      FullMatrix<double> matrix_dofs;

      /**
       * 在自适应细化的情况下， @p matrix_dofs
       * 对每个子单元的投影。
       *
       */
      FullMatrix<double> matrix_dofs_child;

      /**
       * 辅助矩阵，代表存储在每个正交点（第一索引）的数据（第二索引）。
       *
       */
      FullMatrix<double> matrix_quadrature;

      /**
       * parallel::distributed::Triangulation
       * 在注册pack_callback函数时分配给这个对象的句柄。
       *
       */
      unsigned int handle;

      /**
       * 一个指向CellDataStorage类的指针，其数据将被转移。
       *
       */
      CellDataStorage<CellIteratorType, DataType> *data_storage;

      /**
       * 一个指向分布式三角图的指针，单元格数据被附加到该三角图上。
       *
       */
      parallel::distributed::Triangulation<dim> *triangulation;
    };

  } // namespace distributed

} // namespace parallel

#endif

 /*@}*/ 

#ifndef DOXYGEN

// -------------------  inline and template functions ----------------

//--------------------------------------------------------------------
//                         CellDataStorage
//--------------------------------------------------------------------

template <typename CellIteratorType, typename DataType>
template <typename T>
inline void
CellDataStorage<CellIteratorType, DataType>::initialize(
  const CellIteratorType &cell,
  const unsigned int      n_q_points)
{
  static_assert(std::is_base_of<DataType, T>::value,
                "User's T class should be derived from user's DataType class");
  // The first time this method is called, it has to initialize the reference
  // to the triangulation object
  if (!tria)
    tria = &cell->get_triangulation();
  Assert(&cell->get_triangulation() == tria, ExcTriangulationMismatch());

  const auto key = cell->id();
  if (map.find(key) == map.end())
    {
      map[key] = std::vector<std::shared_ptr<DataType>>(n_q_points);
      // we need to initialize one-by-one as the std::vector<>(q, T())
      // will end with a single same T object stored in each element of the
      // vector:
      const auto it = map.find(key);
      for (unsigned int q = 0; q < n_q_points; q++)
        it->second[q] = std::make_shared<T>();
    }
}



template <typename CellIteratorType, typename DataType>
template <typename T>
inline void
CellDataStorage<CellIteratorType, DataType>::initialize(
  const CellIteratorType &cell_start,
  const CellIteratorType &cell_end,
  const unsigned int      number)
{
  for (CellIteratorType it = cell_start; it != cell_end; it++)
    if (it->is_locally_owned())
      initialize<T>(it, number);
}



template <typename CellIteratorType, typename DataType>
inline bool
CellDataStorage<CellIteratorType, DataType>::erase(const CellIteratorType &cell)
{
  const auto key = cell->id();
  const auto it  = map.find(key);
  if (it == map.end())
    return false;
  Assert(&cell->get_triangulation() == tria, ExcTriangulationMismatch());
  for (unsigned int i = 0; i < it->second.size(); i++)
    {
      Assert(
        it->second[i].unique(),
        ExcMessage(
          "Can not erase the cell data multiple objects reference its data."));
    }

  return (map.erase(key) == 1);
}



template <typename CellIteratorType, typename DataType>
inline void
CellDataStorage<CellIteratorType, DataType>::clear()
{
  // Do not call
  // map.clear();
  // as we want to be sure no one uses the stored objects. Loop manually:
  auto it = map.begin();
  while (it != map.end())
    {
      // loop over all objects and see if no one is using them
      for (unsigned int i = 0; i < it->second.size(); i++)
        {
          Assert(
            it->second[i].unique(),
            ExcMessage(
              "Can not erase the cell data, multiple objects reference it."));
        }
      it = map.erase(it);
    }
}



template <typename CellIteratorType, typename DataType>
template <typename T>
inline std::vector<std::shared_ptr<T>>
CellDataStorage<CellIteratorType, DataType>::get_data(
  const CellIteratorType &cell)
{
  static_assert(std::is_base_of<DataType, T>::value,
                "User's T class should be derived from user's DataType class");
  Assert(&cell->get_triangulation() == tria, ExcTriangulationMismatch());

  const auto it = map.find(cell->id());
  Assert(it != map.end(), ExcMessage("Could not find data for the cell"));

  // It would be nice to have a specialized version of this function for
  // T==DataType. However explicit (i.e full) specialization of a member
  // template is only allowed when the enclosing class is also explicitly (i.e
  // fully) specialized. Thus, stick with copying of shared pointers even when
  // the T==DataType:
  std::vector<std::shared_ptr<T>> res(it->second.size());
  for (unsigned int q = 0; q < res.size(); q++)
    {
      res[q] = std::dynamic_pointer_cast<T>(it->second[q]);
      Assert(res[q], ExcCellDataTypeMismatch());
    }
  return res;
}



template <typename CellIteratorType, typename DataType>
template <typename T>
inline std::vector<std::shared_ptr<const T>>
CellDataStorage<CellIteratorType, DataType>::get_data(
  const CellIteratorType &cell) const
{
  static_assert(std::is_base_of<DataType, T>::value,
                "User's T class should be derived from user's DataType class");
  Assert(&cell->get_triangulation() == tria, ExcTriangulationMismatch());

  const auto it = map.find(cell->id());
  Assert(it != map.end(), ExcMessage("Could not find QP data for the cell"));

  // Cast base class to the desired class. This has to be done irrespectively of
  // T==DataType as we need to return shared_ptr<const T> to make sure the user
  // does not modify the content of QP objects
  std::vector<std::shared_ptr<const T>> res(it->second.size());
  for (unsigned int q = 0; q < res.size(); q++)
    {
      res[q] = std::dynamic_pointer_cast<const T>(it->second[q]);
      Assert(res[q], ExcCellDataTypeMismatch());
    }
  return res;
}

template <typename CellIteratorType, typename DataType>
template <typename T>
inline std_cxx17::optional<std::vector<std::shared_ptr<T>>>
CellDataStorage<CellIteratorType, DataType>::try_get_data(
  const CellIteratorType &cell)
{
  static_assert(std::is_base_of<DataType, T>::value,
                "User's T class should be derived from user's DataType class");
  Assert(&cell->get_triangulation() == tria, ExcTriangulationMismatch());

  const auto it = map.find(cell->id());
  if (it != map.end())
    {
      // Cast base class to the desired class. This has to be done
      // irrespectively of T==DataType as we need to return
      // shared_ptr<const T> to make sure the user
      // does not modify the content of QP objects
      std::vector<std::shared_ptr<T>> result(it->second.size());
      for (unsigned int q = 0; q < result.size(); q++)
        {
          result[q] = std::dynamic_pointer_cast<T>(it->second[q]);
          Assert(result[q], ExcCellDataTypeMismatch());
        }
      return {result};
    }
  else
    {
      return {};
    }
}

template <typename CellIteratorType, typename DataType>
template <typename T>
inline std_cxx17::optional<std::vector<std::shared_ptr<const T>>>
CellDataStorage<CellIteratorType, DataType>::try_get_data(
  const CellIteratorType &cell) const
{
  static_assert(std::is_base_of<DataType, T>::value,
                "User's T class should be derived from user's DataType class");
  Assert(&cell->get_triangulation() == tria, ExcTriangulationMismatch());

  const auto it = map.find(cell->id());
  if (it != map.end())
    {
      // Cast base class to the desired class. This has to be done
      // irrespectively of T==DataType as we need to return
      // shared_ptr<const T> to make sure the user
      // does not modify the content of QP objects
      std::vector<std::shared_ptr<const T>> result(it->second.size());
      for (unsigned int q = 0; q < result.size(); q++)
        {
          result[q] = std::dynamic_pointer_cast<const T>(it->second[q]);
          Assert(result[q], ExcCellDataTypeMismatch());
        }
      return {result};
    }
  else
    {
      return {};
    }
}

//--------------------------------------------------------------------
//                    ContinuousQuadratureDataTransfer
//--------------------------------------------------------------------


/*将使用 @p data_storage 存储在 @p cell 中的每个正交点的 @p matrix_data. 类型的单元数据打包到 @p matrix_data 这里 @p matrix_data 是一个矩阵，其第一个索引对应于单元上的不同正交点，而第二个索引代表存储在DataType类中每个正交点的不同值。

* 
*
*/
template <typename CellIteratorType, typename DataType>
inline void
pack_cell_data(const CellIteratorType &                           cell,
               const CellDataStorage<CellIteratorType, DataType> *data_storage,
               FullMatrix<double> &                               matrix_data)
{
  static_assert(
    std::is_base_of<TransferableQuadraturePointData, DataType>::value,
    "User's DataType class should be derived from QPData");

  if (const auto qpd = data_storage->try_get_data(cell))
    {
      const unsigned int m = qpd->size();
      Assert(m > 0, ExcInternalError());
      const unsigned int n = (*qpd)[0]->number_of_values();
      matrix_data.reinit(m, n);

      std::vector<double> single_qp_data(n);
      for (unsigned int q = 0; q < m; ++q)
        {
          (*qpd)[q]->pack_values(single_qp_data);
          AssertDimension(single_qp_data.size(), n);

          for (unsigned int i = 0; i < n; ++i)
            matrix_data(q, i) = single_qp_data[i];
        }
    }
  else
    {
      matrix_data.reinit({0, 0});
    }
}



/* 与上面的打包函数相反。

* 
*
*/
template <typename CellIteratorType, typename DataType>
inline void
unpack_to_cell_data(const CellIteratorType &                     cell,
                    const FullMatrix<double> &                   values_at_qp,
                    CellDataStorage<CellIteratorType, DataType> *data_storage)
{
  static_assert(
    std::is_base_of<TransferableQuadraturePointData, DataType>::value,
    "User's DataType class should be derived from QPData");

  if (const auto qpd = data_storage->try_get_data(cell))
    {
      const unsigned int n = values_at_qp.n();
      AssertDimension((*qpd)[0]->number_of_values(), n);

      std::vector<double> single_qp_data(n);
      AssertDimension(qpd->size(), values_at_qp.m());

      for (unsigned int q = 0; q < qpd->size(); ++q)
        {
          for (unsigned int i = 0; i < n; ++i)
            single_qp_data[i] = values_at_qp(q, i);
          (*qpd)[q]->unpack_values(single_qp_data);
        }
    }
}


#  ifdef DEAL_II_WITH_P4EST

namespace parallel
{
  namespace distributed
  {
    template <int dim, typename DataType>
    inline ContinuousQuadratureDataTransfer<dim, DataType>::
      ContinuousQuadratureDataTransfer(const FiniteElement<dim> &projection_fe_,
                                       const Quadrature<dim> &   lhs_quadrature,
                                       const Quadrature<dim> &   rhs_quadrature)
      : projection_fe(
          std::unique_ptr<const FiniteElement<dim>>(projection_fe_.clone()))
      , data_size_in_bytes(0)
      , n_q_points(rhs_quadrature.size())
      , project_to_fe_matrix(projection_fe->n_dofs_per_cell(), n_q_points)
      , project_to_qp_matrix(n_q_points, projection_fe->n_dofs_per_cell())
      , handle(numbers::invalid_unsigned_int)
      , data_storage(nullptr)
      , triangulation(nullptr)
    {
      Assert(
        projection_fe->n_components() == 1,
        ExcMessage(
          "ContinuousQuadratureDataTransfer requires scalar FiniteElement"));

      FETools::compute_projection_from_quadrature_points_matrix(
        *projection_fe.get(),
        lhs_quadrature,
        rhs_quadrature,
        project_to_fe_matrix);

      FETools::compute_interpolation_to_quadrature_points_matrix(
        *projection_fe.get(), rhs_quadrature, project_to_qp_matrix);
    }



    template <int dim, typename DataType>
    inline void
    ContinuousQuadratureDataTransfer<dim, DataType>::
      prepare_for_coarsening_and_refinement(
        parallel::distributed::Triangulation<dim> &  tr_,
        CellDataStorage<CellIteratorType, DataType> &data_storage_)
    {
      Assert(data_storage == nullptr,
             ExcMessage("This function can be called only once"));
      triangulation = &tr_;
      data_storage  = &data_storage_;

      handle = triangulation->register_data_attach(
        [this](
          const typename parallel::distributed::Triangulation<
            dim>::cell_iterator &cell,
          const typename parallel::distributed::Triangulation<dim>::CellStatus
            status) { return this->pack_function(cell, status); },
         /*returns_variable_size_data=*/ true);
    }



    template <int dim, typename DataType>
    inline void
    ContinuousQuadratureDataTransfer<dim, DataType>::interpolate()
    {
      triangulation->notify_ready_to_unpack(
        handle,
        [this](
          const typename parallel::distributed::Triangulation<
            dim>::cell_iterator &cell,
          const typename parallel::distributed::Triangulation<dim>::CellStatus
            status,
          const boost::iterator_range<std::vector<char>::const_iterator>
            &data_range) { this->unpack_function(cell, status, data_range); });

      // invalidate the pointers
      data_storage  = nullptr;
      triangulation = nullptr;
    }



    template <int dim, typename DataType>
    inline std::vector<char>
    ContinuousQuadratureDataTransfer<dim, DataType>::pack_function(
      const typename parallel::distributed::Triangulation<dim>::cell_iterator
        &cell,
      const typename parallel::distributed::Triangulation<
        dim>::CellStatus  /*status*/ )
    {
      pack_cell_data(cell, data_storage, matrix_quadrature);

      // project to FE
      const unsigned int number_of_values = matrix_quadrature.n();
      matrix_dofs.reinit(project_to_fe_matrix.m(), number_of_values);
      if (number_of_values > 0)
        project_to_fe_matrix.mmult(matrix_dofs, matrix_quadrature);

      return Utilities::pack(matrix_dofs,  /*allow_compression=*/ false);
    }



    template <int dim, typename DataType>
    inline void
    ContinuousQuadratureDataTransfer<dim, DataType>::unpack_function(
      const typename parallel::distributed::Triangulation<dim>::cell_iterator
        &cell,
      const typename parallel::distributed::Triangulation<dim>::CellStatus
        status,
      const boost::iterator_range<std::vector<char>::const_iterator>
        &data_range)
    {
      Assert((status !=
              parallel::distributed::Triangulation<dim, dim>::CELL_COARSEN),
             ExcNotImplemented());
      (void)status;

      matrix_dofs =
        Utilities::unpack<FullMatrix<double>>(data_range.begin(),
                                              data_range.end(),
                                               /*allow_compression=*/ false);
      const unsigned int number_of_values = matrix_dofs.n();
      if (number_of_values == 0)
        return;

      matrix_quadrature.reinit(n_q_points, number_of_values);

      if (cell->has_children())
        {
          // we need to first use prolongation matrix to get dofvalues on child
          // cells based on dofvalues stored in the parent's data_store
          matrix_dofs_child.reinit(projection_fe->n_dofs_per_cell(),
                                   number_of_values);
          for (unsigned int child = 0; child < cell->n_children(); ++child)
            if (cell->child(child)->is_locally_owned())
              {
                projection_fe
                  ->get_prolongation_matrix(child, cell->refinement_case())
                  .mmult(matrix_dofs_child, matrix_dofs);

                // now we do the usual business of evaluating FE on quadrature
                // points:
                project_to_qp_matrix.mmult(matrix_quadrature,
                                           matrix_dofs_child);

                // finally, put back into the map:
                unpack_to_cell_data(cell->child(child),
                                    matrix_quadrature,
                                    data_storage);
              }
        }
      else
        {
          // if there are no children, evaluate FE field at
          // rhs_quadrature points.
          project_to_qp_matrix.mmult(matrix_quadrature, matrix_dofs);

          // finally, put back into the map:
          unpack_to_cell_data(cell, matrix_quadrature, data_storage);
        }
    }

  } // namespace distributed

} // namespace parallel

#  endif // DEAL_II_WITH_P4EST

#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif


