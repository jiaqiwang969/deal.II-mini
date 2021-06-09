//include/deal.II-translator/matrix_free/cuda_matrix_free_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2021 by the deal.II authors
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


#ifndef dealii_cuda_matrix_free_h
#define dealii_cuda_matrix_free_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_COMPILER_CUDA_AWARE

#  include <deal.II/base/cuda_size.h>
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/quadrature.h>
#  include <deal.II/base/tensor.h>

#  include <deal.II/dofs/dof_handler.h>

#  include <deal.II/fe/fe_update_flags.h>
#  include <deal.II/fe/mapping.h>
#  include <deal.II/fe/mapping_q1.h>

#  include <deal.II/grid/filtered_iterator.h>

#  include <deal.II/lac/affine_constraints.h>
#  include <deal.II/lac/cuda_vector.h>
#  include <deal.II/lac/la_parallel_vector.h>


DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  // forward declaration
#  ifndef DOXYGEN
  namespace internal
  {
    template <int dim, typename Number>
    class ReinitHelper;
  }
#  endif

  /**
   * 这个类收集了所有为矩阵自由实现而存储的数据。这个存储方案是针对用相同数据进行的几个循环而定制的，也就是说，通常在同一个网格上做许多矩阵-向量乘积或残差计算。
   * 该类不实现任何涉及有限元基函数的操作，也就是说，关于在单元上进行的操作。对于这些操作，FEEvaluation类被设计为使用该类中收集的数据。
   * 这个类实现了所有单元的循环（cell_loop()）。这个循环的调度方式是，共享自由度的单元不会同时工作，这意味着可以并行地写入向量，而不必明确地同步访问这些向量和矩阵。这个类没有实现任何形状值，它所做的只是缓存相应的数据。要实现有限元操作，请使用类
   * CUDAWrappers::FEEvalutation.
   * 这个类以不同的顺序遍历单元，而不是通常的deal.II的Triangulation类。
   * @note  只支持float和double。
   * @ingroup CUDAWrappers
   *
   */
  template <int dim, typename Number = double>
  class MatrixFree : public Subscriptor
  {
  public:
    using jacobian_type = Tensor<2, dim, Tensor<1, dim, Number>>;
    using point_type    = Point<dim, Number>;
    using CellFilter =
      FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>;

    /**
     * 使用的并行化方案：parallel_in_elem（自由度层面的并行）或
     * parallel_over_elem（单元格层面的并行）。
     *
     */
    enum ParallelizationScheme
    {
      parallel_in_elem,
      parallel_over_elem
    };

    /**
     * 标准化的数据结构，用于向MatrixFree输送额外的数据。
     *
     */
    struct AdditionalData
    {
      /**
       * 构造函数。
       *
       */
      AdditionalData(
        const ParallelizationScheme parallelization_scheme = parallel_in_elem,
        const UpdateFlags           mapping_update_flags   = update_gradients |
                                                 update_JxW_values |
                                                 update_quadrature_points,
        const bool use_coloring                      = false,
        const bool overlap_communication_computation = false)
        : parallelization_scheme(parallelization_scheme)
        , mapping_update_flags(mapping_update_flags)
        , use_coloring(use_coloring)
        , overlap_communication_computation(overlap_communication_computation)
      {
#  ifndef DEAL_II_MPI_WITH_CUDA_SUPPORT
        AssertThrow(
          overlap_communication_computation == false,
          ExcMessage(
            "Overlapping communication and computation requires CUDA-aware MPI."));
#  endif
        if (overlap_communication_computation == true)
          AssertThrow(
            use_coloring == false || overlap_communication_computation == false,
            ExcMessage(
              "Overlapping communication and coloring are incompatible options. Only one of them can be enabled."));
      }
      /**
       * 使用的并行化方案，在自由度或单元上的并行化。
       *
       */
      ParallelizationScheme parallelization_scheme;
      /**
       * 这个标志是用来决定哪些数量应该被缓存。该类可以缓存梯度计算（反雅各布）、雅各布行列式（JxW）、正交点以及Hessians（雅各布导数）的数据。默认情况下，只有梯度和雅各布行列式乘以正交权重JxW的数据被缓存。如果需要二次导数的正交点，必须由这个字段指定。
       *
       */
      UpdateFlags mapping_update_flags;

      /**
       * 如果为真，使用图形着色。否则，使用原子操作。图形着色保证了位数的可重复性，但是在Pascal和较新的架构上会比较慢。
       *
       */
      bool use_coloring;

      /**
       * 将MPI通信与计算重叠。这需要CUDA感知的MPI和use_coloring必须是假的。
       *
       */
      bool overlap_communication_computation;
    };

    /**
     * 传递给内核的结构。它用于将所有必要的信息从CPU传递给GPU。
     *
     */
    struct Data
    {
      /**
       * 指向正交点的指针。
       *
       */
      point_type *q_points;

      /**
       * 将本地向量中的位置映射到全局向量中的位置。
       *
       */
      types::global_dof_index *local_to_global;

      /**
       * 指向反雅各布的指针。
       *
       */
      Number *inv_jacobian;

      /**
       * 指向Jacobian乘以权重的指针。
       *
       */
      Number *JxW;

      /**
       * 相关MatrixFree对象的ID。
       *
       */
      unsigned int id;

      /**
       * 单元的数量。
       *
       */
      unsigned int n_cells;

      /**
       * 填充物的长度。
       *
       */
      unsigned int padding_length;

      /**
       * 行开始（包括填充物）。
       *
       */
      unsigned int row_start;

      /**
       * 决定在一个给定的单元格上设置约束的掩码。
       *
       */
      unsigned int *constraint_mask;

      /**
       * 如果为真，则使用图形着色，我们可以简单地添加到目的地向量中。否则，使用原子操作。
       *
       */
      bool use_coloring;
    };

    /**
     * 默认构造函数。
     *
     */
    MatrixFree();

    /**
     * 解构器。
     *
     */
    ~MatrixFree();

    /**
     * 返回padding的长度。
     *
     */
    unsigned int
    get_padding_length() const;

    /**
     * 提取在单元格上进行循环所需的信息。DoFHandler和AffineConstraints对象描述了自由度的布局，DoFHandler和映射描述了从单元到实数单元的转换，DoFHandler底层的有限元与正交公式一起描述了局部操作。这个函数需要一个IteratorFilters对象（谓词）来循环处理活动单元的一个子集。当使用MPI时，该谓词应过滤掉非本地拥有的单元。
     *
     */
    template <typename IteratorFiltersType>
    void
    reinit(const Mapping<dim> &             mapping,
           const DoFHandler<dim> &          dof_handler,
           const AffineConstraints<Number> &constraints,
           const Quadrature<1> &            quad,
           const IteratorFiltersType &      iterator_filter,
           const AdditionalData &           additional_data = AdditionalData());

    /**
     * 与上述相同，使用 Iterators::LocallyOwnedCell() 作为谓词。
     *
     */
    void
    reinit(const Mapping<dim> &             mapping,
           const DoFHandler<dim> &          dof_handler,
           const AffineConstraints<Number> &constraints,
           const Quadrature<1> &            quad,
           const AdditionalData &           additional_data = AdditionalData());

    /**
     * 初始化数据结构。与上述相同，但使用Q1映射。
     *
     */
    void
    reinit(const DoFHandler<dim> &          dof_handler,
           const AffineConstraints<Number> &constraints,
           const Quadrature<1> &            quad,
           const AdditionalData &           AdditionalData = AdditionalData());

    /**
     * 返回与 @p color. 相关的数据结构。
     *
     */
    Data
    get_data(unsigned int color) const;

    // clang-format off
    /**
     * 这个方法在所有单元格上运行循环，并在每个元素上并行应用局部操作。
     * @p func 是一个漏斗，它被应用于每个颜色。          @p
     * func 需要定义 /code __device__ void operator()( const unsigned int
     * cell, const typename  CUDAWrappers::MatrixFree<dim,
     * Number>::Datagpu_data,   CUDAWrappers::SharedData<dim,  Number>
     * shared_data, const Number src, Number dst) const; static const unsigned
     * int n_dofs_1d; static const unsigned int n_local_dofs; static const
     * unsigned int n_q_points; \endcode
     *
     */
    // clang-format on
    template <typename Functor, typename VectorType>
    void
    cell_loop(const Functor &   func,
              const VectorType &src,
              VectorType &      dst) const;

    /**
     * 这个方法在所有单元格上运行循环，并在每个元素上并行应用局部操作。这个函数与cell_loop()非常相似，但它使用了一个更简单的函数。
     * @p func  需要定义 /code __device__ void operator()( const unsigned
     * int cell, const typename  CUDAWrappers::MatrixFree<dim,
     * Number>::Datagpu_data);  static const unsigned int n_dofs_1d; static
     * const unsigned int n_local_dofs; static const unsigned int n_q_points;
     * \endcode
     *
     */
    template <typename Functor>
    void
    evaluate_coefficients(Functor func) const;

    /**
     * 从 @p src 复制约束项的值到 @p dst.
     * 这是用来施加零迪里希特边界条件。
     *
     */
    template <typename VectorType>
    void
    copy_constrained_values(const VectorType &src, VectorType &dst) const;

    /**
     * 将 @p dst 中对应于约束值的条目设置为 @p val.
     * 这个函数的主要目的是将cell_loop()中使用的源向量的约束条目设置为零。
     *
     */
    template <typename VectorType>
    void
    set_constrained_values(const Number val, VectorType &dst) const;

    /**
     * 初始化一个序列向量。其大小与DoFHandler对象中的自由度数相对应。
     *
     */
    void
    initialize_dof_vector(
      LinearAlgebra::CUDAWrappers::Vector<Number> &vec) const;

    /**
     * 初始化一个分布式向量。本地元素对应于本地拥有的自由度，鬼魂元素对应于（额外的）本地相关自由度。
     *
     */
    void
    initialize_dof_vector(
      LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA> &vec) const;

    /**
     * 返回本地拥有的活动单元的彩色图。
     *
     */
    const std::vector<std::vector<CellFilter>> &
    get_colored_graph() const;

    /**
     * 返回代表本地拥有的数据的分区器，以及单元格循环需要访问的幽灵索引。分区器是由各自字段给出的本地拥有的道夫和幽灵道夫构建的。如果你想知道这些对象的具体信息，你可以用各自的访问函数来查询它们。如果你只是想初始化一个（平行）向量，你通常应该更喜欢这种数据结构，因为数据交换信息可以从一个向量重复使用到另一个向量。
     *
     */
    const std::shared_ptr<const Utilities::MPI::Partitioner> &
    get_vector_partitioner() const;

    /**
     * 释放所有分配的内存。
     *
     */
    void
    free();

    /**
     * 返回DoFHandler。
     *
     */
    const DoFHandler<dim> &
    get_dof_handler() const;

    /**
     * 返回这个类的内存消耗的近似值，单位是字节。
     *
     */
    std::size_t
    memory_consumption() const;

  private:
    /**
     * 初始化数据结构。
     *
     */
    template <typename IteratorFiltersType>
    void
    internal_reinit(const Mapping<dim> &             mapping,
                    const DoFHandler<dim> &          dof_handler,
                    const AffineConstraints<Number> &constraints,
                    const Quadrature<1> &            quad,
                    const IteratorFiltersType &      iterator_filter,
                    std::shared_ptr<const MPI_Comm>  comm,
                    const AdditionalData             additional_data);

    /**
     * 帮助函数。循环所有的单元格，并在每个元素上并行地应用函数。这个函数在不使用MPI时使用。
     *
     */
    template <typename Functor, typename VectorType>
    void
    serial_cell_loop(const Functor &   func,
                     const VectorType &src,
                     VectorType &      dst) const;

    /**
     * 辅助函数。在所有的单元格上循环，并在每个元素上平行地应用函数。这个函数在使用MPI的时候使用。
     *
     */
    template <typename Functor>
    void
    distributed_cell_loop(
      const Functor &                                                      func,
      const LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA> &src,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA> &dst) const;

    /**
     * 这个函数不应该被调用。调用它将导致一个内部错误。这个函数的存在只是因为cell_loop需要distributed_cell_loop()的存在，因为
     * LinearAlgebra::CUDAWrappers::Vector.
     *
     */
    template <typename Functor>
    void
    distributed_cell_loop(
      const Functor &                                    func,
      const LinearAlgebra::CUDAWrappers::Vector<Number> &src,
      LinearAlgebra::CUDAWrappers::Vector<Number> &      dst) const;

    /**
     * 帮助函数。将 @p src 的约束项的值复制到 @p dst.
     * 这个函数在不使用MPI时使用。
     *
     */
    template <typename VectorType>
    void
    serial_copy_constrained_values(const VectorType &src,
                                   VectorType &      dst) const;

    /**
     * 帮助功能。将 @p src 的受限项的值复制到 @p dst.
     * 中，该函数在使用MPI时使用。
     *
     */
    void
    distributed_copy_constrained_values(
      const LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA> &src,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA> &dst) const;

    /**
     * 这个函数永远不应该被调用。调用它将导致一个内部错误。这个函数的存在只是因为copy_constrained_values需要distributed_copy_constrained_values()存在
     * LinearAlgebra::CUDAWrappers::Vector.  。
     *
     */
    void
    distributed_copy_constrained_values(
      const LinearAlgebra::CUDAWrappers::Vector<Number> &src,
      LinearAlgebra::CUDAWrappers::Vector<Number> &      dst) const;

    /**
     * 帮助功能。将 @p dst 的约束项设置为 @p val.
     * ，该函数在不使用MPI时使用。
     *
     */
    template <typename VectorType>
    void
    serial_set_constrained_values(const Number val, VectorType &dst) const;

    /**
     * 帮助功能。将 @p dst 的受限条目设置为 @p val.
     * 当使用MPI时，该函数被使用。
     *
     */
    void
    distributed_set_constrained_values(
      const Number                                                   val,
      LinearAlgebra::distributed::Vector<Number, MemorySpace::CUDA> &dst) const;

    /**
     * 这个函数永远不应该被调用。调用它将导致一个内部错误。这个函数的存在只是因为set_constrained_values需要distributed_set_constrained_values()存在
     * LinearAlgebra::CUDAWrappers::Vector.  。
     *
     */
    void
    distributed_set_constrained_values(
      const Number                                 val,
      LinearAlgebra::CUDAWrappers::Vector<Number> &dst) const;

    /**
     * 与该对象相关的唯一ID。
     *
     */
    int my_id;

    /**
     * 使用的并行化方案，在自由度或单元上的并行化。
     *
     */
    ParallelizationScheme parallelization_scheme;

    /**
     * 如果为真，使用图形着色。否则，使用原子操作。图形着色可以确保位数的可重复性，但是在Pascal和新的架构上会比较慢。
     *
     */
    bool use_coloring;

    /**
     * 将MPI通信与计算重叠。这需要CUDA感知的MPI，并且use_coloring必须为假。
     *
     */
    bool overlap_communication_computation;

    /**
     * 总自由度数。
     *
     */
    types::global_dof_index n_dofs;

    /**
     * 所用有限元的度数。
     *
     */
    unsigned int fe_degree;

    /**
     * 每个单元的自由度数。
     *
     */
    unsigned int dofs_per_cell;

    /**
     * 受限自由度的数量。
     *
     */
    unsigned int n_constrained_dofs;

    /**
     * 每个单元的正交点的数量。
     *
     */
    unsigned int q_points_per_cell;

    /**
     * 图形着色算法产生的颜色数量。
     *
     */
    unsigned int n_colors;

    /**
     * 每种颜色的单元格数量。
     *
     */
    std::vector<unsigned int> n_cells;

    /**
     * 与每种颜色的单元格相关的正交点的指针矢量。
     *
     */
    std::vector<point_type *> q_points;

    /**
     * 将本地向量中的位置映射到全局向量中的位置。
     *
     */
    std::vector<types::global_dof_index *> local_to_global;

    /**
     * 与每种颜色的单元格相关的反雅各布系数的指针向量。
     *
     */
    std::vector<Number *> inv_jacobian;

    /**
     * 指向每个颜色的单元格相关的雅各布系数乘以权重的矢量。
     *
     */
    std::vector<Number *> JxW;

    /**
     * 指向受限自由度的指针。
     *
     */
    types::global_dof_index *constrained_dofs;

    /**
     * 决定在一个给定单元上设置约束的掩码。
     *
     */
    std::vector<unsigned int *> constraint_mask;

    /**
     * 与不同颜色相关的网格尺寸。网格尺寸用于启动CUDA内核。
     *
     */
    std::vector<dim3> grid_dim;

    /**
     * 与不同颜色相关的块尺寸。块的尺寸用于启动CUDA内核。
     *
     */
    std::vector<dim3> block_dim;

    /**
     * 用于cell_loop中的分布式矢量的共享指针。当不使用MPI时，该指针为空。
     *
     */
    std::shared_ptr<const Utilities::MPI::Partitioner> partitioner;

    /**
     * 每块单元格（由函数cell_per_block_shmem()决定）。
     *
     */
    unsigned int cells_per_block;

    /**
     * 用于启动CUDA内核的网格尺寸 in_constrained_values-operations.
     *
     */
    dim3 constraint_grid_dim;

    /**
     * 用于启动CUDA内核的块状尺寸_受限值-操作。
     *
     */
    dim3 constraint_block_dim;

    /**
     * 填充物的长度（大于或等于线程数的最接近2的幂）。
     *
     */
    unsigned int padding_length;

    /**
     * 每个颜色的行开始。
     *
     */
    std::vector<unsigned int> row_start;

    /**
     * 指向与该对象相关的DoFHandler的指针。
     *
     */
    const DoFHandler<dim> *dof_handler;

    /**
     * 本地拥有的活动单元的彩色图谱。
     *
     */
    std::vector<std::vector<CellFilter>> graph;

    friend class internal::ReinitHelper<dim, Number>;
  };



  // TODO find a better place to put these things

  /**
   * 用于将共享内存传递给一般用户函数的结构。
   *
   */
  template <int dim, typename Number>
  struct SharedData
  {
    /**
     * 构造函数。
     *
     */
    __device__
    SharedData(Number *vd, Number *gq[dim])
      : values(vd)
    {
      for (int d = 0; d < dim; ++d)
        gradients[d] = gq[d];
    }

    /**
     * 用于dof和quad值的共享内存。
     *
     */
    Number *values;

    /**
     * 参考坐标系中计算的梯度的共享内存。
     * 每个方向的梯度以数组结构的形式保存，也就是说，首先，X方向的所有梯度都会被保存。
     *
     */
    Number *gradients[dim];
  };



  // This function determines the number of cells per block, possibly at compile
  // time (by virtue of being 'constexpr')
  // TODO this function should be rewritten using meta-programming
  __host__ __device__ constexpr unsigned int
           cells_per_block_shmem(int dim, int fe_degree)
  {
     /* clang-format off */ 
    // We are limiting the number of threads according to the
    // following formulas:
    //  - in 2D: `threads = cells * (k+1)^d <= 4*CUDAWrappers::warp_size`
    //  - in 3D: `threads = cells * (k+1)^d <= 2*CUDAWrappers::warp_size`
    return dim==2 ? (fe_degree==1 ? CUDAWrappers::warp_size :    // 128
                     fe_degree==2 ? CUDAWrappers::warp_size/4 :  //  72
                     fe_degree==3 ? CUDAWrappers::warp_size/8 :  //  64
                     fe_degree==4 ? CUDAWrappers::warp_size/8 :  // 100
                     1) :
           dim==3 ? (fe_degree==1 ? CUDAWrappers::warp_size/4 :  //  64
                     fe_degree==2 ? CUDAWrappers::warp_size/16 : //  54
                     1) : 1;
     /* clang-format on */ 
  }


   /*----------------------- Helper functions ---------------------------------*/ 
  /**
   * 计算给定线程的本地单元中的正交点索引。      @relates
   * CUDAWrappers::MatrixFree 。
   *
   */
  template <int dim>
  __device__ inline unsigned int
  q_point_id_in_cell(const unsigned int n_q_points_1d)
  {
    return (dim == 1 ?
              threadIdx.x % n_q_points_1d :
              dim == 2 ?
              threadIdx.x % n_q_points_1d + n_q_points_1d * threadIdx.y :
              threadIdx.x % n_q_points_1d +
                  n_q_points_1d * (threadIdx.y + n_q_points_1d * threadIdx.z));
  }



  /**
   * 返回一个给定线程的本地正交点索引。该索引只对一个给定的MPI进程是唯一的。
   * @relates   CUDAWrappers::MatrixFree
   *
   */
  template <int dim, typename Number>
  __device__ inline unsigned int
  local_q_point_id(
    const unsigned int                                          cell,
    const typename CUDAWrappers::MatrixFree<dim, Number>::Data *data,
    const unsigned int                                          n_q_points_1d,
    const unsigned int                                          n_q_points)
  {
    return (data->row_start / data->padding_length + cell) * n_q_points +
           q_point_id_in_cell<dim>(n_q_points_1d);
  }



  /**
   * 返回与给定线程相关的正交点。      @relates
   * CUDAWrappers::MatrixFree  *返回与指定线程相关的正交点。
   *
   */
  template <int dim, typename Number>
  __device__ inline typename CUDAWrappers::MatrixFree<dim, Number>::point_type &
  get_quadrature_point(
    const unsigned int                                          cell,
    const typename CUDAWrappers::MatrixFree<dim, Number>::Data *data,
    const unsigned int                                          n_q_points_1d)
  {
    return *(data->q_points + data->padding_length * cell +
             q_point_id_in_cell<dim>(n_q_points_1d));
  }

  /**
   * 传递给内核的结构。它用于将所有必要的信息从CPU传递给GPU。
   *
   */
  template <int dim, typename Number>
  struct DataHost
  {
    /**
     * 正交点的矢量。
     *
     */
    std::vector<Point<dim, Number>> q_points;

    /**
     * 将本地向量中的位置映射到全局向量中的位置。
     *
     */
    std::vector<types::global_dof_index> local_to_global;

    /**
     * 反雅各布系数的矢量。
     *
     */
    std::vector<Number> inv_jacobian;

    /**
     * 雅各布系数乘以权重的向量。
     *
     */
    std::vector<Number> JxW;

    /**
     * 相关MatrixFree对象的ID。
     *
     */
    unsigned int id;

    /**
     * 单元的数量。
     *
     */
    unsigned int n_cells;

    /**
     * 填充物的长度。
     *
     */
    unsigned int padding_length;

    /**
     * 行开始（包括填充物）。
     *
     */
    unsigned int row_start;

    /**
     * 决定在一个给定的单元格上设置约束的掩码。
     *
     */
    std::vector<unsigned int> constraint_mask;

    /**
     * 如果为真，则使用图形着色，我们可以简单地添加到目的地向量中。否则，使用原子操作。
     *
     */
    bool use_coloring;
  };



  /**
   * 将 @p data 从设备复制到设备。  @p update_flags 应该与
   * MatrixFree::AdditionalData.  @relates  CUDAWrappers::MatrixFree
   * 中使用的相同。
   *
   */
  template <int dim, typename Number>
  DataHost<dim, Number>
  copy_mf_data_to_host(
    const typename dealii::CUDAWrappers::MatrixFree<dim, Number>::Data &data,
    const UpdateFlags &update_flags)
  {
    DataHost<dim, Number> data_host;

    data_host.id             = data.id;
    data_host.n_cells        = data.n_cells;
    data_host.padding_length = data.padding_length;
    data_host.row_start      = data.row_start;
    data_host.use_coloring   = data.use_coloring;

    const unsigned int n_elements =
      data_host.n_cells * data_host.padding_length;
    if (update_flags & update_quadrature_points)
      {
        data_host.q_points.resize(n_elements);
        Utilities::CUDA::copy_to_host(data.q_points, data_host.q_points);
      }

    data_host.local_to_global.resize(n_elements);
    Utilities::CUDA::copy_to_host(data.local_to_global,
                                  data_host.local_to_global);

    if (update_flags & update_gradients)
      {
        data_host.inv_jacobian.resize(n_elements * dim * dim);
        Utilities::CUDA::copy_to_host(data.inv_jacobian,
                                      data_host.inv_jacobian);
      }

    if (update_flags & update_JxW_values)
      {
        data_host.JxW.resize(n_elements);
        Utilities::CUDA::copy_to_host(data.JxW, data_host.JxW);
      }

    data_host.constraint_mask.resize(data_host.n_cells);
    Utilities::CUDA::copy_to_host(data.constraint_mask,
                                  data_host.constraint_mask);

    return data_host;
  }



  /**
   * 这个函数是local_q_point_id()的主机版本。      @relates
   * CUDAWrappers::MatrixFree
   *
   */
  template <int dim, typename Number>
  inline unsigned int
  local_q_point_id_host(const unsigned int           cell,
                        const DataHost<dim, Number> &data,
                        const unsigned int           n_q_points,
                        const unsigned int           i)
  {
    return (data.row_start / data.padding_length + cell) * n_q_points + i;
  }



  /**
   * 该函数是主机版的get_quadrature_point()。它假定MatrixFree<dim,
   * Number>::Data
   * 中的数据已经用copy_mf_data_to_host()复制到主机上了。
   * @relates   CUDAWrappers::MatrixFree 。
   *
   */
  template <int dim, typename Number>
  inline Point<dim, Number>
  get_quadrature_point_host(const unsigned int           cell,
                            const DataHost<dim, Number> &data,
                            const unsigned int           i)
  {
    return data.q_points[data.padding_length * cell + i];
  }


   /*----------------------- Inline functions ---------------------------------*/ 

#  ifndef DOXYGEN

  template <int dim, typename Number>
  inline const std::vector<std::vector<
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>>> &
  MatrixFree<dim, Number>::get_colored_graph() const
  {
    return graph;
  }



  template <int dim, typename Number>
  inline const std::shared_ptr<const Utilities::MPI::Partitioner> &
  MatrixFree<dim, Number>::get_vector_partitioner() const
  {
    return partitioner;
  }



  template <int dim, typename Number>
  inline const DoFHandler<dim> &
  MatrixFree<dim, Number>::get_dof_handler() const
  {
    Assert(dof_handler != nullptr, ExcNotInitialized());

    return *dof_handler;
  }

#  endif

} // namespace CUDAWrappers

DEAL_II_NAMESPACE_CLOSE

#endif
#endif


