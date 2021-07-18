//include/deal.II-translator/base/process_grid_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

#ifndef dealii_process_grid_h
#define dealii_process_grid_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>

#ifdef DEAL_II_WITH_SCALAPACK

DEAL_II_NAMESPACE_OPEN


// Forward declaration of class ScaLAPACKMatrix for ProcessGrid
#  ifndef DOXYGEN
template <typename NumberType>
class ScaLAPACKMatrix;
#  endif

namespace Utilities
{
  namespace MPI
  {
    /**
     * 一个负责设置二维处理器网格的类。
     * 例如，一个有5个进程的MPI通信器可以被安排成一个2x2的网格，其中第5个处理器是不活动的。
     * @code
     *    |   0     |   1
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
     * -----|
     *
     * ------- |-----
     * 0    |   P0    |  P1
     *    |         |
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
     * -----|
     *
     * ------- |-----
     * 1    |   P2    |  P3
     * @endcode
     * 这个类的共享指针被提供给ScaLAPACKMatrix矩阵，以执行块循环分布。
     * 请注意，这个类允许设置一个进程网格，它的MPI核数比通信器中的总核数少。
     * 目前，唯一可以使用ProcessGrid对象的地方是与ScaLAPACKMatrix对象相连接。
     *
     */
    class ProcessGrid
    {
    public:
      // Declare class ScaLAPACK as friend to provide access to private members.
      template <typename NumberType>
      friend class dealii::ScaLAPACKMatrix;

      /**
       * 给定 @p mpi_communicator. 的具有 @p n_rows 和 @p n_columns
       * 的进程网格的构造函数，行和列的乘积应该小于或等于
       * @p mpi_communicator. 中的总核数。
       *
       */
      ProcessGrid(const MPI_Comm &   mpi_communicator,
                  const unsigned int n_rows,
                  const unsigned int n_columns);

      /**
       * 在这种情况下，根据 @p n_rows_matrix,   @p n_columns_matrix,
       * @p row_block_size 和 @p column_block_size.
       * 中提供的目标矩阵的尺寸和块循环分布，启发式地选择进程网格，可以利用的最大MPI核心数是
       * $\min\{\frac{M}{MB}\frac{N}{NB}, Np\}$  ，其中 $M,N$
       * ]是矩阵尺寸， $MB,NB$ 是块大小， $Np$ 是 @p
       * mpi_communicator. 中的进程数。
       * 这个函数然后创建一个二维处理器网格，假设进程行
       * $p$ 和列 $q$ 的数量比例等于矩阵尺寸 $M$ 和 $N$
       * 的比例。            例如，一个方形矩阵 $640x640$
       * 的块大小为 $32$ ， @p mpi_communicator
       * 有11个核心，将导致 $3x3$ 的过程网格。
       *
       */
      ProcessGrid(const MPI_Comm &   mpi_communicator,
                  const unsigned int n_rows_matrix,
                  const unsigned int n_columns_matrix,
                  const unsigned int row_block_size,
                  const unsigned int column_block_size);

      /**
       * 销毁器。
       *
       */
      ~ProcessGrid();

      /**
       * 返回进程网格中的行数。
       *
       */
      unsigned int
      get_process_grid_rows() const;

      /**
       * 返回进程网格中的列数。
       *
       */
      unsigned int
      get_process_grid_columns() const;

      /**
       * 返回该进程在进程网格中的行数。
       * 对于不活动的进程，它是负数。
       *
       */
      int
      get_this_process_row() const;

      /**
       * 返回该进程在进程网格中的列。
       * 对于不活动的进程是负数。
       *
       */
      int
      get_this_process_column() const;

      /**
       * 将 @p count 的值从等级为0的进程的 @p value
       * 处开始存储，然后发送给不在进程网格内的进程。
       *
       */
      template <typename NumberType>
      void
      send_to_inactive(NumberType *value, const int count = 1) const;

      /**
       * 如果进程在网格内活动，则返回 <code>true</code> 。
       *
       */
      bool
      is_process_active() const;

    private:
      /**
       * 一个私有的构造函数，它将网格尺寸作为一个
       * <code>std::pair</code> 。
       *
       */
      ProcessGrid(const MPI_Comm &                             mpi_communicator,
                  const std::pair<unsigned int, unsigned int> &grid_dimensions);

      /**
       * 一个与所有进程（活动和非活动）的MPI通信器。
       *
       */
      MPI_Comm mpi_communicator;

      /**
       * 一个带有不活动进程和等级为零的进程的MPI通信器。
       *
       */
      MPI_Comm mpi_communicator_inactive_with_root;

      /**
       * BLACS上下文。这等同于MPI通信器，被ScaLAPACK使用。
       *
       */
      int blacs_context;

      /**
       * 这个MPI进程的等级。
       *
       */
      const unsigned int this_mpi_process;

      /**
       * MPI进程的总数量。
       *
       */
      const unsigned int n_mpi_processes;

      /**
       * 进程网格中的行数。
       *
       */
      int n_process_rows;

      /**
       * 进程网格中的列数。
       *
       */
      int n_process_columns;

      /**
       * 这个进程在网格中的行。
       * 对于不活动的进程，它是负数。
       *
       */
      int this_process_row;

      /**
       * 这个进程在网格中的列。
       * 对于不活动的过程，它是负数。
       *
       */
      int this_process_column;

      /**
       * 一个标志，对于2D进程网格中的进程来说是真的。
       *
       */
      bool mpi_process_is_active;
    };

     /*--------------------- Inline functions --------------------------------*/ 

#  ifndef DOXYGEN

    inline unsigned int
    ProcessGrid::get_process_grid_rows() const
    {
      return n_process_rows;
    }



    inline unsigned int
    ProcessGrid::get_process_grid_columns() const
    {
      return n_process_columns;
    }



    inline int
    ProcessGrid::get_this_process_row() const
    {
      return this_process_row;
    }



    inline int
    ProcessGrid::get_this_process_column() const
    {
      return this_process_column;
    }



    inline bool
    ProcessGrid::is_process_active() const
    {
      return mpi_process_is_active;
    }


#  endif // ifndef DOXYGEN

  } // end of namespace MPI

} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SCALAPACK

#endif


