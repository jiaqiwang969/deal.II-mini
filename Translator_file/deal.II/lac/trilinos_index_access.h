//include/deal.II-translator/lac/trilinos_index_access_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

#ifndef dealii_trilinos_index_access_h
#define dealii_trilinos_index_access_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_TRILINOS

#  include <deal.II/base/types.h>

#  include <Epetra_BlockMap.h>
#  include <Epetra_CrsGraph.h>
#  include <Epetra_CrsMatrix.h>
#  include <Epetra_MultiVector.h>

DEAL_II_NAMESPACE_OPEN

namespace TrilinosWrappers
{
  /**
   * 一个辅助函数，用于查询Epetra_BlockMap对象的大小，并调用32位或64位函数来获得地图中全局元素的数量。
   *
   */
  inline TrilinosWrappers::types::int_type
  n_global_elements(const Epetra_BlockMap &map)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return map.NumGlobalElements64();
#  else
    return map.NumGlobalElements();
#  endif
  }

  /**
   * 一个辅助函数，通过调用32位或64位函数找到调用处理器上的最小全局索引值。
   *
   */
  inline TrilinosWrappers::types::int_type
  min_my_gid(const Epetra_BlockMap &map)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return map.MinMyGID64();
#  else
    return map.MinMyGID();
#  endif
  }

  /**
   * 一个辅助函数，通过调用32位或64位函数，在调用的处理器上找到最大的全局索引值。
   *
   */
  inline TrilinosWrappers::types::int_type
  max_my_gid(const Epetra_BlockMap &map)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return map.MaxMyGID64();
#  else
    return map.MaxMyGID();
#  endif
  }

  /**
   * 一个辅助函数，通过调用32位或64位函数，将本地索引转换为全局索引。
   *
   */
  inline TrilinosWrappers::types::int_type
  global_index(const Epetra_BlockMap &               map,
               const dealii::types::global_dof_index i)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return map.GID64(i);
#  else
    return map.GID(i);
#  endif
  }

  /**
   * 一个辅助函数，通过调用32位或64位函数，返回一个指向包含分配给当前进程的全局索引的数组的指针。
   *
   */
  inline TrilinosWrappers::types::int_type *
  my_global_elements(const Epetra_BlockMap &map)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return map.MyGlobalElements64();
#  else
    return map.MyGlobalElements();
#  endif
  }

  /**
   * 一个辅助函数，通过调用32位或64位函数，找到全局的行数。
   *
   */
  inline TrilinosWrappers::types::int_type
  n_global_rows(const Epetra_CrsGraph &graph)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return graph.NumGlobalRows64();
#  else
    return graph.NumGlobalRows();
#  endif
  }

  /**
   * 一个辅助函数，通过调用32位或64位函数来找到全局的列数。
   *
   */
  inline TrilinosWrappers::types::int_type
  n_global_cols(const Epetra_CrsGraph &graph)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return graph.NumGlobalCols64();
#  else
    return graph.NumGlobalCols();
#  endif
  }

  /**
   * 一个辅助函数，通过调用32位或64位函数来查找全局条目数。
   *
   */
  inline TrilinosWrappers::types::int_type
  n_global_entries(const Epetra_CrsGraph &graph)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return graph.NumGlobalEntries64();
#  else
    return graph.NumGlobalEntries();
#  endif
  }

  /**
   * 一个辅助函数，通过调用32位或64位函数来查找全局行索引。
   *
   */
  inline TrilinosWrappers::types::int_type
  global_row_index(const Epetra_CrsMatrix &              matrix,
                   const dealii::types::global_dof_index i)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return matrix.GRID64(i);
#  else
    return matrix.GRID(i);
#  endif
  }

  /**
   * 一个通过调用32位或64位函数找到全局列索引的辅助函数。
   *
   */
  inline TrilinosWrappers::types::int_type
  global_column_index(const Epetra_CrsMatrix &              matrix,
                      const dealii::types::global_dof_index i)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return matrix.GCID64(i);
#  else
    return matrix.GCID(i);
#  endif
  }

  /**
   * 一个辅助函数，通过调用32位或64位函数来查找向量的全局长度。
   *
   */
  inline TrilinosWrappers::types::int_type
  global_length(const Epetra_MultiVector &vector)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return vector.GlobalLength64();
#  else
    return vector.GlobalLength();
#  endif
  }

  /**
   * 一个辅助函数，通过调用32位或64位函数找到全局的行数。
   *
   */
  inline TrilinosWrappers::types::int_type
  n_global_rows(const Epetra_RowMatrix &matrix)
  {
#  ifdef DEAL_II_WITH_64BIT_INDICES
    return matrix.NumGlobalRows64();
#  else
    return matrix.NumGlobalRows();
#  endif
  }
} // namespace TrilinosWrappers

DEAL_II_NAMESPACE_CLOSE
#endif // DEAL_II_WITH_TRILINOS
#endif // dealii_trilinos_index_access_h


