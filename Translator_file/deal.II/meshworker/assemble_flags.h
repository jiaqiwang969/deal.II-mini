//include/deal.II-translator/meshworker/assemble_flags_0.txt
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

#ifndef dealii_mesh_worker_assemble_flags_h
#define dealii_mesh_worker_assemble_flags_h


#include <deal.II/base/config.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup MeshWorker */ 
 /*@{*/ 

namespace MeshWorker
{
  /**
   * 给予Mesh_loop()函数的枚举类型，告诉该函数哪些元素需要被组装起来。
   * 你可以用位法或 <code>operator|(AssembleFlags,AssembleFlags)</code>
   * 串联的方式选择多个标志。
   *
   */
  enum AssembleFlags
  {
    /**
     * 什么都不做。
     *
     */
    assemble_nothing = 0,
    /**
     * 在本地拥有的单元上进行组合。
     *
     */
    assemble_own_cells = 0x0001,
    /**
     * 集合在幽灵细胞上。
     *
     */
    assemble_ghost_cells = 0x0002,
    /**
     * 在两个本地拥有的单元格之间的内部面上进行集合，每个面只访问一次。
     *
     */
    assemble_own_interior_faces_once = 0x0004,
    /**
     * 在两个本地拥有的单元格之间的内部面上集合，对每个内部面访问两次，从相邻的两个单元格中各访问一次。
     *
     */
    assemble_own_interior_faces_both = 0x0008,
    /**
     * 在一个本地拥有的单元和一个幽灵单元之间的面进行组装，确保只有一个进程将组装这些面（从较细的一面或具有较低mpi等级的进程）。
     *
     */
    assemble_ghost_faces_once = 0x0010,
    /**
     * 在一个本地拥有的单元和一个幽灵单元之间的面进行组装。两个进程都将组装这些面。请注意，它们永远不会在一个进程中从两边进行装配。
     *
     */
    assemble_ghost_faces_both = 0x0020,
    /**
     * 在本地拥有的单元格的边界面上进行组装。
     *
     */
    assemble_boundary_faces = 0x0040,

    /**
     * 默认情况下，我们会在面积分之前装配单元积分。如果指定了这个标志，单元格将在面和边界之后被装配。
     *
     */
    cells_after_faces = 0x0080,

    /**
     * 结合标志来决定是否对单元进行任何工作。
     *
     */
    work_on_cells = assemble_own_cells | assemble_ghost_cells,

    /**
     * 结合标志来决定是否对面进行任何工作。
     *
     */
    work_on_faces = assemble_own_interior_faces_once |
                    assemble_own_interior_faces_both |
                    assemble_ghost_faces_once | assemble_ghost_faces_both,

    /**
     * 结合标志来确定是否对边界面做了任何工作。
     *
     */
    work_on_boundary = assemble_boundary_faces,
  };


  /**
   * 输出运算符，将组合标志作为一组or'd文本值输出。
   * @ref AssembleFlags
   *
   */
  template <class StreamType>
  inline StreamType &
  operator<<(StreamType &s, AssembleFlags u)
  {
    s << " AssembleFlags";
    if (u & assemble_own_cells)
      s << "|own_cells";
    if (u & assemble_own_interior_faces_once)
      s << "|own_faces_once";
    if (u & assemble_own_interior_faces_both)
      s << "|own_faces_both";
    if (u & assemble_ghost_cells)
      s << "|ghost_cells";
    if (u & assemble_ghost_faces_once)
      s << "|ghost_faces_once";
    if (u & assemble_ghost_faces_both)
      s << "|ghost_faces_both";
    if (u & assemble_boundary_faces)
      s << "|boundary_faces";
    return s;
  }


  /**
   * 全局运算符，它返回一个对象，其中所有的位都被设置为第一或第二个参数中的设置。这个操作符的存在是因为如果它不存在，那么bit-or <tt>操作符|</tt>的结果将是一个整数，当我们试图将其分配给AssembleFlags类型的对象时，会引发编译器警告。
   * @ref AssembleFlags
   *
   */
  inline AssembleFlags
  operator|(AssembleFlags f1, AssembleFlags f2)
  {
    return static_cast<AssembleFlags>(static_cast<unsigned int>(f1) |
                                      static_cast<unsigned int>(f2));
  }



  /**
   * 全局运算符，将第二个参数的位也设置在第一个参数中。
   * @ref AssembleFlags
   *
   */
  inline AssembleFlags &
  operator|=(AssembleFlags &f1, AssembleFlags f2)
  {
    f1 = f1 | f2;
    return f1;
  }


  /**
   * 全局操作符，它返回一个对象，其中所有位都被设置在第一个和第二个参数中。这个操作符的存在是因为如果它不存在，那么比特和<tt>操作符&</tt>的结果将是一个整数，当我们试图将其分配给AssembleFlags类型的对象时，会引发编译器警告。
   * @ref AssembleFlags
   *
   */
  inline AssembleFlags operator&(AssembleFlags f1, AssembleFlags f2)
  {
    return static_cast<AssembleFlags>(static_cast<unsigned int>(f1) &
                                      static_cast<unsigned int>(f2));
  }


  /**
   * 全局操作符，如果第一个参数中的所有位没有在第二个参数中设置，则清除它们。
   * @ref AssembleFlags
   *
   */
  inline AssembleFlags &
  operator&=(AssembleFlags &f1, AssembleFlags f2)
  {
    f1 = f1 & f2;
    return f1;
  }
} // namespace MeshWorker

 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


