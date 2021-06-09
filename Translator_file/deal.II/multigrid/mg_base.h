//include/deal.II-translator/multigrid/mg_base_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2021 by the deal.II authors
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

#ifndef dealii_mg_base_h
#define dealii_mg_base_h

/* 该文件包含Multigrid使用的一些抽象基类。

* 
*
*/

#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/vector.h>


DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup mg */ 
 /*@{*/ 


/**
 * 多级矩阵基类。这个类设置了多级算法所需的接口。它与实际的矩阵类型没有关系，只接受向量类作为模板参数。
 * 通常，派生类 mg::Matrix,
 * 对矩阵的MGLevelObject进行操作，对于应用来说已经足够。
 *
 *
 */
template <typename VectorType>
class MGMatrixBase : public Subscriptor
{
public:
  /*虚拟解构器。 
*
*/
  virtual ~MGMatrixBase() override = default;

  /**
   * 矩阵-向量-乘法在一定水平上。
   *
   */
  virtual void
  vmult(const unsigned int level,
        VectorType &       dst,
        const VectorType & src) const = 0;

  /**
   * 在某一层次上增加矩阵-向量-乘法。
   *
   */
  virtual void
  vmult_add(const unsigned int level,
            VectorType &       dst,
            const VectorType & src) const = 0;

  /**
   * 在某一层次上进行矩阵-向量-乘法的转置。
   *
   */
  virtual void
  Tvmult(const unsigned int level,
         VectorType &       dst,
         const VectorType & src) const = 0;

  /**
   * 在某一层次上增加转置的矩阵-向量-乘法。
   *
   */
  virtual void
  Tvmult_add(const unsigned int level,
             VectorType &       dst,
             const VectorType & src) const = 0;

  /**
   * 返回存储矩阵的最小级别。
   *
   */
  virtual unsigned int
  get_minlevel() const = 0;

  /**
   * 返回存储矩阵的最小级别。
   *
   */
  virtual unsigned int
  get_maxlevel() const = 0;
};


/**
 * 粗略网格求解器的基类。
 * 它定义了虚拟括号运算符，是多网格方法使用的接口。任何实现都将由派生类来完成。
 *
 *
 */
template <typename VectorType>
class MGCoarseGridBase : public Subscriptor
{
public:
  /**
   * 虚拟解构器。
   *
   */
  virtual ~MGCoarseGridBase() override = default;

  /**
   * 解算器。
   *
   */
  virtual void
  operator()(const unsigned int level,
             VectorType &       dst,
             const VectorType & src) const = 0;
};


/**
 * 基类用于声明在多网格背景下实现向量延长和限制的具体类所需的操作。这个类是抽象的，没有这些操作的实现。
 * 有几个派生类，反映了细格离散化和多层次实现的向量类型和编号是独立的。
 * 如果你对一个单一的PDE或对你的完整方程组使用多网格，你将和多网格一起使用MGTransferPrebuilt。在精细网格上使用的矢量类型以及用于多级操作的矢量类型可以是Vector或BlockVector。在这两种情况下，MGTransferPrebuilt将对解决方案的所有组成部分进行操作。
 *
 *
 * @todo  更新以下文档，因为它没有反映结构的最新变化。
 * 对于混合系统，可能需要只对单个组件或某些组件做多网格。类MGTransferSelect和MGTransferBlock处理这些情况。
 * MGTransferSelect用于对单一组件使用多重网格（在矢量对象上），可能使用<tt>mg_target_component</tt>进行分组。
 * MGTransferBlock类处理你的多栅格方法对BlockVector对象进行操作的情况。这些对象可以包含完整系统的全部或连续的块集。由于大多数平滑器不能对块结构进行操作，因此不清楚这种情况是否真的有用。因此，需要时将提供这种情况的测试实现。
 *
 *
 */
template <typename VectorType>
class MGTransferBase : public Subscriptor
{
public:
  /**
   * 解构器。在这里不做任何事情，但无论如何需要声明为虚拟。
   *
   */
  virtual ~MGTransferBase() override = default;

  /**
   * 将一个向量从<tt>to_level-1</tt>级延长到<tt>to_level</tt>级。之前<tt>dst</tt>的内容被覆盖。
   * @arg
   * src是一个向量，其元素数与所涉及的较粗层次上的自由度相同。
   * @arg  dst有多少个元素，就有多少个更细层次的自由度。
   *
   */
  virtual void
  prolongate(const unsigned int to_level,
             VectorType &       dst,
             const VectorType & src) const = 0;

  /**
   * 将一个向量从<tt>to_level-1</tt>级延长到<tt>to_level</tt>级，加到<tt>dst</tt>的前一个内容。
   * @arg
   * src是一个向量，其元素数与所涉及的较粗层次上的自由度相同。
   * @arg  dst有多少个元素，就有多少个更细层次的自由度。
   *
   */
  virtual void
  prolongate_and_add(const unsigned int to_level,
                     VectorType &       dst,
                     const VectorType & src) const;

  /**
   * 限制一个从<tt>from_level</tt>层到<tt>from_level-1</tt>层的向量，并将这个限制添加到<tt>dst</tt>。如果<tt>from_level</tt>层的单元所覆盖的区域小于<tt>from_level-1</tt>层的区域（局部细化），那么<tt>dst</tt>中的一些自由度是有效的，将不会被改变。对于其他自由度，限制的结果将被添加。
   * @arg
   * src是一个向量，其元素数量与较细层次上的自由度相同
   * @arg dst的元素数量与较粗层次上的自由度相同。
   *
   */
  virtual void
  restrict_and_add(const unsigned int from_level,
                   VectorType &       dst,
                   const VectorType & src) const = 0;
};



/**
 * 多网格平滑器的基类。除了定义多网格方法所使用的接口外，什么都不做。
 * 平滑器接口提供了两个方法，一个smooth()方法和一个apply()方法。出于效率的考虑，多网格预处理程序接口对这两种方法进行了区分。在进入预处理程序操作时，向量
 * @p u 需要被设置为零，平滑的开始是在 @p rhs
 * 向量上简单应用平滑器。这个方法是由这个类的apply()方法提供的。这与首先将
 * @p u
 * 设置为零，然后调用smooth()是一样的，但是对于许多类来说，单独的apply()接口更有效率，因为它可以跳过一个矩阵-向量乘积。
 * 在多网格预处理接口中，apply()方法被用于预平滑操作，因为解向量中以前的内容需要为新进来的残差所覆盖。另一方面，所有后续的操作都需要平滑向量中已经存在的内容
 * @p u 给定的右手边，这是由smooth()完成的。
 *
 *
 */
template <typename VectorType>
class MGSmootherBase : public Subscriptor
{
public:
  /**
   * 虚拟解构器。
   *
   */
  virtual ~MGSmootherBase() override = default;

  /**
   * 释放矩阵。
   *
   */
  virtual void
  clear() = 0;

  /**
   * 平滑函数，给定右手边的向量 @p rhs. ，对 @p u
   * 中的内容进行平滑处理，这是多网格方法中使用的函数。
   *
   */
  virtual void
  smooth(const unsigned int level,
         VectorType &       u,
         const VectorType & rhs) const = 0;

  /**
   * 与smooth()函数相反，这个函数应用了平滑的动作，覆盖了向量u中以前的内容，这个函数必须等同于以下代码
   * @code
   * u = 0;
   * smooth(level, u, rhs);
   * @endcode
   * 但通常可以比前者更有效地实现。如果一个特定的平滑器没有覆盖apply()方法，就会使用这里描述的默认实现。
   * 在多网格预处理接口中，apply()方法被用于预平滑操作，因为解向量中以前的内容需要被覆盖以获得新的传入残差。另一方面，所有后续的操作都需要平滑已经存在于向量
   * @p u 中的内容，给定的右手边，这是由smooth()完成的。
   *
   */
  virtual void
  apply(const unsigned int level, VectorType &u, const VectorType &rhs) const;
};

 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


