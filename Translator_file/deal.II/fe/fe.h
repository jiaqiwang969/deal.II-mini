//include/deal.II-translator/fe/fe_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_fe_h
#define dealii_fe_h

#include <deal.II/base/config.h>

#include <deal.II/fe/block_mask.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_base.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/lac/full_matrix.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
class FEValuesBase;
template <int dim, int spacedim>
class FEValues;
template <int dim, int spacedim>
class FEFaceValues;
template <int dim, int spacedim>
class FESubfaceValues;
template <int dim, int spacedim>
class FESystem;

/**
 * 这是任意尺寸的有限元的基类。它在成员变量和公共成员函数方面都声明了接口，通过这些接口可以访问有限元的具体实现的属性。这个接口通常由一些变量和函数组组成，可以大致划分为以下几个部分。
 *
 *
 *
 * - 关于有限元的基本信息，如每个顶点、边或单元的自由度数。这类数据被存储在FiniteElementData基类中。(尽管 FiniteElement::get_name() 的成员函数也属于这一类)。
 *
 *
 *
 * - 对参考单元 $[0,1]^d$ 上的形状函数及其导数的描述，如果一个元素确实是通过从参考单元到实际单元的形状函数映射来定义的。
 *
 *
 *
 * - 描述元素的形状函数与父单元或子单元（限制或延长）或相邻单元（用于悬挂节点约束），以及与定义在同一单元上的其他有限元空间（例如，当进行 $p$ 细化时）的关系的矩阵（和访问它们的函数）。
 *
 *
 *
 * - 描述单个形状函数属性的%函数，例如，一个 @ref vector_valued "矢量值有限元 "的形状函数的哪些 @ref GlossComponent "矢量分量 "是非零的，或者一个元素是否是 @ref GlossPrimitive "原始的"。
 *
 *
 *
 * - 对于插值型元素，如常见的 $Q_p$ 拉格朗日元素，描述其 @ref GlossSupport "支持点 "位置的数据。
 *
 *
 *
 * - 定义FEValues类接口的%函数，几乎总是用来从用户代码中访问有限元形状函数。
 * 下面几节更详细地讨论了许多这些概念，并概述了有限元的具体实现可以提供完整描述有限元空间所需细节的策略。
 * 一般来说，派生类提供这种信息的方式有三种。
 *
 *
 *
 * - 一些通常容易计算的字段，这些字段由本类的构造函数（或FiniteElementData基类的构造函数）初始化，因此派生类在调用本类构造函数的过程中必须进行计算。具体来说，上面提到的关于形状函数的基本信息和部分描述性信息就是这种情况。
 *
 *
 *
 * - 一些在库中广泛使用的常见矩阵，本类为这些矩阵提供了派生类的构造函数需要填充的受保护成员变量。在这个类中提供这些矩阵的目的是：（i）它们经常被使用，以及（ii）它们的计算成本很高。因此，只计算一次，而不是每次使用时都计算，是有意义的。在大多数情况下，当前类的构造函数已经将它们设置为正确的大小，因此派生类只需要填充它们。这方面的例子包括将一个单元的形状函数与相邻的、子单元和父单元的形状函数联系起来的矩阵。
 *
 *
 *
 * - 不常见的信息，或者取决于特定输入参数的信息，需要由派生类来实现。对于这些，这个基类只声明抽象的虚拟成员函数，然后派生类必须实现它们。这一类的例子包括计算参考单元上形状函数的值和导数的函数，对于这些函数不可能用表格来表示，因为有无限多的点可能要对它们进行评估。在某些情况下，派生类可能会选择干脆不实现<i>all</i>可能的接口（或者可能没有<i>yet</i>完整的实现）；对于不常见的函数，这时往往有一个派生类可以重载的成员函数来描述一个特定的特征是否被实现。一个例子是，一个元素是否实现了在 $hp$ 有限元背景下使用它的必要信息（见 @ref hp "HP-有限元支持"）。
 *
 *  <h3>Nomenclature</h3>
 * 有限元类必须定义大量的描述有限元空间的不同属性。下面的小节描述了一些将在下面的文档中使用的命名法。
 * <h4>Components and blocks</h4>
 * 另一方面，如果至少有一个形状函数在一个以上的向量分量中是非零的，那么我们把整个元素称为
 * "非原始"。然后可以用 FiniteElement::get_nonzero_components()
 * 来确定形状函数的哪些向量分量是非零的。形状函数的非零分量的数量由
 * FiniteElement::n_components(). 返回
 * 形状函数是否为非原始的可以通过
 * FiniteElement::is_primitive(). 查询。
 * 很多时候，人们可能想把线性系统分成几块，以便它们能反映出底层算子的结构。这通常不是基于向量分量，而是基于使用 @ref GlossBlock "块"
 * ，然后将结果用于BlockVector、BlockSparseMatrix、BlockMatrixArray等类型对象的子结构。如果你使用非原始元素，你不能通过
 * FiniteElement::system_to_component_index().
 * 确定块数，相反，你可以使用
 * FiniteElement::system_to_block_index(). 有限元素的块数可以通过
 * FiniteElement::n_blocks(). 确定。
 * 为了更好地说明这些概念，让我们考虑下面这个多元素系统的例子
 *
 * @code
 * FESystem<dim> fe_basis(FE_Q<dim>(2), dim, FE_Q<dim>(1),1);
 * @endcode
 * 与  <code>dim=2</code>
 * 。产生的有限元有3个分量：两个来自二次元，一个来自线性元。例如，如果这个系统被用来离散流体动力学问题，那么我们可以认为前两个分量代表矢量值的速度场，而最后一个分量则对应于标量压力场。如果不对自由度（DoF）进行重新编号，这个有限元将产生以下的局部DoF分布。
*  @image html fe_system_example.png DoF indices
 * 使用两个函数 FiniteElement::system_to_component_index() 和
 * FiniteElement::system_to_base_index() ，可以得到每个自由度 "i
 * "的以下信息。
 *
 * @code
 * const unsigned int component =
 * fe_basis.system_to_component_index(i).first;
 * const unsigned int within_base =
 * fe_basis.system_to_component_index(i).second;
 * const unsigned int base =
 * fe_basis.system_to_base_index(i).first.first;
 * const unsigned int multiplicity =
 * fe_basis.system_to_base_index(i).first.second;
 * const unsigned int within_base_  =
 * fe_basis.system_to_base_index(i).second; // same as above
 * @endcode
 * 这将导致。 | DoF | Component | Base element | Shape function within
 * base | Multiplicity | :----: | :--------: | :----------: |
 * :------------------------: | :----------: | | 0 | 0 | 0 | 0 | 0 | | 1 | 1 |
 * 0 | 0 | 1 | | 2 | 2 | 1 | 0 | 0 | | 3 | 0 | 0 | 1 | 0 | | 4 | 1 | 0 | 1 | 1
 * | | 5 | 2 | 1 | 1 | 0 | | 6 | 0 | 0 | 2 | 0 | | 7 | 1 | 0 | 2 | 1 | | 8 | 2
 * | 1 | 2 | 0 | | 9 | 0 | 0 | 3 | 0 | | 10 | 1 | 0 | 3 | 1 | |     11 | 2 | 1
 * | 3 | 0 | | 12 | 0 | 0 | 4 | 0 | | 13 | 1 | 0 | 4 | 1 | | 14 | 0 | 0 | 5 |
 * 0 | | 15 | 1 | 0 | 5 | 1 | | 16 | 0 | 0 | 6 | 0 | | 17 | 1 | 0 | 6 | 1 | |
 * 18 | 0 | 0 | 7 | 0 | | 19 | 1 | 0 | 7 | 1 | | 20 | 0 | 0 | 8 | 0 | | 21 | 1
 * | 0 | 8 | 1 |
 * 我们看到的是：这个元素上总共有22个自由度，分量从0到2不等。
 * 每个自由度对应于用于构建FESystem的两个基本元素之一。
 * $\mathbb Q_2$  或  $\mathbb Q_1$
 * 。由于FE_Q是原始元素，我们对二次元总共有9个不同的标量值形状函数，对线性元素有4个。最后，对于对应于第一基元的DoFs，多重性为0或1，意味着我们对速度场的
 * $x$ 和 $y$ 分量使用相同的标量值 $\mathbb Q_2 \otimes \mathbb
 * Q_2$  。对于对应于第二基元的DoFs倍率为零。 <h4>Support
 * points</h4>
 * 有限元经常通过定义一个多项式空间和一组对偶函数来定义。如果这些函数涉及到点的评价，那么这个元素就是
 * "内插
 * "的，可以通过在这些点上评价一个任意的（但足够平滑的）函数来内插到有限元空间。我们称这些点为
 * "支持点"。
 * 一个典型的代码片断是这样做的。
 *
 * @code
 * Quadrature<dim> dummy_quadrature (fe.get_unit_support_points());
 * FEValues<dim>   fe_values (mapping, fe, dummy_quadrature,
 *                          update_quadrature_points);
 * fe_values.reinit (cell);
 * Point<dim> mapped_point = fe_values.quadrature_point (i);
 * @endcode
 *
 * 或者，这些点可以被逐一转换。
 *
 * @code
 * const vector<Point<dim> > &unit_points =
 * fe.get_unit_support_points();
 *
 * Point<dim> mapped_point =
 * mapping.transform_unit_to_real_cell (cell, unit_points[i]);
 * @endcode
 *
 *
 * @note
 * 有限元的get_unit_support_points()函数的实现以与形状函数相同的顺序返回这些点。因此，上面访问的正交点也是以这种方式排序的。形状函数的顺序通常记录在各个有限元类的类文件中。
 *
 *  <h3>Implementing finite element spaces in derived classes</h3>
 * 下面的章节为在派生类中实现具体的有限元空间提供了一些更多的指导。这包括取决于你想提供的东西的维度的信息，然后是帮助在具体案例中生成信息的工具列表。
 * 需要注意的是，有一些中间类可以做很多完整描述有限元空间所需的工作。例如，FE_Poly、FE_PolyTensor和FE_PolyFace类在本质上构建了一个完整的有限元空间，如果你只向它们提供一个抽象的描述，你想在这个空间上构建一个元素。使用这些中间类通常会使实现有限元描述变得非常简单。
 * 一般来说，如果你想实现一个元素，你可能想先看看其他类似元素的实现。由于有限元接口的许多更复杂的部分与它们如何与映射、正交和FEValues类相互作用有关，你也会想通过 @ref
 * FE_vs_Mapping_vs_FEValues 文档模块来阅读。
 *
 *  <h4>Interpolation matrices in one dimension</h4>
 * 在一个空间维度上（即对于 <code>dim==1</code> 和
 * <code>spacedim</code>
 * 的任何值），实现当前基类接口的有限元类只需要设置#限制和#延长矩阵，这两个矩阵分别描述了一个单元上的有限元空间与它的父单元的空间以及它的子单元上的空间之间的插值。当前类的构造函数在一维中预设#interface_constraints矩阵（用于描述不同细化水平的单元之间的界面上的悬挂节点约束）的大小为零，因为1d中没有悬挂节点。
 * <h4>Interpolation matrices in two dimensions</h4>
 * 除了上面讨论的1D的字段外，如果有限元的自由度位于边或顶点上，还需要一个约束矩阵来描述悬挂节点的约束。这些约束由一个
 * $m\times n$
 * -矩阵#interface_constraints表示，其中<i>m</i>是没有角顶点的精炼边上的自由度数量（中间顶点上的那些自由度加上两条线上的自由度），<i>n</i>是未精炼边的自由度（两条顶点上的自由度加上线上的自由度）。因此，该矩阵是一个矩形的。#interface_constraints矩阵的
 * $m\times n$
 * 大小也可以通过interface_constraints_size()函数访问。
 * dofs在未精炼侧的矩阵索引上的映射如下：让 $d_v$
 * 是一个顶点上的dofs数， $d_l$ 是一条线上的dofs数，那么
 * $n=0...d_v-1$ 指的是未精炼线的零顶点上的dofs，
 * $n=d_v...2d_v-1$ 指的是顶点一上的， $n=2d_v...2d_v+d_l-1$
 * 是线上的。 类似地， $m=0...d_v-1$
 * 指的是精炼边中间顶点上的道夫（子线零的顶点一，子线一的顶点零），
 * $m=d_v...d_v+d_l-1$ 指的是子线零上的道夫，
 * $m=d_v+d_l...d_v+2d_l-1$ 是指子线一上的道夫。
 * 请注意，我们不需要为精炼线末端顶点的道夫保留空间，因为这些道夫必须一对一地映射到未精炼线顶点的相应道夫。
 * 通过这种结构，子面的自由度被限制在父面的自由度上。这样提供的信息通常被
 * DoFTools::make_hanging_node_constraints() 函数所消耗。
 *
 *
 * @note  这些矩阵所描述的悬挂节点约束只与相邻（但不同的细化）单元上使用相同的有限元空间的情况有关。一个面的不同侧面的有限元空间不同的情况，即 $hp$ 情况（见 @ref hp "hp-有限元支持"
 * ）由单独的函数处理。参见
 * FiniteElement::get_face_interpolation_matrix() 和
 * FiniteElement::get_subface_interpolation_matrix() 函数。
 *
 *  <h4>Interpolation matrices in three dimensions</h4>
 * 对于接口约束，3d情况与2d情况类似。母面上的指数 $n$
 * 的编号很明显，并保持四边形上自由度的通常编号。
 * 精炼面内部的自由度的索引 $m$ 的编号如下：让 $d_v$ 和
 * $d_l$ 同上， $d_q$
 * 为每个四边形（因此也是每个面）的自由度数量，那么
 * $m=0...d_v-1$ 表示中心顶点上的自由度， $m=d_v...5d_v-1$
 * ]为四边形边界线中心顶点上的自由度， $m=5d_v..5d_v+4*d_l-1$
 * 为连接中心顶点和母面外边界的四条线上的自由度，
 * $m=5d_v+4*d_l...5d_v+4*d_l+8*d_l-1$
 * 为四边形周围小线上的自由度，
 * $m=5d_v+12*d_l...5d_v+12*d_l+4*d_q-1$
 * 为四个子面的自由度。注意四边形边界处的线的方向，如下图所示。
 * 十二条线和四个子面的顺序可以从下面的草图中提取，其中描绘了不同自由度组的总体顺序。
 *
 * @verbatim
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
 * --15--4--16--*
 *  |      |      |
 *  10 19  6  20  12
 *  |      |      |
 *  1--7---0--8---2
 *  |      |      |
 *  9  17  5  18  11
 *  |      |      |
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
 * --13--3--14--*
 * @endverbatim
 * 顶点和线的编号，以及线内子面的编号与《三角形》中描述的一致。因此，这种编号是根据面的不同，分别从外部和内部看到的。
 * 三维的情况下，对于想要实现约束矩阵的派生类来说，有一些陷阱可以利用。考虑一下下面的情况。
 *
 * @verbatim
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
 * -------*
 *       /       /|
 *      /       / |
 *     /       /  |
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
 * -------*   |
 *    |       |
 *
 *
 *
 *
 *
 * -------*
 *    |       |  /       /|
 *    |   1   | /       / |
 *    |       |/       /  |
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
 * -------*-------*   |
 *    |       |       |
 *    |       |       |  /
 *    |   2   |   3   | /
 *    |       |       |/
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
 * -------*-------*
 * @endverbatim
 * 现在假设我们想细化单元格2。我们最终会有两个带有悬挂节点的面，即单元格1和2之间的面，以及单元格2和3之间的面。必须对这两个面的自由度进行约束。问题是，现在有一条边（单元格2的右上方）是两个面的一部分。因此，这条边上的悬挂节点被约束了两次，一次来自两个面。为了有意义，这些约束当然必须是一致的：两个面都必须将边上的悬空节点约束到粗边上的相同节点上（而且只约束到边上，因为这样就不能约束到面的其他部分的节点），而且它们必须以相同的权重这样做。这有时是很棘手的，因为边上的节点可能有不同的局部数字。
 * 对于约束矩阵来说，这意味着：如果一个面的一条边上的自由度被同一条边上的其他一些节点以某种权重约束，那么这些权重必须与其他三条边上的约束节点相对于这些边上的相应节点的权重完全相同。如果不是这样，你会在AffineConstraints类中遇到麻烦，该类是约束信息的主要消费者：虽然该类能够处理多次输入的约束（如上述情况所需），但它坚持认为权重是完全相同的。
 * 使用这个方案，子面的自由度与父面的自由度相制约，这些自由度包含在父面的边上；有可能其中一些反过来制约自己，导致更长的制约链，AffineConstraints类最终将不得不整理出来。上面描述的约束被 DoFTools::make_hanging_node_constraints()
 * 函数所使用，该函数构造了一个AffineConstraints对象）。然而，这对于FiniteElement和派生类来说并不重要，因为它们只对一个单元和它的近邻进行局部操作，并没有看到更大的画面。
 * @ref hp_paper 详细说明了在实践中是如何处理这种链的。
 *
 *  <h4>Helper functions</h4>
 * 构建有限元和计算上述矩阵往往是一项繁琐的工作，尤其是在必须对几个维度进行计算的情况下。通过使用上面提到的中间类（如FE_Poly，FE_PolyTensor等），可以避免大部分的工作。其他任务可以通过命名空间FETools中的一些函数来实现自动化。
 * <h5>Computing the correct basis from a set of linearly independent
 * functions</h5>
 * 首先，计算任意阶数和维度的形状函数的基础可能已经很困难了。另一方面，如果给出 @ref GlossNodes "节点值"
 * ，那么节点函数和基函数之间的二元关系就定义了基。因此，形状函数空间可以从一组线性独立的函数中定义，这样，实际的有限元基是由它们的线性组合计算出来的。这些组合的系数是由节点值的对偶性决定的，并形成一个矩阵。
 * 使用这个矩阵可以分两步构建形状函数的基础。  <ol>
 * <li>
 * 使用任意基<i>w<sub>j</sub></i>定义形状函数的空间，并计算应用于这些基函数的节点函数<i>N<sub>i</sub></i>的矩阵<i>M</i>，从而使其条目为<i>m<sub>ij</sub>
 * = N<sub>i</sub>(w<sub>j</sub>)</i>。
 * <li>  计算有限元形状函数空间的基<i>v<sub>j</sub></i>，将<i>M<sup>-1</sup></i>应用于基<i>w<sub>j</sub></i>。  </ol>
 * 矩阵<i>M</i>可以用 FETools::compute_node_matrix(). 计算。这个函数依赖于#generalized_support_points和 FiniteElement::convert_generalized_support_point_values_to_dof_values() 的存在（更多信息见 @ref GlossGeneralizedSupport "关于广义支持点的术语条目"
 * ）。有了这个，我们就可以在派生自FiniteElement的类的构造函数中使用下面这段代码来计算
 * $M$ 矩阵。
 *
 * @code
 * FullMatrix<double> M(this->n_dofs_per_cell(), this->n_dofs_per_cell());
 * FETools::compute_node_matrix(M,this);
 * this->inverse_node_matrix.reinit(this->n_dofs_per_cell(),
 * this->n_dofs_per_cell()); this->inverse_node_matrix.invert(M);
 * @endcode
 * 不要忘记确保在这之前初始化#unit_support_points或#generalized_support_points!
 * <h5>Computing prolongation matrices</h5>
 * 一旦你有了形状函数，你就可以定义矩阵，将数据从一个单元转移到它的子单元，或者反过来。这当然是多网格中常见的操作，但在网格细化后将解从一个网格内插到另一个网格时也会用到，以及在一些误差估计器的定义中。
 * 为了定义延长矩阵，即那些描述有限元场从一个单元转移到其子单元的矩阵，有限元的实现可以手工填写#延长阵列，或者可以调用
 * FETools::compute_embedding_matrices(). 。
 * 在后一种情况下，所需要的是以下一段代码。
 *
 * @code
 * for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
 * this->prolongation[c].reinit (this->n_dofs_per_cell(),
 *                               this->n_dofs_per_cell());
 * FETools::compute_embedding_matrices (*this, this->prolongation);
 * @endcode
 * 在这个例子中，延长几乎都是通过嵌入来实现的，也就是说，子单元上的函数的节点值可能与父单元上的函数的节点值不同，但作为
 * $\mathbf x\in{\mathbb R}^\text{spacedim}$
 * 的函数，子单元上的有限元场与父单元上的相同。
 *
 *  <h5>Computing restriction matrices</h5>
 * 相反的操作，将定义在子单元上的有限元函数限制在父单元上，通常是通过将子单元上的有限元函数插值到父单元的节点值来实现。在deal.II中，限制操作被实现为对一个单元的子单元的循环，每个单元都对该子单元上的未知数向量应用一个矩阵（这些矩阵存储在#restriction中，可通过get_restriction_matrix()访问）。然后需要实现的操作结果是出乎意料的困难，但描述起来很有启发，因为它也定义了#restriction_is_additive_flags数组的含义（通过restriction_is_additive()函数访问）。
 * 举个具体的例子，假设我们在1d中使用一个 $Q_1$
 * 元素，在每个父细胞和子细胞上的自由度（局部和全局）编号如下。
 *
 * @code
 * meshes:
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
 *
 *
 *
 * -------*
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
 *
 *
 *
 *
 *
 * ---*---*
 * local DoF numbers:  0       1                        0  1|0  1
 * global DoF numbers: 0       1                        0   1   2
 * @endcode
 * 然后我们希望限制操作将子单元0上的第2个自由度的值作为父单元上第2个自由度的值，并将子单元1上的第1个自由度的值作为父单元上第1个自由度的值。理想情况下，我们希望这样写：@f[
 * U^\text{coarse}|_\text{parent}
 * = \sum_{\text{child}=0}^1 R_\text{child} U^\text{fine}|_\text{child}
 * @f]，其中 $U^\text{fine}|_\text{child=0}=(U^\text{fine}_0,U^\text{fine}_1)^T$ 和 $U^\text{fine}|_\text{child=1}=(U^\text{fine}_1,U^\text{fine}_2)^T$  。像这样写所要求的操作，在这里可以选择@f[
 * R_0 = \left(\begin{matrix}1 & 0 \\ 0 & 0\end{matrix}\right), \qquad\qquad
 * R_1 = \left(\begin{matrix}0 & 0 \\ 0 & 1\end{matrix}\right).
 * @f]然而，如果我们去找一个 $Q_2$
 * 元素的自由度如下，这种方法已经失败了。
 *
 * @code
 * meshes:
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
 *
 *
 *
 * -------*
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
 *
 *
 *
 *
 *
 * ----*----*
 * local DoF numbers:  0   2   1                        0 2 1|0 2 1
 * global DoF numbers: 0   2   1                        0 2  1  4 3
 * @endcode
 * 像上面那样把事情写成矩阵运算的总和不容易成功，因为我们必须把非零值加到
 * $U^\text{coarse}_2$ 中两次，每个子元素一次。
 * 因此，限制通常被实现为<i>concatenation</i>操作。也就是说，我们首先计算每个孩子的单独限制，@f[
 * \tilde U^\text{coarse}_\text{child}
 * = R_\text{child} U^\text{fine}|_\text{child},
 * @f]，然后用以下代码计算 $U^\text{coarse}|_\text{parent}$ 的值。
 *
 * @code
 * for (unsigned int child=0; child<cell->n_children(); ++child)
 * for (unsigned int i=0; i<dofs_per_cell; ++i)
 *   if (U_tilde_coarse[child][i] != 0)
 *     U_coarse_on_parent[i] = U_tilde_coarse[child][i];
 * @endcode
 * 换句话说， $\tilde
 * U^\text{coarse}_\text{child}$ 的每个非零元素<i>overwrites</i>，而不是添加到 $U^\text{coarse}|_\text{parent}$ 的相应元素。这通常也意味着来自两个不同单元的限制矩阵应该在它们都想触及的粗大自由度的数值上达成一致（否则结果将取决于我们在子代上循环的顺序，这将是不合理的，因为子代的顺序是一个本来就很随意的约定）。例如，在上面的例子中，限制矩阵将是@f[
 * R_0 = \left(\begin{matrix}1 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 1 & 0
 * \end{matrix}\right), \qquad\qquad R_1 = \left(\begin{matrix}0 & 0 & 0 \\ 0 &
 * 1 & 0 \\ 1 & 0 & 0 \end{matrix}\right),
 * @f]，兼容条件是 $R_{0,21}=R_{1,20}$ ，因为它们都表明 $U^\text{coarse}|_\text{parent,2}$ 应该被设置为 $U^\text{fine}|_\text{child=0,1}$ 和 $U^\text{fine}|_\text{child=1,0}$ 的1倍。
 * 不幸的是，并不是所有的有限元都允许以这种方式写限制操作。例如，对于片状常数FE_DGQ(0)元素，父单元上的有限元域的值不能通过子单元的内插来确定。相反，唯一合理的选择是将其作为子单元之间的<i>average</i>值
 *
 * --所以我们又回到了和的操作，而不是连接的操作。进一步的思考表明，限制是否应该是加法的，是单个形状函数的属性，而不是整个有限元的属性。因此，
 * FiniteElement::restriction_is_additive()
 * 函数返回一个特定的形状函数是否应该通过连接（返回值为
 * @p false) 或通过加法（返回值为 @p true),
 * ，然后整体操作的正确代码如下（事实上，在
 * DoFAccessor::get_interpolated_dof_values()): 中实现了这一点
 *
 * @code
 * for (unsigned int child=0; child<cell->n_children(); ++child)
 * for (unsigned int i=0; i<dofs_per_cell; ++i)
 *   if (fe.restriction_is_additive(i) == true)
 *     U_coarse_on_parent[i] += U_tilde_coarse[child][i];
 *   else
 *     if (U_tilde_coarse[child][i] != 0)
 *       U_coarse_on_parent[i] = U_tilde_coarse[child][i];
 * @endcode
 *
 *  <h5>Computing #interface_constraints</h5> 约束矩阵可以用
 * FETools::compute_face_embedding_matrices().
 * 半自动地计算，这个函数为一个面的每个子节点分别计算粗网格函数对细网格函数的表示。这些矩阵必须被卷积成一个单一的矩形约束矩阵，消除共同顶点和边缘以及粗网格顶点上的自由度。关于这个编号的细节，见上面的讨论。
 *
 *
 * @ingroup febase fe
 *
 *
 */
template <int dim, int spacedim = dim>
class FiniteElement : public Subscriptor, public FiniteElementData<dim>
{
public:
  /**
   * 图像空间的维度，对应于三角法。
   *
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * 一个基类，用于派生有限元类可能希望存储的内部数据。
   * 该类的使用方法如下。每当FEValues(或FEFaceValues或FESubfaceValues)对象被初始化时，它要求它所关联的有限元创建一个从这里的当前类派生的对象。这是通过每个派生类的
   * FiniteElement::get_data()
   * 函数完成的。然后这个对象作为一个常量对象被传递给
   * FiniteElement::fill_fe_values(),  FiniteElement::fill_fe_face_values(),
   * 和 FiniteElement::fill_fe_subface_values()
   * 函数。这些对象的意图是使有限元类可以在开始时（在调用
   * FiniteElement::get_data()
   * 函数时）预先计算一次信息，然后可以在随后访问的每个单元上使用。这方面的一个例子是参考单元正交点的形状函数值，无论访问哪个单元都是一样的，因此可以在开始时计算一次并在以后重复使用。
   * 因为只有派生类才能知道他们可以预先计算什么，所以每个派生类如果想存储在开始时计算过一次的信息，需要从这个类派生出自己的InternalData类，并通过其get_data()函数返回派生类型的对象。
   *
   */
  class InternalDataBase
  {
  public:
    /**
     * 构造函数。设置update_flags为 @p update_default ， @p
     * first_cell 为 @p true. 。
     *
     */
    InternalDataBase();

    /**
     * 解构器。做成虚拟的，以允许多态性。
     *
     */
    virtual ~InternalDataBase() = default;

    /**
     * 禁止复制构造。
     *
     */
    InternalDataBase(const InternalDataBase &) = delete;

    /**
     * 一组更新标志，指定FiniteElement接口的实现需要在每个单元或面计算的信息种类，即在
     * FiniteElement::fill_fe_values() 和朋友中。        这组标志被
     * FiniteElement::get_data(),  FiniteElement::get_face_data(), 或
     * FiniteElement::get_subface_data(),
     * 的实现保存在这里，是传递给那些需要对每个单元进行重新计算的函数的更新标志的子集。(对应于在调用
     * FiniteElement::get_data()
     * 时已经可以一次性计算的信息的标志子集。
     *
     * - 或该接口的实现
     *
     * - 不需要存储在这里，因为它已经被处理过了)。
     *
     */
    UpdateFlags update_each;

    /**
     * 返回这个对象的内存消耗估计值（以字节为单位）。
     *
     */
    virtual std::size_t
    memory_consumption() const;
  };

public:
  /**
   * 构造函数：初始化这个基类的所有有限元素的字段。      @param[in]  fe_data 一个存储关于要构造的元素的识别（通常是积分）信息的对象。特别是，这个对象将包含诸如每个单元（以及每个顶点、线等）的自由度数量、矢量分量的数量等数据。  这个参数被用来初始化当前正在构建的对象的基类。    @param[in]  restriction_is_additive_flags 一个大小为 <code>dofs_per_cell</code> 的向量（或大小为1，见下文），对于每个形状函数来说，它说明该形状函数是否是加性的。这些标志的含义在这个类的一般文档中关于限制矩阵的部分有描述。    @param[in]  nonzero_components 一个大小为 <code>dofs_per_cell</code> （或大小为1，见下文）的向量，为每个形状函数提供一个ComponentMask（大小为 <code>fe_data.n_components()</code> ），指示这个形状函数在哪些向量分量中是非零的（在将形状函数映射到实数单元之后）。对于 "原始 "
   * 形状函数，这个分量掩码将有一个条目（关于原始元素的更多信息见
   * @ref GlossPrimitive
   * ）。另一方面，对于诸如Raviart-Thomas或Nedelec元素，形状函数在一个以上的向量分量中是非零的（在映射到实数单元之后），给定的分量掩码将包含一个以上的条目。(对于这两个元素，事实上所有的条目都将被设置，但是如果你将FE_RaviartThomas和FE_Nedelec一起耦合到一个FES系统中，情况就不是这样了。)
   * @pre   <code>restriction_is_additive_flags.size() == dofs_per_cell</code>
   * ，或  <code>restriction_is_additive_flags.size() == 1</code>
   * 。在后一种情况下，数组被简单地解释为大小为
   * <code>dofs_per_cell</code>
   * ，其中每个元素的值与给出的单个元素相同。      @pre
   * <code>nonzero_components.size() == dofs_per_cell</code>  , 或
   * <code>nonzero_components.size() == 1</code>
   * 。在后一种情况下，数组被简单地解释为大小为
   * <code>dofs_per_cell</code>
   * ，其中每个元素等于给定的单个元素中提供的组件掩码。
   *
   */
  FiniteElement(const FiniteElementData<dim> &    fe_data,
                const std::vector<bool> &         restriction_is_additive_flags,
                const std::vector<ComponentMask> &nonzero_components);

  /**
   * 移动构造器。
   *
   */
  FiniteElement(FiniteElement<dim, spacedim> &&) = default; // NOLINT

  /**
   * 复制构造函数。
   *
   */
  FiniteElement(const FiniteElement<dim, spacedim> &) = default;

  /**
   * 虚拟解构器。确保这个类的指针被正确删除。
   *
   */
  virtual ~FiniteElement() override = default;

  /**
   * 创建以该类为基础元素的FESystem的信息，并且具有多重性
   * @p multiplicity.
   * 特别是，这个函数的返回类型可以用于FESystem对象的构造函数中。
   * 这个函数调用clone()，因此创建了一个当前对象的副本。
   *
   */
  std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>
  operator^(const unsigned int multiplicity) const;

  /**
   * 一种虚拟的拷贝构造函数，该函数返回有限元对象的拷贝。派生类需要在这个基类中覆盖这里的函数，并返回一个与派生类相同类型的对象。
   * 库中的一些地方，例如FESystem以及 hp::FECollection
   * 类的构造函数，需要在不知道其确切类型的情况下对有限元进行复制。他们通过这个函数来做。
   *
   */
  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const = 0;

  /**
   * 返回一个唯一标识有限元的字符串。一般的惯例是，这是类的名称，然后是角括号中的尺寸，以及括号中的多项式程度和其他必要的东西。例如，<tt>FE_Q<2>(3)</tt>是2d中立方体元素的返回值。
   * 元素系统有自己的命名规则，见FESystem类。
   *
   */
  virtual std::string
  get_name() const = 0;

  /**
   * 如果给定的参数等于零，该操作符返回对当前对象的引用。虽然这看起来不是特别有用，但在编写与::DoFHandler和hp-版本
   * hp::DoFHandler,
   * 一起工作的代码时，它是有帮助的，因为这样就可以写出这样的代码。
   * @code
   * dofs_per_cell =
   * dof_handler->get_fe()[cell->active_fe_index()].n_dofs_per_cell();
   * @endcode
   * 这段代码在这两种情况下都不能用，因为
   * DoFHandler::get_fe() 返回一个有限元，而 hp::DoFHandler::get_fe()
   * 返回一个有限元的集合，不提供 <code>dofs_per_cell</code>
   * 成员变量：人们首先要选择对哪个有限元进行处理，这是用操作符[]来完成。幸运的是，
   * <code>cell-@>active_fe_index()</code>
   * 也适用于非hp类，在这种情况下只需返回0。本操作符[]接受这个零参数，返回其集合中索引为零的有限元（当然，无论如何，它只由本有限元组成）。
   * @deprecated 随着 DoFHandler::get_fe(int) 和 hp::DoFHandler
   * 类的废弃，这个操作符不再有任何用途。
   *
   */
  DEAL_II_DEPRECATED const FiniteElement<dim, spacedim> &
                           operator[](const unsigned int fe_index) const;

  /**
   * @name  形状函数访问 
     * @{ 
   *
   */

  /**
   * 返回 @p ith 形状函数在 @p p.  @p p
   * 点的值，该点是参考元素上的一个点。如果有限元是矢量值的，那么返回这个形状函数的矢量值的唯一非零分量的值。如果形状函数有一个以上的非零分量（我们用非原始的术语来指代），那么实现这个函数的派生类应该抛出一个ExcShapeFunctionNotPrimitive类型的异常。在这种情况下，请使用
   * shape_value_component() 函数。
   * 如果所考虑的FiniteElement的形状函数取决于实空间中的单元的形状，也就是说，如果形状函数不是通过从参考单元的映射来定义的，那么这个函数的实现应该抛出一个ExcUnitShapeValuesDoNotExist类型的异常。一些不符合要求的元素就是这样定义的，FE_DGPNonparametric类就是如此，这只是一个例子。
   * 这个虚拟函数的默认实现正是这样做的，也就是说，它只是抛出一个ExcUnitShapeValuesDoNotExist类型的异常。
   *
   */
  virtual double
  shape_value(const unsigned int i, const Point<dim> &p) const;

  /**
   * 就像shape_value()一样，但是当形状函数有一个以上的非零向量分量时，这个函数将被调用。在这种情况下，这个函数应该返回
   * @p component-th 形状函数在 @p p. 点的 @p ith 向量分量的值。
   *
   */
  virtual double
  shape_value_component(const unsigned int i,
                        const Point<dim> & p,
                        const unsigned int component) const;

  /**
   * 返回 @p ith 形状函数在点 @p p. 的梯度 @p p
   * 是参考元素上的一个点，同样，梯度是单元格上关于单元格坐标的梯度。如果有限元是矢量值的，那么返回这个形状函数的矢量值的唯一非零分量的值。如果形状函数有一个以上的非零分量（我们用非原始的术语来指代），那么实现这个函数的派生类应该抛出一个ExcShapeFunctionNotPrimitive类型的异常。在这种情况下，请使用
   * shape_grad_component() 函数。
   * 如果所考虑的FiniteElement的形状函数取决于实空间中的单元的形状，也就是说，如果形状函数不是通过从参考单元的映射来定义的，那么这个函数的实现应该抛出一个ExcUnitShapeValuesDoNotExist类型的异常。一些不符合要求的元素就是这样定义的，FE_DGPNonparametric类就是如此，这只是一个例子。
   * 这个虚拟函数的默认实现正是这样做的，也就是说，它只是抛出一个ExcUnitShapeValuesDoNotExist类型的异常。
   *
   */
  virtual Tensor<1, dim>
  shape_grad(const unsigned int i, const Point<dim> &p) const;

  /**
   * 就像shape_grad()一样，但是当形状函数有一个以上的非零向量分量时，这个函数将被调用。在这种情况下，这个函数应该返回
   * @p component-th 形状函数的 @p ith 向量分量在 @p p.
   * 点的梯度。
   *
   */
  virtual Tensor<1, dim>
  shape_grad_component(const unsigned int i,
                       const Point<dim> & p,
                       const unsigned int component) const;

  /**
   * 返回 @p ith 形状函数在单元格上 @p p
   * 点的二次导数的张量。该导数是单元格上相对于单元格坐标的导数。如果有限元是矢量值的，那么返回这个形状函数的矢量值的唯一非零分量的值。如果形状函数有一个以上的非零分量（我们用非基元一词来指代），那么实现这个函数的派生类应该抛出一个ExcShapeFunctionNotPrimitive类型的异常。在这种情况下，请使用
   * shape_grad_grad_component() 函数。
   * 如果所考虑的FiniteElement的形状函数取决于实空间中的单元的形状，也就是说，如果形状函数不是通过从参考单元的映射来定义的，那么这个函数的实现应该抛出一个ExcUnitShapeValuesDoNotExist类型的异常。一些不符合要求的元素就是这样定义的，FE_DGPNonparametric类就是如此，这只是一个例子。
   * 这个虚拟函数的默认实现正是这样做的，也就是说，它只是抛出一个ExcUnitShapeValuesDoNotExist类型的异常。
   *
   */
  virtual Tensor<2, dim>
  shape_grad_grad(const unsigned int i, const Point<dim> &p) const;

  /**
   * 就像shape_grad_grad()一样，但是当形状函数有一个以上的非零向量分量时，这个函数将被调用。在这种情况下，这个函数应该返回
   * @p component-th 形状函数的 @p ith 向量分量在 @p p.
   * 点的梯度。
   *
   */
  virtual Tensor<2, dim>
  shape_grad_grad_component(const unsigned int i,
                            const Point<dim> & p,
                            const unsigned int component) const;

  /**
   * 返回 @p ith 形状函数在单元格上 @p p
   * 点的三阶导数的张量。该导数是单元格上相对于单元格坐标的导数。如果有限元是矢量值的，那么返回这个形状函数的矢量值的唯一非零分量的值。如果形状函数有一个以上的非零分量（我们用非基元一词来指代），那么实现这个函数的派生类应该抛出一个ExcShapeFunctionNotPrimitive类型的异常。在这种情况下，请使用
   * shape_3rd_derivative_component() 函数。
   * 如果所考虑的FiniteElement的形状函数取决于实空间中的单元的形状，也就是说，如果形状函数不是通过从参考单元的映射来定义的，那么这个函数的实现应该抛出一个ExcUnitShapeValuesDoNotExist类型的异常。一些不符合要求的元素就是这样定义的，FE_DGPNonparametric类就是如此，这只是一个例子。
   * 这个虚拟函数的默认实现正是这样做的，也就是说，它只是抛出一个ExcUnitShapeValuesDoNotExist类型的异常。
   *
   */
  virtual Tensor<3, dim>
  shape_3rd_derivative(const unsigned int i, const Point<dim> &p) const;

  /**
   * 就像shape_3rd_derivative()一样，但这个函数将在形状函数有一个以上的非零矢量分量时被调用。在这种情况下，这个函数应该返回
   * @p component- 形状函数的第 @p ith 个向量分量在点 @p p.
   * 的梯度。
   *
   */
  virtual Tensor<3, dim>
  shape_3rd_derivative_component(const unsigned int i,
                                 const Point<dim> & p,
                                 const unsigned int component) const;

  /**
   * 返回 @p ith 形状函数在单元格上 @p p
   * 点的第四导数的张量。该导数是单元格上相对于单元格坐标的导数。如果有限元是矢量值的，那么返回这个形状函数的矢量值的唯一非零分量的值。如果形状函数有一个以上的非零分量（我们用非基元一词来指代），那么实现这个函数的派生类应该抛出一个ExcShapeFunctionNotPrimitive类型的异常。在这种情况下，请使用
   * shape_4th_derivative_component() 函数。
   * 如果所考虑的FiniteElement的形状函数取决于实空间中的单元的形状，也就是说，如果形状函数不是通过从参考单元的映射来定义的，那么这个函数的实现应该抛出一个ExcUnitShapeValuesDoNotExist类型的异常。一些不符合要求的元素就是这样定义的，FE_DGPNonparametric类就是如此，这只是一个例子。
   * 这个虚拟函数的默认实现正是这样做的，也就是说，它只是抛出一个ExcUnitShapeValuesDoNotExist类型的异常。
   *
   */
  virtual Tensor<4, dim>
  shape_4th_derivative(const unsigned int i, const Point<dim> &p) const;

  /**
   * 就像shape_4th_derivative()一样，但这个函数将在形状函数有一个以上的非零向量分量时被调用。在这种情况下，这个函数应该返回
   * @p component- 形状函数的第 @p ith 个向量分量在 @p p.
   * 点的梯度。
   *
   */
  virtual Tensor<4, dim>
  shape_4th_derivative_component(const unsigned int i,
                                 const Point<dim> & p,
                                 const unsigned int component) const;
  /**
   * 如果形状函数 @p shape_index 在面 @p face_index.
   * 的某处有非零函数值，该函数通常用于确定由面积分产生的一些矩阵元素是否可以被假定为零，因此可以从积分中省略。
   * 在这个基类中提供了一个默认的实现，它总是返回  @p
   * true.  这是一个安全的方法。
   *
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const;

  //@}
  /**
   * @name  转移和约束矩阵  @{
   *
   */

  /**
   * 返回描述从给定的 @p child （由给定的 @p refinement_case)
   * 获得的有限元场到父单元的限制矩阵。返回的矩阵的解释取决于limittion_is_additive()为每个形状函数返回的内容。
   * 行和列指数分别与粗网格和细网格空间相关，与相关运算符的定义一致。
   * 如果投影矩阵没有在派生的有限元类中实现，这个函数会以
   * FiniteElement::ExcProjectionVoid.
   * 类型的异常中止，你可以通过首先调用restriction_is_implemented()或isotropic_restriction_is_implemented()函数检查是否会发生这种情况。
   *
   */
  virtual const FullMatrix<double> &
  get_restriction_matrix(const unsigned int         child,
                         const RefinementCase<dim> &refinement_case =
                           RefinementCase<dim>::isotropic_refinement) const;

  /**
   * 网格间的延长/嵌入矩阵。
   * 从粗网格空间到细网格空间（这两个空间都被识别为定义在父单元和子单元上的函数）的身份运算符与一个矩阵
   * @p P
   * 相关联，该矩阵以其节点值映射这些函数的相应表示。这里返回该矩阵
   * @p P_i 对单个子单元的限制。    矩阵 @p P
   * 是串联的，而不是单元格矩阵 @p
   * P_i的总和。也就是说，如果同一个非零条目<tt>j,k</tt>存在于两个不同的子矩阵
   * @p P_i,
   * 中，该值在两个矩阵中应该是相同的，它只被复制到矩阵
   * @p P 一次。
   * 行和列指数分别与细格和粗格空间相关，与相关运算符的定义一致。
   * 这些矩阵被组装多层次方法的延长矩阵的程序所使用。
   * 在使用这个矩阵阵列组装单元间的转移矩阵时，延长矩阵中的零元素被丢弃，不会填满转移矩阵。
   * 如果延长矩阵没有在派生的有限元类中实现，这个函数将以
   * FiniteElement::ExcEmbeddingVoid.
   * 类型的异常中止。你可以通过首先调用prolongation_is_implemented()或isotropic_prolongation_is_implemented()函数检查是否会发生这种情况。
   *
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix(const unsigned int         child,
                          const RefinementCase<dim> &refinement_case =
                            RefinementCase<dim>::isotropic_refinement) const;

  /**
   * 返回这个元素是否实现了它的延长矩阵。该返回值还表明调用get_prolongation_matrix()函数是否会产生错误。
   * 请注意，只有在各向同性和所有各向异性的细化情况下，该函数才会返回
   * <code>true</code>
   * 的延长矩阵。如果你只对各向同性细化的延长矩阵感兴趣，请使用
   * isotropic_prolongation_is_implemented 函数来代替。
   * 这个函数在这里主要是为了让我们编写更有效的测试程序，我们在各种奇怪的元素上运行，对于这些测试，我们只需要排除某些测试，以防某些东西没有实现。在一般情况下，它在应用中可能不会有很大的帮助，因为如果需要这些功能而它们没有被实现，那就没有什么办法了。这个函数可以用来检查对<tt>get_prolongation_matrix()</tt>的调用是否会成功；然而，人们仍然需要应对这个函数所表达的信息的缺乏。
   *
   */
  bool
  prolongation_is_implemented() const;

  /**
   * 返回这个元素是否为各向同性的孩子实现了其延长矩阵。该返回值也表明调用
   * @p get_prolongation_matrix 函数是否会产生错误。
   * 这个函数在这里主要是为了让我们编写更有效的测试程序，我们在各种奇怪的元素上运行，对于这些元素，我们只需要排除某些测试，以防某些东西没有实现。在一般情况下，它在应用中可能不会有很大的帮助，因为如果需要这些功能而它们没有被实现，那就没有什么办法了。这个函数可以用来检查对<tt>get_prolongation_matrix()</tt>的调用是否会成功；然而，人们仍然需要应对这个函数所表达的信息的缺乏。
   *
   */
  bool
  isotropic_prolongation_is_implemented() const;

  /**
   * 返回这个元素是否实现了它的限制矩阵。该返回值也表明调用get_restriction_matrix()函数是否会产生错误。
   * 注意，只有在各向同性和所有各向异性细化情况下的限制矩阵被实现时，该函数才会返回
   * <code>true</code>
   * 。如果你只对各向同性细化的限制矩阵感兴趣，请使用
   * isotropic_restriction_is_implemented() 函数来代替。
   * 这个函数在这里主要是为了让我们编写更有效的测试程序，我们在各种奇怪的元素上运行，对于这些测试，我们只需要排除某些测试，以防某些东西没有被实现。在一般情况下，它在应用中可能不会有很大的帮助，因为如果需要这些功能而它们没有被实现，那就没有什么办法了。这个函数可以用来检查对<tt>get_restriction_matrix()</tt>的调用是否会成功；然而，人们仍然需要应对这个函数所表达的信息的缺乏。
   *
   */
  bool
  restriction_is_implemented() const;

  /**
   * 返回这个元素是否实现了其各向同性的子女的限制矩阵。该返回值还表明调用get_restriction_matrix()函数是否会产生错误。
   * 这个函数在这里主要是为了让我们编写更有效的测试程序，我们在各种奇怪的元素上运行，对于这些元素，我们只需要排除某些测试，以防某些东西没有实现。在一般情况下，它在应用中可能不会有很大的帮助，因为如果需要这些功能而它们没有被实现，那就没有什么办法了。这个函数可以用来检查对<tt>get_restriction_matrix()</tt>的调用是否会成功；然而，人们仍然需要应对这个函数所表达的信息的缺乏。
   *
   */
  bool
  isotropic_restriction_is_implemented() const;


  /**
   * 访问#restriction_is_additive_flags字段。更多信息请参见一般类文档中关于限制矩阵的讨论。
   * 索引必须在零和这个元素的形状函数的数量之间。
   *
   */
  bool
  restriction_is_additive(const unsigned int index) const;

  /**
   * 返回一个对矩阵的只读引用，该矩阵描述了精炼单元和非精炼单元之间界面上的约束。
   * 一些有限元不（还）实现悬挂节点约束。如果是这种情况，那么这个函数将产生一个异常，因为不能产生有用的返回值。如果你应该有办法接受这种情况，那么你可能想使用constraints_are_implemented()函数来预先检查这个函数是否会成功或产生异常。
   *
   */
  const FullMatrix<double> &
  constraints(const dealii::internal::SubfaceCase<dim> &subface_case =
                dealii::internal::SubfaceCase<dim>::case_isotropic) const;

  /**
   * 返回这个元素是否实现了它的悬挂节点约束。该返回值也表明对约束()函数的调用是否会产生一个错误。
   * 这个函数在这里主要是为了让我们编写更有效的测试程序，我们运行在各种奇怪的元素上，对于这些元素，我们只需要排除某些测试，以防悬挂节点约束没有实现。一般来说，它在应用中可能不会有很大的帮助，因为如果需要悬挂节点约束而它们没有被实现的话，就没有什么可以做的。这个函数可以用来检查对<tt>constraints()</tt>的调用是否会成功；然而，人们仍然需要应对这个函数所表达的信息的缺乏。
   *
   */
  bool
  constraints_are_implemented(
    const dealii::internal::SubfaceCase<dim> &subface_case =
      dealii::internal::SubfaceCase<dim>::case_isotropic) const;


  /**
   * 返回这个元素是否以新的方式实现了它的悬挂节点约束，这必须用于使元素
   * "hp-compatible"。
   * 这意味着，该元素正确实现了get_face_interpolation_matrix和get_subface_interpolation_matrix方法。因此，返回值也表明对get_face_interpolation_matrix()方法和get_subface_interpolation_matrix()方法的调用是否会产生一个错误。
   * 目前这个函数的主要目的是让make_hanging_node_constraints方法决定是否可以使用新的程序，这些程序应该在hp-framework中工作，或者是否应该使用旧的经过充分验证但不具备hp能力的函数。
   * 一旦过渡到计算接口约束的新方案，这个函数将是多余的，可能会消失。
   * 派生类应该相应地实现这个函数。默认的假设是有限元不提供具有hp能力的面插值，因此默认的实现会返回
   * @p false.  。
   *
   */
  virtual bool
  hp_constraints_are_implemented() const;


  /**
   * 返回从给定的有限元到现在的有限元的内插矩阵。矩阵的大小是#dofs_per_cell乘以<tt>source.#dofs_per_cell</tt>。
   * 衍生元素将不得不实现这个函数。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现给定元素的插值，那么他们必须抛出一个ExcInterpolationNotImplemented类型的异常。
   *
   */
  virtual void
  get_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                           FullMatrix<double> &                matrix) const;
  //@}

  /**
   * @name 支持hp的函数  
     * @{ 
   *
   */


  /**
   * 返回从一个元素的面插值到相邻元素的面的矩阵。
   * 矩阵的大小是<tt>source.#dofs_per_face</tt>乘以<tt>this->#dofs_per_face</tt>。
   * 衍生元素将不得不实现这个函数。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现给定元素的插值，那么他们必须抛出一个ExcInterpolationNotImplemented类型的异常。
   *
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double> &                matrix,
                                const unsigned int face_no = 0) const;


  /**
   * 返回从一个元素的面内插到邻近元素的子面的矩阵。
   * 矩阵的大小是<tt>source.#dofs_per_face</tt>乘以<tt>this->#dofs_per_face</tt>。
   * 衍生元素将不得不实现这个函数。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现给定元素的插值，那么他们必须抛出一个ExcInterpolationNotImplemented类型的异常。
   *
   */
  virtual void
  get_subface_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                   const unsigned int                  subface,
                                   FullMatrix<double> &                matrix,
                                   const unsigned int face_no = 0) const;
  //@}


  /**
   * @name  支持HP的函数-  
     * @{ 
   *
   */

  /**
   * 如果在一个顶点上，有几个有限元处于活动状态，hp代码首先为这些FEs的每个自由度分配不同的全局索引。然后调用这个函数来找出其中哪些应该得到相同的值，从而可以得到相同的全局自由度指数。
   * 因此，该函数返回当前有限元对象的自由度与 @p fe_other,
   * 的自由度之间的相同性列表，后者是对代表在该特定顶点上活动的其他有限元之一的有限元对象的引用。该函数计算两个有限元对象的哪些自由度是相等的，这两个自由度的编号都在零和两个有限元的n_dofs_per_vertex()的相应值之间。每一对的第一个索引表示本元素的一个顶点自由度，而第二个是另一个有限元素的相应索引。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(const FiniteElement<dim, spacedim> &fe_other) const;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理线上自由度。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(const FiniteElement<dim, spacedim> &fe_other) const;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理四边形上的自由度。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int                  face_no = 0) const;

  /**
   * 返回这个元素是否支配另一个作为参数的元素  @p fe_other,  是否相反，是否两者都不支配，或者是否两者都可以支配。 @p codim 参数描述了被调查的子空间的二维度，并指定其受此比较。例如，如果`codim==0`，那么这个函数就会比较哪个元素在细胞水平上占优势。如果`codim==1`，那么元素在面进行比较，也就是说，比较发生在限制在面的两个有限元素的函数空间之间。较大的`codim'值也相应地起作用。    关于支配的定义，见 FiniteElementDomination::Domination ，特别是 @ref hp_paper "hp-论文"
   * 。
   *
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int                  codim = 0) const;

  //@}

  /**
   * 比较运算符。
   * 在当前类中的实现是检查当前对象和作为参数给定的对象之间的以下信息是否相等，按照这个顺序。
   *
   *
   *
   *
   *
   * - 当前对象和给定对象的动态类型（即最派生类的类型）。
   *
   *
   *
   *
   *
   *
   * - 由get_name()返回的名称。
   *
   *
   *
   *
   *
   * - 如同FiniteElementData中的所有字段。
   *
   *
   *
   *
   * - 约束矩阵。    这涵盖了大多数元素可以不同的情况，但也有派生元素不同的情况，对于这些情况，当前函数仍然返回  @p true.  对于这些情况，派生类应该重载这个函数。
   * @note  这个操作符特别不检查当前类的下列成员变量。
   *
   *
   *
   *
   *
   *
   * - 限制矩阵。
   *
   *
   * - 此对象和参数的延长矩阵。   这是因为这些成员变量可能只在派生类的要求下被初始化，而不是立即可用。   因此，比较这些成员不仅成本高，因为这些一般都是大数组，而且其计算也可能很昂贵。另一方面，派生类的这些数组对于两个对象来说可能是不同的，即使上面的列表比较起来是相等的，也可能想要实现自己的operator==()。
   *
   */
  virtual bool
  operator==(const FiniteElement<dim, spacedim> &fe) const;

  /**
   * 非等价比较运算符。以平等比较运算符的方式定义。
   *
   */
  bool
  operator!=(const FiniteElement<dim, spacedim> &) const;

  /**
   * @name 索引计算  @{
   *
   */
  /**
   * 从这个有限元内的形状函数的索引中计算出这个形状函数对应的向量分量和索引。    如果元素是标量，那么分量总是零，这个分量内的索引等于总索引。    如果引用的形状函数有一个以上的非零分量，那么它不能与一个矢量分量相关联，并且将引发一个ExcShapeFunctionNotPrimitive类型的异常。    注意，如果元素是由其他（基）元素组成的，并且一个基元素有一个以上的分量，但是它的所有形状函数都是原始的（即在只有一个分量中是非零的），那么这个映射包含有效的信息。然而，这个元素的形状函数在一个分量中的索引（即这个数组中各自条目的第二个数字）并不表示各自的形状函数在基数元素中的索引（因为它有一个以上的向量分量）。关于这个信息，请参考#system_to_base_table字段和system_to_base_index()函数。    关于这个函数的典型使用方法，见上面的类描述。    在 step-8 和 @ref step_20 "  step-20 "
   * 教程程序以及 @ref vector_valued
   * 模块中，对该函数的使用进行了广泛的解释。
   *
   */
  std::pair<unsigned int, unsigned int>
  system_to_component_index(const unsigned int index) const;

  /**
   * 计算给定矢量元件和索引的形状函数。
   * 如果该元素是标量，那么该分量必须为零，并且该分量内的索引等于总索引。
   * 这是与system_to_component_index()函数相反的操作。
   *
   */
  unsigned int
  component_to_system_index(const unsigned int component,
                            const unsigned int index) const;

  /**
   * 与system_to_component_index()函数相同，但对形状函数和它们在面上的索引进行操作。因此允许的索引范围是0...#dofs_per_face。
   * 你在应用程序中很少需要这个函数，因为几乎所有的应用程序代码都只需要处理单元格的索引，而不是面的索引。这个函数主要是为了在库内使用。
   *
   */
  std::pair<unsigned int, unsigned int>
  face_system_to_component_index(const unsigned int index,
                                 const unsigned int face_no = 0) const;

  /**
   * 对于在3D中具有非标准面孔方向的面孔，面孔（四边形）上的道夫必须被移位，以便与正确的形状函数相结合。给出一个四边形上的局部dof
   * @p index
   * ，如果该面具有非标准的面朝向、面朝上或面朝下的旋转，则返回局部索引。在二维和一维中，不需要
   * permutation，因此会抛出一个异常。
   *
   */
  unsigned int
  adjust_quad_dof_index_for_face_orientation(const unsigned int index,
                                             const unsigned int face_no,
                                             const bool face_orientation,
                                             const bool face_flip,
                                             const bool face_rotation) const;

  /**
   * 给定一个面的指数的自然排序，返回单元格上相同自由度的指数。
   * 为了解释这个概念，考虑这样的情况：我们想知道一个面的自由度，例如作为FESystem元素的一部分，是否是原始的。不幸的是，FiniteElement类中的is_primitive()函数需要一个单元格索引，所以我们需要找到对应于当前面的索引的形状函数的单元格索引。
   * 这个函数就是这样做的。
   * 实现这一点的代码将看起来像这样。
   * @code
   * for (i=0; i<dofs_per_face; ++i)
   * if (fe.is_primitive(fe.face_to_cell_index(i, some_face_no)))
   * ... do whatever
   * @endcode
   * 这个函数需要额外的参数，这些参数考虑到实际的面可以是相对于所考虑的单元格的标准排序，或者可以是翻转的，定向的，等等。
   * @param  face_dof_index
   * 一个面的自由度的索引。这个指数必须在零和每个面的自由度之间。
   * @param  face
   * 这个自由度所在的面的编号。这个数字必须在零和
   * GeometryInfo::faces_per_cell. 之间  @param  face_orientation
   * 描述面的方向的一部分。见  @ref GlossFaceOrientation  。
   * @param  face_flip 对脸部方向的描述的一部分。参见  @ref
   * GlossFaceOrientation  。    @param  face_rotation
   * 描述脸部方向的一部分。参见  @ref GlossFaceOrientation  。
   * @return
   * 这个自由度在整个单元上的自由度集合中的索引。返回值将介于0和dofs_per_cell之间。
   * @note
   * 这个函数存在于这个类中，因为这是它首次实现的地方。然而，在不知道我们有什么元素的情况下，它不能真正在最一般的情况下工作。原因是当一个面被翻转或旋转时，我们还需要知道我们是否需要交换这个面上的自由度，或者它们是否可以免于这样。
   * 对于这一点，请考虑2d中 $Q_3$
   * 元素的情况。如果face_flip是真的，那么我们就需要以相反的顺序考虑边缘上的两个自由度。另一方面，如果这个元素是
   * $Q_1^2$
   * ，那么由于这个边上的两个自由度属于不同的向量分量，它们不应该被反向考虑。所有这些表明，如果每条线或四边有一个以上的自由度，这个函数就不能工作，在这些情况下，函数会抛出一个异常，指出这个功能需要由一个知道自由度实际代表的派生类提供。
   *
   */
  virtual unsigned int
  face_to_cell_index(const unsigned int face_dof_index,
                     const unsigned int face,
                     const bool         face_orientation = true,
                     const bool         face_flip        = false,
                     const bool         face_rotation    = false) const;

  /**
   * 对于在三维中具有非标准线方向的线条，线条上的自由度必须被替换，以便与正确的形状函数相结合。给出一条线上的局部道夫
   * @p index
   * ，如果该线有非标准的线方向，则返回局部索引。在二维和一维中，不需要进行置换，所以直接返回给定的索引。
   *
   */
  unsigned int
  adjust_line_dof_index_for_line_orientation(const unsigned int index,
                                             const bool line_orientation) const;

  /**
   * 返回该有限元的哪些向量分量中 @p
   * 第i个形状函数是非零的。返回数组的长度等于这个元素的向量分量的数量。
   * 对于大多数有限元空间，这个函数的结果将是一个恰好有一个元素是
   * @p true,
   * 的向量，因为对于大多数空间，各个向量分量是独立的。在这种情况下，带有单个0的分量也是system_to_component_index()返回的第一个元素。
   * 只有对于那些耦合分量的空间，例如为了使形状函数无发散，才会有一个以上的
   * @p true 条目。 这种情况下的元素被称为非原始元素（见
   * @ref GlossPrimitive  ）。
   *
   */
  const ComponentMask &
  get_nonzero_components(const unsigned int i) const;

  /**
   * 返回 @p ith
   * 形状函数在多少个向量分量中是非零的。这个值等于get_nonzero_components()函数结果中等于
   * @p true 的条目数量。
   * 对于大多数有限元空间来说，该结果将等于1。只有对于那些矢量值形状函数将各个分量结合在一起的解析空间，例如，为了使其无发散，它才不等于1。
   *
   */
  unsigned int
  n_nonzero_components(const unsigned int i) const;

  /**
   * 返回整个有限元是否是原始的，即其所有形状函数都是原始的。如果有限元是标量的，那么情况总是这样的。
   * 由于这是一个极其常见的操作，其结果被缓存并由该函数返回。
   *
   */
  bool
  is_primitive() const;

  /**
   * 返回 @p ith
   * 形状函数是否是原始的，即形状函数只在一个向量分量中非零。那么，非原始的形状函数将是那些无发散的Ansatz空间的形状函数，其中各个向量分量是耦合的。
   * 当且仅当<tt>n_nonzero_components(i)</tt>的结果等于1时，该函数的结果为
   * @p true 。
   *
   */
  bool
  is_primitive(const unsigned int i) const;

  /**
   * 混合离散化中的基元数量。
   * 请注意，即使是矢量值的有限元，组件的数量也不需要与基元的数量相一致，因为它们可能被重复使用。例如，如果你通过使用接受一个有限元和一个倍数的构造函数来创建一个具有三个相同的有限元类的FES系统，那么基元的数量仍然是一个，尽管有限元的分量数量等于倍数。
   *
   */
  unsigned int
  n_base_elements() const;

  /**
   * 对基元对象的访问。如果元素是原子性的，那么
   * <code>base_element(0)</code> 就是 @p this. 。
   *
   */
  virtual const FiniteElement<dim, spacedim> &
  base_element(const unsigned int index) const;

  /**
   * 这个索引表示基础元素 @p index
   * 在一个组成元素中使用的频率。如果该元素是原子性的，那么结果总是等于1。更多细节见n_base_elements()函数的文档。
   *
   */
  unsigned int
  element_multiplicity(const unsigned int index) const;

  /**
   * 返回一个包含的有限元素的引用，该元素与给定的ComponentMask所选择的组件相匹配
   * @p mask.
   * 对于一个任意嵌套的FESystem，该函数返回与给定掩码相匹配的最内层的FiniteElement。如果
   * @p mask
   * 不完全匹配其中一个包含的有限元，则该方法失败。如果当前对象是一个FESystem，该方法是最有用的，因为在所有其他情况下，返回值只能是
   * @p this 。
   * 请注意，如果掩码与之匹配，但不与任何包含的对象匹配，则返回的对象可以是一个FESystem。
   * 让我们用一个有7个组件的FESystem  @p fe
   * 来说明这个函数。
   * @code
   * FESystem<2> fe_velocity(FE_Q<2>(2), 2);
   * FE_Q<2> fe_pressure(1);
   * FE_DGP<2> fe_dg(0);
   * FE_BDM<2> fe_nonprim(1);
   * FESystem<2> fe(fe_velocity, 1, fe_pressure, 1, fe_dg, 2, fe_nonprim, 1);
   * @endcode
   * 下表列出了你可以使用的所有可能的组件掩码。
   * <table> <tr> <th>ComponentMask</th> <th>Result</th> <th>Description</th>
   * </tr> <tr> <td><code>[true,true,true,true,true,true,true]</code></td>
   * <td><code>FESystem<2>[FESystem<2>[FE_Q<2>(2)^2]-FE_Q<2>(1)-FE_DGP<2>(0)^2-FE_BDM<2>(1)]</code></td>
   * <td>@p fe itself, the whole @p FESystem</td> </tr> <tr>
   * <td><code>[true,true,false,false,false,false,false]</code></td>
   * <td><code>FESystem<2>[FE_Q<2>(2)^2]</code></td> <td>just the @p
   * fe_velocity</td> </tr> <tr>
   * <td><code>[true,false,false,false,false,false,false]</code></td>
   * <td><code>FE_Q<2>(2)</code></td> <td>The first component in @p
   * fe_velocity</td> </tr> <tr>
   * <td><code>[false,true,false,false,false,false,false]</code></td>
   * <td><code>FE_Q<2>(2)</code></td> <td>The second component in @p
   * fe_velocity</td> </tr> <tr>
   * <td><code>[false,false,true,false,false,false,false]</code></td>
   * <td><code>FE_Q<2>(1)</code></td> <td>@p fe_pressure</td> </tr> <tr>
   * <td><code>[false,false,false,true,false,false,false]</code></td>
   * <td><code>FE_DGP<2>(0)</code></td> <td>first copy of @p fe_dg</td> </tr>
   * <tr> <td><code>[false,false,false,false,true,false,false]</code></td>
   * <td><code>FE_DGP<2>(0)</code></td> <td>second copy of @p fe_dg</td> </tr>
   * <tr> <td><code>[false,false,false,false,false,true,true]</code></td>
   * <td><code>FE_BDM<2>(1)</code></td> <td>both components of @p
   * fe_nonprim</td> </tr> </table>
   *
   */
  const FiniteElement<dim, spacedim> &
  get_sub_fe(const ComponentMask &mask) const;

  /**
   * 返回一个与组件 @p n_selected_components
   * 相匹配的所含有限元的引用，该组件从索引为 @p
   * first_component.
   * 的组件开始。更多细节见上面的其他get_sub_fe()函数。
   *
   */
  virtual const FiniteElement<dim, spacedim> &
  get_sub_fe(const unsigned int first_component,
             const unsigned int n_selected_components) const;

  /**
   * 为形状函数 @p index
   * 返回它所属的基本元素，这个基本元素的副本数量（介于0和这个元素的倍数之间），以及这个形状函数在这个基本元素中的索引。
   * 如果该元素不是由其他元素组成的，那么基数和实例总是零，而索引等于形状函数的编号。
   * 如果元素是由其他元素的单个实例组成的（即所有的倍数都是1），这些实例都是标量的，那么这个元素中的基值和dof索引就等于#system_to_component_table。只有在该元素由其他元素组成，并且其中至少有一个元素本身是矢量值的情况下，它才会有所不同。
   * 关于这个函数的典型使用方法，请看上面的类文件中的一个例子。
   * 与system_to_component_index()函数相比，这个函数在矢量值（即非原始）形状函数的情况下也能返回有效值。
   *
   */
  std::pair<std::pair<unsigned int, unsigned int>, unsigned int>
  system_to_base_index(const unsigned int index) const;

  /**
   * 与system_to_base_index()函数相同，但是针对位于面的自由度。因此允许的指数范围是0...#dofs_per_face。
   * 你在应用程序中很少需要这个函数，因为几乎所有的应用程序代码都只需要处理单元格的索引，而不是面的索引。这个函数主要是在库内使用的。
   *
   */
  std::pair<std::pair<unsigned int, unsigned int>, unsigned int>
  face_system_to_base_index(const unsigned int index,
                            const unsigned int face_no = 0) const;

  /**
   * 给定一个基本元素数，返回它将生成的BlockVector的第一个块。
   *
   */
  types::global_dof_index
  first_block_of_base(const unsigned int b) const;

  /**
   * 对于每个向量组件，返回哪个基元实现了这个组件以及这个基元中的哪个向量组件。这个信息只对由几个子元素组成的矢量值有限元素有意义。在这种情况下，人们可能想获得关于实现某个矢量分量的元素的信息，这可以用这个函数和
   * FESystem::base_element() 函数来完成。
   * 如果这是一个标量有限元，那么返回值总是等于一对零。
   *
   */
  std::pair<unsigned int, unsigned int>
  component_to_base_index(const unsigned int component) const;


  /**
   * 返回此块的基元和基元的副本数量。
   *
   */
  std::pair<unsigned int, unsigned int>
  block_to_base_index(const unsigned int block) const;

  /**
   * 这个形状函数的向量块和块内的索引。
   *
   */
  std::pair<unsigned int, types::global_dof_index>
  system_to_block_index(const unsigned int component) const;

  /**
   * 这个组件的向量块。
   *
   */
  unsigned int
  component_to_block_index(const unsigned int component) const;

  //@}

  /**
   * @name  分量和块状矩阵  
     * @{ 
   *
   */

  /**
   * 返回一个分量掩码，其元素数量与此对象的向量分量相同，并且其中正好有一个与给定参数相对应的分量是真的。更多信息见 @ref GlossComponentMask "术语表"
   * 。      @param  标量
   * 一个代表该有限元的单一标量向量分量的对象。
   * @return
   * 一个分量掩码，在所有分量中都是假的，除了与参数相对应的那一个。
   *
   */
  ComponentMask
  component_mask(const FEValuesExtractors::Scalar &scalar) const;

  /**
   * 返回一个分量掩码，其元素数量与此对象的向量分量相同，并且其中与给定参数对应的 <code>dim</code> 分量为真。更多信息见 @ref GlossComponentMask "术语表"
   * 。      @param  矢量
   * 一个表示该有限元的暗淡矢量分量的对象。    @return
   * 一个分量掩码，在所有分量中都是假的，除了与参数相对应的分量。
   *
   */
  ComponentMask
  component_mask(const FEValuesExtractors::Vector &vector) const;

  /**
   * 返回一个分量掩码，其元素数与此对象的向量分量相同，其中与给定参数对应的 <code>dim*(dim+1)/2</code> 分量为真。更多信息见 @ref GlossComponentMask  "术语表"
   * 。      @param  sym_tensor
   * 一个表示该有限元的dim*(dim+1)/2分量的对象，这些分量共同被解释为形成一个对称张量。
   * @return
   * 一个分量掩码，在所有分量中都是假的，除了与参数相对应的那些。
   *
   */
  ComponentMask
  component_mask(
    const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const;

  /**
   *
   */
  ComponentMask
  component_mask(const BlockMask &block_mask) const;

  /**
   * 返回一个块掩码，其元素数与此对象的块数相同，并且其中正好有一个与给定参数相对应的成分是真的。更多信息见 @ref GlossBlockMask "术语表"
   * 。
   * @note
   * 这个函数只有在参数所引用的标量包含一个完整的块时才会成功。换句话说，例如，如果你传递了一个单
   * $x$
   * 速度的提取器，并且这个对象代表一个FE_RaviartThomas对象，那么你选择的单标量对象是一个更大的块的一部分，因此没有代表它的块屏蔽。然后，该函数将产生一个异常。
   * @param  标量
   * 一个代表该有限元的单一标量矢量分量的对象。
   * @return
   * 一个分量掩码，在所有分量中都是假的，除了与参数相对应的那一个。
   *
   */
  BlockMask
  block_mask(const FEValuesExtractors::Scalar &scalar) const;

  /**
   * 返回一个分量掩码，其元素数与此对象的向量分量一样多，并且其中对应于给定参数的 <code>dim</code> 分量为真。更多信息见 @ref GlossBlockMask  "术语表"
   * 。
   * @note  同样的注意事项适用于上述函数的版本。
   * 作为参数传递的提取器对象必须使其对应于完整的块，并且不会分割此元素的块。
   * @param  矢量 一个表示该有限元的dim矢量成分的对象。
   * @return
   * 一个分量掩码，除了与参数对应的分量外，所有分量都是假的。
   *
   */
  BlockMask
  block_mask(const FEValuesExtractors::Vector &vector) const;

  /**
   * 返回一个分量掩码，其元素数与此对象的向量分量相同，其中与给定参数对应的 <code>dim*(dim+1)/2</code> 分量为真。更多信息见 @ref GlossBlockMask "术语表"
   * 。
   * @note  同样的注意事项适用于上述函数的版本。
   * 作为参数传递的提取器对象必须使其对应于完整的块，并且不分割此元素的块。
   * @param  sym_tensor
   * 一个代表该有限元的dim*(dim+1)/2组件的对象，这些组件共同被解释为形成一个对称张量。
   * @return
   * 一个分量掩码，在所有分量中都是假的，除了与参数相对应的那些分量。
   *
   */
  BlockMask
  block_mask(const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const;

  /**
   * @note
   * 这个函数只有在参数所引用的组件包含完整的块时才会成功。换句话说，例如，如果你为单个
   * $x$
   * 速度传递了一个组件掩码，而这个对象代表一个FE_RaviartThomas对象，那么你选择的单个组件是一个更大的块的一部分，因此没有代表它的块掩码。该函数就会产生一个异常。
   * @param  component_mask 选择有限元单个组件的掩码  @return
   * 选择与输入参数的选定块相对应的那些块的掩码。
   *
   */
  BlockMask
  block_mask(const ComponentMask &component_mask) const;

  /**
   * 返回一个元素的常数模式列表。结果表中的行数取决于使用的元素。对于标准元素，该表的行数与元素中的元件和dofs_per_cell列的数量相同。对于有限元的每个分量，返回表中的行包含了该元上常数函数1的基础表示。然而，有一些标量元素存在不止一个常数模式，例如FE_Q_DG0元素。
   * 为了将恒定模式与元素中的实际分量相匹配，返回的数据结构也会返回一个向量，其分量与元素上的恒定模式一样多，其中包含分量编号。
   *
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const;

  //@}

  /**
   * @name 支持点和插值 
     * @{ 
   *
   */

  /**
   * 如果派生的有限元定义了支持点，则返回单元格上试验函数的支持点。
   * 允许某种插值操作的有限元通常有支持点。另一方面，通过例如面的力矩或导数来定义其自由度的元素没有支持点。在这种情况下，返回的域是空的。
   * 如果有限元定义了支持点，那么它们的数量就等于该元的自由度数量。
   * 数组中的点的顺序与<tt>cell->get_dof_indices</tt>函数返回的顺序一致。
   * 关于支持点的详细信息，请参见类的文档。
   * @note
   * 有限元素对该函数的实现以与形状函数相同的顺序返回这些点。形状函数的顺序通常记录在各个有限元类的类文件中。特别是，形状函数（以及随之而来的在该类的类文件中讨论的映射的正交点）将首先遍历那些位于顶点上的形状函数，然后是线，然后是四边形，等等。
   * @note
   * 如果这个元素实现了支持点，那么它将为每个形状函数返回一个这样的点。由于多个形状函数可能被定义在同一位置，这里返回的支持点可能是重复的。一个例子是一个
   * <code>FESystem(FE_Q(1),3)</code>
   * 类型的元素，每个支持点在返回的数组中会出现三次。
   *
   */
  const std::vector<Point<dim>> &
  get_unit_support_points() const;

  /**
   * 返回一个有限元是否有定义的支持点。如果结果为真，那么调用get_unit_support_points()就会产生一个非空数组。
   * 如果一个元素不是由插值形状函数定义的，例如由四边形上的P元素定义的，结果可能为假。通常只有当元素通过要求在某一点上为一，在所有与其他形状函数相关的点上为零来构造其形状函数时，它才会为真。
   * 在组成元素中（即对于FESystem类），如果所有的基础元素都有定义的支持点，那么结果将是真的。FE_Nothing是FES系统中的一个特例，因为它有0个支持点，has_support_points()为假，但是在其他元素中包含FE_Nothing的FES系统将返回true。
   *
   */
  bool
  has_support_points() const;

  /**
   * 返回 @p indexth
   * 形状函数的支持点的位置。如果它不存在，引发一个异常。
   * 默认实现只是从你从get_unit_support_points()得到的数组中返回相应的元素，但是派生元素可以重载这个函数。特别要注意的是，FESystem类重载了它，这样它就可以返回单个基础元素的支持点，如果不是所有基础元素都定义了支持点的话。这样，即使get_unit_support_points()只返回一个空数组，你仍然可以要求获得某些支持点。
   *
   */
  virtual Point<dim>
  unit_support_point(const unsigned int index) const;

  /**
   * 如果派生有限元定义了一些支持点，则返回单元面上的试探函数的支持点。
   * 允许某种插值操作的有限元通常有支持点。另一方面，通过例如面的矩或导数来定义自由度的元素没有支持点。在这种情况下，返回的字段是空的
   * 请注意，有支持点的元素不一定在面上有一些支持点，即使插值点在面上有物理位置。例如，不连续元素的插值点在顶点上，对于更高程度的元素也在面上，但它们没有被定义在面上，因为在这种情况下，来自面的两边（或来自顶点的所有相邻元素）的自由度将相互识别，这不是我们想要的）。因此在逻辑上，这些自由度被定义为属于单元，而不是属于面或顶点。
   * 在这种情况下，返回的元素的长度将为零。
   * 如果有限元定义了支持点，那么它们的数量就等于面的自由度数量（#dofs_per_face）。数组中的点的顺序与<tt>cell->face(face)->get_dof_indices</tt>函数返回的顺序一致。
   * 关于支持点的详细信息，请参见类的文档。
   *
   */
  const std::vector<Point<dim - 1>> &
  get_unit_face_support_points(const unsigned int face_no = 0) const;

  /**
   * 返回一个有限元是否在面上定义了支持点。如果结果为真，那么调用get_unit_face_support_points()就会产生一个非空的向量。
   * 更多信息请参见has_support_points()函数的文档。
   *
   */
  bool
  has_face_support_points(const unsigned int face_no = 0) const;

  /**
   * 与unit_support_point()函数相对应的函数，但针对面。更多信息见那里。
   *
   */
  virtual Point<dim - 1>
  unit_face_support_point(const unsigned int index,
                          const unsigned int face_no = 0) const;

  /**
   * 返回一个广义支持点的向量。
   * @note 该函数返回的向量总是唯一*支持点的最小集合。这与get_unit_support_points()的行为相反，后者对于一个有许多（拉格朗日）基元的FESystem来说，返回一个重复的单位支持点列表。    更多信息请参见 @ref GlossGeneralizedSupport "广义支持点的词汇表条目"
   * 。
   *
   */
  const std::vector<Point<dim>> &
  get_generalized_support_points() const;

  /**
   * 返回一个有限元是否有定义的广义支持点。如果结果为真，那么调用get_generalized_support_points()就会产生一个非空的向量。    更多信息请参见 @ref GlossGeneralizedSupport "广义支持点的词汇条"
   * 。
   *
   */
  bool
  has_generalized_support_points() const;

  /**
   * 对于一个给定的自由度，返回它是否与顶点、直线、四边形或六边形有逻辑联系。
   * 例如，对于连续有限元来说，这与自由度的支持点所处的最低维度对象相吻合。举个例子，对于3D的
   * $Q_1$
   * 元素，每个自由度都由一个形状函数定义，我们通过使用位于单元顶点的支持点进行内插得到。这些点的支持当然延伸到与这个顶点相连的所有边，以及相邻的面和单元内部，但我们说逻辑上自由度是与顶点相关的，因为这是它所关联的最低维度的对象。同样，对于3D中的
   * $Q_2$
   * 元素，支持点位于边缘中点的自由度将从这个函数中得到
   * GeometryPrimitive::line
   * 的值，而那些位于3D中面的中心的自由度将返回
   * GeometryPrimitive::quad. 。
   * 为了使之更加正式，这个函数返回的对象种类代表了该对象，因此，对应于自由度的形状函数的支持，（即该函数
   * "生存
   * "的领域的那一部分）是所有共享该对象的单元的联盟。回到上面的例子，对于3D中的
   * $Q_2$
   * ，支持点位于边缘中点的形状函数在所有共享边缘的单元上都有支持，而不仅仅是共享相邻面的单元，因此函数将返回
   * GeometryPrimitive::line.  另一方面，对于 $DGQ_2$
   * 类型的不连续元素，与内插多项式相关的自由度，其支持点物理上位于一个单元的边界线上，但只在一个单元上非零。因此，它在逻辑上与该单元的内部相关联（即，与2d中的
   * GeometryPrimitive::quad 和3d中的 GeometryPrimitive::hex 相关）。
   * @param[in]  cell_dof_index
   * 一个形状函数或自由度的索引。这个索引必须在
   * <code>[0,dofs_per_cell)</code>  的范围内。
   * @note
   * 这个函数返回的对象的整数值等于它所描述的对象的维度，因此可以在通用编程范式中使用。例如，如果一个自由度与一个顶点相关联，那么这个函数返回
   * GeometryPrimitive::vertex, ，其数值为0（顶点的维度）。
   *
   */
  GeometryPrimitive
  get_associated_geometry_primitive(const unsigned int cell_dof_index) const;


  /**
   * 给定参考单元的（广义）支持点的函数值 $f(\mathbf x)$ ，然后该函数计算元素的节点值是什么，即 $\Psi_i[f]$  ，其中 $\Psi_i$ 是元素的节点函数（另见 @ref GlossNodes "节点值或节点函数"
   * ）。  然后，值 $\Psi_i[f]$ 是<i>interpolates</i>给定函数
   * $f(x)$ 的有限元函数的形状函数的扩展系数，即， $
   * f_h(\mathbf x) = \sum_i \Psi_i[f] \varphi_i(\mathbf x)
   * $ 是 $f$ 与当前元素的有限元插值。  这里描述的操作，例如在 FETools::compute_node_matrix() 函数中使用。    更详细地说，让我们假设当前元素的广义支持点（见 @ref GlossGeneralizedSupport "本词汇表条目"
   * ）是 $\hat{\mathbf x}_i$ ，与当前元素相关的节点函数是
   * $\Psi_i[\cdot]$
   * 。然后，该元素基于广义支持点的事实意味着，如果我们将
   * $\Psi_i$ 应用于（可能是矢量值的）有限元函数 $\varphi$
   * ，结果必须具有 $\Psi_i[\varphi] = f_i(\varphi(\hat{\mathbf x}_i))$
   * 的形式。
   *
   * -换句话说，节点函数 $\Psi_i$ 应用于 $\varphi$ <i>only</i>的值取决于<i>values of $\varphi$ at $\hat{\mathbf x}_i$</i>，而不是取决于其他地方的值，或 $\varphi$ 的积分，或任何其他类型的信息。     $f_i$ 的确切形式取决于元素。例如，对于标量 @ref GlossLagrange "拉格朗日元素"
   * ，我们有，事实上 $\Psi_i[\varphi] = \varphi(\hat{\mathbf x}_i)$
   * 。如果你通过FESystem对象组合多个标量拉格朗日元素，那么
   * $\Psi_i[\varphi] = \varphi(\hat{\mathbf x}_i)_{c(i)}$  其中 $c(i)$ 是
   * FiniteElement::system_to_component_index()
   * 函数的返回值的第一个分量的结果。因此，在这两种情况下，
   * $f_i$
   * 只是身份（在标量情况下）或选择其参数的特定向量分量的函数。
   * 另一方面，对于Raviart-Thomas元素，人们会认为 $f_i(\mathbf
   * y) = \mathbf y \cdot \mathbf n_i$ ，其中 $\mathbf n_i$
   * 是定义形状函数的面的法向量。
   * 鉴于所有这些，这个函数的作用如下。如果你输入一个在所有广义支持点的函数
   * $\varphi$
   * 的值的列表（其中每个值实际上是一个值的向量，其分量与元素的分量一样多），那么这个函数返回一个通过对这些值应用节点函数得到的值的向量。换句话说，如果你传入
   * $\{\varphi(\hat{\mathbf x}_i)\}_{i=0}^{N-1}$
   * 那么你将得到一个向量  $\{\Psi[\varphi]\}_{i=0}^{N-1}$  其中
   * $N$  等于  @p dofs_per_cell.   @param[in]  support_point_values
   * 一个大小为  @p dofs_per_cell  的数组（等于
   * get_generalized_support_points()
   * 函数将返回的点数），每个元素是一个向量，其条目数量与该元素的向量成分相同。这个数组应该包含一个函数在当前元素的广义支持点的值。
   * @param[out]  nodal_values 一个大小为 @p dofs_per_cell
   * 的数组，包含应用于给定函数的元素的节点函数值。
   * @note
   * 只有对于具有琐碎的MappingKind的元素，在实数单元上调用这个函数（转换后的）值是安全的。对于所有其他元素（例如符合H(curl)或H(div)的元素），向量值必须首先被转换到参考单元。
   * @note
   * 鉴于该函数应该做什么，该函数显然只能对实际实现（广义）支持点的元素起作用。没有广义支持点的元素
   *
   * - 例如，节点函数评估积分或函数矩的元素（如FE_Q_Hierarchical）。
   *
   * - 一般来说不能理解这个函数所需的操作。因此，他们可能不会实现它。
   *
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const;

  //@}

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   * 这个函数是虚拟的，因为有限元对象通常是通过指向其基类的指针来访问的，而不是类本身。
   *
   */
  virtual std::size_t
  memory_consumption() const;

  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException1(ExcShapeFunctionNotPrimitive,
                 int,
                 << "The shape function with index " << arg1
                 << " is not primitive, i.e. it is vector-valued and "
                 << "has more than one non-zero vector component. This "
                 << "function cannot be called for these shape functions. "
                 << "Maybe you want to use the same function with the "
                 << "_component suffix?");
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcFENotPrimitive);
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(
    ExcUnitShapeValuesDoNotExist,
    "You are trying to access the values or derivatives of shape functions "
    "on the reference cell of an element that does not define its shape "
    "functions through mapping from the reference cell. Consequently, "
    "you cannot ask for shape function values or derivatives on the "
    "reference cell.");

  /**
   * 试图访问非拉格朗日的有限元的支持点。
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(ExcFEHasNoSupportPoints,
                   "You are trying to access the support points of a finite "
                   "element that either has no support points at all, or for "
                   "which the corresponding tables have not been implemented.");

  /**
   * 试图访问一个没有实现这些矩阵的有限元的嵌入矩阵。
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(ExcEmbeddingVoid,
                   "You are trying to access the matrices that describe how "
                   "to embed a finite element function on one cell into the "
                   "finite element space on one of its children (i.e., the "
                   "'embedding' or 'prolongation' matrices). However, the "
                   "current finite element can either not define this sort of "
                   "operation, or it has not yet been implemented.");

  /**
   * 试图访问一个没有实现这些矩阵的有限元的限制矩阵。
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(ExcProjectionVoid,
                   "You are trying to access the matrices that describe how "
                   "to restrict a finite element function from the children "
                   "of one cell to the finite element space defined on their "
                   "parent (i.e., the 'restriction' or 'projection' matrices). "
                   "However, the current finite element can either not define "
                   "this sort of operation, or it has not yet been "
                   "implemented.");

  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException2(ExcWrongInterfaceMatrixSize,
                 int,
                 int,
                 << "The interface matrix has a size of " << arg1 << "x" << arg2
                 << ", which is not reasonable for the current element "
                    "in the present dimension.");
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcInterpolationNotImplemented);

protected:
  /**
   * 将限制矩阵和延长矩阵的向量重设为正确的大小。对于每个细化案例，除了
   * RefinementCase::no_refinement,
   * 和该细化案例的每个子案例，都会分配一个限制和延长矩阵的空间，关于实际的向量大小，请看限制和延长向量的文档。
   * @param  isotropic_restriction_only
   * 只有各向同性细化所需的限制矩阵被重新引用到合适的大小。
   * @param  isotropic_prolongation_only
   * 只有各向同性细化所需的延长矩阵被重新引用到正确的大小。
   *
   */
  void
  reinit_restriction_and_prolongation_matrices(
    const bool isotropic_restriction_only  = false,
    const bool isotropic_prolongation_only = false);

  /**
   * 投影矩阵的矢量。参见上面的get_restriction_matrix()。构造函数将这些矩阵初始化为零维，这可以由实现它们的派生类来改变。
   * 注意， <code>restriction[refinement_case-1][child]</code>
   * 包括了子类 <code>child</code> 对RefinementCase
   * <code>refinement_case</code>. Here, we use <code>refinement_case-1</code>
   * 的限制矩阵，而不是 <code>refinement_case</code> ，因为对于
   * RefinementCase::no_refinement(=0) 没有限制矩阵可用。
   *
   */
  std::vector<std::vector<FullMatrix<double>>> restriction;

  /**
   * 嵌入矩阵的矢量。见上面的<tt>get_prolongation_matrix()</tt>。构造函数将这些矩阵初始化为零维，这可以由实现它们的派生类来改变。
   * 注意， <code>prolongation[refinement_case-1][child]</code>
   * 包括RefinementCase  <code>refinement_case</code>  的子
   * <code>child</code> 的延长矩阵。这里，我们使用
   * <code>refinement_case-1</code> instead of <code>refinement_case</code>
   * ，因为对于 RefinementCase::no_refinement(=0)
   * ，没有可用的延长矩阵。
   *
   */
  std::vector<std::vector<FullMatrix<double>>> prolongation;

  /**
   * 如果线连接两个单元，其中一个单元被精炼过一次，则指定单元界面两边的道夫所依据的约束。
   * 进一步的细节见派生类的一般描述。
   * 这个字段在一维中显然是无用的，在那里有一个零尺寸。
   *
   */
  FullMatrix<double> interface_constraints;

  /**
   * 单元上的支持点的列表，如果有限元有任何支持点的话。构造函数让这个字段为空，派生类可以写入一些内容。
   * 允许某种插值操作的有限元通常有支持点。另一方面，通过例如面的力矩或导数来定义自由度的元素没有支持点。在这种情况下，这个区域仍然是空的。
   *
   */
  std::vector<Point<dim>> unit_support_points;

  /**
   * 对面的情况也一样。参见get_unit_face_support_points()函数的描述，以讨论什么有助于面的支持点。
   *
   */
  std::vector<std::vector<Point<dim - 1>>> unit_face_support_points;

  /**
   * 用于非拉格朗日元素的插值函数的支持点。
   *
   */
  std::vector<Point<dim>> generalized_support_points;

  /**
   * 用于非拉格朗日元素插值函数的面支持点。
   *
   */
  std::vector<std::vector<Point<dim - 1>>> generalized_face_support_points;

  /**
   * 对于在三维中具有非标准面的方向的面，面（四边形）上的道夫必须被移位，以便与正确的形状函数相结合。给定一个四边形上的局部目标
   * @p index ，如果该面具有非标准的面朝向，即 <code>old_index
   * + shift = new_index</code>
   * ，则返回局部索引的移动。在二维和一维中，不需要进行置换，所以这个向量是空的。在三维中，它的大小为
   * <code> #dofs_per_quad 8 </code>
   * ，其中8是方向的数量，一个脸可以在（face_orientation、face_flip和face_rotation三个bool标志的所有组合）。
   * 这个类的构造函数将这个表填上零，也就是说，根本就没有排列组合。派生的有限元类必须用正确的值来填充这个表。
   *
   */
  std::vector<Table<2, int>> adjust_quad_dof_index_for_face_orientation_table;

  /**
   * 对于三维中具有非标准线方向的线条，线条上的道夫必须被置换，以便与正确的形状函数相结合。给出一条线上的局部道夫
   * @p index ，如果该线有非标准的线方向，即 <code>old_index +
   * shift = new_index</code>
   * ，则返回局部索引的移动。在二维和一维中，不需要进行置换，所以这个向量是空的。在三维中，它的大小为#dofs_per_line。
   * 这个类的构造函数用零来填充这个表，也就是说，根本就没有permutation。派生的有限元类必须用正确的值填充这个向量。
   *
   */
  std::vector<int> adjust_line_dof_index_for_line_orientation_table;

  /**
   * 存储system_to_component_index()将返回的内容。
   *
   */
  std::vector<std::pair<unsigned int, unsigned int>> system_to_component_table;

  /**
   * 面上的线性方向和分量方向之间的映射。这在构造函数中是用默认值填充的，但如果有必要，派生类将不得不覆盖这些信息。
   * 我们所说的分量是指向量分量，而不是指基本元素。因此，只有当一个形状函数只在一个分量中非零时，这些信息才有意义。
   *
   */
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>>
    face_system_to_component_table;

  /**
   * 对于每个形状函数，存储它属于哪个基元和这个基元的哪个实例（如果其倍数大于1），以及它在这个基元中的索引。如果该元素不是由其他元素组成的，那么基元素和实例总是零，而索引等于形状函数的编号。如果该元素是由其他元素的单一实例组成的（即所有的倍数都是1），这些实例都是标量的，那么该元素中的基值和道夫指数就等于#system_to_component_table。只有在该元素由其他元素组成，并且其中至少有一个元素本身是矢量值的情况下，它才会有所不同。
   * 与#system_to_component_table相比，这个数组在向量值（即非原始）形状函数的情况下也有有效值。
   *
   */
  std::vector<std::pair<std::pair<unsigned int, unsigned int>, unsigned int>>
    system_to_base_table;

  /**
   * 同样，对于面的指数也是如此。
   *
   */
  std::vector<
    std::vector<std::pair<std::pair<unsigned int, unsigned int>, unsigned int>>>
    face_system_to_base_table;

  /**
   * 对于每个基数元素，存储由基数产生的块数和它将产生的块向量中的第一个块。
   *
   */
  BlockIndices base_to_block_indices;

  /**
   * 基数元素建立一个组件。
   * 对于每个分量编号<tt>c</tt>，其条目有如下含义。  <dl>
   * <dt><tt>table[c].first.first</tt></dt>  <dd>
   * <tt>c</tt>的基础元素的编号。这是你可以传递给base_element()的索引。
   * </dd>  <dt><tt>table[c].first.second</tt></dt>  <dd>
   * <tt>c</tt>的基础元素的组成部分。这个值在0和这个基础元素的n_components()之间。
   * </dd>  <dt><tt>table[c].second</tt></dt>  <dd>
   * 包含<tt>c</tt>的基础元素的倍数索引。这个值在0和这个基础元素的
   * element_multiplicity() 之间。 </dd>   </dl>
   * 这个变量被这个类的构造函数设置为正确的大小，但是需要被派生类初始化，除非它的大小是1，而且唯一的条目是0，对于标量元素就是这样。
   * 在这种情况下，由基类初始化就可以了。
   * @note 此表由 FETools::Compositing::build_cell_tables(). 填写。
   *
   */
  std::vector<std::pair<std::pair<unsigned int, unsigned int>, unsigned int>>
    component_to_base_table;

  /**
   * 一个决定限制矩阵是连接还是相加的标志。更多信息请参见通用类文档中关于限制矩阵的讨论。
   *
   */
  const std::vector<bool> restriction_is_additive_flags;

  /**
   * 对于每个形状函数，给出一个bools的向量（大小等于这个有限元素所具有的向量分量的数量），表明这些形状函数中的每个分量是非零的。
   * 对于原始元素，只有一个非零分量。
   *
   */
  const std::vector<ComponentMask> nonzero_components;

  /**
   * 这个数组持有#nonzero_components元素各自条目中的多少个值是非零的。因此，这个数组是一个捷径，可以让我们更快地获得这些信息，而不是在每次请求获得这些信息时都要计算非零条目。该字段在该类的构造函数中被初始化。
   *
   */
  const std::vector<unsigned int> n_nonzero_components_table;

  /**
   * 存储所有形状函数是否是原始的。由于找出这个是一个非常常见的操作，我们缓存了这个结果，即在构造函数中计算出这个值，以便更简单地访问。
   *
   */
  const bool cached_primitivity;

  /**
   * 返回接口约束矩阵的大小。由于每个派生的有限元类在初始化它们的尺寸时都需要这个，所以它被放在这个函数中，以避免每次都要重新计算这些矩阵的与尺寸有关的尺寸。
   * 注意，有些元素没有实现某些多项式程度的接口约束。在这种情况下，这个函数仍然返回这些矩阵实现时应该有的大小，但实际的矩阵是空的。
   *
   */
  TableIndices<2>
  interface_constraints_size() const;

  /**
   * 给出每个形状函数的非零分量的模式，为每个条目计算每个形状函数有多少分量是非零的。
   * 这个函数在这个类的构造函数中使用。
   *
   */
  static std::vector<unsigned int>
  compute_n_nonzero_components(
    const std::vector<ComponentMask> &nonzero_components);

  /**
   * 给定一组更新标志，计算哪些其他数量<i>also</i>需要被计算以满足给定标志的请求。
   * 然后返回原始标志集和刚刚计算的标志的组合。
   * 例如，如果 @p update_flags
   * 包含update_gradients，一个有限元类通常需要计算雅各布矩阵的逆值，以便将参考单元上的形状函数的梯度旋转到实际单元上。然后，它不仅会返回update_gradients，还会返回update_covariant_transformation，这个标志使映射类产生雅各布矩阵的逆。
   * 关于这个函数和FEValues之间互动的广泛讨论可以在 @ref
   * FE_vs_Mapping_vs_FEValues 文档模块中找到。      @see  UpdateFlags
   *
   */
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const = 0;

  /**
   * 创建一个内部数据对象，并返回一个指针，然后该函数的调用者将拥有该指针的所有权。然后，每次在具体单元上评估有限元形状函数及其导数时，这个对象将被传递给 FiniteElement::fill_fe_values()
   * 。因此，这里创建的对象被派生类用作评估形状函数时的抓取对象，以及存储可以预先计算一次并在每个单元上重复使用的信息（例如，在参考单元上评估形状函数的值和梯度，以便以后在将这些值转换到具体单元时重复使用）。
   * 在为给定的映射和有限元对象初始化FEValues对象的过程中，这个函数是第一个被调用的。返回的对象随后将被传递给
   * FiniteElement::fill_fe_values()
   * 的具体单元，它本身将把其输出放到
   * internal::FEValuesImplementation::FiniteElementRelatedData.
   * 类型的对象中。由于可能有数据已经可以在参考单元上以其<i>final</i>的形式计算，这个函数也接收对
   * internal::FEValuesImplementation::FiniteElementRelatedData
   * 对象的引用作为其最后参数。当与该函数返回的InternalDataBase对象一起使用时，这个输出参数保证总是同一个参数。换句话说，在返回的对象和
   * @p output_data
   * 对象中，从头数据和最终数据的细分是如下的。如果数据可以在参考单元上以以后在具体单元上需要的确切形式预先计算出来，那么这个函数应该已经把它放到了
   * @p output_data
   * 对象中。一个例子是通常拉格朗日元素的正交点的形状函数值，它在具体单元上与参考单元上的相同。另一方面，如果某些数据可以预先计算，以便在具体单元上进行计算<i>cheaper</i>，那么它应该被放入返回的对象中，以便以后在派生类的实现中重新使用
   * FiniteElement::fill_fe_values().  ]
   * 一个例子是拉格朗日元素的形状函数在参考单元上的梯度：为了计算形状函数在具体单元上的梯度，我们必须将参考单元上的梯度乘以映射的雅各布反值；因此，我们不能在调用当前函数时已经计算出具体单元上的梯度，但我们至少可以预先计算出参考单元上的梯度，并将其保存在返回的对象中。
   * 关于这个函数和 FEValues 之间的互动的广泛讨论可以在
   * @ref FE_vs_Mapping_vs_FEValues
   * 文档模块中找到。也请参见InternalDataBase类的文档。
   * @param[in]  update_flags
   * 一组UpdateFlags值，描述FEValues对象要求有限元计算哪种信息。这组标志也可以包括有限元不能计算的信息，例如，与映射产生的数据有关的标志。这个函数的实现需要在返回的对象中设置所有的数据字段，这些字段是产生这些标志所指定的有限元相关数据所必需的，并且可能已经预先计算了这些信息的一部分，如上文所述。元素可能希望在
   * InternalDataBase::update_each
   * 中存储这些更新标志（或这些标志的子集），以便在
   * FiniteElement::fill_fe_values()
   * 被调用时知道他们应该计算什么  @param[in]  映射
   * 对用于计算形状函数值和导数的映射的引用。
   * @param[in]  正交 对描述形状函数应被评估的对象的引用。
   * @param[out]  output_data
   * 对对象的引用，FEValues将与这里返回的对象一起使用，
   * FiniteElement::fill_fe_values()
   * 的实现将在这里放置所需的信息。
   * 这允许当前函数已经预先计算了可以在参考单元上计算的信息，如上所述。FEValues保证这个输出对象和当前函数返回的对象将总是一起使用。
   * @return
   * 一个指向从InternalDataBase派生的类型的对象的指针，派生类可以用它来存储可以预先计算的抓取数据，或者用于抓取数组，然后只需要分配一次。
   * 调用网站承担这个对象的所有权，并在不再需要时将其删除。
   *
   */
  virtual std::unique_ptr<InternalDataBase>
  get_data(const UpdateFlags             update_flags,
           const Mapping<dim, spacedim> &mapping,
           const Quadrature<dim> &       quadrature,
           dealii::internal::FEValuesImplementation::
             FiniteElementRelatedData<dim, spacedim> &output_data) const = 0;

  /**
   * 像get_data()一样，但是返回一个对象，该对象以后将被用于评估单元格面上正交点的形状函数信息。然后该对象将被用于调用
   * FiniteElement::fill_fe_face_values().
   * 的实现，更多信息请参见get_data()的文档。
   * 这个函数的默认实现将面的正交转换为具有适当正交点位置的单元格正交，并以此调用上述必须在派生类中实现的get_data()函数。
   * @param[in]  update_flags
   * 一组UpdateFlags值，描述FEValues对象要求有限元计算哪种信息。这组标志也可以包括有限元不能计算的信息，例如，与映射产生的数据有关的标志。这个函数的实现需要在返回的对象中设置所有的数据字段，这些字段是产生这些标志所指定的有限元相关数据所必需的，并且可能已经预先计算了这些信息的一部分，如上文所述。元素可能希望在
   * InternalDataBase::update_each
   * 中存储这些更新标志（或这些标志的子集），以便在调用
   * FiniteElement::fill_fe_face_values() 时知道他们应该计算什么
   * @param[in]  映射
   * 对用于计算形状函数的值和导数的映射的引用。
   * @param[in]  正交
   * 指对描述形状函数应被评估的对象的引用。    @param[out]
   * output_data
   * 对FEValues将与此处返回的对象一起使用的对象的引用，
   * FiniteElement::fill_fe_face_values()
   * 的实现将在此处放置所要求的信息。这允许当前函数已经预先计算了可以在参考单元上计算的信息，如上所述。FEValues保证这个输出对象和当前函数返回的对象将总是一起使用。
   * @return
   * 一个指向从InternalDataBase派生的类型的对象的指针，派生类可以用它来存储可以预先计算的抓取数据，或者用于抓取数组，然后只需要分配一次。
   * 调用网站承担这个对象的所有权，并在不再需要时将其删除。
   *
   */
  virtual std::unique_ptr<InternalDataBase>
  get_face_data(const UpdateFlags               update_flags,
                const Mapping<dim, spacedim> &  mapping,
                const hp::QCollection<dim - 1> &quadrature,
                dealii::internal::FEValuesImplementation::
                  FiniteElementRelatedData<dim, spacedim> &output_data) const;

  /**
   * @deprecated  使用带有 hp::QCollection 参数的版本。
   *
   */
  virtual std::unique_ptr<InternalDataBase>
  get_face_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim - 1> &   quadrature,
    internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
      &output_data) const;

  /**
   * 与get_data()类似，但返回一个对象，该对象以后将用于评估单元格面子上的正交点的形状函数信息。该对象随后将被用于调用
   * FiniteElement::fill_fe_subface_values().
   * 的实现，更多信息请参见get_data()的文档。
   * 该函数的默认实现将面的正交转换为具有适当正交点位置的单元格正交，并以此调用上述必须在派生类中实现的get_data()函数。
   * @param[in]  update_flags
   * 一组UpdateFlags值，描述FEValues对象要求有限元计算哪种信息。这组标志也可以包括有限元不能计算的信息，例如，与映射产生的数据有关的标志。这个函数的实现需要在返回的对象中设置所有的数据字段，这些字段是产生这些标志所指定的有限元相关数据所必需的，并且可能已经预先计算了这些信息的一部分，如上文所述。元素可能希望在
   * InternalDataBase::update_each
   * 中存储这些更新标志（或这些标志的子集），以便在
   * FiniteElement::fill_fe_subface_values()
   * 被调用时知道他们应该计算什么  @param[in]  映射
   * 对用于计算形状函数的值和导数的映射的引用。
   * @param[in]  quadrature
   * 对描述形状函数应被评估的对象的引用。    @param[out]
   * output_data
   * 对对象的引用，FEValues将与这里返回的对象一起使用，
   * FiniteElement::fill_fe_subface_values()
   * 的实现将在这里放置所需的信息。这允许当前函数已经预先计算了可以在参考单元上计算的信息，如上所述。FEValues保证这个输出对象和当前函数返回的对象将总是一起使用。
   * @return
   * 一个指向从InternalDataBase派生的类型的对象的指针，派生类可以用它来存储可以预先计算的抓取数据，或者用于抓取数组，然后只需要分配一次。
   * 调用网站承担这个对象的所有权，并在不再需要时将其删除。
   *
   */
  virtual std::unique_ptr<InternalDataBase>
  get_subface_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim - 1> &   quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const;

  /**
   * 计算由第一个参数表示的单元上的形状函数信息。派生类将不得不根据它们所代表的元素种类来实现这个函数。它是由 FEValues::reinit().
   * 调用的。从概念上讲，这个函数在由这个函数的正交参数所描述的那些正交点的映射位置上评估形状函数及其导数。在许多情况下，计算形状函数的导数（在某些情况下也计算形状函数的值）需要利用从参考单元到实际单元的映射；这一信息可以从调用此函数之前为当前单元填充的
   * @p mapping_data
   * 对象中获取，或者通过调用与当前单元对应的 @p
   * mapping_internal 对象的Mapping对象的成员函数。
   * 这个函数计算出来的信息被用来填充这个函数的输出参数的各种成员变量。该结构中的哪些成员变量应该被填充，由存储在传递给该函数的对象的
   * FiniteElement::InternalDataBase::update_each
   * 字段中的更新标志决定。这些标志通常由
   * FiniteElement::get_data(),  FiniteElement::get_face_date() 和
   * FiniteElement::get_subface_data()
   * （或者更具体地说，这些函数在派生类中的实现）设置。
   * 关于这个函数和FEValues之间的互动的广泛讨论可以在 @ref
   * FE_vs_Mapping_vs_FEValues 文档模块中找到。      @param[in]  cell
   * 三角形中的单元格，该函数要计算从参考单元格到三角形的映射。
   * @param[in]  cell_similarity
   * 作为第一个参数的单元格是否是最近一次调用此函数的单元格的简单翻译、旋转等。这个信息是通过匹配前一个单元和当前单元之间的顶点（由三角结构存储）简单计算出来的。这里传递的值可能被这个函数的实现所修改，然后应该被返回（见关于这个函数的返回值的讨论）。
   * @param[in]  quadrature
   * 对当前评估中使用的正交公式的引用。这个正交对象与创建
   * @p internal_data
   * 对象时使用的对象相同。然后，当前对象负责在此对象所代表的正交点的映射位置评估形状函数。
   * @param[in]  mapping
   * 对用于从参考单元映射到当前单元的映射对象的引用。在调用当前函数之前，这个对象被用来计算
   * @p mapping_data 对象中的信息。它也是通过 Mapping::get_data().
   * 创建 @p
   * mapping_internal对象的映射对象。你最需要对这个映射对象的引用来调用
   * Mapping::transform()
   * ，将梯度和高导数从参考单元转换到当前单元。
   * @param[in]  mapping_internal 一个特定于映射对象的对象。
   * 映射选择在其中存储的内容与当前函数无关，但你可能必须将这个对象的引用传递给映射类的某些函数（例如，如果你需要从当前函数中调用它们，那么
   * Mapping::transform()) 。    @param[in]  mapping_data
   * Mapping::fill_fe_values()
   * 函数将对应于当前单元格的映射信息写入的输出对象。这包括，例如，可能与当前函数相关的映射的Jacobian，以及
   * FEValues::reinit() 从映射中要求的其他信息。    @param[in]
   * fe_internal
   * 一个对先前由get_data()创建的对象的引用，可用于存储映射在参考单元上可以计算的一次信息。参见
   * FiniteElement::InternalDataBase
   * 类的文档，以了解对这些对象的用途的广泛描述。
   * @param[out]  output_data
   * 对成员变量应被计算的对象的引用。并非所有这个参数的成员都需要被填充；哪些成员需要被填充是由存储在
   * @p fe_internal 对象中的更新标志决定的。
   * @note  FEValues确保这个函数总是用同一对 @p fe_internal 和 @p
   * output_data
   * 对象调用。换句话说，如果这个函数的实现知道它在之前的调用中已经把一个数据写入了输出参数，那么在以后的调用中，如果实现知道这是同一个值，就没有必要再把它复制到那里。
   *
   */
  virtual void
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                     mapping_data,
    const InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const = 0;

  /**
   * 这个函数等同于 FiniteElement::fill_fe_values(),
   * ，但用于单元格的面。关于其目的的广泛讨论见那里。
   * 它被 FEFaceValues::reinit().  @param[in]
   * 单元格所调用，该函数要计算从参考单元格到的映射的三角形单元格。
   * @param[in]  face_no
   * 我们当前考虑的面的编号，在前一个参数所指定的单元格的面中索引。
   * @param[in]  quadrature
   * 对当前评估中使用的正交公式的引用。这个正交对象与创建
   * @p internal_data
   * 对象时使用的对象相同。然后，当前对象负责在此对象所代表的正交点的映射位置评估形状函数。
   * @param[in]  mapping
   * 对用于从参考单元映射到当前单元的映射对象的引用。在调用当前函数之前，该对象被用来计算
   * @p mapping_data 对象中的信息。它也是通过 Mapping::get_data().
   * 创建 @p
   * mapping_internal对象的映射对象。你最需要对这个映射对象的引用来调用
   * Mapping::transform()
   * ，将梯度和高导数从参考单元转换到当前单元。
   * @param[in]  mapping_internal 一个特定于映射对象的对象。
   * 映射选择在其中存储什么与当前函数无关，但你可能必须将这个对象的引用传递给映射类的某些函数（例如，
   * Mapping::transform())  如果你需要从当前函数中调用它们。
   * @param[in]  mapping_data  Mapping::fill_fe_values()
   * 函数将对应于当前单元的映射信息写入的输出对象。这包括，例如，可能与当前函数相关的映射的雅各布，以及
   * FEValues::reinit()  从映射中要求的其他信息。    @param[in]
   * fe_internal
   * 一个对先前由get_data()创建的对象的引用，可用于存储映射在参考单元上可以计算的信息。参见
   * FiniteElement::InternalDataBase
   * 类的文档，了解这些对象的用途的广泛描述。
   * @param[out]  output_data
   * 对成员变量应被计算的对象的引用。并非所有这个参数的成员都需要被填充；哪些成员需要被填充是由存储在
   * @p fe_internal 对象内的更新标志决定的。
   *
   */
  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const hp::QCollection<dim - 1> &                            quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                     mapping_data,
    const InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const;

  /**
   * @deprecated  使用带有 hp::QCollection 参数的版本。
   *
   */
  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &                     mapping_data,
    const InternalDataBase &fe_internal,
    internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
      &output_data) const;

  /**
   * 这个函数等同于 FiniteElement::fill_fe_values(),
   * ，但用于单元格面的子代。关于其目的的广泛讨论见那里。它被
   * FESubfaceValues::reinit().  @param[in]
   * 单元格所调用，该函数要为其计算从参考单元格到的映射。
   * @param[in]  face_no
   * 我们当前考虑的面的编号，在前一个参数所指定的单元格的面中索引。
   * @param[in]  sub_no
   * 我们当前考虑的子面的编号，即一个面的子的编号，在前面参数指定的面的子中索引。
   * @param[in]  quadrature
   * 对当前评估中使用的正交公式的引用。这个正交对象与创建
   * @p internal_data
   * 对象时使用的对象相同。然后，当前对象负责在此对象所代表的正交点的映射位置评估形状函数。
   * @param[in]  mapping
   * 对用于从参考单元映射到当前单元的映射对象的引用。在调用当前函数之前，这个对象被用来计算
   * @p mapping_data 对象中的信息。它也是通过 Mapping::get_data().
   * 创建 @p
   * mapping_internal对象的映射对象。你最需要这个映射对象的引用来调用
   * Mapping::transform()
   * ，将梯度和高导数从参考单元转换到当前单元。
   * @param[in]  mapping_internal 一个特定于映射对象的对象。
   * 映射选择在其中存储什么与当前函数无关，但你可能必须将这个对象的引用传递给映射类的某些函数（例如，
   * Mapping::transform())  如果你需要从当前函数中调用它们。
   * @param[in]  mapping_data  Mapping::fill_fe_values()
   * 函数将对应于当前单元格的映射信息写入的输出对象。这包括，例如，可能与当前函数相关的映射的Jacobians，以及
   * FEValues::reinit() 从映射中要求的其他信息。    @param[in]
   * fe_internal
   * 一个对先前由get_data()创建的对象的引用，可用于存储映射在参考单元上可以计算的一次信息。参见
   * FiniteElement::InternalDataBase
   * 类的文档，了解这些对象的用途的广泛描述。
   * @param[out]  output_data
   * 对成员变量应被计算的对象的引用。并非所有这个参数的成员都需要被填充；哪些成员需要被填充是由存储在
   * @p fe_internal 对象内的更新标志决定的。
   *
   */
  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          sub_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                     mapping_data,
    const InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const = 0;

  friend class InternalDataBase;
  friend class FEValuesBase<dim, spacedim>;
  friend class FEValues<dim, spacedim>;
  friend class FEFaceValues<dim, spacedim>;
  friend class FESubfaceValues<dim, spacedim>;
  friend class FESystem<dim, spacedim>;

  // explicitly check for sensible template arguments, but not on windows
  // because MSVC creates bogus warnings during normal compilation
#ifndef DEAL_II_MSVC
  static_assert(dim <= spacedim,
                "The dimension <dim> of a FiniteElement must be less than or "
                "equal to the space dimension <spacedim> in which it lives.");
#endif
};


//----------------------------------------------------------------------//


template <int dim, int spacedim>
inline const FiniteElement<dim, spacedim> &FiniteElement<dim, spacedim>::
                                           operator[](const unsigned int fe_index) const
{
  (void)fe_index;
  Assert(fe_index == 0,
         ExcMessage("A fe_index of zero is the only index allowed here"));
  return *this;
}



template <int dim, int spacedim>
inline std::pair<unsigned int, unsigned int>
FiniteElement<dim, spacedim>::system_to_component_index(
  const unsigned int index) const
{
  AssertIndexRange(index, system_to_component_table.size());
  Assert(is_primitive(index),
         (typename FiniteElement<dim, spacedim>::ExcShapeFunctionNotPrimitive(
           index)));
  return system_to_component_table[index];
}



template <int dim, int spacedim>
inline unsigned int
FiniteElement<dim, spacedim>::n_base_elements() const
{
  return base_to_block_indices.size();
}



template <int dim, int spacedim>
inline unsigned int
FiniteElement<dim, spacedim>::element_multiplicity(
  const unsigned int index) const
{
  return static_cast<unsigned int>(base_to_block_indices.block_size(index));
}



template <int dim, int spacedim>
inline unsigned int
FiniteElement<dim, spacedim>::component_to_system_index(
  const unsigned int component,
  const unsigned int index) const
{
  AssertIndexRange(component, this->n_components());
  const std::vector<std::pair<unsigned int, unsigned int>>::const_iterator it =
    std::find(system_to_component_table.begin(),
              system_to_component_table.end(),
              std::pair<unsigned int, unsigned int>(component, index));

  Assert(it != system_to_component_table.end(),
         ExcMessage("You are asking for the number of the shape function "
                    "within a system element that corresponds to vector "
                    "component " +
                    Utilities::int_to_string(component) +
                    " and within this to "
                    "index " +
                    Utilities::int_to_string(index) +
                    ". But no such "
                    "shape function exists."));
  return std::distance(system_to_component_table.begin(), it);
}



template <int dim, int spacedim>
inline std::pair<unsigned int, unsigned int>
FiniteElement<dim, spacedim>::face_system_to_component_index(
  const unsigned int index,
  const unsigned int face_no) const
{
  AssertIndexRange(
    index,
    face_system_to_component_table[this->n_unique_faces() == 1 ? 0 : face_no]
      .size());

  // in debug mode, check whether the
  // function is primitive, since
  // otherwise the result may have no
  // meaning
  //
  // since the primitivity tables are
  // all geared towards cell dof
  // indices, rather than face dof
  // indices, we have to work a
  // little bit...
  //
  // in 1d, the face index is equal
  // to the cell index
  Assert(is_primitive(this->face_to_cell_index(index, face_no)),
         (typename FiniteElement<dim, spacedim>::ExcShapeFunctionNotPrimitive(
           index)));

  return face_system_to_component_table[this->n_unique_faces() == 1 ?
                                          0 :
                                          face_no][index];
}



template <int dim, int spacedim>
inline std::pair<std::pair<unsigned int, unsigned int>, unsigned int>
FiniteElement<dim, spacedim>::system_to_base_index(
  const unsigned int index) const
{
  AssertIndexRange(index, system_to_base_table.size());
  return system_to_base_table[index];
}



template <int dim, int spacedim>
inline std::pair<std::pair<unsigned int, unsigned int>, unsigned int>
FiniteElement<dim, spacedim>::face_system_to_base_index(
  const unsigned int index,
  const unsigned int face_no) const
{
  AssertIndexRange(
    index,
    face_system_to_base_table[this->n_unique_faces() == 1 ? 0 : face_no]
      .size());
  return face_system_to_base_table[this->n_unique_faces() == 1 ? 0 : face_no]
                                  [index];
}



template <int dim, int spacedim>
inline types::global_dof_index
FiniteElement<dim, spacedim>::first_block_of_base(
  const unsigned int index) const
{
  return base_to_block_indices.block_start(index);
}



template <int dim, int spacedim>
inline std::pair<unsigned int, unsigned int>
FiniteElement<dim, spacedim>::component_to_base_index(
  const unsigned int index) const
{
  AssertIndexRange(index, component_to_base_table.size());

  return component_to_base_table[index].first;
}



template <int dim, int spacedim>
inline std::pair<unsigned int, unsigned int>
FiniteElement<dim, spacedim>::block_to_base_index(
  const unsigned int index) const
{
  return base_to_block_indices.global_to_local(index);
}



template <int dim, int spacedim>
inline std::pair<unsigned int, types::global_dof_index>
FiniteElement<dim, spacedim>::system_to_block_index(
  const unsigned int index) const
{
  AssertIndexRange(index, this->n_dofs_per_cell());
  // The block is computed simply as
  // first block of this base plus
  // the index within the base blocks
  return std::pair<unsigned int, types::global_dof_index>(
    first_block_of_base(system_to_base_table[index].first.first) +
      system_to_base_table[index].first.second,
    system_to_base_table[index].second);
}



template <int dim, int spacedim>
inline bool
FiniteElement<dim, spacedim>::restriction_is_additive(
  const unsigned int index) const
{
  AssertIndexRange(index, this->n_dofs_per_cell());
  return restriction_is_additive_flags[index];
}



template <int dim, int spacedim>
inline const ComponentMask &
FiniteElement<dim, spacedim>::get_nonzero_components(const unsigned int i) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  return nonzero_components[i];
}



template <int dim, int spacedim>
inline unsigned int
FiniteElement<dim, spacedim>::n_nonzero_components(const unsigned int i) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  return n_nonzero_components_table[i];
}



template <int dim, int spacedim>
inline bool
FiniteElement<dim, spacedim>::is_primitive() const
{
  return cached_primitivity;
}



template <int dim, int spacedim>
inline bool
FiniteElement<dim, spacedim>::is_primitive(const unsigned int i) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());

  // return primitivity of a shape
  // function by checking whether it
  // has more than one non-zero
  // component or not. we could cache
  // this value in an array of bools,
  // but accessing a bit-vector (as
  // std::vector<bool> is) is
  // probably more expensive than
  // just comparing against 1
  //
  // for good measure, short circuit the test
  // if the entire FE is primitive
  return (is_primitive() || (n_nonzero_components_table[i] == 1));
}



template <int dim, int spacedim>
inline GeometryPrimitive
FiniteElement<dim, spacedim>::get_associated_geometry_primitive(
  const unsigned int cell_dof_index) const
{
  AssertIndexRange(cell_dof_index, this->n_dofs_per_cell());

  // just go through the usual cases, taking into account how DoFs
  // are enumerated on the reference cell
  if (cell_dof_index < this->get_first_line_index())
    return GeometryPrimitive::vertex;
  else if (cell_dof_index < this->get_first_quad_index(0))
    return GeometryPrimitive::line;
  else if (cell_dof_index < this->get_first_hex_index())
    return GeometryPrimitive::quad;
  else
    return GeometryPrimitive::hex;
}



DEAL_II_NAMESPACE_CLOSE

#endif


