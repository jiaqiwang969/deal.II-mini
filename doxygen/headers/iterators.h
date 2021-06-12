//include/deal.II-translator/A-headers/iterators_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2021 by the deal.II authors
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

/**
 * @defgroup Iterators Iterators on mesh-like containers   
 * @{ 
 * deal.II有几个类，在概念上被理解为meshes。除了明显的Triangulation之外，这些类还有，例如，DoFHandler和
 * hp::DoFHandler.
 * 所有这些类都定义了一组迭代器，允许用户遍历整个网格，即构成网格的单元格、面、边等的集合，或其中的一部分。这些迭代器在某种意义上都是由TriaIteratorclass派生出来的。
 * 基本上，TriaIterator的模板签名是
 * @code
 * TriaIterator<Accessor>
 * @endcode
 *
 * 从概念上讲，这个类型代表了类似于指向由
 * <code>Accessor</code> 类所代表的对象的指针的东西。
 * 通常情况下，你不会直接使用实际的类名，而是采用网状类所提供的一个类型定义，如<code>typename
 * Triangulation::cell_iterator</code>.
 * 在进入这个问题之前，让我们先讨论一下迭代器的概念，然后再深入研究访问器的作用。
 * 在C++中，迭代器就像指针一样，使用<tt>operator
 * ++</tt>增加到最近的元素，使用<tt>operator减少到之前的元素。
 *
 * - </tt>。我们也可以用加法运算符<tt>it=it+n</tt>来跳转<tt>n</tt>元素，并相应地将一些元素向后移动。此外，按照标准模板库的传统，网格提供了成员函数<tt>begin()</tt>和<tt>end()</tt>，分别提供一个集合的第一个元素和一个超前的迭代器。由于有许多不同的迭代器可用，实际上有一整个系列的这类函数，比如<tt>begin_active()</tt>, <tt>begin_face()</tt>，等等。
 * 就C++标准中定义的迭代器概念而言，deal.II网格迭代器是双向迭代器：它们可以被递增和递减，但是像<tt>it=it+n</tt>这样的操作需要的计算时间与<tt>n</tt>成正比，因为它被实现为一个<tt>n</tt>单独单位的递增序列。请注意，这与下一个更专业的迭代器概念，即随机访问迭代器形成了对比，对任意对象的访问只需要恒定时间，而不是线性时间。
 *
 *
 * @section  IteratorsAndSets 迭代器作为指向对象集合的指针
 * 如上所述，deal.II中的迭代器可以被认为是对构成网格的对象进行整体迭代。(这些对象是线、四边形和六边形，并由作为模板参数的Accessor类的类型表示给迭代器)。这表明我们可以将三角形视为单元格和其他对象的集合，这些单元格和其他对象通过某种数据结构联系在一起，就像链接列表是以线性方式连接对象的数据结构一样。
 * 交易二中的三角形确实可以用这种方式来考虑。特别是，它们使用规则树之林的计算概念来存储它们的数据。这可以理解为如下几点。考虑到粗网格的单元有根；然后，如果这些粗网格的一个单元被细化，它将有2<sup>dim</sup>子，而这些子又可以，但不一定要有自己的2<sup>dim</sup>子，等等。这意味着，粗网中的每个细胞都可以被认为是二叉树（1d）、四叉树（2d）或八叉树（3d）的根。这些从粗网格单元中产生的树的集合构成了完整描述三角形的森林，包括所有的活动和非活动单元。特别是，活动单元是树中那些没有后代的终端节点，即没有进一步细化的单元。相应地，非活动单元对应于树中没有后裔的节点，即被进一步细化的单元。
 * 三角形包含线（每个线可能有2个孩子）、四边形（每个四边形可能有4个孩子）和六边形（每个六边形没有或有8个孩子）的森林。根据维度的不同，这些对象也被称为单元或面。
 * 迭代器在这种森林的元素上循环。通常的迭代器在森林的所有节点上循环，而活动迭代器以相同的顺序在元素上循环，但跳过所有非活动条目，因此只访问终端节点（即活动单元、面等）。遍历森林中的元素有很多方法，例如广度优先或深度优先。根据用于存储森林的数据结构的类型，有些方法比其他方法更有效。目前，在deal.II中，迭代器遍历森林的方式是广度优先。也就是说，迭代器首先访问粗网格的所有元素（单元格、面等），然后再转到直属层的所有元素，即粗网格对象的直属子节点；之后是粗网格的子节点，以此类推。然而，必须注意的是，程序不应该依赖这种特定的树形遍历顺序：这被认为是一个实现细节，可以在不同的版本之间改变，即使我们认为这在目前是一个不太可能的选择。
 *
 *
 * @section  迭代器的区别 不同种类的迭代器
 * 迭代器有两个属性：它们指向什么（即Accessor模板参数的类型），以及它们所迭代的集合的确切定义。一般来说，迭代器总是被声明为
 * @code
 * KindIterator<Accessor>
 * @endcode
 *
 * 这里，<tt>Kind</tt>决定了一个访问器需要有什么属性才能被这个迭代器访问（或者省略）。比如说。
 * @code
 * Iterator<Accessor>
 * @endcode
 * 迭代所有构成网格的Accessor类型的对象（例如所有的单元格，无论它们是否被进一步细化并有子代），而
 * @code
 * ActiveIterator<Accessor>
 * @endcode
 * 因此，主动迭代器操作的是普通迭代器所操作的对象的一个子集，即那些拥有主动属性的对象。请注意，这与我们操作的对象的种类无关：所有有效的访问器类都必须为迭代器类提供一个方法来找出它们是否是活动的。
 * （为了完整起见，让我们提到还有第三种迭代器。"rawiterators
 * "也会遍历三角形中未使用的对象，但出于效率的考虑还是会被分配。用户代码不应该使用rawiterators，它们只为库的内部目的而存在）。)
 * 这是通过使用一个FilteredIterator&lt;BaseIterator&gt;类型的对象实现的，其中BaseIterator通常是上面讨论的标准迭代器之一。
 * FilteredIterator在其构造函数中得到了一个额外的谓词，并将跳过所有该谓词评估为<tt>false</tt>的对象。在命名空间IteratorFilters中可以找到已经实现的谓词的集合。
 *
 *
 * @subsection  迭代器循环 迭代对象
 * 所有相同类型的迭代器和迭代相同类型的几何对象都以相同的顺序遍历网格。以这个代码为例。
 * @code
 * Triangulation<dim> tria;
 * DoFHandler<dim>    dof1(tria);
 * DoFHandler<dim>    dof2(tria);
 * ...
 * typename Trianguation<dim>::cell_iterator ti  = tria.begin();
 * typename DoFHandler<dim>::cell_iterator   di1 = dof1.begin();
 * typename DoFHandler<dim>::cell_iterator   di2 = dof2.begin();
 * ...
 * while (ti != tria.end())
 * {
 *  // do something
 *  ++ti;
 *  ++di1;
 *  ++di2;
 * }
 * @endcode
 *
 * 这里，所有的迭代器总是指向同一个网格单元，即使<tt>DoFHandler</tt>和<tt>Triangulation</tt>是非常不同的类，即使DoFHandler处理的是不同的有限元：它们都以相同的顺序访问单元，区别只在于Accessor。如上所述，迭代器遍历对象森林的顺序实际上是明确的，但是应用程序不应该假定任何这样的顺序，而应该考虑这是库的实现细节。
 * 与上面的例子相对应，在下面的片段中，迭代器遍历活动对象的顺序对所有迭代器都是一样的，与前面的例子不同的是，这里我们只考虑活动单元。
 * @code
 * typename Trianguation<dim>::active_cell_iterator ti  = tria.begin_active();
 * typename DoFHandler<dim>::active_cell_iterator   di1 = dof1.begin_active();
 * typename DoFHandler<dim>::active_cell_iterator   di2 = dof2.begin_active();
 * ...
 * while (ti != tria.end())
 * {
 *  // do something
 *  ++ti;
 *  ++di1;
 *  ++di2;
 * }
 * @endcode
 *
 *
 *   @section  IteratorsAccessors 访问器
 * 迭代器就像指针一样：它们可以被递增和递减，但它们实际上是很笨的。它们的神奇之处在于它们指向一些有用的对象，在这里是指访问器。对于指针来说，它们指向的是一个存储了一些数据的实际对象。另一方面，deal.II迭代器，当被解除引用时，并不返回对一个实际对象的引用，而是返回一个知道如何获取代表单元格的数据的对象。一般来说，这个对象本身并不存储单元格的顶点或其邻居是什么。然而，它知道如何从Triangulation类为描述网格而设置的数组、表格和列表中获取这类信息。
 * 访问表征一个单元的数据总是通过Accessor完成的，即表达式
 * <code>i-&gt;xxx()</code>
 * 允许访问这个Accessor的<b>all</b>属性。你可以从迭代器中查询的属性的例子有
 * @code
 * cell->vertex(1);
 * line->child(0);
 * hex->face(3);
 * cell->at_boundary();
 * face->boundary_id();
 * @endcode
 *
 * 由于对迭代器的解引用产生了访问器对象，这些调用是Tomember函数
 * <code>Accessor::vertex()</code>  ,  <code>Accessor::child()</code>
 * 等。这些函数依次从存储这些数据的各种数据结构中找出相关的数据。这实际上是如何做到的，以及使用什么数据结构，并不是交易.II中应用的作者真正关心的。特别是，通过隐藏实际的数据结构，我们能够以一种有效的方式存储数据，而不一定是以一种使应用程序编写者容易访问或理解的方式。
 *
 *
 * @section  IteratorsTypedefs 访问器的种类
 * 根据你要访问的数据种类，有不同的访问器类。
 *
 * - TriaAccessor类为你提供数据，识别构成三角形的单元、面、线、四边形和六边形的几何属性，以及父子关系。
 *
 * - CellAccessor类是从TriaAccessor类派生出来的，用于对象具有完整维度的情况，即是一个单元，而不是例如一个单元的边界线。在这种情况下，关于一个网格的拓扑连接的额外信息可以从一个访问器中获得，比如请求指向一个单元的邻居的迭代器。
 *
 * - DoFAccessor类可以让你访问与单元格、面等相关的自由度信息；它对DoFHandler和 hp::DoFHandler 对象都是如此。请注意，DoFAccessor类派生于TriaAccessor或CellAccessor（取决于DoFAccessor是否指向全维度的对象），因此能够提供比其基类更高的信息集。此外，DoFAccessor类有两种风格，一种是访问细胞层面的自由度，另一种是访问活动细胞的活动度。
 *
 * - DoFCellAccessor类与DoFCellAccessor的目的和关系与CellAccessor与TriaAccessor的目的和关系相同。
 * 除了查找成员文档，你通常不需要处理上面列出的实际类名。相反，我们使用由网格类Triangulation、DoFHandler和
 * hp::DoFHandler,
 * 提供的类型定义，以及生成此类对象的函数。 <table
 * border=1> <tr> <th>Class</th> <th>cell_iterator type</th> <th>function
 * call</th> </tr> <tr> <th>Triangulation</th> <td>typename
 * Triangulation::cell_iterator</td>   <td>Triangulation::begin()</td>  </tr>
 * <tr> <th>DoFHandler</th> <td>typename  DoFHandler::cell_iterator</td>
 * <td>DoFHandler::begin()</td>  </tr> <tr>  <th>hp::DoFHandler</th>
 * <td>typename  hp::DoFHandler::cell_iterator</td>
 * <td>hp::DoFHandler::begin()</td>  </tr></table>
 * Triangulation类支持用<tt>typename  Triangulation::face_iterator</tt>,
 * 在单元格面之间进行迭代，该类型由
 * Triangulation::begin_face().  返回。 活动迭代器有以下属性。
 * <table border=1> <tr> <th>Class</th> <th>cell_iterator type</th>
 * <th>function call</th> </tr> <tr> <th>Triangulation</th> <td>typename
 * Triangulation::active_cell_iterator</td>
 * <td>Triangulation::begin_active()</td>  </tr> <tr> <th>DoFHandler</th>
 * <td>typename  DoFHandler::active_cell_iterator</td>
 * <td>DoFHandler::begin_active()</td>  </tr> <tr>  <th>hp::DoFHandler</th>
 * <td>typename  hp::DoFHandler::active_cell_iterator</td>
 * <td>hp::DoFHandler::begin_active()</td>  </tr></table>
 * Triangulation类也支持用<tt>typename
 * Triangulation::active_face_iterator</tt>,
 * 遍历活动单元面，这是由 Triangulation::begin_active_face().
 * 返回的类型。
 * 除了这些作用于单元和面的类型和调用（取决于维度的逻辑概念：一个单元在2D中是一个四边形，但在3D中是一个六面体），还有相应的类型和调用，如
 * <code>begin_active_quad()</code> or <code>end_quad()</code>
 * ，作用于独立于维度的几何对象line, quad, and
 * hex。这些调用，就像上面那些调用一样，以活动和非活动的形式存在。
 * 所有Mesh类中的类型定义都是在 <code>begin_active_quad()</code>
 * or <code>end_quad()</code> 中说明的。
 *
 * -  internal::Triangulation::Iterators<1,spacedim>,   internal::Triangulation::Iterators<2,spacedim>,  和  internal::Triangulation::Iterators<3,spacedim>  三角形迭代器的类。
 *
 * - 用于DoFHandler和 hp::DoFHandler 迭代器的<a class="el"
 * href="structinternal_1_1DoFHandler_1_1Iterators_3_01DoFHandlerType_3_011_00_01spacedim_01_4_00_01lda_01_4.html">internal::DoFHandler::Iterators&lt;DoFHandlerType&lt;1,spacedim&gt;,
 * lda&gt;</a>、<a class="el"
 * href="structinternal_1_1Triangulation_1_1Iterators_3_012_00_01spacedim_01_4.html">internal::DoFHandler::Iterators&lt;DoFHandlerType&lt;1,spacedim&gt;,
 * lda&gt;</a>、<a class="el"
 * href="structinternal_1_1Triangulation_1_1Iterators_3_013_00_01spacedim_01_4.html">internal::DoFHandler::Iterators&lt;DoFHandlerType&lt;1,spacedim&gt;,
 * lda&gt;</a>类。 @section  IteratorAccessorInternals
 * 迭代器和访问器的内部结构
 * 迭代器，就像指针一样，就像它们指向一个实际的对象一样，但实际上它们所做的只是在被引用时返回一个访问器。访问器对象包含状态，也就是说，它知道它所代表的对象，例如，通过存储它属于哪个三角形，以及单元格中的级别和索引。因此，它能够访问与它所代表的单元（或面，或边）相对应的数据。
 * 有一个过去-末端指针的表示，由TriaAccessor类中的成员变量
 * <code>present_level</code> and <code>present_index</code>
 * 的特殊值表示。如果 <code>present_level</code> @> =0 and
 * <code>present_index</code>   @>  =0，那么该对象是有效的；如果
 * <code>present_level</code>==-1 and <code>present_index</code>
 * ==-1，那么该迭代器指向了终点；在所有其他情况下，该迭代器被认为是无效的。你可以通过调用
 * TriaAccessorBase::state()  函数来检查。
 * 当向后运行时，过端迭代器也可用于比较迭代器与开始前的值。指向向量两端的迭代器之间是没有区别的。
 * 单元的存储是基于层次结构的，因此上面提到的结构是有用的。然而，面不是按层次组织的，低维度的对象的访问器没有
 * <code>present_level</code> 成员变量。
 *
 * @ingroup grid
 *
 */

//@}


/**
 *    @defgroup Accessors Accessor classes of the mesh iterators
 * @ingroup Iterators
 *
 */


