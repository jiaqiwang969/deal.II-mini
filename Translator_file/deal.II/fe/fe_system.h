//include/deal.II-translator/fe/fe_system_0.txt
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

#ifndef dealii_fe_system_h
#  define dealii_fe_system_h


 /*----------------------------   fe_system.h     ---------------------------*/ 


#  include <deal.II/base/config.h>

#  include <deal.II/base/thread_management.h>

#  include <deal.II/fe/fe.h>
#  include <deal.II/fe/fe_tools.h>

#  include <memory>
#  include <type_traits>
#  include <utility>
#  include <vector>


DEAL_II_NAMESPACE_OPEN

// Forward declaration
#  ifndef DOXYGEN
template <int dim, int spacedim>
class FE_Enriched;
#  endif

/**
 * 这个类提供了一个接口，将几个元素组合成一个，矢量值的元素。作为例子，考虑用于解决斯托克斯和纳维-斯托克斯方程的泰勒-霍德元素。在这里，速度（其分量与域的维数
 * $d$ 一样多）用 $Q_2$ 元素离散，压力用 $Q_1$
 * 元素。在数学上，耦合问题的有限元空间通常被写成 $V_h
 * = Q_2^d \times Q_1$ ，其中指数化被理解为空间的张量乘积
 *
 * - 即在2d中，我们有 $V_h=Q_2\times Q_2\times Q_1$
 *
 * - 而张量积导致向量，其中向量值函数空间的每个分量对应于 $Q_2$ 或 $Q_1$ 空间中的一个标量函数。使用FESystem类，这个空间的创建方法是
 *
 * @code
 * FESystem<dim> taylor_hood_fe (FE_Q<dim>(2)^dim,   // velocity components
 *                               FE_Q<dim>(1));      // pressure component
 * @endcode
 * 在这里创建这个元素相当于在FESystem构造函数的参数列表的第一行中对
 * $Q_2$
 * 元素进行张量乘法，然后通过另一个张量乘法与第二行中的元素进行串联。例如，在
 * step-22 的教程程序中就使用了这种构造。 同样， step-8
 * 解决了一个弹性方程，我们需要解决一个固体物体的位移。如果域是
 * $d$ -维的，位移又有 $d$
 * 分量，因此，组合有限元是用以下方法创建的
 *
 * @code
 * FESystem<dim> displacement_fe (FE_Q<dim>(1)^dim);
 * @endcode
 * 现在组合元素的每个（矢量）分量都对应于一个 $Q_1$
 * 空间。
 * 对外界来说，FESystem对象看起来就像一个普通的有限元对象，它们只是碰巧由其他几个可能是不同类型的有限元组成。这些
 * "基元 "本身可以有多个分量，特别是也可以是矢量值的。
 *
 * - 例如，如果其中一个基元是一个FESystem本身（也见下文）。在命名空间 FETools::Compositing, 的文档中给出了一个使用 "张量乘积 "策略的例子。
 * %矢量值元素在一些教程程序中讨论，例如 step-8
 * ,  step-20  ,  step-21  ,  step-22  ，特别是在 @ref vector_valued
 * 模块中。 @dealiiVideoLecture{19,20}   <h3>FESystem, components and
 * blocks</h3>
 * 一个FESystem，除了最微不足道的情况，会产生一个有多个分量的向量值的有限元。分量的数量n_components()对应于PDE系统中解函数的维度，相应地也对应于你的PDE系统的方程数量。例如， step-20 中涉及的混合拉普拉斯系统在 $d$ 空间维度上有 $d+1$ 个分量：标量压力和速度矢量的 $d$ 分量。同样， step-8 中涉及的弹性方程在 $d$ 空间维度上有 $d$ 分量。一般来说，FES系统元素的分量数是所有基础元素的分量累积数乘以其倍数。在 @ref GlossComponent 的 "元件术语条目 "
 * 中也给出了更多关于元件的信息。
 * 虽然从偏微分方程的角度来看，分量的概念很重要，但在有限元方面看起来有点不同，因为不仅是FESystem，还有像FE_RaviartThomas这样的矢量值元素，都有几个分量。这里需要的概念是一个 @ref GlossBlock "块"
 * 。每个块包含了与FESystem的单个基元相关的自由度集合，其中具有倍数的基元要计算多次。这些块通常使用
 * DoFHandler::block_info().
 * 中的信息来处理，一个FESystem对象的块数只是所有基元的倍数之和，由n_blocks()给出。
 * 例如，用于三维斯托克斯问题的Taylor-Hood元素的FESystem可以用以下代码建立
 *
 * @code
 * const FE_Q<3> u(2);
 * const FE_Q<3> p(1);
 * FESystem<3> sys1(u,3, p,1);
 * @endcode
 * 或者更简洁地通过
 *
 * @code
 * FESystem<3> sys1(FE_Q<3>(2),3,
 *                FE_Q<3>(1),1);
 * @endcode
 * 甚至更短（模仿数学符号，我们正在处理一个 $Q_2^3 \times
 * Q_1$ 元素）。
 *
 * @code
 * FESystem<3> sys1(FE_Q<3>(2)^3,
 *                FE_Q<3>(1));
 * @endcode
 *
 * 这个例子创建了一个FES系统 @p sys1
 * ，有四个元件，三个是速度元件，一个是压力元件，还有四个块，每个速度元件的自由度和压力都在一个单独的块中。由于第一个基元重复了三次，所以块的数量为四。
 * 另一方面，Taylor-Hood元素也可以用以下方法构造
 *
 *
 * @code
 * FESystem<3> U(u,3);
 * FESystem<3> sys2(U, p);
 * @endcode
 *
 * 这里创建的FES系统 @p sys2
 * 具有相同的四个部件，但自由度只分布在两个块中。第一个块有
 * @p U,
 * 的所有速度自由度，而第二个块包含压力自由度。请注意，虽然
 * @p U 本身有3个块，但FES系统 @p sys2 并不试图将 @p U
 * 分割成其基本元素，而是将其视为自己的一个块。通过像
 * @p sys2,
 * 那样先将所有速度封锁在一个系统中，我们实现了相同的块结构，如果我们不使用
 * $Q_2^3$
 * 元素来表示速度，而是使用矢量值的基本元素，例如使用达尔西定律的混合离散化，就会产生相同的块结构。
 *
 *
 * @code
 * FE_RaviartThomas<3> u(1);
 * FE_DGQ<3> p(1);
 * FESystem<3> sys3(u, p);
 * @endcode
 *
 * 这个例子也产生了一个有四个元件的系统，但只有两个块。
 * 在大多数情况下，组成元素的行为就像它是一个普通元素一样。它只是比大多数
 * "普通
 * "元素有更多的自由度。然而，底层结构在限制、延长和界面约束矩阵中是可见的，它们并不与基础元素的自由度相联系。例如，连续性要求是对子对象的形状函数分别施加的；不同子对象的形状函数之间不存在要求，即在上面的例子中：在一个悬挂的节点上，
 * @p u 速度的各自值只与 @p u
 * 的顶点和这个顶点旁边的大单元上的线耦合，但与这个或其他单元的
 * @p v 和 @p w 没有互动。
 *
 *  <h3>Internal information on numbering of degrees of freedom</h3>
 * 自由度的总体编号如下：对于每个子对象（顶点、线、四边形或六边形），自由度的编号是这样的：我们首先运行所有的子元素，然后再转向这个子对象的下一个自由度或下一个子对象。例如，对于一个在一个空间维度上有三个分量的元素，前两个分量是立方滞后元素，第三个是二次滞后元素，系统<tt>s=(u,v,p)</tt>的排序是。
 * <ul>   <li>  第一个顶点。<tt>u0, v0, p0 = s0, s1, s2</tt>  <li>  第二个顶点。<tt>u1, v1, p1 = s3, s4, s5</tt>  <li>  线上第一个分量。<tt>u2, u3 = s4, s5</tt>  <li>  线上的第二个分量。<tt>v2, v3 = s6, s7</tt>。  <li>  线路上的第三个分量。<tt>p2 = s8</tt>。  </ul>  也就是说，你不应该在你的应用程序中依赖这个编号，因为这些%的内部成员将来可能会改变。而是使用函数system_to_component_index()和component_to_system_index()。
 * 关于模板参数<tt>spacedim</tt>的更多信息，请参见Triangulation的文档。
 *
 *
 * @ingroup febase fe vector_valued
 *
 *
 *
 */
template <int dim, int spacedim = dim>
class FESystem : public FiniteElement<dim, spacedim>
{
public:
  /**
   * 删除默认构造函数，以便在没有提供FiniteElement的情况下，`FESystem(FEPairs
   * &&... fe_pairs)`不会被意外选中。
   *
   */
  FESystem() = delete;

  /**
   * 构造函数。取一个有限元和你想用这个类来组合的元素数量。
   * 这个对象 @p fe
   * 除了创建一个副本，然后被当前对象所拥有之外，实际上没有其他用途。换句话说，用一个临时的有限元对象来调用这个构造函数是完全可以的，就像在这个代码片段中。
   * @code
   * FESystem<dim> fe (FE_Q<dim>(2), 2);
   * @endcode
   * 这里， <code>FE_Q@<dim@>(2)</code> 构造了一个未命名的临时对象，传递给FESystem构造函数来创建一个由两个组件组成的有限元，这两个组件都是二次FE_Q元素。在这一行对应的代码结束时，这个临时对象又被销毁了，但这并不重要，因为FESystem创建了自己的FE_Q对象的副本。    这个构造函数（或其下面的变体）基本上在所有处理矢量值问题的教程程序中都会用到。用例见 step-8  ,  step-20  ,  step-22  等。也可参见  @ref vector_valued  "处理向量值问题 "
   * 模块。      @dealiiVideoLecture{19,20}   @param[in]  fe
   * 将被用来表示该组成元素的组件的有限元素。
   * @param[in]  n_elements
   * 一个整数，表示这个元素应该由多少份 @p fe 组成。
   *
   */
  FESystem(const FiniteElement<dim, spacedim> &fe,
           const unsigned int                  n_elements);

  /**
   * 用于具有两个基本元素的混合离散的构造函数。
   * 参见上面的另一个构造函数，以了解组成元素的一般想法。
   *
   */
  FESystem(const FiniteElement<dim, spacedim> &fe1,
           const unsigned int                  n1,
           const FiniteElement<dim, spacedim> &fe2,
           const unsigned int                  n2);

  /**
   * 用于具有三个基本元素的混合离散的构造函数。
   * 参见上面的另一个构造函数，以解释组成元素的一般想法。
   *
   */
  FESystem(const FiniteElement<dim, spacedim> &fe1,
           const unsigned int                  n1,
           const FiniteElement<dim, spacedim> &fe2,
           const unsigned int                  n2,
           const FiniteElement<dim, spacedim> &fe3,
           const unsigned int                  n3);

  /**
   * 用于具有四个基本元素的混合离散的构造函数。
   * 参见上述其他构造函数中的第一个，以获得对组成元素的一般概念的解释。
   *
   */
  FESystem(const FiniteElement<dim, spacedim> &fe1,
           const unsigned int                  n1,
           const FiniteElement<dim, spacedim> &fe2,
           const unsigned int                  n2,
           const FiniteElement<dim, spacedim> &fe3,
           const unsigned int                  n3,
           const FiniteElement<dim, spacedim> &fe4,
           const unsigned int                  n4);

  /**
   * 五个基本元素的混合离散的构造函数。
   * 关于组成元素的一般概念，请看上面第一个其他构造函数的解释。
   *
   */
  FESystem(const FiniteElement<dim, spacedim> &fe1,
           const unsigned int                  n1,
           const FiniteElement<dim, spacedim> &fe2,
           const unsigned int                  n2,
           const FiniteElement<dim, spacedim> &fe3,
           const unsigned int                  n3,
           const FiniteElement<dim, spacedim> &fe4,
           const unsigned int                  n4,
           const FiniteElement<dim, spacedim> &fe5,
           const unsigned int                  n5);

  /**
   * 与上述相同，但适用于任何数量的基元。指向基数元素的指针和它们的乘数被作为向量传递给这个构造函数。这些向量的长度被认为是相等的。
   * 如上所述，第一个参数所指向的有限元对象除了在内部创建副本外，实际上并不使用。因此，你可以在调用这个构造函数后立即再次删除这些指针。
   * <h4>How to use this constructor</h4>
   * 使用这个构造函数有时会有点尴尬，因为你需要在一个地方传递两个向量，而这个地方可能无法直接构造这样一个向量
   *
   * - 例如，在一个有FESystem成员变量的类的成员初始化列表中。例如，如果你的主类看起来像这样。
   * @code
   * template <int dim>
   * class MySimulator {
   * public:
   *   MySimulator (const unsigned int polynomial_degree);
   * private:
   *   FESystem<dim> fe;
   * };
   *
   * template <int dim>
   * MySimulator<dim>::MySimulator (const unsigned int polynomial_degree)
   *   :
   *   fe (...)  // what to pass here???
   * {}
   * @endcode
   * 使用C++11语言标准（或更高版本），你可以这样做来创建一个具有四个基元和乘数1、2、3、4的元素。
   * @code
   * template <int dim>
   * MySimulator<dim>::MySimulator (const unsigned int polynomial_degree)
   *   :
   *   fe (std::vector<const FiniteElement<dim>*> { new FE_Q<dim>(1),
   *                                                new FE_Q<dim>(2),
   *                                                new FE_Q<dim>(3),
   *                                                new FE_Q<dim>(4) },
   *       std::vector<unsigned int> { 1, 2, 3, 4 })
   * {}
   * @endcode
   * 这就在原地创建了两个向量，并使用大括号中的初始化器列表对它们进行初始化
   * <code>{ ... }</code>  。
   * 这段代码有一个问题：它产生了四个内存泄漏，因为上面的第一个向量是用指向
   * <code>new</code> 分配的元素的指针创建的，但从未销毁。
   * 解决其中第二个问题的方法是创建两个可以创建向量的静态成员函数。下面是一个例子。
   * @code
   * template <int dim>
   * class MySimulator {
   * public:
   *   MySimulator (const unsigned int polynomial_degree);
   *
   * private:
   *   FESystem<dim> fe;
   *
   *   static std::vector<const FiniteElement<dim>*>
   *   create_fe_list (const unsigned int polynomial_degree);
   *
   *   static std::vector<unsigned int>
   *   create_fe_multiplicities ();
   * };
   *
   * template <int dim>
   * std::vector<const FiniteElement<dim>*>
   * MySimulator<dim>::create_fe_list (const unsigned int polynomial_degree)
   * {
   *   std::vector<const FiniteElement<dim>*> fe_list;
   *   fe_list.push_back (new FE_Q<dim>(1));
   *   fe_list.push_back (new FE_Q<dim>(2));
   *   fe_list.push_back (new FE_Q<dim>(3));
   *   fe_list.push_back (new FE_Q<dim>(4));
   *   return fe_list;
   * }
   *
   * template <int dim>
   * std::vector<unsigned int>
   * MySimulator<dim>::create_fe_multiplicities ()
   * {
   *   std::vector<unsigned int> multiplicities;
   *   multiplicities.push_back (1);
   *   multiplicities.push_back (2);
   *   multiplicities.push_back (3);
   *   multiplicities.push_back (4);
   *   return multiplicities;
   * }
   *
   * template <int dim>
   * MySimulator<dim>::MySimulator (const unsigned int polynomial_degree)
   *   :
   *   fe (create_fe_list (polynomial_degree),
   *       create_fe_multiplicities ())
   * {}
   * @endcode
   * 这样做的方法是，我们有两个静态成员函数来创建必要的向量，以传递给成员变量的构造函数
   * <code>fe</code>  。它们需要是静态的，因为它们是在
   * <code>MySimulator</code>  的构造函数中被调用的，此时
   * <code>*this</code>
   * 对象还没有完全构造好，因此，常规的成员函数还不能被调用。
   * 不过上面的代码还没有解决内存泄漏的问题：
   * <code>create_fe_list()</code>
   * 函数创建了一个指针向量，但没有任何东西破坏这些指针。这就是解决方案。
   * @code
   * template <int dim>
   * class MySimulator
   * {
   * public:
   * MySimulator (const unsigned int polynomial_degree);
   *
   * private:
   * FESystem<dim> fe;
   *
   * struct VectorElementDestroyer
   * {
   *   const std::vector<const FiniteElement<dim>*> data;
   *
   *   VectorElementDestroyer(
   *     const std::vector<const FiniteElement<dim>*> &pointers);
   *
   *    // destructor to delete the pointers
   *   ~VectorElementDestroyer ();
   *
   *   const std::vector<const FiniteElement<dim>*> & get_data () const;
   * };
   *
   * static std::vector<const FiniteElement<dim>*>
   * create_fe_list (const unsigned int polynomial_degree);
   *
   * static std::vector<unsigned int>
   * create_fe_multiplicities ();
   * };
   *
   * template <int dim>
   * MySimulator<dim>::VectorElementDestroyer::
   * VectorElementDestroyer(
   * const std::vector<const FiniteElement<dim>*> &pointers)
   * :
   * data(pointers)
   * {}
   *
   * template <int dim>
   * MySimulator<dim>::VectorElementDestroyer::
   * ~VectorElementDestroyer ()
   * {
   * for (unsigned int i=0; i<data.size(); ++i)
   *   delete data[i];
   * }
   *
   * template <int dim>
   * const std::vector<const FiniteElement<dim>*> &
   * MySimulator<dim>::VectorElementDestroyer::
   * get_data () const
   * {
   * return data;
   * }
   *
   * template <int dim>
   * MySimulator<dim>::MySimulator (const unsigned int polynomial_degree)
   * :
   * fe (VectorElementDestroyer(create_fe_list (polynomial_degree)).get_data(),
   *   create_fe_multiplicities ())
   * {}
   * @endcode
   * 换句话说，我们从 <code>create_fe_list()</code>
   * 中收到的向量被打包到一个 <code>VectorElementDestroyer</code>
   * 类型的临时对象中；然后我们立即从这个临时对象中获得向量，将其传递给
   * <code>fe</code>; and finally, the <code>VectorElementDestroyer</code>
   * 的构造器，析构器在整个表达式的最后被调用（在
   * <code>fe</code>
   * 的构造器完成后），并销毁了临时向量的元素。瞧：既不短也不优雅，但它是有效的!
   *
   */
  FESystem(const std::vector<const FiniteElement<dim, spacedim> *> &fes,
           const std::vector<unsigned int> &multiplicities);

#  if !defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1900
  /**
   * 构造函数接受任意数量的参数，类型为
   * <code>std::pair<std::unique_ptr<FiniteElement<dim,  spacedim>>, unsigned
   * int></code>。与 FiniteElement::operator^,
   * 相结合，可以构造如下的FESystem对象。
   * @code
   * FiniteElementType1<dim,spacedim> fe_1;
   * FiniteElementType1<dim,spacedim> fe_2;
   * FESystem<dim,spacedim> fe_system ( fe_1^dim, fe_2 );
   * @endcode
   * `fe_1`和`fe_2`对象除了创建一个副本，然后被当前对象所拥有之外，实际上没有其他用途。换句话说，用一个临时的有限元对象来调用这个构造函数是完全可以的，就像这个代码片断一样。
   * @code
   * FESystem<dim> fe (FE_Q<dim>(2)^2);
   * @endcode
   * 这里， <code>FE_Q@<dim@>(2)</code>
   * 构造了一个未命名的临时对象，传递给FESystem构造函数来创建一个由两个组件组成的有限元，这两个组件都是二次FE_Q元素。在这一行对应的代码结束时，这个临时对象又被销毁了，但这并不重要，因为FESystem创建了自己的FE_Q对象的副本。
   * 作为一种快捷方式，这个构造函数也允许调用
   * @code
   * FESystem<dim> fe (FE_Q<dim>(2)^dim, FE_Q<dim>(1));
   * @endcode
   * 而不是更明确的
   * @code
   * FESystem<dim> fe (FE_Q<dim>(2)^dim, FE_Q<dim>(1)^1);
   * @endcode
   * 换句话说，如果没有通过指数化操作明确指定一个元素的倍数，那么就假定它是1（正如人们所期望的）。
   * @warning
   * 这个功能在19.0版之前的英特尔编译器上是不可用的。对于18.0之前的英特尔编译器，定义这个构造函数会导致内部编译器错误。
   *
   */
  template <
    class... FEPairs,
    typename = typename enable_if_all<
      (std::is_same<typename std::decay<FEPairs>::type,
                    std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
                              unsigned int>>::value ||
       std::is_base_of<FiniteElement<dim, spacedim>,
                       typename std::decay<FEPairs>::type>::value)...>::type>
  FESystem(FEPairs &&... fe_pairs);

  /**
   * 与上述相同，允许采用以下语法。
   * @code
   * FiniteElementType1<dim,spacedim> fe_1;
   * FiniteElementType1<dim,spacedim> fe_2;
   * FESystem<dim,spacedim> fe_system = { fe_1^dim, fe_2^1 };
   * @endcode
   * @warning
   * 这个功能对19.0版以前的英特尔编译器是不可用的。构造函数只是没有被选中进行重载解析。
   *
   */
  FESystem(
    const std::initializer_list<
      std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>>
      &fe_systems);
#  endif

  /**
   * 复制构造函数。该构造函数被删除，即不允许复制FESystem对象。
   *
   */
  FESystem(const FESystem<dim, spacedim> &) = delete;

  /**
   * 移动构造函数。
   *
   */
  FESystem(FESystem<dim, spacedim> &&other_fe_system) noexcept
    : FiniteElement<dim, spacedim>(std::move(other_fe_system))
  {
    base_elements = std::move(other_fe_system.base_elements);
    generalized_support_points_index_table =
      std::move(other_fe_system.generalized_support_points_index_table);
  }

  /**
   * 解构器。
   *
   */
  virtual ~FESystem() override = default;

  /**
   * 返回一个唯一标识一个有限元素的字符串。该元素返回一个字符串，该字符串由基础元素返回的字符串
   * @p name1...@p
   * nameN组成。从这些中，我们创建一个序列<tt>FESystem<dim>[name1^m1-name2^m2-...-nameN^mN]</tt>，其中
   * @p mi
   * 是基元的倍率。如果一个倍数等于1，那么上标就省略了。
   *
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  // make variant with ComponentMask also available:
  using FiniteElement<dim, spacedim>::get_sub_fe;

  /**
   * @copydoc   FiniteElement<dim,spacedim>::get_sub_fe() 。
   *
   */
  virtual const FiniteElement<dim, spacedim> &
  get_sub_fe(const unsigned int first_component,
             const unsigned int n_selected_components) const override;

  /**
   * 返回 @p ith 点的形状函数值  @p p.   @p p
   * 是参考元素上的一个点。由于这个有限元总是矢量值的，我们返回这个形状函数的矢量值的唯一非零分量的值。如果形状函数有一个以上的非零分量（我们用非原始分量来指代），那么抛出一个
   * @p ExcShapeFunctionNotPrimitive. 类型的异常 如果 @p FiniteElement
   * 的形状值（对应于 @p ith
   * 形状函数）取决于实空间中的单元的形状，就会抛出 @p
   * ExcUnitShapeValuesDoNotExist 。
   *
   */
  virtual double
  shape_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 返回 @p componentth 形状函数的 @p ith 矢量分量在 @p p.
   * 点的值，关于这个函数的语义，请看FiniteElement基类的更多信息。
   * 由于这个元素一般都是矢量值，所以它把这些值的计算转交给基元素。
   *
   */
  virtual double
  shape_value_component(const unsigned int i,
                        const Point<dim> & p,
                        const unsigned int component) const override;

  /**
   * 返回 @p ith 形状函数在 @p p. 点的梯度 @p p
   * 是参考元素上的一个点，同样，梯度是单元格上关于单元格坐标的梯度。由于这个有限元总是矢量值的，我们返回这个形状函数的矢量值的唯一非零分量的值。如果形状函数有一个以上的非零分量（我们用非原始分量一词来指代），那么就抛出一个类型为
   * @p  ExcShapeFunctionNotPrimitive的异常。    如果 @p FiniteElement
   * 的形状值（对应于 @p ith
   * 形状函数）取决于实空间中的单元格形状，则抛出一个
   * @p ExcUnitShapeValuesDoNotExist 。
   *
   */
  virtual Tensor<1, dim>
  shape_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 返回 @p componentth 形状函数的 @p ith 矢量分量在 @p p.
   * 点的梯度，关于这个函数的语义，请看FiniteElement基类的更多信息。
   * 由于这个元素一般都是矢量值，它把这些值的计算转给了基元素。
   *
   */
  virtual Tensor<1, dim>
  shape_grad_component(const unsigned int i,
                       const Point<dim> & p,
                       const unsigned int component) const override;

  /**
   * 返回 @p ith 形状函数在单元格上 @p p
   * 点的二次导数的张量。该导数是单元格上相对于单元格坐标的导数。由于这个有限元总是矢量值的，我们返回这个形状函数的矢量值的唯一非零分量的值。如果形状函数有一个以上的非零分量（我们用非原始分量一词来指代），那么抛出一个类型为
   * @p  ExcShapeFunctionNotPrimitive的异常。    如果 @p FiniteElement
   * 的形状值（对应于 @p ith
   * 形状函数）取决于实空间中的单元格形状，则抛出 @p
   * ExcUnitShapeValuesDoNotExist 。
   *
   */
  virtual Tensor<2, dim>
  shape_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * 返回 @p componentth 形状函数的 @p ith 矢量分量在 @p p.
   * 点的二阶导数，关于这个函数的语义，请看FiniteElement基类的更多信息。
   * 由于这个元素一般都是矢量值，它把这些值的计算转给了基元素。
   *
   */
  virtual Tensor<2, dim>
  shape_grad_grad_component(const unsigned int i,
                            const Point<dim> & p,
                            const unsigned int component) const override;

  /**
   * 返回 @p ith 形状函数在单元格上 @p p
   * 点的三次导数的张量。这些导数是单元格上相对于单元格坐标的导数。由于这个有限元总是矢量值的，我们返回这个形状函数的矢量值的唯一非零分量的值。如果形状函数有一个以上的非零分量（我们用非原始分量一词来指代），那么抛出一个
   * @p  ExcShapeFunctionNotPrimitive类型的异常。    如果 @p
   * FiniteElement 的形状值（对应于 @p ith
   * 形状函数）取决于实空间中的单元格形状，则抛出 @p
   * ExcUnitShapeValuesDoNotExist 。
   *
   */
  virtual Tensor<3, dim>
  shape_3rd_derivative(const unsigned int i,
                       const Point<dim> & p) const override;

  /**
   * 返回 @p componentth 形状函数的 @p ith 矢量分量在点 @p p.
   * 处的三阶导数
   * 关于此函数的语义，请参见FiniteElement基类。
   * 由于这个元素一般都是矢量值，它把这些值的计算转给了基元素。
   *
   */
  virtual Tensor<3, dim>
  shape_3rd_derivative_component(const unsigned int i,
                                 const Point<dim> & p,
                                 const unsigned int component) const override;

  /**
   * 返回 @p ith 形状函数在单元格上 @p p
   * 点的第四导数的张量。这些导数是单元格上关于单元格坐标的导数。由于这个有限元总是矢量值的，我们返回这个形状函数的矢量值的唯一非零分量的值。如果形状函数有一个以上的非零分量（我们用非原始分量一词来指代），那么抛出一个
   * @p  ExcShapeFunctionNotPrimitive类型的异常。    如果 @p
   * FiniteElement 的形状值（对应于 @p ith
   * 形状函数）取决于实空间中的单元格的形状，则抛出 @p
   * ExcUnitShapeValuesDoNotExist 。
   *
   */
  virtual Tensor<4, dim>
  shape_4th_derivative(const unsigned int i,
                       const Point<dim> & p) const override;

  /**
   * 返回 @p componentth 形状函数的 @p ith 矢量分量在点 @p p.
   * 处的四次导数，关于此函数的语义，请参见FiniteElement基类。
   * 因为这个元素一般都是矢量值，所以它把这些值的计算转交给基元素。
   *
   */
  virtual Tensor<4, dim>
  shape_4th_derivative_component(const unsigned int i,
                                 const Point<dim> & p,
                                 const unsigned int component) const override;

  /**
   * 返回从给定的有限元内插到现在的矩阵。然后矩阵的大小是
   * @p dofs_per_cell 乘以<tt>source.n_dofs_per_cell()</tt>。
   * 如果源元素和目的元素都是 @p FESystem
   * 元素，有相同数量的基本元素，有相同的元素倍数，并且这些基本元素也实现了它们的
   * @p
   * get_interpolation_matrix函数，这些矩阵就可以使用。否则，会抛出一个
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented
   * 类型的异常。
   *
   */
  virtual void
  get_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                           FullMatrix<double> &matrix) const override;

  /**
   * 访问一个合成元素。索引需要小于基数元素的数量。请注意，如果乘数大于1，基础元素的数量可能反过来小于系统元素的组件数量。
   *
   */
  virtual const FiniteElement<dim, spacedim> &
  base_element(const unsigned int index) const override;

  /**
   * 如果形状函数 @p shape_index 在面 @p face_index.
   * 的某处有非零函数值，该函数返回 @p true, 。
   *
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  /**
   * 从精细网格空间投射到粗略网格空间。重写FiniteElement中的相应方法，实现懒人评估（在请求时初始化）。
   * 如果这个投影运算符与一个矩阵 @p P,
   * 相关联，那么这里将返回这个矩阵 @p P_i
   * 对一个子单元的限制。    矩阵 @p P 是单元格矩阵 @p
   * P_i的连接或相加，取决于
   * FiniteElement::restriction_is_additive().
   * 的值，这区分了插值（连接）和标量积（相加）方面的投影。
   * 行和列指数分别与粗网格和细网格空间有关，与相关运算符的定义一致。
   * 如果投影矩阵没有在派生的有限元类中实现，这个函数会以
   * FiniteElement::ExcProjectionVoid.
   * 类型的异常中止，你可以通过首先调用restriction_is_implemented()或isotropic_restriction_is_implemented()函数检查是否会发生这种情况。
   *
   */
  virtual const FullMatrix<double> &
  get_restriction_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * 网格间的嵌入矩阵。重写FiniteElement中的相应方法，实现懒人评估（查询时初始化）。
   * 从粗网格空间到细网格空间的身份运算符与一个矩阵 @p
   * P. 相关联，该矩阵 @p P_i
   * 对单个子单元的限制在这里被返回。    矩阵 @p P
   * 是串联的，而不是单元格矩阵 @p
   * P_i的总和。也就是说，如果同一个非零条目<tt>j,k</tt>存在于两个不同的子矩阵
   * @p P_i,
   * 中，其值在两个矩阵中应该是相同的，它只被复制到矩阵
   * @p P 中一次。
   * 行和列指数分别与细格和粗格空间相关，与相关运算符的定义一致。
   * 这些矩阵被组装多级方法的延长矩阵的程序所使用。
   * 在使用这个矩阵阵列组装单元间的转移矩阵时，延长矩阵中的零元素被丢弃，不会填满转移矩阵。
   * 如果延长矩阵没有在一个基本的有限元类中实现，这个函数会以
   * FiniteElement::ExcEmbeddingVoid.
   * 类型的异常中止。你可以通过首先调用prolongation_is_implemented()或isotropic_prolongation_is_implemented()函数来检查是否会发生。
   *
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * 给出一个面的指数自然排序中的指数，返回单元格上相同自由度的指数。
   * 为了解释这个概念，考虑这样的情况：我们想知道一个面的自由度，例如作为FESystem元素的一部分，是否是原始的。不幸的是，FiniteElement类中的is_primitive()函数需要一个单元格索引，所以我们需要找到对应于当前面的索引的形状函数的单元格索引。
   * 这个函数可以做到这一点。
   * 实现这一点的代码将看起来像这样。
   * @code
   * for (i=0; i<dofs_per_face; ++i)
   * if (fe.is_primitive(fe.face_to_cell_index(i, some_face_no)))
   * ... do whatever
   * @endcode
   * 这个函数需要额外的参数，以考虑到实际的面可以是相对于所考虑的单元格的标准排序，或者可以是翻转的，定向的，等等。
   * @param  face_dof_index
   * 一个面的自由度的索引。这个指数必须在零和每个面的自由度之间。
   * @param  face
   * 这个自由度所在的面的编号。这个数字必须介于零和
   * GeometryInfo::faces_per_cell.   @param  face_orientation
   * 描述面的方向的一个部分。见  @ref GlossFaceOrientation  。
   * @param  face_flip 对脸部方向的描述的一部分。参见  @ref
   * GlossFaceOrientation  。    @param  face_rotation
   * 描述脸部方向的一部分。见  @ref GlossFaceOrientation  。
   * @return
   * 这个自由度在整个单元上的自由度集合中的索引。返回值将介于0和dofs_per_cell之间。
   *
   */
  virtual unsigned int
  face_to_cell_index(const unsigned int face_dof_index,
                     const unsigned int face,
                     const bool         face_orientation = true,
                     const bool         face_flip        = false,
                     const bool         face_rotation = false) const override;

  /**
   * 在基类中实现相应的函数。
   *
   */
  virtual Point<dim>
  unit_support_point(const unsigned int index) const override;

  /**
   * 在基类中实现相应的函数。
   *
   */
  virtual Point<dim - 1>
  unit_face_support_point(const unsigned int index,
                          const unsigned int face_no = 0) const override;

  /**
   * 返回一个元素的常量模式列表。返回表有多少行，就有多少个元素中的元件和dofs_per_cell列。对于有限元的每个分量，返回表中的行包含该元上常数函数1的基础表示。将每个基元的常数模式串联起来。
   *
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  /**
   * @name 支持hp的函数  
     * @{ 
   *
   */

  /**
   * 返回该元素是否以新的方式实现其悬挂的节点约束，这必须用于使元素
   * "hp-兼容"。    当且仅当其所有基础元素都返回 @p true
   * 时，此函数才会返回 @p true 。
   *
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * 返回从一个元素的面插值到邻近元素的面的矩阵。
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。
   * 这个元素的基础元素将不得不实现这个功能。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现给定元素的插值，那么他们必须抛出一个类型为
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented,
   * 的异常，这个异常将从这个元素传播出去。
   *
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double> &                matrix,
                                const unsigned int face_no = 0) const override;


  /**
   * 返回从一个元素的面插值到邻近元素的子面的矩阵。
   * 矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。
   * 这个元素的基础元素将不得不实现这个功能。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现给定元素的插值，那么他们必须抛出一个类型为
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented,
   * 的异常，这个异常将从这个元素传播出去。
   *
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source,
    const unsigned int                  subface,
    FullMatrix<double> &                matrix,
    const unsigned int                  face_no = 0) const override;

  /**
   * 如果在一个顶点上，有几个有限元被激活，hp-code首先为这些FEs的自由度分配不同的全局索引。然后调用这个函数来找出其中哪些应该得到相同的值，从而可以得到相同的全局自由度指数。
   * 因此，该函数返回当前有限元对象的自由度与 @p fe_other,
   * 的自由度之间的相同性列表，后者是对代表在该特定顶点上活动的其他有限元之一的有限元对象的引用。该函数计算两个有限元对象的哪些自由度是等价的，两个自由度的编号都在零和两个有限元的n_dofs_per_vertex()的相应值之间。每一对的第一个索引表示本元素的一个顶点自由度，而第二个是另一个有限元素的相应索引。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理线上自由度。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * 与hp_vertex_dof_indices()相同，只是该函数处理四边形上的自由度。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  /**
   * @copydoc   FiniteElement::compare_for_domination() .
   *
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

  //@}

  /**
   * FiniteElement::convert_generalized_support_point_values_to_dof_values()
   * 函数的实现。    这个函数简单地调用
   * FiniteElement::convert_generalized_support_point_values_to_dof_values
   * 的基本元素，并将所有内容重新组合到输出参数中。如果一个基元是非插值的，那么相应的dof值将用
   * "信号 "NaN来代替。
   * 如果FES系统的基本元素没有一个是插值的，则该函数失败。
   *
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              dof_values) const override;

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。
   * 这个函数是虚拟的，因为有限元对象通常是通过指向其基类的指针来访问的，而不是类本身。
   *
   */
  virtual std::size_t
  memory_consumption() const override;

protected:
  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim> &       quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  using FiniteElement<dim, spacedim>::get_face_data;

  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_face_data(
    const UpdateFlags               update_flags,
    const Mapping<dim, spacedim> &  mapping,
    const hp::QCollection<dim - 1> &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_subface_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim - 1> &   quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  using FiniteElement<dim, spacedim>::fill_fe_face_values;

  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const hp::QCollection<dim - 1> &                            quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

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
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  /**
   * 为三个<tt>fill_fe*_values</tt>函数做工作。
   * 调用（除其他外）<tt>fill_fe_([sub]face)_values</tt>的基础元素。如果<tt>face_no==invalid_face_no</tt>和<tt>sub_no==invalid_face_no</tt>，调用
   * @p fill_fe_values ；调用 @p fill_fe_face_values  ]
   * 如果<tt>face_no==invalid_face_no</tt>和<tt>sub_no!=invalid_face_no</tt>；如果<tt>face_no!=invalid_face_no</tt>和<tt>sub_no!
   *
   */
  template <int dim_1>
  void
  compute_fill(
    const Mapping<dim, spacedim> &                              mapping,
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          sub_no,
    const hp::QCollection<dim_1> &                              quadrature,
    const CellSimilarity::Similarity                            cell_similarity,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_data,
    const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &mapping_data,
    internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
      &output_data) const;

private:
  /**
   * 表示一个给定的面或子面编号无效的值。
   *
   */
  static const unsigned int invalid_face_number = numbers::invalid_unsigned_int;

  /**
   * 指向底层有限元对象的指针。
   * 这个对象包含一个指向混合离散化的每个贡献元素的指针和它的倍率。它是由构造函数创建的，之后是常量。
   *
   */
  std::vector<std::pair<std::unique_ptr<const FiniteElement<dim, spacedim>>,
                        unsigned int>>
    base_elements;

  /**
   * 一个索引表，将基础元素的广义支持点映射到FE系统的广义支持点的矢量。
   * 它成立的原因是
   * @code
   * auto n = generalized_support_points_index_table[i][j];
   * generalized_support_points[n] ==
   *         base_elements[i].generalized_support_points[j];
   * @endcode
   * 对于每个基元（以i为索引）和基元的每个g.s.点（以j为索引）。
   *
   */
  std::vector<std::vector<std::size_t>> generalized_support_points_index_table;

  /**
   * 这个函数是简单地从构造函数中挑出来的，因为它有几个。它设置了系统的索引表，以及
   * @p 限制和 @p prolongation 矩阵。
   *
   */
  void
  initialize(const std::vector<const FiniteElement<dim, spacedim> *> &fes,
             const std::vector<unsigned int> &multiplicities);

  /**
   * 由 @p initialize. 使用。
   *
   */
  void
  build_interface_constraints();

  /**
   * 一个计算hp_vertex_dof_identities()、hp_line_dof_identities()或hp_quad_dof_identities()的函数，这取决于模板参数的值。
   *
   */
  template <int structdim>
  std::vector<std::pair<unsigned int, unsigned int>>
  hp_object_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                           const unsigned int face_no = 0) const;

  /**
   * 通常情况下。独立于细胞的数据字段。
   * 然而，在这里，这个类本身并不存储数据，而只是指向每个基本元素的
   * @p InternalData 对象的指针。
   *
   */
  class InternalData : public FiniteElement<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * 构造函数。由 @p get_data 函数调用。设置 @p base_fe_datas
     * 向量的大小为 @p n_base_elements. 。
     *
     */
    InternalData(const unsigned int n_base_elements);

    /**
     * 销毁器。删除所有 @p InternalDatas ，其指针由 @p
     * base_fe_datas 向量存储。
     *
     */
    ~InternalData() override;

    /**
     * 对 @p base_noth基元的 @p InternalData 的指针给予写权限。
     *
     */
    void
    set_fe_data(
      const unsigned int base_no,
      std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>);

    /**
     * 给予对 @p InternalData 的 @p
     * 基数元素的指针的读访问权。
     *
     */
    typename FiniteElement<dim, spacedim>::InternalDataBase &
    get_fe_data(const unsigned int base_no) const;

    /**
     * 当调用 FiniteElement::fill_fe_values()
     * 和类似函数时，给读访问指向 <code>base_no</code>
     * 第1个基元的对象的指针，该对象将写入其输出。
     *
     */
    internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &
    get_fe_output_object(const unsigned int base_no) const;

  private:
    /**
     * 指向每个基元的 @p InternalData 对象的指针。它们被 @p
     * set_ 和 @p get_fe_data 函数所访问。
     * 这个向量的大小由InternalData构造函数设置为 @p
     * n_base_elements 。 它由 @p get_data 函数填充。
     * 请注意，由于基类的每个实例的数据必然是相同的，我们只需要有多少个基类元素就有多少个这样的对象，而不考虑它们的多重性。
     *
     */
    typename std::vector<
      std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>>
      base_fe_datas;

    /**
     * 一个对象的集合，当我们对它们调用
     * FiniteElement::fill_fe_values()
     * 和相关函数时，基元将把它们的输出写入其中。
     * 这个向量的大小由InternalData构造函数设置为 @p
     * n_base_elements 。
     *
     */
    mutable std::vector<
      internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>>
      base_fe_output_objects;
  };

  /**
   * 用于保护限制和嵌入矩阵的初始化的互斥器。
   *
   */
  mutable std::mutex mutex;

  friend class FE_Enriched<dim, spacedim>;
};

//------------------------variadic template constructor------------------------

#  ifndef DOXYGEN
namespace internal
{
  namespace FESystemImplementation
  {
    template <int dim, int spacedim>
    unsigned int
    count_nonzeros(
      const std::initializer_list<
        std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>>
        &fe_systems)
    {
      return std::count_if(
        fe_systems.begin(),
        fe_systems.end(),
        [](const std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
                           unsigned int> &fe_system) {
          return fe_system.second > 0;
        });
    }



    template <int dim, int spacedim>
    std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>
    promote_to_fe_pair(const FiniteElement<dim, spacedim> &fe)
    {
      return std::make_pair(std::move(fe.clone()), 1u);
    }



    template <int dim, int spacedim>
    auto
    promote_to_fe_pair(std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
                                 unsigned int> &&p)
      -> decltype(
        std::forward<std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
                               unsigned int>>(p))
    {
      return std::forward<
        std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>>(
        p);
    }
  } // namespace FESystemImplementation
} // namespace internal



#    if !defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1900
// We are just forwarding/delegating to the constructor taking a
// std::initializer_list. If we decide to remove the deprecated constructors, we
// might just use the variadic constructor with a suitable static_assert instead
// of the std::enable_if.
template <int dim, int spacedim>
template <class... FEPairs, typename>
FESystem<dim, spacedim>::FESystem(FEPairs &&... fe_pairs)
  : FESystem<dim, spacedim>(
      {internal::FESystemImplementation::promote_to_fe_pair<dim, spacedim>(
        std::forward<FEPairs>(fe_pairs))...})
{}



template <int dim, int spacedim>
FESystem<dim, spacedim>::FESystem(
  const std::initializer_list<
    std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>>
    &fe_systems)
  : FiniteElement<dim, spacedim>(
      FETools::Compositing::multiply_dof_numbers<dim, spacedim>(fe_systems),
      FETools::Compositing::compute_restriction_is_additive_flags<dim,
                                                                  spacedim>(
        fe_systems),
      FETools::Compositing::compute_nonzero_components<dim, spacedim>(
        fe_systems))
  , base_elements(internal::FESystemImplementation::count_nonzeros(fe_systems))
{
  std::vector<const FiniteElement<dim, spacedim> *> fes;
  std::vector<unsigned int>                         multiplicities;

  const auto extract =
    [&fes, &multiplicities](
      const std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
                      unsigned int> &fe_system) {
      fes.push_back(fe_system.first.get());
      multiplicities.push_back(fe_system.second);
    };

  for (const auto &p : fe_systems)
    extract(p);

  initialize(fes, multiplicities);
}
#    endif

#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
 /*----------------------------  fe_system.h  ---------------------------*/ 


