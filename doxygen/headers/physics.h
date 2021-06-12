//include/deal.II-translator/A-headers/physics_0.txt
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

/**
 *    @defgroup physics Physics
 *
 * @brief A module dedicated to the implementation of functions and
 * 与连续体物理学、物理场和材料有关的课程。
 *
 */

/**
 * 命名空间和实用程序的集合，以协助定义、构建和操作与物理场和材料有关的数据。
 *
 */
namespace Physics
{

  /**
   * 减少张量的顺序的符号，有效地将它们存储在某种一致的压缩存储模式中。一个例子是将
   * $3\times 3$
   * 等级2的对称张量的6个独立分量存储为一个有6个分量的向量，然后将等级4的对称
   * $3\times 3 \times 3\times 3$
   * 张量的36个独立元素（当应用于一个对称等级2张量时，会产生另一个对称等级2张量）表示为一个
   * $6 \times 6$ 矩阵。
   * 尽管这种表示张量的方法最常与四阶弹性张量的有效存储联系在一起，但随着它的泛化，它具有更广泛的适用性。这种表示方法在物理学、材料科学和有限元文献中也很常见。
   * 张量符号有几种变化，每一种结构都略有不同。各种形式的张量符号之间的主要区别是对压缩张量的各种元素所规定的权重。
   * 这个<a href="https://en.wikipedia.org/wiki/Voigt_notation">wikipedia
   * article</a>对这个主题有一些进一步的一般见解。
   * @ingroup physics   @author  Jean-Paul Pelteret, 2017
   *
   */
  namespace Notation
  {

  }

  /**
   * 一组操作，以协助张量从参考到空间配置的转换，反之亦然。
   * 这些类型的转换通常用于以第二种配置重新表达在一种配置中测量或计算的数量。
   * <h3>Notation</h3> 我们将对坐标 $\mathbf{X}, \mathbf{x}$ 、变换
   * $\varphi$ 、微分算子 $\nabla_{0}$ 和变形梯度 $\mathbf{F}$
   * 使用与命名空间 Physics::Elasticity. 相同的符号。
   * 作为符号的另一点，我们将遵循Holzapfel（2007）的做法，将前推变换表示为
   * $\chi\left(\bullet\right)$ ，后拉变换表示为
   * $\chi^{-1}\left(\bullet\right)$  。  我们还将使用注释
   * $\left(\bullet\right)^{\sharp}$ 来表示张量 $\left(\bullet\right)$
   * 是反变量张量， $\left(\bullet\right)^{\flat}$
   * 表示它是协变的。换句话说，这些索引实际上并不改变张量，它们只是表明某个张量的<i>kind</i>对象。
   *
   * @note
   * 对于这些变换，除非另有说明，我们将严格假定变换后的张量的所有指数来自一个坐标系；也就是说，它们不是多点张量（如弹性中的皮奥拉应力）。
   *
   * @ingroup physics
   * @author  Jean-Paul Pelteret, Andrew McBride, 2016
   *
   */
  namespace Transformations
  {
  }

  /**
   * 这个命名空间提供了一个符合（非线性）弹性中使用的标准符号的定义集合。
   * <h3>Notation</h3> 这个符号的参考文献包括：。
   * @code{.bib}
   *   @Book{Holzapfel2007a,
   *      title =     {Nonlinear solid mechanics. A Continuum Approach for Engineering},
   *      publisher = {John Wiley \& Sons Ltd.},
   *      year =      {2007},
   *      author =    {Holzapfel, G. A.},
   *      address =   {West Sussex, England},
   *      note =      {ISBN: 0-471-82304-X}
   *    }
   *    @Book{Wriggers2008a,
   *      title =     {Nonlinear finite element methods},
   *      publisher = {Springer Berlin Heidelberg},
   *      year =      {2008},
   *      author =    {Wriggers, P.},
   *      volume =    {4},
   *      address =   {Berlin, Germany},
   *      note =      {ISBN: 978-3-540-71000-4},
   *      doi =       {10.1007/978-3-540-71001-1}
   *    }
   * @endcode
   * 为方便起见，我们将预先定义一些常用的参考张量和操作。   考虑到参考（材料）配置中的位置向量 $\mathbf{X}$ ，点 $\mathbf{X}$ 通过非线性图@f[
   * \mathbf{x}
   * \dealcoloneq \boldsymbol{\varphi} \left( \mathbf{X} \right)
   *  = \mathbf{X} + \mathbf{u}(\mathbf{X}) \, ,
   * @f]被转换为当前（空间）配置中的点 $\mathbf{x}$ ，其中 $\mathbf{u}(\mathbf{X})$ 代表位移向量。   由此我们可以计算出变形梯度张量为@f[
   * \mathbf{F} \dealcoloneq \mathbf{I} + \nabla_{0}\mathbf{u} \, ,
   * @f]，其中微分算子 $\nabla_{0}$ 被定义为 $\frac{\partial}{\partial \mathbf{X}}$ ， $\mathbf{I}$ 是身份张量。     最后，两个普通张量算子由 $\cdot$ 和 $:$ 算子表示。它们分别代表对内部张量指数的单缩和双缩。   向量和二阶张量用粗体字突出，而四阶张量则用卡列字体表示。     我们可以把四阶张量看作是将二阶张量（矩阵）映射到自己身上的线性运算符，其方式与矩阵将向量映射到向量上一样。   为了给已实现的类成员和函数提供一些背景，考虑对具有特殊属性的张量进行以下基本操作。     如果我们把一般的二阶张量表示为 $\mathbf{A}$ ，那么一般的四阶单位张量 $\mathcal{I}$ 和 $\overline{\mathcal{I}}$ 由@f[
   * \mathbf{A} = \mathcal{I}:\mathbf{A} \qquad \text{and} \qquad \mathbf{A}^T
   * = \overline{\mathcal{I}}:\mathbf{A} \, ,
   * @f]定义，或者用表记法表示为@f[
   * I_{ijkl} = \delta_{ik}\delta_{jl} \qquad \text{and} \qquad \overline
   * I_{ijkl} = \delta_{il}\delta_{jk}
   * @f]，克朗克三角采用其共同定义。   请注意， $\mathcal{I} \neq \overline{\mathcal{I}}^T$  。     然后我们用@f[
   * \mathcal{S} \dealcoloneq \dfrac{1}{2}[\mathcal{I} +
   * \overline{\mathcal{I}}] \qquad \text{and} \qquad \mathcal{W} \dealcoloneq
   * \dfrac{1}{2}[\mathcal{I}
   *
   * - \overline{\mathcal{I}}] \, ,
   * @f]定义对称和偏斜对称的四阶单位张量，这样@f[
   * \mathcal{S}:\mathbf{A} = \dfrac{1}{2}[\mathbf{A} + \mathbf{A}^T] \qquad
   * \text{and} \qquad \mathcal{W}:\mathbf{A} = \dfrac{1}{2}[\mathbf{A}
   *
   * - \mathbf{A}^T] \, .
   * @f] identity_tensor()返回的四阶对称张量是 $\mathcal{S}$  。
   * @author  Jean-Paul Pelteret, Andrew McBride, 2016年
   * @ingroup physics
   *
   */
  namespace Elasticity
  {
  }

}


