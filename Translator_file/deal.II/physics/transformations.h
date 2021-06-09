//include/deal.II-translator/physics/transformations_0.txt
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

#ifndef dealii_transformations_h
#define dealii_transformations_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

DEAL_II_NAMESPACE_OPEN


namespace Physics
{
  namespace Transformations
  {
    /**
     * 以旋转角度和旋转轴定义的变换函数和张量。
     *
     */
    namespace Rotations
    {
      /**
       * @name 旋转矩阵
       *
       */
      //@{

      /**
       * 返回2维欧氏空间的旋转矩阵，即@f[
       * \mathbf{R} \dealcoloneq \left[ \begin{array}{cc}
       * cos(\theta) &
       *
       * -sin(\theta) \\
       * sin(\theta) & cos(\theta)
       * \end{array}\right]
       * @f]，其中 $\theta$ 是以弧度给出的旋转角度。特别是，这描述了一个向量相对于<a href="http://mathworld.wolfram.com/RotationMatrix.html">fixed
       * set of right-handed axes</a>的逆时针旋转。
       * @param[in] 角度
       * 以弧度为单位的旋转角度（关于Z轴）。
       *
       */
      template <typename Number>
      Tensor<2, 2, Number>
      rotation_matrix_2d(const Number &angle);


      /**
       * 返回3维欧几里得空间的旋转矩阵。用罗德里格斯的旋转公式最简洁地说明，这个函数返回相当于@f[
       * \mathbf{R} \dealcoloneq cos(\theta)\mathbf{I} + sin(\theta)\mathbf{W}
       *            + (1-cos(\theta))\mathbf{u}\otimes\mathbf{u}
       * @f]，其中 $\mathbf{u}$ 是轴向量（一个轴向向量）， $\theta$ 是以弧度给出的旋转角度， $\mathbf{I}$ 是身份张量， $\mathbf{W}$ 是 $\mathbf{u}$ 的斜对称张量。              @dealiiWriggersA{374,9.194}  这提出了罗德里格斯的旋转公式，但在这个函数中使用的实现在这个<a
       * href="https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle">wikipedia
       * link</a>中描述。特别是，这描述了一个向量<a
       * href="http://mathworld.wolfram.com/RotationMatrix.html">in a plane
       * with its normal</a>的逆时针旋转。由旋转的 @p axis
       * 定义。在<a
       * href="https://www.gamedev.net/resources/_/technical/math-and-physics/do-we-really-need-quaternions-r1199">this
       * link</a>中讨论了另一种实现方式，但与罗德里格斯的旋转公式不一致（符号上），因为它描述的是坐标系的旋转。
       * @param[in]  axis 一个定义旋转轴的单位向量  @param[in]
       * angle 旋转角度，单位为弧度
       *
       */
      template <typename Number>
      Tensor<2, 3, Number>
      rotation_matrix_3d(const Point<3, Number> &axis, const Number &angle);

      //@}

    } // namespace Rotations

    /**
     * 以一组禁忌基数定义的张量的变换。等级1和等级2的忌变张量 $\left(\bullet\right)^{\sharp} = \mathbf{T}$ （以及它的空间对应物 $\mathbf{t}$ ）通常满足关系@f[
     *  \int_{V_{0}} \nabla_{0} \cdot \mathbf{T} \; dV
     *    = \int_{\partial V_{0}} \mathbf{T} \cdot \mathbf{N} \; dA
     *    = \int_{\partial V_{t}} \mathbf{T} \cdot \mathbf{n} \; da
     *    = \int_{V_{t}} \nabla \cdot \mathbf{t} \; dv
     * @f]，其中 $V_{0}$ 和 $V_{t}$ 分别是参考和空间配置中的控制体积，其表面 $\partial
     * V_{0}$ 和 $\partial V_{t}$ 具有向外的法线 $\mathbf{N}$ 和
     * $\mathbf{n}$  。
     *
     */
    namespace Contravariant
    {
      /**
       * @name  向前推操作
       *
       */
      //@{

      /**
       * 返回对一个禁忌向量进行前推变换的结果，即 @f[
       * \chi\left(\bullet\right)^{\sharp}
       *  \dealcoloneq \mathbf{F} \cdot \left(\bullet\right)^{\sharp}
       * @f]  @param[in]  V 被操作的（参考）向量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi\left( \mathbf{V} \right)$
       *
       */
      template <int dim, typename Number>
      Tensor<1, dim, Number>
      push_forward(const Tensor<1, dim, Number> &V,
                   const Tensor<2, dim, Number> &F);

      /**
       * 返回对第2级禁忌张量的前推变换结果，即 @f[
       * \chi\left(\bullet\right)^{\sharp}
       *  \dealcoloneq \mathbf{F} \cdot \left(\bullet\right)^{\sharp} \cdot
       * \mathbf{F}^{T}
       * @f]  @param[in]  T 要操作的（参考）第2级张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi\left( \mathbf{T} \right)$  。
       *
       */
      template <int dim, typename Number>
      Tensor<2, dim, Number>
      push_forward(const Tensor<2, dim, Number> &T,
                   const Tensor<2, dim, Number> &F);

      /**
       * 返回对等级2的禁忌对称张量进行前推变换的结果，即 @f[
       * \chi\left(\bullet\right)^{\sharp}
       *  \dealcoloneq \mathbf{F} \cdot \left(\bullet\right)^{\sharp} \cdot
       * \mathbf{F}^{T}
       * @f]  @param[in]  T 要操作的（参考）等级2的对称张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi\left( \mathbf{T} \right)$  。
       *
       */
      template <int dim, typename Number>
      SymmetricTensor<2, dim, Number>
      push_forward(const SymmetricTensor<2, dim, Number> &T,
                   const Tensor<2, dim, Number> &         F);

      /**
       * 返回对一个等级4的禁忌张量进行前推变换的结果，即（用索引符号表示）。      @f[
       * \left[ \chi\left(\bullet\right)^{\sharp} \right]_{ijkl}
       *  \dealcoloneq F_{iI} F_{jJ}
       *  \left(\bullet\right)^{\sharp}_{IJKL} F_{kK} F_{lL}
       * @f]  @param[in]  H 要操作的（参考）秩-4张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi\left( \mathbf{H} \right)$
       *
       */
      template <int dim, typename Number>
      Tensor<4, dim, Number>
      push_forward(const Tensor<4, dim, Number> &H,
                   const Tensor<2, dim, Number> &F);

      /**
       * 返回对等级4的禁忌对称张量进行推前变换的结果，即（用索引符号表示）。      @f[
       * \left[ \chi\left(\bullet\right)^{\sharp} \right]_{ijkl}
       *  \dealcoloneq F_{iI} F_{jJ}
       *  \left(\bullet\right)^{\sharp}_{IJKL} F_{kK} F_{lL}
       * @f]  @param[in]  H 要操作的（参考）等级-4对称张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi\left( \mathbf{H} \right)$
       *
       */
      template <int dim, typename Number>
      SymmetricTensor<4, dim, Number>
      push_forward(const SymmetricTensor<4, dim, Number> &H,
                   const Tensor<2, dim, Number> &         F);

      //@}

      /**
       * @name  拉回操作
       *
       */
      //@{

      /**
       * 返回对一个禁忌向量进行回拉变换的结果，即 @f[
       * \chi^{-1}\left(\bullet\right)^{\sharp}
       *  \dealcoloneq \mathbf{F}^{-1} \cdot \left(\bullet\right)^{\sharp}
       * @f]  @param[in]  v 被操作的（空间）向量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi^{-1}\left( \mathbf{v} \right)$
       *
       */
      template <int dim, typename Number>
      Tensor<1, dim, Number>
      pull_back(const Tensor<1, dim, Number> &v,
                const Tensor<2, dim, Number> &F);

      /**
       * 返回2级禁忌张量的回拉变换结果，即 @f[
       * \chi^{-1}\left(\bullet\right)^{\sharp}
       *  \dealcoloneq \mathbf{F}^{-1} \cdot \left(\bullet\right)^{\sharp}
       *  \cdot \mathbf{F}^{-T}
       * @f]  @param[in]  t 要操作的（空间）张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi^{-1}\left( \mathbf{t} \right)$
       * 。
       *
       */
      template <int dim, typename Number>
      Tensor<2, dim, Number>
      pull_back(const Tensor<2, dim, Number> &t,
                const Tensor<2, dim, Number> &F);

      /**
       * 返回2级禁忌对称张量的回拉变换结果，即 @f[
       * \chi^{-1}\left(\bullet\right)^{\sharp}
       *  \dealcoloneq \mathbf{F}^{-1} \cdot \left(\bullet\right)^{\sharp}
       *  \cdot \mathbf{F}^{-T}
       * @f]  @param[in]  t 要操作的（空间）对称张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi^{-1}\left( \mathbf{t} \right)$
       * 。
       *
       */
      template <int dim, typename Number>
      SymmetricTensor<2, dim, Number>
      pull_back(const SymmetricTensor<2, dim, Number> &t,
                const Tensor<2, dim, Number> &         F);

      /**
       * 返回对一个等级为4的禁忌张量进行回拉变换的结果，即（用索引符号表示）。      @f[
       * \left[ \chi^{-1}\left(\bullet\right)^{\sharp} \right]_{IJKL}
       *  \dealcoloneq F^{-1}_{Ii} F^{-1}_{Jj}
       * \left(\bullet\right)^{\sharp}_{ijkl} F^{-1}_{Kk} F^{-1}_{Ll}
       * @f]  @param[in]  h 要操作的（空间）张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi^{-1}\left( \mathbf{h} \right)$
       *
       */
      template <int dim, typename Number>
      Tensor<4, dim, Number>
      pull_back(const Tensor<4, dim, Number> &h,
                const Tensor<2, dim, Number> &F);

      /**
       * 返回对等级4的禁忌对称张量进行回拉变换的结果，即（用索引符号表示）。      @f[
       * \left[ \chi^{-1}\left(\bullet\right)^{\sharp} \right]_{IJKL}
       *  \dealcoloneq F^{-1}_{Ii} F^{-1}_{Jj}
       *  \left(\bullet\right)^{\sharp}_{ijkl} F^{-1}_{Kk} F^{-1}_{Ll}
       * @f]  @param[in]  h 要操作的（空间）对称张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi^{-1}\left( \mathbf{h} \right)$
       *
       */
      template <int dim, typename Number>
      SymmetricTensor<4, dim, Number>
      pull_back(const SymmetricTensor<4, dim, Number> &h,
                const Tensor<2, dim, Number> &         F);

      //@}
    } // namespace Contravariant

    /**
     * 以一组协变基向量定义的张量的变换。等级1和等级2协变张量 $\left(\bullet\right)^{\flat} = \mathbf{T}$ （及其空间对应物 $\mathbf{t}$ ）通常满足关系@f[
     *  \int_{\partial V_{0}} \left[ \nabla_{0} \times \mathbf{T} \right]
     * \cdot \mathbf{N} \; dA = \oint_{\partial A_{0}} \mathbf{T} \cdot
     * \mathbf{L} \; dL = \oint_{\partial A_{t}} \mathbf{t} \cdot \mathbf{l} \;
     * dl = \int_{\partial V_{t}} \left[ \nabla \times \mathbf{t} \right] \cdot
     * \mathbf{n} \; da
     * @f]，其中控制面 $\partial V_{0}$ 和 $\partial V_{t}$  ] 的法线 $\mathbf{N}$ 和 $\mathbf{n}$ 被曲线 $\partial A_{0}$ 和 $\partial A_{t}$ 所约束，这些曲线分别与线长 $\mathbf{L}$ 和 $\mathbf{l}$ 有关。
     *
     */
    namespace Covariant
    {
      /**
       * @name  推进操作
       *
       */
      //@{

      /**
       * 返回共变向量上的前推变换结果，即 @f[
       * \chi\left(\bullet\right)^{\flat}
       *  \dealcoloneq \mathbf{F}^{-T} \cdot \left(\bullet\right)^{\flat}
       * @f]  @param[in]  V 被操作的（参考）向量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi\left( \mathbf{V} \right)$
       *
       */
      template <int dim, typename Number>
      Tensor<1, dim, Number>
      push_forward(const Tensor<1, dim, Number> &V,
                   const Tensor<2, dim, Number> &F);

      /**
       * 返回对等级2共变张量的前推变换结果，即 @f[
       * \chi\left(\bullet\right)^{\flat}
       *  \dealcoloneq \mathbf{F}^{-T} \cdot \left(\bullet\right)^{\flat}
       *  \cdot \mathbf{F}^{-1}
       * @f]  @param[in]  T 要操作的（参考）等级2张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi\left( \mathbf{T} \right)$  。
       *
       */
      template <int dim, typename Number>
      Tensor<2, dim, Number>
      push_forward(const Tensor<2, dim, Number> &T,
                   const Tensor<2, dim, Number> &F);

      /**
       * 返回对等级2共变对称张量的前推变换结果，即 @f[
       * \chi\left(\bullet\right)^{\flat}
       *  \dealcoloneq \mathbf{F}^{-T} \cdot \left(\bullet\right)^{\flat}
       *  \cdot \mathbf{F}^{-1}
       * @f]  @param[in]  T 要操作的（参考）等级2对称张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi\left( \mathbf{T} \right)$  。
       *
       */
      template <int dim, typename Number>
      SymmetricTensor<2, dim, Number>
      push_forward(const SymmetricTensor<2, dim, Number> &T,
                   const Tensor<2, dim, Number> &         F);

      /**
       * 返回对等级4共变张量的前推变换的结果，即（用索引符号表示）。      @f[
       * \left[ \chi\left(\bullet\right)^{\flat} \right]_{ijkl}
       *  \dealcoloneq F^{-T}_{iI} F^{-T}_{jJ}
       *  \left(\bullet\right)^{\flat}_{IJKL} F^{-T}_{kK} F^{-T}_{lL}
       * @f]  @param[in]  H 要操作的（参考）等级-4张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi\left( \mathbf{H} \right)$
       *
       */
      template <int dim, typename Number>
      Tensor<4, dim, Number>
      push_forward(const Tensor<4, dim, Number> &H,
                   const Tensor<2, dim, Number> &F);

      /**
       * 返回对等级4共变对称张量的前推变换的结果，即（用索引符号表示）。      @f[
       * \left[ \chi\left(\bullet\right)^{\flat} \right]_{ijkl}
       *  \dealcoloneq F^{-T}_{iI} F^{-T}_{jJ}
       *  \left(\bullet\right)^{\flat}_{IJKL} F^{-T}_{kK} F^{-T}_{lL}
       * @f]  @param[in]  H 要操作的（参考）等级-4对称张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi\left( \mathbf{H} \right)$ .
       *
       */
      template <int dim, typename Number>
      SymmetricTensor<4, dim, Number>
      push_forward(const SymmetricTensor<4, dim, Number> &H,
                   const Tensor<2, dim, Number> &         F);

      //@}

      /**
       * @name  拉回操作
       *
       */
      //@{

      /**
       * 返回共变向量上的回拉变换结果，即 @f[
       * \chi^{-1}\left(\bullet\right)^{\flat}
       *  \dealcoloneq \mathbf{F}^{T} \cdot \left(\bullet\right)^{\flat}
       * @f]  @param[in]  v 被操作的（空间）向量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi^{-1}\left( \mathbf{v} \right)$
       *
       */
      template <int dim, typename Number>
      Tensor<1, dim, Number>
      pull_back(const Tensor<1, dim, Number> &v,
                const Tensor<2, dim, Number> &F);

      /**
       * 返回2级协变张量的回拉变换结果，即 @f[
       * \chi^{-1}\left(\bullet\right)^{\flat}
       *  \dealcoloneq \mathbf{F}^{T} \cdot \left(\bullet\right)^{\flat} \cdot
       * \mathbf{F}
       * @f]  @param[in]  t 要操作的（空间）张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi^{-1}\left( \mathbf{t} \right)$
       *
       */
      template <int dim, typename Number>
      Tensor<2, dim, Number>
      pull_back(const Tensor<2, dim, Number> &t,
                const Tensor<2, dim, Number> &F);

      /**
       * 返回对等级2共变对称张量的回拉变换结果，即 @f[
       * \chi^{-1}\left(\bullet\right)^{\flat}
       *  \dealcoloneq \mathbf{F}^{T} \cdot \left(\bullet\right)^{\flat}
       *  \cdot \mathbf{F}
       * @f]  @param[in]  t 要操作的（空间）对称张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi^{-1}\left( \mathbf{t} \right)$
       * 。
       *
       */
      template <int dim, typename Number>
      SymmetricTensor<2, dim, Number>
      pull_back(const SymmetricTensor<2, dim, Number> &t,
                const Tensor<2, dim, Number> &         F);

      /**
       * 返回对一个等级4的禁忌张量进行回拉变换的结果，即（用索引符号表示）。      @f[
       * \left[ \chi^{-1}\left(\bullet\right)^{\flat} \right]_{IJKL}
       * \dealcoloneq F^{T}_{Ii} F^{T}_{Jj}
       * \left(\bullet\right)^{\flat}_{ijkl} F^{T}_{Kk} F^{T}_{Ll}
       * @f]  @param[in]  h 要操作的（空间）张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi^{-1}\left( \mathbf{h} \right)$
       *
       */
      template <int dim, typename Number>
      Tensor<4, dim, Number>
      pull_back(const Tensor<4, dim, Number> &h,
                const Tensor<2, dim, Number> &F);

      /**
       * 返回对等级4的禁忌对称张量进行回拉变换的结果，即（用索引符号表示）。      @f[
       * \left[ \chi^{-1}\left(\bullet\right)^{\flat} \right]_{IJKL}
       * \dealcoloneq F^{T}_{Ii} F^{T}_{Jj}
       * \left(\bullet\right)^{\flat}_{ijkl} F^{T}_{Kk} F^{T}_{Ll}
       * @f]  @param[in]  h 要操作的（空间）对称张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\chi^{-1}\left( \mathbf{h} \right)$
       * 。
       *
       */
      template <int dim, typename Number>
      SymmetricTensor<4, dim, Number>
      pull_back(const SymmetricTensor<4, dim, Number> &h,
                const Tensor<2, dim, Number> &         F);

      //@}
    } // namespace Covariant

    /**
     * 张量的变换是以一组禁忌基向量来定义的，并以与映射相关的体积变化的逆数来缩放。
     *
     */
    namespace Piola
    {
      /**
       * @name  前推操作
       *
       */
      //@{

      /**
       * 返回对一个禁忌向量进行前推变换的结果，即 @f[
       * \textrm{det} \mathbf{F}^{-1} \; \chi\left(\bullet\right)^{\sharp}
       * \dealcoloneq \frac{1}{\textrm{det} \mathbf{F}} \; \mathbf{F} \cdot
       * \left(\bullet\right)^{\sharp}
       * @f]  @param[in]  V 被操作的（参考）向量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\frac{1}{\textrm{det} \mathbf{F}} \;
       * \chi\left( \mathbf{V} \right)$
       *
       */
      template <int dim, typename Number>
      Tensor<1, dim, Number>
      push_forward(const Tensor<1, dim, Number> &V,
                   const Tensor<2, dim, Number> &F);

      /**
       * 返回对第2级禁忌张量的前推变换结果，即 @f[
       * \textrm{det} \mathbf{F}^{-1} \; \chi\left(\bullet\right)^{\sharp}
       *  \dealcoloneq \frac{1}{\textrm{det} \mathbf{F}} \; \mathbf{F} \cdot
       * \left(\bullet\right)^{\sharp} \cdot \mathbf{F}^{T}
       * @f]  @param[in]  T 被操作的（参考）第2级张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\frac{1}{\textrm{det} \mathbf{F}} \;
       * \chi\left( \mathbf{T} \right)$
       *
       */
      template <int dim, typename Number>
      Tensor<2, dim, Number>
      push_forward(const Tensor<2, dim, Number> &T,
                   const Tensor<2, dim, Number> &F);

      /**
       * 返回对等级2的禁忌对称张量进行前推变换的结果，即 @f[
       * \textrm{det} \mathbf{F}^{-1} \; \chi\left(\bullet\right)^{\sharp}
       *  \dealcoloneq \frac{1}{\textrm{det} \mathbf{F}} \; \mathbf{F} \cdot
       * \left(\bullet\right)^{\sharp} \cdot \mathbf{F}^{T}
       * @f]  @param[in]  T 要操作的（参考）等级2的对称张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\frac{1}{\textrm{det} \mathbf{F}} \;
       * \chi\left( \mathbf{T} \right)$  。
       *
       */
      template <int dim, typename Number>
      SymmetricTensor<2, dim, Number>
      push_forward(const SymmetricTensor<2, dim, Number> &T,
                   const Tensor<2, dim, Number> &         F);

      /**
       * 返回对一个等级4的禁忌张量进行前推变换的结果，即（用索引符号表示）。      @f[
       * \textrm{det} \mathbf{F}^{-1} \; \left[
       * \chi\left(\bullet\right)^{\sharp} \right]_{ijkl}
       *  \dealcoloneq \frac{1}{\textrm{det} \mathbf{F}} \; F_{iI} F_{jJ}
       * \left(\bullet\right)^{\sharp}_{IJKL} F_{kK} F_{lL}
       * @f]  @param[in]  H 被操作的（参考）等级-4张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\frac{1}{\textrm{det} \mathbf{F}} \;
       * \chi\left( \mathbf{H} \right)$
       *
       */
      template <int dim, typename Number>
      Tensor<4, dim, Number>
      push_forward(const Tensor<4, dim, Number> &H,
                   const Tensor<2, dim, Number> &F);

      /**
       * 返回对等级4的禁忌对称张量进行推前变换的结果，即（用索引符号表示）。      @f[
       * \textrm{det} \mathbf{F}^{-1} \; \left[
       * \chi\left(\bullet\right)^{\sharp} \right]_{ijkl}
       *  \dealcoloneq \frac{1}{\textrm{det} \mathbf{F}} \; F_{iI} F_{jJ}
       * \left(\bullet\right)^{\sharp}_{IJKL} F_{kK} F_{lL}
       * @f]  @param[in]  H 要操作的（参考）等级4对称张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\frac{1}{\textrm{det} \mathbf{F}} \;
       * \chi\left( \mathbf{H} \right)$
       *
       */
      template <int dim, typename Number>
      SymmetricTensor<4, dim, Number>
      push_forward(const SymmetricTensor<4, dim, Number> &H,
                   const Tensor<2, dim, Number> &         F);

      //@}

      /**
       * @name  拉回操作
       *
       */
      //@{

      /**
       * 返回对一个禁忌向量进行回拉变换的结果，即 @f[
       * \textrm{det} \mathbf{F} \; \chi^{-1}\left(\bullet\right)^{\sharp}
       *  \dealcoloneq \textrm{det} \mathbf{F} \; \mathbf{F}^{-1} \cdot
       * \left(\bullet\right)^{\sharp}
       * @f]  @param[in]  v 要操作的（空间）向量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\textrm{det} \mathbf{F} \;
       * \chi^{-1}\left( \mathbf{v} \right)$
       *
       */
      template <int dim, typename Number>
      Tensor<1, dim, Number>
      pull_back(const Tensor<1, dim, Number> &v,
                const Tensor<2, dim, Number> &F);

      /**
       * 返回2级禁忌张量的回拉变换结果，即 @f[
       * \textrm{det} \mathbf{F} \; \chi^{-1}\left(\bullet\right)^{\sharp}
       *  \dealcoloneq \textrm{det} \mathbf{F} \; \mathbf{F}^{-1} \cdot
       * \left(\bullet\right)^{\sharp} \cdot \mathbf{F}^{-T}
       * @f]  @param[in]  t 要操作的（空间）张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\textrm{det} \mathbf{F} \;
       * \chi^{-1}\left( \mathbf{t} \right)$ 。
       *
       */
      template <int dim, typename Number>
      Tensor<2, dim, Number>
      pull_back(const Tensor<2, dim, Number> &t,
                const Tensor<2, dim, Number> &F);

      /**
       * 返回2级禁忌对称张量的回拉变换结果，即 @f[
       * \textrm{det} \mathbf{F} \; \chi^{-1}\left(\bullet\right)^{\sharp}
       *  \dealcoloneq \textrm{det} \mathbf{F} \; \mathbf{F}^{-1} \cdot
       * \left(\bullet\right)^{\sharp} \cdot \mathbf{F}^{-T}
       * @f]  @param[in]  t 要操作的（空间）对称张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\textrm{det} \mathbf{F} \;
       * \chi^{-1}\left( \mathbf{t} \right)$  。
       *
       */
      template <int dim, typename Number>
      SymmetricTensor<2, dim, Number>
      pull_back(const SymmetricTensor<2, dim, Number> &t,
                const Tensor<2, dim, Number> &         F);

      /**
       * 返回对一个等级4的禁忌张量进行回拉变换的结果，即（用索引符号表示）。      @f[
       * \textrm{det} \mathbf{F} \; \left[
       * \chi^{-1}\left(\bullet\right)^{\sharp} \right]_{IJKL}
       *  \dealcoloneq \textrm{det} \mathbf{F} \; F^{-1}_{Ii} F^{-1}_{Jj}
       * \left(\bullet\right)^{\sharp}_{ijkl} F^{-1}_{Kk} F^{-1}_{Ll}
       * @f]  @param[in]  h 要操作的（空间）张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\textrm{det} \mathbf{F} \;
       * \chi^{-1}\left( \mathbf{h} \right)$  。
       *
       */
      template <int dim, typename Number>
      Tensor<4, dim, Number>
      pull_back(const Tensor<4, dim, Number> &h,
                const Tensor<2, dim, Number> &F);

      /**
       * 返回对等级4的禁忌对称张量进行回拉变换的结果，即（用索引符号表示）。      @f[
       * \textrm{det} \mathbf{F} \; \left[
       * \chi^{-1}\left(\bullet\right)^{\sharp} \right]_{IJKL}
       *  \dealcoloneq \textrm{det} \mathbf{F} \; F^{-1}_{Ii} F^{-1}_{Jj}
       * \left(\bullet\right)^{\sharp}_{ijkl} F^{-1}_{Kk} F^{-1}_{Ll}
       * @f]  @param[in]  h 要操作的（空间）对称张量  @param[in]  F 变形梯度张量  $\mathbf{F} \left(
       * \mathbf{X} \right)$   @return   $\textrm{det} \mathbf{F} \;
       * \chi^{-1}\left( \mathbf{h} \right)$  。
       *
       */
      template <int dim, typename Number>
      SymmetricTensor<4, dim, Number>
      pull_back(const SymmetricTensor<4, dim, Number> &h,
                const Tensor<2, dim, Number> &         F);

      //@}
    } // namespace Piola

    /**
     * @name  特殊操作
     *
     */
    //@{

    /**
     * 返回在非线性变换图 $\mathbf{x} = \boldsymbol{\varphi} \left( \mathbf{X} \right)$ 下，应用南森公式对材料表面积元素 $d\mathbf{A}$ 进行变换的结果。        返回的结果是由参考和空间表面元素之间的面积比缩放的空间法线，即 @f[
     * \mathbf{n} \frac{da}{dA}
     * \dealcoloneq \textrm{det} \mathbf{F} \, \mathbf{F}^{-T} \cdot \mathbf{N}
     * = \textrm{cof} \mathbf{F} \cdot \mathbf{N} \, .
     * @f]  @param[in]  N 参考法线单位向量  $\mathbf{N}$   @param[in]  F 变形梯度张量  $\mathbf{F} \left(
     * \mathbf{X} \right)$   @return  缩放的空间法向量  $\mathbf{n}
     * \frac{da}{dA}$   @dealiiHolzapfelA{75,2.55}   @dealiiWriggersA{23,3.11}
     *
     */
    template <int dim, typename Number>
    Tensor<1, dim, Number>
    nansons_formula(const Tensor<1, dim, Number> &N,
                    const Tensor<2, dim, Number> &F);

    //@}

    /**
     * @name  基准转换
     *
     */
    //@{

    /**
     * 返回一个改变了基数的向量，即 @f[
     * \mathbf{V}^{\prime} \dealcoloneq \mathbf{B} \cdot \mathbf{V}
     * @f]  @param[in]  V 要转换的向量  $\mathbf{V}$   @param[in]  B 转换矩阵  $\mathbf{B}$   @return   $\mathbf{V}^{\prime}$
     *
     */
    template <int dim, typename Number>
    Tensor<1, dim, Number>
    basis_transformation(const Tensor<1, dim, Number> &V,
                         const Tensor<2, dim, Number> &B);

    /**
     * 返回一个改变了基础的第2级张量，即 @f[
     * \mathbf{T}^{\prime} \dealcoloneq \mathbf{B} \cdot \mathbf{T} \cdot
     * \mathbf{B}^{T}
     * @f]  @param[in]  T 要转换的张量  $\mathbf{T}$   @param[in]  B 转换矩阵  $\mathbf{B}$   @return   $\mathbf{T}^{\prime}$  。
     *
     */
    template <int dim, typename Number>
    Tensor<2, dim, Number>
    basis_transformation(const Tensor<2, dim, Number> &T,
                         const Tensor<2, dim, Number> &B);

    /**
     * 返回一个改变了基础的对称秩-2张量，即 @f[
     * \mathbf{T}^{\prime} \dealcoloneq \mathbf{B} \cdot \mathbf{T} \cdot
     * \mathbf{B}^{T}
     * @f]  @param[in]  T 要转换的张量  $\mathbf{T}$   @param[in]  B 转换矩阵  $\mathbf{B}$   @return   $\mathbf{T}^{\prime}$  。
     *
     */
    template <int dim, typename Number>
    SymmetricTensor<2, dim, Number>
    basis_transformation(const SymmetricTensor<2, dim, Number> &T,
                         const Tensor<2, dim, Number> &         B);

    /**
     * 返回一个改变了基础的等级4张量，即（用索引符号表示）。    @f[
     * H_{ijkl}^{\prime} \dealcoloneq B_{iI} B_{jJ} H_{IJKL} B_{kK} B_{lL}
     * @f]  @param[in]  H 要转换的张量  $\mathbf{T}$   @param[in]  B 转换矩阵  $\mathbf{B}$   @return   $\mathbf{H}^{\prime}$ .
     *
     */
    template <int dim, typename Number>
    Tensor<4, dim, Number>
    basis_transformation(const Tensor<4, dim, Number> &H,
                         const Tensor<2, dim, Number> &B);

    /**
     * 返回一个改变了基础的对称秩-4张量，即（用索引符号表示）。    @f[
     * H_{ijkl}^{\prime} \dealcoloneq B_{iI} B_{jJ} H_{IJKL} B_{kK} B_{lL}
     * @f]  @param[in]  H 要转换的张量  $\mathbf{T}$   @param[in]  B 转换矩阵  $\mathbf{B}$   @return   $\mathbf{H}^{\prime}$  。
     *
     */
    template <int dim, typename Number>
    SymmetricTensor<4, dim, Number>
    basis_transformation(const SymmetricTensor<4, dim, Number> &H,
                         const Tensor<2, dim, Number> &         B);

    //@}

  } // namespace Transformations
} // namespace Physics



#ifndef DOXYGEN



template <typename Number>
Tensor<2, 2, Number>
Physics::Transformations::Rotations::rotation_matrix_2d(const Number &angle)
{
  const Number rotation[2][2] = {{std::cos(angle), -std::sin(angle)},
                                 {std::sin(angle), std::cos(angle)}};
  return Tensor<2, 2>(rotation);
}



template <typename Number>
Tensor<2, 3, Number>
Physics::Transformations::Rotations::rotation_matrix_3d(
  const Point<3, Number> &axis,
  const Number &          angle)
{
  Assert(std::abs(axis.norm() - 1.0) < 1e-9,
         ExcMessage("The supplied axial vector is not a unit vector."));
  const Number c              = std::cos(angle);
  const Number s              = std::sin(angle);
  const Number t              = 1. - c;
  const Number rotation[3][3] = {{t * axis[0] * axis[0] + c,
                                  t * axis[0] * axis[1] - s * axis[2],
                                  t * axis[0] * axis[2] + s * axis[1]},
                                 {t * axis[0] * axis[1] + s * axis[2],
                                  t * axis[1] * axis[1] + c,
                                  t * axis[1] * axis[2] - s * axis[0]},
                                 {t * axis[0] * axis[2] - s * axis[1],
                                  t * axis[1] * axis[2] + s * axis[0],
                                  t * axis[2] * axis[2] + c}};
  return Tensor<2, 3, Number>(rotation);
}



template <int dim, typename Number>
inline Tensor<1, dim, Number>
Physics::Transformations::Contravariant::push_forward(
  const Tensor<1, dim, Number> &V,
  const Tensor<2, dim, Number> &F)
{
  return Physics::Transformations::basis_transformation(V, F);
}



template <int dim, typename Number>
inline Tensor<2, dim, Number>
Physics::Transformations::Contravariant::push_forward(
  const Tensor<2, dim, Number> &T,
  const Tensor<2, dim, Number> &F)
{
  return Physics::Transformations::basis_transformation(T, F);
}



template <int dim, typename Number>
inline SymmetricTensor<2, dim, Number>
Physics::Transformations::Contravariant::push_forward(
  const SymmetricTensor<2, dim, Number> &T,
  const Tensor<2, dim, Number> &         F)
{
  return Physics::Transformations::basis_transformation(T, F);
}



template <int dim, typename Number>
inline Tensor<4, dim, Number>
Physics::Transformations::Contravariant::push_forward(
  const Tensor<4, dim, Number> &H,
  const Tensor<2, dim, Number> &F)
{
  return Physics::Transformations::basis_transformation(H, F);
}



template <int dim, typename Number>
inline SymmetricTensor<4, dim, Number>
Physics::Transformations::Contravariant::push_forward(
  const SymmetricTensor<4, dim, Number> &H,
  const Tensor<2, dim, Number> &         F)
{
  return Physics::Transformations::basis_transformation(H, F);
}



template <int dim, typename Number>
inline Tensor<1, dim, Number>
Physics::Transformations::Contravariant::pull_back(
  const Tensor<1, dim, Number> &v,
  const Tensor<2, dim, Number> &F)
{
  return Physics::Transformations::basis_transformation(v, invert(F));
}



template <int dim, typename Number>
inline Tensor<2, dim, Number>
Physics::Transformations::Contravariant::pull_back(
  const Tensor<2, dim, Number> &t,
  const Tensor<2, dim, Number> &F)
{
  return Physics::Transformations::basis_transformation(t, invert(F));
}



template <int dim, typename Number>
inline SymmetricTensor<2, dim, Number>
Physics::Transformations::Contravariant::pull_back(
  const SymmetricTensor<2, dim, Number> &t,
  const Tensor<2, dim, Number> &         F)
{
  return Physics::Transformations::basis_transformation(t, invert(F));
}



template <int dim, typename Number>
inline Tensor<4, dim, Number>
Physics::Transformations::Contravariant::pull_back(
  const Tensor<4, dim, Number> &h,
  const Tensor<2, dim, Number> &F)
{
  return Physics::Transformations::basis_transformation(h, invert(F));
}



template <int dim, typename Number>
inline SymmetricTensor<4, dim, Number>
Physics::Transformations::Contravariant::pull_back(
  const SymmetricTensor<4, dim, Number> &h,
  const Tensor<2, dim, Number> &         F)
{
  return Physics::Transformations::basis_transformation(h, invert(F));
}



template <int dim, typename Number>
inline Tensor<1, dim, Number>
Physics::Transformations::Covariant::push_forward(
  const Tensor<1, dim, Number> &V,
  const Tensor<2, dim, Number> &F)
{
  return Physics::Transformations::basis_transformation(V,
                                                        transpose(invert(F)));
}



template <int dim, typename Number>
inline Tensor<2, dim, Number>
Physics::Transformations::Covariant::push_forward(
  const Tensor<2, dim, Number> &T,
  const Tensor<2, dim, Number> &F)
{
  return Physics::Transformations::basis_transformation(T,
                                                        transpose(invert(F)));
}



template <int dim, typename Number>
inline SymmetricTensor<2, dim, Number>
Physics::Transformations::Covariant::push_forward(
  const SymmetricTensor<2, dim, Number> &T,
  const Tensor<2, dim, Number> &         F)
{
  return Physics::Transformations::basis_transformation(T,
                                                        transpose(invert(F)));
}



template <int dim, typename Number>
inline Tensor<4, dim, Number>
Physics::Transformations::Covariant::push_forward(
  const Tensor<4, dim, Number> &H,
  const Tensor<2, dim, Number> &F)
{
  return Physics::Transformations::basis_transformation(H,
                                                        transpose(invert(F)));
}



template <int dim, typename Number>
inline SymmetricTensor<4, dim, Number>
Physics::Transformations::Covariant::push_forward(
  const SymmetricTensor<4, dim, Number> &H,
  const Tensor<2, dim, Number> &         F)
{
  return Physics::Transformations::basis_transformation(H,
                                                        transpose(invert(F)));
}



template <int dim, typename Number>
inline Tensor<1, dim, Number>
Physics::Transformations::Covariant::pull_back(const Tensor<1, dim, Number> &v,
                                               const Tensor<2, dim, Number> &F)
{
  return Physics::Transformations::basis_transformation(v, transpose(F));
}



template <int dim, typename Number>
inline Tensor<2, dim, Number>
Physics::Transformations::Covariant::pull_back(const Tensor<2, dim, Number> &t,
                                               const Tensor<2, dim, Number> &F)
{
  return Physics::Transformations::basis_transformation(t, transpose(F));
}



template <int dim, typename Number>
inline SymmetricTensor<2, dim, Number>
Physics::Transformations::Covariant::pull_back(
  const SymmetricTensor<2, dim, Number> &t,
  const Tensor<2, dim, Number> &         F)
{
  return Physics::Transformations::basis_transformation(t, transpose(F));
}



template <int dim, typename Number>
inline Tensor<4, dim, Number>
Physics::Transformations::Covariant::pull_back(const Tensor<4, dim, Number> &h,
                                               const Tensor<2, dim, Number> &F)
{
  return Physics::Transformations::basis_transformation(h, transpose(F));
}



template <int dim, typename Number>
inline SymmetricTensor<4, dim, Number>
Physics::Transformations::Covariant::pull_back(
  const SymmetricTensor<4, dim, Number> &h,
  const Tensor<2, dim, Number> &         F)
{
  return Physics::Transformations::basis_transformation(h, transpose(F));
}



template <int dim, typename Number>
inline Tensor<1, dim, Number>
Physics::Transformations::Piola::push_forward(const Tensor<1, dim, Number> &V,
                                              const Tensor<2, dim, Number> &F)
{
  return Number(1.0 / determinant(F)) * Contravariant::push_forward(V, F);
}



template <int dim, typename Number>
inline Tensor<2, dim, Number>
Physics::Transformations::Piola::push_forward(const Tensor<2, dim, Number> &T,
                                              const Tensor<2, dim, Number> &F)
{
  return Number(1.0 / determinant(F)) * Contravariant::push_forward(T, F);
}



template <int dim, typename Number>
inline SymmetricTensor<2, dim, Number>
Physics::Transformations::Piola::push_forward(
  const SymmetricTensor<2, dim, Number> &T,
  const Tensor<2, dim, Number> &         F)
{
  return Number(1.0 / determinant(F)) * Contravariant::push_forward(T, F);
}



template <int dim, typename Number>
inline Tensor<4, dim, Number>
Physics::Transformations::Piola::push_forward(const Tensor<4, dim, Number> &H,
                                              const Tensor<2, dim, Number> &F)
{
  return Number(1.0 / determinant(F)) * Contravariant::push_forward(H, F);
}



template <int dim, typename Number>
inline SymmetricTensor<4, dim, Number>
Physics::Transformations::Piola::push_forward(
  const SymmetricTensor<4, dim, Number> &H,
  const Tensor<2, dim, Number> &         F)
{
  return Number(1.0 / determinant(F)) * Contravariant::push_forward(H, F);
}



template <int dim, typename Number>
inline Tensor<1, dim, Number>
Physics::Transformations::Piola::pull_back(const Tensor<1, dim, Number> &v,
                                           const Tensor<2, dim, Number> &F)
{
  return Number(determinant(F)) * Contravariant::pull_back(v, F);
}



template <int dim, typename Number>
inline Tensor<2, dim, Number>
Physics::Transformations::Piola::pull_back(const Tensor<2, dim, Number> &t,
                                           const Tensor<2, dim, Number> &F)
{
  return Number(determinant(F)) * Contravariant::pull_back(t, F);
}



template <int dim, typename Number>
inline SymmetricTensor<2, dim, Number>
Physics::Transformations::Piola::pull_back(
  const SymmetricTensor<2, dim, Number> &t,
  const Tensor<2, dim, Number> &         F)
{
  return Number(determinant(F)) * Contravariant::pull_back(t, F);
}



template <int dim, typename Number>
inline Tensor<4, dim, Number>
Physics::Transformations::Piola::pull_back(const Tensor<4, dim, Number> &h,
                                           const Tensor<2, dim, Number> &F)
{
  return Number(determinant(F)) * Contravariant::pull_back(h, F);
}



template <int dim, typename Number>
inline SymmetricTensor<4, dim, Number>
Physics::Transformations::Piola::pull_back(
  const SymmetricTensor<4, dim, Number> &h,
  const Tensor<2, dim, Number> &         F)
{
  return Number(determinant(F)) * Contravariant::pull_back(h, F);
}



template <int dim, typename Number>
inline Tensor<1, dim, Number>
Physics::Transformations::nansons_formula(const Tensor<1, dim, Number> &N,
                                          const Tensor<2, dim, Number> &F)
{
  return cofactor(F) * N;
}


template <int dim, typename Number>
inline Tensor<1, dim, Number>
Physics::Transformations::basis_transformation(const Tensor<1, dim, Number> &V,
                                               const Tensor<2, dim, Number> &B)
{
  return contract<1, 0>(B, V);
}



template <int dim, typename Number>
inline Tensor<2, dim, Number>
Physics::Transformations::basis_transformation(const Tensor<2, dim, Number> &T,
                                               const Tensor<2, dim, Number> &B)
{
  return contract<1, 0>(B, contract<1, 1>(T, B));
}



template <int dim, typename Number>
inline SymmetricTensor<2, dim, Number>
Physics::Transformations::basis_transformation(
  const SymmetricTensor<2, dim, Number> &T,
  const Tensor<2, dim, Number> &         B)
{
  Tensor<2, dim, Number> tmp_1;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int J = 0; J < dim; ++J)
      // Loop over I but complex.h defines a macro I, so use I_ instead
      for (unsigned int I_ = 0; I_ < dim; ++I_)
        tmp_1[i][J] += B[i][I_] * T[I_][J];

  SymmetricTensor<2, dim, Number> out;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      for (unsigned int J = 0; J < dim; ++J)
        out[i][j] += B[j][J] * tmp_1[i][J];

  return out;
}



template <int dim, typename Number>
inline Tensor<4, dim, Number>
Physics::Transformations::basis_transformation(const Tensor<4, dim, Number> &H,
                                               const Tensor<2, dim, Number> &B)
{
  // This contraction order and indexing might look a bit dubious, so a
  // quick explanation as to what's going on is probably in order:
  //
  // When the contract() function operates on the inner indices, the
  // result has the inner index and outer index transposed, i.e.
  // contract<2,1>(H,F) implies
  // T_{IJLk} = (H_{IJMN} F_{mM}) \delta_{mL} \delta_{Nk}
  // rather than T_{IJkL} (the desired result).
  // So, in effect, contraction of the 3rd (inner) index with F as the
  // second argument results in its transposition with respect to its
  // adjacent neighbor. This is due to the position of the argument F,
  // leading to the free index being on the right hand side of the result.
  // However, given that we can do two transformations from the LHS of H
  // and two from the right we can undo the otherwise erroneous
  // swapping of the outer indices upon application of the second
  // sets of contractions.
  //
  // Note: Its significantly quicker (in 3d) to push forward
  // each index individually
  return contract<1, 1>(
    B, contract<1, 1>(B, contract<2, 1>(contract<2, 1>(H, B), B)));
}



template <int dim, typename Number>
inline SymmetricTensor<4, dim, Number>
Physics::Transformations::basis_transformation(
  const SymmetricTensor<4, dim, Number> &H,
  const Tensor<2, dim, Number> &         B)
{
  // The first and last transformation operations respectively
  // break and recover the symmetry properties of the tensors.
  // We also want to perform a minimal number of operations here
  // and avoid some complications related to the transposition of
  // tensor indices when contracting inner indices using the contract()
  // function. (For an explanation of the contraction operations,
  // please see the note in the equivalent function for standard
  // Tensors.) So what we'll do here is manually perform the first
  // and last contractions that break/recover the tensor symmetries
  // on the inner indices, and use the contract() function only on
  // the outer indices.
  //
  // Note: Its significantly quicker (in 3d) to push forward
  // each index individually

  // Push forward (inner) index 1
  Tensor<4, dim, Number> tmp;
  // Loop over I but complex.h defines a macro I, so use I_ instead
  for (unsigned int I_ = 0; I_ < dim; ++I_)
    for (unsigned int j = 0; j < dim; ++j)
      for (unsigned int K = 0; K < dim; ++K)
        for (unsigned int L = 0; L < dim; ++L)
          for (unsigned int J = 0; J < dim; ++J)
            tmp[I_][j][K][L] += B[j][J] * H[I_][J][K][L];

  // Push forward (outer) indices 0 and 3
  tmp = contract<1, 0>(B, contract<3, 1>(tmp, B));

  // Push forward (inner) index 2
  SymmetricTensor<4, dim, Number> out;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = i; j < dim; ++j)
      for (unsigned int k = 0; k < dim; ++k)
        for (unsigned int l = k; l < dim; ++l)
          for (unsigned int K = 0; K < dim; ++K)
            out[i][j][k][l] += B[k][K] * tmp[i][j][K][l];

  return out;
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif




