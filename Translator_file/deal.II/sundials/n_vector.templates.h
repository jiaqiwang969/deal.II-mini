//include/deal.II-translator/sundials/n_vector.templates_0.txt
//-----------------------------------------------------------
//
//    Copyright (C) 2020 - 2021 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

#ifndef dealii_sundials_n_vector_templates_h
#define dealii_sundials_n_vector_templates_h

#include <deal.II/base/config.h>

#include <deal.II/sundials/n_vector.h>

#ifdef DEAL_II_WITH_SUNDIALS

#  include <deal.II/base/exceptions.h>

#  include <deal.II/lac/block_vector.h>
#  include <deal.II/lac/la_parallel_block_vector.h>
#  include <deal.II/lac/la_parallel_vector.h>
#  include <deal.II/lac/la_vector.h>
#  include <deal.II/lac/petsc_block_vector.h>
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/trilinos_parallel_block_vector.h>
#  include <deal.II/lac/trilinos_vector.h>
#  include <deal.II/lac/vector_memory.h>

#  if DEAL_II_SUNDIALS_VERSION_LT(5, 0, 0)
#    include <deal.II/sundials/sundials_backport.h>
#  endif

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
  namespace internal
  {
    /**
     * 一个内部类，用于存储一个向量的指针，并在必要时管理内存。这个类的对象被用作SUNDIALS
     * N_Vector模块中的`内容'字段。此外，这个类有一个标志，用来存储存储的向量是否应该被视为常量。当对该类的非常量对象调用get()时，会检查该标志，如果该向量实际上是常量，则会抛出一个异常。因此，我们保留了一种
     * "运行时常态正确性"，尽管由于SUNDIALS
     * N_Vector不支持常态，所以静态常态正确性已经丧失。
     *
     */
    template <typename VectorType>
    class NVectorContent
    {
    public:
      /**
       * 用现有的 @p vector.  @param
       * 向量创建一个非所有权的内容
       * 要包裹在这个对象中的基础向量。
       *
       */
      NVectorContent(VectorType *vector);

      /**
       * 用一个现有的const创建一个非所有权的内容  @p vector.
       * 如果使用这个构造函数，只允许通过get()const方法访问。
       * @param  vector 在这个对象中要包裹的基础向量。
       *
       */
      NVectorContent(const VectorType *vector);

      /**
       * 分配一个新的（非const）向量，包裹在一个新的内容对象中。当这个对象被销毁时，这个向量将被自动取消分配。
       * @note 这个构造函数是为SUNDIALS的N_VClone()调用准备的。
       *
       */
      NVectorContent();

      /**
       * 对存储向量的非const访问。只有在使用了不同于NVectorContent(const
       * VectorTypevector)的构造函数时才允许。        @return
       *
       */
      VectorType *
      get();

      /**
       * 对存储向量的常量访问。总是允许的。
       *
       */
      const VectorType *
      get() const;

    private:
      using PointerType =
        std::unique_ptr<VectorType, std::function<void(VectorType *)>>;

      /**
       * 如果调用NVectorContent()，可能用于分配存储的向量内存。
       *
       */
      GrowingVectorMemory<VectorType> mem;

      /**
       * 实际存储的向量内容。
       *
       */
      PointerType vector;

      /**
       * 存储的指针是否被视为常数的标志。如果在构造函数中传递的指针确实是常数，那么它将被抛开，但这个标志将被设置为真。然后，对指针的访问必须检查这个标志是否被正确设置。
       *
       */
      const bool is_const;
    };

    /**
     * 帮助创建一个具有所有操作和给定 @p content.   @param
     * 内容的向量，以附加到N_Vector。      @return
     * 一个新的N_Vector
     *
     */
    template <typename VectorType>
    N_Vector
    create_nvector(NVectorContent<VectorType> *content);

    /**
     * 帮助创建一个具有所有操作但没有内容的空向量。
     * @return  一个新的N_Vector
     *
     */
    template <typename VectorType>
    N_Vector
    create_empty_nvector();

    /**
     * SUNDIALS
     * N_Vector文档指定的所有操作的集合。这些函数被附加到通用的N_Vector结构中。
     *
     */
    namespace NVectorOperations
    {
      N_Vector_ID
      get_vector_id(N_Vector v);

      N_Vector
      clone_empty(N_Vector w);

      template <typename VectorType>
      N_Vector
      clone(N_Vector w);

      template <typename VectorType>
      void
      destroy(N_Vector v);

      template <typename VectorType>
      sunindextype
      get_global_length(N_Vector v);

      template <typename VectorType>
      void
      linear_sum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);

      template <typename VectorType>
      realtype
      dot_product(N_Vector x, N_Vector y);

      template <typename VectorType>
      realtype
      weighted_l2_norm(N_Vector x, N_Vector y);

      template <typename VectorType>
      realtype
      l1_norm(N_Vector x);

      template <typename VectorType>
      void
      elementwise_product(N_Vector x, N_Vector y, N_Vector z);

      template <
        typename VectorType,
        typename std::enable_if_t<!IsBlockVector<VectorType>::value, int> = 0>
      void
      elementwise_div(N_Vector x, N_Vector y, N_Vector z);

      template <
        typename VectorType,
        typename std::enable_if_t<IsBlockVector<VectorType>::value, int> = 0>
      void
      elementwise_div(N_Vector x, N_Vector y, N_Vector z);

      template <
        typename VectorType,
        typename std::enable_if_t<!IsBlockVector<VectorType>::value, int> = 0>
      void
      elementwise_inv(N_Vector x, N_Vector z);

      template <
        typename VectorType,
        typename std::enable_if_t<IsBlockVector<VectorType>::value, int> = 0>
      void
      elementwise_inv(N_Vector x, N_Vector z);

      template <
        typename VectorType,
        typename std::enable_if_t<!IsBlockVector<VectorType>::value, int> = 0>
      void
      elementwise_abs(N_Vector x, N_Vector z);

      template <
        typename VectorType,
        typename std::enable_if_t<IsBlockVector<VectorType>::value, int> = 0>
      void
      elementwise_abs(N_Vector x, N_Vector z);

      template <typename VectorType>
      realtype
      weighted_rms_norm(N_Vector x, N_Vector w);

      template <typename VectorType>
      realtype
      weighted_rms_norm_mask(N_Vector x, N_Vector w, N_Vector mask);

      template <typename VectorType>
      realtype
      max_norm(N_Vector x);

      template <
        typename VectorType,
        typename std::enable_if_t<is_serial_vector<VectorType>::value, int> = 0>
      realtype
      min_element(N_Vector x);

      template <
        typename VectorType,
        typename std::enable_if_t<!is_serial_vector<VectorType>::value &&
                                    !IsBlockVector<VectorType>::value,
                                  int> = 0>
      realtype
      min_element(N_Vector x);

      template <
        typename VectorType,
        typename std::enable_if_t<!is_serial_vector<VectorType>::value &&
                                    IsBlockVector<VectorType>::value,
                                  int> = 0>
      realtype
      min_element(N_Vector x);

      template <typename VectorType>
      void
      scale(realtype c, N_Vector x, N_Vector z);

      template <typename VectorType>
      void
      set_constant(realtype c, N_Vector v);

      template <typename VectorType>
      void
      add_constant(N_Vector x, realtype b, N_Vector z);

      template <
        typename VectorType,
        typename std::enable_if_t<!IsBlockVector<VectorType>::value, int> = 0>
      MPI_Comm
      get_communicator(N_Vector v);

      template <
        typename VectorType,
        typename std::enable_if_t<IsBlockVector<VectorType>::value, int> = 0>
      MPI_Comm
      get_communicator(N_Vector v);

      /**
       * Sundials喜欢一个void*，但我们想在内部使用上述函数的安全类型。
       *
       */
      template <
        typename VectorType,
        typename std::enable_if_t<is_serial_vector<VectorType>::value, int> = 0>
      inline void *
      get_communicator_as_void_ptr(N_Vector v);

      template <typename VectorType,
                typename std::enable_if_t<!is_serial_vector<VectorType>::value,
                                          int> = 0>
      inline void *
      get_communicator_as_void_ptr(N_Vector v);

    } // namespace NVectorOperations
  }   // namespace internal
} // namespace SUNDIALS



template <typename VectorType>
SUNDIALS::internal::NVectorContent<VectorType>::NVectorContent()
  : vector(typename VectorMemory<VectorType>::Pointer(mem))
  , is_const(false)
{}



template <typename VectorType>
SUNDIALS::internal::NVectorContent<VectorType>::NVectorContent(
  VectorType *vector)
  : vector(vector, [](VectorType *) {  /* not owning memory -> don't free*/  })
  , is_const(false)
{}



template <typename VectorType>
SUNDIALS::internal::NVectorContent<VectorType>::NVectorContent(
  const VectorType *vector)
  : vector(const_cast<VectorType *>(vector),
           [](VectorType *) {  /* not owning memory -> don't free*/  })
  , is_const(true)
{}



template <typename VectorType>
VectorType *
SUNDIALS::internal::NVectorContent<VectorType>::get()
{
  AssertThrow(
    !is_const,
    ExcMessage("Tried to access a constant vector content as non-const."
               " This most likely happened because a vector that is passed to a"
               " NVectorView() call should not be const.\n"
               "Alternatively, if you tried to access the vector, use"
               " unwrap_nvector_const() for this vector."));
  return vector.get();
}



template <typename VectorType>
const VectorType *
SUNDIALS::internal::NVectorContent<VectorType>::get() const
{
  return vector.get();
}



template <typename VectorType>
SUNDIALS::internal::NVectorView<VectorType>
SUNDIALS::internal::make_nvector_view(VectorType &vector)
{
  return NVectorView<VectorType>(vector);
}



template <typename VectorType>
VectorType *
SUNDIALS::internal::unwrap_nvector(N_Vector v)
{
  Assert(v != nullptr, ExcInternalError());
  Assert(v->content != nullptr, ExcInternalError());
  auto *pContent = reinterpret_cast<NVectorContent<VectorType> *>(v->content);
  return pContent->get();
}



template <typename VectorType>
const VectorType *
SUNDIALS::internal::unwrap_nvector_const(N_Vector v)
{
  Assert(v != nullptr, ExcInternalError());
  Assert(v->content != nullptr, ExcInternalError());
  auto *pContent =
    reinterpret_cast<const NVectorContent<VectorType> *>(v->content);
  return pContent->get();
}



template <typename VectorType>
SUNDIALS::internal::NVectorView<VectorType>::NVectorView(VectorType &vector)
  : vector_ptr(
      create_nvector(
        new NVectorContent<typename std::remove_const<VectorType>::type>(
          &vector)),
      [](N_Vector v) { N_VDestroy(v); })
{}



template <typename VectorType>
SUNDIALS::internal::NVectorView<VectorType>::operator N_Vector() const
{
  Assert(vector_ptr != nullptr, ExcNotInitialized());
  return vector_ptr.get();
}



template <typename VectorType>
N_Vector SUNDIALS::internal::NVectorView<VectorType>::operator->() const
{
  Assert(vector_ptr != nullptr, ExcNotInitialized());
  return vector_ptr.get();
}



template <typename VectorType>
N_Vector
SUNDIALS::internal::create_nvector(NVectorContent<VectorType> *content)
{
  // Create an N_Vector with operators attached and empty content
  N_Vector v = create_empty_nvector<VectorType>();
  Assert(v != nullptr, ExcInternalError());

  v->content = content;
  Assert(v->content != nullptr, ExcInternalError());
  return (v);
}



N_Vector_ID SUNDIALS::internal::NVectorOperations::get_vector_id(N_Vector)
{
  return SUNDIALS_NVEC_CUSTOM;
}



N_Vector
SUNDIALS::internal::NVectorOperations::clone_empty(N_Vector w)
{
  Assert(w != nullptr, ExcInternalError());
  N_Vector v = N_VNewEmpty();
  Assert(v != nullptr, ExcInternalError());

  int status = N_VCopyOps(w, v);
  Assert(status == 0, ExcInternalError());
  (void)status;

  return v;
}



template <typename VectorType>
N_Vector
SUNDIALS::internal::NVectorOperations::clone(N_Vector w)
{
  N_Vector v = clone_empty(w);

  // the corresponding delete is called in destroy()
  auto  cloned   = new NVectorContent<VectorType>();
  auto *w_dealii = unwrap_nvector_const<VectorType>(w);

  // reinit the cloned vector based on the layout of the source vector
  cloned->get()->reinit(*w_dealii);
  v->content = cloned;

  return v;
}



template <typename VectorType>
void
SUNDIALS::internal::NVectorOperations::destroy(N_Vector v)
{
  // support destroying a nullptr because SUNDIALS vectors also do it
  if (v == nullptr)
    return;

  if (v->content != nullptr)
    {
      auto *content =
        reinterpret_cast<NVectorContent<VectorType> *>(v->content);
      // the NVectorContent knows if it owns the memory and will free
      // correctly
      delete content;
      v->content = nullptr;
    }

  N_VFreeEmpty(v);
}



template <typename VectorType,
          std::enable_if_t<IsBlockVector<VectorType>::value, int>>
MPI_Comm
SUNDIALS::internal::NVectorOperations::get_communicator(N_Vector v)
{
  return unwrap_nvector_const<VectorType>(v)->block(0).get_mpi_communicator();
}



template <typename VectorType,
          std::enable_if_t<!IsBlockVector<VectorType>::value, int>>
MPI_Comm
SUNDIALS::internal::NVectorOperations::get_communicator(N_Vector v)
{
  return unwrap_nvector_const<VectorType>(v)->get_mpi_communicator();
}



template <typename VectorType,
          typename std::enable_if_t<is_serial_vector<VectorType>::value, int>>
void *
  SUNDIALS::internal::NVectorOperations::get_communicator_as_void_ptr(N_Vector)
{
  // required by SUNDIALS: MPI-unaware vectors should return the nullptr as comm
  return nullptr;
}



template <typename VectorType,
          typename std::enable_if_t<!is_serial_vector<VectorType>::value, int>>
void *
SUNDIALS::internal::NVectorOperations::get_communicator_as_void_ptr(N_Vector v)
{
#  ifndef DEAL_II_WITH_MPI
  (void)v;
  return nullptr;
#  else
  return get_communicator<VectorType>(v);
#  endif
}



template <typename VectorType>
sunindextype
SUNDIALS::internal::NVectorOperations::get_global_length(N_Vector v)
{
  return unwrap_nvector_const<VectorType>(v)->size();
}



template <typename VectorType>
void
SUNDIALS::internal::NVectorOperations::linear_sum(realtype a,
                                                  N_Vector x,
                                                  realtype b,
                                                  N_Vector y,
                                                  N_Vector z)
{
  auto *x_dealii = unwrap_nvector_const<VectorType>(x);
  auto *y_dealii = unwrap_nvector_const<VectorType>(y);
  auto *z_dealii = unwrap_nvector<VectorType>(z);

  if (z_dealii == x_dealii)
    z_dealii->sadd(a, b, *y_dealii);
  else if (z_dealii == y_dealii)
    z_dealii->sadd(b, a, *x_dealii);
  else
    {
      *z_dealii = 0;
      z_dealii->add(a, *x_dealii, b, *y_dealii);
    }
}



template <typename VectorType>
realtype
SUNDIALS::internal::NVectorOperations::dot_product(N_Vector x, N_Vector y)
{
  return *unwrap_nvector_const<VectorType>(x) *
         *unwrap_nvector_const<VectorType>(y);
}



template <typename VectorType>
realtype
SUNDIALS::internal::NVectorOperations::weighted_l2_norm(N_Vector x, N_Vector w)
{
  // TODO copy can be avoided by a custom kernel
  VectorType tmp      = *unwrap_nvector_const<VectorType>(x);
  auto *     w_dealii = unwrap_nvector_const<VectorType>(w);
  tmp.scale(*w_dealii);
  return tmp.l2_norm();
}



template <typename VectorType>
realtype
SUNDIALS::internal::NVectorOperations::l1_norm(N_Vector x)
{
  return unwrap_nvector_const<VectorType>(x)->l1_norm();
}



template <typename VectorType>
void
SUNDIALS::internal::NVectorOperations::set_constant(realtype c, N_Vector v)
{
  auto *v_dealii = unwrap_nvector<VectorType>(v);
  *v_dealii      = c;
}



template <typename VectorType>
void
SUNDIALS::internal::NVectorOperations::add_constant(N_Vector x,
                                                    realtype b,
                                                    N_Vector z)
{
  auto *x_dealii = unwrap_nvector_const<VectorType>(x);
  auto *z_dealii = unwrap_nvector<VectorType>(z);

  if (x_dealii == z_dealii)
    z_dealii->add(b);
  else
    {
      *z_dealii = *x_dealii;
      z_dealii->add(b);
    }
}



template <typename VectorType>
void
SUNDIALS::internal::NVectorOperations::elementwise_product(N_Vector x,
                                                           N_Vector y,
                                                           N_Vector z)
{
  auto *x_dealii = unwrap_nvector_const<VectorType>(x);
  auto *y_dealii = unwrap_nvector_const<VectorType>(y);
  auto *z_dealii = unwrap_nvector<VectorType>(z);
  if (z_dealii == x_dealii)
    z_dealii->scale(*y_dealii);
  else if (z_dealii == y_dealii)
    z_dealii->scale(*x_dealii);
  else
    {
      *z_dealii = *y_dealii;
      z_dealii->scale(*x_dealii);
    }
}



template <typename VectorType>
realtype
SUNDIALS::internal::NVectorOperations::weighted_rms_norm(N_Vector x, N_Vector w)
{
  // TODO copy can be avoided by a custom kernel
  VectorType tmp      = *unwrap_nvector_const<VectorType>(x);
  auto *     w_dealii = unwrap_nvector_const<VectorType>(w);
  const auto n        = tmp.size();
  tmp.scale(*w_dealii);
  return tmp.l2_norm() / std::sqrt(n);
}



template <typename VectorType>
realtype
SUNDIALS::internal::NVectorOperations::weighted_rms_norm_mask(N_Vector x,
                                                              N_Vector w,
                                                              N_Vector mask)
{
  // TODO copy can be avoided by a custom kernel
  VectorType tmp         = *unwrap_nvector_const<VectorType>(x);
  auto *     w_dealii    = unwrap_nvector_const<VectorType>(w);
  auto *     mask_dealii = unwrap_nvector_const<VectorType>(mask);
  const auto n           = tmp.size();
  tmp.scale(*w_dealii);
  tmp.scale(*mask_dealii);
  return tmp.l2_norm() / std::sqrt(n);
}



template <typename VectorType>
realtype
SUNDIALS::internal::NVectorOperations::max_norm(N_Vector x)
{
  return unwrap_nvector_const<VectorType>(x)->linfty_norm();
}



template <typename VectorType,
          typename std::enable_if_t<is_serial_vector<VectorType>::value, int>>
realtype
SUNDIALS::internal::NVectorOperations::min_element(N_Vector x)
{
  auto *vector = unwrap_nvector_const<VectorType>(x);
  return *std::min_element(vector->begin(), vector->end());
}



template <typename VectorType,
          typename std::enable_if_t<!is_serial_vector<VectorType>::value &&
                                      !IsBlockVector<VectorType>::value,
                                    int>>
realtype
SUNDIALS::internal::NVectorOperations::min_element(N_Vector x)
{
  auto *vector = unwrap_nvector_const<VectorType>(x);


  const auto indexed_less_than = [&](const IndexSet::size_type idxa,
                                     const IndexSet::size_type idxb) {
    return (*vector)[idxa] < (*vector)[idxb];
  };

  auto local_elements = vector->locally_owned_elements();

  const auto local_min = *std::min_element(local_elements.begin(),
                                           local_elements.end(),
                                           indexed_less_than);
  return Utilities::MPI::min((*vector)[local_min],
                             get_communicator<VectorType>(x));
}



template <typename VectorType,
          typename std::enable_if_t<!is_serial_vector<VectorType>::value &&
                                      IsBlockVector<VectorType>::value,
                                    int>>
realtype
SUNDIALS::internal::NVectorOperations::min_element(N_Vector x)
{
  auto *vector = unwrap_nvector_const<VectorType>(x);

  // initialize local minimum to the largest possible value
  auto proc_local_min =
    std::numeric_limits<typename VectorType::value_type>::max();

  for (unsigned i = 0; i < vector->n_blocks(); ++i)
    {
      const auto indexed_less_than = [&](const IndexSet::size_type idxa,
                                         const IndexSet::size_type idxb) {
        return vector->block(i)[idxa] < vector->block(i)[idxb];
      };

      auto local_elements = vector->block(i).locally_owned_elements();

      const auto block_local_min_element =
        std::min_element(local_elements.begin(),
                         local_elements.end(),
                         indexed_less_than);

      // guard against empty blocks on this processor
      if (block_local_min_element != local_elements.end())
        proc_local_min =
          std::min(proc_local_min, vector->block(i)[*block_local_min_element]);
    }

  return Utilities::MPI::min(proc_local_min, get_communicator<VectorType>(x));
}



template <typename VectorType>
void
SUNDIALS::internal::NVectorOperations::scale(realtype c, N_Vector x, N_Vector z)
{
  auto *x_dealii = unwrap_nvector_const<VectorType>(x);
  auto *z_dealii = unwrap_nvector<VectorType>(z);

  if (x_dealii == z_dealii)
    (*z_dealii) *= c;
  else
    z_dealii->sadd(0.0, c, *x_dealii);
}



template <typename VectorType,
          std::enable_if_t<!IsBlockVector<VectorType>::value, int>>
void
SUNDIALS::internal::NVectorOperations::elementwise_div(N_Vector x,
                                                       N_Vector y,
                                                       N_Vector z)
{
  auto *x_dealii = unwrap_nvector_const<VectorType>(x);
  auto *y_dealii = unwrap_nvector_const<VectorType>(y);
  auto *z_dealii = unwrap_nvector<VectorType>(z);

  AssertDimension(x_dealii->size(), z_dealii->size());
  AssertDimension(x_dealii->size(), y_dealii->size());

  auto x_ele = x_dealii->locally_owned_elements();
  for (const auto idx : x_ele)
    {
      (*z_dealii)[idx] = (*x_dealii)[idx] / (*y_dealii)[idx];
    }
  z_dealii->compress(VectorOperation::insert);
}



template <typename VectorType,
          std::enable_if_t<IsBlockVector<VectorType>::value, int>>
void
SUNDIALS::internal::NVectorOperations::elementwise_div(N_Vector x,
                                                       N_Vector y,
                                                       N_Vector z)
{
  auto *x_dealii = unwrap_nvector_const<VectorType>(x);
  auto *y_dealii = unwrap_nvector_const<VectorType>(y);
  auto *z_dealii = unwrap_nvector<VectorType>(z);

  AssertDimension(x_dealii->size(), z_dealii->size());
  AssertDimension(x_dealii->size(), y_dealii->size());
  AssertDimension(x_dealii->n_blocks(), z_dealii->n_blocks());
  AssertDimension(x_dealii->n_blocks(), y_dealii->n_blocks());

  for (unsigned i = 0; i < x_dealii->n_blocks(); ++i)
    {
      auto x_ele = x_dealii->block(i).locally_owned_elements();
      for (const auto idx : x_ele)
        {
          z_dealii->block(i)[idx] =
            x_dealii->block(i)[idx] / y_dealii->block(i)[idx];
        }
    }
  z_dealii->compress(VectorOperation::insert);
}



template <typename VectorType,
          std::enable_if_t<!IsBlockVector<VectorType>::value, int>>
void
SUNDIALS::internal::NVectorOperations::elementwise_inv(N_Vector x, N_Vector z)
{
  auto *x_dealii = unwrap_nvector_const<VectorType>(x);
  auto *z_dealii = unwrap_nvector<VectorType>(z);

  AssertDimension(x_dealii->size(), z_dealii->size());

  auto x_ele = x_dealii->locally_owned_elements();
  for (const auto idx : x_ele)
    {
      (*z_dealii)[idx] = 1.0 / (*x_dealii)[idx];
    }
  z_dealii->compress(VectorOperation::insert);
}



template <typename VectorType,
          std::enable_if_t<IsBlockVector<VectorType>::value, int>>
void
SUNDIALS::internal::NVectorOperations::elementwise_inv(N_Vector x, N_Vector z)
{
  auto *x_dealii = unwrap_nvector_const<VectorType>(x);
  auto *z_dealii = unwrap_nvector<VectorType>(z);

  AssertDimension(x_dealii->size(), z_dealii->size());
  AssertDimension(x_dealii->n_blocks(), z_dealii->n_blocks());

  for (unsigned i = 0; i < x_dealii->n_blocks(); ++i)
    {
      auto x_ele = x_dealii->block(i).locally_owned_elements();
      for (const auto idx : x_ele)
        {
          z_dealii->block(i)[idx] = 1.0 / x_dealii->block(i)[idx];
        }
    }

  z_dealii->compress(VectorOperation::insert);
}



template <typename VectorType,
          std::enable_if_t<!IsBlockVector<VectorType>::value, int>>
void
SUNDIALS::internal::NVectorOperations::elementwise_abs(N_Vector x, N_Vector z)
{
  auto *x_dealii = unwrap_nvector_const<VectorType>(x);
  auto *z_dealii = unwrap_nvector<VectorType>(z);

  AssertDimension(x_dealii->size(), z_dealii->size());

  auto x_ele = x_dealii->locally_owned_elements();
  for (const auto idx : x_ele)
    {
      (*z_dealii)[idx] = std::fabs((*x_dealii)[idx]);
    }

  z_dealii->compress(VectorOperation::insert);
}



template <typename VectorType,
          std::enable_if_t<IsBlockVector<VectorType>::value, int>>
void
SUNDIALS::internal::NVectorOperations::elementwise_abs(N_Vector x, N_Vector z)
{
  auto *x_dealii = unwrap_nvector_const<VectorType>(x);
  auto *z_dealii = unwrap_nvector<VectorType>(z);

  AssertDimension(x_dealii->size(), z_dealii->size());
  AssertDimension(x_dealii->n_blocks(), z_dealii->n_blocks());

  for (unsigned i = 0; i < x_dealii->n_blocks(); ++i)
    {
      auto x_ele = x_dealii->block(i).locally_owned_elements();
      for (const auto idx : x_ele)
        {
          z_dealii->block(i)[idx] = std::fabs(x_dealii->block(i)[idx]);
        }
    }

  z_dealii->compress(VectorOperation::insert);
}



template <typename VectorType>
N_Vector
SUNDIALS::internal::create_empty_nvector()
{
  N_Vector v = N_VNewEmpty();
  Assert(v != nullptr, ExcInternalError());

   /* constructors, destructors, and utility operations */ 
  v->ops->nvgetvectorid = NVectorOperations::get_vector_id;
  v->ops->nvclone       = NVectorOperations::clone<VectorType>;
  v->ops->nvcloneempty  = NVectorOperations::clone_empty;
  v->ops->nvdestroy     = NVectorOperations::destroy<VectorType>;
  //  v->ops->nvspace           = undef;
#  if DEAL_II_SUNDIALS_VERSION_GTE(5, 0, 0)
  v->ops->nvgetcommunicator =
    NVectorOperations::get_communicator_as_void_ptr<VectorType>;
  v->ops->nvgetlength = NVectorOperations::get_global_length<VectorType>;
#  endif

   /* standard vector operations */ 
  v->ops->nvlinearsum = NVectorOperations::linear_sum<VectorType>;
  v->ops->nvconst     = NVectorOperations::set_constant<VectorType>;
  v->ops->nvprod      = NVectorOperations::elementwise_product<VectorType>;
  v->ops->nvdiv       = NVectorOperations::elementwise_div<VectorType>;
  v->ops->nvscale     = NVectorOperations::scale<VectorType>;
  v->ops->nvabs       = NVectorOperations::elementwise_abs<VectorType>;
  v->ops->nvinv       = NVectorOperations::elementwise_inv<VectorType>;
  v->ops->nvaddconst  = NVectorOperations::add_constant<VectorType>;
  v->ops->nvdotprod   = NVectorOperations::dot_product<VectorType>;
  v->ops->nvmaxnorm   = NVectorOperations::max_norm<VectorType>;
  v->ops->nvwrmsnorm  = NVectorOperations::weighted_rms_norm<VectorType>;
  v->ops->nvmin       = NVectorOperations::min_element<VectorType>;
  v->ops->nvwl2norm   = NVectorOperations::weighted_l2_norm<VectorType>;
  v->ops->nvl1norm    = NVectorOperations::l1_norm<VectorType>;
  v->ops->nvwrmsnormmask =
    NVectorOperations::weighted_rms_norm_mask<VectorType>;
  //  v->ops->nvcompare      = undef;
  //  v->ops->nvinvtest      = undef;
  //  v->ops->nvconstrmask   = undef;
  //  v->ops->nvminquotient  = undef;

   /* fused and vector array operations are disabled (NULL) by default */ 

   /* local reduction operations */ 
  //  v->ops->nvdotprodlocal     = undef;
  //  v->ops->nvmaxnormlocal     = undef;
  //  v->ops->nvminlocal         = undef;
  //  v->ops->nvl1normlocal      = undef;
  //  v->ops->nvinvtestlocal     = undef;
  //  v->ops->nvconstrmasklocal  = undef;
  //  v->ops->nvminquotientlocal = undef;
  //  v->ops->nvwsqrsumlocal     = undef;
  //  v->ops->nvwsqrsummasklocal = undef;
  return (v);
}

DEAL_II_NAMESPACE_CLOSE

#endif
#endif


