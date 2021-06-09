//include/deal.II-translator/distributed/p4est_wrappers_0.txt
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

#ifndef dealii_p4est_wrappers_h
#define dealii_p4est_wrappers_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>

#ifdef DEAL_II_WITH_P4EST
#  include <p4est_bits.h>
#  include <p4est_communication.h>
#  include <p4est_extended.h>
#  include <p4est_ghost.h>
#  include <p4est_iterate.h>
#  include <p4est_vtk.h>
#  include <p8est_bits.h>
#  include <p8est_communication.h>
#  include <p8est_extended.h>
#  include <p8est_ghost.h>
#  include <p8est_iterate.h>
#  include <p8est_vtk.h>

#  include <map>
#  include <set>


DEAL_II_NAMESPACE_OPEN

// Forward declaration
#  ifndef DOXYGEN
namespace parallel
{
  namespace distributed
  {
    template <int dim, int spacedim>
    class Triangulation;
  }
} // namespace parallel
#  endif

namespace internal
{
  namespace p4est
  {
    /**
     * 一个结构，其明确的特殊化包含了对相关p4est_*和p8est_*类型的别名。使用这个结构，比如说通过说
     * <tt>types<dim>::connectivity</tt>
     * ，我们可以以独立于维度的方式编写代码，根据模板参数，可以引用p4est_connectivity_t或者p8est_connectivity_t。
     *
     */
    template <int>
    struct types;

    // these struct mimics p4est for 1D
    template <>
    struct types<1>
    {
      // id of a quadrant is an integeger
      using quadrant = int;

      // maximum number of children
      static const int max_n_child_indices_bits = 27;

      // number of bits the data type of id has
      static const int n_bits = std::numeric_limits<quadrant>::digits;
    };

    template <>
    struct types<2>
    {
      using connectivity     = p4est_connectivity_t;
      using forest           = p4est_t;
      using tree             = p4est_tree_t;
      using quadrant         = p4est_quadrant_t;
      using topidx           = p4est_topidx_t;
      using locidx           = p4est_locidx_t;
      using gloidx           = p4est_gloidx_t;
      using balance_type     = p4est_connect_type_t;
      using ghost            = p4est_ghost_t;
      using transfer_context = p4est_transfer_context_t;
    };

    template <>
    struct types<3>
    {
      using connectivity     = p8est_connectivity_t;
      using forest           = p8est_t;
      using tree             = p8est_tree_t;
      using quadrant         = p8est_quadrant_t;
      using topidx           = p4est_topidx_t;
      using locidx           = p4est_locidx_t;
      using gloidx           = p4est_gloidx_t;
      using balance_type     = p8est_connect_type_t;
      using ghost            = p8est_ghost_t;
      using transfer_context = p8est_transfer_context_t;
    };



    /**
     * 一个结构，其明确的专业化包含指向相关p4est_*和p8est_*函数的指针。使用这个结构，例如说
     * functions<dim>::quadrant_compare,
     * 我们可以以独立于维度的方式编写代码，可以调用p4est_quadrant_compare或者p8est_quadrant_compare，取决于模板参数。
     *
     */
    template <int dim>
    struct functions;

    template <>
    struct functions<2>
    {
      static int (&quadrant_compare)(const void *v1, const void *v2);

      static void (&quadrant_childrenv)(const types<2>::quadrant *q,
                                        types<2>::quadrant        c[]);

      static int (&quadrant_overlaps_tree)(types<2>::tree *          tree,
                                           const types<2>::quadrant *q);

      static void (&quadrant_set_morton)(types<2>::quadrant *quadrant,
                                         int                 level,
                                         uint64_t            id);

      static int (&quadrant_is_equal)(const types<2>::quadrant *q1,
                                      const types<2>::quadrant *q2);

      static int (&quadrant_is_sibling)(const types<2>::quadrant *q1,
                                        const types<2>::quadrant *q2);

      static int (&quadrant_is_ancestor)(const types<2>::quadrant *q1,
                                         const types<2>::quadrant *q2);

      static int (&quadrant_ancestor_id)(const types<2>::quadrant *q,
                                         int                       level);

      static int (&comm_find_owner)(types<2>::forest *        p4est,
                                    const types<2>::locidx    which_tree,
                                    const types<2>::quadrant *q,
                                    const int                 guess);

      static types<2>::connectivity *(&connectivity_new)(
        types<2>::topidx num_vertices,
        types<2>::topidx num_trees,
        types<2>::topidx num_corners,
        types<2>::topidx num_vtt);

      static types<2>::connectivity *(&connectivity_new_copy)(
        types<2>::topidx        num_vertices,
        types<2>::topidx        num_trees,
        types<2>::topidx        num_corners,
        const double *          vertices,
        const types<2>::topidx *ttv,
        const types<2>::topidx *ttt,
        const int8_t *          ttf,
        const types<2>::topidx *ttc,
        const types<2>::topidx *coff,
        const types<2>::topidx *ctt,
        const int8_t *          ctc);

      static void (&connectivity_join_faces)(types<2>::connectivity *conn,
                                             types<2>::topidx        tree_left,
                                             types<2>::topidx        tree_right,
                                             int                     face_left,
                                             int                     face_right,
                                             int orientation);



      static void (&connectivity_destroy)(p4est_connectivity_t *connectivity);

      static types<2>::forest *(&new_forest)(
        MPI_Comm                mpicomm,
        types<2>::connectivity *connectivity,
        types<2>::locidx        min_quadrants,
        int                     min_level,
        int                     fill_uniform,
        std::size_t             data_size,
        p4est_init_t            init_fn,
        void *                  user_pointer);

      static types<2>::forest *(&copy_forest)(types<2>::forest *input,
                                              int               copy_data);

      static void (&destroy)(types<2>::forest *p4est);

      static void (&refine)(types<2>::forest *p4est,
                            int               refine_recursive,
                            p4est_refine_t    refine_fn,
                            p4est_init_t      init_fn);

      static void (&coarsen)(types<2>::forest *p4est,
                             int               coarsen_recursive,
                             p4est_coarsen_t   coarsen_fn,
                             p4est_init_t      init_fn);

      static void (&balance)(types<2>::forest *     p4est,
                             types<2>::balance_type btype,
                             p4est_init_t           init_fn);

      static types<2>::gloidx (&partition)(types<2>::forest *p4est,
                                           int partition_for_coarsening,
                                           p4est_weight_t weight_fn);

      static void (&save)(const char *      filename,
                          types<2>::forest *p4est,
                          int               save_data);

      static types<2>::forest *(&load_ext)(const char *filename,
                                           MPI_Comm    mpicomm,
                                           std::size_t data_size,
                                           int         load_data,
                                           int         autopartition,
                                           int         broadcasthead,
                                           void *      user_pointer,
                                           types<2>::connectivity **p4est);

      static int (&connectivity_save)(const char *            filename,
                                      types<2>::connectivity *connectivity);

      static int (&connectivity_is_valid)(types<2>::connectivity *connectivity);

      static types<2>::connectivity *(&connectivity_load)(const char * filename,
                                                          std::size_t *length);

      static unsigned int (&checksum)(types<2>::forest *p4est);

      static void (&vtk_write_file)(types<2>::forest *p4est,
                                    p4est_geometry_t *,
                                    const char *baseName);

      static types<2>::ghost *(&ghost_new)(types<2>::forest *     p4est,
                                           types<2>::balance_type btype);

      static void (&ghost_destroy)(types<2>::ghost *ghost);

      static void (&reset_data)(types<2>::forest *p4est,
                                std::size_t       data_size,
                                p4est_init_t      init_fn,
                                void *            user_pointer);

      static std::size_t (&forest_memory_used)(types<2>::forest *p4est);

      static std::size_t (&connectivity_memory_used)(
        types<2>::connectivity *p4est);

      template <int spacedim>
      static void
        iterate(dealii::internal::p4est::types<2>::forest *parallel_forest,
                dealii::internal::p4est::types<2>::ghost * parallel_ghost,
                void *                                     user_data);

      static constexpr unsigned int max_level = P4EST_MAXLEVEL;

      static void (&transfer_fixed)(const types<2>::gloidx *dest_gfq,
                                    const types<2>::gloidx *src_gfq,
                                    MPI_Comm                mpicomm,
                                    int                     tag,
                                    void *                  dest_data,
                                    const void *            src_data,
                                    std::size_t             data_size);

      static types<2>::transfer_context *(&transfer_fixed_begin)(
        const types<2>::gloidx *dest_gfq,
        const types<2>::gloidx *src_gfq,
        MPI_Comm                mpicomm,
        int                     tag,
        void *                  dest_data,
        const void *            src_data,
        std::size_t             data_size);

      static void (&transfer_fixed_end)(types<2>::transfer_context *tc);

      static void (&transfer_custom)(const types<2>::gloidx *dest_gfq,
                                     const types<2>::gloidx *src_gfq,
                                     MPI_Comm                mpicomm,
                                     int                     tag,
                                     void *                  dest_data,
                                     const int *             dest_sizes,
                                     const void *            src_data,
                                     const int *             src_sizes);

      static types<2>::transfer_context *(&transfer_custom_begin)(
        const types<2>::gloidx *dest_gfq,
        const types<2>::gloidx *src_gfq,
        MPI_Comm                mpicomm,
        int                     tag,
        void *                  dest_data,
        const int *             dest_sizes,
        const void *            src_data,
        const int *             src_sizes);

      static void (&transfer_custom_end)(types<2>::transfer_context *tc);
    };


    template <>
    struct functions<3>
    {
      static int (&quadrant_compare)(const void *v1, const void *v2);

      static void (&quadrant_childrenv)(const types<3>::quadrant *q,
                                        types<3>::quadrant        c[]);

      static int (&quadrant_overlaps_tree)(types<3>::tree *          tree,
                                           const types<3>::quadrant *q);

      static void (&quadrant_set_morton)(types<3>::quadrant *quadrant,
                                         int                 level,
                                         uint64_t            id);

      static int (&quadrant_is_equal)(const types<3>::quadrant *q1,
                                      const types<3>::quadrant *q2);

      static int (&quadrant_is_sibling)(const types<3>::quadrant *q1,
                                        const types<3>::quadrant *q2);

      static int (&quadrant_is_ancestor)(const types<3>::quadrant *q1,
                                         const types<3>::quadrant *q2);
      static int (&quadrant_ancestor_id)(const types<3>::quadrant *q,
                                         int                       level);

      static int (&comm_find_owner)(types<3>::forest *        p4est,
                                    const types<3>::locidx    which_tree,
                                    const types<3>::quadrant *q,
                                    const int                 guess);

      static types<3>::connectivity *(&connectivity_new)(
        types<3>::topidx num_vertices,
        types<3>::topidx num_trees,
        types<3>::topidx num_edges,
        types<3>::topidx num_ett,
        types<3>::topidx num_corners,
        types<3>::topidx num_ctt);

      static types<3>::connectivity *(&connectivity_new_copy)(
        types<3>::topidx        num_vertices,
        types<3>::topidx        num_trees,
        types<3>::topidx        num_edges,
        types<3>::topidx        num_corners,
        const double *          vertices,
        const types<3>::topidx *ttv,
        const types<3>::topidx *ttt,
        const int8_t *          ttf,
        const types<3>::topidx *tte,
        const types<3>::topidx *eoff,
        const types<3>::topidx *ett,
        const int8_t *          ete,
        const types<3>::topidx *ttc,
        const types<3>::topidx *coff,
        const types<3>::topidx *ctt,
        const int8_t *          ctc);

      static void (&connectivity_join_faces)(types<3>::connectivity *conn,
                                             types<3>::topidx        tree_left,
                                             types<3>::topidx        tree_right,
                                             int                     face_left,
                                             int                     face_right,
                                             int orientation);

      static void (&connectivity_destroy)(p8est_connectivity_t *connectivity);

      static types<3>::forest *(&new_forest)(
        MPI_Comm                mpicomm,
        types<3>::connectivity *connectivity,
        types<3>::locidx        min_quadrants,
        int                     min_level,
        int                     fill_uniform,
        std::size_t             data_size,
        p8est_init_t            init_fn,
        void *                  user_pointer);

      static types<3>::forest *(&copy_forest)(types<3>::forest *input,
                                              int               copy_data);

      static void (&destroy)(types<3>::forest *p8est);

      static void (&refine)(types<3>::forest *p8est,
                            int               refine_recursive,
                            p8est_refine_t    refine_fn,
                            p8est_init_t      init_fn);

      static void (&coarsen)(types<3>::forest *p8est,
                             int               coarsen_recursive,
                             p8est_coarsen_t   coarsen_fn,
                             p8est_init_t      init_fn);

      static void (&balance)(types<3>::forest *     p8est,
                             types<3>::balance_type btype,
                             p8est_init_t           init_fn);

      static types<3>::gloidx (&partition)(types<3>::forest *p8est,
                                           int partition_for_coarsening,
                                           p8est_weight_t weight_fn);

      static void (&save)(const char *      filename,
                          types<3>::forest *p4est,
                          int               save_data);

      static types<3>::forest *(&load_ext)(const char *filename,
                                           MPI_Comm    mpicomm,
                                           std::size_t data_size,
                                           int         load_data,
                                           int         autopartition,
                                           int         broadcasthead,
                                           void *      user_pointer,
                                           types<3>::connectivity **p4est);

      static int (&connectivity_save)(const char *            filename,
                                      types<3>::connectivity *connectivity);

      static int (&connectivity_is_valid)(types<3>::connectivity *connectivity);

      static types<3>::connectivity *(&connectivity_load)(const char * filename,
                                                          std::size_t *length);

      static unsigned int (&checksum)(types<3>::forest *p8est);

      static void (&vtk_write_file)(types<3>::forest *p8est,
                                    p8est_geometry_t *,
                                    const char *baseName);
      static types<3>::ghost *(&ghost_new)(types<3>::forest *     p4est,
                                           types<3>::balance_type btype);

      static void (&ghost_destroy)(types<3>::ghost *ghost);

      static void (&reset_data)(types<3>::forest *p4est,
                                std::size_t       data_size,
                                p8est_init_t      init_fn,
                                void *            user_pointer);

      static std::size_t (&forest_memory_used)(types<3>::forest *p4est);

      static std::size_t (&connectivity_memory_used)(
        types<3>::connectivity *p4est);

      static constexpr unsigned int max_level = P8EST_MAXLEVEL;

      static void (&transfer_fixed)(const types<3>::gloidx *dest_gfq,
                                    const types<3>::gloidx *src_gfq,
                                    MPI_Comm                mpicomm,
                                    int                     tag,
                                    void *                  dest_data,
                                    const void *            src_data,
                                    std::size_t             data_size);

      static types<3>::transfer_context *(&transfer_fixed_begin)(
        const types<3>::gloidx *dest_gfq,
        const types<3>::gloidx *src_gfq,
        MPI_Comm                mpicomm,
        int                     tag,
        void *                  dest_data,
        const void *            src_data,
        std::size_t             data_size);

      static void (&transfer_fixed_end)(types<3>::transfer_context *tc);

      static void (&transfer_custom)(const types<3>::gloidx *dest_gfq,
                                     const types<3>::gloidx *src_gfq,
                                     MPI_Comm                mpicomm,
                                     int                     tag,
                                     void *                  dest_data,
                                     const int *             dest_sizes,
                                     const void *            src_data,
                                     const int *             src_sizes);

      static types<3>::transfer_context *(&transfer_custom_begin)(
        const types<3>::gloidx *dest_gfq,
        const types<3>::gloidx *src_gfq,
        MPI_Comm                mpicomm,
        int                     tag,
        void *                  dest_data,
        const int *             dest_sizes,
        const void *            src_data,
        const int *             src_sizes);

      static void (&transfer_custom_end)(types<3>::transfer_context *tc);
    };



    /**
     * 该结构模板化了p4est迭代结构和函数原型，用于执行面、边和角的回调函数，这些函数需要本地邻域信息，即相邻的单元格。
     *
     */
    template <int dim>
    struct iter;

    template <>
    struct iter<2>
    {
      using corner_info = p4est_iter_corner_info_t;
      using corner_side = p4est_iter_corner_side_t;
      using corner_iter = p4est_iter_corner_t;
      using face_info   = p4est_iter_face_info_t;
      using face_side   = p4est_iter_face_side_t;
      using face_iter   = p4est_iter_face_t;
    };

    template <>
    struct iter<3>
    {
      using corner_info = p8est_iter_corner_info_t;
      using corner_side = p8est_iter_corner_side_t;
      using corner_iter = p8est_iter_corner_t;
      using edge_info   = p8est_iter_edge_info_t;
      using edge_side   = p8est_iter_edge_side_t;
      using edge_iter   = p8est_iter_edge_t;
      using face_info   = p8est_iter_face_info_t;
      using face_side   = p8est_iter_face_side_t;
      using face_iter   = p8est_iter_face_t;
    };



    /**
     * 初始化单元格p4est_cell的
     * GeometryInfo<dim>::max_children_per_cell 子女。
     *
     */
    template <int dim>
    void
    init_quadrant_children(
      const typename types<dim>::quadrant &p4est_cell,
      typename types<dim>::quadrant (
        &p4est_children)[dealii::GeometryInfo<dim>::max_children_per_cell]);



    /**
     * 初始化象限以代表粗略的单元。
     *
     */
    template <int dim>
    void
    init_coarse_quadrant(typename types<dim>::quadrant &quad);



    /**
     * 返回q1和q2是否相等
     *
     */
    template <int dim>
    bool
    quadrant_is_equal(const typename types<dim>::quadrant &q1,
                      const typename types<dim>::quadrant &q2);



    /**
     * 返回q1是否是q2的祖先
     *
     */
    template <int dim>
    bool
    quadrant_is_ancestor(const typename types<dim>::quadrant &q1,
                         const typename types<dim>::quadrant &q2);



    /**
     * 返回粗略的单元格的子代是否被存储在本地
     *
     */
    template <int dim>
    bool
    tree_exists_locally(const typename types<dim>::forest *parallel_forest,
                        const typename types<dim>::topidx  coarse_grid_cell);


    /**
     * 深度复制一个p4est连接对象。
     *
     */
    template <int dim>
    typename types<dim>::connectivity *
    copy_connectivity(const typename types<dim>::connectivity *connectivity);

#  ifndef DOXYGEN
    template <>
    typename types<2>::connectivity *
    copy_connectivity<2>(const typename types<2>::connectivity *connectivity);

    template <>
    typename types<3>::connectivity *
    copy_connectivity<3>(const typename types<3>::connectivity *connectivity);
#  endif
  } // namespace p4est
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_P4EST

#endif // dealii_p4est_wrappers_h


