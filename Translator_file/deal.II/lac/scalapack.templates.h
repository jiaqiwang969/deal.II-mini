//include/deal.II-translator/lac/scalapack.templates_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2019 by the deal.II authors
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

#ifndef dealii_scalapack_templates_h
#define dealii_scalapack_templates_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SCALAPACK

#  include <deal.II/base/mpi.h>
#  include <deal.II/base/mpi.templates.h>

#  ifdef DEAL_II_HAVE_FP_EXCEPTIONS
#    include <cfenv>
#  endif

// useful examples:
// https://stackoverflow.com/questions/14147705/cholesky-decomposition-scalapack-error/14203864
// http://icl.cs.utk.edu/lapack-forum/viewtopic.php?t=139   // second post by
// Julien Langou
// https://andyspiros.wordpress.com/2011/07/08/an-example-of-blacs-with-c/
// http://qboxcode.org/trac/browser/qb/tags/rel1_63_4/src/Matrix.C
// https://gitlab.phys.ethz.ch/lwossnig/lecture/blob/a534f562dfb2ad5c564abe5c2356d5d956fb7218/examples/mpi/scalapack.cpp
// https://github.com/elemental/Elemental/blob/master/src/core/imports/scalapack.cpp
// https://scicomp.stackexchange.com/questions/7766/performance-optimization-or-tuning-possible-for-scalapack-gemm
//
// info:
// http://www.netlib.org/scalapack/slug/index.html       // User guide
// http://www.netlib.org/scalapack/slug/node135.html // How to Measure Errors

extern "C"
{
   /* Basic Linear Algebra Communication Subprograms (BLACS) declarations */ 
  // https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dinitb.htm#dinitb

  /**
   * 确定有多少进程可用以及当前的进程等级。https://www.ibm.com/support/knowledgecenter/en/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dbpnf.htm。
   *
   */
  void
  Cblacs_pinfo(int *rank, int *nprocs);

  /**
   * 根据输入的 @p what 和 @p icontxt. 返回 @p val
   * 中的内部BLACS值。最常见的用途是检索默认的系统上下文（
   * @p what =0， @p icontxt
   * 被忽略），用于BLACS_GRIDINIT或BLACS_GRIDMAP。
   * https://www.ibm.com/support/knowledgecenter/en/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dbget.htm
   *
   */
  void
  Cblacs_get(int icontxt, int what, int *val);

  /**
   * 以行为主或列为主的顺序将进程映射到进程网格中。每个进程的输入参数必须是相同的。
   * 在返回时， @p context
   * 是BLACS上下文的整数句柄，而在进入时，它是一个系统上下文，用于创建BLACS上下文。https://www.ibm.com/support/knowledgecenter/en/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dbint.htm
   *
   */
  void
  Cblacs_gridinit(int *       context,
                  const char *order,
                  int         grid_height,
                  int         grid_width);

  /**
   * 返回进程的行和列索引。https://www.ibm.com/support/knowledgecenter/en/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dbinfo.htm
   *
   */
  void
  Cblacs_gridinfo(int  context,
                  int *grid_height,
                  int *grid_width,
                  int *grid_row,
                  int *grid_col);

  /**
   * 给出系统进程号，返回BLACS的进程网格中的行和列坐标。
   *
   */
  void
  Cblacs_pcoord(int ictxt, int pnum, int *prow, int *pcol);

  /**
   * 释放一个BLACS上下文。
   *
   */
  void
  Cblacs_gridexit(int context);

  /**
   * 该例程将暂停执行指定范围内的所有进程，直到它们全部调用该例程。
   *
   */
  void
  Cblacs_barrier(int, const char *);

  /**
   * 释放所有BLACS上下文并释放所有分配的内存。
   *
   */
  void
  Cblacs_exit(int error_code);

  /**
   * 从进程 @prsrc,  @p csrc
   * 接收一个信息，将其转化为一般的矩形矩阵。https://software.intel.com/en-us/mkl-developer-reference-c-gerv2d。
   *
   */
  void
  Cdgerv2d(int context, int M, int N, double *A, int lda, int rsrc, int csrc);
  void
  Csgerv2d(int context, int M, int N, float *A, int lda, int rsrc, int csrc);

  /**
   * 将一般矩形矩阵A发送到进程网格中的目标进程 @p rdest
   * @p cdest
   * 。https://software.intel.com/en-us/mkl-developer-reference-c-2018-beta-gesd2d
   *
   */
  void
  Cdgesd2d(int context, int M, int N, double *A, int lda, int rdest, int cdest);
  void
  Csgesd2d(int context, int M, int N, float *A, int lda, int rdest, int cdest);

  /**
   * 从MPI获取BLACS上下文  @p comm.  。
   *
   */
  int
  Csys2blacs_handle(MPI_Comm comm);

  /**
   * 计算每个进程拥有多少行和列（NUMber of Rows Or
   * Columns）。https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_dnumy.htm
   *
   */
  int
  numroc_(const int *n,
          const int *nb,
          const int *iproc,
          const int *isproc,
          const int *nprocs);

  /**
   * 计算N乘N实数对称正定分布矩阵sub( A
   * )的Cholesky分解，表示A(IA:IA+N-1,
   * JA:JA+N-1)。http://www.netlib.org/scalapack/explore-html/d5/d9e/pdpotrf_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lpotrf.htm
   *
   */
  void
  pdpotrf_(const char *UPLO,
           const int * N,
           double *    A,
           const int * IA,
           const int * JA,
           const int * DESCA,
           int *       INFO);
  void
  pspotrf_(const char *UPLO,
           const int * N,
           float *     A,
           const int * IA,
           const int * JA,
           const int * DESCA,
           int *       INFO);

  /**
   * 计算一般分布式矩阵sub( A
   * )的LU因式分解，使用部分枢轴和行互换。http://www.netlib.org/scalapack/explore-html/df/dfe/pdgetrf_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/en/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lgetrf.htm
   *
   */
  void
  pdgetrf_(const int *m,
           const int *n,
           double *   A,
           const int *IA,
           const int *JA,
           const int *DESCA,
           int *      ipiv,
           int *      INFO);
  void
  psgetrf_(const int *m,
           const int *n,
           float *    A,
           const int *IA,
           const int *JA,
           const int *DESCA,
           int *      ipiv,
           int *      INFO);

  /**
   * 使用PDPOTRF计算的Cholesky因式分解sub( A ) =
   * U**T*U或L*L**T，计算实对称正定分布矩阵sub( A ) =
   * A(IA:IA+N-1,JA:JA+N-1)的逆。http://www.netlib.org/scalapack/explore-html/d2/d44/pdpotri_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lpotri.htm
   * https://software.intel.com/en-us/mkl-developer-reference-c-p-potri
   *
   */
  void
  pdpotri_(const char *UPLO,
           const int * N,
           double *    A,
           const int * IA,
           const int * JA,
           const int * DESCA,
           int *       INFO);
  void
  pspotri_(const char *UPLO,
           const int * N,
           float *     A,
           const int * IA,
           const int * JA,
           const int * DESCA,
           int *       INFO);

  /**
   * PDGETRI使用PDGETRF计算的LU因子化来计算分布式矩阵的逆。该方法反转U，然后通过求解InvA的系统InvA*L=inv(U)，计算sub(
   * A )=A(IA:IA+N-1,JA:JA+N-1)的反，表示InvA。
   * http://www.netlib.org/scalapack/explore-html/d3/df3/pdgetri_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lgetri.htm
   *
   */
  void
  pdgetri_(const int *N,
           double *   A,
           const int *IA,
           const int *JA,
           const int *DESCA,
           const int *ipiv,
           double *   work,
           int *      lwork,
           int *      iwork,
           int *      liwork,
           int *      info);
  void
  psgetri_(const int *N,
           float *    A,
           const int *IA,
           const int *JA,
           const int *DESCA,
           const int *ipiv,
           float *    work,
           int *      lwork,
           int *      iwork,
           int *      liwork,
           int *      info);


  /**
   * PDTRTRI计算上或下三角分布矩阵sub( A ) =
   * A(IA:IA+N-1,JA:JA+N-1)的逆。http://www.netlib.org/scalapack/explore-html/d9/dc0/pdtrtri_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lpdtri.htm
   * https://software.intel.com/en-us/mkl-developer-reference-c-p-trtri
   *
   */
  void
  pdtrtri_(const char *UPLO,
           const char *DIAG,
           const int * N,
           double *    A,
           const int * IA,
           const int * JA,
           const int * DESCA,
           int *       INFO);
  void
  pstrtri_(const char *UPLO,
           const char *DIAG,
           const int * N,
           float *     A,
           const int * IA,
           const int * JA,
           const int * DESCA,
           int *       INFO);

  /**
   * 使用Cholesky分解法估计实数对称正定分布矩阵的条件数的倒数(l1-norm)。
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lpocon.htm#lpocon
   * http://www.netlib.org/scalapack/explore-html/d4/df7/pdpocon_8f.html
   * https://software.intel.com/en-us/mkl-developer-reference-fortran-pocon
   *
   */
  void
  pdpocon_(const char *  uplo,
           const int *   N,
           const double *A,
           const int *   IA,
           const int *   JA,
           const int *   DESCA,
           const double *ANORM,
           double *      RCOND,
           double *      WORK,
           const int *   LWORK,
           int *         IWORK,
           const int *   LIWORK,
           int *         INFO);
  void
  pspocon_(const char * uplo,
           const int *  N,
           const float *A,
           const int *  IA,
           const int *  JA,
           const int *  DESCA,
           const float *ANORM,
           float *      RCOND,
           float *      WORK,
           const int *  LWORK,
           int *        IWORK,
           const int *  LIWORK,
           int *        INFO);

  /**
   * 实对称矩阵的准则
   * http://www.netlib.org/scalapack/explore-html/dd/d12/pdlansy_8f_source.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_pdlansy.htm#pdlansy
   *
   */
  double
  pdlansy_(const char *  norm,
           const char *  uplo,
           const int *   N,
           const double *A,
           const int *   IA,
           const int *   JA,
           const int *   DESCA,
           double *      work);
  float
  pslansy_(const char * norm,
           const char * uplo,
           const int *  N,
           const float *A,
           const int *  IA,
           const int *  JA,
           const int *  DESCA,
           float *      work);

  /**
   * 计算两个正整数 @p M 和 @p N. 的最小公倍数（LCM）
   * 事实上，该例程计算最大公除数（GCD），并利用M*N=GCD*LCM的事实。http://www.netlib.org/scalapack/explore-html/d0/d9b/ilcm_8f_source.html
   *
   */
  int
  ilcm_(const int *M, const int *N);

  /**
   * 返回两个整数的除法的上限。http://www.netlib.org/scalapack/explore-html/df/d07/iceil_8f_source.html
   *
   */
  int
  iceil_(const int *i1, const int *i2);

  /**
   * 用8个输入参数初始化描述符向量
   *
   */
  void
  descinit_(int *      desc,
            const int *m,
            const int *n,
            const int *mb,
            const int *nb,
            const int *irsrc,
            const int *icsrc,
            const int *ictxt,
            const int *lld,
            int *      info);

  /**
   * 计算由 @p iproc. 所示进程的本地索引 @p indxloc
   * 指向的分布式矩阵条目的全局索引  @param  indxloc
   * 分布式矩阵条目的本地索引。    @param  nb
   * 块大小，分布式矩阵被分割成的块的大小。    @param
   * iproc 要确定其本地数组行或列的进程的坐标  @param
   * isrcproc 拥有分布式矩阵第一行/列的进程的坐标  @param
   * nprocs 分布式矩阵的总进程数
   *
   */
  int
  indxl2g_(const int *indxloc,
           const int *nb,
           const int *iproc,
           const int *isrcproc,
           const int *nprocs);

  /**
   * 计算一个实数线性方程组的解
   *
   */
  void
  pdgesv_(const int *n,
          const int *nrhs,
          double *   A,
          const int *ia,
          const int *ja,
          const int *desca,
          int *      ipiv,
          double *   B,
          const int *ib,
          const int *jb,
          const int *descb,
          int *      info);
  void
  psgesv_(const int *n,
          const int *nrhs,
          float *    A,
          const int *ia,
          const int *ja,
          const int *desca,
          int *      ipiv,
          float *    B,
          const int *ib,
          const int *jb,
          const int *descb,
          int *      info);

  /**
   * 执行矩阵-矩阵操作中的一个。
   * @f{align*}{
   * \mathrm{sub}(C) &\dealcoloneq \alpha op(\mathrm{sub}(A))op(\mathrm{sub}(B))
   *                          + \beta \mathrm{sub}(C), \\
   * \mathrm{sub}(C) &\dealcoloneq \alpha op(\mathrm{sub}(A))op(\mathrm{sub}(B))
   *                          + beta sub(C),
   * @f}
   * 其中 $\mathrm{sub}(C)$ 表示C(IC:IC+M-1,JC:JC+N-1)，并且， $op(X)$
   * 是 $op(X) = X$ 或 $op(X) = X^T$ 之一。
   *
   */
  void
  pdgemm_(const char *  transa,
          const char *  transb,
          const int *   m,
          const int *   n,
          const int *   k,
          const double *alpha,
          const double *A,
          const int *   IA,
          const int *   JA,
          const int *   DESCA,
          const double *B,
          const int *   IB,
          const int *   JB,
          const int *   DESCB,
          const double *beta,
          double *      C,
          const int *   IC,
          const int *   JC,
          const int *   DESCC);
  void
  psgemm_(const char * transa,
          const char * transb,
          const int *  m,
          const int *  n,
          const int *  k,
          const float *alpha,
          const float *A,
          const int *  IA,
          const int *  JA,
          const int *  DESCA,
          const float *B,
          const int *  IB,
          const int *  JB,
          const int *  DESCB,
          const float *beta,
          float *      C,
          const int *  IC,
          const int *  JC,
          const int *  DESCC);

  /**
   * 返回一准则的值，或弗罗本纽斯准则，或无穷大准则，或分布式矩阵的最大绝对值元素。
   *
   */
  double
  pdlange_(char const *  norm,
           const int *   m,
           const int *   n,
           const double *A,
           const int *   ia,
           const int *   ja,
           const int *   desca,
           double *      work);
  float
  pslange_(const char * norm,
           const int *  m,
           const int *  n,
           const float *A,
           const int *  ia,
           const int *  ja,
           const int *  desca,
           float *      work);

  /**
   * 计算拥有由全局索引指定的分布式矩阵的条目的过程坐标
   *
   */
  int
  indxg2p_(const int *glob,
           const int *nb,
           const int *iproc,
           const int *isproc,
           const int *nprocs);

  /**
   * 通过调用推荐的ScaLAPACK程序序列，计算实数对称矩阵A的所有特征值和可选的特征向量。在目前的形式下，该例程假定是一个同质系统，并且不检查不同进程中的特征值或特征向量的一致性。正因为如此，异质系统有可能在没有任何错误信息的情况下返回不正确的结果。http://www.netlib.org/scalapack/explore-html/d0/d1a/pdsyev_8f.html
   * https://www.ibm.com/support/knowledgecenter/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lsyev.htm#lsyev
   *
   */
  void
  pdsyev_(const char *jobz,
          const char *uplo,
          const int * m,
          double *    A,
          const int * ia,
          const int * ja,
          int *       desca,
          double *    w,
          double *    z,
          const int * iz,
          const int * jz,
          int *       descz,
          double *    work,
          const int * lwork,
          int *       info);
  void
  pssyev_(const char *jobz,
          const char *uplo,
          const int * m,
          float *     A,
          const int * ia,
          const int * ja,
          int *       desca,
          float *     w,
          float *     z,
          const int * iz,
          const int * jz,
          int *       descz,
          float *     work,
          const int * lwork,
          int *       info);

  /**
   * 将一个分布式矩阵A的全部或部分复制到另一个分布式矩阵B。不进行通信，pdlacpy执行本地复制
   * $\mathrm{sub}(A) \dealcoloneq \mathrm{sub}(B)$  ，其中
   * $\mathrm{sub}(A)$  表示  $A(ia:ia+m-1, ja:ja+n-1)$  ，
   * $\mathrm{sub}(B)$  表示  $B(ib:ib+m-1, jb:jb+n-1)$  。
   *
   */
  void
  pdlacpy_(const char *  uplo,
           const int *   m,
           const int *   n,
           const double *A,
           const int *   ia,
           const int *   ja,
           const int *   desca,
           double *      B,
           const int *   ib,
           const int *   jb,
           const int *   descb);
  void
  pslacpy_(const char * uplo,
           const int *  m,
           const int *  n,
           const float *A,
           const int *  ia,
           const int *  ja,
           const int *  desca,
           float *      B,
           const int *  ib,
           const int *  jb,
           const int *  descb);

  /**
   * 将一般矩形分布式矩阵 @p A
   * 的内容复制到另一个分布式矩阵 @p B
   * 不要求矩阵A和B具有相同的进程网格或块大小，例如，将矩阵从一维进程网格复制到二维进程网格
   * @p ictxt 的上下文至少是上下文A和B中所有进程的联合。
   *
   */
  void
  pdgemr2d_(const int *   m,
            const int *   n,
            const double *A,
            const int *   ia,
            const int *   ja,
            const int *   desca,
            double *      B,
            const int *   ib,
            const int *   jb,
            const int *   descb,
            const int *   ictxt);
  void
  psgemr2d_(const int *  m,
            const int *  n,
            const float *A,
            const int *  ia,
            const int *  ja,
            const int *  desca,
            float *      B,
            const int *  ib,
            const int *  jb,
            const int *  descb,
            const int *  ictxt);

  /**
   * 确定机器精度的辅助程序
   *
   */
  double
  pdlamch_(const int *ictxt, const char *cmach);
  float
  pslamch_(const int *ictxt, const char *cmach);


  /**
   * psyevx计算一个实数对称矩阵A的选定的特征值和可选的特征向量。可以通过指定所需特征值的数值范围或指数范围来选择特征值/向量。
   *
   */
  void
  pdsyevx_(const char *  jobz,
           const char *  range,
           const char *  uplo,
           const int *   n,
           double *      A,
           const int *   ia,
           const int *   ja,
           const int *   desca,
           const double *VL,
           const double *VU,
           const int *   il,
           const int *   iu,
           const double *abstol,
           const int *   m,
           const int *   nz,
           double *      w,
           double *      orfac,
           double *      Z,
           const int *   iz,
           const int *   jz,
           const int *   descz,
           double *      work,
           int *         lwork,
           int *         iwork,
           int *         liwork,
           int *         ifail,
           int *         iclustr,
           double *      gap,
           int *         info);
  void
  pssyevx_(const char * jobz,
           const char * range,
           const char * uplo,
           const int *  n,
           float *      A,
           const int *  ia,
           const int *  ja,
           const int *  desca,
           const float *VL,
           const float *VU,
           const int *  il,
           const int *  iu,
           const float *abstol,
           const int *  m,
           const int *  nz,
           float *      w,
           float *      orfac,
           float *      Z,
           const int *  iz,
           const int *  jz,
           const int *  descz,
           float *      work,
           int *        lwork,
           int *        iwork,
           int *        liwork,
           int *        ifail,
           int *        iclustr,
           float *      gap,
           int *        info);

  /* PDGESVD计算M乘N矩阵A的奇异值分解（SVD），可以选择计算左和/或右奇异向量  
*
*/
  void
  pdgesvd_(const char *jobu,
           const char *jobvt,
           const int * m,
           const int * n,
           double *    A,
           const int * ia,
           const int * ja,
           const int * desca,
           double *    S,
           double *    U,
           const int * iu,
           const int * ju,
           const int * descu,
           double *    VT,
           const int * ivt,
           const int * jvt,
           const int * descvt,
           double *    work,
           int *       lwork,
           int *       info);
  void
  psgesvd_(const char *jobu,
           const char *jobvt,
           const int * m,
           const int * n,
           float *     A,
           const int * ia,
           const int * ja,
           const int * desca,
           float *     S,
           float *     U,
           const int * iu,
           const int * ju,
           const int * descu,
           float *     VT,
           const int * ivt,
           const int * jvt,
           const int * descvt,
           float *     work,
           int *       lwork,
           int *       info);

  /* P_GELS使用A的QR或LQ因子化，解决涉及M乘N矩阵A或其转置的过定或欠定实数线性系统，假定A具有全秩。 
*
*/
  void
  pdgels_(const char *trans,
          const int * m,
          const int * n,
          const int * nrhs,
          double *    A,
          const int * ia,
          const int * ja,
          const int * desca,
          double *    B,
          const int * ib,
          const int * jb,
          const int * descb,
          double *    work,
          int *       lwork,
          int *       info);
  void
  psgels_(const char *trans,
          const int * m,
          const int * n,
          const int * nrhs,
          float *     A,
          const int * ia,
          const int * ja,
          const int * desca,
          float *     B,
          const int * ib,
          const int * jb,
          const int * descb,
          float *     work,
          int *       lwork,
          int *       info);

  /* 进行矩阵求和。  @f{equation*}{
   C \dealcoloneq \beta C + \alpha op(A),
   @f}
* 其中 $op(A)$ 表示 $op(A) = A$ 或 $op(A)=A^T$  。 
*
*/
  void
  pdgeadd_(const char *  transa,
           const int *   m,
           const int *   n,
           const double *alpha,
           const double *A,
           const int *   IA,
           const int *   JA,
           const int *   DESCA,
           const double *beta,
           double *      C,
           const int *   IC,
           const int *   JC,
           const int *   DESCC);
  void
  psgeadd_(const char * transa,
           const int *  m,
           const int *  n,
           const float *alpha,
           const float *A,
           const int *  IA,
           const int *  JA,
           const int *  DESCA,
           const float *beta,
           float *      C,
           const int *  IC,
           const int *  JC,
           const int *  DESCC);

  /**
   * 对一个矩阵进行转置的程序。  C = beta C + alpha A^T
   *
   */
  void
  pdtran_(const int *   m,
          const int *   n,
          const double *alpha,
          const double *A,
          const int *   IA,
          const int *   JA,
          const int *   DESCA,
          const double *beta,
          double *      C,
          const int *   IC,
          const int *   JC,
          const int *   DESCC);
  void
  pstran_(const int *  m,
          const int *  n,
          const float *alpha,
          const float *A,
          const int *  IA,
          const int *  JA,
          const int *  DESCA,
          const float *beta,
          float *      C,
          const int *  IC,
          const int *  JC,
          const int *  DESCC);

  /**
   * psyevr使用MRR算法的并行实现，计算实数对称矩阵A的选定特征值和可选的特征向量。特征值/向量可以通过指定所需特征值的数值范围或指数范围来选择。
   *
   */
  void
  pdsyevr_(const char *  jobz,
           const char *  range,
           const char *  uplo,
           const int *   n,
           double *      A,
           const int *   IA,
           const int *   JA,
           const int *   DESCA,
           const double *VL,
           const double *VU,
           const int *   IL,
           const int *   IU,
           int *         m,
           int *         nz,
           double *      w,
           double *      Z,
           const int *   IZ,
           const int *   JZ,
           const int *   DESCZ,
           double *      work,
           int *         lwork,
           int *         iwork,
           int *         liwork,
           int *         info);
  void
  pssyevr_(const char * jobz,
           const char * range,
           const char * uplo,
           const int *  n,
           float *      A,
           const int *  IA,
           const int *  JA,
           const int *  DESCA,
           const float *VL,
           const float *VU,
           const int *  IL,
           const int *  IU,
           int *        m,
           int *        nz,
           float *      w,
           float *      Z,
           const int *  IZ,
           const int *  JZ,
           const int *  DESCZ,
           float *      work,
           int *        lwork,
           int *        iwork,
           int *        liwork,
           int *        info);
}



/* 在下面的内容中，我们为ScaLAPACK例程提供了模板包装器，其他数字类型的包装器可以在将来添加。

* 
*
*/
template <typename number>
inline void
Cgerv2d(int  /*context*/ ,
        int  /*M*/ ,
        int  /*N*/ ,
        number *  /*A*/ ,
        int  /*lda*/ ,
        int  /*rsrc*/ ,
        int  /*csrc*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
Cgerv2d(int context, int M, int N, double *A, int lda, int rsrc, int csrc)
{
  Cdgerv2d(context, M, N, A, lda, rsrc, csrc);
}

inline void
Cgerv2d(int context, int M, int N, float *A, int lda, int rsrc, int csrc)
{
  Csgerv2d(context, M, N, A, lda, rsrc, csrc);
}


template <typename number>
inline void
Cgesd2d(int  /*context*/ ,
        int  /*M*/ ,
        int  /*N*/ ,
        number *  /*A*/ ,
        int  /*lda*/ ,
        int  /*rdest*/ ,
        int  /*cdest*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
Cgesd2d(int context, int M, int N, double *A, int lda, int rdest, int cdest)
{
  Cdgesd2d(context, M, N, A, lda, rdest, cdest);
}

inline void
Cgesd2d(int context, int M, int N, float *A, int lda, int rdest, int cdest)
{
  Csgesd2d(context, M, N, A, lda, rdest, cdest);
}


template <typename number>
inline void
ppotrf(const char *  /*UPLO*/ ,
       const int *  /*N*/ ,
       number *  /*A*/ ,
       const int *  /*IA*/ ,
       const int *  /*JA*/ ,
       const int *  /*DESCA*/ ,
       int *  /*INFO*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
ppotrf(const char *UPLO,
       const int * N,
       double *    A,
       const int * IA,
       const int * JA,
       const int * DESCA,
       int *       INFO)
{
  pdpotrf_(UPLO, N, A, IA, JA, DESCA, INFO);
}

inline void
ppotrf(const char *UPLO,
       const int * N,
       float *     A,
       const int * IA,
       const int * JA,
       const int * DESCA,
       int *       INFO)
{
  pspotrf_(UPLO, N, A, IA, JA, DESCA, INFO);
}


template <typename number>
inline void
pgetrf(const int *  /*m*/ ,
       const int *  /*n*/ ,
       number *  /*A*/ ,
       const int *  /*IA*/ ,
       const int *  /*JA*/ ,
       const int *  /*DESCA*/ ,
       int *  /*ipiv*/ ,
       int *  /*INFO*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
pgetrf(const int *m,
       const int *n,
       double *   A,
       const int *IA,
       const int *JA,
       const int *DESCA,
       int *      ipiv,
       int *      INFO)
{
  pdgetrf_(m, n, A, IA, JA, DESCA, ipiv, INFO);
}

inline void
pgetrf(const int *m,
       const int *n,
       float *    A,
       const int *IA,
       const int *JA,
       const int *DESCA,
       int *      ipiv,
       int *      INFO)
{
  psgetrf_(m, n, A, IA, JA, DESCA, ipiv, INFO);
}


template <typename number>
inline void
ppotri(const char *  /*UPLO*/ ,
       const int *  /*N*/ ,
       number *  /*A*/ ,
       const int *  /*IA*/ ,
       const int *  /*JA*/ ,
       const int *  /*DESCA*/ ,
       int *  /*INFO*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
ppotri(const char *UPLO,
       const int * N,
       double *    A,
       const int * IA,
       const int * JA,
       const int * DESCA,
       int *       INFO)
{
  pdpotri_(UPLO, N, A, IA, JA, DESCA, INFO);
}

inline void
ppotri(const char *UPLO,
       const int * N,
       float *     A,
       const int * IA,
       const int * JA,
       const int * DESCA,
       int *       INFO)
{
  pspotri_(UPLO, N, A, IA, JA, DESCA, INFO);
}


template <typename number>
inline void
pgetri(const int *  /*N*/ ,
       number *  /*A*/ ,
       const int *  /*IA*/ ,
       const int *  /*JA*/ ,
       const int *  /*DESCA*/ ,
       const int *  /*ipiv*/ ,
       number *  /*work*/ ,
       int *  /*lwork*/ ,
       int *  /*iwork*/ ,
       int *  /*liwork*/ ,
       int *  /*info*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
pgetri(const int *N,
       double *   A,
       const int *IA,
       const int *JA,
       const int *DESCA,
       const int *ipiv,
       double *   work,
       int *      lwork,
       int *      iwork,
       int *      liwork,
       int *      info)
{
  pdgetri_(N, A, IA, JA, DESCA, ipiv, work, lwork, iwork, liwork, info);
}

inline void
pgetri(const int *N,
       float *    A,
       const int *IA,
       const int *JA,
       const int *DESCA,
       const int *ipiv,
       float *    work,
       int *      lwork,
       int *      iwork,
       int *      liwork,
       int *      info)
{
  psgetri_(N, A, IA, JA, DESCA, ipiv, work, lwork, iwork, liwork, info);
}

template <typename number>
inline void
ptrtri(const char *  /*UPLO*/ ,
       const char *  /*DIAG*/ ,
       const int *  /*N*/ ,
       number *  /*A*/ ,
       const int *  /*IA*/ ,
       const int *  /*JA*/ ,
       const int *  /*DESCA*/ ,
       int *  /*INFO*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
ptrtri(const char *UPLO,
       const char *DIAG,
       const int * N,
       double *    A,
       const int * IA,
       const int * JA,
       const int * DESCA,
       int *       INFO)
{
  pdtrtri_(UPLO, DIAG, N, A, IA, JA, DESCA, INFO);
}

inline void
ptrtri(const char *UPLO,
       const char *DIAG,
       const int * N,
       float *     A,
       const int * IA,
       const int * JA,
       const int * DESCA,
       int *       INFO)
{
  pstrtri_(UPLO, DIAG, N, A, IA, JA, DESCA, INFO);
}

template <typename number>
inline void
ppocon(const char *  /*uplo*/ ,
       const int *  /*N*/ ,
       const number *  /*A*/ ,
       const int *  /*IA*/ ,
       const int *  /*JA*/ ,
       const int *  /*DESCA*/ ,
       const number *  /*ANORM*/ ,
       number *  /*RCOND*/ ,
       number *  /*WORK*/ ,
       const int *  /*LWORK*/ ,
       int *  /*IWORK*/ ,
       const int *  /*LIWORK*/ ,
       int *  /*INFO*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
ppocon(const char *  uplo,
       const int *   N,
       const double *A,
       const int *   IA,
       const int *   JA,
       const int *   DESCA,
       const double *ANORM,
       double *      RCOND,
       double *      WORK,
       const int *   LWORK,
       int *         IWORK,
       const int *   LIWORK,
       int *         INFO)
{
  pdpocon_(
    uplo, N, A, IA, JA, DESCA, ANORM, RCOND, WORK, LWORK, IWORK, LIWORK, INFO);
}

inline void
ppocon(const char * uplo,
       const int *  N,
       const float *A,
       const int *  IA,
       const int *  JA,
       const int *  DESCA,
       const float *ANORM,
       float *      RCOND,
       float *      WORK,
       const int *  LWORK,
       int *        IWORK,
       const int *  LIWORK,
       int *        INFO)
{
  pspocon_(
    uplo, N, A, IA, JA, DESCA, ANORM, RCOND, WORK, LWORK, IWORK, LIWORK, INFO);
}


template <typename number>
inline number
plansy(const char *  /*norm*/ ,
       const char *  /*uplo*/ ,
       const int *  /*N*/ ,
       const number *  /*A*/ ,
       const int *  /*IA*/ ,
       const int *  /*JA*/ ,
       const int *  /*DESCA*/ ,
       number *  /*work*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline double
plansy(const char *  norm,
       const char *  uplo,
       const int *   N,
       const double *A,
       const int *   IA,
       const int *   JA,
       const int *   DESCA,
       double *      work)
{
  return pdlansy_(norm, uplo, N, A, IA, JA, DESCA, work);
}

inline float
plansy(const char * norm,
       const char * uplo,
       const int *  N,
       const float *A,
       const int *  IA,
       const int *  JA,
       const int *  DESCA,
       float *      work)
{
  return pslansy_(norm, uplo, N, A, IA, JA, DESCA, work);
}


template <typename number>
inline void
pgesv(const int *  /*n*/ ,
      const int *  /*nrhs*/ ,
      number *  /*A*/ ,
      const int *  /*ia*/ ,
      const int *  /*ja*/ ,
      const int *  /*desca*/ ,
      int *  /*ipiv*/ ,
      number *  /*B*/ ,
      const int *  /*ib*/ ,
      const int *  /*jb*/ ,
      const int *  /*descb*/ ,
      int *  /*info*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
pgesv(const int *n,
      const int *nrhs,
      double *   A,
      const int *ia,
      const int *ja,
      const int *desca,
      int *      ipiv,
      double *   B,
      const int *ib,
      const int *jb,
      const int *descb,
      int *      info)
{
  pdgesv_(n, nrhs, A, ia, ja, desca, ipiv, B, ib, jb, descb, info);
}

inline void
pgesv(const int *n,
      const int *nrhs,
      float *    A,
      const int *ia,
      const int *ja,
      const int *desca,
      int *      ipiv,
      float *    B,
      const int *ib,
      const int *jb,
      const int *descb,
      int *      info)
{
  psgesv_(n, nrhs, A, ia, ja, desca, ipiv, B, ib, jb, descb, info);
}


template <typename number>
inline void
pgemm(const char *  /*transa*/ ,
      const char *  /*transb*/ ,
      const int *  /*m*/ ,
      const int *  /*n*/ ,
      const int *  /*k*/ ,
      const number *  /*alpha*/ ,
      number *  /*A*/ ,
      const int *  /*IA*/ ,
      const int *  /*JA*/ ,
      const int *  /*DESCA*/ ,
      number *  /*B*/ ,
      const int *  /*IB*/ ,
      const int *  /*JB*/ ,
      const int *  /*DESCB*/ ,
      const number *  /*beta*/ ,
      number *  /*C*/ ,
      const int *  /*IC*/ ,
      const int *  /*JC*/ ,
      const int *  /*DESCC*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
pgemm(const char *  transa,
      const char *  transb,
      const int *   m,
      const int *   n,
      const int *   k,
      const double *alpha,
      const double *A,
      const int *   IA,
      const int *   JA,
      const int *   DESCA,
      const double *B,
      const int *   IB,
      const int *   JB,
      const int *   DESCB,
      const double *beta,
      double *      C,
      const int *   IC,
      const int *   JC,
      const int *   DESCC)
{
  pdgemm_(transa,
          transb,
          m,
          n,
          k,
          alpha,
          A,
          IA,
          JA,
          DESCA,
          B,
          IB,
          JB,
          DESCB,
          beta,
          C,
          IC,
          JC,
          DESCC);
}

inline void
pgemm(const char * transa,
      const char * transb,
      const int *  m,
      const int *  n,
      const int *  k,
      const float *alpha,
      const float *A,
      const int *  IA,
      const int *  JA,
      const int *  DESCA,
      const float *B,
      const int *  IB,
      const int *  JB,
      const int *  DESCB,
      const float *beta,
      float *      C,
      const int *  IC,
      const int *  JC,
      const int *  DESCC)
{
  psgemm_(transa,
          transb,
          m,
          n,
          k,
          alpha,
          A,
          IA,
          JA,
          DESCA,
          B,
          IB,
          JB,
          DESCB,
          beta,
          C,
          IC,
          JC,
          DESCC);
}


template <typename number>
inline number
plange(const char *  /*norm*/ ,
       const int *  /*m*/ ,
       const int *  /*n*/ ,
       const number *  /*A*/ ,
       const int *  /*ia*/ ,
       const int *  /*ja*/ ,
       const int *  /*desca*/ ,
       number *  /*work*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline double
plange(const char *  norm,
       const int *   m,
       const int *   n,
       const double *A,
       const int *   ia,
       const int *   ja,
       const int *   desca,
       double *      work)
{
  return pdlange_(norm, m, n, A, ia, ja, desca, work);
}

inline float
plange(const char * norm,
       const int *  m,
       const int *  n,
       const float *A,
       const int *  ia,
       const int *  ja,
       const int *  desca,
       float *      work)
{
  return pslange_(norm, m, n, A, ia, ja, desca, work);
}


template <typename number>
inline void
psyev(const char *  /*jobz*/ ,
      const char *  /*uplo*/ ,
      const int *  /*m*/ ,
      number *  /*A*/ ,
      const int *  /*ia*/ ,
      const int *  /*ja*/ ,
      int *  /*desca*/ ,
      number *  /*w*/ ,
      number *  /*z*/ ,
      const int *  /*iz*/ ,
      const int *  /*jz*/ ,
      int *  /*descz*/ ,
      number *  /*work*/ ,
      const int *  /*lwork*/ ,
      int *  /*info*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
psyev(const char *jobz,
      const char *uplo,
      const int * m,
      double *    A,
      const int * ia,
      const int * ja,
      int *       desca,
      double *    w,
      double *    z,
      const int * iz,
      const int * jz,
      int *       descz,
      double *    work,
      const int * lwork,
      int *       info)
{
  pdsyev_(
    jobz, uplo, m, A, ia, ja, desca, w, z, iz, jz, descz, work, lwork, info);
}

inline void
psyev(const char *jobz,
      const char *uplo,
      const int * m,
      float *     A,
      const int * ia,
      const int * ja,
      int *       desca,
      float *     w,
      float *     z,
      const int * iz,
      const int * jz,
      int *       descz,
      float *     work,
      const int * lwork,
      int *       info)
{
  pssyev_(
    jobz, uplo, m, A, ia, ja, desca, w, z, iz, jz, descz, work, lwork, info);
}


template <typename number>
inline void
placpy(const char *  /*uplo*/ ,
       const int *  /*m*/ ,
       const int *  /*n*/ ,
       const number *  /*A*/ ,
       const int *  /*ia*/ ,
       const int *  /*ja*/ ,
       const int *  /*desca*/ ,
       number *  /*B*/ ,
       const int *  /*ib*/ ,
       const int *  /*jb*/ ,
       const int *  /*descb*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
placpy(const char *  uplo,
       const int *   m,
       const int *   n,
       const double *A,
       const int *   ia,
       const int *   ja,
       const int *   desca,
       double *      B,
       const int *   ib,
       const int *   jb,
       const int *   descb)
{
  pdlacpy_(uplo, m, n, A, ia, ja, desca, B, ib, jb, descb);
}

inline void
placpy(const char * uplo,
       const int *  m,
       const int *  n,
       const float *A,
       const int *  ia,
       const int *  ja,
       const int *  desca,
       float *      B,
       const int *  ib,
       const int *  jb,
       const int *  descb)
{
  pslacpy_(uplo, m, n, A, ia, ja, desca, B, ib, jb, descb);
}


template <typename number>
inline void
pgemr2d(const int *  /*m*/ ,
        const int *  /*n*/ ,
        const number *  /*A*/ ,
        const int *  /*ia*/ ,
        const int *  /*ja*/ ,
        const int *  /*desca*/ ,
        number *  /*B*/ ,
        const int *  /*ib*/ ,
        const int *  /*jb*/ ,
        const int *  /*descb*/ ,
        const int *  /*ictxt*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
pgemr2d(const int *   m,
        const int *   n,
        const double *A,
        const int *   ia,
        const int *   ja,
        const int *   desca,
        double *      B,
        const int *   ib,
        const int *   jb,
        const int *   descb,
        const int *   ictxt)
{
  pdgemr2d_(m, n, A, ia, ja, desca, B, ib, jb, descb, ictxt);
}

inline void
pgemr2d(const int *  m,
        const int *  n,
        const float *A,
        const int *  ia,
        const int *  ja,
        const int *  desca,
        float *      B,
        const int *  ib,
        const int *  jb,
        const int *  descb,
        const int *  ictxt)
{
  psgemr2d_(m, n, A, ia, ja, desca, B, ib, jb, descb, ictxt);
}


template <typename number>
inline void
plamch(const int *  /*ictxt*/, const char * /*cmach*/, number & /*val*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
plamch(const int *ictxt, const char *cmach, double &val)
{
  val = pdlamch_(ictxt, cmach);
}

inline void
plamch(const int *ictxt, const char *cmach, float &val)
{
  val = pslamch_(ictxt, cmach);
}


template <typename number>
inline void
psyevx(const char *  /*jobz*/ ,
       const char *  /*range*/ ,
       const char *  /*uplo*/ ,
       const int *  /*n*/ ,
       number *  /*A*/ ,
       const int *  /*ia*/ ,
       const int *  /*ja*/ ,
       const int *  /*desca*/ ,
       number *  /*VL*/ ,
       number *  /*VU*/ ,
       const int *  /*il*/ ,
       const int *  /*iu*/ ,
       number *  /*abstol*/ ,
       const int *  /*m*/ ,
       const int *  /*nz*/ ,
       number *  /*w*/ ,
       number *  /*orfac*/ ,
       number *  /*Z*/ ,
       const int *  /*iz*/ ,
       const int *  /*jz*/ ,
       const int *  /*descz*/ ,
       number *  /*work*/ ,
       int *  /*lwork*/ ,
       int *  /*iwork*/ ,
       int *  /*liwork*/ ,
       int *  /*ifail*/ ,
       int *  /*iclustr*/ ,
       number *  /*gap*/ ,
       int *  /*info*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
psyevx(const char *jobz,
       const char *range,
       const char *uplo,
       const int * n,
       double *    A,
       const int * ia,
       const int * ja,
       const int * desca,
       double *    VL,
       double *    VU,
       const int * il,
       const int * iu,
       double *    abstol,
       const int * m,
       const int * nz,
       double *    w,
       double *    orfac,
       double *    Z,
       const int * iz,
       const int * jz,
       const int * descz,
       double *    work,
       int *       lwork,
       int *       iwork,
       int *       liwork,
       int *       ifail,
       int *       iclustr,
       double *    gap,
       int *       info)
{
  pdsyevx_(jobz,
           range,
           uplo,
           n,
           A,
           ia,
           ja,
           desca,
           VL,
           VU,
           il,
           iu,
           abstol,
           m,
           nz,
           w,
           orfac,
           Z,
           iz,
           jz,
           descz,
           work,
           lwork,
           iwork,
           liwork,
           ifail,
           iclustr,
           gap,
           info);
}

inline void
psyevx(const char *jobz,
       const char *range,
       const char *uplo,
       const int * n,
       float *     A,
       const int * ia,
       const int * ja,
       const int * desca,
       float *     VL,
       float *     VU,
       const int * il,
       const int * iu,
       float *     abstol,
       const int * m,
       const int * nz,
       float *     w,
       float *     orfac,
       float *     Z,
       const int * iz,
       const int * jz,
       const int * descz,
       float *     work,
       int *       lwork,
       int *       iwork,
       int *       liwork,
       int *       ifail,
       int *       iclustr,
       float *     gap,
       int *       info)
{
  pssyevx_(jobz,
           range,
           uplo,
           n,
           A,
           ia,
           ja,
           desca,
           VL,
           VU,
           il,
           iu,
           abstol,
           m,
           nz,
           w,
           orfac,
           Z,
           iz,
           jz,
           descz,
           work,
           lwork,
           iwork,
           liwork,
           ifail,
           iclustr,
           gap,
           info);
}


template <typename number>
inline void
pgesvd(const char *  /*jobu*/ ,
       const char *  /*jobvt*/ ,
       const int *  /*m*/ ,
       const int *  /*n*/ ,
       number *  /*A*/ ,
       const int *  /*ia*/ ,
       const int *  /*ja*/ ,
       const int *  /*desca*/ ,
       number *  /*S*/ ,
       number *  /*U*/ ,
       const int *  /*iu*/ ,
       const int *  /*ju*/ ,
       const int *  /*descu*/ ,
       number *  /*VT*/ ,
       const int *  /*ivt*/ ,
       const int *  /*jvt*/ ,
       const int *  /*descvt*/ ,
       number *  /*work*/ ,
       int *  /*lwork*/ ,
       int *  /*info*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
pgesvd(const char *jobu,
       const char *jobvt,
       const int * m,
       const int * n,
       double *    A,
       const int * ia,
       const int * ja,
       const int * desca,
       double *    S,
       double *    U,
       const int * iu,
       const int * ju,
       const int * descu,
       double *    VT,
       const int * ivt,
       const int * jvt,
       const int * descvt,
       double *    work,
       int *       lwork,
       int *       info)
{
  pdgesvd_(jobu,
           jobvt,
           m,
           n,
           A,
           ia,
           ja,
           desca,
           S,
           U,
           iu,
           ju,
           descu,
           VT,
           ivt,
           jvt,
           descvt,
           work,
           lwork,
           info);
}

inline void
pgesvd(const char *jobu,
       const char *jobvt,
       const int * m,
       const int * n,
       float *     A,
       const int * ia,
       const int * ja,
       const int * desca,
       float *     S,
       float *     U,
       const int * iu,
       const int * ju,
       const int * descu,
       float *     VT,
       const int * ivt,
       const int * jvt,
       const int * descvt,
       float *     work,
       int *       lwork,
       int *       info)
{
  psgesvd_(jobu,
           jobvt,
           m,
           n,
           A,
           ia,
           ja,
           desca,
           S,
           U,
           iu,
           ju,
           descu,
           VT,
           ivt,
           jvt,
           descvt,
           work,
           lwork,
           info);
}


template <typename number>
inline void
pgels(const char *  /*trans*/ ,
      const int *  /*m*/ ,
      const int *  /*n*/ ,
      const int *  /*nrhs*/ ,
      number *  /*A*/ ,
      const int *  /*ia*/ ,
      const int *  /*ja*/ ,
      const int *  /*desca*/ ,
      number *  /*B*/ ,
      const int *  /*ib*/ ,
      const int *  /*jb*/ ,
      const int *  /*descb*/ ,
      number *  /*work*/ ,
      int *  /*lwork*/ ,
      int *  /*info*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
pgels(const char *trans,
      const int * m,
      const int * n,
      const int * nrhs,
      double *    A,
      const int * ia,
      const int * ja,
      const int * desca,
      double *    B,
      const int * ib,
      const int * jb,
      const int * descb,
      double *    work,
      int *       lwork,
      int *       info)
{
  pdgels_(
    trans, m, n, nrhs, A, ia, ja, desca, B, ib, jb, descb, work, lwork, info);
}

inline void
pgels(const char *trans,
      const int * m,
      const int * n,
      const int * nrhs,
      float *     A,
      const int * ia,
      const int * ja,
      const int * desca,
      float *     B,
      const int * ib,
      const int * jb,
      const int * descb,
      float *     work,
      int *       lwork,
      int *       info)
{
  psgels_(
    trans, m, n, nrhs, A, ia, ja, desca, B, ib, jb, descb, work, lwork, info);
}


template <typename number>
inline void
pgeadd(const char *  /*transa*/ ,
       const int *  /*m*/ ,
       const int *  /*n*/ ,
       const number *  /*alpha*/ ,
       const number *  /*A*/ ,
       const int *  /*IA*/ ,
       const int *  /*JA*/ ,
       const int *  /*DESCA*/ ,
       const number *  /*beta*/ ,
       number *  /*C*/ ,
       const int *  /*IC*/ ,
       const int *  /*JC*/ ,
       const int *  /*DESCC*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
pgeadd(const char *  transa,
       const int *   m,
       const int *   n,
       const double *alpha,
       const double *A,
       const int *   IA,
       const int *   JA,
       const int *   DESCA,
       const double *beta,
       double *      C,
       const int *   IC,
       const int *   JC,
       const int *   DESCC)
{
  pdgeadd_(transa, m, n, alpha, A, IA, JA, DESCA, beta, C, IC, JC, DESCC);
}

inline void
pgeadd(const char * transa,
       const int *  m,
       const int *  n,
       const float *alpha,
       const float *A,
       const int *  IA,
       const int *  JA,
       const int *  DESCA,
       const float *beta,
       float *      C,
       const int *  IC,
       const int *  JC,
       const int *  DESCC)
{
  psgeadd_(transa, m, n, alpha, A, IA, JA, DESCA, beta, C, IC, JC, DESCC);
}


template <typename number>
inline void
ptran(const int *  /*m*/ ,
      const int *  /*n*/ ,
      const number *  /*alpha*/ ,
      const number *  /*A*/ ,
      const int *  /*IA*/ ,
      const int *  /*JA*/ ,
      const int *  /*DESCA*/ ,
      const number *  /*beta*/ ,
      number *  /*C*/ ,
      const int *  /*IC*/ ,
      const int *  /*JC*/ ,
      const int *  /*DESCC*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
ptran(const int *   m,
      const int *   n,
      const double *alpha,
      const double *A,
      const int *   IA,
      const int *   JA,
      const int *   DESCA,
      const double *beta,
      double *      C,
      const int *   IC,
      const int *   JC,
      const int *   DESCC)
{
  pdtran_(m, n, alpha, A, IA, JA, DESCA, beta, C, IC, JC, DESCC);
}

inline void
ptran(const int *  m,
      const int *  n,
      const float *alpha,
      const float *A,
      const int *  IA,
      const int *  JA,
      const int *  DESCA,
      const float *beta,
      float *      C,
      const int *  IC,
      const int *  JC,
      const int *  DESCC)
{
  pstran_(m, n, alpha, A, IA, JA, DESCA, beta, C, IC, JC, DESCC);
}


template <typename number>
inline void
psyevr(const char *  /*jobz*/ ,
       const char *  /*range*/ ,
       const char *  /*uplo*/ ,
       const int *  /*n*/ ,
       number *  /*A*/ ,
       const int *  /*IA*/ ,
       const int *  /*JA*/ ,
       const int *  /*DESCA*/ ,
       const number *  /*VL*/ ,
       const number *  /*VU*/ ,
       const int *  /*IL*/ ,
       const int *  /*IU*/ ,
       int *  /*m*/ ,
       int *  /*nz*/ ,
       number *  /*w*/ ,
       number *  /*Z*/ ,
       const int *  /*IZ*/ ,
       const int *  /*JZ*/ ,
       const int *  /*DESCZ*/ ,
       number *  /*work*/ ,
       int *  /*lwork*/ ,
       int *  /*iwork*/ ,
       int *  /*liwork*/ ,
       int *  /*info*/ )
{
  Assert(false, dealii::ExcNotImplemented());
}

inline void
psyevr(const char *  jobz,
       const char *  range,
       const char *  uplo,
       const int *   n,
       double *      A,
       const int *   IA,
       const int *   JA,
       const int *   DESCA,
       const double *VL,
       const double *VU,
       const int *   IL,
       const int *   IU,
       int *         m,
       int *         nz,
       double *      w,
       double *      Z,
       const int *   IZ,
       const int *   JZ,
       const int *   DESCZ,
       double *      work,
       int *         lwork,
       int *         iwork,
       int *         liwork,
       int *         info)
{
  /* Netlib ScaLAPACK在调用pdsyevr时进行浮点测试（例如除以零），导致浮点异常被抛出（至少在调试模式下）。因此，我们把对pdsyevr的调用包装成以下代码，以抑制异常。 
*
*/
#  ifdef DEAL_II_HAVE_FP_EXCEPTIONS
  fenv_t fp_exceptions;
  feholdexcept(&fp_exceptions);
#  endif

  pdsyevr_(jobz,
           range,
           uplo,
           n,
           A,
           IA,
           JA,
           DESCA,
           VL,
           VU,
           IL,
           IU,
           m,
           nz,
           w,
           Z,
           IZ,
           JZ,
           DESCZ,
           work,
           lwork,
           iwork,
           liwork,
           info);

#  ifdef DEAL_II_HAVE_FP_EXCEPTIONS
  fesetenv(&fp_exceptions);
#  endif
}

inline void
psyevr(const char * jobz,
       const char * range,
       const char * uplo,
       const int *  n,
       float *      A,
       const int *  IA,
       const int *  JA,
       const int *  DESCA,
       const float *VL,
       const float *VU,
       const int *  IL,
       const int *  IU,
       int *        m,
       int *        nz,
       float *      w,
       float *      Z,
       const int *  IZ,
       const int *  JZ,
       const int *  DESCZ,
       float *      work,
       int *        lwork,
       int *        iwork,
       int *        liwork,
       int *        info)
{
  /* Netlib ScaLAPACK在调用pssyevr时进行浮点测试（例如除以零），导致浮点异常被抛出（至少在调试模式下）。因此，我们把对pssyevr的调用包装成以下代码，以抑制异常。 
*
*/
#  ifdef DEAL_II_HAVE_FP_EXCEPTIONS
  fenv_t fp_exceptions;
  feholdexcept(&fp_exceptions);
#  endif

  pssyevr_(jobz,
           range,
           uplo,
           n,
           A,
           IA,
           JA,
           DESCA,
           VL,
           VU,
           IL,
           IU,
           m,
           nz,
           w,
           Z,
           IZ,
           JZ,
           DESCZ,
           work,
           lwork,
           iwork,
           liwork,
           info);

#  ifdef DEAL_II_HAVE_FP_EXCEPTIONS
  fesetenv(&fp_exceptions);
#  endif
}

#endif // DEAL_II_WITH_SCALAPACK

#endif // dealii_scalapack_templates_h


