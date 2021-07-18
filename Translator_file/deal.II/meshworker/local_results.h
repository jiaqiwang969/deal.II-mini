//include/deal.II-translator/meshworker/local_results_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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


#ifndef dealii_mesh_worker_local_results_h
#define dealii_mesh_worker_local_results_h

#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/matrix_block.h>

#include <deal.II/meshworker/vector_selector.h>

#include <functional>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
class BlockIndices;
#endif

/**
 * 一个用于网格循环的函数和类的集合，它是每个有限元程序中无处不在的部分。
 * 这个命名空间的主力是loop()函数，它实现了对所有网格单元的完全通用循环。由于对loop()的调用由于其通用性而容易出错，对于许多应用来说，最好从
 * MeshWorker::LocalIntegrator 派生出一个类，并使用不太通用的
 * integration_loop()。
 * loop()依赖于作为参数交给它的某些对象。这些对象有两种类型，
 * @p info 对象，如 DoFInfo 和 IntegrationInfo，以及 LocalWorker 和
 * IntegrationWorker 等工作对象。
 * 工作者对象通常做两个不同的工作：首先，他们计算一个单元或面对全局操作的本地贡献。其次，它们将这个局部贡献集合到全局结果中，无论是函数、形式还是双线性形式。第一项工作是针对被解决的问题的，而第二项工作是通用的，只取决于数据结构。因此，工人组装成全局数据的基类在命名空间Assembler中提供。
 * <h3>Template argument types</h3>
 * 函数loop()和cell_action()需要一些参数，这些参数是模板参数。让我们在这里列出这些类的最低要求并描述它们的属性。
 * <h4>ITERATOR</h4>
 * 任何具有<tt>operator++()</tt>并指向TriaAccessor或派生类的对象。
 * <h4>DOFINFO</h4>
 * 关于一个实现的例子，请参考类模板DoFInfo。为了与cell_action()和loop()协同工作，DOFINFO需要遵循以下接口。
 *
 * @code
 * class DOFINFO
 * {
 * private:
 *   DOFINFO();
 *   DOFINFO(const DOFINFO&);
 *   DOFINFO& operator=(const DOFINFO&);
 *
 * public:
 *   template <class CellIt>
 *   void reinit(const CellIt& c);
 *
 *   template <class CellIt, class FaceIt>
 *   void reinit(const CellIt& c, const FaceIt& f, const unsigned int n);
 *
 *   template <class CellIt, class FaceIt>
 *   void reinit(const CellIt& c, const FaceIt& f, const unsigned int n,
 *               const unsigned int s);
 *
 * friend template class DoFInfoBox<int dim, DOFINFO>;
 * };
 * @endcode
 *
 * 这三个私有函数是由DoFInfoBox调用的，其他地方应该不需要。很明显，它们可以被公开，然后末尾的friend声明可能会被遗漏。
 * 此外，你将需要至少一个公共构造函数。此外，DOFINFO还相当无用：需要与INTEGRATIONINFO和ASSEMBLER接口的函数。
 * DOFINFO对象被聚集在一个DoFInfoBox中。在这些对象中，我们存储了对每个单元及其面的局部操作的结果。一旦所有这些信息都被收集起来，一个ASSEMBLER被用来将其组装成全局数据。
 * <h4>INFOBOX</h4>
 * 这种类型在IntegrationInfoBox中得到了体现。它在INFO对象中收集单元格和面的动作的输入数据（见下文）。它为loop()和cell_action()提供以下接口。
 *
 *
 * @code
 * class INFOBOX
 * {
 * public:
 *   template <int dim, class DOFINFO>
 *   void post_cell(const DoFInfoBox<dim, DOFINFO>&);
 *
 *   template <int dim, class DOFINFO>
 *   void post_faces(const DoFInfoBox<dim, DOFINFO>&);
 *
 *   INFO cell;
 *   INFO boundary;
 *   INFO face;
 *   INFO subface;
 *   INFO neighbor;
 * };
 * @endcode
 *
 * 这个类的主要目的是收集五个INFO对象，它们包含了每个单元或面所使用的临时数据。对这些对象的要求列在下面。在这里，我们只注意到需要有这5个对象的名字在上面列出。
 * 两个函数模板是在cell_action()中调用的回调函数。第一个是在工作面之前调用的，第二个是在工作面之后调用的。
 * <h4>INFO</h4>
 * 关于这些对象的例子，请参见IntegrationInfo。它们包含了每个单元或面所需的临时数据，以计算出结果。MeshWorker只使用接口
 *
 * @code
 * class INFO
 * {
 * public:
 *   void reinit(const DOFINFO& i);
 * };
 * @endcode
 *
 * <h3>Simplified interfaces</h3>
 * 由于loop()是相当普遍的，所以有一个专门的integration_loop()，它是一个具有简化接口的loop()的包装器。
 * integration_loop()函数循环从IntegrationInfoBox对象中获取大部分需要传递给loop()的信息。它的用途在
 * step-12
 * 中作了解释，但简而言之，它需要在单元格、内部或边界面上做局部积分的函数，它还需要一个对象（称为
 * "装配器"）将这些局部贡献复制到全局矩阵和右手对象中。
 * 在我们能够运行积分循环之前，我们必须在我们的IntegrationWorker和assembler对象中初始化几个数据结构。例如，我们必须决定正交规则，或者我们可能需要比默认更新标志更多的东西。
 *
 *
 * @ingroup MeshWorker
 *
 * @ingroup Integrators
 *
 */
namespace MeshWorker
{
  /**
   * 提供剪贴簿的类，以填充局部整合的结果。根据网格工作者循环所执行的任务，局部结果可以是不同的类型。它们可以是标量，也可以是等于积分中使用的自由度数的向量，或者是同样大小的方阵。所有这些都有一个共同点，那就是它们是在一个单元或面上进行局部积分的结果。哪种对象是一个操作的结果，是由使用它们的汇编器决定的。也是装配者决定<i>how many</i>每种对象的产生（例如，一个装配者可以同时创建质量和刚度矩阵的局部贡献），以及将局部结果的数组设置为所需的大小。    该类的接口允许通过以下函数访问所有这些信息。      <ol>   <li>  标量：n_values()返回该类对象存储的标量数量，通过value()函数访问它们。      <li>  向量：n_vectors()返回该类对象存储的向量数量（每个向量的长度等于该单元上发生积分的自由度数量）。  向量是由vector()函数访问的。      <li>  矩阵：n_matrices()返回存储的矩阵数量，每个矩阵的维数等于每个单元的自由度数。矩阵由matrix()访问，第二个参数是<tt>false</tt>。这些矩阵在同一个单元中耦合自由度。对于跨面的通量，还有一组同样大小的矩阵，这些矩阵的维度与两个单元的自由度有关。这些矩阵可以通过matrix()访问，使用第二个参数<tt>true</tt>。    </ol>  本地矩阵由 @p info 对象的reinit()初始化，然后由Assembler类组装成全局系统。
   * @ingroup MeshWorker
   *
   */
  template <typename number>
  class LocalResults
  {
  public:
    /**
     * 当前对象所存储的标量值的数量。        这个数字被
     * Assembler::CellsAndFaces 设置为非零值。
     *
     */
    unsigned int
    n_values() const;

    /**
     * 当前对象所存储的向量的数量。        这个数字被
     * Assembler::ResidualSimple 和
     * Assembler::ResidualLocalBlocksToGlobalBlocks. 设置为非零值。
     *
     */
    unsigned int
    n_vectors() const;

    /**
     * 当前对象存储的矩阵数量。
     *
     */
    unsigned int
    n_matrices() const;

    /**
     * quadrature_values()中正交点的数量。
     *
     */
    unsigned int
    n_quadrature_points() const;

    /**
     * quadrature_values()中每个正交点的值的数量。
     *
     */
    unsigned int
    n_quadrature_values() const;

    /**
     * 对该类存储的第i个标量的读写访问。
     *
     */
    number &
    value(const unsigned int i);

    /**
     * 对该类存储的第`i`个标量的读取访问。
     *
     */
    number
    value(const unsigned int i) const;

    /**
     * 对该类所存储的第i个向量的读写访问。
     *
     */
    BlockVector<number> &
    vector(const unsigned int i);

    /**
     * 对该类存储的第`i`个向量的读写权限
     *
     */
    const BlockVector<number> &
    vector(const unsigned int i) const;

    /**
     * 对该类存储的第`i`个矩阵的读写访问。
     * 关于第二个参数的解释，请看当前类本身的文档。
     *
     */
    MatrixBlock<FullMatrix<number>> &
    matrix(const unsigned int i, const bool external = false);

    /**
     * 读取该类所存储的第`i`个矩阵的访问权限。
     * 关于第二个参数的解释，请看当前类本身的文档。
     *
     */
    const MatrixBlock<FullMatrix<number>> &
    matrix(const unsigned int i, const bool external = false) const;

    /**
     * 访问正交点数据的向量#quadrature_data，其组织方式是每个点都有一个向量，每个分量都包含一个条目。
     *
     */
    Table<2, number> &
    quadrature_values();

    /**
     * 访问正交点<i>i</i>的第<i>k</i>个值。
     *
     */
    number &
    quadrature_value(const unsigned int k, const unsigned int i);

    /**
     * 读取正交点<i>k</i>的<i>i</i>的值。
     *
     */
    number
    quadrature_value(const unsigned int k, const unsigned int i) const;

    /**
     * 用标量值初始化向量。
     * @note 这个函数通常只由汇编程序调用。
     *
     */
    void
    initialize_numbers(const unsigned int n);

    /**
     * 用向量值初始化向量。
     * @note  这个函数通常只由汇编器调用。
     *
     */
    void
    initialize_vectors(const unsigned int n);

    /**
     * 分配 @p n
     * 本地矩阵。此外，将它们的块行和列坐标设置为零。矩阵本身的大小由reinit()来调整。
     * @note  这个函数通常只由汇编器调用。
     *
     */
    void
    initialize_matrices(const unsigned int n, bool both);

    /**
     * 为 @p matrices.
     * 中的每个全局矩阵分配一个本地矩阵，另外，设置它们的块行和列坐标。矩阵本身的大小由reinit()来调整。
     * @note 这个函数通常只由汇编器调用。
     *
     */
    template <typename MatrixType>
    void
    initialize_matrices(const MatrixBlockVector<MatrixType> &matrices,
                        bool                                 both);

    /**
     * 为 @p
     * 矩阵中的每个全局级对象分配一个本地矩阵。另外，设置它们的块行和块列坐标。矩阵本身的大小由reinit()来调整。
     * @note 这个函数通常只由汇编器调用。
     *
     */
    template <typename MatrixType>
    void
    initialize_matrices(const MGMatrixBlockVector<MatrixType> &matrices,
                        bool                                   both);

    /**
     * 将正交值初始化为<tt>nv</tt>值在<tt>np</tt>正交点。
     *
     */
    void
    initialize_quadrature(const unsigned int np, const unsigned int nv);

    /**
     * 为新单元重新初始化矩阵。不调整存储在此对象中的任何数据向量的大小，但调整#R中的向量以及#M1和#M2中的矩阵的大小为hp，并将它们设置为零。
     *
     */
    void
    reinit(const BlockIndices &local_sizes);

    template <class StreamType>
    void
    print_debug(StreamType &os) const;

    /**
     * 这个对象所使用的内存。
     *
     */
    std::size_t
    memory_consumption() const;

  private:
    /**
     * 本地的数字，在一个单元或一个面上计算。
     *
     */
    std::vector<number> J;

    /**
     * 本地向量。这个字段是公开的，所以本地积分器可以写到它。
     *
     */
    std::vector<BlockVector<number>> R;

    /**
     * 耦合单元本身或面的第一个单元的自由度的局部矩阵。
     *
     */
    std::vector<MatrixBlock<FullMatrix<number>>> M1;

    /**
     * 耦合单元上的测试函数和其他单元上的试验函数的局部矩阵。
     * 只在内部面使用。
     *
     */
    std::vector<MatrixBlock<FullMatrix<number>>> M2;

    /**
     * 用于写入补丁数据的正交点的值。
     *
     */
    Table<2, number> quadrature_data;
  };

  //----------------------------------------------------------------------//

  template <typename number>
  inline void
  LocalResults<number>::initialize_numbers(const unsigned int n)
  {
    J.resize(n);
  }


  template <typename number>
  inline void
  LocalResults<number>::initialize_vectors(const unsigned int n)
  {
    R.resize(n);
  }


  template <typename number>
  template <typename MatrixType>
  inline void
  LocalResults<number>::initialize_matrices(
    const MatrixBlockVector<MatrixType> &matrices,
    bool                                 both)
  {
    M1.resize(matrices.size());
    if (both)
      M2.resize(matrices.size());
    for (unsigned int i = 0; i < matrices.size(); ++i)
      {
        const unsigned int row = matrices.block(i).row;
        const unsigned int col = matrices.block(i).column;

        M1[i].row    = row;
        M1[i].column = col;
        if (both)
          {
            M2[i].row    = row;
            M2[i].column = col;
          }
      }
  }


  template <typename number>
  template <typename MatrixType>
  inline void
  LocalResults<number>::initialize_matrices(
    const MGMatrixBlockVector<MatrixType> &matrices,
    bool                                   both)
  {
    M1.resize(matrices.size());
    if (both)
      M2.resize(matrices.size());
    for (unsigned int i = 0; i < matrices.size(); ++i)
      {
        const MGLevelObject<MatrixBlock<MatrixType>> &o = matrices.block(i);
        const unsigned int row                          = o[o.min_level()].row;
        const unsigned int col = o[o.min_level()].column;

        M1[i].row    = row;
        M1[i].column = col;
        if (both)
          {
            M2[i].row    = row;
            M2[i].column = col;
          }
      }
  }


  template <typename number>
  inline void
  LocalResults<number>::initialize_matrices(const unsigned int n,
                                            const bool         both)
  {
    M1.resize(n);
    if (both)
      M2.resize(n);
    for (unsigned int i = 0; i < n; ++i)
      {
        M1[i].row    = 0;
        M1[i].column = 0;
        if (both)
          {
            M2[i].row    = 0;
            M2[i].column = 0;
          }
      }
  }


  template <typename number>
  inline void
  LocalResults<number>::initialize_quadrature(const unsigned int np,
                                              const unsigned int nv)
  {
    quadrature_data.reinit(np, nv);
  }


  template <typename number>
  inline unsigned int
  LocalResults<number>::n_values() const
  {
    return J.size();
  }


  template <typename number>
  inline unsigned int
  LocalResults<number>::n_vectors() const
  {
    return R.size();
  }


  template <typename number>
  inline unsigned int
  LocalResults<number>::n_matrices() const
  {
    return M1.size();
  }


  template <typename number>
  inline unsigned int
  LocalResults<number>::n_quadrature_points() const
  {
    return quadrature_data.n_rows();
  }


  template <typename number>
  inline unsigned int
  LocalResults<number>::n_quadrature_values() const
  {
    return quadrature_data.n_cols();
  }


  template <typename number>
  inline number &
  LocalResults<number>::value(const unsigned int i)
  {
    AssertIndexRange(i, J.size());
    return J[i];
  }


  template <typename number>
  inline BlockVector<number> &
  LocalResults<number>::vector(const unsigned int i)
  {
    AssertIndexRange(i, R.size());
    return R[i];
  }


  template <typename number>
  inline MatrixBlock<FullMatrix<number>> &
  LocalResults<number>::matrix(const unsigned int i, const bool external)
  {
    if (external)
      {
        AssertIndexRange(i, M2.size());
        return M2[i];
      }
    AssertIndexRange(i, M1.size());
    return M1[i];
  }


  template <typename number>
  inline number &
  LocalResults<number>::quadrature_value(const unsigned int k,
                                         const unsigned int i)
  {
    return quadrature_data(k, i);
  }


  template <typename number>
  inline Table<2, number> &
  LocalResults<number>::quadrature_values()
  {
    return quadrature_data;
  }


  template <typename number>
  inline number
  LocalResults<number>::value(const unsigned int i) const
  {
    AssertIndexRange(i, J.size());
    return J[i];
  }


  template <typename number>
  inline const BlockVector<number> &
  LocalResults<number>::vector(const unsigned int i) const
  {
    AssertIndexRange(i, R.size());
    return R[i];
  }


  template <typename number>
  inline const MatrixBlock<FullMatrix<number>> &
  LocalResults<number>::matrix(const unsigned int i, const bool external) const
  {
    if (external)
      {
        AssertIndexRange(i, M2.size());
        return M2[i];
      }
    AssertIndexRange(i, M1.size());
    return M1[i];
  }


  template <typename number>
  inline number
  LocalResults<number>::quadrature_value(const unsigned int k,
                                         const unsigned int i) const
  {
    return quadrature_data(k, i);
  }


  template <typename number>
  template <class StreamType>
  void
  LocalResults<number>::print_debug(StreamType &os) const
  {
    os << "J: " << J.size() << std::endl;
    os << "R: " << R.size() << std::endl;
    for (unsigned int i = 0; i < R.size(); ++i)
      {
        os << "  " << R[i].n_blocks() << " -";
        for (unsigned int j = 0; j < R[i].n_blocks(); ++j)
          os << ' ' << R[i].block(j).size();
        os << std::endl;
      }
    os << "M: " << M1.size() << " face " << M2.size() << std::endl;
    for (unsigned int i = 0; i < M1.size(); ++i)
      {
        os << "  " << M1[i].row << "," << M1[i].column << " "
           << M1[i].matrix.m() << 'x' << M1[i].matrix.n();
        if (i < M2.size())
          os << " face " << M2[i].row << "," << M2[i].column << " "
             << M2[i].matrix.m() << 'x' << M2[i].matrix.n();
        os << std::endl;
      }
  }

} // namespace MeshWorker


DEAL_II_NAMESPACE_CLOSE

#endif


