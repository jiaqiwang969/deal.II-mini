//include/deal.II-translator/base/hdf5_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2021 by the deal.II authors
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

#ifndef dealii_hdf5_h
#define dealii_hdf5_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_HDF5

#  include <deal.II/base/array_view.h>

#  include <deal.II/lac/full_matrix.h>

#  include <hdf5.h>

#  include <numeric>

DEAL_II_NAMESPACE_OPEN

// It is necessary to turn clang-format off in order to maintain the Doxygen
// links because they are longer than 80 characters
// clang-format off
/**
 * 包含deal.II的HDF5接口的命名空间。
 * [层次数据格式（HDF）]（https://www.hdfgroup.org/）是一种跨平台和高I/O性能的格式，旨在存储大量的数据。它支持串行和MPI
 * I/O访问。这组类提供了一个与[HDF5库](https://www.hdfgroup.org/downloads/hdf5/)的接口。
 * 教程 step-62 展示了如何使用deal.II的HDF5接口。 #
 * 组、数据集和属性
 * 一个HDF5文件被组织在[组](https://bitbucket.hdfgroup.org/pages/HDFFV/hdf5doc/master/browse/html/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Groups/HDF5_Groups.htm)和[数据集](https://bitbucket.hdfgroup.org/pages/HDFFV/hdf5doc/master/browse/html/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Datasets/HDF5_Datasets.htm)中。组可以包含数据集和其他组。数据集是由数据元素的集合组成的对象。数据集等同于张量和矩阵。此外，属性可以被附加到根文件、组或数据集。一个[HDF5属性]（https://bitbucket.hdfgroup.org/pages/HDFFV/hdf5doc/master/browse/html/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Attributes/HDF5_Attributes.htm）是一个小型的元数据。方法
 * HDF5Object::get_attribute() 和 HDF5Object::set_attribute()
 * 可以用来获取和设置属性。 一个例子显示如下
 *
 * @code
 * HDF5::File data_file(filename, HDF5::File::FileAccessMode::create);
 * double double_attribute = 2.2;
 * data_file.set_attribute("double_attribute", double_attribute);
 * auto group = data_file.create_group("group");
 * group.set_attribute("simulation_type", std::string("elastic_equation"));
 * auto dataset = group.create_dataset<double>("dataset_name", dimensions);
 * dataset.set_attribute("complex_double_attribute",
 *                     std::complex<double>(2,2.3));
 * @endcode
 *
 * # MPI I/O
 * 一个HDF5文件可以用串行（一个单一进程）或MPI支持（几个进程访问同一个HDF5文件）来打开/创建。
 * File::File(const   std::string  &, const
 * FileAccessMode）为串行操作打开/创建一个HDF5文件。
 * File::File(const   std::string  &, const FileAccessMode, const MPI_Comm &)
 * 使用MPI并行地创建或打开一个HDF5文件。修改文件结构的HDF5调用总是集体进行的，而数据集中的原始数据的写入和读取可以独立进行，也可以集体进行。[集体访问通常更快](https://www.hdfgroup.org/2015/08/parallel-io-with-hdf5/)，因为它允许MPI进行优化。在deal.II的HDF5接口中，为了最大限度地提高性能，所有的调用都被设置为集体调用。这意味着所有的MPI进程都必须对每一次调用做出贡献，即使他们没有数据需要写入。MPI
 * HDF5要求deal.II和HDF5已经被编译为MPI支持。 ##
 * 写一个并行的hyperslab
 * Hyperslab是数据集的一部分。一个hyperslab可以是一个数据集中连续的点的集合，也可以是一个数据集中有规律的点或块的模式。Hyperslabs等同于python
 * numpy和h5py
 * [slices](http://docs.h5py.org/en/latest/high/dataset.html#reading-writing-data)。参见HDF5用户指南中的<a
 * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Dataspaces/HDF5_Dataspaces_and_Partial_I_O.htm?rhtocid=7.2#TOC_7_4_Dataspaces_and_Databc-6">Dataspaces
 * and Data Transfer</a>部分。也可参见<a
 * href="https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab">H5Sselect_hyperslab
 * definition</a>。
 * 下面的例子显示了如何编写一个简单的矩形超文本。偏移量定义了原始数据集中超文本的原点。超文本的尺寸是`hyperslab_dimensions
 * ={2,
 * 5}`。注意，每个进程可以写一个不同尺寸的超文本。如果一个进程根本不写任何数据，该进程应该调用函数
 * DataSet::write_none()
 * ，因为该操作是集体的*，所有MPI进程都必须为该调用作出贡献，即使它们没有数据可写。
 *
 * @code
 * std::vector<hsize_t> dataset_dimensions = {50, 30};
 * auto dataset = group.create_dataset<double>("name", dataset_dimensions);
 * if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 * {
 *   // hyperslab_data can be std::vector, FullMatrix or Vector
 *   FullMatrix<double> hyperslab_data = {...};
 *   std::vector<hsize_t> hyperslab_offset     = {1, 2};
 *   std::vector<hsize_t> hyperslab_dimensions = {2, 3};
 *   dataset.write_hyperslab(hyperslab_data,
 *                           hyperslab_offset,
 *                           hyperslab_dimensions);
 * }
 * else
 * {
 *   dataset.write_none<double>();
 * }
 * @endcode
 *
 * 函数 DataSet::write_hyperslab(const  Container &,const
 * std::vector<hsize_t>  &, const  std::vector<hsize_t>
 * &)用于写简单的超板，函数 DataSet::write_hyperslab(const
 * Container &,const  std::vector<hsize_t>  &, const  std::vector<hsize_t>  &,
 * const  std::vector<hsize_t>  &, const  std::vector<hsize_t>  &, const
 * std::vector<hsize_t>  &)用于写复杂超板。 ##并行写入无序数据
 * 下面的例子显示了如何写入一个选择的数据。请注意，每个进程可以写入不同数量的数据。如果一个进程根本不写任何数据，该进程应该调用函数
 * DataSet::write_none()
 * ，因为该操作是collective*，所有的MPI进程都必须为该调用作出贡献，即使他们没有数据可写。一个更详细的例子可以在
 * step-62 中找到。
 *
 * @code
 * std::vector<hsize_t> dataset_dimensions = {50, 30};
 * auto dataset = group.create_dataset<double>("name", dataset_dimensions);
 *
 * if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 * {
 *   std::vector<hsize_t> coordinates = {0,
 *                                       0, // first point
 *                                       0,
 *                                       2, // second point
 *                                       3,
 *                                       4, // third point
 *                                       25,
 *                                       12}; // fourth point
 *   std::vector<double>  data        = {2, 3, 5, 6};
 *   dataset.write_selection(data, coordinates);
 * }
 * else if (Utilities::MPI::this_mpi_process(mpi_communicator) == 1)
 * {
 *   std::vector<hsize_t> coordinates = {5,
 *                                       0, // first point
 *                                       0,
 *                                       4, // second point
 *                                       5,
 *                                       4, // third point
 *                                       26,
 *                                       12}; // fourth point
 *   std::vector<double>  data        = {9, 4, 7, 6};
 *   dataset.write_selection(data, coordinates);
 * }
 * else
 * {
 *   dataset.write_none<double>();
 * }
 * @endcode
 *
 * ## 查询HDF5在最后一次并行I/O调用中使用的I/O模式
 * 在deal.II的HDF5
 * C++接口中，默认的访问模式是集体访问，这通常会更快，因为它允许MPI做更多的优化。在某些情况下，比如有类型转换时，HDF5库可以决定做独立的I/O而不是集体I/O，即使用户要求集体I/O。参见下面的[文章]（https://www.hdfgroup.org/2015/08/parallel-io-with-hdf5/）。在需要最大性能的情况下，确保所有的MPI读/写操作都是集体的，这很重要。HDF5库提供了API例程，可以在读/写I/O操作之后使用，以查询I/O模式。如果
 * DataSet::query_io_mode
 * 为True，那么在每次读/写操作之后，deal.II的HDF5接口都会调用例程[H5Pget_mpio_actual_io_mode()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioActualIoMode)和[H5Pget_mpio_no_collective_cause()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause)
 * 。结果存储在  DataSet::io_mode,   DataSet::local_no_collective_cause
 * 和  DataSet::get_global_no_collective_cause.
 * 我们建议只在调试模式下查询I/O模式，因为它需要调用额外的HDF5程序。
 * 以下代码可用于查询I/O方式。
 *
 * @code
 * auto dataset = group.create_dataset<double>("name", dimensions);
 * #ifdef DEBUG
 * dataset.set_query_io_mode(true);
 * #endif
 *
 * if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 * {
 *   dataset.write(data);
 * }
 * else
 * {
 *   dataset.write_none<double>();
 * }
 *
 * if(dataset.get_query_io_mode()){
 * pcout << "IO mode: " << dataset.io_mode() << std::endl;
 * pcout << "Local no collective cause: "
 *       << dataset.local_no_collective_cause() << std::endl;
 * pcout << "Global no collective cause: "
 *       << dataset.get_global_no_collective_cause() <<
 * std::endl;
 * }
 * @endcode
 *
 * 如果写操作是集体的，那么输出应该是
 *
 * @code
 * IO mode: H5D_MPIO_CONTIGUOUS_COLLECTIVE
 * Local no collective cause: H5D_MPIO_COLLECTIVE
 * Global no collective cause: H5D_MPIO_COLLECTIVE
 * @endcode
 * 参见 DataSet::get_io_mode(),   DataSet::get_local_no_collective_cause()
 * 和  DataSet::get_global_no_collective_cause()
 * 所有可能的返回代码。 # HDF5数据集和hyperslabs的等级
 * deal.II的HDF5接口可以用来向任何特定等级的数据集和hyperslabs写入/读取数据。`FullMatrix`只能用于向等级为2的数据集和超文本写入/读取数据。另一方面，
 * `std::vector`
 * 和`Vector`可以用来向等级为1、2、3和更高的数据集和hyperslab写/读数据，数据是按照C和C++矩阵中常用的[row-major
 * order]（https://en.wikipedia.org/wiki/Row-_and_column-major_order）组织的。我们可以用
 * std::vector 重写上一节的代码
 *
 * @code
 * // Dataset of rank 2. dim_0 = 50, dim_1 = 30
 * std::vector<hsize_t> dataset_dimensions = {50, 30};
 * auto dataset = group.create_dataset<double>("name", dataset_dimensions);
 * if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
 * {
 *   // hyperslab_data can be std::vector, FullMatrix or Vector
 *   std::vector<double> hyperslab_data = {0,1,2,3,4,5};
 *   // hyperslab of rank 2. dim_0 = 2 and dim_1 = 3
 *   std::vector<hsize_t> hyperslab_offset     = {1, 2};
 *   std::vector<hsize_t> hyperslab_dimensions = {2, 3};
 *   dataset.write_hyperslab(hyperslab_data,
 *                           hyperslab_offset,
 *                           hyperslab_dimensions);
 * }
 * else
 * {
 *   dataset.write_none<double>();
 * }
 * @endcode
 * 前面的代码写出了以下的超简约矩阵
 *
 * @code
 * 0 1
 * 2 3
 * 4 5
 * @endcode
 *
 * # 数据类型 属性数据类型可以是float, `double`,
 * `std::complex<float>`,   `std::complex<double>`,  `int`, `unsigned int`,
 * `bool`和  `std::string`.   HDF5Object::get_attribute()  和
 * HDF5Object::set_attribute()  可以使用所有这些数据类型。
 * 数据集数据类型可以是`float`，`double`，  `std::complex<float>`,
 * `std::complex<double>`,  `int`和`unsigned int`。  DataSet::read(),
 * DataSet::write(),   DataSet::read_selection(),
 * 等可以与所有这些数据类型一起使用。注意，数据集数据类型不能是`bool'，原因是不能假设
 * `std::vector<bool>` 以连续的方式存储元素。
 *
 *  ## 复数和HDF5 在HDF5文件中没有正式的HDF5格式来存储
 * `std::complex` 数字。但是事实上*的标准是将 `std::complex`
 * 数字存储在一个复合类型中，其中`r`对应实部，`i`对应虚部。在这个接口中，我们定义了两个复合类型，一个用于
 * `std::complex<double>` ，对应于`(double,double)`，另一个用于
 * `std::complex<float>`
 * ，对应于`（float,float）`。这两种类型分别对应于python/numpy/h5py的类型：`complex128`和`complex64`。这意味着本接口生成的文件将被python/numpy/h5py正确读取，同时本接口能够读取python/numpy/h5py生成的文件。
 * # 与python脚本交换数据
 * HDF5格式可以用来与python脚本交换数据。字符串被存储为HDF5可变长度的UTF-8字符串，复数，如上所述，被存储为HDF5复合数据类型，与[h5py](https://www.h5py.org/)和[numpy](http://www.numpy.org/)兼容。
 * 下面的python脚本写了一个deal.II模拟的参数。~~~~~~~~~~~~~{.py}
 * h5_file = h5py.File('simulation.hdf5','w') data =
 * h5_file.create_group('data') data.attrs['nb_frequency_points'] = 50 # int
 * data.attrs['rho'] = 2300.5 # double data.attrs['save_vtk_files'] = True #
 * bool data.attrs['simulation_type'] = 'elastic_equation' # utf8 string
 * ~~~~~~~~~~~~~ 用MPI HDF5进行C++ deal.II仿真。
 *
 * @code
 * HDF5::File data_file("simulation.hdf5",
 *                    HDF5::File::FileAccessMode::open,
 *                    MPI_COMM_WORLD);
 * HDF5::Group data = data_file.open_group("data");
 *
 * auto nb_frequency_points = data.get_attribute<int>("nb_frequency_points");
 * auto rho = data.get_attribute<double>("rho");
 * auto save_vtk_files = data.get_attribute<bool>("save_vtk_files");
 * auto simulation_type = data.get_attribute<std::string>("simulation_type");
 *
 * std::vector<std::complex<double>> displacement = {...};
 *
 * data.write_dataset("displacement", displacement);
 *
 * // Write the simulation metadata
 * data.set_attribute("active_cells", triangulation.n_active_cells());
 * @endcode
 *
 * 用python读取模拟结果。~~~~~~~~~~~~~{.py} h5_file =
 * h5py.File('simulation.hdf5','r+') data = h5_file['data'] displacement =
 * data['displacement'] # complex128 dtype active_cells =
 * data.attrs['deges_of_freedom']) ~~~~~~~~~~~~~ # HDF5和线程安全
 * 默认情况下，HDF5不是线程安全的。HDF5库可以被配置为线程安全的，参见[HDF5文档](https://support.hdfgroup.org/HDF5/faq/threadsafe.html)。线程安全的HDF5版本将API序列化，但不提供任何级别的并发性。为了实现HDF5的高并行性能，我们建议将HDF5与MPI一起使用。
 *
 *
 */
// clang-format on
namespace HDF5
{
  /**
   * HDF5对象的基类。
   *
   */
  class HDF5Object
  {
  protected:
    /**
     * 构造函数。  @p string_name  是HDF5对象的名称。如果 @p mpi
     * 为True，则使用MPI I/O。
     *
     */
    HDF5Object(const std::string &name, const bool mpi);

  public:
    /**
     * 读取一个属性。  @p T  可以是`float`, `double`,
     * `std::complex<float>`,   `std::complex<double>`,  `int`, `unsigned
     * int`, `bool` 或  `std::string`.  注意  `std::string`
     * 的编码是UTF8，以便与python3兼容。
     * 数据类型转换在读或写时进行，是自动的。参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>部分。
     *
     */
    template <typename T>
    T
    get_attribute(const std::string &attr_name) const;

    /**
     * 写入一个属性。  @p T  可以是`float`, `double`,
     * `std::complex<float>`,   `std::complex<double>`,  `int`, `unsigned
     * int`, `bool` 或  `std::string`.  注意，为了与python3兼容，
     * `std::string`  的编码为UTF8。
     * 数据类型转换在读或写时进行，是自动的。参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>部分。
     *
     */
    template <typename T>
    void
    set_attribute(const std::string &attr_name, const T value);

    /**
     * 返回对象的#名称。在文件的情况下，#name对应的是文件名。在Group和DataSet的情况下，#name对应于HDF5文件中的对象名称。
     *
     */
    std::string
    get_name() const;

  protected:
    /**
     * HDF5Oject的名称。在文件的情况下， @p name
     * 对应于文件名。在Group和DataSet的情况下， @p name
     * 对应于HDF5文件中的对象名称。
     *
     */
    const std::string name;

    /**
     * 文件、组和数据集对象的HDF5标识符。 `std::shared_ptr<>`
     * 指针允许对象被复制。例如，程序的几个部分可以共享和访问同一个组；当所有访问该组的函数关闭时，该组的HDF5资源将被自动释放。
     *
     */
    std::shared_ptr<hid_t> hdf5_reference;

    /**
     * 如果为真则使用并行HDF5，如果为假则使用串行HDF5。
     *
     */
    const bool mpi;
  };

  /**
   * 这个类实现了一个HDF5数据集。
   *
   */
  class DataSet : public HDF5Object
  {
    friend class Group;

  protected:
    /**
     * 打开数据集。这是一个内部构造函数。应该使用函数
     * Group::open_dataset() 来打开一个数据集。
     *
     */
    DataSet(const std::string &name, const hid_t &parent_group_id, bool mpi);

    /**
     * 创建数据集。这是一个内部构造函数。应该使用函数
     * Group::create_dataset() 来创建一个数据集。
     *
     */
    DataSet(const std::string &           name,
            const hid_t &                 parent_group_id,
            const std::vector<hsize_t> &  dimensions,
            const std::shared_ptr<hid_t> &t_type,
            const bool                    mpi);

  public:
    /**
     * 读取数据集的所有数据。
     * 数据类型转换在读取操作时进行，是自动的。参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>部分。        容器
     * "可以是 `std::vector<float>`,   `std::vector<double>`,
     * `std::vector<std::complex<float>>`,
     * `std::vector<std::complex<double>>`,   `std::vector<int>`,
     * `std::vector<unsigned  int>`、`Vector<float>`、`Vector<double>`、
     * `Vector<std::complex<float>>`,   `Vector<std::complex<double>>`,
     * `FullMatrix<float>`、`FullMatrix<double>`、
     * `FullMatrix<std::complex<float>>`  或
     * `FullMatrix<std::complex<double>>`.
     *
     */
    template <typename Container>
    Container
    read();

    /**
     * 读取数据集的一个子集的数据。
     * 数据类型转换在读取操作时进行，并且是自动的。参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>部分。
     * 选定的元素可以是分散的，在数据集中采取任何形状。
     * 例如，在一个等级为4的数据集的情况下，选择3个点将由一个3乘4的数组来描述。请注意，索引是基于零的。要选择(1,1,1,1)、(14,6,12,18)和(8,22,30,22)这几个点，点选择阵列将如下。
     * @code
     *  0  0  0  0
     *
     * 13  5 11 17
     *
     *  7 21 29 21
     * @endcode
     * <a
     * href="https://support.hdfgroup.org/newsletters/newsletter140.html">Parallel
     * HDF5 supports collective I/O on point selections.</a>
     * 数据类型转换在读操作时进行，是自动的。参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>部分。
     *
     */
    template <typename Container>
    Container
    read_selection(const std::vector<hsize_t> &coordinates);

    // clang-format off
    /**
     * 从数据集中读取一个超文本。参数总结如下。
     *
     *
     *
     *
     *
     * -  @p offset:  超级板块的起始位置。
     *
     *
     *
     *
     *
     *
     * -  @p count:  沿着每个维度要选择的元素数量。        当读取一个超文本时，HDF5也允许提供 "stride "和 "block "参数（参见[HDF5文档](https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab)）。    这些参数不会被当前函数使用，而是被设置为 "nullptr"。然而，这些参数可以在函数read_hyperslab(const  std::vector<hsize_t>  &, const  std::vector<hsize_t>  &, const  std::vector<hsize_t>  &, const  std::vector<hsize_t>  &, const  std::vector<hsize_t>  &) 中使用。        参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Dataspaces/HDF5_Dataspaces_and_Partial_I_O.htm?rhtocid=7.2#TOC_7_4_Dataspaces_and_Databc-6">Dataspaces
     * and Data Transfer</a>部分。也可以参见<a
     * href="https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab">H5Sselect_hyperslab
     * definition</a>。
     * 数据类型转换是在读或写的时候进行的，是自动的。参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>部分。        容器
     * "可以是 `std::vector<float>`,   `std::vector<double>`,
     * `std::vector<std::complex<float>>`,
     * `std::vector<std::complex<double>>`,   `std::vector<int>`,
     * `std::vector<unsigned  int>`、`Vector<float>`、`Vector<double>`、
     * `Vector<std::complex<float>>`,   `Vector<std::complex<double>>`,
     * `FullMatrix<float>`、`FullMatrix<double>`、
     * `FullMatrix<std::complex<float>>`  或
     * `FullMatrix<std::complex<double>>`.
     *
     */
    // clang-format on
    template <typename Container>
    Container
    read_hyperslab(const std::vector<hsize_t> &offset,
                   const std::vector<hsize_t> &count);

    /**
     * 向数据集写入一个数据超文本。参数总结如下。
     *
     *
     *
     *
     *
     * -  @p dataset_dimensions:  数据存储块的尺寸。
     *
     *
     *
     *
     *
     * -  @p offset:  超级板块的起始位置。
     *
     *
     *
     *
     *
     *
     * -  @p stride:  分隔每个要选择的元素或块的元素数量。
     *
     *
     *
     *
     *
     *
     * -  @p count:  沿着每个维度要选择的元素或块的数量。
     *
     *
     *
     *
     *
     *
     * -  @p block:  从数据空间选择的块的大小。        参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Dataspaces/HDF5_Dataspaces_and_Partial_I_O.htm?rhtocid=7.2#TOC_7_4_Dataspaces_and_Databc-6">Dataspaces
     * and Data Transfer</a>部分。也可参见<a
     * href="https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab">H5Sselect_hyperslab
     * definition</a>。
     * 数据类型转换是在读或写的时候进行的，而且是自动的。参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>部分。        容器
     * "可以是 `std::vector<float>`,   `std::vector<double>`,
     * `std::vector<std::complex<float>>`,
     * `std::vector<std::complex<double>>`,   `std::vector<int>`,
     * `std::vector<unsigned  int>`、`Vector<float>`、`Vector<double>`、
     * `Vector<std::complex<float>>`,   `Vector<std::complex<double>>`,
     * `FullMatrix<float>`、`FullMatrix<double>`、
     * `FullMatrix<std::complex<float>>`  或
     * `FullMatrix<std::complex<double>>`.
     *
     */
    template <typename Container>
    Container
    read_hyperslab(const std::vector<hsize_t> &data_dimensions,
                   const std::vector<hsize_t> &offset,
                   const std::vector<hsize_t> &stride,
                   const std::vector<hsize_t> &count,
                   const std::vector<hsize_t> &block);

    /**
     * 这个函数不读取任何数据，但它可以为集体读取调用作出贡献。
     * @p number  可以是`float`, `double`,  `std::complex<float>`,
     * `std::complex<double>`,  `int`或`unsigned int`。
     * 数据类型的转换是在读或写时进行的，并且是自动的。参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>部分。
     *
     */
    template <typename number>
    void
    read_none();

    /**
     * 写入数据集中的数据。  @p number
     * 可以是`float`，`double`， `std::complex<float>`,
     * `std::complex<double>`,  `int`或`unsigned int`。
     * 数据类型的转换是在读或写时进行的，并且是自动的。参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>部分。        容器
     * "可以是 `std::vector<float>`,   `std::vector<double>`,
     * `std::vector<std::complex<float>>`,
     * `std::vector<std::complex<double>>`,   `std::vector<int>`,
     * `std::vector<unsigned  int>`、`Vector<float>`、`Vector<double>`、
     * `Vector<std::complex<float>>`,   `Vector<std::complex<double>>`,
     * `FullMatrix<float>`、`FullMatrix<double>`、
     * `FullMatrix<std::complex<float>>`  或
     * `FullMatrix<std::complex<double>>`.
     *
     */
    template <typename Container>
    void
    write(const Container &data);

    /**
     * 将数据写入数据集的一个子集。  @p number
     * 可以是`float`, `double`,  `std::complex<float>`,
     * `std::complex<double>`,  `int`或`unsigned int`。
     * 选择的元素可以是分散的，在数据集中采取任何形状。
     * 例如，在一个等级为4的数据集的情况下，3个点的选择将由一个3乘4的数组来描述。请注意，索引是基于零的。为了选择(1,1,1,1)、(14,6,12,18)和(8,22,30,22)这几个点，点选择阵列将如下。
     * @code
     *  0  0  0  0
     *
     * 13  5 11 17
     *
     *  7 21 29 21
     * @endcode
     * <a
     * href="https://support.hdfgroup.org/newsletters/newsletter140.html">Parallel
     * HDF5 supports collective I/O on point selections.</a>
     * 数据类型转换在读或写时进行，是自动的。参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>部分。
     *
     */
    template <typename Container>
    void
    write_selection(const Container &           data,
                    const std::vector<hsize_t> &coordinates);

    // clang-format off
    /**
     * 向数据集写入一个数据超文本。参数总结如下。
     *
     *
     *
     *
     *
     * -  @p offset:  超文本的起始位置。
     *
     *
     *
     *
     *
     *
     * -  @p count:  沿着每个维度要选择的元素数量。        在编写超文本时，HDF5还允许提供 "stride "和 "block "参数（见[HDF5文档](https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab)）。    这些参数不会被当前函数使用，而是被设置为 "nullptr"。但是这些参数可以在函数write_hyperslab(const Container &data, const  std::vector<hsize_t>  &data_dimensions, const  std::vector<hsize_t>  &offset, const  std::vector<hsize_t>  &stride, const  std::vector<hsize_t>  &count, const  std::vector<hsize_t>  &block) 中使用。        参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Dataspaces/HDF5_Dataspaces_and_Partial_I_O.htm?rhtocid=7.2#TOC_7_4_Dataspaces_and_Databc-6">Dataspaces
     * and Data Transfer</a>部分。也可参见<a
     * href="https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab">H5Sselect_hyperslab
     * definition</a>。
     * 数据类型转换是在读或写的时候进行的，而且是自动的。参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>部分。
     *
     */
    // clang-format on
    template <typename Container>
    void
    write_hyperslab(const Container &           data,
                    const std::vector<hsize_t> &offset,
                    const std::vector<hsize_t> &count);

    /**
     * 向数据集写入一个数据超文本。参数总结如下。
     *
     *
     *
     *
     *
     * -  @p dataset_dimensions:  数据存储块的尺寸。
     *
     *
     *
     *
     *
     *
     * -  @p offset:  超级板块的起始位置。
     *
     *
     *
     *
     *
     *
     * -  @p stride:  分隔每个要选择的元素或块的元素数量。
     *
     *
     *
     *
     *
     *
     * -  @p count:  沿着每个维度要选择的元素或块数。
     *
     *
     *
     *
     *
     *
     * -  @p block:  从数据空间中选择的块的大小。        参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/HDF5_Users_Guide/Dataspaces/HDF5_Dataspaces_and_Partial_I_O.htm?rhtocid=7.2#TOC_7_4_Dataspaces_and_Databc-6">Dataspaces
     * and Data Transfer</a>部分。也可参见<a
     * href="https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-SelectHyperslab">H5Sselect_hyperslab
     * definition</a>。
     * 数据类型转换是在读或写的时候进行的，而且是自动的。参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>部分。        容器
     * "可以是 `std::vector<float>`,   `std::vector<double>`,
     * `std::vector<std::complex<float>>`,
     * `std::vector<std::complex<double>>`,   `std::vector<int>`,
     * `std::vector<unsigned  int>`、`Vector<float>`、`Vector<double>`、
     * `Vector<std::complex<float>>`,   `Vector<std::complex<double>>`,
     * `FullMatrix<float>`、`FullMatrix<double>`、
     * `FullMatrix<std::complex<float>>`  或
     * `FullMatrix<std::complex<double>>`.
     *
     */
    template <typename Container>
    void
    write_hyperslab(const Container &           data,
                    const std::vector<hsize_t> &data_dimensions,
                    const std::vector<hsize_t> &offset,
                    const std::vector<hsize_t> &stride,
                    const std::vector<hsize_t> &count,
                    const std::vector<hsize_t> &block);

    /**
     * 这个函数不写任何数据，但它可以为集体写调用作出贡献。在MPI集体写调用的情况下，如果一个进程根本不写任何数据，该进程应该调用这个函数，因为该操作是集体的*，所有的MPI进程都必须为该调用作出贡献，即使他们没有数据可写。
     * @p number  可以是 "float"、"double"、 `std::complex<float>`,
     * `std::complex<double>`,  `int'或`unsigned int'。
     * 数据类型的转换是在读或写时进行的，并且是自动的。参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>部分。
     * 如何使用这个函数的例子可以在  step-62  中找到。
     *
     */
    template <typename number>
    void
    write_none();

    /**
     * 该函数返回boolean query_io_mode。
     * 在必须实现最大性能的情况下，确保所有的MPI读/写操作是集体的，这一点很重要。HDF5库提供了API例程，可以在读/写I/O操作之后使用，以查询I/O模式。如果query_io_mode设置为true，那么在每次读/写操作之后，deal.II的HDF5接口都会调用例程[H5Pget_mpio_actual_io_mode()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioActualIoMode)和[H5Pget_mpio_no_collective_cause()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause)
     * 。
     * 结果存储在io_mode、local_no_collective_cause和global_no_collective_cause。我们建议只在Debug模式下查询I/O模式，因为这需要调用额外的HDF5例程。
     *
     */
    bool
    get_query_io_mode() const;

    /**
     * 这个函数设置布尔查询_io_mode。
     *
     */
    void
    set_query_io_mode(const bool new_query_io_mode);

    /**
     * 该函数返回最后一次并行I/O调用时使用的I/O模式。参见<a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioActualIoMode">H5Pget_mpio_actual_io_mode</a>。
     * 返回值是一个 `std::string` ，可以是Value | Meaning。
     *
     *
     *
     *
     *
     *
     * - ---------------------------- |
     *
     * - ----- H5D_MPIO_NO_COLLECTIVE | 没有执行集体I/O。没有要求进行集体I/O，或者在这个数据集上不可能进行集体I/O。    H5D_MPIO_CHUNK_INDEPENDENT | HDF5执行了分块集体优化方案，每个分块被独立访问。    H5D_MPIO_CHUNK_COLLECTIVE | HDF5执行大块集体优化方案，每个大块被集体访问。    H5D_MPIO_CHUNK_MIXED | HDF5执行大块集体优化方案，有些大块被独立访问，有些被集体访问。    H5D_MPIO_CONTIGUOUS_COLLECTIVE | 对一个连续的数据集进行了集体I/O。
     *
     */
    std::string
    get_io_mode();

    /**
     * 该函数返回最后一次并行I/O调用时使用的I/O模式。参见<a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioActualIoMode">H5Pget_mpio_actual_io_mode</a>。
     * 返回类型为`H5D_mpio_actual_io_mode_t`，对应于H5Pget_mpio_actual_io_mode的返回值。
     * 返回值可以是Value | Meaning
     *
     *
     *
     *
     *
     *
     * - ---------------------------- |
     *
     * - ----- H5D_MPIO_NO_COLLECTIVE | 没有执行集体I/O。没有要求进行集体I/O，或者在这个数据集上不可能进行集体I/O。    H5D_MPIO_CHUNK_INDEPENDENT | HDF5进行了分块集体优化，每个分块被独立访问。    H5D_MPIO_CHUNK_COLLECTIVE | HDF5执行大块集体优化，每个大块被集体访问。    H5D_MPIO_CHUNK_MIXED | HDF5执行大块集体优化，有些大块被独立访问，有些被集体访问。    H5D_MPIO_CONTIGUOUS_COLLECTIVE | 对一个连续的数据集进行了集体I/O。
     *
     */
    H5D_mpio_actual_io_mode_t
    get_io_mode_as_hdf5_type();

    /**
     * 这个函数返回在最后一次并行I/O调用中破坏集体I/O的本地原因。见<a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">H5Pget_mpio_no_collective_cause</a>。
     * 返回值是一个字符串，可以是 Value | Meaning
     *
     *
     *
     *
     *
     *
     * - ---------------------------------------- |
     *
     * - ----- H5D_MPIO_COLLECTIVE | 集体I/O已成功执行。    H5D_MPIO_SET_INDEPENDENT | 因为要求独立的I/O，所以没有执行集体I/O。    H5D_MPIO_DATATYPE_CONVERSION | 因为需要进行数据类型转换，所以没有执行集体I/O。    H5D_MPIO_DATA_TRANSFORMS | 因为需要应用数据转换，所以没有执行集体I/O。    H5D_MPIO_SET_MPIPOSIX | 因为选择的文件驱动是MPI-POSIX，所以没有执行集体I/O。    H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES | 因为其中一个数据空间既不简单也不标量，所以没有执行集体I/O。    H5D_MPIO_POINT_SELECTIONS | 由于其中一个数据空间中存在点选择，所以没有执行集体I/O。    H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET | 由于数据集既不是连续的也不是分块的，所以没有执行集体I/O。    H5D_MPIO_FILTERS | 因为需要应用过滤器而没有执行集体I/O。
     *
     */
    std::string
    get_local_no_collective_cause();

    /**
     * 这个函数返回在最后一次并行I/O调用中破坏集体I/O的本地原因。参见<a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">H5Pget_mpio_no_collective_cause</a>。
     * 返回类型为`uint32_t`，对应于[H5Pget_mpio_no_collective_cause](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause)返回的值。
     * 返回值可以是Value | Meaning
     *
     *
     *
     *
     *
     *
     * - ---------------------------------------- |
     *
     * - ----- H5D_MPIO_COLLECTIVE | 集体I/O已成功执行。    H5D_MPIO_SET_INDEPENDENT | 因为要求独立的I/O，所以没有执行集体I/O。    H5D_MPIO_DATATYPE_CONVERSION | 因为需要进行数据类型转换，所以没有执行集体I/O。    H5D_MPIO_DATA_TRANSFORMS | 因为需要应用数据转换，所以没有执行集体I/O。    H5D_MPIO_SET_MPIPOSIX | 因为选择的文件驱动是MPI-POSIX，所以没有执行集体I/O。    H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES | 因为其中一个数据空间既不简单也不标量，所以没有执行集体I/O。    H5D_MPIO_POINT_SELECTIONS | 由于其中一个数据空间中存在点选择，所以没有执行集体I/O。    H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET | 由于数据集既不是连续的也不是分块的，所以没有执行集体I/O。    H5D_MPIO_FILTERS | 因为需要应用过滤器而没有执行集体I/O。
     *
     */
    uint32_t
    get_local_no_collective_cause_as_hdf5_type();

    /**
     * 这个函数检索在最后一次并行I/O调用中破坏集体I/O的全局原因。见<a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">H5Pget_mpio_no_collective_cause</a>。
     * 返回值是一个 std::string ，可以是Value | Meaning。
     *
     *
     *
     *
     *
     *
     * - ---------------------------------------- |
     *
     * - ----- H5D_MPIO_COLLECTIVE | 集体I/O已成功执行。    H5D_MPIO_SET_INDEPENDENT | 因为要求独立的I/O，所以没有执行集体I/O。    H5D_MPIO_DATATYPE_CONVERSION | 因为需要进行数据类型转换，所以没有执行集体I/O。    H5D_MPIO_DATA_TRANSFORMS | 因为需要应用数据转换，所以没有执行集体I/O。    H5D_MPIO_SET_MPIPOSIX | 因为选择的文件驱动是MPI-POSIX，所以没有执行集体I/O。    H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES | 因为其中一个数据空间既不简单也不标量，所以没有执行集体I/O。    H5D_MPIO_POINT_SELECTIONS | 由于其中一个数据空间中存在点选择，所以没有执行集体I/O。    H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET | 由于数据集既不是连续的也不是分块的，所以没有执行集体I/O。    H5D_MPIO_FILTERS | 因为需要应用过滤器而没有执行集体I/O。
     *
     */
    std::string
    get_global_no_collective_cause();

    /**
     * 这个函数返回在最后一次并行I/O调用中破坏集体I/O的全局原因。参见<a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">H5Pget_mpio_no_collective_cause</a>。
     * 返回类型为`uint32_t`，与H5Pget_mpio_no_collective_cause返回的值对应。
     * 返回值可以是 Value | Meaning
     *
     *
     *
     *
     *
     *
     * - ---------------------------------------- |
     *
     * - ----- H5D_MPIO_COLLECTIVE | 集体I/O已成功执行。    H5D_MPIO_SET_INDEPENDENT | 因为要求独立的I/O，所以没有执行集体I/O。    H5D_MPIO_DATATYPE_CONVERSION | 因为需要进行数据类型转换，所以没有执行集体I/O。    H5D_MPIO_DATA_TRANSFORMS | 因为需要应用数据转换，所以没有执行集体I/O。    H5D_MPIO_SET_MPIPOSIX | 因为选择的文件驱动是MPI-POSIX，所以没有执行集体I/O。    H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES | 因为其中一个数据空间既不简单也不标量，所以没有执行集体I/O。    H5D_MPIO_POINT_SELECTIONS | 由于其中一个数据空间中存在点选择，所以没有执行集体I/O。    H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET | 由于数据集既不是连续的也不是分块的，所以没有执行集体I/O。    H5D_MPIO_FILTERS | 因为需要应用过滤器，所以没有执行集体I/O。
     *
     */
    uint32_t
    get_global_no_collective_cause_as_hdf5_type();

    /**
     * 这个函数返回数据集的尺寸。向量dimensions是一个大小为rank的一维数组，指定数据集的每个维度的大小。
     *
     */
    std::vector<hsize_t>
    get_dimensions() const;

    /**
     * 该函数返回数据集的总元素数。
     *
     */
    unsigned int
    get_size() const;

    /**
     * 此函数返回数据集的等级。
     *
     */
    unsigned int
    get_rank() const;

  private:
    /**
     * 数据集的等级
     *
     */
    unsigned int rank;

    /**
     * 向量`dimensions`是一个大小为rank的一维数组，指定数据集的每个维度的大小。
     *
     */
    std::vector<hsize_t> dimensions;

    /**
     * HDF5数据空间标识符。
     *
     */
    std::shared_ptr<hid_t> dataspace;

    /**
     * 数据集的总元素数。
     *
     */
    unsigned int size;

    /**
     * 如果query_io_mode设置为true，那么在每次读/写操作之后，deal.II的HDF5接口都会调用例程[H5Pget_mpio_actual_io_mode()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioActualIoMode)和[H5Pget_mpio_no_collective_cause()](https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause)
     * 。
     * 结果存储在io_mode、local_no_collective_cause和global_no_collective_cause。
     *
     */
    bool query_io_mode;

    /**
     * 在最后一次并行I/O调用时执行的I/O模式。
     *
     */
    H5D_mpio_actual_io_mode_t io_mode;

    /**
     * 在最后一次并行I/O调用中破坏集体I/O的本地原因。见<a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">H5Pget_mpio_no_collective_cause</a>。
     *
     */
    uint32_t local_no_collective_cause;

    /**
     * 在最后一次并行I/O调用中破坏集体I/O的全局原因。见<a
     * href="https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause">H5Pget_mpio_no_collective_cause</a>。
     *
     */
    uint32_t global_no_collective_cause;
  };

  /**
   * 该类实现了一个HDF5组
   *
   */
  class Group : public HDF5Object
  {
  protected:
    /**
     * 组的访问模式
     *
     */
    enum class GroupAccessMode
    {
      /**
       * 打开一个现有的组
       *
       */
      open,
      /**
       * 创建一个新的群组
       *
       */
      create
    };
    /**
     * 这个构造函数创建或打开一个组，取决于 @p mode.
     * 的值，该组将被放在组内 @p parent_group. 的参数 @p mpi
     * 定义了I/O操作是串行还是并行。这是一个内部构造函数，应该使用当前类的函数open_group()和create_group()来打开或创建一个组。
     *
     */
    Group(const std::string &   name,
          const Group &         parent_group,
          const bool            mpi,
          const GroupAccessMode mode);

    /**
     * 文件使用的内部构造函数。该构造函数设置HDF5Group的受保护常量成员。
     * @p name  和  @p mpi.  它不会创建或打开一个组。
     *
     */
    Group(const std::string &name, const bool mpi);

  public:
    /**
     * 打开当前组或文件的一个子组。
     *
     */
    Group
    open_group(const std::string &name) const;

    /**
     * 在当前组或文件中创建一个子组。
     *
     */
    Group
    create_group(const std::string &name) const;

    /**
     * 打开一个数据集。
     *
     */
    DataSet
    open_dataset(const std::string &name) const;

    /**
     * 创建一个数据集。  @p number  可以是`float`, `double`,
     * `std::complex<float>`,   `std::complex<double>`,  `int`或`unsigned
     * int`。
     * 数据类型的转换是在读或写时进行的，并且是自动的。参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>部分。
     *
     */
    template <typename number>
    DataSet
    create_dataset(const std::string &         name,
                   const std::vector<hsize_t> &dimensions) const;

    /**
     * 创建并向数据集写入数据。  @p number
     * 可以是`float`、`double`、 `std::complex<float>`,
     * `std::complex<double>`,  `int`或`unsigned int`。
     * 数据类型的转换是在读或写时进行的，并且是自动的。参见HDF5用户指南中的<a
     * href="https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FDatatypes%2FHDF5_Datatypes.htm%23TOC_6_10_Data_Transferbc-26&rhtocid=6.5_2">Data
     * Transfer: Datatype Conversion and Selection</a>部分。        容器
     * "可以是 `std::vector<float>`,   `std::vector<double>`,
     * `std::vector<std::complex<float>>`,
     * `std::vector<std::complex<double>>`,   `std::vector<int>`,
     * `std::vector<unsigned  int>`、`Vector<float>`、`Vector<double>`、
     * `Vector<std::complex<float>>`,   `Vector<std::complex<double>>`,
     * `FullMatrix<float>`、`FullMatrix<double>`、
     * `FullMatrix<std::complex<float>>`  或
     * `FullMatrix<std::complex<double>>`.
     *
     */
    template <typename Container>
    void
    write_dataset(const std::string &name, const Container &data) const;
  };

  /**
   * 该类实现了一个HDF5文件
   *
   */
  class File : public Group
  {
  public:
    /**
     * 文件访问模式
     *
     */
    enum class FileAccessMode
    {
      /**
       * 读/写，文件必须存在
       *
       */
      open,
      /**
       * 创建文件，如果存在则截断
       *
       */
      create
    };

    /**
     * 创建或打开一个HDF5文件进行串行操作。这个调用不需要MPI支持。它创建或打开一个HDF5文件，取决于
     * @p mode. 的值。
     *
     */
    File(const std::string &name, const FileAccessMode mode);

    /**
     * 使用MPI并行地创建或打开一个HDF5文件。这需要deal.II和HDF5在编译时支持MPI。它创建或打开一个HDF5文件，取决于
     * @p mode.   @p mpi_communicator
     * 的值，定义了参与此调用的进程；`MPI_COMM_WORLD`是MPI通信器的一个通用值。
     *
     */
    File(const std::string &  name,
         const FileAccessMode mode,
         const MPI_Comm &     mpi_communicator);

  private:
    /**
     * 授权的内部构造函数。    File(const  std::string  &, const
     * MPI_Comm &, const Mode); 和 File(const  std::string  &, const Mode)
     * 应该被用来打开或创建HDF5文件。
     *
     */
    File(const std::string &  name,
         const FileAccessMode mode,
         const bool           mpi,
         const MPI_Comm &     mpi_communicator);
  };

  namespace internal
  {
    /**
     * 该函数返回与C++类型对应的HDF5数据类型。    在
     * std::complex 类型的情况下，HDF5处理程序使用
     * `std::shared_ptr`. 的析构器自动释放 `std::shared_ptr` 而不是
     * `std::unique_ptr` ，因为 `std::shared_ptr`
     * 的析构器不需要在模板参数中定义。另一方面，
     * `std::unique`
     * 的析构器必须在模板参数中定义。像`H5T_NATIVE_DOUBLE`这样的本地类型不需要析构器，但是像
     * std::complex<double>
     * 这样的复合类型需要一个析构器来释放HDF5资源。
     *
     */
    template <typename number>
    std::shared_ptr<hid_t>
    get_hdf5_datatype();

    /**
     * 返回`data`的尺寸。对于一个 std::vector 该函数返回
     * `std::vector<hsize_t>{vector_size}`.
     * 几个HDF5函数，如H5Screate_simple()需要一个一维数组，指定容器的每个维度的大小，见：https://support.hdfgroup.org/HDF5/doc1.8/RM/RM_H5S.html#Dataspace-CreateSimple
     *
     */
    template <typename number>
    std::vector<hsize_t>
    get_container_dimensions(const std::vector<number> &data);

    /**
     * 返回`data`的尺寸。对于一个矢量，这个函数返回
     * `std::vector<hsize_t>{vector_size}`. 。
     *
     */
    template <typename number>
    std::vector<hsize_t>
    get_container_dimensions(const Vector<number> &data);

    /**
     * 返回`data`的尺寸。对于FullMatrix，该函数返回
     * `std::vector<hsize_t>{rows, 列}`。
     *
     */
    template <typename number>
    std::vector<hsize_t>
    get_container_dimensions(const FullMatrix<number> &data);

    /**
     * 这个函数返回容器的总尺寸。对于一个 std::vector
     * ，该函数返回`int(vector_size)`。
     *
     */
    template <typename number>
    unsigned int
    get_container_size(const std::vector<number> &data);

    /**
     * 此函数返回容器的总大小。对于一个向量，该函数返回`int(vector_size)`。
     *
     */
    template <typename number>
    unsigned int
    get_container_size(const Vector<number> &data);

    /**
     * 此函数返回容器的总大小。对于FullMatrix，该函数返回`int(rows*columns)`。
     *
     */
    template <typename number>
    unsigned int
    get_container_size(const FullMatrix<number> &data);

    /**
     * 这个函数初始化并返回一个 std::vector,
     * Vector或FullMatrix类型的容器。该函数不设置容器中元素的值。该容器可以存储HDF5数据集或HDF5选择的数据。维度参数保存HDF5数据集或选择的维度。
     * 在 std::vector,
     * 的情况下，向量的大小将是由维度给出的总大小。例如，在一个等级为3的数据集的情况下，尺寸为
     * `std::vector<hsize_t>{dim_0,dim_1,dim_2}`. ，返回的 std::vector
     * 的大小将是`dim_0*dim_1*dim_2`。        如果是 dealii::Vector,
     * ，返回的 dealii::Vector 的大小也将是`dim_0*dim_1*dim_2'。
     * 一个FullMatrix只能存储等级为2的HDF5数据集的数据。FullMatrix的大小将是FullMatrix(dim_0,dim_2)
     *
     */
    template <typename Container>
    typename std::enable_if<
      std::is_same<Container,
                   std::vector<typename Container::value_type>>::value,
      Container>::type
    initialize_container(const std::vector<hsize_t> &dimensions);

    /**
     * 同上。
     *
     */
    template <typename Container>
    typename std::enable_if<
      std::is_same<Container, Vector<typename Container::value_type>>::value,
      Container>::type
    initialize_container(const std::vector<hsize_t> &dimensions);

    /**
     * 同上。
     *
     */
    template <typename Container>
    typename std::enable_if<
      std::is_same<Container,
                   FullMatrix<typename Container::value_type>>::value,
      Container>::type
    initialize_container(const std::vector<hsize_t> &dimensions);

    /**
     * 这个辅助函数设置了DataSet的读写操作的属性列表。必须为MPI驱动创建一个属性列表。对于串行驱动，可以使用默认的H5P_DEFAULT。此外，H5Pset_dxpl_mpio被用来设置MPI模式为集体模式。
     *
     */
    inline void
    set_plist(hid_t &plist, const bool mpi);

    /**
     * 这个辅助函数释放了DataSet的读写操作的属性列表处理程序。对于串行版本，不需要释放属性列表处理程序，因为已经使用了H5P_DEFAULT。如果query_io_mode为True，那么H5Pget_mpio_actual_io_mode和H5Pget_mpio_no_collective_cause被用来检查该操作是否已经被集合。
     *
     */
    inline void
    release_plist(hid_t &                    plist,
                  H5D_mpio_actual_io_mode_t &io_mode,
                  uint32_t &                 local_no_collective_cause,
                  uint32_t &                 global_no_collective_cause,
                  const bool                 mpi,
                  const bool                 query_io_mode);

    /**
     * 将HDF5 no_collective_cause代码转换成人类可读的字符串。
     *
     */
    inline std::string
    no_collective_cause_to_string(const uint32_t no_collective_cause);
  } // namespace internal



  // definitions

  namespace internal
  {
    template <typename number>
    std::shared_ptr<hid_t>
    get_hdf5_datatype()
    {
      static_assert(std::is_same<number, float>::value ||
                      std::is_same<number, double>::value ||
                      std::is_same<number, int>::value ||
                      std::is_same<number, bool>::value ||
                      std::is_same<number, unsigned int>::value ||
                      std::is_same<number, std::complex<float>>::value ||
                      std::is_same<number, std::complex<double>>::value,
                    "The data type you are trying to get the HDF5 tag for "
                    "is not supported by this function.");

      if (std::is_same<number, float>::value)
        {
          return std::make_shared<hid_t>(H5T_NATIVE_FLOAT);
        }
      else if (std::is_same<number, double>::value)
        {
          return std::make_shared<hid_t>(H5T_NATIVE_DOUBLE);
        }
      else if (std::is_same<number, int>::value)
        {
          return std::make_shared<hid_t>(H5T_NATIVE_INT);
        }
      else if (std::is_same<number, unsigned int>::value)
        {
          return std::make_shared<hid_t>(H5T_NATIVE_UINT);
        }
      else if (std::is_same<number, bool>::value)
        {
          return std::make_shared<hid_t>(H5T_NATIVE_HBOOL);
        }
      else if (std::is_same<number, std::complex<float>>::value)
        {
          std::shared_ptr<hid_t> t_type =
            std::shared_ptr<hid_t>(new hid_t, [](hid_t *pointer) {
              // Release the HDF5 resource
              const herr_t ret = H5Tclose(*pointer);
              AssertNothrow(ret >= 0, ExcInternalError());
              (void)ret;
              delete pointer;
            });

          *t_type = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<float>));
          //  The C++ standards committee agreed to mandate that the storage
          //  format used for the std::complex type be binary-compatible with
          //  the C99 type, i.e. an array T[2] with consecutive real [0] and
          //  imaginary [1] parts.
          herr_t ret = H5Tinsert(*t_type, "r", 0, H5T_NATIVE_FLOAT);
          Assert(ret >= 0, ExcInternalError());
          ret = H5Tinsert(*t_type, "i", sizeof(float), H5T_NATIVE_FLOAT);
          Assert(ret >= 0, ExcInternalError());
          (void)ret;
          return t_type;
        }
      else if (std::is_same<number, std::complex<double>>::value)
        {
          std::shared_ptr<hid_t> t_type =
            std::shared_ptr<hid_t>(new hid_t, [](hid_t *pointer) {
              // Release the HDF5 resource
              const herr_t ret = H5Tclose(*pointer);
              AssertNothrow(ret >= 0, ExcInternalError());
              (void)ret;
              delete pointer;
            });
          *t_type = H5Tcreate(H5T_COMPOUND, sizeof(std::complex<double>));
          //  The C++ standards committee agreed to mandate that the storage
          //  format used for the std::complex type be binary-compatible with
          //  the C99 type, i.e. an array T[2] with consecutive real [0] and
          //  imaginary [1] parts.
          herr_t ret = H5Tinsert(*t_type, "r", 0, H5T_NATIVE_DOUBLE);
          Assert(ret >= 0, ExcInternalError());
          ret = H5Tinsert(*t_type, "i", sizeof(double), H5T_NATIVE_DOUBLE);
          Assert(ret >= 0, ExcInternalError());
          (void)ret;
          return t_type;
        }

      // The function should not reach this point
      Assert(false, ExcInternalError());
      return {};
    }



    template <typename number>
    std::vector<hsize_t>
    get_container_dimensions(const std::vector<number> &data)
    {
      std::vector<hsize_t> dimensions = {data.size()};
      return dimensions;
    }



    template <typename number>
    std::vector<hsize_t>
    get_container_dimensions(const Vector<number> &data)
    {
      std::vector<hsize_t> dimensions = {data.size()};
      return dimensions;
    }



    template <typename number>
    std::vector<hsize_t>
    get_container_dimensions(const FullMatrix<number> &data)
    {
      std::vector<hsize_t> dimensions = {data.m(), data.n()};
      return dimensions;
    }



    template <typename number>
    unsigned int
    get_container_size(const std::vector<number> &data)
    {
      return static_cast<unsigned int>(data.size());
    }



    template <typename number>
    unsigned int
    get_container_size(const Vector<number> &data)
    {
      return static_cast<unsigned int>(data.size());
    }



    template <typename number>
    unsigned int
    get_container_size(const FullMatrix<number> &data)
    {
      return static_cast<unsigned int>(data.m() * data.n());
    }



    template <typename Container>
    typename std::enable_if<
      std::is_same<Container,
                   std::vector<typename Container::value_type>>::value,
      Container>::type
    initialize_container(const std::vector<hsize_t> &dimensions)
    {
      return Container(std::accumulate(
        dimensions.begin(), dimensions.end(), 1, std::multiplies<int>()));
    }



    template <typename Container>
    typename std::enable_if<
      std::is_same<Container, Vector<typename Container::value_type>>::value,
      Container>::type
    initialize_container(const std::vector<hsize_t> &dimensions)
    {
      return Container(std::accumulate(
        dimensions.begin(), dimensions.end(), 1, std::multiplies<int>()));
    }



    template <typename Container>
    typename std::enable_if<
      std::is_same<Container,
                   FullMatrix<typename Container::value_type>>::value,
      Container>::type
    initialize_container(const std::vector<hsize_t> &dimensions)
    {
      // If the rank is higher than 2, then remove single-dimensional entries
      // from the shape defined by dimensions. This is equivalent to the squeeze
      // function of python/numpy. For example the following code would convert
      // the vector {1,3,1,2} to {3,2}
      std::vector<hsize_t> squeezed_dimensions;

      if (dimensions.size() > 2)
        {
          for (const auto &dimension : dimensions)
            {
              if (dimension > 1)
                squeezed_dimensions.push_back(dimension);
            }
        }
      else
        {
          squeezed_dimensions = dimensions;
        }

      AssertDimension(squeezed_dimensions.size(), 2);
      return Container(squeezed_dimensions[0], squeezed_dimensions[1]);
    }


    inline void
    set_plist(hid_t &plist, const bool mpi)
    {
      if (mpi)
        {
#  ifdef DEAL_II_WITH_MPI
          plist = H5Pcreate(H5P_DATASET_XFER);
          Assert(plist >= 0, ExcInternalError());
          const herr_t ret = H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
          (void)ret;
          Assert(ret >= 0, ExcInternalError());
#  else
          AssertThrow(false, ExcNotImplemented());
#  endif
        }
      else
        {
          plist = H5P_DEFAULT;
        }

      (void)plist;
      (void)mpi;
    }


    inline void
    release_plist(hid_t &                    plist,
                  H5D_mpio_actual_io_mode_t &io_mode,
                  uint32_t &                 local_no_collective_cause,
                  uint32_t &                 global_no_collective_cause,
                  const bool                 mpi,
                  const bool                 query_io_mode)
    {
      if (mpi)
        {
#  ifdef DEAL_II_WITH_MPI
          herr_t ret;
          (void)ret;
          if (query_io_mode)
            {
              ret = H5Pget_mpio_actual_io_mode(plist, &io_mode);
              Assert(ret >= 0, ExcInternalError());
              ret =
                H5Pget_mpio_no_collective_cause(plist,
                                                &local_no_collective_cause,
                                                &global_no_collective_cause);
              Assert(ret >= 0, ExcInternalError());
            }
          ret = H5Pclose(plist);
          Assert(ret >= 0, ExcInternalError());
#  else
          AssertThrow(false, ExcNotImplemented());
#  endif
        }

      (void)plist;
      (void)io_mode;
      (void)local_no_collective_cause;
      (void)global_no_collective_cause;
      (void)mpi;
      (void)query_io_mode;
    }


    inline std::string
    no_collective_cause_to_string(const uint32_t no_collective_cause)
    {
      std::string message;

      auto append_to_message = [&message](const char *p) {
        if (message.length() > 0)
          message += ", ";
        message += p;
      };

      // The first is not a bitmask comparison, the rest are bitmask
      // comparisons.
      // https://support.hdfgroup.org/HDF5/doc/RM/RM_H5P.html#Property-GetMpioNoCollectiveCause
      // See H5Ppublic.h
      // Hex codes are used because the HDF5 Group can deprecate some of the
      // enum codes. For example the enum code H5D_MPIO_FILTERS is not defined
      // in 1.10.2 because it is possible to use compressed datasets with the
      // MPI/IO driver.

      // H5D_MPIO_COLLECTIVE
      if (no_collective_cause == 0x00)
        {
          append_to_message("H5D_MPIO_COLLECTIVE");
        }
      // H5D_MPIO_SET_INDEPENDENT
      if ((no_collective_cause & 0x01) == 0x01)
        {
          append_to_message("H5D_MPIO_SET_INDEPENDENT");
        }
      // H5D_MPIO_DATATYPE_CONVERSION
      if ((no_collective_cause & 0x02) == 0x02)
        {
          append_to_message("H5D_MPIO_DATATYPE_CONVERSION");
        }
      // H5D_MPIO_DATA_TRANSFORMS
      if ((no_collective_cause & 0x04) == 0x04)
        {
          append_to_message("H5D_MPIO_DATA_TRANSFORMS");
        }
      // H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES
      if ((no_collective_cause & 0x10) == 0x10)
        {
          append_to_message("H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES");
        }
      // H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET
      if ((no_collective_cause & 0x20) == 0x20)
        {
          append_to_message("H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET");
        }
      // H5D_MPIO_FILTERS
      if ((no_collective_cause & 0x40) == 0x40)
        {
          append_to_message("H5D_MPIO_FILTERS");
        }
      return message;
    }
  } // namespace internal


  template <typename T>
  T
  HDF5Object::get_attribute(const std::string &attr_name) const
  {
    const std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<T>();
    T                            value;
    hid_t                        attr;
    herr_t                       ret;

    attr = H5Aopen(*hdf5_reference, attr_name.data(), H5P_DEFAULT);
    Assert(attr >= 0, ExcMessage("Error at H5Aopen"));
    (void)ret;
    ret = H5Aread(attr, *t_type, &value);
    Assert(ret >= 0, ExcMessage("Error at H5Aread"));
    (void)ret;
    ret = H5Aclose(attr);
    Assert(ret >= 0, ExcMessage("Error at H5Aclose"));

    return value;
  }



  template <>
  inline std::string
  HDF5Object::get_attribute(const std::string &attr_name) const
  {
    // Reads a UTF8 variable string
    //
    // code inspired from
    // https://support.hdfgroup.org/ftp/HDF5/examples/misc-examples/vlstratt.c
    //
    // In the case of a variable length string the user does not have to reserve
    // memory for string_out. H5Aread will reserve the memory and the
    // user has to free the memory.
    //
    // Todo:
    // - Use H5Dvlen_reclaim instead of free

    char * string_out;
    hid_t  attr;
    hid_t  type;
    herr_t ret;

     /* Create a datatype to refer to. */ 
    type = H5Tcopy(H5T_C_S1);
    Assert(type >= 0, ExcInternalError());

    // Python strings are encoded in UTF8
    ret = H5Tset_cset(type, H5T_CSET_UTF8);
    Assert(type >= 0, ExcInternalError());

    ret = H5Tset_size(type, H5T_VARIABLE);
    Assert(ret >= 0, ExcInternalError());

    attr = H5Aopen(*hdf5_reference, attr_name.data(), H5P_DEFAULT);
    Assert(attr >= 0, ExcInternalError());

    ret = H5Aread(attr, type, &string_out);
    Assert(ret >= 0, ExcInternalError());

    std::string string_value(string_out);
    // The memory of the variable length string has to be freed.
    // H5Dvlen_reclaim could be also used
    free(string_out);
    ret = H5Tclose(type);
    Assert(ret >= 0, ExcInternalError());

    ret = H5Aclose(attr);
    Assert(ret >= 0, ExcInternalError());

    (void)ret;
    return string_value;
  }



  template <typename T>
  void
  HDF5Object::set_attribute(const std::string &attr_name, const T value)
  {
    hid_t  attr;
    hid_t  aid;
    herr_t ret;

    const std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<T>();


    /* 创建标量属性。   
*
*/
    aid = H5Screate(H5S_SCALAR);
    Assert(aid >= 0, ExcMessage("Error at H5Screate"));
    attr = H5Acreate2(*hdf5_reference,
                      attr_name.data(),
                      *t_type,
                      aid,
                      H5P_DEFAULT,
                      H5P_DEFAULT);
    Assert(attr >= 0, ExcMessage("Error at H5Acreate2"));

    /* 编写标量属性。   
*
*/
    ret = H5Awrite(attr, *t_type, &value);
    Assert(ret >= 0, ExcMessage("Error at H5Awrite"));

    ret = H5Sclose(aid);
    Assert(ret >= 0, ExcMessage("Error at H5Sclose"));
    ret = H5Aclose(attr);
    Assert(ret >= 0, ExcMessage("Error at H5Aclose"));

    (void)ret;
  }



  template <>
  inline void
  HDF5Object::set_attribute(const std::string &attr_name,
                            const std::string  value) // NOLINT
  {
    // Writes a UTF8 variable string
    //
    // code inspired from
    // https://support.hdfgroup.org/ftp/HDF5/examples/misc-examples/vlstratt.c

    hid_t  attr;
    hid_t  aid;
    hid_t  t_type;
    herr_t ret;

     /* Create a datatype to refer to. */ 
    t_type = H5Tcopy(H5T_C_S1);
    Assert(t_type >= 0, ExcInternalError());

    // Python strings are encoded in UTF8
    ret = H5Tset_cset(t_type, H5T_CSET_UTF8);
    Assert(t_type >= 0, ExcInternalError());

    ret = H5Tset_size(t_type, H5T_VARIABLE);
    Assert(ret >= 0, ExcInternalError());

    /* 创建标量属性。   
*
*/
    aid = H5Screate(H5S_SCALAR);
    Assert(aid >= 0, ExcMessage("Error at H5Screate"));
    attr = H5Acreate2(
      *hdf5_reference, attr_name.data(), t_type, aid, H5P_DEFAULT, H5P_DEFAULT);
    Assert(attr >= 0, ExcMessage("Error at H5Acreate2"));

    /* 编写标量属性。    在大多数情况下，H5Awrite和H5Dwrite需要一个指向数据的指针。    但是在可变长度的字符串的特殊情况下，H5Awrite取的是字符串的指针的地址。   
*
*/
    const char *c_string_value = value.c_str();
    ret                        = H5Awrite(attr, t_type, &c_string_value);
    Assert(ret >= 0, ExcInternalError());

    ret = H5Sclose(aid);
    Assert(ret >= 0, ExcMessage("Error at H5Sclose"));
    ret = H5Aclose(attr);
    Assert(ret >= 0, ExcMessage("Error at H5Aclose"));

    (void)ret;
  }



  template <typename Container>
  Container
  DataSet::read()
  {
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    hid_t  plist;
    herr_t ret;

    Container data = internal::initialize_container<Container>(dimensions);

    internal::set_plist(plist, mpi);

    ret = H5Dread(*hdf5_reference,
                  *t_type,
                  H5S_ALL,
                  H5S_ALL,
                  plist,
                  make_array_view(data).data());
    Assert(ret >= 0, ExcInternalError());


    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    (void)ret;
    return data;
  }



  template <typename Container>
  Container
  DataSet::read_selection(const std::vector<hsize_t> &coordinates)
  {
    Assert(coordinates.size() % rank == 0,
           ExcMessage(
             "The dimension of coordinates has to be divisible by the rank"));
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    hid_t  plist;
    hid_t  memory_dataspace;
    herr_t ret;

    std::vector<hsize_t> data_dimensions{
      static_cast<hsize_t>(coordinates.size() / rank)};

    Container data = internal::initialize_container<Container>(data_dimensions);

    memory_dataspace = H5Screate_simple(1, data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_elements(*dataspace,
                             H5S_SELECT_SET,
                             data.size(),
                             coordinates.data());
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_elements"));

    internal::set_plist(plist, mpi);

    ret = H5Dread(*hdf5_reference,
                  *t_type,
                  memory_dataspace,
                  *dataspace,
                  plist,
                  make_array_view(data).data());
    Assert(ret >= 0, ExcMessage("Error at H5Dread"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5SClose"));

    (void)ret;
    return data;
  }



  template <typename Container>
  Container
  DataSet::read_hyperslab(const std::vector<hsize_t> &offset,
                          const std::vector<hsize_t> &count)
  {
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    hid_t  plist;
    hid_t  memory_dataspace;
    herr_t ret;

    // In this particular overload of read_hyperslab the data_dimensions are
    // the same as count
    std::vector<hsize_t> data_dimensions = count;

    Container data = internal::initialize_container<Container>(data_dimensions);

    memory_dataspace =
      H5Screate_simple(data_dimensions.size(), data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_hyperslab(*dataspace,
                              H5S_SELECT_SET,
                              offset.data(),
                              nullptr,
                              count.data(),
                              nullptr);
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_hyperslab"));

    internal::set_plist(plist, mpi);

    ret = H5Dread(*hdf5_reference,
                  *t_type,
                  memory_dataspace,
                  *dataspace,
                  plist,
                  make_array_view(data).data());
    Assert(ret >= 0, ExcMessage("Error at H5Dread"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5SClose"));

    (void)ret;
    return data;
  }



  template <typename Container>
  Container
  DataSet::read_hyperslab(const std::vector<hsize_t> &data_dimensions,
                          const std::vector<hsize_t> &offset,
                          const std::vector<hsize_t> &stride,
                          const std::vector<hsize_t> &count,
                          const std::vector<hsize_t> &block)
  {
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    hid_t  plist;
    hid_t  memory_dataspace;
    herr_t ret;

    Container data = internal::initialize_container<Container>(data_dimensions);

    memory_dataspace =
      H5Screate_simple(data_dimensions.size(), data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_hyperslab(*dataspace,
                              H5S_SELECT_SET,
                              offset.data(),
                              stride.data(),
                              count.data(),
                              block.data());
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_hyperslab"));

    internal::set_plist(plist, mpi);

    ret = H5Dread(*hdf5_reference,
                  *t_type,
                  memory_dataspace,
                  *dataspace,
                  plist,
                  make_array_view(data).data());
    Assert(ret >= 0, ExcMessage("Error at H5Dread"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5SClose"));

    (void)ret;
    return data;
  }



  template <typename number>
  void
  DataSet::read_none()
  {
    const std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<number>();
    const std::vector<hsize_t>   data_dimensions = {0};

    hid_t  memory_dataspace;
    hid_t  plist;
    herr_t ret;

    memory_dataspace = H5Screate_simple(1, data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_none(*dataspace);
    Assert(ret >= 0, ExcMessage("H5Sselect_none"));

    internal::set_plist(plist, mpi);

    // The pointer of data can safely be nullptr, see the discussion at the HDF5
    // forum:
    // https://forum.hdfgroup.org/t/parallel-i-o-does-not-support-filters-yet/884/17
    ret = H5Dread(
      *hdf5_reference, *t_type, memory_dataspace, *dataspace, plist, nullptr);
    Assert(ret >= 0, ExcMessage("Error at H5Dread"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5SClose"));

    (void)ret;
  }



  template <typename Container>
  void
  DataSet::write(const Container &data)
  {
    AssertDimension(size, internal::get_container_size(data));
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    hid_t  plist;
    herr_t ret;

    internal::set_plist(plist, mpi);

    ret = H5Dwrite(*hdf5_reference,
                   *t_type,
                   H5S_ALL,
                   H5S_ALL,
                   plist,
                   make_array_view(data).data());
    Assert(ret >= 0, ExcMessage("Error at H5Dwrite"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    (void)ret;
  }



  template <typename Container>
  void
  DataSet::write_selection(const Container &           data,
                           const std::vector<hsize_t> &coordinates)
  {
    AssertDimension(coordinates.size(), data.size() * rank);
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    const std::vector<hsize_t> data_dimensions =
      internal::get_container_dimensions(data);

    hid_t  memory_dataspace;
    hid_t  plist;
    herr_t ret;


    memory_dataspace = H5Screate_simple(1, data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_elements(*dataspace,
                             H5S_SELECT_SET,
                             data.size(),
                             coordinates.data());
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_elements"));

    internal::set_plist(plist, mpi);

    ret = H5Dwrite(*hdf5_reference,
                   *t_type,
                   memory_dataspace,
                   *dataspace,
                   plist,
                   make_array_view(data).data());
    Assert(ret >= 0, ExcMessage("Error at H5Dwrite"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5SClose"));

    (void)ret;
  }



  template <typename Container>
  void
  DataSet::write_hyperslab(const Container &           data,
                           const std::vector<hsize_t> &offset,
                           const std::vector<hsize_t> &count)
  {
    AssertDimension(std::accumulate(count.begin(),
                                    count.end(),
                                    1,
                                    std::multiplies<unsigned int>()),
                    internal::get_container_size(data));
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();
    // In this particular overload of write_hyperslab the data_dimensions are
    // the same as count
    const std::vector<hsize_t> &data_dimensions = count;

    hid_t  memory_dataspace;
    hid_t  plist;
    herr_t ret;

    memory_dataspace =
      H5Screate_simple(data_dimensions.size(), data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_hyperslab(*dataspace,
                              H5S_SELECT_SET,
                              offset.data(),
                              nullptr,
                              count.data(),
                              nullptr);
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_hyperslab"));

    internal::set_plist(plist, mpi);

    ret = H5Dwrite(*hdf5_reference,
                   *t_type,
                   memory_dataspace,
                   *dataspace,
                   plist,
                   make_array_view(data).data());
    Assert(ret >= 0, ExcMessage("Error at H5Dwrite"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5Sclose"));

    (void)ret;
  }



  template <typename Container>
  void
  DataSet::write_hyperslab(const Container &           data,
                           const std::vector<hsize_t> &data_dimensions,
                           const std::vector<hsize_t> &offset,
                           const std::vector<hsize_t> &stride,
                           const std::vector<hsize_t> &count,
                           const std::vector<hsize_t> &block)
  {
    const std::shared_ptr<hid_t> t_type =
      internal::get_hdf5_datatype<typename Container::value_type>();

    hid_t  memory_dataspace;
    hid_t  plist;
    herr_t ret;

    memory_dataspace =
      H5Screate_simple(data_dimensions.size(), data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_hyperslab(*dataspace,
                              H5S_SELECT_SET,
                              offset.data(),
                              stride.data(),
                              count.data(),
                              block.data());
    Assert(ret >= 0, ExcMessage("Error at H5Sselect_hyperslab"));

    internal::set_plist(plist, mpi);

    ret = H5Dwrite(*hdf5_reference,
                   *t_type,
                   memory_dataspace,
                   *dataspace,
                   plist,
                   make_array_view(data).data());
    Assert(ret >= 0, ExcMessage("Error at H5Dwrite"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5Sclose"));

    (void)ret;
  }



  template <typename number>
  void
  DataSet::write_none()
  {
    std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<number>();
    std::vector<hsize_t>   data_dimensions = {0};

    hid_t  memory_dataspace;
    hid_t  plist;
    herr_t ret;

    memory_dataspace = H5Screate_simple(1, data_dimensions.data(), nullptr);
    Assert(memory_dataspace >= 0, ExcMessage("Error at H5Screate_simple"));
    ret = H5Sselect_none(*dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5PSselect_none"));

    internal::set_plist(plist, mpi);

    // The pointer of data can safely be nullptr, see the discussion at the HDF5
    // forum:
    // https://forum.hdfgroup.org/t/parallel-i-o-does-not-support-filters-yet/884/17
    ret = H5Dwrite(
      *hdf5_reference, *t_type, memory_dataspace, *dataspace, plist, nullptr);
    Assert(ret >= 0, ExcMessage("Error at H5Dwrite"));

    internal::release_plist(plist,
                            io_mode,
                            local_no_collective_cause,
                            global_no_collective_cause,
                            mpi,
                            query_io_mode);

    ret = H5Sclose(memory_dataspace);
    Assert(ret >= 0, ExcMessage("Error at H5Sclose"));

    (void)ret;
  }



  template <typename number>
  DataSet
  Group::create_dataset(const std::string &         name,
                        const std::vector<hsize_t> &dimensions) const
  {
    std::shared_ptr<hid_t> t_type = internal::get_hdf5_datatype<number>();
    return {name, *hdf5_reference, dimensions, t_type, mpi};
  }



  template <typename Container>
  void
  Group::write_dataset(const std::string &name, const Container &data) const
  {
    std::vector<hsize_t> dimensions = internal::get_container_dimensions(data);
    auto                 dataset =
      create_dataset<typename Container::value_type>(name, dimensions);
    dataset.write(data);
  }
} // namespace HDF5

DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_WITH_HDF5

#endif // dealii_hdf5_h


