// B.10. The Stokes Equation Using Taylor–Hood
// Elements
// File: functions-stokes.cc
// Equation: Stokes equation
// Discretization: Taylor–Hood finite elements
// Grid: Uniform quadrilateral grid, implemented with YaspGrid
// Solver: GMRes
// Distributed: no
// Discussed in: Chapter 10.8

#include <config.h>

#include <array>
#include <vector>

#include <dune/common/indices.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/taylorhoodbasis.hh>
#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>

using namespace Dune;

// Compute the stiffness matrix for a single element
// { local_assembler_signature_begin }
template <class LocalView>
void assembleElementStiffnessMatrix(const LocalView &localView,
                                    Matrix<double> &elementMatrix)
// { local_assembler_signature_end }
{
    // Get the grid element from the local FE basis view
    // { local_assembler_get_element_information_begin }
    using Element = typename LocalView::Element;
    const Element element = localView.element();

    constexpr int dim = Element::dimension;
    auto geometry = element.geometry();
    // { local_assembler_get_element_information_end }

    // Set all matrix entries to zero
    // { initialize_element_matrix_begin }
    elementMatrix.setSize(localView.size(), localView.size());
    elementMatrix = 0; // Fills the entire matrix with zeros
    // { initialize_element_matrix_end }

    // Get set of shape functions for this element
    // { get_local_fe_begin }
    using namespace Indices;
    const auto &velocityLocalFiniteElement = localView.tree().child(_0, 0).finiteElement();
    const auto &pressureLocalFiniteElement = localView.tree().child(_1).finiteElement();
    // { get_local_fe_end }

    // Get a quadrature rule
    // { begin_quad_loop_begin }
    int order = 2 * (dim * velocityLocalFiniteElement.localBasis().order() - 1);
    const auto &quad = QuadratureRules<double, dim>::rule(element.type(), order);

    // Loop over all quadrature points
    for (const auto &quadPoint : quad)
    {
        // { begin_quad_loop_end }
        // { quad_loop_preamble_begin }
        // The transposed inverse Jacobian of the map from the
        // reference element to the element
        const auto jacobianInverseTransposed = geometry.jacobianInverseTransposed(quadPoint.position());

        // The multiplicative factor in the integral transformation formula
        const auto integrationElement = geometry.integrationElement(quadPoint.position());
        // { quad_loop_preamble_end }

        ///////////////////////////////////////////////////////////////////////
        // Velocity--velocity coupling
        ///////////////////////////////////////////////////////////////////////

        // The gradients of the shape functions on the reference element
        // { velocity_gradients_begin }
        std::vector<FieldMatrix<double, 1, dim>> referenceGradients;
        velocityLocalFiniteElement.localBasis().evaluateJacobian(
            quadPoint.position(),
            referenceGradients);

        // Compute the shape function gradients on the grid element
        std::vector<FieldVector<double, dim>> gradients(referenceGradients.size());
        for (size_t i = 0; i < gradients.size(); i++)
            jacobianInverseTransposed.mv(referenceGradients[i][0], gradients[i]);
        // { velocity_gradients_end }

        // Compute the actual matrix entries
        // { velocity_velocity_coupling_begin }
        for (size_t i = 0; i < velocityLocalFiniteElement.size(); i++)
            for (size_t j = 0; j < velocityLocalFiniteElement.size(); j++)
                for (size_t k = 0; k < dim; k++)
                {
                    size_t row = localView.tree().child(_0, k).localIndex(i);
                    size_t col = localView.tree().child(_0, k).localIndex(j);
                    elementMatrix[row][col] += (gradients[i] * gradients[j]) * quadPoint.weight() * integrationElement;
                }
        // { velocity_velocity_coupling_end }

        ///////////////////////////////////////////////////////////////////////
        // Velocity--pressure coupling
        ///////////////////////////////////////////////////////////////////////

        // The values of the pressure shape functions
        // { pressure_values_begin }
        std::vector<FieldVector<double, 1>> pressureValues;
        pressureLocalFiniteElement
            .localBasis()
            .evaluateFunction(quadPoint.position(), pressureValues);
        // { pressure_values_end }

        // Compute the actual matrix entries
        // { velocity_pressure_coupling_begin }
        for (size_t i = 0; i < velocityLocalFiniteElement.size(); i++)
            for (size_t j = 0; j < pressureLocalFiniteElement.size(); j++)
                for (size_t k = 0; k < dim; k++)
                {
                    size_t vIndex = localView.tree().child(_0, k).localIndex(i);
                    size_t pIndex = localView.tree().child(_1).localIndex(j);

                    auto value = gradients[i][k] * pressureValues[j] * quadPoint.weight() * integrationElement;
                    elementMatrix[vIndex][pIndex] += value;
                    elementMatrix[pIndex][vIndex] += value;
                }
        // { velocity_pressure_coupling_end }
    }
}

// Set the occupation pattern of the stiffness matrix
template <class Basis, class Matrix>
void setOccupationPattern(const Basis &basis, Matrix &matrix)
{
    enum
    {
        dim = Basis::GridView::dimension
    };

    // MatrixIndexSets store the occupation pattern of a sparse matrix.
    // They are not particularly efficient, but simple to use.
    std::array<std::array<MatrixIndexSet, 2>, 2> nb;

    // Set sizes of the 2x2 submatrices
    for (size_t i = 0; i < 2; i++)
        for (size_t j = 0; j < 2; j++)
            nb[i][j].resize(basis.size({i}), basis.size({j}));

    // A view on the FE basis on a single element
    auto localView = basis.localView();

    // Loop over all leaf elements
    for (const auto &element : elements(basis.gridView()))
    {
        // Bind the local view to the current element
        localView.bind(element);

        // Add element stiffness matrix onto the global stiffness matrix
        for (size_t i = 0; i < localView.size(); i++)
        {
            // Global index of the i-th local degree of freedom of the current element
            auto row = localView.index(i);

            for (size_t j = 0; j < localView.size(); j++)
            {
                // Global index of the j-th local degree of freedom of the current element
                auto col = localView.index(j);

                nb[row[0]][col[0]].add(row[1], col[1]);
            }
        }
    }

    // Give the matrix the occupation pattern we want.
    using namespace Indices;
    nb[0][0].exportIdx(matrix[_0][_0]);
    nb[0][1].exportIdx(matrix[_0][_1]);
    nb[1][0].exportIdx(matrix[_1][_0]);
    nb[1][1].exportIdx(matrix[_1][_1]);
}

// { matrixentry_begin }
template <class Matrix, class MultiIndex>
decltype(auto) matrixEntry(
    Matrix &matrix, const MultiIndex &row, const MultiIndex &col)
{
    using namespace Indices;
    if ((row[0] == 0) && (col[0] == 0))
        return matrix[_0][_0][row[1]][col[1]][row[2]][col[2]];
    if ((row[0] == 0) && (col[0] == 1))
        return matrix[_0][_1][row[1]][col[1]][row[2]][0];
    if ((row[0] == 1) && (col[0] == 0))
        return matrix[_1][_0][row[1]][col[1]][0][col[2]];
    return matrix[_1][_1][row[1]][col[1]];
}
// { matrixentry_end }

// Assemble the Laplace stiffness matrix on the given grid view
// { global_assembler_signature_begin }
template <class Basis, class Matrix>
void assembleStokesMatrix(const Basis &basis, Matrix &matrix)
// { global_assembler_signature_end }
{
    // { setup_matrix_pattern_begin }
    // Set matrix size and occupation pattern
    setOccupationPattern(basis, matrix);

    // Set all entries to zero
    matrix = 0;
    // { setup_matrix_pattern_end }

    // A view on the FE basis on a single element
    // { get_localview_begin }
    auto localView = basis.localView();
    // { get_localview_end }

    // A loop over all elements of the grid
    // { element_loop_and_bind_begin }
    for (const auto &element : elements(basis.gridView()))
    {
        // Bind the local FE basis view to the current element
        localView.bind(element);
        // { element_loop_and_bind_end }

        // Now let’s get the element stiffness matrix
        // A dense matrix is used for the element stiffness matrix
        // { setup_element_stiffness_begin }
        Dune::Matrix<double> elementMatrix;
        assembleElementStiffnessMatrix(localView, elementMatrix);
        // { setup_element_stiffness_end }

        // Add element stiffness matrix onto the global stiffness matrix
        // { accumulate_global_matrix_begin }
        for (size_t i = 0; i < elementMatrix.N(); i++)
        {
            // The global index of the i-th local degree of freedom
            // of the current element
            auto row = localView.index(i);

            for (size_t j = 0; j < elementMatrix.M(); j++)
            {
                // The global index of the j-th local degree of freedom
                // of the current element
                auto col = localView.index(j);
                matrixEntry(matrix, row, col) += elementMatrix[i][j];
            }
        }
        // { accumulate_global_matrix_end }
    }
}

// { main_begin }
int main(int argc, char *argv[])
{
    // Set up MPI, if available
    MPIHelper::instance(argc, argv);
    // { mpi_setup_end }

    ///////////////////////////////////
    // Generate the grid
    ///////////////////////////////////

    // { grid_setup_begin }
    constexpr int dim = 2;
    using Grid = YaspGrid<dim>;
    FieldVector<double, dim> upperRight = {1, 1};
    std::array<int, dim> nElements = {4, 4};
    Grid grid(upperRight, nElements);

    using GridView = typename Grid::LeafGridView;
    GridView gridView = grid.leafGridView();
    // { grid_setup_end }

    /////////////////////////////////////////////////////////
    // Choose a finite element space
    /////////////////////////////////////////////////////////

    // { function_space_basis_begin }
    using namespace Functions::BasisFactory;

    constexpr std::size_t p = 1; // Pressure order for Taylor-Hood

    auto taylorHoodBasis = makeBasis(
        gridView,
        composite(
            power<dim>(
                lagrange<p + 1>(),
                blockedInterleaved()),
            lagrange<p>()));
    // { function_space_basis_end }

    /////////////////////////////////////////////////////////
    // Stiffness matrix and right hand side vector
    /////////////////////////////////////////////////////////

    // { linear_algebra_setup_begin }
    using VelocityVector = BlockVector<FieldVector<double, dim>>;
    using PressureVector = BlockVector<double>;
    using Vector = MultiTypeBlockVector<VelocityVector, PressureVector>;

    using Matrix00 = BCRSMatrix<FieldMatrix<double, dim, dim>>;
    using Matrix01 = BCRSMatrix<FieldMatrix<double, dim, 1>>;
    using Matrix10 = BCRSMatrix<FieldMatrix<double, 1, dim>>;
    using Matrix11 = BCRSMatrix<double>;
    using MatrixRow0 = MultiTypeBlockVector<Matrix00, Matrix01>;
    using MatrixRow1 = MultiTypeBlockVector<Matrix10, Matrix11>;
    using Matrix = MultiTypeBlockMatrix<MatrixRow0, MatrixRow1>;
    // { linear_algebra_setup_end }

    /////////////////////////////////////////////////////////
    // Assemble the system
    /////////////////////////////////////////////////////////

    // { rhs_assembly_begin }
    Vector rhs;

    auto rhsBackend = Functions::istlVectorBackend(rhs);

    rhsBackend.resize(taylorHoodBasis);
    rhs = 0;
    // { rhs_assembly_end }

    // { matrix_assembly_begin }
    Matrix stiffnessMatrix;
    assembleStokesMatrix(taylorHoodBasis, stiffnessMatrix);
    // { matrix_assembly_end }

    /////////////////////////////////////////////////////////
    // Set Dirichlet values.
    // Only velocity components have Dirichlet boundary values
    /////////////////////////////////////////////////////////

    // { initialize_boundary_dofs_vector_begin }
    using VelocityBitVector = std::vector<std::array<char, dim>>;
    using PressureBitVector = std::vector<char>;
    using BitVector = TupleVector<VelocityBitVector, PressureBitVector>;

    BitVector isBoundary;

    auto isBoundaryBackend = Functions::istlVectorBackend(isBoundary);
    isBoundaryBackend.resize(taylorHoodBasis);

    using namespace Indices;
    for (auto &&b0i : isBoundary[_0])
        for (std::size_t j = 0; j < b0i.size(); ++j)
            b0i[j] = false;
    std::fill(isBoundary[_1].begin(), isBoundary[_1].end(), false);
    // { initialize_boundary_dofs_vector_end }

    // { determine_boundary_dofs_begin }
    Functions::forEachBoundaryDOF(
        Functions::subspaceBasis(taylorHoodBasis, _0),
        [&](auto &&index)
        {
            isBoundaryBackend[index] = true;
        });
    // { determine_boundary_dofs_end }

    // { interpolate_dirichlet_values_begin }
    using Coordinate = GridView::Codim<0>::Geometry::GlobalCoordinate;
    using VelocityRange = FieldVector<double, dim>;
    auto &&g = [](Coordinate x)
    {
        return VelocityRange{0.0, (x[0] < 1e-8) ? 1.0 : 0.0};
    };

    Functions::interpolate(Functions::subspaceBasis(taylorHoodBasis, _0),
                           rhs,
                           g,
                           isBoundary);
    // { interpolate_dirichlet_values_end }

    ////////////////////////////////////////////
    // Modify Dirichlet rows
    ////////////////////////////////////////////

    // Loop over the matrix rows
    // { set_dirichlet_matrix_begin }
    auto localView = taylorHoodBasis.localView();
    for (const auto &element : elements(gridView))
    {
        localView.bind(element);
        for (size_t i = 0; i < localView.size(); ++i)
        {
            auto row = localView.index(i);
            // If row corresponds to a boundary entry,
            // modify it to be an identity matrix row.
            if (isBoundaryBackend[row])
                for (size_t j = 0; j < localView.size(); ++j)
                {
                    auto col = localView.index(j);
                    matrixEntry(stiffnessMatrix, row, col) = (i == j) ? 1 : 0;
                }
        }
    }
    // { set_dirichlet_matrix_end }

    ////////////////////////////
    // Compute solution
    ////////////////////////////
    // { stokes_solve_begin }
    // Initial iterate: Start from the rhs vector,
    // that way the Dirichlet entries are already correct.
    Vector x = rhs;

    // Turn the matrix into a linear operator
    MatrixAdapter<Matrix, Vector, Vector> stiffnessOperator(stiffnessMatrix);

    // Fancy (but only) way to not have a preconditioner at all
    Richardson<Vector, Vector> preconditioner(1.0);

    // Construct the iterative solver
    RestartedGMResSolver<Vector> solver(
        stiffnessOperator, // Operator to invert
        preconditioner,    // Preconditioner
        1e-10,             // Desired residual reduction factor
        500,               // Number of iterations between restarts,
        // here: no restarting
        500, // Maximum number of iterations
        2);  // Verbosity of the solver

    // Object storing some statistics about the solving process
    InverseOperatorResult statistics;

    // Solve!
    solver.apply(x, rhs, statistics);
    // { stokes_solve_end }

    ////////////////////////////////////////////////////////////////////////////
    // Make a discrete function from the FE basis and the coefficient vector
    ////////////////////////////////////////////////////////////////////////////

    // { make_result_functions_begin }
    using VelocityRange = FieldVector<double, dim>;
    using PressureRange = double;

    auto velocityFunction = Functions::makeDiscreteGlobalBasisFunction<VelocityRange>(
        Functions::subspaceBasis(taylorHoodBasis, _0), x);
    auto pressureFunction = Functions::makeDiscreteGlobalBasisFunction<PressureRange>(
        Functions::subspaceBasis(taylorHoodBasis, _1), x);
    // { make_result_functions_end }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Write result to VTK file
    // We need to subsample, because the dune-grid VTKWriter cannot natively display
    // second-order functions
    //////////////////////////////////////////////////////////////////////////////////////////////
    // { vtk_output_begin }
    SubsamplingVTKWriter<GridView> vtkWriter(
        gridView,
        refinementLevels(2));
    vtkWriter.addVertexData(
        velocityFunction,
        VTK::FieldInfo("velocity", VTK::FieldInfo::Type::vector, dim));
    vtkWriter.addVertexData(
        pressureFunction,
        VTK::FieldInfo("pressure", VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.write("stokes-taylorhood-result");
    // { vtk_output_end }
}