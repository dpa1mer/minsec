#include "mex_helper.h"

#include <Eigen/Core>

#include "geometrycentral/utilities/eigen_interop_helpers.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/direction_fields.h"

// Adapted from Geometry Central tutorial

namespace data = matlab::data;

using namespace geometrycentral;
using namespace geometrycentral::surface;

class MexFunction : public MexHelper {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        scoped_redirect_cout mycout_redirect(*this, matlabPtr, factory);

        if (inputs.size() != 2) {
            displayError("Two inputs required!");
        }

        if (inputs[0].getType() != data::ArrayType::STRUCT) {
            displayError("First input must be struct.");
        }

        if (inputs[1].getType() != data::ArrayType::DOUBLE) {
            displayError("Second input must be double.");
        }

        // import the triangle mesh from MATLAB
        const data::StructArray meshDataArray = inputs[0];
        const data::Struct meshData = meshDataArray[0];
        data::TypedArray<double> verts = std::move(meshData["verts"]);
        data::TypedArray<double> faces = std::move(meshData["faces"]);
        const double nb = meshData["nb"][0];

        // Set up mesh
        data::ArrayDimensions inVertsDims = verts.getDimensions();
        data::buffer_ptr_t<double> vertsPtr = verts.release();
        Eigen::Map<const Eigen::MatrixXd> vertsEigen(vertsPtr.get(), inVertsDims[0], inVertsDims[1]);

        data::ArrayDimensions inFacesDims = faces.getDimensions();
        Eigen::MatrixXi facesEigen(inFacesDims[0], inFacesDims[1]);
        for (int j = 0; j < inFacesDims[1]; ++j) {
            for (int i = 0; i < inFacesDims[0]; ++i) {
                facesEigen(i, j) = faces[i][j] - 1;
            }
        }

        std::unique_ptr<ManifoldSurfaceMesh> mesh;
        std::unique_ptr<VertexPositionGeometry> geometry;
        std::tie(mesh, geometry) = makeManifoldSurfaceMeshAndGeometry(vertsEigen, facesEigen);

        const double degree = inputs[1][0];

        // Compute a RoSy field
        displayMessage("Computing RoSy field...\n");
        VertexData<Vector2> powerField;
        if (nb == 0) {
            powerField = computeSmoothestVertexDirectionField(*geometry, degree);
        } else {
            powerField = computeSmoothestBoundaryAlignedVertexDirectionField(*geometry, degree);
        }
        VertexData<Vector2> intrinsicField(*mesh);
        for (Vertex v : mesh->vertices()) {
            intrinsicField[v] = powerField[v].pow(1.0 / degree);
        }

        displayMessage("Converting to extrinsic field...\n");
        geometry->requireVertexTangentBasis();
        Eigen::Matrix<double, 2, Eigen::Dynamic> intrinsicFieldEigen = EigenMap<double, 2>(intrinsicField).transpose();
        Eigen::Matrix<double, 6, Eigen::Dynamic> intrinsicBasis = EigenMap<double, 6>(geometry->vertexTangentBasis).transpose();
        Eigen::Matrix<double, 3, Eigen::Dynamic> extrinsicField(3, intrinsicFieldEigen.cols());
        for (int j = 0; j < intrinsicFieldEigen.cols(); ++j) {
            extrinsicField.col(j) = intrinsicBasis.col(j).reshaped(3, 2) * intrinsicFieldEigen.col(j);
        }

        displayMessage("Passing field back to MATLAB...\n");
        outputs[0] = factory.createArray({static_cast<size_t>(extrinsicField.rows()), static_cast<size_t>(extrinsicField.cols())},
                        extrinsicField.reshaped().cbegin(),
                        extrinsicField.reshaped().cend());
    }
};
