#include "mex_helper.h"

#include <directional/TriMesh.h>
#include <directional/IntrinsicVertexTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>

// Adapted from Directional tutorial 301

namespace data = matlab::data;

class MexFunction : public MexHelper {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        scoped_redirect_cout mycout_redirect(*this, matlabPtr, factory);

        if (inputs.size() != 3) {
            displayError("Three inputs required!");
        }

        if (inputs[0].getType() != data::ArrayType::STRUCT) {
            displayError("First input must be struct.");
        }

        if (inputs[1].getType() != data::ArrayType::DOUBLE) {
            displayError("Second input must be double.");
        }

        if (inputs[2].getType() != data::ArrayType::DOUBLE) {
            displayError("Third input must be double array.");
        }

        // import the triangle mesh from MATLAB
        const data::StructArray meshDataArray = inputs[0];
        const data::Struct meshData = meshDataArray[0];
        data::TypedArray<double> verts = std::move(meshData["verts"]);
        data::TypedArray<double> faces = std::move(meshData["faces"]);
        data::TypedArray<double> bdryIdx = std::move(meshData["bdryIdx"]);

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

        data::ArrayDimensions bdryIdxDims = bdryIdx.getDimensions();
        if (bdryIdxDims[1] != 1) {
            displayError("bdryIdx must be an index vector.");
        }
        Eigen::VectorXi bdryIdxEigen(bdryIdxDims[0], bdryIdxDims[1]);
        for (int i = 0; i < bdryIdxDims[0]; ++i) {
            bdryIdxEigen(i) = bdryIdx[i] - 1;
        }

        directional::TriMesh mesh;
        directional::IntrinsicVertexTangentBundle vtb;
        mesh.set_mesh(vertsEigen, facesEigen);
        vtb.init(mesh);

        const double degree = inputs[1][0];

        // import the boundary values from MATLAB
        data::TypedArray<double> bdryVals = std::move(inputs[2]);
        data::ArrayDimensions bdryValsDims = bdryVals.getDimensions();
        data::buffer_ptr_t<double> bdryValsPtr = bdryVals.release();
        Eigen::Map<const Eigen::MatrixXd> bdryValsMapped(bdryValsPtr.get(), bdryValsDims[0], bdryValsDims[1]);

        // Compute a RoSy field
        displayMessage("Computing RoSy field...\n");
        directional::CartesianField rosyField;
        rosyField.init(vtb, directional::fieldTypeEnum::POWER_FIELD, degree);
        directional::power_field(vtb, bdryIdxEigen, bdryValsMapped, Eigen::VectorXd::Constant(bdryIdxEigen.size(), -1.0), degree, rosyField);

        displayMessage("Converting to raw field...\n");
        directional::CartesianField rawField;
        directional::power_to_raw(rosyField, degree, rawField, false);

        displayMessage("Passing field back to MATLAB...\n");
        outputs[0] = factory.createArray({static_cast<size_t>(rawField.extField.rows()), static_cast<size_t>(rawField.extField.cols())},
                        rawField.extField.reshaped().cbegin(),
                        rawField.extField.reshaped().cend());
    }
};
