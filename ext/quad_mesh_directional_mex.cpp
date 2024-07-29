#include "mex_helper.h"

#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/IntrinsicVertexTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/curl_matching.h>
#include <directional/power_to_raw.h>
#include <directional/combing.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/branched_isolines.h>
#include <directional/mesh_function_isolines.h>
#include <directional/setup_mesh_function_isolines.h>

// Adapted from Directional tutorial 505

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
            displayError("Second input must be double array.");
        }

        if (inputs[2].getType() != data::ArrayType::STRUCT) {
            displayError("Third input must be structure array.");
        }

        const data::StructArray paramsArray = inputs[2];
        const data::Struct params = paramsArray[0];

        // import the triangle mesh from MATLAB
        const data::StructArray meshDataArray = inputs[0];
        const data::Struct meshData = meshDataArray[0];
        data::TypedArray<double> verts = std::move(meshData["verts"]);
        data::TypedArray<double> faces = std::move(meshData["faces"]);
        data::TypedArray<double> nf = meshData["nf"];

        if (inputs[1].getDimensions()[0] != static_cast<size_t>(nf[0])) {
            displayError("Only works with face-based fields!");
        }

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

        directional::TriMesh meshWhole, meshCut;
        directional::IntrinsicFaceTangentBundle ftb;
        meshWhole.set_mesh(vertsEigen, facesEigen);
        ftb.init(meshWhole);

        // import the field values from MATLAB
        data::TypedArray<double> field = std::move(inputs[1]);
        data::ArrayDimensions fieldDims = field.getDimensions();
        data::buffer_ptr_t<double> fieldPtr = field.release();
        Eigen::Map<const Eigen::MatrixXd> fieldMapped(fieldPtr.get(), fieldDims[0], fieldDims[1]);

        directional::CartesianField rawField;
        rawField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 4);
        rawField.set_extrinsic_field(fieldMapped);
        directional::principal_matching(rawField);

        // Integrate to get a map
        displayMessage("Starting integration...\n");
        directional::CartesianField combedField;
        directional::IntegrationData intData(4);
        directional::setup_integration(rawField, intData, meshCut, combedField);
        
        intData.verbose          = data::TypedArray<bool>(params["verbose"])[0];
        intData.integralSeamless = data::TypedArray<bool>(params["seamless"])[0];
        intData.roundSeams       = data::TypedArray<bool>(params["round_seams"])[0]; // i.e., rather than rounding singularities
        intData.lengthRatio      = data::TypedArray<double>(params["length_ratio"])[0];
        // intData.localInjectivity=true;
    
        Eigen::MatrixXd NFunction, NCornerFunction;
        directional::integrate(combedField,  intData, meshCut, NFunction, NCornerFunction);
        
        // set up polygonal mesh data from integration data
        directional::MeshFunctionIsolinesData mfiData;
        directional::setup_mesh_function_isolines(meshCut, intData, mfiData);
        
        // generate the quad mesh
        displayMessage("Generating mesh...\n");
        Eigen::MatrixXd VPolyMesh;
        Eigen::VectorXi DPolyMesh;
        Eigen::MatrixXi FPolyMesh;
        directional::mesh_function_isolines(meshWhole, mfiData, data::TypedArray<bool>(params["verbose"])[0], VPolyMesh, DPolyMesh, FPolyMesh);
        
        // Convert to 1-indexed
        FPolyMesh.array() += 1;

        displayMessage("Passing mesh back to MATLAB...\n");
        outputs[0] = factory.createArray({static_cast<size_t>(VPolyMesh.rows()), static_cast<size_t>(VPolyMesh.cols())},
                        VPolyMesh.reshaped().cbegin(),
                        VPolyMesh.reshaped().cend());
        outputs[1] = factory.createArray({static_cast<size_t>(DPolyMesh.rows()), 1},
                        DPolyMesh.cbegin(), DPolyMesh.cend());
        outputs[2] = factory.createArray({static_cast<size_t>(FPolyMesh.rows()), static_cast<size_t>(FPolyMesh.cols())},
                        FPolyMesh.reshaped().cbegin(), FPolyMesh.reshaped().cend());
    }
};
