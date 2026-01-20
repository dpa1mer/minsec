function [z, abs_z, err] = mbo_intrinsic(meshData, tauMult, degree)

if nargin < 2
    tauMult = 0.5;
end

if nargin < 3
    degree = 1;
end

nv = meshData.nv;
nf = meshData.nf;

[G_intr, M_intr] = intrinsic_ops(meshData);
V2F = covariant_ops(meshData, degree);

M_cov = V2F' * M_intr * V2F;
G_cov = G_intr * V2F;
Area = kron(spdiags(meshData.areas, 0, nf, nf), speye(2));
L_cov = G_cov' * Area * G_cov;

%% Set up boundary conditions
[~, z_B] = field_bdry_conds(meshData, degree);

Bdry = sparse(1:meshData.nb, meshData.bdryIdx, 1, meshData.nb, meshData.nv);
Int = sparse(1:meshData.ni, meshData.intIdx, 1, meshData.ni, meshData.nv);

%% Compute cross field by MBO

if tauMult <= 0
    % Use "Globally Optimal Direction Fields" method
    A = L_cov;
else
    lambda = eigs(L_cov, M_cov, 2, 'smallestabs');
    tau = tauMult/lambda(2);
    
    A = M_cov + tau * L_cov;
end

A_IB = Int * (A * Bdry');
A_II = Int * (A * Int');
A_II = 0.5 * (A_II + A_II'); % for good measure
cholA_II = decomposition(A_II, 'chol', 'lower');

if tauMult <= 0
    if meshData.nb > 0
        z_I = cholA_II \ (-A_IB * z_B);
        z = Int' * z_I + Bdry' * z_B;
    else
        [z, ~] = eigs(L_cov, M_cov, 1, 'smallestabs');
    end
    abs_z = abs(z);
    z = z ./ abs_z;
else
    z = randn(nv, 1) + 1i * randn(nv, 1);
    z(meshData.bdryIdx) = z_B;
    for j = 1:10000
        z_old = z;
    
        % Diffuse
        z_I = cholA_II \ (Int * (M_cov * z) - A_IB * z_B);
        z = Int' * z_I + Bdry' * z_B;
    
        % Project
        abs_z = abs(z);
        z = z ./ abs_z;
    
        err(j) = norm(z - z_old);
        if err(j) < 1e-4
            break;
        end
    end
end

end