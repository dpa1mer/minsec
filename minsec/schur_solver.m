% Schur complement solver for systems of the form
%   [A B'; B 0] * [x; y] = [u; v]
% where A is positive definite and B is "thin"
function solver = schur_solver(A, B, gpuflag)

solver.A = A;
solver.B = B;

solver.gpuflag = gpuflag;

% factor A
solver.A_chol = decomposition(solver.A, 'chol', 'lower');

if ~isempty(solver.B)
    solver.A_comp = solver.B * (solver.A_chol \ full(solver.B'));
    solver.A_comp_chol = decomposition(solver.A_comp, 'chol', 'lower');

    % if solver.gpuflag
    %     solver.B = gpuArray(solver.B);
    % end
end

% if solver.gpuflag
%     solver.A = gpuArray(solver.A);
% end

solver.solve = @solve;

function [x, y] = solve(u, v, warmstart)
    if nargin < 3
        warmstart = [];
    end

    % Currently, it seems to be faster to solve the system on the CPU
    % because MATLAB does not support sparse pre-factorization on the GPU.
    % So we just pass the data back to the CPU to do the solve.
    if solver.gpuflag
        u = gather(u);
        v = gather(v);
    end

    if ~isempty(solver.B)
        y = solver.A_comp_chol \ (solver.B * solve_A(u, warmstart) - v);
        x = solve_A(u - solver.B' * y, warmstart);
    else
        x = solve_A(u, warmstart);
        y = [];
    end

    if solver.gpuflag
        x = gpuArray(x);
        y = gpuArray(y);
    end
end

function x = solve_A(u, warmstart)
    % if solver.gpuflag
    %     x = pcg(solver.A, u, 1e-6, 2000, [], [], warmstart);
    % else
    x = solver.A_chol \ u;
    % end
end

end