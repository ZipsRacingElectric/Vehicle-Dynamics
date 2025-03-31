clc; clear all;

% Show information about the solver, including required parameters.
help csolve

% Attach random problem data to the params structure.
params.A = randn(4, 4);
params.P = randn(4, 4);
params.l = randn(4, 1);
params.q = randn(4, 1);
params.u = randn(4, 1);

% Exercise the high-speed solver.
[vars, status] = csolve(params);  % solve, saving results.

% Check convergence, and display the optimal variable value.
% Note: we are providing random data, so we can expect some infesiable
% solutions
if ~status.converged, error 'failed to converge'; end
vars.x