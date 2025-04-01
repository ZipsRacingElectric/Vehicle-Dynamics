% csolve  Solves a custom quadratic program very rapidly.
%
% [vars, status] = csolve(params, settings)
%
% solves the convex optimization problem
%
%   minimize((1/2)*quad_form(x, P) + q'*x)
%   subject to
%     l <= A*x
%     A*x <= u
%
% with variables
%        x   4 x 1
%
% and parameters
%        A   4 x 4
%        P   4 x 4    PSD
%        l   4 x 1
%        q   4 x 1
%        u   4 x 1
%
% Note:
%   - Check status.converged, which will be 1 if optimization succeeded.
%   - You don't have to specify settings if you don't want to.
%   - To hide output, use settings.verbose = 0.
%   - To change iterations, use settings.max_iters = 20.
%   - You may wish to compare with cvxsolve to check the solver is correct.
%
% Specify params.A, ..., params.u, then run
%   [vars, status] = csolve(params, settings)
% Produced by CVXGEN, 2025-03-31 14:19:24 -0400.
% CVXGEN is Copyright (C) 2006-2017 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2017 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: csolve.m.
% Description: Help file for the Matlab solver interface.
