%% generate_C_code.m
% This script loads the symbolic QP matrices P and q from QP_matrices.mat,
% automatically determines all free symbolic parameters used in these matrices,
% and generates a C source file that defines a function to compute P and q.
% The generated C function will have a signature with parameters corresponding to
% all free symbols in P and q.

clc; clear all;

%% Load the QP matrices
load('QP_matrices.mat', 'P', 'q');

%% Automatically extract all free symbols from P and q.
% symvar returns a sorted list of unique symbols in the expression.
paramsP = symvar(P);
paramsQ = symvar(q);
allParams = unique([paramsP, paramsQ]);  % Combine and remove duplicates

% Convert the symbols to strings for use in the C function signature.
paramNames = arrayfun(@char, allParams, 'UniformOutput', false);

%% Create the function signature string.
% The function will have the form:
% void computeQP(double <param1>, double <param2>, ..., double* P, double* q)
% where P and q are output arrays.
paramStr = '';
for k = 1:length(paramNames)
    % Append each parameter in the format "double param,".
    paramStr = [paramStr, sprintf('double %s, ', paramNames{k})];
end
% Append the output arrays.
paramStr = [paramStr, 'double* P, double* q'];

%% Get dimensions for P and q.
[nP_rows, nP_cols] = size(P);
nq = numel(q);  % assuming q is a column vector

%% Generate C code for computing P.
c_code_str = '';
c_code_str = [c_code_str, sprintf('  /* Compute P matrix elements */%s', newline)];
for i = 1:nP_rows
    for j = 1:nP_cols
        expr = simplify(P(i,j));
        % Use ccode to generate C code for this expression.
        tmp = ccode(expr);
        % Extract the expression after the '='. ccode returns something like "t0 = <expr>;"
        eqIdx = strfind(tmp, '=');
        if ~isempty(eqIdx)
            c_expr = strtrim(tmp(eqIdx+1:end));
            % Remove trailing semicolon if present.
            if c_expr(end)==';'
                c_expr(end) = [];
            end
        else
            c_expr = tmp;
        end
        % Compute a linear index for the element (assuming row-major order).
        idx = (i-1)*nP_cols + (j-1);
        line_str = sprintf('  P[%d] = %s;\n', idx, c_expr);
        c_code_str = [c_code_str, line_str];
    end
end

%% Generate C code for computing q.
c_code_str = [c_code_str, sprintf('  /* Compute q vector elements */%s', newline)];
for i = 1:nq
    expr = simplify(q(i));
    tmp = ccode(expr);
    eqIdx = strfind(tmp, '=');
    if ~isempty(eqIdx)
        c_expr = strtrim(tmp(eqIdx+1:end));
        if c_expr(end)==';'
            c_expr(end) = [];
        end
    else
        c_expr = tmp;
    end
    line_str = sprintf('  q[%d] = %s;\n', i-1, c_expr);
    c_code_str = [c_code_str, line_str];
end

%% Assemble the complete C function code.
c_function_code = sprintf(['/* Generated C code for computing QP matrices\n',...
    '   Free parameters: %s\n',...
    '   Outputs: P (matrix) and q (vector) stored in arrays\n*/\n',...
    '#include <math.h>\n\n',...
    'void computeQP(%s) {\n',...
    '%s',...
    '}\n'], strjoin(paramNames, ', '), paramStr, c_code_str);

%% Write the generated C code to a file.
fid = fopen('qp_matrix.c', 'w');
if fid == -1
    error('Could not open file for writing.');
end
fprintf(fid, '%s', c_function_code);
fclose(fid);

disp('C code has been generated in qp_matrix.c');
