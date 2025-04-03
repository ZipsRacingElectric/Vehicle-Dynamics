% Script to take the D2704 RBF Model and generate a C file which declares
% these as lookup tables in flash memory

clc; clear all;

load('../tire_data/d2704/d2704_7in_rbf_full_model.mat');


%% Reduce Tire Model Lookup Table Data
% --- Helper function to reduce symmetrical breakpoint array while keeping zero ---
reduce_sym_array = @(arr) arr( ...
    ismember(arr, 0) | ...               % keep zero
    mod(1:length(arr), 2) == 1);         % pick every other element (odd indices)

% Reduce sa_vals and fz_vals while preserving zero
reduced_sa_vals = reduce_sym_array(sa_vals);
reduced_fz_vals = reduce_sym_array(fz_vals);

% Find the indices of the new reduced values
[~, sa_idx] = ismember(reduced_sa_vals, sa_vals);
[~, fz_idx] = ismember(reduced_fz_vals, fz_vals);

% Truncate and cast grid data
Fx_reduced = int16(fix(Fx_grid(:, sa_idx, fz_idx)));
Fy_reduced = int16(fix(Fy_grid(:, sa_idx, fz_idx)));

% Display info
fprintf('Original sa_vals: %d → Reduced: %d\n', numel(sa_vals), numel(reduced_sa_vals));
fprintf('Original fz_vals: %d → Reduced: %d\n', numel(fz_vals), numel(reduced_fz_vals));

%% Create C data file
% Open C file
fid = fopen('../tire_data/d2704/tire_data.c', 'w');

% Write file header
fprintf(fid, '// Auto-generated lookup table\n');
fprintf(fid, '#include <stdint.h>\n\n');

% Write breakpoint arrays
write_float_array(fid, 'sl_breakpoints', sl_vals);
write_float_array(fid, 'sa_breakpoints', reduced_sa_vals);
write_float_array(fid, 'fz_breakpoints', reduced_fz_vals);

% Write lookup tables
write_int16_3d_array(fid, 'fx', Fx_reduced);
write_int16_3d_array(fid, 'fy', Fy_reduced);

fclose(fid);
fprintf('\nC file "../tire_data/d2704/tire_data.c" has been generated.\n');

%% Estimate data sizes
% --- Data Type Sizes (in bytes)
bytes_per_int16 = 2;
bytes_per_float = 4;

% --- Dimensions
[num_sl, num_sa, num_fz] = size(Fx_reduced);

% --- Size Calculations
Fx_size_bytes = num_sl * num_sa * num_fz * bytes_per_int16;
Fy_size_bytes = Fx_size_bytes;  % Same size as Fx

sl_vals_size_bytes = numel(sl_vals) * bytes_per_float;
sa_vals_size_bytes = numel(reduced_sa_vals) * bytes_per_float;
fz_vals_size_bytes = numel(reduced_fz_vals) * bytes_per_float;

% --- Print estimated sizes
fprintf('\n=== Memory Footprint Estimate ===\n');
fprintf('Fx:        %d x %d x %d int16_t = %d bytes (%.2f KB)\n', ...
    num_sl, num_sa, num_fz, Fx_size_bytes, Fx_size_bytes / 1024);
fprintf('Fy:        %d x %d x %d int16_t = %d bytes (%.2f KB)\n', ...
    num_sl, num_sa, num_fz, Fy_size_bytes, Fy_size_bytes / 1024);
fprintf('sl_vals:   %d floats = %d bytes (%.2f KB)\n', ...
    numel(sl_vals), sl_vals_size_bytes, sl_vals_size_bytes / 1024);
fprintf('sa_vals:   %d floats = %d bytes (%.2f KB)\n', ...
    numel(reduced_sa_vals), sa_vals_size_bytes, sa_vals_size_bytes / 1024);
fprintf('fz_vals:   %d floats = %d bytes (%.2f KB)\n', ...
    numel(reduced_fz_vals), fz_vals_size_bytes, fz_vals_size_bytes / 1024);

total_size = Fx_size_bytes + Fy_size_bytes + sl_vals_size_bytes + sa_vals_size_bytes + fz_vals_size_bytes;
fprintf('---------------------------------------\n');
fprintf('Total Estimated Size: %d bytes (%.2f KB)\n\n', total_size, total_size / 1024);

%% Helper Functions
% --- Helper function to write float arrays ---
function write_float_array(fid, name, data)
    fprintf(fid, 'const float %s[%d] = {', name, numel(data));

    for i = 1:numel(data)
        if i < numel(data)
            fprintf(fid, '%.6ff, ', data(i));
        else
            fprintf(fid, '%.6ff', data(i));  % no comma after last element
        end
    end

    fprintf(fid, '};\n\n');
end

% --- Helper function to write 3D int16_t array ---
function write_int16_3d_array(fid, name, data)
    [dim1, dim2, dim3] = size(data);
    fprintf(fid, 'const int16_t %s[%d][%d][%d] = {\n', name, dim1, dim2, dim3);
    for i = 1:dim1
        fprintf(fid, '  {\n');
        for j = 1:dim2
            row = sprintf('%d, ', data(i, j, :));
            row = row(1:end-2);  % remove trailing comma
            fprintf(fid, '    { %s },\n', row);
        end
        fprintf(fid, '  },\n');
    end
    fprintf(fid, '};\n\n');
end
