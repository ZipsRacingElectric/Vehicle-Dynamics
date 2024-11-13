% Battery parameters for ZR25 design.
% 9/26/24

% TODO: load data from CSV so it can later be replaced with test data

%% Individual Cell Data (COSMX 13Ah cell) 

% rows: SOC increasing columns: temps increasing
cell_OCV = [
    3.49, 3.5, 3.51;
    3.55, 3.57, 3.56;
    3.62, 3.63, 3.64;
    3.71, 3.71, 3.72;
    3.91, 3.93, 3.94;
    4.07, 4.08, 4.08;
    4.19, 4.19, 4.19];

% rows: SOC increasing columns: temps increasing
cell_DC_IR = [
    .0117, .0085, .009;
    .011, .0085, .009;
    .0114, .0087, .0092;
    .0107, .0082, .0088;
    .0107, .0083, .0091;
    .0113, .0085, .0089;
    .0116, .0085, .0089];

cell_nominal_capacity = 13;

%% battery pack configuration (it assumes all segments are identical)
segment_parallel_cells = 1;
segment_series_cells = 24;
pack_parallel_segments = 1;
pack_series_segments = 6;

%% Battery Pack Data
% Compile data into parameters that can be used by a single Simscape
% table-based battery block:

pack_SOC = Simulink.Parameter([0, .1, .25, .5, .75, .9, 1]);                              
pack_temps = Simulink.Parameter([5, 20, 40]); 

pack_OCV = Simulink.Parameter(cell_OCV * segment_series_cells * pack_series_segments);
pack_DC_IR = Simulink.Parameter((((cell_DC_IR / segment_parallel_cells) * segment_series_cells) / pack_parallel_segments) * pack_series_segments);
pack_nominal_capacity = Simulink.Parameter(cell_nominal_capacity * segment_parallel_cells * pack_parallel_segments);

%% Battery Pack Design Data
segment_bus_bar_resistance = 0.01 * pack_series_segments / pack_parallel_segments;
cell_interconnect_resistance = (0.01 * (segment_series_cells + 1) / segment_parallel_cells) * pack_series_segments / pack_parallel_segments;
segment_connector_resistance = 0.01 * pack_series_segments / pack_parallel_segments;
maintence_plug_cable_resistance = 0.01 * pack_series_segments / pack_parallel_segments;
isolation_relay_contact_resistance = 0.01;
hvd_contact_resistance = 0.01;
ts_cable_resistance =  (0.0171E-6) * 2 / (pi * (8.37/1000/1000)); % R = p*l/a

ts_equivalent_resistance = Simulink.Parameter(segment_bus_bar_resistance + cell_interconnect_resistance + segment_connector_resistance + maintence_plug_cable_resistance + isolation_relay_contact_resistance + hvd_contact_resistance + ts_cable_resistance);
