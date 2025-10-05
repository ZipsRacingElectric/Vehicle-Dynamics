% test function to upload data from motec, mattGPT, other sources to the
% data pond (influxDB). At the moment, we are just using some test CSV 
% files. 
% 
% Not sure on:
% how to handle security, esp. with uploading to github. but, no worry for
% now, data retention policy is 30 days unless we pay. Long term, we should
% have cloud backup and slower, long term backup, as well as faster local
% server for at the track. 


% Written by: Ben Model, Abigail Tucker, May 2025



%% Import test data
%filepath = '..\DATA_POND\SampleDataForDataBucket.csv';

%CleanWaterCols = string({"Time","lapDistance","rpm","a_y_","a_x_","dampFL","dampFR","dampRL","dampRR","SteeringDeg","throttle_","VehicleSpeed"});

%testData = Data_Brita(filepath, CleanWaterCols);


%% Matlab's InfluxDB Example
% changed from a locally hosted jawn to their cloud jawn
        % only needs to run once per account, can move to setup script or
        % something? or should edit to be in current path if env set to
        % VehicleDynamics/MATLAB like : "..\\DATA_POND" or something, IDK
addpath(genpath('C:\Users\ATuck\OneDrive - The University of Akron\Documents\GitHub\Vehicle-Dynamics\MATLAB\DATA_POND'))
savepath


% Authenticate with authentication token
    % needs to be run once? or multiple times? or movey to setuppy scripty?
getSecret("influxdbToken3")

% === CONFIGURATION ===
clusterURL = "https://us-east-1-1.aws.cloud2.influxdata.com";  % Your InfluxDB Cloud cluster
organization = "VehicleDynamics";  % subject to change as we shuffle around accounts and stuff still

% === CONNECT using saved token ===
conn = influxdb("hostURL", clusterURL, ...
                "authToken", getSecret("influxdbToken3"), ...
                "org", organization);



fluxQuery = ['from(bucket: "DATA_POND_TEST") ' ...
             '|> range(start: -30m) ' ...
             '|> filter(fn: (r) => r._measurement == "telemetry")'];




% Run query using method query()
data = conn.queryData(fluxQuery);

disp(data)


% TODO: create conditonal check if we actually are connected or not.
% print this if we are, or if not give some info or stuff. Or yell at them
        % angree frog
disp("✅ Connected to InfluxDB Cloud securely without exposing token.");

