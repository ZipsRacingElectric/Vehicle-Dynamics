% test function to upload data from motec, mattGPT, other sources to the
% data pond (influxDB). At the moment, we are just using some test CSV 
% files. 
% 
% Not sure on:
% how to handle security, esp. with uploading to github. but, no worry for
% now, data retention policy is 30 days unless we pay. Long term, we should
% have cloud backup and slower, long term backup, as well as faster local
% server for at the track. 


% Written by: Ben Model, May 1 2025



%% Import test data
filepath = '..\DATA_POND\SampleDataForDataBucket.csv';

CleanWaterCols = string({"Time","lapDistance","rpm","a_y_","a_x_","dampFL","dampFR","dampRL","dampRR","SteeringDeg","throttle_","VehicleSpeed"});

testData = Data_Brita(filepath, CleanWaterCols);


%% Matlab's InfluxDB Example
% changed from a locally hosted jawn to their cloud jawn

% Authenticate with authentication token
setSecret("influxdbToken");
conn = influxdb("hostURL","http://influxdb:8086",...
"authToken", getSecret("influxdbToken"),"org","Mathworks");
    
% Authenticate with username and password
setSecret("usernameinfluxdb");
setSecret("passwordinfluxdb");


clusterURL = "https://us-east-1-1.aws.cloud2.influxdata.com";
VD_organizationID = "a1cdb447f71147e4";
UserID = "853b37c49d5133ba";


conn = influxdb("hostURL","http://influxdb:8086",...
    "username",getSecret("usernameinfluxdb"),"password",getSecret("passwordinfluxdb"),"org","Mathworks");

% Connect to default server localhost:8086
conn = influxdb("authToken",getSecret("influxdbToken"),...
    "org","mathworks");









