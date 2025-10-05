# Getting Started with MATLAB® Interface for InfluxDB

## Description

This interface allows users to import data from influxDB v2 server to MATLAB® and allows users to export data from MATLAB® to influxDB v2 server.
## System Requirements

- MATLAB R2025a or later
- Database Toolbox R2025a or later

## Features

Users can retrieve influxDB data directly from MATLAB. This software is written based on InfluxDB v2 API documentation, which can be found here:
https://docs.influxdata.com/influxdb/v2/reference/api/

This software supports authenticating using authentication token and username, password. InfluxDB server provides such authentication information. 

## Create a influxdb connection object
```MATLAB
% Authenticate with authentication token
setSecret("influxdbToken");
conn = influxdb("hostURL","http://influxdb:8086",...
"authToken", getSecret("influxdbToken"),"org","Mathworks");
    
% Authenticate with username and password
setSecret("usernameinfluxdb");
setSecret("passwordinfluxdb");
conn = influxdb("hostURL","http://influxdb:8086",...
    "username",getSecret("usernameinfluxdb"),"password",getSecret("passwordinfluxdb"),"org","Mathworks");

% Connect to default server localhost:8086
conn = influxdb("authToken",getSecret("influxdbToken"),...
    "org","mathworks");
```

## Write data to an existing bucket
Here's the mapping between MATLAB data types and InfluxDB data types when exporting data to InfluxDB server:
|MATLAB data types|InfluxDB data types|
|-----------------|-------------------|
|chars, string|string|
|datetime|unix timestamp (precision ns)|
|int8,int16,int32,int64|signed 64-bit integer|
|logical|boolean|
|single, double|float|
|uint8,uint16,uint32,uint64|unsigned 64-bit integer|

```MATLAB
MeasurementTime = datetime({'2015-12-18 08:03:05';'2015-12-18 10:03:17';'2015-12-18 12:03:13'});
Temp = [37.3;39.1;42.3];
Pressure = [30.1;30.03;29.9];
WindSpeed = [13.4;6.5;7.3];
ID = uint64([1;2;3]);
Description=repmat("Weather Data",3,1);
Locations = ["New York";"Boston";"New York"];
TT = timetable(MeasurementTime, ID, Temp,Pressure,WindSpeed,Locations,Description);
 
bucketName = "write-test";
createBucket(conn,bucketName);
measurementName = "basic-test";
tag = "Locations";
writeData(conn,TT,bucketName,measurementName,"Tags",tag);
```
## Query Data from influxDB server
Here's the mapping between MATLAB data types and InfluxDB data types when querying data from InfluxDB server:
|InfluxDB Data types|MATLAB Data Types|
|---|---|
|unsigned 64-bit integer|uint64|
|unix timestamp (precision ns)|datetime|
|string|string|
|signed 64-bit integer|int64|
|float|double|
|boolean|logical|

```MATLAB
>> query = ['from(bucket:"write-test") |> range(start: 0) |> filter(fn: (r) => r._measurement == "basic-test")', ...
'|> pivot(rowKey: ["_time"], columnKey: ["_field"], valueColumn: "_value")'];
 
 
>> T = queryData(conn, query)
T =

  1×1 cell array

    {1000×8 timetable}
>> T{1}
ans = 
   1000×8 timetable

           _time            _measurement    cityTags    countryTags    charTypedCol    doubleTypedCol    intTypedCol    stringTypedCol    uintTypedCol
    ____________________    ____________    ________    ___________    ____________    ______________    ___________    ______________    ____________

    01-Jan-2000 00:00:03    "basic-test"    "Boston"      "US"            "test"               3               3            "test"               3    
    01-Jan-2000 00:00:04    "basic-test"    "Boston"      "US"            "test"               4               4            "test"               4    
    01-Jan-2000 00:00:05    "basic-test"    "Boston"      "US"            "test"               5               5            "test"               5    
    01-Jan-2000 00:00:08    "basic-test"    "Boston"      "US"            "test"               8               8            "test"               8    
```

## Bucket Management
```MATLAB
% Create a bucket named "write-test"
createBucket(conn,"write-test");

% Check the information of a bucket
info = bucketInfo(conn,"write-test");

% Create a bucket with description and retention rule
createBucket(conn,"write-test","Description","A bucket for testing writeData function","retentionRules",struct("everyseconds",4e5,"type","expire"));

% Delete the bucket "write-test"
deleteBucket(conn,"write-test");

% List all the buckets within an organization
buckets = listBuckets(conn,"mathworks");

```

## Check if the server is alive
```MATLAB
>> [status, message] = healthCheck(conn)

status =
 
    "pass"
 
 
message =
 
    "ready for queries and writes"
```
