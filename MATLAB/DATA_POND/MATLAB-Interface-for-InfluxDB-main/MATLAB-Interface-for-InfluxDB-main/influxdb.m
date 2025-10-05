classdef influxdb < handle
    % INFLUXDB http connection

    % conn = influxdb(NAME1,VALUE1,...NAMEN,VALUEN) creates a
    % connection object using Name-Value pairs to set connection details

    % Name-Value Pairs:
    % ----------------
    % hostURL - InfluxDB server API url. Default value - "http://localhost:8086/".
    % authToken - Authentication token for requests to InfluxDB server.
    % organization - Organization name (used as default in class methods,
    % unless specified).
    % username - User name. Used together with password for InfluxDB server authentication.
    % password - Password. used together with username for InfluxDB server authentication.
    %
    % influxdb class supports either authenticating using authentication
    % token or username and password.
    %
    % Methods:
    % --------
    % createBucket - Create a bucket of given name and retention policy(if any)
    % deleteBucket - Delete a bucket of given name
    % listBucket - List all the buckets in the database
    % bucketInfo - Get information about a bucket
    % healthCheck - Check if the influxDB server is up running
    % queryData - Read data from the influxDB server
    % writeData - Export data from MATLAB to the influxDB server
    %
    % Example:
    % --------
    % conn = influxdb("hostURL","http://influxdb:8086/","authToken",token,"org",org);
    %
    % conn = influxdb("username",username,"password",dbpassword,"org",org);

    % Copyright 2024 - 2025 The MathWorks, Inc.


    properties
        Org string = "";
    end

    properties(SetAccess = private)
        HostURL string;
    end

    properties(Access=private)
        HttpURI matlab.net.URI;
        CookieInfo matlab.net.http.CookieInfo = matlab.net.http.CookieInfo.empty();
        AuthToken string = "";
    end

    methods
        function this = influxdb(varargin)
            % Constructs the influxdb object with the specified hostname.
            % The default hosturl is localhost at port 8086

            p = inputParser();
            p.addParameter("hostURL","http://localhost:8086/",@(x)validateattributes(x,["char","string"],"scalartext"));
            p.addParameter("organization","",@(x)validateattributes(x,["char","string"],"scalartext"));

            p.addParameter("authToken","",@(x)validateattributes(x,["char","string"],"scalartext"));

            p.addParameter("username","",@(x)validateattributes(x,["char","string"],"scalartext"));
            p.addParameter("password","",@(x)validateattributes(x,["char","string"],"scalartext"));

            p.parse(varargin{:});
            this.HostURL = p.Results.hostURL;
            this.HttpURI = matlab.net.URI(this.HostURL);
            this.HttpURI.Path = "api/v2";

            authToken = p.Results.authToken;
            this.Org = p.Results.organization;

            usrname = p.Results.username;
            pwd = p.Results.password;

            if strlength(authToken)>0 ...
                    && (strlength(usrname)>0||strlength(pwd)>0)
                warning("database:influxDB:TwoAuthenticationMethodsWarning","The constructor takes either username and password, or authentication token as input parameters. Only authentication token will be used for authentication.")
            end

            if strlength(authToken)>0
                this.AuthToken = authToken;
            elseif strlength(usrname)>0
                % Use the Basic scheme for credential verification
                credentials = matlab.net.http.Credentials("Username",usrname,"Password",pwd,"Scheme","Basic");
                options = matlab.net.http.HTTPOptions("Credentials",credentials);
                request = matlab.net.http.RequestMessage('POST',[],matlab.net.http.MessageBody('boiler plater'));
                try
                    response = request.send(matlab.net.URI(strcat(this.HttpURI.EncodedURI,"/signin")),options);
                catch me
                    error(me.identifier,me.message)
                end
                setCookieFields = response.getFields("Set-Cookie");

                if isempty(setCookieFields)
                    error("database:influxDB:ServerError","Influxdb server:%s",response.Body.Data.message)
                else
                    this.CookieInfo = setCookieFields.convert;
                end

            else
                error("database:influxDB:AuthenticationError","The constructor takes either username and password, or authentication token as input parameters. Please input authentication information.")
            end
        end

        function createBucket(this, bucket, varargin)
            % createBucket method creates a bucket with given information
            % createBucket(conn, bucket)
            % creates a bucket with the given name bucket
            %
            % createBucket(conn, bucket, NAME1, VALUE1, ..., NAMEN, VALUEN)
            % creates a bucket with the given name bucket and other options
            %
            % Input Arguments:
            % ----------------
            % conn - influxdb connection object
            % bucket - bucket name

            % Name-value pairs:
            % -----------------
            % Description - bucket description
            % RetentionRules - a struct with retention policy information
            %
            % Examples:
            % ---------
            % createBucket(conn,"test-bucket")
            %
            % createBucket(conn,"test-bucket","Description","A bucket for
            % testing","retentionRules",struct("everySeconds",400000,"type","expire"))

            parser = inputParser;
            parser.addRequired("bucket",@(x)validateattributes(x,["char","string"],"scalartext"));
            parser.addParameter("org",this.Org,@(x)validateattributes(x,["char","string"],"scalartext"));
            parser.addParameter("description","",@(x)validateattributes(x,["char","string"],"scalartext"));
            parser.addParameter("retentionRules",struct,@(x)validateattributes(x,"struct","vector"));

            parser.parse(bucket,varargin{:})

            bucketName = parser.Results.bucket;
            orgName = parser.Results.org;
            bucketDescription = parser.Results.description;
            retentionRules = parser.Results.retentionRules;

            if isempty(fields(retentionRules))
                % empty retentionRules
                bucketInfo = struct('orgID', this.getOrgID(orgName), 'name', bucketName,...
                    'description', bucketDescription);
            else
                bucketInfo = containers.Map({'orgID','name','description','retentionRules'},...
                    {this.getOrgID(orgName), bucketName, bucketDescription,{retentionRules}});
            end

            data = matlab.net.http.MessageBody();
            data.Payload = jsonencode(bucketInfo);

            if isempty(this.CookieInfo) % AuthenToken is used
                request = matlab.net.http.RequestMessage( 'POST', ...
                    [matlab.net.http.HeaderField( 'Authorization', sprintf('Token %s', this.AuthToken)), ...
                    matlab.net.http.HeaderField('Content-type', 'application/json')],...
                    data);
            else % Username and password is used
                request = matlab.net.http.RequestMessage( 'POST', ...
                    [matlab.net.http.HeaderField('Content-type', 'application/json'),...
                    matlab.net.http.field.CookieField(this.CookieInfo.Cookie)],...
                    data);
            end
            try
                response = request.send(sprintf('%s/buckets', this.HttpURI.EncodedURI));
            catch me
                error(me.identifier,me.message)
            end

            if isstruct(response.Body.Data) && isfield(response.Body.Data,'message')
               error("database:influxDB:ServerError","Influxdb server:%s",response.Body.Data.message)
            end
        end

        function info = bucketInfo(this, bucket)
            % This method returns a struct with the information of the
            % targeting bucket
            %
            % info = bucketInfo(conn,bucket);
            %
            % Input Arguments:
            % ----------------
            % conn - influxdb connection object
            % bucket - bucket name

            % Outputs:
            % --------
            % info - a struct with bucket information of createt time,
            % updated time and retentionRules
            %
            % Examples:
            % ---------
            % info = bucketInfo(conn,"test-bucket");
            
            validateattributes(bucket,["char","string"],"scalartext");

            if isempty(this.CookieInfo) % AuthenToken is used
                request = matlab.net.http.RequestMessage( 'GET', ...
                    [matlab.net.http.HeaderField('Authorization', sprintf('Token %s', this.AuthToken)), ...
                    matlab.net.http.HeaderField('Content-type', 'application/json')]);
            else % Username and password is used
                request = matlab.net.http.RequestMessage( 'GET', ...
                    [matlab.net.http.HeaderField('Content-type', 'application/json'),...
                    matlab.net.http.field.CookieField(this.CookieInfo.Cookie)]);
            end
            bucketID = getBucketID(this,bucket);
            
            try
                response = request.send(sprintf('%s/buckets/%s', this.HttpURI.EncodedURI ,bucketID));
            catch me
                error(me.identifier,me.message)
            end

            if isstruct(response.Body.Data) && isfield(response.Body.Data,'message')
               error("database:influxDB:ServerError","Influxdb server:%s",response.Body.Data.message)
            end

            info = response.Body.Data;
            info = rmfield(info,["links","id","orgID","labels","type"]);
            if isfield(info,"createdAt")
                info.createdAt = timeString2datetime(this,info.createdAt);
            end
            if isfield(info,"updatedAt")
                info.updatedAt = timeString2datetime(this,info.updatedAt);
            end

        end

        function deleteBucket(this, bucketName)
            % deleteBucket method deletes an existing bucket of given name.
            %
            % deleteBucket(conn,bucket);
            %
            % Input Arguments:
            % ----------------
            % conn - influxdb connection object
            % bucket - bucket name
            %
            % Examples:
            % ---------
            % deleteBucket(conn,"test-bucket")

            % Get BucketID
            bucketID = getBucketID(this, bucketName);

            % Delete the bucket with the identifying bucketID
            if isempty(this.CookieInfo) % AuthenToken is used
                request = matlab.net.http.RequestMessage( 'DELETE', ...
                    [matlab.net.http.HeaderField( 'Authorization', sprintf('Token %s', this.AuthToken)), ...
                    matlab.net.http.HeaderField('Content-type', 'application/json')]);
            else % Username and password is used
                request = matlab.net.http.RequestMessage( 'DELETE', ...
                    [matlab.net.http.HeaderField('Content-type', 'application/json'),...
                    matlab.net.http.field.CookieField(this.CookieInfo.Cookie)]);
            end

            try
                response = request.send(sprintf('%s/buckets/%s',this.HttpURI.EncodedURI, bucketID));
            catch me
                error(me.identifier,me.message)
            end

            if isstruct(response.Body.Data)&&isfield(response.Body.Data,'message')
                error("database:influxDB:ServerError","Influxdb server:%s",response.Body.Data.message)
            end
        end

        function writeData(this, data, bucket, measurement, varargin)
            % writeData method writes data to an existing bucket.
            %
            % writeData(conn, data, bucket, measurement)
            % exports data from timetable data to influxDB server.
            %
            % writeData(this, data, bucket, measurement, NAME1, VALUE1,
            % ..., NAMEN, VALUEN)
            %
            % Input Arguments:
            % ----------------
            % conn - influxdb connection object
            % data - timetable data
            % bucket - Bucket name
            % measurement - Measurement name

            % Name-value pairs:
            % -----------------
            % Tags - A string array with tags names. The tags names must be
            % one of column names in the timedata variable data.
            % Org - Organization name.
            % BatchSize - The row size of data to be sent to the sever for
            % each request.
            %
            % Examples:
            % ---------
            % writeData(conn,TT,"write-test","measurement")
            % writeData(conn,TT,"write-test","measurement","Tags",["countryTags","cityTags"])


            p = inputParser;
            p.addRequired("data",@(x)validateattributes(x,"timetable","2d"));
            p.addRequired("bucket",@(x)validateattributes(x,["char","string"],"scalartext"));
            p.addRequired("measurement",@(x)validateattributes(x,["char","string"],"scalartext"));
            p.addParameter("Tags",string.empty,@(x)validateattributes(x,["char","string"],"vector"));
            p.addParameter("Org",this.Org,@(x)validateattributes(x,["char","string"],"scalartext"));
            p.addParameter("Debug",false,@(x)validateattributes(x,"logical","scalar"));
            p.addParameter("BatchSize",1e6,@(x)validateattributes(x,"numeric",{"scalar","finite","nonnegative"}));

            p.parse(data,bucket,measurement,varargin{:})

            data = p.Results.data;
            bucket = p.Results.bucket;
            measurement = p.Results.measurement;
            tags = p.Results.Tags;
            org = p.Results.Org;
            debug = p.Results.Debug;
            batchSize = p.Results.BatchSize;

            stringProvider = internal.influxDBStringProvider(data, tags, measurement,batchSize,debug);

            if isempty(this.CookieInfo) % AuthenToken is used
                request = matlab.net.http.RequestMessage( 'POST', ...
                    [matlab.net.http.HeaderField( 'Authorization', sprintf('Token %s', this.AuthToken)), ...
                    matlab.net.http.HeaderField('Content-Type', 'text/plain; charset=utf-8')],...
                    stringProvider);
            else % Username and password is used
                request = matlab.net.http.RequestMessage( 'POST', ...
                    [matlab.net.http.field.CookieField(this.CookieInfo.Cookie),...
                    matlab.net.http.HeaderField('Content-Type', 'text/plain; charset=utf-8')],...
                    stringProvider);
            end
            
            try
                response = request.send(sprintf('%s/write?org=%s&bucket=%s&precision=ns', ...
                    this.HttpURI.EncodedURI, org, bucket));
            catch me
                error(me.identifier,me.message)
            end

            if isstruct(response.Body.Data)&&isfield(response.Body.Data,'message')
                error("database:influxDB:ServerError","Influxdb server:%s",response.Body.Data.message)
            end
            
        end

        function dataTable = queryData(this, query,varargin)
            % queryData reads data from influxDB server. By default, the
            % columns _start, _stop, result and table are not included in
            % the output timetable/table.
            %
            % TT = queryData(conn,query)
            %
            % Input arugments:
            % ----------------
            % conn - influxdb connection object
            % query - a valid flux query
            %
            % Name-value arguments:
            % ---------------------
            % ReturnAllColumns - By default it's false. When set to true, the
            % columns _start, _stop, result and table will be included in
            % the output timetable/table.

            % Outputs:
            % --------
            % dataTable - If there's a time column returned, then the
            % output is a timetable; otherwise, the output is a table.

            % example query: 'from(bucket:"connection_test") |> range(start: 0) |> filter(fn: (r) => r._measurement == "weather")'
            %                'from(bucket:"connection_test") |> range(start: 0) |> filter(fn: (r) => r._measurement == "weather" and r._field== "Temp") |> mean()'
            %                'from(bucket:"connection_test") |> range(start: 0) |> filter(fn: (r) => r._measurement == "weather") |> pivot(rowKey:["_time"],columnKey:["_field"],valueColumn:"_value")'

            % fetch meta data
            p = inputParser;
            p.addRequired("query",@(x)validateattributes(x,["char","string"],"scalartext"))
            p.addParameter("ReturnAllColumns",false,@(x)validateattributes(x,"logical","scalar"));
            p.addParameter("Debug",false,@(x)validateattributes(x,"logical","scalar"));
            p.parse(query,varargin{:})
            
            Query = p.Results.query;
            returnAllColumns = p.Results.ReturnAllColumns;
            debug = p.Results.Debug;

            metaData.query = Query;
            metaData.type = 'flux';
            metaData.dialect.annotations = ["datatype" "group"];
            metaData.dialect.delimiter = ",";
            metaData.dialect.dateTimeFormat = "RFC3339";

            metaDataMessageBody = matlab.net.http.MessageBody();
            metaDataMessageBody.Payload = jsonencode(metaData);
            

            if isempty(this.CookieInfo) % AuthenToken is used
                metaDataRequest = matlab.net.http.RequestMessage( 'POST', ...
                [matlab.net.http.HeaderField('Authorization', sprintf('Token %s', this.AuthToken)), ...
                matlab.net.http.HeaderField( 'Content-type', 'application/json')], ...
                metaDataMessageBody);
            else % Username and password is used
                metaDataRequest = matlab.net.http.RequestMessage( 'POST', ...
                [matlab.net.http.field.CookieField(this.CookieInfo.Cookie),...
                matlab.net.http.HeaderField( 'Content-type', 'application/json')], ...
                metaDataMessageBody);
            end

            stringConsumer = internal.influxDBStringConsumer(returnAllColumns,debug);

            try
                response = metaDataRequest.send(sprintf('%s/query?org=%s', this.HttpURI.EncodedURI,this.Org), [], stringConsumer);
            catch me
                error(me.identifier,me.message)
            end
            if isstruct(response.Body.Data) && isfield(response.Body.Data,'message')
                error("database:influxDB:ServerError","Influxdb server:%s",response.Body.Data.message)
            end
            dataTable = lineProtocol2table(this,stringConsumer.dataList,stringConsumer.metadataList,stringConsumer.tableCount,returnAllColumns);

        end

        function bucketNames = listBuckets(this,organization)
            % listBuckets method lists all buckets that have been created in
            % the influxDB instance. The list returned might depend on the
            % permissions assigned to the Authtoken provided.
            %
            % bucketNames = listBuckets(conn);
            %
            % bucketNames = listBuckets(conn,org);
            %
            % Input arugments:
            % ----------------
            % conn - influxdb connection object
            % org - organization name
            %
            % Outputs:
            % --------
            % bucketNames - An array of strings of bucket names

            % Get Bucket list
            if nargin==1
                org = this.Org;
            else
                validateattributes(organization,["char","string"],"scalartext")
                org = organization;
            end
            myConsumer = matlab.net.http.io.JSONConsumer;
            if isempty(this.CookieInfo)
                request = matlab.net.http.RequestMessage( 'GET', ...
                    [matlab.net.http.HeaderField( 'Authorization', sprintf('Token %s', this.AuthToken))]);
            else
                request = matlab.net.http.RequestMessage('GET',...
                    matlab.net.http.field.CookieField(this.CookieInfo.Cookie));
            end

            try
                response = request.send(sprintf('%s/buckets?org=%s', this.HttpURI.EncodedURI, org), [], myConsumer);
            catch me
                error(me.identifier,me.message)
            end

            if isstruct(response.Body.Data) && isfield(response.Body.Data,'message')
                error("database:influxDB:ServerError","Influxdb server:%s",response.Body.Data.message)
            end

            buckets = response.Body.Data.buckets;
            nBuckets = size(buckets,1);
            bucketNames = strings(nBuckets,1);
            for i = 1:nBuckets
                bucketNames(i) = buckets{i}.name;
            end

        end




        function [status,message] = healthCheck(this)
            % healthCheck method checks if the influxdb instance at the specified port is
            % healthy.
            %
            % [status, msg] = healthCheck(conn);
            %
            % Input arugments:
            % ----------------
            % conn - influxdb connection object
            %
            % Outputs:
            % --------
            % status - Status of the server. Possible values are "pass" or
            % "fail".
            % message - Extra information of the server status.

            if isempty(this.CookieInfo) % AuthenToken is used
                request = matlab.net.http.RequestMessage( 'GET', ...
                    matlab.net.http.HeaderField( 'Authorization', sprintf('Token %s', this.AuthToken)));
            else % Username and password is used
                request = matlab.net.http.RequestMessage( 'GET', ...
                    matlab.net.http.field.CookieField(this.CookieInfo.Cookie));
            end
            uri = matlab.net.URI(this.HostURL);
            uri.EncodedPath = "/";

            try
                response = send(request, sprintf('%shealth',uri.EncodedURI));
            catch me
                error(me.identifier, me.message)
            end
            status = string(response.Body.Data.status);
            message = string(response.Body.Data.message);
        end
    end
    methods(Access = private)
        function orgID = getOrgID(this, org)
            if isempty(this.CookieInfo)
                request = matlab.net.http.RequestMessage( 'GET', ...
                    [matlab.net.http.HeaderField( 'Authorization', sprintf('Token %s', this.AuthToken)), ...
                    matlab.net.http.HeaderField('org', org)]);
            else
                request = matlab.net.http.RequestMessage( 'GET', ...
                    [matlab.net.http.HeaderField('org', org),...
                    matlab.net.http.field.CookieField(this.CookieInfo.Cookie)]);
            end

            try
                response = request.send(sprintf('%s/orgs', this.HttpURI.EncodedURI));
            catch me
                error(me.identifier,me.message)
            end

            if isstruct(response.Body.Data)&&isfield(response.Body.Data,'message')
                error("database:influxDB:ServerError",response.Body.Data.message)
            end

            orgID = response.Body.Data.orgs.id;
        end

        function bucketID = getBucketID(this, bucketName)
            % Get ID of specified bucketName.

            % Next steps:
            % - check if bucketName is valid and exists

            if isempty(this.CookieInfo) % AuthenToken is used
                request = matlab.net.http.RequestMessage( 'GET', ...
                    [matlab.net.http.HeaderField( 'Authorization', sprintf('Token %s', this.AuthToken)), ...
                    matlab.net.http.HeaderField( 'Accept', 'application/json' ),...
                    matlab.net.http.HeaderField('Content-type', 'application/json')]);
            else % Username and password is used
                request = matlab.net.http.RequestMessage( 'GET', ...
                    [matlab.net.http.HeaderField('Content-type', 'application/json'),...
                    matlab.net.http.HeaderField( 'Accept', 'application/json' ),...
                    matlab.net.http.field.CookieField(this.CookieInfo.Cookie)]);
            end

            myConsumer = matlab.net.http.io.JSONConsumer;

            try
                response = request.send(sprintf('%s/buckets?', this.HttpURI.EncodedURI),[],myConsumer);
            catch me
                error(me.identifier,me.message)
            end

            if ~isfield(response.Body.Data,'buckets') && isfield(response.Body.Data,'message')
                error("database:influxDB:ServerError",response.Body.Data.message)
            end
            
            if isa(response.Body.Data.buckets, 'struct')
                for iBkt = 1:numel(response.Body.Data.buckets)
                    bucketIDIdx(iBkt) = strcmpi(response.Body.Data.buckets(iBkt).name, bucketName); %#ok<AGROW>
                end
                selectedBucket = response.Body.Data.buckets(bucketIDIdx);
                if isempty(selectedBucket)
                    error("database:influxDB:BucketNotFound","The bucket %s doesn't exist.", bucketName)
                end
                bucketID = selectedBucket.id;
            else
                bucketIDIdx = cellfun(@(x)strcmpi(x.name, bucketName), response.Body.Data.buckets);
                selectedBucket = response.Body.Data.buckets(bucketIDIdx);
                if isempty(selectedBucket)
                    error("database:influxDB:BucketNotFound","The bucket %s doesn't exist.", bucketName)
                end
                bucketID = selectedBucket{1}.id;
            end
        end

        function T = lineProtocol2table(this,dataList,metadataList,tableCount,returnAllColumns)
            % helper function for converting the response of the a read query to a
            % MATLAB table format.
            dataNode = dataList.headNode; % headNode is a dummy node, no actual data stored
            metadataNode = metadataList.headNode; % headNode is a dummy node, no actual data stored
            % dataList.reset();

            % Process the very first real node
           if isempty(metadataNode.nextNode)
                error("database:influxDB:EmptyDataSet","Empty Data Set.");
            end
            dataNode = dataNode.nextNode;
            
            T = cell(tableCount,1);
            data = dataNode.getData();

            for i=1:tableCount
                metadataNode = metadataNode.nextNode;
                metadata = metadataNode.getData();
                typeArrayMod = metadata.dataType;
                for j=1:numel(typeArrayMod)
                    switch typeArrayMod(j)
                        case "dateTime:RFC3339"
                            typeArrayMod(j) = "string";
                        case "long"
                            typeArrayMod(j) = "int64";
                        case "unsignedLong"
                            typeArrayMod(j) = "uint64";
                        case "boolean"
                            typeArrayMod(j) = "logical";
                    end
                end

                if returnAllColumns
                    validColumns = ones(1,numel(metadata.columnName),"logical");
                else
                    validColumns = ~(strcmp(metadata.columnName,"_start")|strcmp(metadata.columnName,"_stop"));
                end
                columnNames = metadata.columnName(validColumns);
                typeArrayMod = typeArrayMod(validColumns);
                typeArray = metadata.dataType(validColumns);

                curTable = table('Size',[metadata.rowCount,size(typeArrayMod,2)],...
                'VariableTypes',typeArrayMod,'VariableNames',columnNames);
                curRow = 0;
                while metadata.rowCount > size(data,2)
                    metadata.rowCount = metadata.rowCount - size(data,2);
                    data = split(data',",",2);
                    if returnAllColumns
                        data = data(:,2:end);
                    else
                        data = data(:,4:end);
                    end
                    data = data(:,validColumns);
                    startRow = curRow + 1;
                    curRow = curRow + size(data,1);

                    for colN = 1:numel(columnNames)
                        switch typeArray(colN)
                            case {"string","dateTime:RFC3339"}
                                curTable.(columnNames(colN))(startRow:curRow) = data(:,colN);
                            case "double"
                                curTable.(columnNames(colN))(startRow:curRow) = double(data(:,colN));
                            case "long"
                                curTable.(columnNames(colN))(startRow:curRow) = int64(double(data(:,colN)));
                            case "unsignedLong"
                                curTable.(columnNames(colN))(startRow:curRow) = uint64(double(data(:,colN)));
                            case "boolean"
                                curTable.(columnNames(colN))(startRow:curRow) = strcmpi(data(:,colN),"true");
                        end
                    end
                    dataNode = dataNode.nextNode;
                    data = dataNode.getData();
                end
                if metadata.rowCount > 0
                    startRow = curRow + 1;
                    curRow = curRow + metadata.rowCount;
                    tempdata = data(1:metadata.rowCount);
                    tempdata = split(tempdata',",",2);
                    if returnAllColumns
                        tempdata = tempdata(:,(2:end));
                    else
                        tempdata = tempdata(:,(4:end));
                    end
                    tempdata = tempdata(:,validColumns);
                    for colN = 1:numel(columnNames)
                        switch typeArray(colN)
                            case {"string","dateTime:RFC3339"}
                                curTable.(columnNames(colN))(startRow:curRow) = tempdata(1:metadata.rowCount,colN);
                            case "double"
                                curTable.(columnNames(colN))(startRow:curRow) = double(tempdata(1:metadata.rowCount,colN));
                            case "long"
                                curTable.(columnNames(colN))(startRow:curRow) = int64(double(tempdata(1:metadata.rowCount,colN)));
                            case "unsignedLong"
                                curTable.(columnNames(colN))(startRow:curRow) = uint64(double(tempdata(1:metadata.rowCount,colN)));
                            case "boolean"
                                curTable.(columnNames(colN))(startRow:curRow) = strcmpi(tempdata(1:metadata.rowCount,colN),"true");
                        end
                    end
                    data = data(metadata.rowCount+1:end);
                    metadata.rowCount = 0;
                end
                if isempty(data)
                    dataNode = dataNode.nextNode;
                    if ~isempty(dataNode)
                        data = dataNode.getData();
                    end
                end

                for j=1:numel(columnNames)
                    if strcmp(typeArray(j),"dateTime:RFC3339")
                        curTable.(columnNames(j)) = timeString2datetime(this,curTable.(columnNames(j)));
                    end
                end
    
                % Check if the there's a column "_time". If true, return a
                % timetable, else return a table
                isTT = any(strcmp(columnNames,"_time"),2);
                if isTT
                    curTable = table2timetable(curTable,'RowTimes','_time');
                end

                T{i} = curTable;

            end
            

        end

        function datetimeArray = timeString2datetime(~,timeStringArr)
            timeStringArr = string(timeStringArr);
            tempData = replace(timeStringArr,["T","Z"],[" ",""]);
            if contains(tempData(1),".")
                datetimeArray = datetime(tempData,"inputFormat","yyyy-MM-dd HH:mm:ss.SSSSSSSSS");
            else
                datetimeArray = datetime(tempData,"inputFormat","yyyy-MM-dd HH:mm:ss");
            end
        end

    end
end
