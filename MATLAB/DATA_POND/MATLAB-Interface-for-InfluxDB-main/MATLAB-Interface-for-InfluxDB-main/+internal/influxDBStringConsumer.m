classdef influxDBStringConsumer < matlab.net.http.io.StringConsumer
    % Copyright 2024 - 2025 The MathWorks, Inc.
    properties
        dataList
        buffer
        metadataList
        tableCount
    end

    properties(Access = private)
        ReturnAllColumns
        Debug
    end

    methods (Access = protected)
        function extractString(obj,strData)
            keyword = ',_result';
            if ~isempty(strData)
                if obj.Debug
                    disp(strData)
                end
                if isempty(obj.buffer)
                    obj.buffer = '';
                end
                % Find the last occurence of _result. Combine the string before
                % the last _result and the old buffer. Add it to dataList.
                % Update the buffer to the string after _result (including
                % _result)
                
                indx = findTheLastMatch(obj,strData,keyword);

                oldLen = strlength(obj.buffer);
                if indx>0
                    len = oldLen + indx-1;
                else
                    len = oldLen + strlength(strData);
                end
                % if len>0 % len is the new combined buffer length
                % Copy the characters in the old buffer and strData
                % (all the data before the last _result) into newChar
                newChar = char(nan(1,len));
                newChar(1:oldLen) = obj.buffer;
                if indx>0
                    newChar((oldLen+1):end) = strData(1:indx-1);
                    obj.buffer = strData(indx:end);
                    strArray = strtrim(strsplit(string(newChar),'\n'));
                    strArray = strArray(strlength(strArray)>0);
                    % strArray has to contain complete metadata, or no meta
                    % data (all or nothing)
                    obj.addStrData(strArray);
                else
                    newChar((oldLen+1):end) = strData;
                    obj.buffer = newChar;
                end
            else
                if strlength(obj.buffer)>0
                    % Add old buffer contents to dataList
                    strArray = strtrim(strsplit(string(obj.buffer),'\n'));
                    strArray = strArray(strlength(strArray)>0);
                    if ~isempty(strArray)
                        obj.addStrData(strArray);
                    end
                end
            end
            obj.CurrentLength = strlength(obj.buffer);
        end
    end

    methods
        function obj = influxDBStringConsumer(returnAllColumns,debug)
            obj.dataList = database.internal.utilities.LinkedList;
            obj.metadataList = database.internal.utilities.LinkedList;
            obj.AppendFcn = @(obj,data)obj.extractString(data);
            obj.ReturnAllColumns = returnAllColumns;
            obj.tableCount = 0;
            obj.Debug = debug;
        end
    end

    methods(Access = private)
        function indx = findTheLastMatch(~,inputStr,pattern)

            start = numel(inputStr)-numel(pattern)+1;
            for i= start:-1:1
                if inputStr(i)==pattern(1) && all(inputStr(i:i+numel(pattern)-1)==pattern)
                    indx = i;
                    return;
                end
            end
            % When no match found, or inputStr is shorter than pattern,
            % return 0;
            indx = 0;
        end

        function addStrData(obj,strArray)
            % strArray has to contain complete metadata, or no meta
            % data (all or nothing)
            metaDataIndicator = startsWith(strArray,",result");
            metaDataIndx = find(metaDataIndicator);
            metaDataGroupNum = numel(metaDataIndx);
            for i = 1:metaDataGroupNum
                metaDataArray = strArray((metaDataIndx(i)-2):metaDataIndx(i));
                metaDataArray = split(metaDataArray',",",2);
                metaData = struct();
                if obj.ReturnAllColumns
                    metaData.dataType = metaDataArray(1,2:end);
                    metaData.columnName = metaDataArray(3,2:end);
                else
                    metaData.dataType = metaDataArray(1,4:end);
                    metaData.columnName = metaDataArray(3,4:end);
                end
                metaData.rowCount = 0;
                if obj.metadataList.isempty()
                    obj.metadataList.add(metaData);
                elseif i == 1
                    obj.metadataList.lastNode.data.rowCount = obj.metadataList.lastNode.data.rowCount + metaDataIndx(i) - 3;
                    obj.metadataList.add(metaData);
                else
                    obj.metadataList.lastNode.data.rowCount = obj.metadataList.lastNode.data.rowCount + metaDataIndx(i) - metaDataIndx(i-1) - 3;
                    obj.metadataList.add(metaData);
                end
                obj.tableCount = obj.tableCount + 1;
                % change the metaDataIndicator to indicate all the
                % metadata rows, not just the column names rows
                metaDataIndicator(metaDataIndx(i)-2)=true;
                metaDataIndicator(metaDataIndx(i)-1)=true;

            end
            if ~isempty(metaDataIndx)
                obj.metadataList.lastNode.data.rowCount = obj.metadataList.lastNode.data.rowCount + size(strArray,2)-metaDataIndx(end);
            else
                obj.metadataList.lastNode.data.rowCount = obj.metadataList.lastNode.data.rowCount + size(strArray,2);
            end
            
            strArray = strArray(~metaDataIndicator);
            if strlength(strArray)~=0
                %strArray could be empty when there's only a row of data
                %and it is stored in buffer. See g3432359
                obj.dataList.add(strArray);
            end
        end
    end
end
