classdef influxDBStringProvider < matlab.net.http.io.StringProvider
    
    %   Copyright 2024 - 2025 The MathWorks, Inc. 
    
    properties(Access=private)
        measurementName (1,1) string
        tagNames (1,:) string
        fieldNames (1,:) string
        batch (1,1) double
        startRow (1,1) double = 1
        endRow (1,1) double
        cellData
        dataHeight
        debug
    end

    methods
        function obj = influxDBStringProvider(dataIn, tagNames, measurementName,batch,debug)

            % Add escape character '\' in front of the special characters
            % Special characters in the measurement name,
            % tag names and values, field names and values:
            % measurement - comma, space
            % tag key - comma, equals sign, space
            % tag value - comma, equals sign, space
            % field key - comma, equals sign, space
            % field value - double quote, backslashs

            pattern = '([ ,=])';
            fieldValuePattern = '(["\\])';
            obj.measurementName = regexprep(measurementName,pattern, '\\$1');
            obj.debug = debug;

            obj.batch = batch;
            obj.endRow = batch;
            obj.dataHeight = height(dataIn);
            timeVar = dataIn.Properties.DimensionNames{1};
            
            variableNames = dataIn.Properties.VariableNames;
            if ~all(ismember(tagNames,variableNames))
                error("database:influxDB:TagsNotFound","tags %s should be part of the following: %s",strjoin(tagNames," "),strjoin(variableNames," "))
            end
            fieldNames = setdiff(variableNames,tagNames);
            obj.cellData = cell(numel(variableNames)+1,1);
            % Assign time array into the 1st cell of cellData
            obj.cellData{1} = convertTo(dataIn.(timeVar),"epochtime",'TicksPerSecond',1e9);
            
            for i=1:length(tagNames)
                % Assign tag arrays into the cells of cellData
                modifiedArray = regexprep(dataIn.(tagNames(i)), pattern, '\\$1');
                obj.cellData{i+1} = modifiedArray;
            end

            for i=1:length(fieldNames)
                array = dataIn.(fieldNames(i));
                if strcmpi(class(array(1)),"string")||...
                        (strcmpi(class(array(1)),"cell")&& strcmpi(class(array{1}),"char"))
                    array = regexprep(array, fieldValuePattern, '\\$1');
                end
                    obj.cellData{length(tagNames)+1+i} = array;
               
            end
            
            obj.tagNames = regexprep(tagNames, pattern, '\\$1');
            obj.fieldNames = regexprep(fieldNames, pattern, '\\$1');

            obj.ForceChunked = true;

        end

        function [data, stop] = getData(obj, ~)
            obj.Data = generateNextBufferOfData(obj);
            if obj.debug
                obj.show
            end
            if isempty(obj.Data) || strlength(obj.Data)==0
                stop = true;
                data = [];
            else
                data = unicode2native(char(obj.Data));
                stop = false;
            end
        end
    end

    methods(Access=private)

        function strData = generateNextBufferOfData(obj)
            if obj.startRow > obj.dataHeight
                strData = "";
                return
            end
            strEnd = min(obj.dataHeight,obj.endRow);
            strData = internal.convert2lineprotocol(obj.measurementName,obj.fieldNames,...
                                                obj.tagNames,obj.cellData,...
                                                obj.startRow-1,strEnd);
            obj.startRow = obj.endRow + 1;
            obj.endRow = obj.endRow + obj.batch;
        end
    end
end
