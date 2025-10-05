% Data Brita
% Filters data from data pond, returns table of desired data columns. Pass
% CleanWaterCols as string array.
% Abigail Tucker, 04/24/25 

% Example: FOR NOW UNTIL DATA POND IS ONLINE!!!!
% CSVs ARE TO BE ELIMINATED!!! DEAD!!! KILLED!!!! GONE!!!☢️🍄

%filepath = '..DATA_POND\SampleDataForDataBucket.csv'
%CleanWaterCols = string({'a_y_', 'a_x_'})
%CleanWater = DataBrita(filepath, CleanWaterCols)



function CleanWater = DataBrita(filepath, CleanWaterCols)

opts = detectImportOptions(filepath);
PondWater = readtable(filepath)

PondWaterCols = string(PondWater.Properties.VariableNames)

indexVector = not(any(CleanWaterCols(:)==PondWaterCols));


% clean water ISNT in pond water 🐸, kick out
CleanWater = removevars(PondWater, indexVector);

if length(CleanWaterCols) ~= length(CleanWater.Properties.VariableNames)
    for i=1:length(CleanWaterCols)
        if not(any(PondWaterCols(:)==CleanWaterCols(i)))
            % not in pond water 🐸
            fprintf('Column Name %s was not found\n',CleanWaterCols(i))
        end
    end
    
    error('🚫🐸 some of your variable names are missing, bro')


    end
end
