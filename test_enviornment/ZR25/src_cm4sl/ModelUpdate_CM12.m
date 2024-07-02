% Updates CM4SL Simulink models of Type CM11 to Type CM12
clear all;
close all;

addpath('ModelUpdate'); % additional scripts and functions for update
SetGlobal;

ValidProd       = {'CarMaker'}; % script only for CM4SL update
postfix_updated = '_updated';   % appendix for filename of updated model
postfix_backup  = '_backup';    % appendix for filename of backup model
publish_report  = 1;            % create html report of all changes (slow!)
ReferenceModel  =  'generic_src'; % Do not modify!

nConvert = 0;

fprintf('=====================================\n');
fprintf('Starting CM4SL/TM4SL models update script\n\n');
fprintf('Executing cmenv.m to ensure the necessary paths are included\n\n');
cmenv;

if(~isempty(bdroot))
    if strcmp(get_param(bdroot,'Dirty'),'on')
        fprintf('Aborting update script due to unsaved simulink models.\n');
        fprintf('Ensure all Simulink models are saved and closed before executing script\n');
        fprintf('=====================================\n');
        return
    end
end

fprintf('\n');
fprintf('Closing all open CM4SL/TM4SL models\n\n');
bdclose('all');

% select models
FileListAll = [dir('*.slx')' dir('*.mdl')']; % scan for all simulink models with extension .mdl and .slx
FileList    = FileListAll(:);         % select single models if desired (default: update all)

for iFile=1:numel(FileList)

    Mdl_FileName = FileList(iFile).name;
    [Mdl_FilePath, Mdl_BaseName, Mdl_FileExt] = fileparts(Mdl_FileName);

    % try/catch used in case some models can not be updated successfully
    try
        % load model which will be updated
        trg_mdl = Mdl_BaseName;
        load_system(trg_mdl);

        % check if model is valid for update
        if (IsValidForUpdate(Mdl_FileName, ValidProd)==1)

            if (numel(find_system(bdroot,'Name','CarMaker'))==1)
                XM = 'CM';
                XMaker = 'CarMaker';
            elseif (numel(find_system(bdroot,'Name','TruckMaker'))==1)
                XM = 'TM';
                XMaker = 'TruckMaker';
            end

            % Saving a backup of the simulation model if not yet existing
            save_backup;
           
            % Load source model (contains new S-Functions)
            src_mdl = ReferenceModel;
            load_system(src_mdl);

            % Store all IPG S-Functions
            SFun_list_src = FindAllSFun(src_mdl);
            SFun_list_trg = FindAllSFun(trg_mdl);

            % Annotations at first level of model
            update_Annotations

            %%% Update subsystems and blocks %%%
            % VehicleControl


            % IPG Vehicle/CarAndTrailer/


            % IPG Vehicle/Brake


            % IPG Vehicle/Powertrain
            update_PT_TCU_Example_in;
            update_PT_TCU_Example_out;

            update_PT_Control_Example_in;
            update_PT_Control_Example_out;

            update_PT_BCU_Example_in;
            update_PT_BCU_Example_out;


            % Saving the updated simulation model
            save_system(trg_mdl, strcat(trg_mdl, postfix_updated, Mdl_FileExt));

            % Save html report which contains all changes
            publish_changes;

            fprintf('Updated: %s\n', Mdl_FileName);
            nConvert = nConvert + 1;

            bdclose('all');
        end
    catch ME
        % Print error message and continue with next model in FileList
        fprintf('Failed to Update: %s.\nError Message was:\n%s', Mdl_FileName);
        warning(getReport(ME))
    end
end

fprintf('\n');
fprintf('Finishing CM4SL/TM4SL models update script\n');
fprintf('%d file updated\n', nConvert);
fprintf('=====================================\n');

bdclose('all');

