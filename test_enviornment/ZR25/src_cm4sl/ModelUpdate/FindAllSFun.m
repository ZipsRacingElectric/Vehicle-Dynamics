function [SFun_list] = FindAllSFun(mdl)
% Find all IPG SFuns in src and trg mdl
% Store information in list list for later substitutions
SFun_list.FullName = [];

IPG_SFuns = {'TM_Sfun',  'CM_Sfun',  'MM_Sfun',...
             'cm_ptcontrol_in','cm_ptcontrol_out',...
             'cm_pttransmcu_in','cm_pttransmcu_out',...
             'cm_ptbatterycu_in','cm_ptbatterycu_out'
    };

% Find S-Functions
for (i=1:numel(IPG_SFuns))
    tmp  = find_system(mdl, 'IncludeCommented', 'on', 'FunctionName', IPG_SFuns{i});
    SFun_list.FullName = [SFun_list.FullName; tmp];
end

% Get additional information
for (i=1:numel(SFun_list.FullName))
    h = SFun_list.FullName{i};

    SFun_list.parameters{i,1}   = get_param(h, 'parameters');
    SFun_list.Handle{i,1}       = get_param(h, 'Handle');
    SFun_list.FunctionName{i,1} = get_param(h, 'FunctionName');
    SFun_list.Name{i,1}         = get_param(h, 'Name');
end

end

