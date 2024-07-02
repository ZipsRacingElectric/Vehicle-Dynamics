function [SFun_matches] = FindSFun(sfunparam)
% Find system and location in src and trg model
% Search for sfunparam in FunctionNames and parameters
SetGlobal;
srcBlck = {};
trgBlck = {};

sfunname  = sfunparam;
sfunparam = ['''',sfunparam, '''']; % add quotes for strcmp

for (i=1:numel(SFun_list_src.FullName))
    % CM S-Fun and example S-Fun are identified differently
    if (strcmp(SFun_list_src.parameters{i}, sfunparam) || strcmp(SFun_list_src.FunctionName{i}, sfunname))
        srcBlck{end+1} = SFun_list_src.FullName{i};
    end
end
srcSys  = get_param(srcBlck, 'Parent');


for (i=1:numel(SFun_list_trg.FullName))
    if (strcmp(SFun_list_trg.parameters{i}, sfunparam) || strcmp(SFun_list_trg.FunctionName{i}, sfunname))
        trgBlck{end+1} = SFun_list_trg.FullName{i};
    end
end
trgSys    = get_param(trgBlck, 'Parent');
BlockType = get_param( trgBlck, 'BlockType');

SFun_matches.srcBlck   = srcBlck;
SFun_matches.srcSys    = srcSys; 
SFun_matches.trgBlck   = trgBlck; 
SFun_matches.trgSys    = trgSys; 
SFun_matches.BlockType = BlockType;
end

