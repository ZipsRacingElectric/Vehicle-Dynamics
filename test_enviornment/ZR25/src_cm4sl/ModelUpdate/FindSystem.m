function [System_matches] = FindSystem(blck)
% find system and location in src and trg model
SetGlobal;

srcBlck = find_system(src_mdl, 'name', blck);
srcSys  = get_param(srcBlck, 'Parent');

trgBlck = find_system(trg_mdl, 'name', blck);
trgSys  = get_param(trgBlck, 'Parent');

BlockType = get_param( trgBlck, 'BlockType');

System_matches.srcBlck   = srcBlck;
System_matches.srcSys    = srcSys; 
System_matches.trgBlck   = trgBlck; 
System_matches.trgSys    = trgSys; 
System_matches.BlockType = BlockType;
end

