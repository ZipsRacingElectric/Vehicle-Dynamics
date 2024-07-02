% Update Trailer
subs = [];
opts = [];


sfunparam    = 'cm_ptcontrol_out';
ModSys       = 'ignored';
ModifiedPort = -1; % replace the S-Function itself
PortKind     = 'ignored';

subs_add = FindSubsFromSFun (sfunparam, ModSys, ModifiedPort, PortKind);
subs     = [subs, subs_add];

opts.OldPortNumsIn  = [1:8, -ones(1,4), 9:14, -1, -1, 15];
opts.OldPortNumsOut = [];
% opts.FontSizeLabel  = 2;
opts.AddTerms       = 1;

ReplaceAndReconnect(subs, opts);

clear opts subs;