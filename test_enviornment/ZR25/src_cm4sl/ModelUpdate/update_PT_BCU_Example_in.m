subs = []; % used to store componets that will be replaced/reconnected
opts = [];


sfunparam    = 'cm_ptbatterycu_in';
ModSys       = 'CreateBus BCU_Cfg BattLV';
ModifiedPort = 1; % [FL SpringFrc, ...]
PortKind     = 'Outport';

subs_add = FindSubsFromSFun (sfunparam, ModSys, ModifiedPort, PortKind);
subs     = [subs, subs_add];

sfunparam    = 'cm_ptbatterycu_in';
ModSys       = 'CreateBus BCU_Cfg BattHV';
ModifiedPort = 2; % [FL SpringFrc, ...]
PortKind     = 'Outport';

subs_add = FindSubsFromSFun (sfunparam, ModSys, ModifiedPort, PortKind);
subs     = [subs, subs_add];

sfunparam    = 'cm_ptbatterycu_in';
ModSys       = 'CreateBus BCU BattLV';
ModifiedPort = 4; % [FL SpringFrc, ...]
PortKind     = 'Outport';

subs_add = FindSubsFromSFun (sfunparam, ModSys, ModifiedPort, PortKind);
subs     = [subs, subs_add];

sfunparam    = 'cm_ptbatterycu_in';
ModSys       = 'CreateBus BCU BattHV';
ModifiedPort = 5; % [FL SpringFrc, ...]
PortKind     = 'Outport';

subs_add = FindSubsFromSFun (sfunparam, ModSys, ModifiedPort, PortKind);
subs     = [subs, subs_add];



opts.OldPortNumsIn  = [1]; % Numbers of old ports on new block
opts.OldPortNumsOut = [1];
% opts.FontSizeLabel  = 2;
opts.AddTerms       = 1;

ReplaceAndReconnect(subs, opts);

clear opts subs;