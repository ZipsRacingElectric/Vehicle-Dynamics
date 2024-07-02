function [subs] = FindSubsFromSFun (sfunparam, ModSys, ModifiedPort, PortKind)
% Search for systems attached to S-Function 'sfunparam' at Port 'ModifiedPort'.
% 'ModSys' is the block name in the src system.
% 'PortKind' = Inport | Outport.
% 'ModifiedPort' = -1 stores the S-Fun itself in 'subs' for later substitution.
subs = [];


SFun_matches = FindSFun(sfunparam);

srcBlck   = SFun_matches.srcBlck;
srcSys    = SFun_matches.srcSys;
trgBlck   = SFun_matches.trgBlck;
trgSys    = SFun_matches.trgSys;
BlockType = SFun_matches.BlockType;

if (ModifiedPort == -1)
    % replace the S-Function(s) itself
    for (j=1:numel(trgBlck))
        subs(end+1).blck   = get_param(trgBlck{j},'Name');
        subs(end).trgSys   = trgSys{j};
        subs(end).blck_src = get_param(srcBlck{1},'Name');
        subs(end).srcSys   = srcSys{1};
    end
else
    for (j=1:numel(trgBlck))
        % get all ports on trg block
        ports.old = get_param(trgBlck{j},'PortHandles');

        % get all lines attached to trg block
        h     = get_param(trgBlck{j},'LineHandles');
        line = SaveConnections(h);

        if (strcmp(PortKind, 'Inport'))
            for (i=1:numel(line.in))
                if (line.in(i).DstPortHandle == ports.old.Inport(ModifiedPort))
                    subs(end+1).blck   = get_param(line.in(i).SrcBlockHandle,'Name');
                    subs(end).trgSys   = get_param(line.in(i).SrcBlockHandle,'Parent');
                    subs(end).blck_src = ModSys; % name from generic_truck in TM12
                    subs(end).srcSys   = srcSys{1};
                end
            end
        elseif (strcmp(PortKind, 'Outport'))
            for (i=1:numel(line.out))
                if (line.out(i).SrcPortHandle==ports.old.Outport(ModifiedPort))
                    subs(end+1).blck   = get_param(line.out(i).DstBlockHandle,'Name');
                    subs(end).trgSys   = get_param(line.out(i).DstBlockHandle,'Parent');
                    subs(end).blck_src = ModSys; % name from generic_truck in TM12
                    subs(end).srcSys   = srcSys{1};
                end
            end
        else
            % Bug/missuse
            warning("No valid PortKind (Inport | Outport). Got ''%s'' instead.", PortKind);
        end
    end
end
