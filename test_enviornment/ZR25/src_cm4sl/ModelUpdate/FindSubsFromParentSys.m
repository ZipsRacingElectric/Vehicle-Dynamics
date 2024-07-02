function [subs] = FindSubsFromParentSys (blck, sys, ModSys, ModifiedPort, PortKind)
% Search for a system which lays in a given parent system
% Is used if systems are to be mofified that are not directly connected to
% unique S-Functions and may exist more than once at differnt places
subs = [];

% Find all blocks with matching name
System_matches = FindSystem(blck);

srcBlck   = System_matches.srcBlck;
srcSys    = System_matches.srcSys;
trgBlck   = System_matches.trgBlck;
trgSys    = System_matches.trgSys;
BlockType = System_matches.BlockType;

% find system with given parent in src system
for (j=1:numel(srcBlck))
    if (strcmp(get_param(srcSys{j},'Name'),sys))
        srcSys_static = srcSys{j};
    end
end

for (j=1:numel(trgBlck))
    % Only update instance in current sys
    if (strcmp(get_param(trgSys{j},'Name'),sys))
        % get all ports on trg block
        ports.old = get_param(trgBlck{j},'PortHandles');

        % get all lines attached to trg block
        h    = get_param(trgBlck{j},'LineHandles');
        line = SaveConnections(h);

        if (strcmp(PortKind, 'Inport'))
            for (i=1:numel(line.in))
                if (line.in(i).DstPortHandle == ports.old.Inport(ModifiedPort))
                    subs(end+1).blck   = get_param(line.in(i).SrcBlockHandle,'Name');
                    subs(end).trgSys   = get_param(line.in(i).SrcBlockHandle,'Parent');
                    subs(end).blck_src = ModSys; % name from generic_truck in TM12
                    subs(end).srcSys   = srcSys_static;
                end
            end
        elseif (strcmp(PortKind, 'Outport'))
            for (i=1:numel(line.out))
                if (line.out(i).SrcPortHandle==ports.old.Outport(ModifiedPort))
                    subs(end+1).blck   = get_param(line.out(i).DstBlockHandle,'Name');
                    subs(end).trgSys   = get_param(line.out(i).DstBlockHandle,'Parent');
                    subs(end).blck_src = ModSys; % name from generic_truck in TM12
                    subs(end).srcSys   = srcSys_static;
                end
            end
        else
            % Bug/missuse
            warning("No valid PortKind (Inport | Outport). Got ''%s'' instead.", PortKind);
        end
    end
end

end


