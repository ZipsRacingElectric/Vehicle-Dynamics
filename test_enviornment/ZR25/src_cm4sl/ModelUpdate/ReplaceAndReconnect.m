function [] = ReplaceAndReconnect(subs, opts)
%Ssubstitute systems defined in subs and reconnect ports if possible.
% opts.OldPortNumsIn     -> Numbers of old Inports on new block
% opts.OldPortNumsOut  -> Numbers of old Outports on new block
% opts.FontSizeLabel      -> new font size to avoid overlapping labels

for (i=1:numel(subs))
    blck      = subs(i).blck;
    blck_src  = subs(i).blck_src;
    srcSys    = subs(i).srcSys;
    trgSys    = subs(i).trgSys;
    srcBlck   = strcat(srcSys, '/',blck_src);
    trgBlck   = strcat(trgSys,'/',blck);
    BlockType = get_param( trgBlck, 'BlockType');

    % save current ports and connections
    ports.old  = get_param(trgBlck,'PortHandles');
    h          = get_param(trgBlck,'LineHandles');
    lines_old  = SaveConnections(h);

    % clear all lines
    DeleteLines(h);

    % substitute block
    pos       = get_param(trgBlck, 'Position');
    new_block = replace_block(trgBlck, BlockType, srcBlck, 'noprompt');
    set_param(trgBlck,'Position', pos);

    % get all ports from new block
    ports.new = get_param(trgBlck,'PortHandles');

    % get all connections and labels from src block
    h         = get_param(srcBlck,'LineHandles');
    lines_src = SaveConnections(h);

    % Outputs: reconnect new block
    if (isfield(lines_src,'out')) % catch blocks with no outports
        for (iH=1:numel(lines_src.out)) % all ports of new block
            % copy labels from src blcok
            NewSrcPort = ports.new.Outport(iH);
            Label      = lines_src.out(iH).Label;
            set_param(NewSrcPort,'SignalNameFromLabel', Label);

            % reconnect if possible
            if (opts.OldPortNumsOut(iH) ~= -1)
                OldPortNum = opts.OldPortNumsOut(iH);
                OldDstPort = lines_old.out(OldPortNum).DstPortHandle;

                if (numel(opts.OldPortNumsOut) == 1)
                    h = add_line(trgSys, NewSrcPort,OldDstPort, 'autorouting','on');
                else
                    h = add_line(trgSys, NewSrcPort,OldDstPort, 'autorouting','off');
                end

                if (isfield(opts,'FontSizeLabel'))
                    set_param(h,'FontSize', opts.FontSizeLabel);
                end
            else
                h = TerminateUnconnected(trgSys, NewSrcPort, -1, -1);
                if (isfield(opts,'FontSizeLabel'))
                    set_param(h,'FontSize', opts.FontSizeLabel);
                end
            end

        end
    end

    % Inputs: reconnect new block
    if (isfield(lines_src,'in')) % catch blocks with no inports
        for (iH=1:numel(lines_src.in)) % all ports of new block
            % copy labels from src blcok
            NewDstPort = ports.new.Inport(iH); % no labels at Inports!

            % reconnect if possible
            if (opts.OldPortNumsIn(iH) ~= -1)
                OldPortNum = opts.OldPortNumsIn(iH);
                OldSrcPort = lines_old.in(OldPortNum).SrcPortHandle;
                if (numel(opts.OldPortNumsIn) == 1)
                    h = add_line(trgSys,OldSrcPort, NewDstPort, 'autorouting','on');
                else
                    h = add_line(trgSys,OldSrcPort, NewDstPort, 'autorouting','off');
                end
            else
                % addterms will be added later
            end


        end
        % add terminator to unconnetced new ports (if any)
        if (isfield(opts,'AddTerms') && opts.AddTerms == 1)
            addterms(trgSys);
        end
    end
end
end

