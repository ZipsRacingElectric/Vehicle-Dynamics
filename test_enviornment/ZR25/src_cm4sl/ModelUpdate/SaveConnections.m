function line = SaveConnections(h)
% save all connected ports
if (isfield(h,'Outport'))
    for (iH=1:numel(h.Outport))
        line.out(iH).n              = iH;
        line.out(iH).LineHandle     = h.Outport(iH);
        line.out(iH).SrcBlockHandle = get_param(h.Outport(iH), 'SrcBlockHandle');
        line.out(iH).SrcPortHandle  = get_param(h.Outport(iH), 'SrcPortHandle');
        line.out(iH).DstBlockHandle = get_param(h.Outport(iH), 'DstBlockHandle');
        line.out(iH).DstPortHandle  = get_param(h.Outport(iH), 'DstPortHandle');
        line.out(iH).Label          = get_param(line.out(iH).SrcPortHandle, 'SignalNameFromLabel');
    end
end
if (isfield(h,'Inport'))
    for (iH=1:numel(h.Inport))
        line.in(iH).n               = iH;
        line.in(iH).LineHandle      = h.Inport(iH);
        line.in(iH).SrcBlockHandle  = get_param(h.Inport(iH), 'SrcBlockHandle');
        line.in(iH).SrcPortHandle   = get_param(h.Inport(iH), 'SrcPortHandle');
        line.in(iH).DstBlockHandle  = get_param(h.Inport(iH), 'DstBlockHandle');
        line.in(iH).DstPortHandle   = get_param(h.Inport(iH), 'DstPortHandle');
        line.in(iH).Label           = get_param(line.in(iH).DstPortHandle, 'SignalNameFromLabel');
    end
end
end

