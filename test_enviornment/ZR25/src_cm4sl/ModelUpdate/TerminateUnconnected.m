function h = TerminateUnconnected(trgSys, Port, dx, size)
% terminate all unconnected ports

% default values if dx and size are are set to negative values
if (dx<0)
    dx = 20;
end
if (size<0)
    size = 5;
end

posPort = get_param(Port,'Position');
pos(1)  = posPort(1) + dx;
pos(2)  = posPort(2) - size/2;
pos(3)  = posPort(1) + dx + size;
pos(4)  = posPort(2) + size/2;

h   = add_block('simulink/Sinks/Terminator', [trgSys,'/Terminator'], 'Position', pos, 'MakeNameUnique', 'on');
tmp = get_param(h,'PortHandles');
h   = add_line(trgSys, Port, tmp.Inport);
end