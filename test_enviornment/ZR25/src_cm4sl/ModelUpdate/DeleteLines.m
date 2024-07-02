% Delete all lines attached to Inports and Outports
function [] = DeleteLines(h)
    for (iH=1:numel(h.Outport))
        delete_line(h.Outport(iH));
    end
    for (iH=1:numel(h.Inport))
        delete_line(h.Inport(iH));
    end
end

