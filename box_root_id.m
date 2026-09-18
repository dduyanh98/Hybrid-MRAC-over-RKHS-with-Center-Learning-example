function RootID = box_root_id(Boxes, BoxID)
if isfield(Boxes, 'root_id') && ~isempty(Boxes(BoxID).root_id)
    RootID = Boxes(BoxID).root_id;
else
    RootID = BoxID;
end
end
