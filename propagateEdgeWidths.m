function propagateEdgeWidths(NT,edge,w)
% go down the tree, setting any undefined edge widths to the width of the 
% parent edge (if set) or to the default w if not set

if (isnan(NT.rootnode) | isempty(NT.rootnode))
    error('Network has not been set up as a directed tree')
end

if (isnan(w))
    % default width
    w = 1;
end
if (isnan(edge))
    % default edge, start at trunk (longest edge from current node)
    rootedges = NT.nodeedges(NT.rootnode,1:NT.degrees(NT.rootnode));
    [~,ind] = max(NT.edgelens(rootedges));
    edge = rootedges(ind);
end

if (isempty(NT.edgewidth{edge}))
    % set current edge to the desired width
    NT.edgewidth{edge} = makeEdgeWidth(NT,edge,w);
end

% downstream node
downnode = NT.edgenodes(edge,2);
deg = NT.degrees(downnode);
if (deg == 1)
    % terminal node
    return
end

% child edges
childwidth = mean(NT.edgewidth{edge}(:,1));
if (deg>1)
    for dc = 2:deg
        childedge = NT.nodeedges(downnode,dc);
        propagateEdgeWidths(NT,childedge,childwidth);
    end
end

end