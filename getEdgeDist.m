function [dist] = getEdgeDist(NT,trunkedge,dist)

    if (~exist('dist','var'))    
        dist = zeros(1,NT.nedge);        
    end
    if (isempty(dist))
        % initialize all arrays to 0
        dist = zeros(1,NT.nedge);        
    end
    %% calculate graph distance of each edge away from the rootnode
    
    % distance of this edge from its own start
    edgedist = NT.edgelens(trunkedge)/2;
    
    % node from which trunk comes
    rootnode = NT.edgenodes(trunkedge,1);
    if (NT.degrees(rootnode)==1)
        % parent node of whole tree
        dist(trunkedge) = edgedist;
    elseif (NT.degrees(rootnode)==3)        
        parentedge = NT.nodeedges(rootnode,1); % parent edge to this trunk    
        dist(trunkedge) = dist(parentedge) + NT.edgelens(parentedge)/2 + edgedist;
    else
        error('bad degree')
    end
    
    % downstream edges
    junc = NT.edgenodes(trunkedge,2);

    if (NT.degrees(junc) ==3) % junction node        
        edge1 = NT.nodeedges(junc,2);
        edge2 = NT.nodeedges(junc,3);
        
        [dist] = getEdgeDist(NT,edge1,dist);
        [dist] = getEdgeDist(NT,edge2,dist)    ;    
    end
end