function [junctedges,isterminal] = getJunctionEdges(NT)
% get the 3 edges (parent, then 2 daughters) for each junction
% isterminal = true for every junction where both daughters are dead ends

isterminal = true(0,1);
junctedges = [];

jc = 0;
for ec = 1:NT.nedge    
    junc = NT.edgenodes(ec,2);

    if (NT.degrees(junc)==1) % this is a dead-end not a junction
        continue
    end
    
    jc = jc+1;
    
    edge1 = NT.nodeedges(junc,2); n1 = NT.nodenodes(junc,2);
    edge2 = NT.nodeedges(junc,3); n2 = NT.nodenodes(junc,3);   
    isterminal(jc) = (NT.degrees(n1)==1 & NT.degrees(n2)==1);
    
    junctedges(jc,:) = [ec edge1 edge2];
end

end