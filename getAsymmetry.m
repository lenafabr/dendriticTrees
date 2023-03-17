function [asymmetry,nontermjunc] = getAsymmetry(NT,stM, stV)
% calculate asymmetry in mitochondrial volume density
% for network with radii and mitochondrial concentrations precalculated
% nontermjunc = list of nonterminal junctions
% junctions indexed by edge leading up to them

terminaledge = false(1,NT.nedge);
terminaljunc = false(1,NT.nedge);
nontermjunc = [];
avgC = zeros(1,NT.nedge);

for ec = 1:NT.nedge    
    junc = NT.edgenodes(ec,2);

    if (NT.degrees(junc)==1)
        terminaledge(ec) = true;        
        continue
    end
    
    edge1 = NT.nodeedges(junc,2); n1 = NT.nodenodes(junc,2);
    edge2 = NT.nodeedges(junc,3); n2 = NT.nodenodes(junc,3);
    if (NT.degrees(n1)==1 & NT.degrees(n2)==1);
        terminaljunc(ec) = true;
    else
        nontermjunc = [nontermjunc ec];
    end
    
    avgC(edge1) = stM(edge1)/(stV(edge1));
    avgC(edge2) = stM(edge2)/(stV(edge2));
    
    asymmetry(ec) = (avgC(edge1)-avgC(edge2))/(avgC(edge1)+avgC(edge2));
end


end