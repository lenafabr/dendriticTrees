function [stL,stV,stD] = setSubtreeInfo_fromRadii(NT,trunkedge,stL,stV,stD,radii)
% walk down a tree, starting with the trunkedge
% compute information on subtrees based on tree topologies  and
% precalculated or measured radii.
% does assume a scaling exponent alpha
% stL = total length of subtree below each junction
% stV = volume for the subtree below each junction
% stD = depth for the subtree below each junction


% if input stL, stEta, muvals arrays not provided (or empty), initialize to zeros
if (~exist('stL','var'))
    % initialize all arrays to 0
    stL = zeros(1,NT.nedge);    
end
if (~exist('stV','var'))
    % initialize all arrays to 0
    stV = zeros(1,NT.nedge);    
end
if (~exist('stD','var'))
    % initialize all arrays to 0
    stD = zeros(1,NT.nedge);    
end

if (isempty(stL))
    % initialize all arrays to 0
    stL = zeros(1,NT.nedge);    
end
if (isempty(stV))
    % initialize all arrays to 0
    stV = zeros(1,NT.nedge);    
end
if (isempty(stD))
    % initialize all arrays to 0
    stD = zeros(1,NT.nedge);    
end


%trunk length
ell0 = NT.edgelens(trunkedge);

%% junction node below this edge
junc = NT.edgenodes(trunkedge,2);

if (NT.degrees(junc)==1)
    % reached a terminal node
    % this edge is a single-trunk subtree
    stV(trunkedge) = ell0*radii(trunkedge)^2;   
    stD(trunkedge) = ell0;    
    stL(trunkedge) = ell0;

elseif (NT.degrees(junc)==2)
    error('not set up to deal with deg 2 nodes yet')
    
elseif (NT.degrees(junc)==3)
   
    edge1 = NT.nodeedges(junc,2);
    edge2 = NT.nodeedges(junc,3);
    
    % compute subtree properties for the daughter trunks
    % including radii splitting (mu) for all downstream junctions
    [stL,stV,stD] = setSubtreeInfo_fromRadii(NT,edge1,stL,stV,stD,radii);
    [stL,stV,stD] = setSubtreeInfo_fromRadii(NT,edge2,stL,stV,stD,radii);
    
    % tree length
    stL(trunkedge) = stL(edge1)+stL(edge2)+ell0;    
    % tree volume
    stV(trunkedge) = stV(edge1)+stV(edge2)+ell0*radii(trunkedge).^2;
    
    % subtree depth
    stD(trunkedge) = ell0 + (stL(edge1) + stL(edge2))/(stL(edge1)/stD(edge1) + stL(edge2)/stD(edge2));
else
    error('not set up to deal with node degrees > 3')
end



end