function [radii] = setRadiiFromSubtreeInfo(NT,trunkedge,a,rtrunk,muvals,radii)
% walk down a tree, starting with the trunkedge
% given a radius value for the trunk edge (rtrunk),
% set the radii of the rest of the branches.
% assumes the subtree information (splitting ratios, etc) has
% already been calculated with setSubtreeInfo


% if just starting the recursion,
% initialize arrays to 0 
if (~exist('radii','var'))    
    radii = zeros(1,NT.nedge);        
end
if (isempty(radii))
    % initialize all arrays to 0
    radii = zeros(1,NT.nedge);        
end

% set radius for current trunk
radii(trunkedge) = rtrunk;

%% junction node below this edge
junc = NT.edgenodes(trunkedge,2);

if (NT.degrees(junc)==1)
    % no daughter branches to deal with
    return

elseif (NT.degrees(junc)==2)
    error('not set up to deal with deg 2 nodes yet')
    
elseif (NT.degrees(junc)==3)
   
    edge1 = NT.nodeedges(junc,2);
    edge2 = NT.nodeedges(junc,3);
    
    % set radii for daughter branches
    r1 = muvals(trunkedge)^(1/a)*radii(trunkedge);
    r2 = (1-muvals(trunkedge))^(1/a)*radii(trunkedge);
    [radii] = setRadiiFromSubtreeInfo(NT,edge1,a,r1,muvals,radii);
    [radii] = setRadiiFromSubtreeInfo(NT,edge2,a,r2,muvals,radii);    
else
    error('not set up to deal with node degrees > 3')
end



end