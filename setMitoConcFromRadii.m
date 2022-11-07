function [rhoSvals,rhoWvals,stM] = setMitoConcFromRadii(NT,trunkedge,beta,kskw,rhoWtrunk,radii,rhoSvals,rhoWvals,stM)
% walk down a tree, starting with the trunkedge
% using precalculated radii (from setRadiiFromSubtreeInfo)
% set the linear concentrations of mitochondria (rhovals)
% set density of stationary mitos rhoS and walking mitos rhoW
% kskw = ks^*/kw prefactor 
% stopping rate = ks*/r^beta
% a = alpha for parent-daughter radius scaling
% ** assuming walking mitos split in proportion to r^2 **
% rhoSvals= stationary mito linear density in each edge
% rhoWvals = walking mito linear density in each edge
% stM = total mitochondria in subtree stemming from each edge.


% if just starting the recursion,
% initialize arrays to 0 
if (~exist('rhoSvals','var'))    
    rhoSvals = zeros(1,NT.nedge);        
end
if (isempty(rhoSvals))
    % initialize all arrays to 0
    rhoSvals = zeros(1,NT.nedge);        
end
if (~exist('rhoWvals','var'))    
    rhoWvals = zeros(1,NT.nedge);        
end
if (isempty(rhoWvals))
    % initialize all arrays to 0
    rhoWvals = zeros(1,NT.nedge);        
end
if (~exist('stM','var'))    
    stM = zeros(1,NT.nedge);        
end
if (isempty(stM))
    % initialize all arrays to 0
    stM = zeros(1,NT.nedge);        
end


% set rho for current trunk
rhoWvals(trunkedge) = rhoWtrunk;
rhoSvals(trunkedge) = rhoWtrunk*kskw/radii(trunkedge)^beta;



%% junction node below this edge
junc = NT.edgenodes(trunkedge,2);

if (NT.degrees(junc)==1)
    % no daughter branches to deal with
    % just set the total mitos in this branch
    stM(trunkedge) = (rhoWvals(trunkedge)+rhoSvals(trunkedge))*NT.edgelens(trunkedge);
    return

elseif (NT.degrees(junc)==2)
    error('not set up to deal with deg 2 nodes yet')
    
elseif (NT.degrees(junc)==3)
   
    edge1 = NT.nodeedges(junc,2);
    edge2 = NT.nodeedges(junc,3);
    
    % get linear densities for daughter branches
    r1 = radii(edge1); r2 = radii(edge2);
    rhoW1 = r1^2/(r1^2+r2^2)*rhoWvals(trunkedge);
    rhoW2 = r2^2/(r1^2+r2^2)*rhoWvals(trunkedge);
    
    [rhoSvals,rhoWvals,stM] = setMitoConcFromRadii(NT,edge1,beta,kskw,rhoW1,radii,rhoSvals,rhoWvals,stM);
    [rhoSvals,rhoWvals,stM] = setMitoConcFromRadii(NT,edge2,beta,kskw,rhoW2,radii,rhoSvals,rhoWvals,stM);    
    
    % get total mitochondria in full tree
    stM(trunkedge) = stM(edge1) + stM(edge2) + (rhoWvals(trunkedge)+rhoSvals(trunkedge))*NT.edgelens(trunkedge);
    
    
else
    error('not set up to deal with node degrees > 3')
end



end