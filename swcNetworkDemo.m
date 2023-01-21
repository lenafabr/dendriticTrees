% update the path names here to match your system
addpath(genpath('../treestoolbox-master/'))
addpath('../../networktools')

%% load tree data from an SWC file
swcfile = './example/1.VS-brain1_left_220811.swc';

[exTree, name, path] = load_tree(swcfile);


%% convert to network object, including saving edge widths along edges
% NOTEs: 
% - node 1 is assumed to be the tree root in exTree
% - the edge order will get shuffled around in this process
% - all degree 2 nodes are removed and turned into edge paths
% 
% - mergednodes lists, for each edge in NT, what were the original intermediate node
% indices in exTree that got removed
%
% NT.edgewidth keeps width measurements for all intermediate nodes and the
% downstream end-node for each edge
% 1st column is radius, 2nd column is position along edge
[NT,mergednodes] = convertSWC2NT(exTree);

%% Check that network has been properly cleaned.
% set trunk edge
% (longest edge starting at root)
rootedges = NT.nodeedges(NT.rootnode,1:NT.degrees(NT.rootnode));
[~,ind] = max(NT.edgelens(rootedges));
trunkedge = rootedges(ind);

if (NT.degrees(NT.rootnode)>1)
    warning('Root node has degree > 1. May break code!')
end

if (any(NT.degrees>3))
    warning('Network still has high-degree nodes. May break code')
end
if (any(NT.degrees==2))
    warning('Network still has degree 2 nodes. May break code')
end
%% check which edges are missing a width measurement
missingedges = [];
for ec = 1:NT.nedge
    if isempty(NT.edgewidth{ec})
        disp(sprintf('Missing width measurement: %d', ec))
        missingedges = [missingedges ec];
    end
end
if (isempty(missingedges))
    disp('No edge widths are missing.')
end

%% define a single radius for each edge, in um units
umperpx = 1;
% for now, this will be the average of all measured radii on the edge
for ec = 1:NT.nedge %don't need to run if creating synthetic ones
    radii(ec) = mean(NT.edgewidth{ec}(:,1),1)/umperpx;
end

%% calculate subtree statistics from the defined radii
% subtrees are indexed by trunk edge
stL = []; stV = []; stD = [];
[stL,stV,stD] = setSubtreeInfo_fromRadii(NT,trunkedge,stL,stV,stD,radii);

%% for the network with measured edge widths, calculate linear mitochondria concentrations
%Ask Erin what she wants for these values
% default linear concentration of walking mitos at trunk
rhoWtrunk = 1;
% kskw is the ratio of stopping to walking rate for a radius of 1
kskw = 10;
beta = 1.3

trunkedge = NT.nodeedges(NT.rootnode,1);
% rhoSvals = stationary density on each edge
% rhoWvals = mobile density on each edge
% stM = average mito density in the subtree below (and including) each edge
[rhoSvals,rhoWvals,stM] = setMitoConcFromRadii(NT,trunkedge,beta,kskw,rhoWtrunk,radii);
rhovals = rhoSvals + rhoWvals;

%% Calculate asymmetry in mito concentrations below each junction
% calculate asymmetry at each junction
% indexed by edge leading to that junction
[asymmetry,nontermjunc] = getAsymmetry(NT,stM, stV);

%% get rms asymmetry in predicted mito densities:
rmsasymmetry = sqrt(mean(asymmetry.^2))