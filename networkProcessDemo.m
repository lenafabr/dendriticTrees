%% clean up network and get edge widths
% set up file paths to add network tools and gui
% change this to match the path where you put this repository
% network tools repository should be cloned from 
% https://github.com/lenafabr/networktools

addpath('../../networktools')
addpath('../../networktools/gui')

%% if network has been previously processed, load in the workspace
% and skip the network cleaning steps below
load('./example/cleanNetwork_MCFO-HSE-1.mat','NT','origimgfile','img','umperpx')

% ------------------ EDITING AND MEASURING NETWORK ---------------
%% Load image files
% original image file, not segmented, will be used for visualization only
origimgfile = './example/MCFO-HSE-1.tif'; %simpler with no loops

% single-frame image file for a black and white image (can be skeletonized
% already, or merely segmented)
bwimgfile = './example/MCFO-HSE-1_skeleton.tif';

% load the images
img = imread(origimgfile);
bwimg = imread(bwimgfile);

% in case of rgb - use relevant channel
cno = 1;
img = img(:,:,cno);


%% Extract a network object structure from the black and while image
% this step can be  slow for a large image and has not been fully optimized
NT= getNetworkFromBWImage(bwimg);
% make a backup copy of the network before you edit it
NT0 = copy(NT);

%% visualize the network object to make sure it looks ok
% plotopt: specify linewidth, color etc for network skeleton

% plotting options (avoids drawing black lines on black background)
plotopt = struct('plotoverimage',true,'datatipindex',true);

figure(2)
imshow(img,[])
%imshow(img(:,:,:))
hold all
NT.plotNetwork(plotopt)
hold off
%% save existing workspace before doing anything to the network
save('../imgData/mysavefile.mat','NT','img','plotopt')

%% Edit network structure semi-manually using a GUI
% use the GUI to do the following:
% 1) remove loops (hit Show Loops button to visualize)
% 2) merge edges at degree-2 nodes (to do all: set Merge all checkbox, then hit
% Merge at Selected Nodes)
% 3) clear out high degree nodes (hit Color by Degree to see them)
% 4) Set the root node of the tree (this will update edge directions into a
% directed tree)
% 5) Make sure the root node has degree 1!
% 6) Can also measure edge widths at the same time or later. Measurements
% can later be adjusted.
% 7) **Make sure you hit "Update Network" before leaving GUI**

% GUI to run network edit code - may want to save existing NT workspace first!
networkEdit('NT',NT,'img',img,'plotopt',plotopt)

%% after editing is a good time to save your work

% specify in how many microns per pixels for our cell (to make real length units)
umperpx = 0.30145
save('../imgData/mysavefile_cleaned.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% START HERE IF ALREADY HAVE CLEANED UP NETWORK STRUCTURE %%%%
% Network must be directed, have no loops, and no degree 2 nodes.

%% load a workspace and set scale factor for convertion to um units
load('./example/cleanNetwork_MCFO-HSE-1.mat','NT','origimgfile','img','umperpx')


%% Check that network has been properly cleaned.
% set trunk edge 
% (longest edge starting at root)
rootedges = NT.nodeedges(NT.rootnode,1:NT.degrees(NT.rootnode));
[~,ind] = max(NT.edgelens(rootedges));
trunkedge = rootedges(ind);

disp('Running checks on network structure')

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

%% reset any missing edge widths based on the parent edge
% WARNING: only do this if you want to include fake edge widths
% goes down the tree starting with the trunk
% will set edge width of everything to default if trunk width is not set
% these widths can be adjusted in the GUI
defaultWidth = 10;
propagateEdgeWidths(NT,NaN,defaultWidth)

%% define a single radius for each edge, in um units
% for now, this will be the average of all measured radii on the edge
for ec = 1:NT.nedge
    radii(ec) = mean(NT.edgewidth{ec}(:,1),1)/2*umperpx;
end

%% calculate subtree statistics from the defined radii
% subtrees are indexed by trunk edge
stL = []; stV = []; stD = [];
[stL,stV,stD] = setSubtreeInfo_fromRadii(NT,trunkedge,stL,stV,stD,radii)

%% for the network with measured edge widths, calculate linear mitochondria concentrations

% default linear concentration of walking mitos at trunk
rhoWtrunk = 1; 
% beta is the exponent relating stopping rate and trunk width: ks ~ 1/r^beta
beta = 1.3;
% kskw is the ratio of stopping to walking rate for a radius of 1
kskw = 10;

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

% get the average asymmetry
avgasym = sqrt(sum(asymmetry(nontermjunc).^2)/length(nontermjunc))

%% Calculate distal enrichment of mitochondrial density

% first define what we consider distal edges
% distal = all edges with path length > this factor * max path length (from
% root node)
cutoff = 0.75;

% distance of each edge from root node
[edgedist] = getEdgeDist(NT,trunkedge);
edgedist = edgedist*umperpx;
maxdist = max(edgedist);
distind = find(edgedist > cutoff*maxdist);

Mdist = sum(rhovals(distind).*NT.edgelens(distind));
Vdist = sum(NT.edgelens(distind).*radii(distind).^2);

% distal enrichment: average mito volume density in all distal edges put together
% compared to volume density in the trunk
distenrich = (Mdist/Vdist) / (rhovals(trunkedge)/radii(trunkedge).^2)

%% ------------------ Generate synthetic radii if desired -----------------------
% If you want to determine synthetic radii for the tree edges:
clear stL stEta stD muvals
a = 2; % alpha value relating parent and daughter width: r0^a = r1^a + r2^a
rm = 0; % `minimal radius' (from Liao paper) such that r0^a + rm^a = r1^a + r2^a

% how to split radii. Possible values:
% LV = split so that subtree L ~ V; does *not* use rm
% LD = split so that branch radius ~ subtree length / depth
% L = split so that branch radius ~ subtree length
% equal = split equally
splittype = 'LD';


% define radius of trunk
rtrunk = NT.edgewidth{trunkedge}(1,1);


% calculate some basic statistics
[stL,stEta,stD,muvals] = setSubtreeInfo(NT,trunkedge,a,'LV');        

switch(splittype)
    case('LV')
        radii = setRadiiFromSubtreeInfo(NT,trunkedge,a,rtrunk,muvals);   
        stV = stEta.*radii.^2;
    case('LD')
        [radii,stV] = setRadiiWithRm(NT,trunkedge,a,rm,rtrunk,stL./stD);
    case('equal')       
        [radii,stV] = setRadiiWithRm(NT,trunkedge,a,rm,rtrunk,ones(size(stL)));
    case('L')        
        [radii,stV] = setRadiiWithRm(NT,trunkedge,a,rm,rtrunk,stL);
end

% after setting synthetic radii, can go back to "calculate mitochondria
% concentrations", get asymmetry and distal enrichment

%% to visualize the radii in gui, store fake width info in tree
for ec = 1:NT.nedge
    NT.edgewidth{ec} = makeEdgeWidth(NT,ec,radii(ec)/umperpx);
end


% --------------------
%% Plot network with widths and mito densities
% --------------------

figure

% colors to scale between
color2 = [0.2, 0.8,0.8];     
color1 = [0,  0,  0]; 
% scaling factor for edge widths
sclwidth = 0.2;

nedge = NT.nedge;
cmap = colormapinterp([color1;color2],nedge);
maxconc = max(rhovals);
colormap(cmap);

for e = 1:NT.nedge
    if (radii(e)>0)        
        avgWidth = 2*radii(e)*umperpx*sclwidth;    
        conc = rhovals(e);

        % set edge color according to conc
        if (conc>maxconc)
            cinterp = color2;
        else
            cinterp = interp1(linspace(0,maxconc,size(cmap,1)),cmap,conc);       
        end
        savecinterp(e,:) = cinterp;
        
        %plot(NT.edgepath{e}(:,1),NT.edgepath{e}(:,2),'color',(alpha*color2) + ((1-alpha)*color1),'Linewidth',(avgWidth+0.01)/1.5)
        plot(NT.edgepath{e}(:,1),NT.edgepath{e}(:,2),'color',cinterp,'Linewidth',(avgWidth+0.01)*10)
        hold on       
    end
end

colormap(cmap); caxis([0,10])
set(gca,'Visible','off','FontSize',16)
h=colorbar;
set(get(h,'ylabel'),'string','Mito Density');
hold off

        