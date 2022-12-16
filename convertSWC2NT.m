function [NT,mergednodes] = convertSWC2NT(exTree)
% convert an SWC tree structure, extracted from a .swc file using load_tree
% into a network object structure
% all degree 2 nodes get merged away and long edges saved as edgepath
% edge widths are stored for all intermediate nodes and downstream node
% along edge
% NOTEs: 
% assumes node 1 is the parent node
% the order of edges in the resulting network gets mixed around

%% convert to network structure
NT=NetworkObj;
NT.nodepos=[exTree.X, exTree.Y, exTree.Z];
% find reads column first
% if A is symmetric matrix, this will store duplicate edges
[n1,n2,~] = find(exTree.dA);
% sort connectivity based on node order
NT.edgenodes = sortrows([n1,n2]')';
NT.setupNetwork;

% store widths, keeping the node downstream of each edge 
for ec = 1:NT.nedge
    n1 = NT.edgenodes(ec,1);
    n2 = NT.edgenodes(ec,2);
    NT.edgewidth{ec}(1,1) = exTree.D(n2);
    %NT.edgewidth{ec}(2,1) = exTree.D(n2);
    % store contour length
    NT.edgewidth{ec}(1,2) = 0.99*NT.edgelens(ec);
   % NT.edgewidth{ec}(2,2) = 0.99*NT.edgelens(ec);
    % store end points; XY of 1 endpoint, XY of 2nd endpoint
end

% order the tree
NT.rootnode = 1;
isset = false(NT.nedge,1);
wasreversed = false(NT.nedge,1);
[isset,wasreversed] = directedTreeEdges(NT,NT.rootnode,isset,wasreversed);

%% Filter out degree 2 nodes along all edges, keeping width info
mergednodes = NT.mergeAllEdgePaths(true);

% order the tree again just in case
NT.rootnode = 1;
isset = false(NT.nedge,1);
wasreversed = false(NT.nedge,1);
[isset,wasreversed] = directedTreeEdges(NT,NT.rootnode,isset,wasreversed);




