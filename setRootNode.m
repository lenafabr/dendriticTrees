%% find a node with outgoing edges only and set it to be root
function rootnode = setRootNode(NT)
% set the root node of a tree as one that has an outgoing edge but no
% incoming ones
% WARNING: this assumes the tree has edges oriented to point along a
% directed tree
for nc = 1:NT.nnode
    % is there an outgoing edge from this node?
    hasout = any(NT.edgenodes(:,1)==nc);
    % is there an incoming edge to this node?
    hasin = any(NT.edgenodes(:,2)==nc);
    if (hasout & ~ hasin)
        % disp(sprintf('Root node: %d',nc))
        rootnode = nc;
        NT.rootnode = nc;
    end
end

