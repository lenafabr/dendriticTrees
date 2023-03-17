function datatable = getJunctionDataTable(allNetworks,allumperpx, rm)

%% prepare empty table

varnames = {'cell','n','r0^2','r0^2+rm^2', 'r1^2','r2^2','r1^2+r2^2','L1','L2','D1','D2','L1/D1','L2/D2','V1','V2',...
    'r^2_asym','r^2_asym^2','L_asym','L_asym^2','D_asym','D_asym^2','L/D_asym','L/D_asym^2','V_asym','V_asym^2'};

datatable = cell2table(cell(0,length(varnames)), 'VariableNames', varnames);

%%
for nc = 1:length(allNetworks)
    % go through each network
    NT = allNetworks(nc);
    umperpx = allumperpx(nc);
    
    % set trunk edge
    % (longest edge starting at root)
    rootedges = NT.nodeedges(NT.rootnode,1:NT.degrees(NT.rootnode));
    [~,ind] = max(NT.edgelens(rootedges));
    trunkedge = rootedges(ind);
    
    %% define a single radius for each edge, in um units
    % for now, this will be the average of all measured radii on the edge    
    for ec = 1:NT.nedge
        radii(ec) = mean(NT.edgewidth{ec}(:,1),1)/2*umperpx;
    end

    %% calculate subtree statistics from the defined radii
    % subtrees are indexed by trunk edge
    % this gets the subtree length stL, volume stV, and depth stD
    stL = []; stV = []; stD = [];
    [stL,stV,stD] = setSubtreeInfo_fromRadii(NT,trunkedge,stL,stV,stD,radii);
    
    %% get a list of parent & daughter edges for each nonterminal junction
    [junctedges,isterminal] = getJunctionEdges(NT);   
    junctedges = junctedges(~isterminal,:);
    
    %% save info for each nonterminal junction
    startjc = size(datatable,1);
    
    for jcur = 1:size(junctedges,1)
        jc = startjc+jcur; % index in full junction list
        edge0 = junctedges(jcur,1); % parent
        edge1 = junctedges(jcur,2); % daughter 1
        edge2 = junctedges(jcur,3); % daughter 2
        
        datatable.("cell")(jc) = nc;
        datatable.("n")(jc) = jcur;
        
        % radii
        datatable.("r0^2")(jc) = radii(edge0)^2;
        datatable.("r0^2+rm^2")(jc) = radii(edge0)^2+rm^2;
        datatable.("r1^2")(jc) = radii(edge1)^2;
        datatable.("r2^2")(jc) = radii(edge2)^2;
        datatable.("r1^2+r2^2")(jc) = radii(edge1)^2 + radii(edge2)^2;
        
        % subtree sizes
        datatable.("L1")(jc) = stL(edge1);        
        datatable.("L2")(jc) = stL(edge2);
        datatable.("D1")(jc) = stD(edge1);
        datatable.("D2")(jc) = stD(edge2);        
        datatable.("L1/D1")(jc) = stL(edge1)/stD(edge1);
        datatable.("L2/D2")(jc) = stL(edge2)/stD(edge2);
        datatable.("V1")(jc) = stV(edge1);        
        datatable.("V2")(jc) = stV(edge2);        
    end            
end

% get asymmetries for whole table
datatable.("r^2_asym") = (datatable.("r1^2") - datatable.("r2^2"))./(datatable.("r1^2") + datatable.("r2^2"));
datatable.("r^2_asym^2") = datatable.("r^2_asym").^2;
datatable.("L_asym") = (datatable.("L1") - datatable.("L2"))./(datatable.("L1") + datatable.("L2"));
datatable.("L_asym^2") = datatable.("L_asym").^2;
datatable.("D_asym") = (datatable.("D1") - datatable.("D2"))./(datatable.("D1") + datatable.("D2"));
datatable.("D_asym^2") = datatable.("D_asym").^2;
datatable.("L/D_asym") = (datatable.("L1/D1") - datatable.("L2/D2"))./(datatable.("L1/D1") + datatable.("L2/D2"));
datatable.("L/D_asym^2") = datatable.("L/D_asym").^2;
datatable.("V_asym") = (datatable.("V1") - datatable.("V2"))./(datatable.("V1") + datatable.("V2"));
datatable.("V_asym^2") = datatable.("V_asym").^2;