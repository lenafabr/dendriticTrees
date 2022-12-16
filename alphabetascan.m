%% directory name containing multiple .mat files
% each .mat file needs to contain the following saved variables:
% NT = a network object, already arranged to be a directed tree, rootnode set
% umperpx = conversion factor, micrometers per pixel
dirname = '../test/';

%% load in all .mat files within the directory
% and save data
files = dir([dirname '*.mat'])
filenames = {files.name};

clear allNetworks
for fc = 1:length(filenames)
        
    fname = [dirname filenames{fc}];
    disp(fname)
    load(fname);
    
    allNetworks(fc) = NT;    
    allumperpx(fc) = umperpx;
end


%% Go through and calculate asymmetries and enrichment for different alpha & beta values

rm = 0; % `minimal radius' (from Liao paper) such that r0^a + rm^a = r1^a + r2^a
% how to split radii. Possible values:
% LV = split so that subtree L ~ V; does *not* use rm
% LD = split so that branch radius ~ subtree length / depth
% L = split so that branch radius ~ subtree length
% equal = split equally
splittype = 'LD';

% list of alpha values to use
alphalist = linspace(1,3,20);
% list of beta values to use
betalist = linspace(0,3,19);

rtrunk = 3; % trunk radius in um
% ratio of stopped to walking rates
kskw = 1000;
% density of mitos in trunk (arbitrary units)
rhoWtrunk = 1;

% distal edges are defined as > this prefactor times max distance to root
distalprefactor= 0.75;

clear avgasym distenrich
for nc = 1:length(allNetworks)    
    NT = allNetworks(nc);

    trunkedge = NT.nodeedges(NT.rootnode);
        
    %% define which are the distal edges
     % distance of each edge from root node
     [edgedist] = getEdgeDist(NT,trunkedge);
     edgedist = edgedist*allumperpx(nc);
    
     % distal is > 75% of max distance   
     distcutoff(1) = max(edgedist)*distalprefactor;
     distcutoff(2) = inf;
     distind = find(edgedist>distcutoff(1) & edgedist<distcutoff(2));  

     for ac = 1:length(alphalist)
        [nc ac]
        
        a = alphalist(ac);
        %% generate synthetic radii based on this value of alpha
        % for the selected rule of splitting between sisters (splittype)
        
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

        % loop go over different scalings for mitochondrial stopping rate
        for bc = 1:length(betalist)
            beta = betalist(bc);
            % calculate linear mitochondria concentrations
            
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
            avgasym(ac,bc,nc) = sqrt(sum(asymmetry(nontermjunc).^2)/length(nontermjunc));

            %% Calculate distal enrichment of mitochondrial density
            Mdist = sum(rhovals(distind).*NT.edgelens(distind));
            Vdist = sum(NT.edgelens(distind).*radii(distind).^2);
            distenrich(ac,bc,nc) = (Mdist/Vdist) / (rhovals(trunkedge)/radii(trunkedge).^2);
            
        end
     end
end

%% take the average over all networks
allavgasym = squeeze(mean(avgasym,3));
allavgenrich = squeeze(mean(distenrich,3));

%% plot asymmetries
cmap = BBVYWcolormap(100);
colormap(cmap)
pcolor(betalist,alphalist,allavgasym)
shading flat
plot_cleanup(gca,'FontSize',14,'pcolor',true)
xlabel('transport scaling $\beta$')
ylabel('parent-daughter scaling $\alpha$')

title('Average Asymmetry','Interpreter','latex')

% add alpha confidence interval from measured values
amin = 1.99; amax = 2.54;
hold all
plot(betalist,amin*ones(size(betalist)),'g--','LineWidth',2)
plot(betalist,amax*ones(size(betalist)),'g--','LineWidth',2)
hold off

ylim([1.3,3])
caxis([0,0.4])

%% plot distal enrichment
cmap = BBVYWcolormap(100);
colormap(cmap)
pcolor(betalist,alphalist,allavgenrich)
shading flat
plot_cleanup(gca,'FontSize',14,'pcolor',true)
xlabel('transport scaling $\beta$')
ylabel('parent-daughter scaling $\alpha$')
title('Distal enrichment','Interpreter','latex')
ylim([1.3,3])
caxis([1e-1,1e3])
set(gca,'ColorScale','log')

% add alpha confidence interval
amin = 1.99; amax = 2.54;
hold all
plot(betalist,amin*ones(size(betalist)),'g--','LineWidth',2)
plot(betalist,amax*ones(size(betalist)),'g--','LineWidth',2)
hold off

