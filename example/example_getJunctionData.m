
addpath('~/proj/networktools/')
addpath('../../treestoolbox-master');
addpath('../../dendriticTrees_public/')

%% load in all .mat files within current directory
% and save data
dirname = ['./']
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


%% get table of data for each junction
rm = 0; % assume rm=0;
datatable = getJunctionDataTable(allNetworks,allumperpx, rm);

%% output data table
writetable(datatable,'../../examples/example_datatable.csv')

%% plot L and V asymmetries

xvals = datatable.("L_asym");
yvals = datatable.("V_asym");

plot(xvals,yvals,'.')

mdl = fitlm(xvals,yvals)

%% Plot r^2 and L/D asymmetries
figure
xvals = datatable.("L/D_asym");
yvals = datatable.("r^2_asym");

plot(xvals,yvals,'.')
xlabel('L/D asymmetry')
ylabel('r^2 asymmetry')

mdl = fitlm(xvals,yvals)
mdl.plot()

%% Plot r^2 and L asymmetries
figure
xvals = datatable.("L_asym");
yvals = datatable.("r^2_asym");

plot(xvals,yvals,'.')
xlabel('L asymmetry')
ylabel('r^2 asymmetry')

mdl = fitlm(xvals,yvals)
mdl.plot()
%% Plot r^2 and D asymmetries
figure
xvals = datatable.("D_asym");
yvals = datatable.("r^2_asym");

plot(xvals,yvals,'.')
xlabel('D asymmetry')
ylabel('r^2 asymmetry')

mdl = fitlm(xvals,yvals)
mdl.plot()