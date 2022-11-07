function cmap = colormapinterp(nodes,npt)
% create a colormap by interpolating between defined node colors
% nodes is an nx3 array of colors
% npt is the number of points for the colormap

nn = size(nodes,1);

cmap = interp1(1:nn,nodes,linspace(1,nn,npt))

end