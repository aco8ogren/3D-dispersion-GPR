clear; close all;

data_path = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\full BZ output 18-May-2023 18-16-56\checkpoint233_ndgrid.mat";
data = load(data_path);

e = data.EIGENVALUE_DATA;
w = data.WAVEVECTOR_DATA;

% Find index for which the wavevector is [0 0 0]
[~,index] = ismember([0 0 0],w,'rows');

% w(index,:)

figure
plot(squeeze(e(index,1,data.c.struct_idxs)))