clear; close all;

addpath('../3D-dispersion-comsol')

data_fn = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\full BZ output 18-May-2023 18-16-56\checkpoint233.mat";

data = load(data_fn);

vars_to_save = fields(data);

unpack_struct(data);

[WAVEVECTOR_DATA,indices] = sortrows(WAVEVECTOR_DATA,[3 2 1]);
EIGENVALUE_DATA = EIGENVALUE_DATA(indices,:,:);

save(data_fn,vars_to_save{:},'-v7.3')