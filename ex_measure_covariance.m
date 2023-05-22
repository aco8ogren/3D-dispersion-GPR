clear; close all;

addpath('../3D-dispersion-comsol') % So that design_parameters object can be loaded

% Measure the empirical covariance of a 3D dataset, and plot it.

covariance_options = struct();
covariance_options.eig_idxs = [1];
covariance_options.isComputeCovarianceGradient = false;

data = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\IBZ output 16-May-2023 17-15-12\checkpoint38.mat");

WAVEVECTOR_DATA = data.WAVEVECTOR_DATA;
EIGENVALUE_DATA = data.EIGENVALUE_DATA;

[Cs,~,kfcns,~,Vs,~,~,~,X_grid_vec,Y_grid_vec,covariance_info] = get_empirical_covariance(WAVEVECTOR_DATA,EIGENVALUE_DATA,covariance_options);

figure
imagesc(reshape(Cs{1},[prod(data.c.N_wv) prod(data.c.N_wv)]))

figure
idxs = 1:6;
axis_pairs = combvec(idxs,idxs);
axis_pairs(:,axis_pairs(1,:)==axis_pairs(2,:)) = [];

fig = figure;
tlo = tiledlayout('flow');

C = Cs{1};
for i = 1:size(axis_pairs,2)
    for j = 1:length(size(C))
        r = randi(size(C,j));
        indices{j} = r;
    end
    for k = 1:2
        indices{axis_pairs(k,i)} = 1:size(C,axis_pairs(k,i));
    end
    temp = squeeze(C(indices{:}));
    
    nexttile
    imagesc(temp)
    colorbar
    % title(axis_pairs(:,i)')
    xlabel(['axis ' num2str(max(axis_pairs(:,i)))])
    ylabel(['axis ' num2str(min(axis_pairs(:,i)))])
end

title(tlo,'Examining the smoothness of the 6D covariance array along each pair of slice directions')
