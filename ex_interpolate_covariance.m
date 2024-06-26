clear; close all;

addpath('../3D-dispersion-comsol') % So that design_parameters object can be loaded

% Measure the empirical covariance of a 3D dataset, and plot it.

covariance_options = struct();
covariance_options.eig_idxs = [1];
covariance_options.isComputeCovarianceGradient = false;

% data = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\IBZ output 16-May-2023 17-15-12\checkpoint38.mat");
% data = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\output 18-May-2023 18-16-56\checkpoint152 - Copy.mat");
data = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\full BZ output 18-May-2023 18-16-56\checkpoint233_ndgrid.mat");

WAVEVECTOR_DATA = data.WAVEVECTOR_DATA;
EIGENVALUE_DATA = data.EIGENVALUE_DATA(:,:,data.c.struct_idxs);

[Cs,~,kfcns,~,Vs,~,~,~,X_grid_vec,Y_grid_vec,Z_grid_vec,covariance_info] = get_empirical_covariance(WAVEVECTOR_DATA,EIGENVALUE_DATA,covariance_options);

% Look at the covariance matrix reshaped to a 2D array
figure
imagesc(reshape(Cs{1},[prod(data.c.N_wv) prod(data.c.N_wv)]))

% Look at multiple slices of the 6D covariance array long each combination
% of directions
% There is randomness associated with these images because I pick random
% slices each time, but ensure I'm always looking at *every possible* pair
% of directions.
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

clear temp

%% Interpolation of covariance
% Do covariance measurements match?
idx = 1;
C6D_temp = cov(squeeze(EIGENVALUE_DATA(:,idx,:))');
C6D_temp = reshape(C6D_temp,[data.c.N_wv data.c.N_wv]);

C6D = Cs{idx};

any(C6D_temp - C6D,'all')

kfcn = kfcns{idx};
% Generate a high resolution list of wavevectors
% N_wv_h = [2*data.c.N_wv(1) 2*data.c.N_wv(2) 2*data.c.N_wv(3)-1];
N_wv_h = data.c.N_wv;
wv_h = get_IBZ_wavevectors(N_wv_h,data.c.L,'none',1);

figure
tlo = tiledlayout(1,2);
nexttile
scatter3(WAVEVECTOR_DATA(:,1),WAVEVECTOR_DATA(:,2),WAVEVECTOR_DATA(:,3))
title('wavevectors of the original dataset that created the covariance function')
daspect([1 1 1])

nexttile
scatter3(wv_h(:,1),wv_h(:,2),wv_h(:,3))
title('wavevectors in the high resolution space that we are interpolating the covariance to')
daspect([1 1 1])

% Interpolate covariance to the high resolution set of wavevectors
% Scattered method
C_interp_scattered = kfcn(wv_h,wv_h,'scattered');

C6D_interp_scattered = reshape(C_interp_scattered,[data.c.N_wv data.c.N_wv]);

disp('scattered kfcn error')
any(C6D - C6D_interp_scattered,'all')

% Gridded method
C_interp_gridded = kfcn(wv_h,wv_h,'scattered');

C6D_interp_gridded = reshape(C_interp_gridded,[data.c.N_wv data.c.N_wv]);

disp('gridded kfcn error')
any(C6D - C6D_interp_gridded,'all')
% Interpolate the covariance against a few single wavevectors and make sure
% the results look the same.

fig = figure;
tlo = tiledlayout('flow');

% Define the wavevector pairs to look at
wv_idx = 1;
wv_i_temp = WAVEVECTOR_DATA; % I want to compute the covariance between all wavevectors and...
wv_j_temp = WAVEVECTOR_DATA(wv_idx,:); % [-pi/data.c.L -pi/data.c.L 0]; % a single other wavevector 

nexttile
C3D_temp1 = C(:,:,:,wv_idx);
scatter3(WAVEVECTOR_DATA(:,1),WAVEVECTOR_DATA(:,2),WAVEVECTOR_DATA(:,3),[],C3D_temp1(:),'filled')
daspect([1 1 1])
title(['covariance field for all wavevectors with a single other selected wavevector' newline 'wavevector = ' num2str(wv_j_temp,3)])
colorbar

nexttile
C3D_temp2 = kfcn(wv_i_temp,wv_j_temp,'scattered');
scatter3(WAVEVECTOR_DATA(:,1),WAVEVECTOR_DATA(:,2),WAVEVECTOR_DATA(:,3),[],C3D_temp2,'filled')
daspect([1 1 1])
title(['covariance field for all wavevectors with a single other selected wavevector' newline 'wavevector = ' num2str(wv_j_temp,3)])
colorbar

nexttile
C3D_temp_diff_ratio = (C3D_temp1(:)-C3D_temp2(:))./C3D_temp1(:);
scatter3(WAVEVECTOR_DATA(:,1),WAVEVECTOR_DATA(:,2),WAVEVECTOR_DATA(:,3),[],C3D_temp_diff_ratio,'filled')
daspect([1 1 1])
title(['covariance field for all wavevectors with a single other selected wavevector' newline 'wavevector = ' num2str(wv_j_temp,3)])
colorbar

% Interpolate the covariance and make sure it still looks good.


% C_interp = kfcn()

% Idea - compute the interpolated covariance field for the covariance between
% a higher resolution grid and a single wavevector, and check that this
% field is the same for a variety of different wavevectors. That way I can
% actually just plot the entire covariance field since it can actually be
% thought of in 3D.