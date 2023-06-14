clear; close all;
%% Instructions
% =========================================================================
% MUST BE RUN IN MATLAB R2020a (or maybe a later version would be okay too)
% =========================================================================

%% Turn off warnings
warning('off','MATLAB:MKDIR:DirectoryExists')
% warning('off','MATLAB:handle_graphics:Layout:NoPositionSetInTiledChartLayout')

%% Add directories
% addpath('../2D-dispersion-optimization-GPR')
addpath('../3D-dispersion-comsol')

%% Flags
% isSavePlots = false;
% isUseHomemade = true;
% isHighlightFirstStructure = true;
% isMakeBoxPlots = false;
% isMakeQuantilePlots = true;
isPlot = true;
isSaveData = true;

%% Settings
% Iteratable variables
eig_idx = 10;
struct_idx = 1; % which struct_idxs to test in the test set
model_name = 'GPR'; % {'GPR','linear','nearest','cubic','makima','spline'};
N_sample = 20; %round(logspace(log10(5),log10(1000),10));
sigma = [0]; % If not a scalar value or NaN, any non-GPR models will be wastefully run multiple times.
training_set_size = 234; % [50 83 139 232 387 646 1077 1797 2997 5000];

data_path_train = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\IBZ output 18-May-2023 18-16-56\checkpoint233_ndgrid.mat";
data_path_train = char(data_path_train);

sample_order_path = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\IBZ output 18-May-2023 18-16-56\sample_order_data.mat";
sample_order_path = char(sample_order_path);

sample_order_data = load(sample_order_path);
N_sample_max = size(sample_order_data.sample_orders,1); % kind of correct?

regexp_idx = regexp(data_path_train,'\');
data_dir = data_path_train(1:(regexp_idx(end)));
script_start_time = replace(char(datetime),':','-');

% Load training set
data_train = load(data_path_train,'EIGENVALUE_DATA','WAVEVECTOR_DATA');
EIGENVALUE_DATA_TRAIN = data_train.EIGENVALUE_DATA;
WAVEVECTOR_DATA_TRAIN = data_train.WAVEVECTOR_DATA;

[N_wv_train,N_eig_train,N_struct_train] = size(EIGENVALUE_DATA_TRAIN);

data_path_test = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\output 29-May-2023 17-01-51\checkpoint15 - Copy.mat";
data_path_test = char(data_path_test);

% Load test set
data_test = load(data_path_test,'EIGENVALUE_DATA','WAVEVECTOR_DATA');
EIGENVALUE_DATA_TEST = data_test.EIGENVALUE_DATA;
WAVEVECTOR_DATA_TEST = data_test.WAVEVECTOR_DATA;

N_evaluate = [numel(unique(WAVEVECTOR_DATA_TEST(:,1))) numel(unique(WAVEVECTOR_DATA_TEST(:,2))) numel(unique(WAVEVECTOR_DATA_TEST(:,3)))];

[N_wv_test,N_eig_test,N_struct_test] = size(EIGENVALUE_DATA_TEST);

covariance_options.isAllowGPU = false;
covariance_options.isComputeCovarianceGradient = false;

%% Assertions
assert(N_sample <= N_sample_max)

%% Initializations
original_wv_x = unique(sort(WAVEVECTOR_DATA_TRAIN(:,1)));
original_wv_y = unique(sort(WAVEVECTOR_DATA_TRAIN(:,2)));
original_wv_z = unique(sort(WAVEVECTOR_DATA_TRAIN(:,3)));

covariance_options.eig_idxs = eig_idx;
[Cs,C_grads,kfcns,kfcn_grads] = get_empirical_covariance(WAVEVECTOR_DATA_TRAIN(:,:),EIGENVALUE_DATA_TRAIN(:,:,1:training_set_size),covariance_options); %#ok<ASGLU>
model_options = struct();
model_options.kfcn = kfcns{1};
model_options.kfcn_grad = {};
model_options.N_sample = N_sample;

fr = squeeze(EIGENVALUE_DATA_TEST(:,eig_idx,struct_idx));
wv = squeeze(WAVEVECTOR_DATA_TEST(:,:,struct_idx));

model_options.model_name = model_name;
if strcmp(model_name,'GPR')
    model_options.kfcn = kfcns{1};
    model_options.sample_points = sample_order_data.sample_orders(:,:,eig_idx); % the eig_idx_idx is wrong here, should align the eig_idx within this script with the eig_idx in sampling_order_data
    if isnan(model_options.sample_points(N_sample,1))
        error('why is there a nan?')% break
    end
end

model_options.sigma = sigma;
out = interp_model_2D_smart_sample(fr,wv,model_options);
% err_L2(eig_idx_idx,struct_idx_idx,model_idx,N_sample_idx,sigma_idx,training_set_size_idx) = out.e_L2;
% err_H1(eig_idx_idx,struct_idx_idx,model_idx,N_sample_idx,sigma_idx,training_set_size_idx) = out.e_H1;

% err_L2_rel(eig_idx_idx,struct_idx_idx,model_idx,N_sample_idx,sigma_idx,training_set_size_idx) = out.e_L2_rel;
% err_H1_rel(eig_idx_idx,struct_idx_idx,model_idx,N_sample_idx,sigma_idx,training_set_size_idx) = out.e_H1_rel;

disp(['e_L2 = ' num2str(out.e_L2)])
disp(['e_H1 = ' num2str(out.e_H1)])

figure
tlo = tiledlayout(3,5,'TileIndexing','columnmajor');
tlo.TileSpacing = 'compact';
slice_idxs = [1 3 5 7 10];
for i = 1:length(slice_idxs)
    s = slice_idxs(i);
    for j = 1:3
        ind = sub2ind([3,5],j,i);
        ax(i,j) = nexttile(tlo,ind);
        if j == 1
            imagesc(out.V_e(:,:,s))
            title([ 'True - slice ' num2str(s)])
        elseif j == 2
            imagesc(out.V_pred(:,:,s))
            title([ 'Pred - slice ' num2str(s)])
        elseif j == 3
            imagesc(out.V_err(:,:,s))
            title([ 'Error - slice ' num2str(s)])
        end
        daspect([1 1 1])
        colorbar
    end
end

uniform_colorscale(ax(:,1:2))

% if isSaveData
%     q = linspace(0,1,101);
%     Q_L2 = quantile(err_L2,q,2);
%     Q_H1 = quantile(err_H1,q,2);
%     Q_L2_rel = quantile(err_L2_rel,q,2);
%     Q_H1_rel = quantile(err_H1_rel,q,2);
%     save('error_analysis_data',...
%         'q','Q_L2','Q_H1','err_L2','err_H1','Q_L2_rel','Q_H1_rel','err_L2_rel','err_H1_rel',...
%         'eig_idxs','struct_idxs','model_names','num_sample_points','sigmas','training_set_sizes',...
%         'data_path_train','data_path_test')
% end


% Plot the qth error quantile for each model & eigensurface
% q = .95;
% Q_L2_temp = quantile(err_L2,q,2);
% Q_H1_temp = quantile(err_H1,q,2);
% cm = lines(length(eig_idxs));
% mm = {'s','d','^','x','*','o'};
% if isPlot
%     clear p
%     for err_idx = 1:2
%         if err_idx == 1
%             Q = Q_L2_temp;
%             err_name = 'L2';
%         elseif err_idx == 2
%             Q = Q_H1_temp;
%             err_name = 'H1';
%         end
%         for sigma_idx = 1:length(sigmas)
%             sigma = sigmas(sigma_idx);
% 
%             for model_idx = 1:N_model
%                 figure
%                 hold on
%                 for eig_idx_idx = 1:length(eig_idxs)
%                     eig_idx = eig_idxs(eig_idx_idx);
%                     if length(num_sample_points)>1
%                         p_temp = plot(num_sample_points,squeeze(Q(eig_idx_idx,:,model_idx,:,sigma_idx,:)),'color',cm(eig_idx_idx,:),'marker',mm{model_idx});
%                     elseif length(training_set_sizes)>1
%                         p_temp = plot(training_set_sizes,squeeze(Q(eig_idx_idx,:,model_idx,:,sigma_idx,:)),'color',cm(eig_idx_idx,:),'marker',mm{model_idx});
%                     end
%                     p(eig_idx_idx) = p_temp(1);
%                     p(eig_idx_idx).DisplayName = [model_names{model_idx} ' eig\_idx = ' num2str(eig_idx)];
%                 end
%                 set(gca,'yscale','log')
%                 legend(reshape(p,[],1),'location','northeast')
%                 title([err_name ' error quantiles ' regexprep(num2str(q),' +',',') newline model_names{model_idx} ' \sigma = 1e' num2str(log10(sigma))])
%                 xlabel('prod(N\_sample)')
%                 ylabel([err_name ' error'])
%                 grid on
%                 grid minor
%             end
%         end
%     end
% end

