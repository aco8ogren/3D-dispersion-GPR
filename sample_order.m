clear; close all;

warning('off','MATLAB:nearlySingularMatrix');

% data_path_train = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\full BZ output 18-May-2023 18-16-56\checkpoint152 - Copy.mat";
data_path_train = "C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\full BZ output 18-May-2023 18-16-56\checkpoint233_ndgrid.mat";

N_sample = 20;
disp_idxs = 'all'; % typically set to 'all';
% N_evaluate = [31 NaN]; N_evaluate(2) = ceil(N_evaluate(1)/2);
sigma = 0;
rcond_tol = 1e-16;

isPause = true;
isPlot = true;

disp('Loading training set...')
% [WAVEVECTOR_DATA,EIGENVALUE_DATA] = load_dispersion_dataset(data_path_train);
data = load(data_path_train);
disp('done.')

if strcmp(disp_idxs,'all')
    disp_idxs = data.c.struct_idxs;
    % disp_idxs = 1:size(data.EIGENVALUE_DATA,3); 
end

EIGENVALUE_DATA = data.EIGENVALUE_DATA(:,:,disp_idxs);
WAVEVECTOR_DATA = data.WAVEVECTOR_DATA(:,:);
% [~,idxs] = sort(WAVEVECTOR_DATA(:,2));
% EIGENVALUE_DATA = EIGENVALUE_DATA(idxs,:,:);
% WAVEVECTOR_DATA = WAVEVECTOR_DATA(idxs,:);
N_evaluate = data.c.N_wv;

N_eig = size(EIGENVALUE_DATA,2);

N_wv_train(1) = numel(unique(WAVEVECTOR_DATA(:,1)));
N_wv_train(2) = numel(unique(WAVEVECTOR_DATA(:,2)));
N_wv_train(3) = numel(unique(WAVEVECTOR_DATA(:,3)));

covariance_options.isComputeCovarianceGradient = false;
covariance_options.isAllowGPU = false;
% eig_idxs = 1:data.const.N_eig;
eig_idxs = 1;

[X_e,Y_e,Z_e] = ndgrid(linspace(-pi,pi,N_evaluate(1)),linspace(-pi,pi,N_evaluate(2)),linspace(0,pi,N_evaluate(3)));
wv_e = [reshape(X_e,[],1) reshape(Y_e,[],1) reshape(Z_e,[],1)];
% [~,idxs] = sort(wv_e(:,1));
% wv_e = wv_e(idxs,:);

figure
tiledlayout('flow')

sample_orders = nan(N_sample,3,length(eig_idxs));
ranks = nan(length(eig_idxs),N_sample);
rconds = nan(length(eig_idxs),N_sample);

for eig_idx_idx = eig_idxs
    eig_idx = eig_idxs(eig_idx_idx);
    
    wb = waitbar(0,['Finding sampling order for eig\_idx = ' num2str(eig_idx)]);
    
    isBreak = false;
    
    covariance_options.eig_idxs = eig_idx;
    [Cs,C_grads,kfcns,kfcn_grads,Vs,V_grads,vfcns,vfcn_grads,X_grid_vec,Y_grid_vec,Z_grid_vec,covariance_info] = get_empirical_covariance(WAVEVECTOR_DATA,EIGENVALUE_DATA,covariance_options);
    kfcn = kfcns{1};
    vfcn = vfcns{1};
    
    [V,D] = eig(reshape(Cs{1},[prod(N_wv_train) prod(N_wv_train)]));
    covariance_eigenvalues(eig_idx_idx,:) = sort(diag(D),'descend');
    
    nexttile
    title(['eig\_idx = ' num2str(eig_idx)])
    wv_s = [];
    sample_idx_labels = {};
    %     hold on
    
    % V = permute(vfcn(wv_e,'gridded'),[2 1 3]);
    V = vfcn(wv_e,'scattered'); % V does not have sensible array ordering (it's slices are not smooth) % THIS MAY NO LONGER BE TRUE
    variances = reshape(V,[],1); % variances does, however, line up with the values of wv_e in a way that is smooth and correct % THIS MAY NO LONGER BE TRUE

    cla
    hold on
    % imagesc(wv_e(:,1),wv_e(:,2),V)
    scatter3(wv_e(:,1),wv_e(:,2),wv_e(:,3),[],variances,'filled') % this plot shows the smoothness and correctness of the ordering of variances with the ordering of wv_e
    % set(gca,'YDir','normal')
    colorbar
    daspect([1 1 1])
    view(3)
    axis tight
    
    if isPause
        pause
    end
    
    for sample_idx = 1:N_sample
        if ceil(sample_idx/10) == sample_idx/10
            waitbar(sample_idx/N_sample,wb,['Finding sampling order for eig\_idx = ' num2str(eig_idx)]);
        end
        if isBreak
            break
        end
        
        [~,idx] = max(variances);
        wv_s = [wv_s; wv_e(idx,:)];
        
        sample_idx_labels{end + 1} = num2str(sample_idx);
        
        if isPlot
            scatter3(wv_s(:,1),wv_s(:,2),wv_s(:,3),100,'r')
            text(wv_s(:,1),wv_s(:,2),wv_s(:,3), sample_idx_labels,'HorizontalAlignment','center','VerticalAlignment','middle')
        end
        
        if isPause
            pause
        end

        V_temp = vfcn(wv_e,'scattered');
        t1_comp = reshape(V_temp,[],1);
        t2_comp = kfcn(wv_e,wv_s,'scattered');
        t3_comp = kfcn(wv_s,wv_s,'scattered');
        ranks(eig_idx_idx,sample_idx) = rank(t3_comp);
        rconds(eig_idx_idx,sample_idx) = rcond(t3_comp);
        if rconds(eig_idx_idx,sample_idx) < rcond_tol
            isBreak = true;
            break
        end
        t4_comp = permute(t2_comp,[2 1]);
        % t2_comp = permute(t2_comp,[3 2 1]);
        % t5_comp = t3_comp\t4_comp;
        % t5_comp = permute(t5_comp,[1 3 2]);
        % t_subtract_comp = permute(pagemtimes(t2_comp,t5_comp),[3 2 1]);
        % variances_comp = t1_comp - t_subtract_comp;
        % variances = variances_comp;
        % V = reshape(variances,N_evaluate(2),N_evaluate(1),N_evaluate(3));
        
        % ABOVE: This is the method I used to use to only compute the
        % matrix multiplication for the entries along the diagonal. Note
        % that the method BELOW works, but it computes an entire matrix,
        % only to filter out the diagonal elements and throw away the rest
        % of the matrix.

        t5_comp = t3_comp\t4_comp;
        t_subtract_comp = t2_comp*t5_comp;
        % t_subtract_comp = t1_comp*t5_comp;
        variances_comp = t1_comp - diag(t_subtract_comp);
        variances = variances_comp;
        
        if isPlot
            cla
            hold on
            title(['eig\_idx = ' num2str(eig_idx)])
            % imagesc(wv_e(:,1),wv_e(:,2),V)
            % set(gca,'YDir','normal')
            scatter3(wv_e(:,1),wv_e(:,2),wv_e(:,3),[],variances,'filled')
            scatter3(wv_s(:,1),wv_s(:,2),wv_s(:,3),100,'r')
            text(wv_s(:,1),wv_s(:,2),wv_s(:,3),sample_idx_labels,'HorizontalAlignment','center','VerticalAlignment','middle')
            colorbar
            daspect([1 1 1])
            axis tight
        end
        
        if isPause
            pause
        end
        
    end
    close(wb)
    sample_orders(1:size(wv_s,1),:,eig_idx_idx) = wv_s;
    final_Vs{eig_idx_idx} = V;
end

figure
imagesc(isnan(rconds))
title('isnan(rconds)')
xlabel('sample\_idx')
ylabel('eig\_idx\_idx')

figure
plot(1:size(covariance_eigenvalues,2),covariance_eigenvalues,'.')
set(gca,'yscale','log')
ylabel('eigenvalue of covariance matrix')
xlabel('index of eigenvalue')
for eig_idx_temp = 1:size(covariance_eigenvalues,1)
    legend_labels{eig_idx_temp} = ['eig\_idx = ' num2str(eig_idx_temp)];
end
legend(legend_labels,'location','northeastoutside')

vars_to_save = {'sample_orders','eig_idxs','sigma','rconds','ranks','final_Vs','data_path_train','N_evaluate','disp_idxs'};
save('sample_order_data',vars_to_save{:})



