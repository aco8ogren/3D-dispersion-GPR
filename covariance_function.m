function C_interp = covariance_function(wv_i,wv_j,query_format,C_griddedInterpolant,covariance_options)
%     X_grid_vec = original_C_struct.X_grid_vec;
%     Y_grid_vec = original_C_struct.Y_grid_vec;
%     C = original_C_struct.C;
%     gridded_interpolant_function = original_C_struct.gridded_interpolant_function;
%     C_gpu = original_C_struct.C_gpu;
    % query_format can be 'scattered' or 'gridded'
%     M = length(X_grid_vec);
%     N = length(Y_grid_vec);
%     if numel(C)~=M*N*M*N
%         error('error in covariance_function')
%     end
%     C4D = reshape(C,M,N,M,N);
    
    if strcmp(query_format,'scattered')
%         wv = combvec(wv_i',wv_j'); % gives a 4 x N_combinations array
        wv = combvec2(wv_i',wv_j'); % gives a 4 x N_combinations array. Faster, but less general than MatLab's combvec.
%         C_interp_4D = interpn(X_grid_vec,Y_grid_vec,X_grid_vec,Y_grid_vec,C4D,wv(1,:),wv(2,:),wv(3,:),wv(4,:)); % important that query vecs have same orientation on this line --> signals to interpn that they represent a list of scattered points in R^n.
        C_interp_4D = C_griddedInterpolant(wv');
        %         C_interp_4D = myGriddedInterpolant(wv);
        [N_h,~] = size(wv_i);
        [M_h,~] = size(wv_j);
        C_interp = reshape(C_interp_4D,N_h,M_h);
    elseif strcmp(query_format,'scattered - precomputed combvec')
        wv = wv_i.values;
        N_h = wv_i.size(1);
        M_h = wv_i.size(2);
        C_interp_4D = C_griddedInterpolant(wv');
        C_interp = reshape(C_interp_4D,N_h,M_h);
    elseif strcmp(query_format,'gridded')
        % Requires wavevectors to be input sorted by wv_x (for both wv_i and wv_j)
        wv_i_x = sort(unique(wv_i(:,1)));
        wv_i_y = sort(unique(wv_i(:,2)));
        wv_j_x = sort(unique(wv_j(:,1)));
        wv_j_y = sort(unique(wv_j(:,2)));
        if covariance_options.isAllowGPU && length(wv_i_x)*length(wv_i_y)*length(wv_j_x)*length(wv_j_y) > 40^4
            %             C4D_gpu = gpuArray(C4D);
            %             C4D_gpu = reshape(original_C_struct.C_gpu,M,N,M,N);
            C4D_gpu = original_C_struct.C_gpu;
            C_interp_4D = interpn(X_grid_vec,Y_grid_vec,X_grid_vec,Y_grid_vec,C4D_gpu,wv_i_x,wv_i_y',wv_j_x,wv_j_y'); % important that query vecs have different orientation on this line --> signals to interpn that query vecs are grid vecs in R^n.
%             C_interp_4D = gather(C_interp_4D);
        else
%             C_interp_4D = interpn(X_grid_vec,Y_grid_vec,X_grid_vec,Y_grid_vec,C4D,wv_i_x,wv_i_y',wv_j_x,wv_j_y'); % important that query vecs have different orientation on this line --> signals to interpn that query vecs are grid vecs in R^n.
            C_interp_4D = C_griddedInterpolant({wv_i_x,wv_i_y,wv_j_x,wv_j_y});
        end
        C_interp_4D = permute(C_interp_4D,[2 1 4 3]);
        C_interp = reshape(C_interp_4D,size(wv_i_x,1)*size(wv_i_y,1),size(wv_j_x,1)*size(wv_j_y,1));
    else
        error('query_format not recognized')
    end
end