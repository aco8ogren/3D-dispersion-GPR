function V_interp = variance_function(wv,query_format,original_V_struct,variance_options)
    X_grid_vec = original_V_struct.X_grid_vec;
    Y_grid_vec = original_V_struct.Y_grid_vec;
    V = original_V_struct.V;
    %     C_gpu = original_C_struct.C_gpu;
    % query_format can be 'scattered' or 'gridded'
    M = length(X_grid_vec);
    N = length(Y_grid_vec);
    if numel(V)~= M*N
        error('error in variance_function')
    end
    %     C4D = reshape(V,M,N,M,N);
    
    if strcmp(query_format,'scattered')
        %         wv = combvec(wv_i',wv_j'); % gives a 4 x N_combinations array
        %         V_interp = interp2(X_grid_vec,Y_grid_vec,V,wv(1,:),wv(2,:)); % important that query vecs have same orientation on this line --> signals to interpn that they represent a list of scattered points in R^n.
        % interpn and interp2 take arguments in different orders. interpn
        % assumes inputs are in axis order, interp2 flips axis order for
        % geometric reason (axis 2 then axis 1, since x then y)
        V_interp = interpn(X_grid_vec,Y_grid_vec,V,wv(1,:),wv(2,:)); % important that query vecs have same orientation on this line --> signals to interpn that they represent a list of scattered points in R^n.
    elseif strcmp(query_format,'gridded')
        % Requires wavevectors to be input with wv_x varying first and wv_y
        % varying second (for both wv_i and wv_j)
        X_grid_vec_query = sort(unique(wv(:,1)));
        Y_grid_vec_query = sort(unique(wv(:,2)));
        %         V_interp = interp2(X_grid_vec,Y_grid_vec,V,X_grid_vec_query,Y_grid_vec_query'); % important that query vecs have different orientation on this line --> signals to interpn that query vecs are grid vecs in R^n.
        % interpn and interp2 take arguments in different orders. interpn
        % assumes inputs are in axis order, interp2 flips axis order for
        % geometric reason (axis 2 then axis 1, since x then y)
        V_interp = interpn(X_grid_vec,Y_grid_vec,V,X_grid_vec_query,Y_grid_vec_query'); % important that query vecs have different orientation on this line --> signals to interpn that query vecs are grid vecs in R^n.
        
        %         V_interp = V_interp';
        %         wv_i_x = sort(unique(wv_i(:,1)));
        %         wv_i_y = sort(unique(wv_i(:,2)));
        %         wv_j_x = sort(unique(wv_j(:,1)));
        %         wv_j_y = sort(unique(wv_j(:,2)));
        %         if variance_options.isAllowGPU && length(wv_i_x)*length(wv_i_y)*length(wv_j_x)*length(wv_j_y) > 40^4
        %             %             C4D_gpu = gpuArray(C4D);
        %             %             C4D_gpu = reshape(original_C_struct.C_gpu,M,N,M,N);
        %             C4D_gpu = original_V_struct.C_gpu;
        %             C_interp_4D = interpn(X_grid_vec,Y_grid_vec,X_grid_vec,Y_grid_vec,C4D_gpu,wv_i_x,wv_i_y',wv_j_x,wv_j_y'); % important that query vecs have different orientation on this line --> signals to interpn that query vecs are grid vecs in R^n.
        % %             C_interp_4D = gather(C_interp_4D);
        %         else
        %             C_interp_4D = interpn(X_grid_vec,Y_grid_vec,X_grid_vec,Y_grid_vec,C4D,wv_i_x,wv_i_y',wv_j_x,wv_j_y'); % important that query vecs have different orientation on this line --> signals to interpn that query vecs are grid vecs in R^n.
        %         end
        %         C_interp_4D = permute(C_interp_4D,[2 1 4 3]);
        %         C_interp = reshape(C_interp_4D,size(wv_i_x,1)*size(wv_i_y,1),size(wv_j_x,1)*size(wv_j_y,1));
    end
end