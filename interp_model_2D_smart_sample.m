function out = interp_model_2D_smart_sample(fr,wv,model_options)
    % Applies and analyzed error of a 2D interpolation model (including
    % GPR)
    
    N_sample = model_options.N_sample;
       
    a = 1; % Hard coded!
    
    N_evaluate = [numel(unique(wv(:,1))) numel(unique(wv(:,2))) numel(unique(wv(:,3)))];
    
    X = reshape(wv(:,1),N_evaluate(1),N_evaluate(2),N_evaluate(3)); % Does the transpose fix the fact that I do the N_wv components out of order?
    Y = reshape(wv(:,2),N_evaluate(1),N_evaluate(2),N_evaluate(3));
    Z = reshape(wv(:,3),N_evaluate(1),N_evaluate(2),N_evaluate(3));
    V = reshape(fr,N_evaluate(1),N_evaluate(2),N_evaluate(3));
    
    original_domain_X = reshape(X(:,1,1),[],1);
    original_domain_Y = reshape(Y(1,:,1),[],1);
    original_domain_Z = reshape(Z(1,1,:),[],1);

    out.N_points = prod(N_sample);

    % if strcmp(model_options.model_name,'GPR')
    %     isTrimRightBoundary = true;
    % else
    %     isTrimRightBoundary = false;
    % end
    
%     [X_s,Y_s] = get_wavevectors(N_sample,a,struct('isTrimRightBoundary',isTrimRightBoundary,'format','grid'));
    [X_e,Y_e,Z_e] = get_wavevectors(N_evaluate,a,struct('format','grid'));
        
%     wv_s = get_wavevectors(N_sample,a,struct('isTrimRightBoundary',isTrimRightBoundary,'format','list'));
    wv_s = model_options.sample_points(1:model_options.N_sample,:);
    wv_e = get_wavevectors(N_evaluate,a,struct('format','list'));
    
    h_x = X_e(2,1,1) - X_e(1,1,1); h_y = Y_e(1,2,1) - Y_e(1,1,1); h_z = Z_e(1,1,2) - Z_e(1,1,1);
    
    V_s = interpn(X,Y,Z,V,wv_s(:,1),wv_s(:,2),wv_s(:,3));
    V_e = interpn(X,Y,Z,V,X_e,Y_e,Z_e);
    if V_e ~= V
        warning('Evaluation points are being interpolated') % This is actually undesirable, and doesn't really need to be done if I have high resolution 'ground truth' datasets. Maybe interpolation of ground truth should never be done since it requires an interpolation model itself and could be inaccurate.
    end
    
    fr_s = reshape(V_s,[],1);
    fr_e = reshape(V_e,[],1);
      
    x_train = wv_s;
    y_train = fr_s;
    
    if strcmp(model_options.model_name,'GPR')
        model = create_GPR_model3(x_train,y_train,model_options.sigma,model_options.kfcn,model_options.kfcn_grad,'scattered');
        fr_pred = model.pred(wv_e,'scattered')';
        V_pred = reshape(fr_pred,[N_evaluate(1) N_evaluate(2) N_evaluate(3)]);
    else
        warning('using unchecked part of code')
        % I'm not sure linear interpolation model will work here
        [X_s,Y_s] = get_wavevectors(model_options.N_wv,a,struct('isTrimRightBoundary',isTrimRightBoundary,'format','grid'));
        
        V_pred = interp2(X_s,Y_s,reshape(Z_s,size(X_s)),X_e,Y_e,model_options.model_name);
        model = {};
    end
    
    V_err = V_pred - V_e; % in matrix format
    [grad_x,grad_y,grad_z] = gradient(V_err,h_x,h_y,h_z);
    derrdgamma = cat(4,grad_x,grad_y,grad_z);
    
    e_L2 = LP_norm(X_e,Y_e,Z_e,V_err,2);
    e_H1 = H1_norm(X_e,Y_e,Z_e,V_err);

    % Z_err_rel = zeros(size(Z_e));
    % idxs = Z_e~=0;
    % Z_err_rel(idxs) = Z_err(idxs)./Z_e(idxs);
    % e_L2_rel = LP_norm(X_e,Y_e,Z_err_rel,2);
    % e_H1_rel = H1_norm(X_e,Y_e,Z_err_rel);
    
    out.e_L2 = e_L2;
    out.e_H1 = e_H1;
    % out.e_L2_rel = e_L2_rel;
    % out.e_H1_rel = e_H1_rel;
    out.V_err = V_err;
    out.V_pred = V_pred;
    out.X_e = X_e;
    out.Y_e = Y_e;
    out.Z_e = Z_e;
    out.V_e = V_e;
%     out.X_s = X_s;
%     out.Y_s = Y_s;
    out.V_s = V_s;
    out.wv_s = wv_s;
    out.original_domain_X = original_domain_X;
    out.original_domain_Y = original_domain_Y;
    out.original_domain_Z = original_domain_Z;
%     out.covariance = covariance;
    out.model = model;
end