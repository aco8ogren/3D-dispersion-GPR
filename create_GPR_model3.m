function model = create_GPR_model3(x_train,y_train,sigma,kfcn,kfcn_grad,train_format)
    % This function creates a structure. To predict on new points, use
    % y_new = model.pred(x_new)
    % To compute the gradient of the prediction at x_new with respect to
    % the training points, use 
    % d(y_new)/d(y_train) = model.grad(x_new)
    model.kfcn = kfcn;
    model.kfcn_grad = kfcn_grad;
    model.alpha = get_alpha(x_train,y_train,kfcn,sigma,train_format);
    model.pred = @(x_pred,query_format) get_pred(x_pred,x_train,model.alpha,kfcn,query_format);

    % TODO someday or never
    % model.grad_dy = @(x_pred,query_format) get_grad_dy(x_pred,x_train,sigma,kfcn,query_format);
    % model.grad_dx = @(x_pred,query_format) get_grad_dx(x_pred,x_train,model.alpha,kfcn_grad,query_format);
    % model.grad_dxdy = @(x_pred,query_format) get_grad_dxdy(x_pred,x_train,sigma,kfcn,kfcn_grad,query_format);
    
    % This is a fun idea, maybe for later. Could make sample_order script
    % more elegant
%     model.posterior_covariance = @(x_pred,query_format) kfcn(x_pred,x_pred,query_format) - K(x_pred,x_train,query_format)*K(x_train,x_train,train_format)\K(x_train,x_pred,query_format);
%     model.posterior_variance = @(x_pred,query_format)
end