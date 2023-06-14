function alpha = get_alpha(x_train,y_train,kfcn,sigma,train_format)
    alpha = (kfcn(x_train,x_train,train_format) + sigma^2*eye(size(x_train,1),size(x_train,1)))\y_train;
%     alpha = pinv(kfcn(x_train,x_train,train_format) + sigma^2*eye(size(x_train,1),size(x_train,1)),1e-8)*y_train;
end