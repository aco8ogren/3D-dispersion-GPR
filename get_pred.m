function pred = get_pred(x_pred,x_train,alpha,kfcn,query_format)
   pred = kfcn(x_pred,x_train,query_format)*alpha;
end