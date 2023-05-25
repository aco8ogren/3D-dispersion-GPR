clear; close all;

N = 3;

xgv = linspace(0,1,N);
ygv = linspace(0,1,N);
zgv = linspace(0,1,N);

[X,Y,Z] = ndgrid(xgv,ygv,zgv);

L = [X(:) Y(:) Z(:)];

L(1:20,:)

[~,idxs] = sort(L(:,1))

[~,idxs] = sort(L(:,2))

[~,idxs] = sort(L(:,3))