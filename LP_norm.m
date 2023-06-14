function out = LP_norm(X,Y,Z,V,p)
    % Only works for scalar valued functions V defined on uniform grids
    % X,Y,Z
    % with spacings h_x,h_y,h_z
    h_x = X(2,1,1) - X(1,1,1); h_y = Y(1,2,1) - Y(1,1,1); h_z = Z(1,1,2) - Z(1,1,1);
    out = (sum(V.^p,'all')*h_x*h_y*h_z)^(1/p);
end