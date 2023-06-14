function out = H1_norm(X,Y,Z,V)
    % Only works for scalar valued functions V defined on uniform grids
    % X,Y,Z
    % with spacings h_x,h_y,h_z
    h_x = X(2,1,1) - X(1,1,1); h_y = Y(1,2,1) - Y(1,1,1); h_z = Z(1,1,2) - Z(1,1,1);
    [V_x,V_y,V_z] = gradient(V,h_x,h_y,h_z);
    
    out = sqrt(LP_norm(X,Y,Z,V_x,2)^2 + LP_norm(X,Y,Z,V_y,2)^2 + LP_norm(X,Y,Z,V_z,2)^2);
end