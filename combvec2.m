function y = combvec2(wv_i,wv_j)
    %COMBVEC Create all combinations of vectors.
    
    M = size(wv_i,2);
    N = size(wv_j,2);
    
    y = [repmat(wv_i,1,N); repelem(wv_j,1,M)];
end
