clear; close all;

% This script shows how I can reorder data that I've previously generated
% using meshgrid ordering to make it follow ndgrid ordering.

N = 3;

for i = 1:3
    grid_vecs{i} = linspace(0,1,N);
end

[X,Y,Z] = meshgrid(grid_vecs{:});

meshgrid_data = [X(:) Y(:) Z(:)];

[X,Y,Z] = ndgrid(grid_vecs{:});

ndgrid_data = [X(:) Y(:) Z(:)];

disp(['any(meshgrid_data ~= ndgrid_data,''all'') --> ' num2str(any(meshgrid_data ~= ndgrid_data,'all'))])

% Re-sort the meshgrid_data to be the same ordering as the ndgrid data
resorted_meshgrid_data = sortrows(ndgrid_data,[3 2 1]); % I think the column sorting order should always by the descending vector NDIM:1

disp(['any(resorted_meshgrid_data ~= ndgrid_data,''all'') --> ' num2str(any(resorted_meshgrid_data ~= ndgrid_data,'all'))])