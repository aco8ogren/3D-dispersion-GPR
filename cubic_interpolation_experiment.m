clear; close all;

% This script is evidence that even when using cubic interpolation, the
% interpolation always passes exactly through all of the original points 
% with exactly zero error

[x,y] = ndgrid(linspace(0,10,1000),linspace(0,10,1000));

z = sin(x.*y.^2);
z(40:60,40:60) = -1;
z(300:500,300:500) = 10;

figure
tlo = tiledlayout('flow')
nexttile
surf(x,y,z)

gi = griddedInterpolant(x,y,z,'cubic');

z_interp = gi(x,y);

nexttile
surf(x,y,z_interp)

figure
z_diff = z_interp - z;
imagesc(z_diff)
set(gca,'colorscale','log')

any(z_diff,'all')


