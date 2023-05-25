clear; close all;

% data = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\output 18-May-2023 18-16-56\checkpoint152 - Copy.mat");
data = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\full BZ output 18-May-2023 18-16-56\checkpoint233_ndgrid.mat");
% data = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\output 24-May-2023 18-10-17\checkpoint2 - Copy.mat");

eig_idx = 1;

e = data.EIGENVALUE_DATA(:,:,data.c.struct_idxs);
w = data.WAVEVECTOR_DATA(:,:);

idxs = sort(w(:,2));

C = cov(squeeze(e(:,eig_idx,:))');

figure
for i = 1:10 % numel(C(1,:))
    c = C(:,i);
    scatter3(w(:,1),w(:,2),w(:,3),[],c(:),'filled')
    pause(0.01)
end

C6D = reshape(C,[data.c.N_wv data.c.N_wv]);

% Show that C6D(:,:,:,idx) gives the covariance of every wavevector in the
% IBZ with a single other wavevector
% The reason this is convincing to me is because C(:,idx) is FOR SURE the
% covariance of every wavevector in the IBZ with a single other wavevector.
figure
idx = 12;
temp = C6D(:,:,:,idx);
diff = C(:,idx) - temp(:);
scatter3(w(:,1),w(:,2),w(:,3),[],diff(:),'filled')
colorbar

% Create an interpolater and make sure it "goes through" the original data
w_grid_vecs = {unique(w(:,1)),unique(w(:,2)),unique(w(:,3))};
gi = griddedInterpolant({w_grid_vecs{:},w_grid_vecs{:}},C6D,'cubic');

w_query = combvec2(w',w')';
C1D_interp = gi(w_query);

C6D_interp = reshape(C1D_interp,[data.c.N_wv data.c.N_wv]);
% C6D_interp = permute(C6D_interp,[2 1 3 5 4 6]);

max(abs(C6D_interp - C6D),[],'all')

% Plot as an image in 2D
C2D_interp = reshape(C6D_interp,[size(w,1) size(w,1)]);
figure
tiledlayout('flow')

nexttile
imagesc(C2D_interp)
title('interpolated')
colorbar
daspect([1 1 1])

nexttile
imagesc(C)
title('original')
colorbar
daspect([1 1 1])

nexttile
imagesc(C2D_interp - C)
title('difference')
colorbar
daspect([1 1 1])

% Plot in a line
idx = 1;

figure
tiledlayout('flow')

nexttile
temp1 = C6D_interp(:,:,:,idx);
plot(temp1(:))

nexttile
temp2 = C6D(:,:,:,idx);
plot(temp2(:))

nexttile
plot(temp1(:)-temp2(:))


% Plot the ordering of WAVEVECTOR_DATA
figure
for i = 1:size(w,1)
    x = w(i,1);
    y = w(i,2);
    z = w(i,3);
    text(x,y,z,num2str(i))
    hold on
    scatter3(x,y,z,'k.')
end
view(3)
daspect([1 1 1])
xlabel('x')
ylabel('y')
zlabel('z')







