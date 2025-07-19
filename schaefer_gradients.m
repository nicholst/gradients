load('schaefer.mat')

%%
subplot(1,2,1)
imagesc(grmat)
subplot(1,2,2)
doac = 1;
% imagesc(apower(grmat, 2.3))
grmat_top = keep_top(grmat, 0.1);
grmat_grad = normalized_angle(grmat_top, doac);

grmat_grad_sim = normalized_angle(grmat, doac);
imagesc(grmat_grad)

grmat_grad_noac = normalized_angle(grmat_top, 0);

%%
[V, D] = eig(grmat); V=fliplr(V);
plot(V(:,1), V(:,2), '*')

%% Both hemispheres
subplot(1,3,1)
ngrads = 2;
% [coeff, score, ~] = pca(apower(grmat, 2.3));
[coeff, score, ~] = pca(grmat);

% Extract the top `ngrads` components
grads = score(:, 1:ngrads);

plot(grads(:,1), grads(:,2), '*')
title('PCA')

subplot(1,3,2)
ngrads = 2;
% [coeff, score, ~] = pca(apower(grmat, 2.3));
% [coeff, score, ~] = pca(grmat_grad);
[coeff, score, ~] = pca(grmat_top);

% Extract the top `ngrads` components
grads = score(:, 1:ngrads);

plot(grads(:,1), grads(:,2), '*')
title('Gradients - just thresholding')

subplot(1,3,3)
[coeff, score, ~] = pca(grmat_grad);
% [coeff, score, ~] = pca(grmat_grad_sim);

% Extract the top `ngrads` components
grads = score(:, 1:ngrads);

plot(grads(:,1), grads(:,2), '*')
title('Gradients')
BigFont(30)
fullscreen

%% Both hemispheres
subplot(1,3,1)
ngrads = 2;
% [coeff, score, ~] = pca(apower(grmat, 2.3));
[coeff, score, ~] = pca(grmat);

% Extract the top `ngrads` components
grads = score(:, 1:ngrads);

plot(grads(:,1), grads(:,2), '*')
title('PCA')

subplot(1,3,2)
ngrads = 2;
% [coeff, score, ~] = pca(apower(grmat, 2.3));
% [coeff, score, ~] = pca(grmat_grad);
[coeff, score, ~] = pca(grmat_grad_noac);

% Extract the top `ngrads` components
grads = score(:, 1:ngrads);

plot(grads(:,1), grads(:,2), '*')
title('Gradients - no ac')

subplot(1,3,3)
[coeff, score, ~] = pca(grmat_grad);
% [coeff, score, ~] = pca(grmat_grad_sim);

% Extract the top `ngrads` components
grads = score(:, 1:ngrads);

plot(grads(:,1), grads(:,2), '*')
title('Gradients')
BigFont(30)
fullscreen

%% Left hemisphere
grmat_lh = grmat(1:200, 1:200);
subplot(1,2,1)
imagesc(grmat_lh)
subplot(1,2,2)
% imagesc(apower(grmat, 2.3))
grmat_top = keep_top(grmat_lh, 0.1);
grmat_grad_lh = normalized_angle(grmat_top);
imagesc(grmat_grad_lh)

%%
subplot(1,2,1)
ngrads = 2;
% [coeff, score, ~] = pca(apower(grmat, 2.3));
[coeff, score, ~] = pca(grmat_lh);

% Extract the top `ngrads` components
grads = score(:, 1:ngrads);

plot(grads(:,1), grads(:,2), '*')
title('PCA')

subplot(1,2,2)
ngrads = 2;
% [coeff, score, ~] = pca(apower(grmat, 2.3));
[coeff, score, ~] = pca(grmat_grad_lh);

% Extract the top `ngrads` components
grads = score(:, 1:ngrads);

plot(grads(:,1), grads(:,2), '*')
title('Gradients')