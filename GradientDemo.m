N=20;
V=300;
K=5;
rho0=0.05;
% rho0 = 0;
rho=0.5;
% rhoMin=0.5; % minimum non-zero fraction of rho in a block
rhoMin=1;
ngrads=2;
ThresFrac=0.1; % threshold fraction for the gradient matrix

Lab=repmat(1:K,[V/K,1]);Lab=Lab(:);

xx=linspace(1,rhoMin,V/K);
Sig0=toeplitz(xx);
Sig=rho0+(1-rho)*eye(V)+rho*kron(eye(K),Sig0);

Y = mvnrnd(repmat(0,[1,V]),Sig,N);

addnoise = 1;
if addnoise == 1
    Y = Y + randn(size(Y))*0.75;
end

% Y = randn(N,V)

S = corrcoef(Y);
[V,D]=eig(S);V=fliplr(V);D=flip(diag(D));

% subplot(2,2,1);imagesc(Sig,[-0.1,min(1,rho*2)]);axis image;colorbar
subplot(2,3,1);imagesc(S,[-0.1,min(1,rho*2)]);axis image;colorbar
grmat_top = keep_top(S, ThresFrac);
% grmat_top(grmat_top > 0) = grmat_top(grmat_top > 0) - (3/2)*mean(grmat_top(grmat_top > 0));
subplot(2,3,2);imagesc(grmat_top,[-0.1,min(1,rho*2)]);axis image;colorbar
grmat_grad = normalized_angle(grmat_top);
subplot(2,3,3);imagesc(grmat_grad,[-0.1,min(1,rho*2)]);axis image;colorbar

% subplot(2,2,3);plot(V(:,1:3));set(get(gca,'Children'),'LineWidth',2)
% legend({'PC1','PC2','PC3'})
[V,D]=eig(S); V=fliplr(V);
subplot(2,3,4);
[coeff, score, ~] = pca(S);
ngrads = 2;
grads = score(:, 1:ngrads);

% plot(grads(:,1), grads(:,2), '*')
gscatter(grads(:,1), grads(:,2), Lab)
legend off;
title('PCA')
% plot(V(:,1),V(:,2), '*')
% [coeff, score, ~] = pca(S);
% grads = score(:, 1:ngrads);
% 
% plot(grads(:,1), grads(:,2), '*')

subplot(2,3,5)

% grmat_grad = S;
% [V, D] = eig(grmat_grad); V=fliplr(V);
[coeff, score, ~] = pca(grmat_top);
grads = score(:, 1:ngrads);

% plot(grads(:,1), grads(:,2), '*')
gscatter(grads(:,1), grads(:,2), Lab)
legend off;
title('Just thresholding')

subplot(2,3,6)

% grmat_grad = S;
% [V, D] = eig(grmat_grad); V=fliplr(V);
[coeff, score, ~] = pca(grmat_grad);
grads = score(:, 1:ngrads);

% plot(grads(:,1), grads(:,2), '*')
gscatter(grads(:,1), grads(:,2), Lab)
legend off;
title('Gradients')
% subplot(2,2,4); plot(V(:,1),V(:,2), '*')
% xlim

% plot(grmat_grad*V(:,1), grmat_grad*V(:,2), '*')
BigFont(25)

%% On the original matrix 
[coeff, score, ~] = pca(Sig);
grads = score(:, 1:ngrads);
plot(grads(:,1), grads(:,2), '*')

%%
N=20;
V=500;
K=10;
rho0=0.05;
% rho0 = 0;
rho=0.5;
% rhoMin=0.5; % minimum non-zero fraction of rho in a block
rhoMin=1;

Lab=repmat(1:K,[V/K,1]);Lab=Lab(:);

xx=linspace(1,rhoMin,V/K);
Sig0=toeplitz(xx);
Sig=rho0+(1-rho)*eye(V)+rho*kron(eye(K),Sig0);

Y = mvnrnd(repmat(0,[1,V]),Sig,N);

addnoise = 1;
if addnoise == 1
    Y = Y + randn(size(Y))*0.75;
end

% Y = randn(N,V)

S = corrcoef(Y);
[V,D]=eig(S);V=fliplr(V);D=flip(diag(D));

% subplot(2,2,1);imagesc(Sig,[-0.1,min(1,rho*2)]);axis image;colorbar
subplot(2,3,1);imagesc(S,[-0.1,min(1,rho*2)]);axis image;colorbar
grmat_top = keep_top(S, 0.1);
% grmat_top(grmat_top > 0) = grmat_top(grmat_top > 0) - (3/2)*mean(grmat_top(grmat_top > 0));
subplot(2,3,2);imagesc(grmat_top,[-0.1,min(1,rho*2)]);axis image;colorbar
grmat_grad = normalized_angle(grmat_top);
subplot(2,3,3);imagesc(grmat_grad,[-0.1,min(1,rho*2)]);axis image;colorbar

% subplot(2,2,3);plot(V(:,1:3));set(get(gca,'Children'),'LineWidth',2)
% legend({'PC1','PC2','PC3'})
[V,D]=eig(S); V=fliplr(V);
subplot(2,3,4);
[coeff, score, ~] = pca(S);
ngrads = 3;
grads = score(:, 1:ngrads);

% plot(grads(:,1), grads(:,2), '*')
scatter3(grads(:,1), grads(:,2), grads(:,3), 36, Lab, 'filled'); 
legend off;
title('PCA')
% plot(V(:,1),V(:,2), '*')
% [coeff, score, ~] = pca(S);
% grads = score(:, 1:ngrads);
% 
% plot(grads(:,1), grads(:,2), '*')

subplot(2,3,5)

% grmat_grad = S;
% [V, D] = eig(grmat_grad); V=fliplr(V);
[coeff, score, ~] = pca(grmat_top);
grads = score(:, 1:ngrads);

% plot(grads(:,1), grads(:,2), '*')
scatter3(grads(:,1), grads(:,2), grads(:,3), 36, Lab, 'filled'); 
legend off;
title('Just thresholding')

subplot(2,3,6)

% grmat_grad = S;
% [V, D] = eig(grmat_grad); V=fliplr(V);
[coeff, score, ~] = pca(grmat_grad);
grads = score(:, 1:ngrads);

% plot(grads(:,1), grads(:,2), '*')
scatter3(grads(:,1), grads(:,2), grads(:,3), 36, Lab, 'filled'); 
legend off;
title('Gradients')
% subplot(2,2,4); plot(V(:,1),V(:,2), '*')
% xlim
fullscreen
% plot(grmat_grad*V(:,1), grmat_grad*V(:,2), '*')
BigFont(25)
