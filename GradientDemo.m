N=50;
V=300;
K=10;
rho0=0.05;
rho=0.5;
% rhoMin=0.5; % minimum non-zero fraction of rho in a block
rhoMin=1;

xx=linspace(1,rhoMin,V/K);
Sig0=toeplitz(xx);
Sig=rho0+(1-rho)*eye(V)+rho*kron(eye(K),Sig0);

Y = mvnrnd(repmat(0,[1,V]),Sig,N);

addnoise = 1;
if addnoise == 1
    Y = Y + randn(size(Y))*0.75;
end

S = corrcoef(Y);
[V,D]=eig(S);V=fliplr(V);D=flip(diag(D));


% subplot(2,2,1);imagesc(Sig,[-0.1,min(1,rho*2)]);axis image;colorbar
subplot(2,2,1);imagesc(S,[-0.1,min(1,rho*2)]);axis image;colorbar
grmat_top = keep_top(S, 0.1);
grmat_grad = normalized_angle(grmat_top);
subplot(2,2,2);imagesc(grmat_grad,[-0.1,min(1,rho*2)]);axis image;colorbar

% subplot(2,2,3);plot(V(:,1:3));set(get(gca,'Children'),'LineWidth',2)
% legend({'PC1','PC2','PC3'})
[V,D]=eig(S); V=fliplr(V);
subplot(2,2,3);
[coeff, score, ~] = pca(S);
grads = score(:, 1:ngrads);

plot(grads(:,1), grads(:,2), '*')
% plot(V(:,1),V(:,2), '*')
% [coeff, score, ~] = pca(S);
% grads = score(:, 1:ngrads);
% 
% plot(grads(:,1), grads(:,2), '*')

subplot(2,2,4)

% grmat_grad = S;
% [V, D] = eig(grmat_grad); V=fliplr(V);
[coeff, score, ~] = pca(grmat_grad);
grads = score(:, 1:ngrads);

plot(grads(:,1), grads(:,2), '*')
% subplot(2,2,4); plot(V(:,1),V(:,2), '*')
% xlim
fullscreen
% plot(grmat_grad*V(:,1), grmat_grad*V(:,2), '*')


