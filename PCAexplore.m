N=20;
K=8;
V=K*30;
rho0=0.05;
rho=0.5;
% rhoMin=0.5; % minimum non-zero fraction of rho in a block
rhoMin=1;
ngrads=2;
ThresFrac=1/K; % threshold fraction for the gradient matrix
Lab=repmat(1:K,[V/K,1]);Lab=Lab(:);

xx=linspace(1,rhoMin,V/K);
Sig0=toeplitz(xx);
Sig=rho0+(1-rho)*eye(V)+rho*kron(eye(K),Sig0);
addnoise = 0;
flipblocks=0; % flip blocks of the covariance matrix
if flipblocks==1
    for i=1:K
        I=(i-1)*V/K+1:i*V/K;
        Sig(I,I)=Sig(I,I)*(rand(1)>=0.5);
    end
    if addnoise == 1
        Sig = Sig + 0.1*randn(size(Sig));
    end
else
    Y = mvnrnd(repmat(0,[1,V]),Sig,N);
    if addnoise == 1
        Y = Y + randn(size(Y))*0.75;
    end
end

S = corrcoef(Y);

[V,D]=eig(S);V=fliplr(V);D=flip(diag(D));

%subplot(2,2,1);imagesc(Sig,[-0.1,min(1,rho*2)]);axis image;colorbar
subplot(2,2,1);imagesc(S,[-0.1,min(1,rho*2)]);axis image;colorbar
grmat_top = keep_top(S, ThresFrac); 
grmat_grad = normalized_angle(grmat_top);
subplot(2,2,2);imagesc(grmat_grad,[-0.1,min(1,rho*2)]);axis image;colorbar

% subplot(2,2,3);plot(V(:,1:3));set(get(gca,'Children'),'LineWidth',2)
% legend({'PC1','PC2','PC3'})
[V,D]=eig(S); V=fliplr(V);
subplot(2,2,3);
%[coeff, score, ~] = pca(S); % Sample covariance
[coeff, score, ~] = pca(grmat_top); % Thresholded covariance
% [coeff, score, ~] = pca(Sig); % True covariance
grads = score(:, 1:ngrads);

% plot(grads(:,1), grads(:,2), '.','MarkerSize',10)
gscatter(grads(:,1), grads(:,2), Lab)

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

% plot(grads(:,1), grads(:,2), '.','MarkerSize',10)
gscatter(grads(:,1), grads(:,2), Lab)
% subplot(2,2,4); plot(V(:,1),V(:,2), '*')
% xlim

% plot(grmat_grad*V(:,1), grmat_grad*V(:,2), '*')

mean(S(:)<0)
