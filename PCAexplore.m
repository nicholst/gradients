N=500;
K=4;
V=K*30;
rho0=0.01;
rho=0.15;
% rhoMin=0.5; % minimum non-zero fraction of rho in a block
rhoMin=1;
ngrads=2;
ThresFrac=1/K; % threshold fraction for the gradient matrix
Lab=repmat(1:K,[V/K,1]);Lab=Lab(:);
fwhm=0;

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

% Smoothing kernel
if fwhm==0
    Ks = eye(V);
else
    sig=fwhm/sqrt(8*log(2));
    kernlen=floor(fwhm*3/2)*2+1; % must be odd
    xx=(1:kernlen);xx=xx-mean(xx);
    kern=exp(-xx.^2/(2*sig^2));
    kern=kern/sum(kern);
    Ks=convmtx(kern,V);
    Ks=Ks(:,(kernlen-1)/2+1:end-(kernlen-1)/2); % symmetric & square
    es=mvnrnd(repmat(0,[1,V]),Ks*Ks',N);
end

Ys=Y*Ks';
%Ys=Y+es;

S = corrcoef(Y);
Ss = corrcoef(Ys);
% %  Add FastICA_25 - https://research.ics.aalto.fi/ica/fastica/code/dlcode.shtml
% Si = fastica(Y.^6, 'numOfIC', 2)';


[V,D]=eig(S);V=fliplr(V);D=flip(diag(D));

grmat_top = keep_top(S, ThresFrac);       
grmat_grad = normalized_angle(grmat_top);
grmat_tops = keep_top(Ss, ThresFrac);       
grmat_grads = normalized_angle(grmat_tops);


%% Raw

%subplot(2,2,1);imagesc(Sig,[-0.1,min(1,rho*2)]);axis image;colorbar
subplot(3,3,1);imagesc(S,[-0.1,1]);axis image;colorbar

% subplot(2,2,3);plot(V(:,1:3));set(get(gca,'Children'),'LineWidth',2)
% legend({'PC1','PC2','PC3'})

subplot(3,3,4);
[V,D]=eig(S); V=fliplr(V);
[coeff, score, ~] = pca(S); % Sample covariance
% [coeff, score, ~] = pca(Sig); % True covariance
%grads = score(:, 1:ngrads);
%gscatter(grads(:,1), grads(:,2), Lab)
grads = score(:, 1:ngrads+1);
scatter3(grads(:,1), grads(:,2),grads(:,3), 200,Lab,'marker','.')

subplot(3,3,7)
plot(V(:,1:K))

if (0)
%% Top threshold

%subplot(2,2,1);imagesc(Sig,[-0.1,min(1,rho*2)]);axis image;colorbar
subplot(3,3,1);imagesc(grmat_top,[-0.1,1]);axis image;colorbar

% subplot(2,2,3);plot(V(:,1:3));set(get(gca,'Children'),'LineWidth',2)
% legend({'PC1','PC2','PC3'})

subplot(3,3,4);
%[V,D]=eig(S); V=fliplr(V);
%[coeff, score, ~] = pca(S); % Sample covariance
[coeff, score, ~] = pca(grmat_top); % Thresholded covariance
% [coeff, score, ~] = pca(Sig); % True covariance
%grads = score(:, 1:ngrads);
%gscatter(grads(:,1), grads(:,2), Lab)
grad = score(:, 1:ngrads+1);
scatter3(grad(:,1), grad(:,2),grad(:,3), 200,Lab,'marker','.')

subplot(3,3,7)
plot(V(:,1:K))
end


%% Thresholded

subplot(3,3,2);imagesc(grmat_top,[-0.1,1]);axis image;colorbar

subplot(3,3,5);
%[Vs,Ds]=eig(gra); Vs=fliplr(Vs);
%[coeff, score, ~] = pca(S); % Sample covariance
[coeffs, scores, ~] = pca(grmat_top); % Thresholded covariance
% [coeff, score, ~] = pca(Sig); % True covariance
%gradss = scores(:, 1:ngrads);
%gscatter(gradss(:,1), gradss(:,2), Lab)
gradss = scores(:, 1:ngrads+1);
scatter3(gradss(:,1), gradss(:,2), gradss(:,3),200, Lab,'marker','.')

subplot(3,3,8)
plot(Vs(:,1:K))

if (0)

%% Smoothed

subplot(3,3,2);imagesc(Ss,[-0.1,1]);axis image;colorbar

subplot(3,3,5);
[Vs,Ds]=eig(Ss); Vs=fliplr(Vs);
%[coeff, score, ~] = pca(S); % Sample covariance
[coeffs, scores, ~] = pca(grmat_tops); % Thresholded covariance
% [coeff, score, ~] = pca(Sig); % True covariance
%gradss = scores(:, 1:ngrads);
%gscatter(gradss(:,1), gradss(:,2), Lab)
gradss = scores(:, 1:ngrads+1);
scatter3(gradss(:,1), gradss(:,2), gradss(:,3),200, Lab,'marker','.')

subplot(3,3,8)
plot(Vs(:,1:K))

end

%% Gradients
subplot(3,3,3);imagesc(grmat_grads,[-0.1,1]);axis image;colorbar

subplot(3,3,6)
% grmat_grad = S;
% [V, D] = eig(grmat_grad); V=fliplr(V);
[coeffs, scores, ~] = pca(grmat_grads);
gradss = scores(:, 1:ngrads);
% plot(grads(:,1), grads(:,2), '.','MarkerSize',10)
%gscatter(gradss(:,1), gradss(:,2), Lab)
gradss = scores(:, 1:ngrads+1);
scatter3(gradss(:,1), gradss(:,2), gradss(:,3), 200,Lab,'marker','.')
% subplot(2,2,4); plot(V(:,1),V(:,2), '*')
% plot(grmat_grad*V(:,1), grmat_grad*V(:,2), '*')

subplot(3,3,9)
plot(scores(:,1:K))




