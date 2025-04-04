x = [1, 0;0,1]

y = pdist(x, 'cosine')% This is a list of 1 - inner product between the rows.

z = 1 - squareform(pdist(x, 'cosine')); % Each entry of this matrix is the inner product. 

z2 = 1 - acos(z) / pi % This is then the angle between each row. 

% z2 = 1 - acos(z) / (pi/2)

%%
% inner_product = (x/normalizeX2(x)'*(x/normalizeX2(x);
1 - acos((normalizeX2(x))'*(normalizeX2(x)))/pi % this is one 1 minus the angle between each row.

%%
data = randn()

%%

x = 1 - squareform(pdist(x, 'cosine'));

%%
nsim = 1000;
a = mvnrnd(zeros([1,500]), eye(500), nsim);
size(cov(a));
imagesc(cov(a))

S = cov(a);

S_top = keep_top(S, 0.1);
S_grad = normalized_angle(S);
imagesc(S_grad)

subplot(1,2,1)
[coeff, score, ~] = pca(S);

% Extract the top `ngrads` components
grads = score(:, 1:ngrads);

plot(grads(:,1), grads(:,2), '*')
title('PCA')

subplot(1,2,2)
[coeff, score, ~] = pca(S_grad);

% Extract the top `ngrads` components
grads = score(:, 1:ngrads);

plot(grads(:,1), grads(:,2), '*')

title('gradients')
