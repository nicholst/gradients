function S = FCG_ScaleReg(G,Gt)
% Estimate scale so that marginal variance of each gradient of G matches Gt
% G  - nV x nG matrix, to be scaled
% Gt - nV x nG matrix, template/reference
% Scale - Whether to scale; 0 - no scale, return identity, ~0, scale

    nG = size(G,2);
    
    % Match marginal variances of each gradient component after rotation
    S=sqrt(sum(Gt.^2)./(sum(G.^2)+eps));

end
