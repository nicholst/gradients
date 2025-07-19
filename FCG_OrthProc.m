function [R,u,v] = FCGOrthProc(G,Gt,Pow)
% Orthogonal Procrustes registration, one subject to a template
%
% G  - nV x nG matrix, to be moved
% Gt - nV x nG matrix, template/reference
% Pow - power to set weights; 0 no weighting; 2 square distance from origin    

    % Register G to Gt (template)
    % Solve rotation-only Procrustes: R = argmin || G R - Gt ||_F
    if Pow==0
        [u,s,v]=svd(G'*Gt);
    else
        W=(sqrt(sum(G.^2,2)).^Pow).*(sqrt(sum(Gt.^2,2)).^Pow);
        [u,s,v]=svd(G'*(W.*Gt));
    end
    R = u*v';

end
