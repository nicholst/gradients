function [Gr,Gm]=FCG_Reg1(G,Gt)
% Register to 1 template
% FUNCTION [Gr,Gm]=FCG_Reg1(G,Gt)
%
%  G   - nV x nG x nSubj matrix of gradients
%  Gt  - nV x nG matrix, template/reference -- or if scalar, then G(:,:,Gt) is template

if isscalar(Gt)
    iRef=Gt;
    Gt=G(:,:,iRef);
else
    iRef=[];
end


% Standard method...
    
%  Register all to 1st
R   = FCG_OrthProcReg(G,Gt);
Gr  = FCG_ApplyReg(G,R);
% Compute mean, the target
Gt  = mean(Gr,3);
% Register to mean
R   = FCG_OrthProcReg(Gr,Gt);
Gr  = FCG_ApplyReg(Gr,R);
Gm  = mean(Gr,3);

end
