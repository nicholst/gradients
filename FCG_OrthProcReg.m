function Rs = FCG_OrthProcReg(G,Gt)
% FUNCTION Rs = FCG_OrthProcReg(G,Gt)
% Orthogonal Procrustes group registration, all subjects to a template
%   
% G  - nV x nG x nSubj matrix
% Gt - nV x nG matrix, template/reference -- or if scalar, then G(:,:,Gt) is template

    if isscalar(Gt)
        iRef=Gt;
        Gt=G(:,:,iRef);
    else
        iRef=[];
    end

    nV = size(G,1);
    nG = size(G,2);
    nSubj = size(G,3);

    Rs = cell(nSubj,1);

    for i=1:nSubj
        if (i==iRef)
            Rs{i}=eye(nG);
        else
            Rs{i}=FCG_OrthProc(G(:,:,i),Gt,0);
        end
    end
end
