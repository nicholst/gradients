function Gr = FCG_ApplyReg(G,R,varargin)
% FORMAT Gr = FCG_ApplyReg(G,R[,S])
% Apply registration
%
% G   - nV x nG x nSubj matrix
% R   - nSubjx1 cell array of nG x nG transformation matricies, such that 
%       G(:,:,i)*R{i} is aligned to atlas
% S   - nSubj cell array of nG x nG diagonal scaling matricies, such that
%       G(:,:,i)*R{i}.*S{i} is aligned to atlas
%
    

    nSubj=size(G,3);

    if isempty(varargin)
        S=cell(1,nSubj);S(:)={1};
    else
        S=varargin{1};
    end

    Gr=zeros(size(G));
    for i=1:nSubj
        Gr(:,:,i)=G(:,:,i)*R{i}.*S{i};
    end

end
