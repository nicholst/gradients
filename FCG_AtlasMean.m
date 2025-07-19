function [mG,varargout] = FCGAtlasMean(G,R,varargin)
% FORMAT [mG[,sdG]] = FCGAtlasMean(G,R[,S])
% G   - nV x nG x nSubj matrix
% R   - nSubjx1 cell array of nG x nG transformation matricies, such that 
%       G(:,:,i)*R{i} is aligned to atlas
% S   - nSubj cell array of nG x nG diagonal scaling matricies, such that
%       G(:,:,i)*R{i}.*S{i} is aligned to atlas
% 
% mG  - Mean of registered gradients
% sdG - nG x nSubj x 2 matrix: sdGr(:,:,1) stdev before scaling, sdGr(:,:,1) stdev after scaling
          
    nV = size(G,1);
    nG = size(G,2);
    nSubj = size(G,3);

    if isempty(varargin)
        S=cell(1,nSubj);S(:)={1};
    else
        S=varargin{1};
    end
    if nargout==1
        sdGr = [];
    else
        sdGr = zeros(nG,nSubj,2);
    end

    % Mean gradient registered to atlas
    mG=zeros(nV,nG);
    for i=1:nSubj
        Gr = G(:,:,i)*R{i};
        Grs = Gr.*S{i};
        if ~isempty(sdGr)
            sdGr(:,i,1)=std(Gr);
            sdGr(:,i,2)=std(Grs);
        end
        mG=mG+Grs;
    end
    mG=mG/nSubj;
    if ~isempty(sdGr)
        varargout{1}=sdGr;
    end
end
