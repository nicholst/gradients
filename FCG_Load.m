function [G,varargout] = FCG_Load(Dat,varargin)
% [G,ID,AxL]=FCGload(Dat[,Subjs,Gs])
%
%  Dat is 'HCP' or 'SCZ'
%  Subjs is vector of subject indicies to load (or empty for all)
%  Gs is is vector of gradient indicies to load (or empty for all)
%
%  G   - nV x nG x nSubj matrix of gradients
%  ID  - nSubj cell vector of strings, subject identifiers
%  AxL - 2-element list with axis limits
%    

% Parse arguments
 
if nargin<1 || isempty(varargin{1})
    Subjs=[];
else
    Subjs=varargin{1};
end
if nargin<2 || isempty(varargin{2})
    Gs=[];
else
    Gs=varargin{2};
end

% Read data
switch Dat

  case 'HCP'
    Base='/Users/nichols/Proj/SpatMod/2025_Gradients/Data/hcpGrads';
    Subj=dir(Base);
    Subj=Subj(3:end);
    if ~isempty(Subjs)
        Subj=Subj(Subjs);
    end
    nSubj=length(Subj);
    Type='.mapalign.diffmaps.04mm.npy';
    % Type='.mapalign.ses1.diffmap.s04mm.npy']);
    % Type='.pca.ses1.s04mm.npy'
    % Type='.pca.ses2.s04mm.npy'
    
    ID=Subj(1).name;
    a=readNPY(fullfile(Base,ID,[ID Type]));
    if isempty(Gs)
        Gs=1:size(a,1);
    end
    nG=length(Gs);
    nV=size(a,2);
    G=zeros(nV,nG,nSubj);
    for i=1:nSubj
        ID=Subj(i).name;
        tmp=readNPY(fullfile(Base,ID,[ID Type]))';
        G(:,:,i)=tmp(:,Gs);
        if i>1
            G(:,:,i)=G(:,:,i)*diag(sign(sum(G(:,:,i).*G(:,:,1))));
        end
    end

    ID = {Subj.name}';

    AxL={[-8,8],[-8,8]};  % axis limits

  case 'SCZ'

    Base='/Users/nichols/Proj/SpatMod/2025_Gradients/Data/';
    G=h5read(fullfile(Base,'scz_gradients.nc'),'/grads');
    if isempty(Subjs)
        Subjs=1:size(G,3);
    else
        G=G(:,:,Subjs);
    end
    nSubj=size(G,3);
    if ~isempty(Gs)
        G=G(:,Gs,:);
    end
    nG=size(G,2);

    ID=h5read(fullfile(Base,'scz_gradients.nc'),'/subject');
    ID=cellstr(ID(Subjs));   

    AxL={[-4,3],[-2.5,2.5]};  % axis limits
end

if nargout>=1
    varargout{1}=ID;
end
if nargout>=2
    varargout{2}=AxL;
end

end

