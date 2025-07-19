function [Gr,iRef,R,S]=FCG_RegAll(G,Pow,Scale,varargin);
% Compute all possible registrations, select closest, register to closest
% FUNCTION [Gr,iRef,R,S]=FCG_RegAll(G,Pow,Scale[,Gt])
%
%  G   - nV x nG x nSubj matrix of gradients
%  Pow - power to set weights; 0 no weighting; 2 square distance from origin    
%  Scale - 0 no scaling; >0 scaling (nG parmeters), <0 affine (15 parameters)
%  Gt   - empty for all possible; otw nV x nG matrix or scalar (G(:,:,Gt)) to set template
%    
%
%  Gr   -  nV x nG x nSubj matrix of gradients, registered
%  iRef -  Subject index of target; 
%  R    -  nSubj x nSubj cell array of rotation matricies, R{i,j} maps subject
%          i to subject j     
%  S    -  nSubj x nSubj cell array of scaling vectors, S{i,j} scales subject
%          i to subject j
%
%  Outputs satisfy
%    Gr(:,:,i) = G(:,:,i) * R{i,iRef} .* S{i,iRef}
%    

nScaleIter=3;

iRef=[];
if nargin<4 || isempty(varargin{1})
    Gi=[];
else
    if isscalar(varargin{1})
        iRef=varargin{1};
        Gi=G(:,:,iRef);
    else
        iRef=NaN;
        Gi=varargin{1};
    end
end

nG=size(G,2);
nSubj=size(G,3);
    
% Initialise output
R=cell(nSubj,nSubj);R(:,:)={eye(nG)}; % R{i,j} rotation of i to match j
S=cell(nSubj,nSubj);S(:,:)={ones(1,nG)};
DE=zeros(nSubj,nSubj); % Eucledian distance (optimised cost function)

% Estimate All Registrations 
for i = 1:nSubj
    if isempty(iRef)
        fprintf('* -> %d\n',i)
        Gi = G(:,:,i);
    elseif ~isnan(iRef)
        fprintf('* -> %d\n',iRef)
    else
        fprintf('* -> Gt\n')
    end
    for j = i+1:nSubj
        % Register subj j to subj i
        Gj = G(:,:,j);
        if Scale < 0 % affine
            % Solve affine Procrustes: R = argmin_R || G_i R - G_j ||_F
            % ... but note registration no longer symmetric (I think!?)
            R{j,i} = pinv(Gj)*Gi;
            R{i,j} = pinv(Gi)*Gj;
        else
            % Solve rotation-only Procrustes: R = argmin || G_i R - G_j ||_F
            R{j,i} = FCG_OrthProc(Gj,Gi,Pow);
            R{i,j} = R{j,i}';
            if Scale>0
                % Solve rotation-scale Procrustes: R = argmin || G_i R S - G_j ||_F
                % ... but note registration no longer symmetric
                Rj2i=R{j,i};
                Ri2j=R{i,j};
                for k=1:nScaleIter
                    [Sj2i,Si2j]=deal(ones(1,nG));
                    if k>1
                        Rj2i=Rj2i*FCG_OrthProc(Gj*Rj2i,Gi./Sj2i,Pow);
                        Ri2j=Ri2j*FCG_OrthProc(Gi*Ri2j,Gj./Si2j,Pow);
                    end
                    Gj2i = Gj*Rj2i;  
                    Gi2j = Gi*Ri2j;
                    Sj2i = Sj2i.*FCG_ScaleReg(Gj2i,Gi./Sj2i);
                    Si2j = Si2j.*FCG_ScaleReg(Gi2j,Gj./Si2j);
                    Gj2i = Gj.*Sj2i; 
                    Gi2j = Gi.*Si2j;
                end
                R{j,i}=Rj2i;
                R{i,j}=Ri2j;
                S{j,i}=Sj2i;
                S{i,j}=Si2j;
            end
        end
        DE(j,i)=NormEucl(Gj*R{j,i}.*S{j,i} - Gi);
        DE(i,j)=NormEucl(Gi*R{i,j}.*S{i,j} - Gj);
    end
    if ~isempty(iRef)
        R=R(:,1);
        S=cellfun(@(a,b)(a+b)/2,S(:,1),S(1,:)','uniformoutput',false);
        DE=(DE(:,1)+DE(1,:)')/2;
        break
    end
end

if isempty(iRef)
    % Ref: Best Euclidean subject
    if Scale==0
        [~,iRef]=min(sum(DE));
    else
        [~,iRef]=min(sum(DE(:,:,end))+sum(DE(:,:,end),2)');
    end
    Gr  = FCG_ApplyReg(G,R(:,iRef),S(:,iRef));
else
    Gr  = FCG_ApplyReg(G,R,S);
end

end

function E = NormEucl(Err)
% Mean Eucledian norm on registration errors
    E = mean(sqrt(sum(Err.^2,2)));
end
