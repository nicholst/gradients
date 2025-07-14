function FCGreg(Dat,varargin)
% ReadHCP(Dat,[nSubj,nGrad,Pow,Scale])
%
%  Dat is 'HCP' or 'SCZ'
%

% Supress warnings about non-principal logm, due to R \in O(K) instead of SO(K)
warnState = warning('query', 'MATLAB:logm:nonPosRealEig');
warning('off', 'MATLAB:logm:nonPosRealEig');

% Parse arguments
if nargin<1 || isempty(varargin{1})
    nSubj=Inf;
else
    nSubj=varargin{1};
end
if nargin<2 || isempty(varargin{2})
    nG=Inf;
else
    nG=varargin{2};
end
if nargin<3 || isempty(varargin{3})
    Pow=0;
else
    Pow=varargin{3};
end
if nargin<4 || isempty(varargin{4})
    Scale=0;
else
    Scale=varargin{4};
end

% Read data
switch Dat

  case 'HCP'
    Base='/Users/nichols/Proj/SpatMod/2025_Gradients/Data/hcpGrads';
    Subj=dir(Base);
    Subj=Subj(3:end);
    nSubj=min(nSubj,length(Subj));
    Type='.mapalign.diffmaps.04mm.npy';
    % Type='.mapalign.ses1.diffmap.s04mm.npy']);
    % Type='.pca.ses1.s04mm.npy'
    % Type='.pca.ses2.s04mm.npy'
    
    ID=Subj(1).name;
    a=readNPY(fullfile(Base,ID,[ID Type]));
    nG=min(nG,size(a,1));
    G=zeros(size(a,2),nG,nSubj);
    for i=1:nSubj
        ID=Subj(i).name;
        tmp=readNPY(fullfile(Base,ID,[ID Type]))';
        G(:,:,i)=tmp(:,1:nG);
        if i>1
            G(:,:,i)=G(:,:,i)*diag(sign(sum(G(:,:,i).*G(:,:,1))));
        end
    end
    AxL={[-5,7],[-4,4]};  % axis limits

  case 'SCZ'

    Base='/Users/nichols/Proj/SpatMod/2025_Gradients/Data/';
    G=h5read(fullfile(Base,'scz_gradients.nc'),'/grads');
    nSubj=min(nSubj,size(G,3));
    nG=min(nG,size(G,2));
    G=G(:,1:nG,1:nSubj);

    AxL={[-4,3],[-2.5,2.5]};  % axis limits
end

nV=size(G,1);

% Remove crazy subjects
skew=zeros(1,nSubj);
for i=1:nSubj
    skew(i)=max(skewness(G(:,:,i)));
end
I=skew>=4;
if sum(I)>0
    fprintf("Rejecting %d skew>=4 subjects\n",sum(I));
    G=G(:,:,~I); 
    nSubj=size(G,3);
end

Out=sprintf('%s-nS%d-nG%d-P%d-S%d',Dat,nSubj,nG,Pow,Scale);

% Initialise output
[R0,Rm]=deal(cell(nSubj,1));
[R,S]=deal(cell(nSubj,nSubj));
DF=zeros(nSubj,nSubj); % Frobenius distance, orthogonal distance
DL=zeros(nSubj,nSubj); % Lie algebra norm distance (with flip in case det=-1)
DE=zeros(nSubj,nSubj); % Eucledian distance (optimised cost function)

R0{1}=eye(nG);
G1=G(:,1:nG,1);
for i=2:nSubj
    % Register subj i to subj 1
    % Solve rotation-only Procrustes: R = argmin || G_i R - G_1 ||_F
    Gi=G(:,1:nG,i);
    if 1 % Pow==0
        [u,s,v]=svd(Gi'*G1,'econ');
    else
        W=(sqrt(sum(Gi.^2,2)).^Pow).*(sqrt(sum(G1.^2,2)).^Pow);
        [u,s,v]=svd(Gi'*(W.*G1),'econ');
    end
    R0{i} = u*v';
end
mGr0=zeros(nV,nG);
for i=1:nSubj
    mGr0=mGr0+G(:,1:nG,i)*R0{i};
end
mGr0=mGr0/nSubj;

% Estimate Registration

for i=1:nSubj
    % Register subj i to mean
    % Solve rotation-only Procrustes: R = argmin || G_i R - mGm0 ||_F
    Gi=G(:,1:nG,i);
    if 1 % Pow==0
        [u,s,v]=svd(Gi'*mGr0,'econ');
    else
        W=(sqrt(sum(Gi.^2,2)).^Pow).*(sqrt(sum(mGr0.^2,2)).^Pow);
        [u,s,v]=svd(Gi'*(W.*mGr0),'econ');
    end
    Rm{i} = u*v';

    fprintf('* -> %d\n',i)
    R{i,i}=eye(nG);
    S{i,i}=eye(nG);
    for j=i+1:nSubj
        % Register subj j to subj i
        % Solve rotation-only Procrustes: R = argmin || G_i R - G_j ||_F
        Gj=G(:,1:nG,j);
        if Pow==0
            [u,s,v]=svd(Gj'*Gi,'econ');
        else
            W=(sqrt(sum(Gj.^2,2)).^Pow).*(sqrt(sum(Gi.^2,2)).^Pow);
            [u,s,v]=svd(Gj'*(W.*Gi),'econ');
        end
        Rj2i = u*v';
        R{j,i}=Rj2i;
        R{i,j}=Rj2i';
        DF(j,i)=norm(Rj2i - eye(nG), 'fro'); 
        DF(i,j)=DF(j,i);
        if det(Rj2i)>0
            DL(j,i)=norm(logm(Rj2i));
        else
            % Project into SO(n) to avoid logm error on det=-1
            u(:,end)=-u(:,end);
            DL(j,i)=norm(logm(u*v'));
        end
        DL(i,j)=DL(j,i);
        Gj2i=Gj*Rj2i;
        Gi2j=Gi*Rj2i';
        DE(j,i)=sqrt(mean(mean((Gi-Gj2i).^2)));
        DE(i,j)=sqrt(mean(mean((Gj-Gi2j).^2)));
        % For registering Gj to Gi, s.t. Gj*R*S approximates Gi
        if (Scale)
            %% Scale subj j to subj i -- This is least squares optimal, but
            %% badly behaved (shrinks to zero) when match is not good
            % Sj2i=diag(diag(Gj2i'*Gi)./sum(Gj2i.^2)');
            % Si2j=diag(diag(Gi2j'*Gj)./sum(Gi2j.^2)');

            % Match marginal variances of each gradient component after rotation
            Sj2i=diag(sqrt(sum(Gi.^2)./sum(Gj2i.^2)+eps));
            Si2j=diag(sqrt(sum(Gj.^2)./sum(Gi2j.^2)+eps));
        else
            Sj2i=eye(nG);
            Si2j=eye(nG);
        end
        S{j,i}=Sj2i;
        S{i,j}=Si2j;
    end
end

% Apply registration, compute mean gradients (mG)
% Note SD
[sdGm,sdGF,sdGL,sdGE]=deal(nan(nG,nSubj,2)); % 3rd dim: 1 - SD before scaling; 2 - SD after scaling
[sdmGm,sdmGF,sdmGL,sdmGE]=deal(nan(nG,2));   % 2nd dim: 1-  SD of target; 2 - SD of atlas mean

% Ref: mean (after alignment to 1st subject)
mGrm=zeros(nV,nG);
for i=1:nSubj
    tG=squeeze(G(:,1:nG,i))*Rm{i};
    sdGm(:,i,1)=std(tG);
    mGrm=mGrm+tG;
end
mGrm=mGrm/nSubj;
sdmGm(:,1)=std(G(:,1:nG,1));
sdmGm(:,2)=std(mGrm);

% Ref: Best Frobenius subject
[~,iRefF]=min(sum(DF));
mGrF=zeros(nV,nG);
for j=1:nSubj
    tG=G(:,1:nG,j)*R{j,iRefF};
    sdGF(:,j,1)=std(tG);
    tG=tG*S{j,iRefF};
    sdGF(:,j,2)=std(tG);
    mGrF=mGrF+tG;
end
mGrF=mGrF/nSubj;
sdmGF(:,1)=std(G(:,1:nG,iRefF));
sdmGF(:,2)=std(mGrF);
mGrF=mGrF*R{iRefF,1}*S{iRefF,1}; % rotate/scale to subject 1

% Ref: Best Lie subject
[~,iRefL]=min(sum(DL));
mGrL=zeros(nV,nG);
for j=1:nSubj
    tG=G(:,1:nG,j)*R{j,iRefL};
    sdGL(:,j,1)=std(tG);
    tG=tG*S{j,iRefL};
    sdGL(:,j,2)=std(tG);
    mGrL=mGrL+tG;
end
mGrL=mGrL/nSubj;
sdmGL(:,1)=std(G(:,1:nG,iRefL));
sdmGL(:,2)=std(mGrL);
mGrL=mGrL*R{iRefL,1}*S{iRefL,1}; % rotate/scale to subject 1

% Ref: Best Euclidean subject
[~,iRefE]=min(sum(DE));
mGrE=zeros(nV,nG);
for j=1:nSubj
    tG=G(:,1:nG,j)*R{j,iRefE};
    sdGE(:,j,1)=std(tG);
    tG=tG*S{j,iRefE};
    sdGE(:,j,2)=std(tG);
    mGrE=mGrE+tG;
end
mGrE=mGrE/nSubj;
sdmGE(:,1)=std(G(:,1:nG,iRefE));
sdmGE(:,2)=std(mGrE);
mGrE=mGrE*R{iRefE,1}*S{iRefE,1}; % rotate/scale to subject 1

LS1={'LineStyle','-'}; % Line Spec: solid line style 

if (Scale)
    I=1:nSubj;
    g=1;
    subplot(4,1,1)
    plot(I,sdGm(g,:,1),'o',I,sdGm(g,:,2),"diamond");xlabel('Subject');ylabel('SD (o prescaling)')
    abline('h',sdmGm(g,1));abline('h',sdmGm(g,2),LS1{:});title(sprintf('Frob G%d SD (- - target, -- atlas)',g))
    subplot(4,1,2)
    plot(I,sdGF(g,:,1),'o',I,sdGF(g,:,2),"diamond");xlabel('Subject');ylabel('SD (o prescaling)')
    abline('h',sdmGF(g,1));abline('h',sdmGL(g,2),LS1{:});title(sprintf('Frob G%d SD (- - target, -- atlas)',g))
    subplot(4,1,3)
    plot(I,sdGL(g,:,1),'o',I,sdGL(g,:,2),"diamond");xlabel('Subject');ylabel('SD (o prescaling)')
    abline('h',sdmGL(g,1));abline('h',sdmGL(g,2),LS1{:});title(sprintf('Lie G%d SD (- - target, -- atlas)',g))
    subplot(4,1,4)
    plot(I,sdGE(g,:,1),'o',I,sdGE(g,:,2),"diamond");xlabel('Subject');ylabel('SD (o prescaling)')
    abline('h',sdmGE(g,1));abline('h',sdmGE(g,2),LS1{:});title(sprintf('Eucl G%d SD (- - target, -- atlas)',g))

    MyPrint(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-sds'],'v')
end

YF = squareform(DF);ZF = linkage(YF, 'average');IF = optimalleaforder(ZF, YF);
YL = squareform(DL);ZL = linkage(YL, 'average');IL = optimalleaforder(ZL, YL);
YE = squareform(DE);ZE = linkage(YE, 'average');IE = optimalleaforder(ZE, YE);

clf
subplot(3,2,1);imagesc(DF(IF,IF));colorbar;title('Frob. Dist. (hclus ord)');axis image
subplot(3,2,3);imagesc(DL(IL,IL));colorbar;title('Lie alg. Dist. (hclus ord)');axis image
subplot(3,2,5);imagesc(DE(IE,IE));colorbar;title('Eucl. Dist. (hclus ord)');axis image
subplot(3,2,2);plot(sum(DF(:,IF)));
[m,I]=min(sum(DF(:,IF)));abline('h',m,LS1{:});line(I,m,'marker','o')
Is=1:nSubj;Is=Is(IF);title(sprintf('Frob ref: subj %d (- - subj 1)',Is(I)));abline('v',find(IF==1));xlabel('hclus ord')
subplot(3,2,4);plot(sum(DL(:,IL)));
[m,I]=min(sum(DL(:,IL)));abline('h',m,LS1{:});line(I,m,'marker','o');
Is=1:nSubj;Is=Is(IL);title(sprintf('Lie ref: subj %d (- - subj 1)',Is(I)));abline('v',find(IL==1));xlabel('hclus ord')
subplot(3,2,6);plot(sum(DE(:,IE)));
[m,I]=min(sum(DE(:,IE)));abline('h',m,LS1{:});line(I,m,'marker','o');
Is=1:nSubj;Is=Is(IE);title(sprintf('Eucl ref: subj %d (- - subj 1)',Is(I)));abline('v',find(IE==1));xlabel('hclus ord')

MyPrint(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-sim'],'v')

clf
MyPlotMatrix(G(:,1:nG,1),'Mean ref (1)',AxL)
MyPrint(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-tm'],'h')

MyPlotMatrix(G(:,1:nG,iRefF)*R{iRefF,1}*S{iRefF,1},sprintf('Frob ref (%d)',iRefF),AxL)
MyPrint(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-tF'],'h')

MyPlotMatrix(G(:,1:nG,iRefL)*R{iRefF,1}*S{iRefF,1},sprintf('Lie ref (%d)',iRefL),AxL)
MyPrint(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-tL'],'h')

MyPlotMatrix(G(:,1:nG,iRefE)*R{iRefE,1}*S{iRefE,1},sprintf('Eucl ref (%d)',iRefE),AxL)
MyPrint(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-tE'],'h')

clf
MyPlotMatrix(mGrm(:,1:nG),'Reg Mean ref',AxL)
MyPrint(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-Rm'],'h')

MyPlotMatrix(mGrF(:,1:nG),'Reg Frob ref',AxL)
MyPrint(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-RF'],'h')

MyPlotMatrix(mGrL(:,1:nG),'Reg Lie ref',AxL)
MyPrint(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-RL'],'h')

MyPlotMatrix(mGrE(:,1:nG),'Reg Eucl ref',AxL)
MyPrint(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-RE'],'h')
clf

subplot(3,1,1); h=plot(std(mGrm),'color',[0.7 0.7 0.7]);title('sd(G) Reg Mean ref')
xlabel('Gradient Number');ylabel('sd(G)');ThickLine
hold on 
plot(std(mGr0),'o','color',[.7 .7 .7]);uistack(h,'top')
hold off

subplot(3,1,2); h=plot(std(mGrF));title('sd(G) Reg Frob ref')
xlabel('Gradient Number');ylabel('sd(G)');ThickLine
hold on 
plot(std(mGr0),'o','color',[.7 .7 .7]);uistack(h,'top')
hold off

subplot(3,1,3); h=plot(std(mGrL));title('sd(G) Reg Lie ref')
xlabel('Gradient Number');ylabel('sd(G)');ThickLine
hold on 
plot(std(mGr0),'o','color',[.7 .7 .7]);uistack(h,'top')
hold off

SetAllAxLim
MyPrint(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',nSubj,nG,Pow,Scale),[Out '-sd'],'v')

warning(warnState.state, 'MATLAB:logm:nonPosRealEig');

end

function MyPrint(Titl,PDFnm,Orient)
    axes('position',[0 0 1 1],'visible','off');
    if (Orient=='h')
        fullpage('h')
        axes('position',[0 0 1 1],'visible','off');
        text(0.05,0.5,...
             sprintf(Titl),...
                     'rotation',90,...
                     'fontsize',20,...
                     'horizontalalign','center',...
                     'verticalalign','middle')
    else
        fullpage(Orient)
        text(0.5,0.96,...
             sprintf(Titl),...
                     'rotation',0,...
                     'fontsize',20,...
                     'horizontalalign','center',...
                     'verticalalign','top')
    end
    print('-dpdf',[PDFnm '.pdf'])
end

function MyPlotMatrix(G,Titl,Ax)
    nG=size(G,2);
    for g1=1:nG-1
        for g2=g1+1:nG
            i=(nG-1)*(nG-g2)+g1;
            subplot(nG-1,nG-1,i)
            dscatter(G(:,g1),G(:,g2),'marker','.');
            xlabel(sprintf('G%d - sd %g',g1,std(G(:,g1))))
            ylabel(sprintf('G%d - sd %g',g2,std(G(:,g2))))
            title(Titl)
            %axis equal
            set(gca,'xlim',Ax{1},'ylim',Ax{2})
        end
    end
    if nG==3
        g1=3;g2=2;
        subplot(2,2,4)
        dscatter(G(:,g1),G(:,g2),'marker','.');
        xlabel(sprintf('G%d - sd %g',g1,std(G(:,g1))))
        ylabel(sprintf('G%d - sd %g',g2,std(G(:,g2))))
        title(Titl)
        %axis equal
        set(gca,'xlim',Ax{1},'ylim',Ax{2})
    end
    SetAllAxLim
end
