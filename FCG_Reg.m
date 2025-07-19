function FCG_Reg(Dat,varargin)
% ReadHCP(Dat,[Subjs,Gs,Pow,Scale])
%
%  Dat is 'HCP' or 'SCZ'
%  Subjs is vector of subject indicies to load (or empty for all)
%  Gs is is vector of gradient indicies to load (or empty for all)
%  Pow - 1 to weight by power of vox/vertex length, 0 no weighting
%  Scale - 1 to include axis scaling to match template, 0 no scaling    

nScaleIter=3;

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

[G,ID,AxL] = FCG_Load(Dat,Subjs,Gs);
nV=size(G,1);
nG=size(G,2);
nSubj=size(G,3);

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
R=cell(nSubj,nSubj);R(:,:)={eye(nG)}; % R{i,j} rotation of i to match j
S=cell(nSubj,nSubj);S(:,:)={ones(1,nG)};
DF=zeros(nSubj,nSubj); % Frobenius distance, orthogonal distance
DL=zeros(nSubj,nSubj); % Lie algebra norm distance (with flip in case det=-1)
DE=zeros(nSubj,nSubj); % Eucledian distance (optimised cost function)
if Scale
    DEs=zeros(nSubj,nSubj,nScaleIter+1);
end

% Standard method - Register all to 1st, make mean, register to mean
R0   = FCG_OrthProcReg(G,1,0);
% Treat reorientation to subj 1 as a preprocessing step
G    = FCG_ApplyReg(G,R0);
mG0  = mean(G,3);
Rm   = FCG_OrthProcReg(G,mG0,0);

% Estimate All Pairwise Registrations
for i = 1:nSubj
    fprintf('* -> %d\n',i)
    Gi = G(:,:,i);
    for j = i+1:nSubj
        % Register subj j to subj i
        % Solve rotation-only Procrustes: R = argmin || G_i R - G_j ||_F
        Gj = G(:,:,j);
        R{j,i} = FCG_OrthProc(Gj,Gi,Pow);
        R{i,j} = R{j,i}';
        if Scale
            % Solve rotation-scale Procrustes: R = argmin || G_i R S - G_j ||_F
            % ... but note registration no longer symmetric
            DEs(j,i,1) = NormEucl(Gj*R{j,i}-Gi);
            DEs(i,j,1) = NormEucl(Gi*R{i,j}-Gj);
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
                DEs(j,i,k+1) = NormEucl(Gj2i-Gi);
                DEs(i,j,k+1) = NormEucl(Gi2j-Gj);
            end
            R{j,i}=Rj2i;
            R{i,j}=Ri2j;
            S{j,i}=Sj2i;
            S{i,j}=Si2j;
        end
        DE(j,i)=NormEucl(Gj*R{j,i} - Gi);
        DE(i,j)=NormEucl(Gi*R{i,j} - Gj);
        [DF(j,i),DF(j,i)]=deal(ONormFrob(R{j,i}));
        [DL(i,j),DL(j,i)]=deal(ONormLie(R{j,i}));
    end
end

% Apply registration, compute:
%   mG* - mean of registered gradients (nV x nG)
%   sdG* - stdev of registered gradients (nG x nSubj x 2; (:,:,1) sd before scaling, (:,:,2) after
%   sdmG* - stdev of refs (nG x 2; (:,1) sd of target, (:,2) sd of atlas mean

% Ref: mean (after alignment to 1st subject)
[mGm,sdGm] = FCG_AtlasMean(G,Rm);
sdmGm = [std(G(:,:,1)); std(mGm)];

% Ref: Best Frobenius subject
[~,iRefF]=min(sum(DF));
[mGF,sdGF] = FCG_AtlasMean(G,R(:,iRefF));
sdmGF = [std(G(:,:,iRefF)); std(mGF)];
mGF=mGF*R{iRefF,1}.*S{iRefF,1}; % rotate/scale to subject 1

% Ref: Best Lie subject
[~,iRefL]=min(sum(DL));
[mGL,sdGL] = FCG_AtlasMean(G,R(:,iRefL));
sdmGL = [std(G(:,:,iRefL)); std(mGL)];
mGL=mGL*R{iRefL,1}.*S{iRefL,1}; % rotate/scale to subject 1

% Ref: Best Euclidean subject
[~,iRefE]=min(sum(DE));
[mGE,sdGE] = FCG_AtlasMean(G,R(:,iRefE));
sdmGE = [std(G(:,:,iRefE)); std(mGE)];
mGE=mGE*R{iRefE,1}.*S{iRefE,1}; % rotate/scale to subject 1

if (Scale)
    % Ref: Best Euclidean subject with scaling
    [~,iRefEs]=min(sum(DEs(:,:,end))+sum(DEs(:,:,end),2)');
    [mGEs,sdGEs] = FCG_AtlasMean(G,R(:,iRefEs),S);
    sdmGEs = [std(G(:,:,iRefEs)); std(mGEs)];
    mGEs=mGEs*R{iRefEs,1}.*S{iRefEs,1}; % rotate/scale to subject 1
end


%
% Plots
%

%% Interactive plot
%densityScatter(mGE');xlabel('G1');ylabel('G2');zlabel('G3'); axis([-8 8 -4 4 -4 4])

LS1={'LineStyle','-'}; % Line Spec: solid line style 

YF = squareform(DF);ZF = linkage(YF, 'average');IF = optimalleaforder(ZF, YF);
YL = squareform(DL);ZL = linkage(YL, 'average');IL = optimalleaforder(ZL, YL);
YE = squareform(DE);ZE = linkage(YE, 'average');IE = optimalleaforder(ZE, YE);
if Scale
    DEsy = (DEs(:,:,end)+DEs(:,:,end)')/2;
    YEs = squareform(DEsy);ZEs = linkage(YEs, 'average');IEs = optimalleaforder(ZEs, YEs);
    nP=4;
else
    nP=3;
end

clf
subplot(nP,2,1);imagesc(DF(IF,IF));colorbar;title('Frob. Dist. (hclus ord)');axis image
subplot(nP,2,3);imagesc(DL(IL,IL));colorbar;title('Lie alg. Dist. (hclus ord)');axis image
subplot(nP,2,5);imagesc(DE(IE,IE));colorbar;title('Eucl. Dist. (hclus ord)');axis image
if Scale
    subplot(nP,2,7);imagesc(DEs(IEs,IEs));colorbar;title('Eucl. scal Dist. (hclus ord)');axis image
end
subplot(nP,2,2);plot(sum(DF(:,IF)));
[m,I]=min(sum(DF(:,IF)));abline('h',m,LS1{:});line(I,m,'marker','o')
Is=1:nSubj;Is=Is(IF);title(sprintf('Frob ref: subj %d (- - subj 1)',Is(I)));abline('v',find(IF==1));xlabel('hclus ord')
subplot(nP,2,4);plot(sum(DL(:,IL)));
[m,I]=min(sum(DL(:,IL)));abline('h',m,LS1{:});line(I,m,'marker','o');
Is=1:nSubj;Is=Is(IL);title(sprintf('Lie ref: subj %d (- - subj 1)',Is(I)));abline('v',find(IL==1));xlabel('hclus ord')
subplot(nP,2,6);plot(sum(DE(:,IE)));
[m,I]=min(sum(DE(:,IE)));abline('h',m,LS1{:});line(I,m,'marker','o');
Is=1:nSubj;Is=Is(IE);title(sprintf('Eucl ref: subj %d (- - subj 1)',Is(I)));abline('v',find(IE==1));xlabel('hclus ord')
if Scale
    subplot(nP,2,8);plot(sum(DEsy(:,IE)));
    [m,I]=min(sum(DEsy(:,IE)));abline('h',m,LS1{:});line(I,m,'marker','o');
    Is=1:nSubj;Is=Is(IE);title(sprintf('Eucl scal ref: subj %d (- - subj 1)',Is(I)));abline('v',find(IE==1));xlabel('hclus ord')
end

FCG_Print(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-sim'],'v')

if (Scale)
    I=1:nSubj;
    g=1;
    subplot(4,1,1)
    plot(I,sdGm(g,IE,1),'o',I,sdGm(g,IE,2),"diamond");xlabel('Subject');ylabel('SD (o prescaling)')
    abline('h',sdmGm(g,1));abline('h',sdmGm(g,2),LS1{:});title(sprintf('Mean G%d SD (- - target, -- atlas)',g))
    subplot(4,1,2)
    plot(I,sdGF(g,IE,1),'o',I,sdGF(g,IE,2),"diamond");xlabel('Subject');ylabel('SD (o prescaling)')
    abline('h',sdmGF(g,1));abline('h',sdmGL(g,2),LS1{:});title(sprintf('Frob G%d SD (- - target, -- atlas)',g))
    subplot(4,1,3)
    plot(I,sdGL(g,IE,1),'o',I,sdGL(g,IE,2),"diamond");xlabel('Subject');ylabel('SD (o prescaling)')
    abline('h',sdmGL(g,1));abline('h',sdmGL(g,2),LS1{:});title(sprintf('Lie G%d SD (- - target, -- atlas)',g))
    subplot(4,1,4)
    plot(I,sdGE(g,IE,1),'o',I,sdGE(g,IE,2),"diamond");xlabel('Subject');ylabel('SD (o prescaling)')
    abline('h',sdmGE(g,1));abline('h',sdmGE(g,2),LS1{:});title(sprintf('Eucl G%d SD (- - target, -- atlas)',g))
    subplot(5,1,5)
    plot(I,sdGEs(g,IE,1),'o',I,sdGEs(g,IE,2),"diamond");xlabel('Subject');ylabel('SD (o prescaling)')
    abline('h',sdmGEs(g,1));abline('h',sdmGEs(g,2),LS1{:});title(sprintf('Eucl-scaled G%d SD (- - target, -- atlas)',g))

    FCG_Print(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-sds'],'v')
end



clf
FCG_PlotMatrix(G(:,:,1),'Mean ref (1)',AxL)
FCG_Print(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-tm'],'h')

FCG_PlotMatrix(G(:,:,iRefF)*R{iRefF,1}.*S{iRefF,1},sprintf('Frob ref (%d)',iRefF),AxL)
FCG_Print(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-tF'],'h')

FCG_PlotMatrix(G(:,:,iRefL)*R{iRefF,1}.*S{iRefF,1},sprintf('Lie ref (%d)',iRefL),AxL)
FCG_Print(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-tL'],'h')

FCG_PlotMatrix(G(:,:,iRefE)*R{iRefE,1}.*S{iRefE,1},sprintf('Eucl ref (%d)',iRefE),AxL)
FCG_Print(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-tE'],'h')

FCG_PlotMatrix(G(:,:,iRefE)*R{iRefE,1}.*S{iRefE,1},sprintf('Eucl ref (%d)',iRefE),AxL)
FCG_Print(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-tE'],'h')

if Scale
    FCG_PlotMatrix(G(:,:,iRefEs)*R{iRefEs,1}.*S{iRefEs,1},sprintf('Eucl scal ref (%d)',iRefEs),AxL)
    FCG_Print(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-tEs'],'h')
end

clf
FCG_PlotMatrix(mGm(:,1:nG),'Reg Mean ref',AxL)
FCG_Print(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-Rm'],'h')

FCG_PlotMatrix(mGF(:,1:nG),'Reg Frob ref',AxL)
FCG_Print(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-RF'],'h')

FCG_PlotMatrix(mGL(:,1:nG),'Reg Lie ref',AxL)
FCG_Print(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-RL'],'h')

FCG_PlotMatrix(mGE(:,1:nG),'Reg Eucl ref',AxL)
FCG_Print(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-RE'],'h')

if Scale
    FCG_PlotMatrix(mGEs(:,1:nG),'Reg Eucl scal ref',AxL)
    FCG_Print(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',Dat,nSubj,nG,Pow,Scale),[Out '-REs'],'h')
end

clf

subplot(3,1,1); h=plot(std(mGm),'color',[0.7 0.7 0.7]);title('sd(G) Reg Mean ref')
xlabel('Gradient Number');ylabel('sd(G)');ThickLine
hold on 
plot(std(mG0),'o','color',[.7 .7 .7]);uistack(h,'top')
hold off

subplot(3,1,2); h=plot(std(mGF));title('sd(G) Reg Frob ref')
xlabel('Gradient Number');ylabel('sd(G)');ThickLine
hold on 
plot(std(mG0),'o','color',[.7 .7 .7]);uistack(h,'top')
hold off

subplot(3,1,3); h=plot(std(mGL));title('sd(G) Reg Lie ref')
xlabel('Gradient Number');ylabel('sd(G)');ThickLine
hold on 
plot(std(mG0),'o','color',[.7 .7 .7]);uistack(h,'top')
hold off

SetAllAxLim
FCG_Print(sprintf('%s FCGrad  %d Subj  %d Grad  %d Pow  %d Scale',nSubj,nG,Pow,Scale),[Out '-sd'],'v')


end



function F = ONormFrob(R)
% Frobenius norm for orthogonal matricies
    F = norm( R - eye(size(R)), 'fro');
end

function L = ONormLie(R)
% Lie algebra norm for orthogonal matricies

    % Supress warnings about non-principal logm, due to R \in O(K) instead of SO(K)
    warnState = warning('query', 'MATLAB:logm:nonPosRealEig');
    warning('off', 'MATLAB:logm:nonPosRealEig');

    if det(R)>0
        L=norm(logm(R));
    else
        % Project into SO(n) to avoid logm error on det=-1
        [u,~,v]=svd(R);
        u(:,end)=-u(:,end);
        L=norm(logm(u*v'));
    end

    warning(warnState.state, 'MATLAB:logm:nonPosRealEig');

end

function E = NormEucl(Err)
% Mean Eucledian norm on registration errors
    E = mean(sqrt(sum(Err.^2,2)));
end

