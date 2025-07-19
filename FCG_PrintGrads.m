function FCG_PrintGrads(Dat,pdfNm,varargin)
% [G,AxL,ID]=FCG_PrintGrads(Dat,pdfNm[,nSubj,Gs,Reg])
% [G,AxL,ID]=FCG_PrintGrads(G,pdfNm[,ID,Ax,Titl])
% Print gradients
%
%  Dat is 'HCP' or 'SCZ'
%  pdfNm is output basename (no ".pdf")
%  Subjs is vector of subject indicies to load (or empty for all)
%  Gs is is a 2-vector of gradient indicies to load (or empty for [1,2])
%  Reg: 0 no registration;  1 is reg to first
%  Titl: Main title
%
%  If specifying G, must specify
%    ID - nSubj cell-vector
%    Ax - Axis limits    

nCol=6;
nRow=4;

% Parse arguments, 2 different modes
if ~isnumeric(Dat)
    if nargin<3 || isempty(varargin{1})
        Subjs=[];
    else
        Subjs=varargin{1};
    end
    if nargin<4 || isempty(varargin{2})
        Gs=[1,2];
    else
        Gs=varargin{2};
        if length(Gs)~=2
            error('Must specify exactly 2 gradients to plot')
        end
    end
    if nargin<5 || isempty(varargin{3})
        Reg=0;
    else
        Reg=varargin{1};
    end
    [G,ID,Ax]=FCG_Load(Dat,Subjs,Gs);
else
    G=Dat; clear Dat
    if size(G,2)~=2
        error('Only can plot pairs of gradients')
    end
    if nargin<3
        ID={};
    else
        ID=varargin{1};
    end
    if nargin<4
        Ax={};
    else
        Ax=varargin{2};
    end
    if nargin<5
        MainTitle='';
    else
        MainTitle=varargin{3};
    end
    Reg=[];
    Gs=1:size(G,2);
end

nSubj=size(G,3);
Gstr=cellstr(strcat('G',num2str(Gs(:))));
if isempty(ID)
    ID=cellstr(string(1:nSubj));
end
if ~isempty(Reg) && Reg
    MainTitle=sprintf('Registered to first (%s)',ID{1});
elseif ~isempty(Reg) && ~Reg
    MainTitle=sprintf('No registration');
end

g1=Gs(1);g2=Gs(2);

c=1;p=1;
TmpNm=tempname;
Tmps=cell(0);
fprintf('Page ');
for i=1:nSubj
    subplot(nRow,nCol,c);
    
    dscatter(G(:,g1,i),G(:,g2,i))
    xlabel(sprintf('G%d - sd %g',g1,std(G(:,g1,i))))
    ylabel(sprintf('G%d - sd %g',g2,std(G(:,g2,i))))
    grid on
    title(sprintf('%d: "%s"',i,ID{i}))
    if ~isempty(Ax)
        set(gca,'xlim',Ax{1},'ylim',Ax{2})
    end

    c=c+1;
    if c>nRow*nCol || i==nSubj
        Tmps{p}=[TmpNm '_' num2str(p)];
        FCG_Print(MainTitle,Tmps{p},'h')
        clf
        c=1;
        fprintf('%d ',p);
        p=p+1;
    end
end
fprintf('\n');

gsCmd=sprintf('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/ebook -sOutputFile=%s.pdf %s',pdfNm,strjoin(strcat(Tmps,'.pdf')));
if system(gsCmd)
    error("Ghostscript command failed... Is it installed?")
end
cellfun(@delete,strcat(Tmps,'.pdf'));

end


