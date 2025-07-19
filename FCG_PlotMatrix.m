function FCG_PlotMatrix(G,Titl,Ax)
    nG=size(G,2);
    for g1=1:nG-1
        for g2=g1+1:nG
            i=(nG-1)*(nG-g2)+g1;
            subplot(nG-1,nG-1,i)
            dscatter(G(:,g1),G(:,g2),'marker','.');
            xlabel(sprintf('G%d - sd %g',g1,std(G(:,g1))))
            ylabel(sprintf('G%d - sd %g',g2,std(G(:,g2))))
            grid on
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
        grid on
        title(Titl)
        %axis equal
        set(gca,'xlim',Ax{1},'ylim',Ax{2})
    end
    SetAllAxLim
end

