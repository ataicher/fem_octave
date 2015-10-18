phiVec = [1e-2 5e-3 1e-3 5e-4 1e-4 5e-5 1e-5 5e-6];
% phiVec = [1e-2 1e-5];
setPhi = [5];
% setPhi = [3 4 5];
formulation = [1 2 3 4];
% formulation = [4];
nx = 20;
cc = hsv(length(formulation));

    figure, clf
j=1;
for SETPHI = setPhi

    subplot(1,4*length(setPhi),4*(j-1)+(1:4))

    i=1;
    for FORMULATION = formulation
        condJacVals = [];
        for phi_m = phiVec
            [~,~,~,appCtx] = exec(@compactingColumn,0,0,nx,1,SETPHI,phi_m,phiVec(1),FORMULATION);
            condJacVals = [condJacVals appCtx.condJac];
        end
        loglog(1./phiVec,condJacVals,'color',cc(i,:),'LineWidth',3,'Marker','x','MarkerSize', 20)
        i = i+1;
        hold on
    end

    ax = gca;
    set(ax,'Box','on','FontSize',15);
    xlabel('1/\phi','FontSize',29)
    ylabel('condition number','FontSize',29)
    legend('standard', 'expanded','symmetric','scaled')
    
%     % plot porosity
%     subplot(1,4*length(setPhi),4*j)
%     xGrid = -.2:.001:.2;
%     phiVal = zeros(size(xGrid));
%     for i = 1:length(xGrid)
%         phiVal(i) = appCtx.phi(xGrid(i),1);
%     end
%     plot(xGrid, phiVal, 'LineWidth', 4)
%     view(90,90)
%     ax = gca;
%     set(ax,'Box','on','FontSize',15);
%     xlabel('depth')
%     title('porosity','FontSize', 21)
    j=j+1;
end