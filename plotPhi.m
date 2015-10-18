function plotPhi(phi)

xGrid = 0:.1:LX;
yGrid = 0:.1:LY;
phiVal = zeros(length(xGrid),length(yGrid));
for j = 1:length(yGrid)
    for i = 1:length(xGrid)
        phiVal(i,j) = phi(xGrid(i),yGrid(j));
    end
end
[X,Y] = meshgrid(xGrid,yGrid);
surf(X,Y,phiVal')
xlabel('x')
ylabel('y')