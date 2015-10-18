function computeStreamLines(u,appCtx)

RT0VelocityField = 1;
cells = appCtx.mesh.cells;
edges = appCtx.mesh.edges;
nodes = appCtx.mesh.nodes;
edgeCells = appCtx.mesh.edgeCells;
bndryEdges = appCtx.mesh.bndryEdges;
numBndryEdges = appCtx.mesh.numBndryEdges;
numLocalEdges = appCtx.mesh.numLocalEdges;
v_0Ref = appCtx.quad.v_0Ref;
refNodes = appCtx.quad.refNodes;
numComp = appCtx.field(RT0VelocityField).numComp;
lNormals = [0 -1 ; 1 0 ; 0 1 ; -1 0];

x = [2*ones(1,10) ; .11:.2:1.91];
figure, clf
for i = 1:10
    x_in = x(:,i);
    
    for bEPtr = 1:numBndryEdges
        ePtr = bndryEdges(bEPtr);
        n1Ptr = edges(1,ePtr);
        n2Ptr = edges (2,ePtr);
        n1 = nodes(:,n1Ptr);
        n2 = nodes(:,n2Ptr);
        if sign(x_in(1)-n1(1))*sign(x_in(1)-n2(1)) <= 0
            if sign(x_in(2)-n1(2))*sign(x_in(2)-n2(2)) <= 0
                c = edgeCells(1,ePtr);
                if (sign(x_in(2)-n1(2))*sign(x_in(2)-n2(2)) == 0) &&  (sign(x_in(1)-n1(1))*sign(x_in(1)-n2(1)) == 0)
                    error('cant start at a node')
                end
            end
        end
    end
    
    path = x_in;
    iter = 0; maxIter = 100;
    while c > 0 && iter < maxIter
        
        J = appCtx.cellGeometry.J{c};
        invJ = appCtx.cellGeometry.invJ{c};
        detJ = appCtx.cellGeometry.detJ{c};
        v_0 = appCtx.cellGeometry.v_0{c};
        cellDOF = appCtx.field(RT0VelocityField).cellDOF(:,c);
        
        x_in_ref = invJ*(x_in - v_0) + v_0Ref;
        
        uLocal = u(cellDOF);
        uLocalRef = zeros(size(uLocal));
        for lEPtr=1:numLocalEdges
            orientation = sign(cells(lEPtr,c));
            for comp = 1:numComp
                uLocalRef(lEPtr) = uLocalRef(lEPtr) + abs(lNormals(lEPtr,comp))*detJ*invJ(comp,comp)*orientation*uLocal(lEPtr);
            end
        end
        
        A = (uLocalRef(2)+uLocalRef(4))/4; B = (uLocalRef(2)-uLocalRef(4))/4;
        C = (uLocalRef(1)+uLocalRef(3))/4; D = (uLocalRef(3)-uLocalRef(1))/4;
        
        if A == 0
            CInt(1) = x_in_ref(1);
        else
            CInt(1) = x_in_ref(1) + B/A;
        end
        if C == 0
            CInt(2) = x_in_ref(2);
        else
            CInt(2) = x_in_ref(2) + D/C;
        end
        
        if C == 0
            if D == 0
                t(1) = 0;
            else
                t(1) = (-1-CInt(2))/D;
            end
        else
            if (-1+D/C)/CInt(2) > 0
                t(1) = (1/C)*log((-1+D/C)/CInt(2));
            else
                t(1) = 0;
            end
        end
        
        if A == 0
            if B == 0
                t(2) = 0;
            else
                t(2) = (1-CInt(1))/B;
            end
        else
            if (1+B/A)/CInt(1) > 0
                t(2) = (1/A)*log((1+B/A)/CInt(1));
            else
                t(2) = 0;
            end
        end
        
        if C == 0
            if D == 0
                t(3) = 0;
            else
                t(1) = (1-CInt(2))/D;
            end
        else
            if (1+D/C)/CInt(2) > 0
                t(3) = (1/C)*log((1+D/C)/CInt(2));
            else
                t(3) = 0;
            end
        end
        
        if A == 0
            if B == 0
                t(4) = 0;
            else
                t(4) = (-1-CInt(1))/B;
            end
        else
            if (-1+B/A)/CInt(1) > 0
                t(4) = (1/A)*log((-1+B/A)/CInt(1));
            else
                t(4) = 0;
            end
        end
        if isempty(t(t>0))
            break
        end
        tMin = min(t(t>0));
        lEPtr = find(t==tMin);
        ePtr = abs(cells(lEPtr,c));
        c = edgeCells((edgeCells(:,ePtr)~=c),ePtr);
        
        %         if A == 0
        x_in_ref(1) = -B/A + CInt(1)*exp(A*tMin);
        %         if abs(x_in_ref(1) - 1) < 1e-5
        %             x_in_ref(1) = 1;
        %         elseif abs(x_in_ref(1) + 1) < 1e-5
        %             x_in_ref(1) = -1;
        %         end
        x_in_ref(2) = -D/C + CInt(2)*exp(C*tMin);
        %         if abs(x_in_ref(2) - 1) < 1e-5
        %             x_in_ref(2) = 1;
        %         elseif abs(x_in_ref(2) + 1) < 1e-5
        %             x_in_ref(2) = -1;
        %         end
        if ~isempty(find(sum(refNodes-repmat(x_in_ref,1,4))==0))
            error('node is a corner')
        end
        
        x_in = J*(x_in_ref - v_0Ref) + v_0;
        path = [path x_in];
        iter = iter + 1;
        
    end
    
    plot(path(1,:),path(2,:)), hold on
end
