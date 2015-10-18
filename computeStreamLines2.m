function computeStreamLines2(u,appCtx)

dim = appCtx.dim;
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

figure, clf

% starting points for the steam lines
x = [2*ones(1,11) ; .11:.2:1.91 1.95];

for i = 1:length(x)
    x_in = x(:,i);
    
    % find initial cell
    for bEPtr = 1:numBndryEdges
        ePtr = bndryEdges(bEPtr);
        n1Ptr = edges(1,ePtr);
        n2Ptr = edges (2,ePtr);
        n1 = nodes(:,n1Ptr);
        n2 = nodes(:,n2Ptr);
        if sign(x_in(1)-n1(1))*sign(x_in(1)-n2(1)) <= 0
            if sign(x_in(2)-n1(2))*sign(x_in(2)-n2(2)) <= 0
                c = edgeCells(1,ePtr);
                if (sign(x_in(2)-n1(2))*sign(x_in(2)-n2(2)) == 0) && (sign(x_in(1)-n1(1))*sign(x_in(1)-n2(1)) == 0)
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
        
        % calculate position on reference element
        x_in_ref = invJ*(x_in - v_0) + v_0Ref;
        
        % calculate DOF on reference element [-1 1]^2
        % b1 = [0 ; (y-1)/4]
        % b2 = [(x+1)/4 ; 0]
        % b3 = [0 ; (y+1)/4]
        % b4 = [(x-1)/4 ; 0]
        uLocal = u(cellDOF);
        uLocalRef = zeros(size(uLocal));
        for lEPtr=1:numLocalEdges
            orientation = sign(cells(lEPtr,c));
            for comp = 1:numComp
                uLocalRef(lEPtr) = uLocalRef(lEPtr) + abs(lNormals(lEPtr,comp))*detJ*invJ(comp,comp)*orientation*uLocal(lEPtr);
            end
        end
        
        % v_x = A(1)*x + B(1)
        % v_x = A(2)*x + B(2)
        A(1) = (uLocalRef(2)+uLocalRef(4))/4; B(1) = (uLocalRef(2)-uLocalRef(4))/4;
        A(2) = (uLocalRef(1)+uLocalRef(3))/4; B(2) = (uLocalRef(3)-uLocalRef(1))/4;
        
        % find time to exit cell
        % if A(1) = 0
        %   x(t) = B(1)*t + x_in_ref(1)
        %   t(x) = (x - x_in_ref(1))/B(1)
        % if A(1) ~= 0
        %   x(t) = -B(1)/A(1) + (x_in_ref(1) + B(1)/A(1))*exp(-A(1)*t)
        %   t(x) = (-1/A(1))*log((x+B(1)/A(1))/(x_in_ref+B(1)/A(1)))
        % similarly for y
        t = zeros(numLocalEdges,1);
        %x-edges
        if A(1) == 0
            if B(1) ~= 0
                t(2) =  (1 - x_in_ref(1))/B(1);
                t(4) = (-1 - x_in_ref(1))/B(1);
            end
        else
            if abs(x_in_ref(1) + B(1)/A(1)) > 0
                if  x_in_ref(1) ~=  1 && (1 + B(1)/A(1))/(x_in_ref(1) + B(1)/A(1)) > 0
                    t(2) = (1/A(1))*log( (1 + B(1)/A(1))/(x_in_ref(1) + B(1)/A(1)));
                end
                if  x_in_ref(1) ~= -1 && (-1 + B(1)/A(1))/(x_in_ref(1) + B(1)/A(1)) > 0
                    t(4) = (1/A(1))*log((-1 + B(1)/A(1))/(x_in_ref(1) + B(1)/A(1)));
                end
            end
        end
        % y-edges
        if A(2) == 0
            if B(2) ~= 0
                t(3) =  (1 - x_in_ref(2))/B(2);
                t(1) = (-1 - x_in_ref(2))/B(2);
            end
        else
            if abs(x_in_ref(2) + B(2)/A(2)) > 0
                if x_in_ref(2) ~=  1 && (1 + B(2)/A(2))/(x_in_ref(2) + B(2)/A(2)) > 0
                    t(3) = (1/A(2))*log( (1 + B(2)/A(2))/(x_in_ref(2) + B(2)/A(2)));
                end
                if x_in_ref(2) ~= -1 && (-1 + B(2)/A(2))/(x_in_ref(2) + B(2)/A(2)) > 0
                    t(1) = (1/A(2))*log((-1 + B(2)/A(2))/(x_in_ref(2) + B(2)/A(2)));
                end
            end
        end
        
        % check that cell is exited
        if isempty(t(t>0))
            warning('did not exit cell')   
            break
        end
        tMin = min(t(t>0));
        
        % find the edge where it exited
        lEPtr = t==tMin;
        ePtr = abs(cells(lEPtr,c));
        
        % calculate next cell. If at boundary set c = 0
        c = edgeCells((edgeCells(:,ePtr)~=c),ePtr);
        
        % calculate x_in for the next cell
        for d = 1:dim
            if A(d) == 0
                x_in_ref(d) = B(d)*tMin + x_in_ref(d);
            else
                x_in_ref(d) = -B(d)/A(d) + (x_in_ref(d) + B(d)/A(d))*exp(A(d)*tMin);
            end
        end
        x_in = J*(x_in_ref - v_0Ref) + v_0;
        
        % check that the new entry is not a node
        if ~isempty(find(sum(refNodes-repmat(x_in_ref,1,4)) == 0, 1))
            error('node is a corner')
        end
        
        % update path
        path = [path x_in];
        
        % update iteration and test for stuck stream lines
        iter = iter + 1;
        if iter == maxIter
            warning('reached max iterations')
        end
    end
    
    plot(path(1,:),path(2,:)), hold on
end