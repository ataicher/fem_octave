function outputProblemInfo(appCtx)

numCells = appCtx.mesh.numCells;
numEdges = appCtx.mesh.numEdges;
numNodes = appCtx.mesh.numNodes;
numBndryEdges = appCtx.mesh.numBndryEdges;
h = appCtx.h;
prec = appCtx.quad.prec;
DEBUGLEVEL = appCtx.DEBUGLEVEL;
NewtonAbsTol = appCtx.NewtonAbsTol;
NewtonRelTol = appCtx.NewtonRelTol;
NewtonMaxIter = appCtx.NewtonMaxIter;
numFields = appCtx.numFields;
globalSize = appCtx.globalSize;
DEGENERATE = appCtx.DEGENERATE;
if DEGENERATE
    phi = appCtx.degenerate.phi;
    RT0VelocityField = appCtx.degenerate.RT0VelocityField;
    RT0PressureField = appCtx.degenerate.RT0PressureField;
    UNDERINTEGRATE = appCtx.degenerate.UNDERINTEGRATE;
    Theta = appCtx.degenerate.Theta;
end

if DEBUGLEVEL > 1
    fprintf('\nMesh Info:\n');
    fprintf('   Number of cells: %d\n', numCells);
    fprintf('   Number of edges: %d\n', numEdges);
    fprintf('   Number of nodes: %d\n', numNodes);
    fprintf('   Number of boundary edges: %d\n', numBndryEdges);
    fprintf('\nGaussian quadrature precision: %d\n', prec);
    if appCtx.h > 0
        fprintf('\nJacobian finite difference parameter: %d\n', h);
    end
    fprintf('\nNewton method parameters:\n');
    fprintf('   absolute tolerance: %d\n', NewtonAbsTol);
    fprintf('   relative tolerance: %d\n', NewtonRelTol);
    fprintf('   max iterations: %d\n', NewtonMaxIter);
    fprintf('\nTotal number of degrees of freedom: %d\n', globalSize);
    
    if DEGENERATE
        fprintf('\nSolving degenerate problem:\n');
        fprintf('     the variational form for the vector field %d is supplemented with the term: -(u.field{%d}, phi^(-1/2) div (phi^(%d/2) v))\n',RT0VelocityField, RT0PressureField, Theta);
        fprintf('     the variational form for the scalar field %d is supplemented with the term: -(phi^(-1/2) div (phi^(%d/2) u.field{%d} , q)\n',RT0PressureField, Theta, RT0VelocityField);
        if UNDERINTEGRATE
            fprintf('     Under-integration using trapezoidal rule:\n');
            fprintf('          the variational form for the vector field %d is supplemented with the term: (u.field{%d}, v))_Q\n',RT0VelocityField,RT0VelocityField);
        end
    end
            
    for f = 1:numFields
        name = appCtx.field(f).name;
        numComp = appCtx.field(f).numComp;
        if appCtx.field(f).numComp == 1
            fprintf('\nfield %d: %s (scalar valued)\n', f, name);
        else
            fprintf('\nfield %d: %s (vector valued)\n', f, name);
        end
        fprintf('\n   Finite element Basis: %s\n', appCtx.field(f).basisType);
        fprintf('      Number of basis functions: %d\n', appCtx.field(f).numBasisFuncs);
        if appCtx.field(f).NULLSPACE
            fprintf('\n   Non-zero nullspace set\n')
        end
        
        if appCtx.field(f).BNDRYCONDTION
            fprintf('\n   boundary condition:\n')
            numBC = appCtx.field(f).numBC;
            for bc = 1:numBC
                alpha = appCtx.field(f).bndry(bc).alpha;
                beta = appCtx.field(f).bndry(bc).beta;
                eta = appCtx.field(f).bndry(bc).eta;
                
                fprintf('      boundary %d:\n',bc);
                if numComp == 1
                    
                    str = func2str(eta);
                    eta = str(7:end);
                    if alpha ~= 0
                        if beta == 0
                            if alpha == 1
                                fprintf('         (natural condition) = %s\n', eta);
                            else
                                fprintf('         %d(natural condition) = %s\n', alpha, eta);
                            end
                        else
                            if alpha == 1 && beta == 1
                                fprintf('         (natural condition) + p = %s\n', eta);
                            elseif alpha ~= 1 && beta == 1
                                fprintf('         %d(natural condition) + p = %s\n', alpha, eta);
                            elseif alpha == 1 && beta ~= 1
                                fprintf('         (natural condition) + %d*p = %s\n', beta, eta);
                            else
                                fprintf('         %d(natural condition) + %d*p = %s\n', alpha, beta, eta);
                            end
                        end
                    else
                        if beta == 1
                            fprintf('         p = %s\n', eta);
                        else
                            fprintf('         %d*p = %s\n', beta, eta);
                        end
                    end
                    
                elseif numComp == 2
                    
                    str = func2str(eta{1});
                    eta{1} = str(7:end);
                    str = func2str(eta{2});
                    eta{2} = str(7:end);
                    
                    if alpha(1) ~= 0
                        if beta(1) == 0
                            if alpha(1) == 1
                                fprintf('         n^T*(natural condition) = %s\n', eta{1});
                            else
                                fprintf('         %3.fn^T*(natural condition)*n = %s\n', alpha(1), eta{1});
                            end
                        else
                            if alpha(1) == 1 && beta(1) == 1
                                fprintf('         n^T*(natural condition) + u*n = %s\n', eta{1});
                            elseif alpha(1) ~= 1 && beta(1) == 1
                                fprintf('         %3.fn^T*(natural condition) + u*n = %s\n', alpha(1), eta{1});
                            elseif alpha(1) == 1 && beta(1) ~= 1
                                fprintf('         n^T*(natural condition) + %3.f*u*n = %s\n', beta(1), eta{1});
                            else
                                fprintf('         %.3fn^T*(natural condition) + %3.f*u*n = %s\n', alpha(1), beta(1), eta{1});
                            end
                        end
                    else
                        if beta(1) == 1
                            fprintf('         u*n = %s\n', eta{1});
                        else
                            fprintf('         %3f*u*n = %s\n', beta(1), eta{1});
                        end
                    end
                    
                    if alpha(2) ~= 0
                        if beta(2) == 0
                            if alpha(2) == 1
                                fprintf('         tau^T*(natural condition) = %s\n', eta{2});
                            else
                                fprintf('         %4.1ftau^T*(natural condition) = %s\n', alpha(2), eta{2});
                            end
                        else
                            if alpha(2) == 1 && beta(2) == 1
                                fprintf('         tau^T*(natural condition) + u*tau = %s\n', eta{2});
                            elseif alpha(2) ~= 1 && beta(2) == 1
                                fprintf('         %4.1ftau^T*(natural condition) + u*tau = %s\n', alpha(2), eta{2});
                            elseif alpha(2) == 1 && beta(2) ~= 1
                                fprintf('         tau^T*(natural condition) + %4.1f*u*tau = %s\n', beta(2), eta{2});
                            else
                                fprintf('         %4.1ftau^T*(natural condition) + %4.1f*u*tau = %s\n', alpha(2), beta(2), eta{2});
                            end
                        end
                    else
                        if beta(2) == 1
                            fprintf('         u*tau = %s\n', eta{2});
                        else
                            fprintf('         %4.1f*u*tau = %s\n', beta(2), eta{2});
                        end
                    end
                end
            end
        end
        
        %         if EXISTEXACTSOL
        %             fprintf('\n   Exact solution:\n')
        %             for comp = 1:numComp
        %                 fprintf('      %s\n', func2str(appCtx.field(f).uExact{comp}));
        %             end
        %             if EXISTEXACTSOLGRAD
        %                 fprintf('\n   Exact solution gradient:\n')
        %                 for comp = 1:numComp
        %                     for d = 1:appCtx.dim
        %                         fprintf('      %s\n', func2str(appCtx.field(f).gradUExact{comp,d}));
        %                     end
        %                 end
        %             end
        %         end
        %         fprintf('\n   Variational Form:\n')
        %         fprintf('      %s\n', func2str(appCtx.field(f).varForm.v));
        %         fprintf('      %s\n', func2str(appCtx.field(f).varForm.gradV));
        %         if appCtx.h == 0
        %             fprintf('   \nJacobian Form:\n')
        %         end
    end
    fprintf('\n------------------------------------------------------------\n');
    
end