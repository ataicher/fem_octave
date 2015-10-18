switch FORMULATION
    % standard (divide by phi^(2+2*Theta))
    case 1
        % ( phi^(-2-2*Theta)*u_r , v_r ) - ( q_f , div v_r )
        field(1).varForm.v{1} = @(u, gradU, x, y) (U_0/w_0)*phi(x,y)^(-2-2*Theta)*u.field{1}(1);
        field(1).varForm.v{2} = @(u, gradU, x, y) (U_0/w_0)*phi(x,y)^(-2-2*Theta)*u.field{1}(2);
        field(1).varForm.gradV{1,1} = @(u, gradU, x, y) -u.field{2}(1);
        field(1).varForm.gradV{2,2} = @(u, gradU, x, y) -u.field{2}(1);
        
        % -( div u_r, w_f ) - ( (phi/(1-phi))*(q_f-q) , w_f )
        field(2).varForm.v{1} = @(u, gradU, x, y)  -(gradU.field{1}(1,1) + gradU.field{1}(2,2)) - ...
            phi(x,y)*(u.field{2}(1)-u.field{4}(1));
        
        % -( q , divV_m ) + ( 2*(1-phi)*DU_m , DV_m ) - ( (2/3)*(1-phi)*divU_m , divV_m ) - ([(1-phi) ; 0] , v_m)
        field(3).varForm.gradV{1,1} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(1,1) - ((5-2*phi(x,y))/3)*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
        field(3).varForm.gradV{1,2} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
        field(3).varForm.gradV{2,1} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
        field(3).varForm.gradV{2,2} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(2,2) - ((5-2*phi(x,y))/3)*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
        field(3).varForm.v{1} = @(u, gradU, x, y) -(1-phi(x,y));
        
        % -( divU_m, w ) + ( (phi/(1-phi))*(q_f-q) , w )
        field(4).varForm.v{1} = @(u, gradU, x, y) -(gradU.field{3}(1,1) + gradU.field{3}(2,2)) + phi(x,y)*(u.field{2}(1)-u.field{4}(1));
        
        % auxilliary velocity
    case 2
        % ( u_r , v_r ) - ( q_f , div v_r )
        field(1).varForm.v{1} = @(u, gradU, x, y) u.field{5}(1);
        field(1).varForm.v{2} = @(u, gradU, x, y) u.field{5}(2);
        field(1).varForm.gradV{1,1} = @(u, gradU, x, y) -u.field{2}(1);
        field(1).varForm.gradV{2,2} = @(u, gradU, x, y) -u.field{2}(1);
        
        % -( div u_r, w_f ) - ( (phi/(1-phi))*(q_f-q) , w_f )
        field(2).varForm.v{1} = @(u, gradU, x, y)  -(gradU.field{1}(1,1) + gradU.field{1}(2,2)) - ...
            phi(x,y)*(u.field{2}(1)-u.field{4}(1));
        
        % -( q , divV_m ) + ( 2*(1-phi)*DU_m , DV_m ) - ( (2/3)*(1-phi)*divU_m , divV_m ) - ([(1-phi) ; 0] , v_m)
        field(3).varForm.gradV{1,1} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(1,1) - ((5-2*phi(x,y))/3)*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
        field(3).varForm.gradV{1,2} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
        field(3).varForm.gradV{2,1} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
        field(3).varForm.gradV{2,2} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(2,2) - ((5-2*phi(x,y))/3)*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
        field(3).varForm.v{1} = @(u, gradU, x, y) -(1-phi(x,y));
        
        % -( divU_m, w ) + ( (phi/(1-phi))*(q_f-q) , w )
        field(4).varForm.v{1} = @(u, gradU, x, y) -(gradU.field{3}(1,1) + gradU.field{3}(2,2)) + phi(x,y)*(u.field{2}(1)-u.field{4}(1));
        
        %(phi^(2*(1+Theta))*su_r, sv_r) - (u_r,sv_r)
        field(5).varForm.v{1} = @(u, gradU, x, y) phi(x,y)^(2*(1+Theta))*u.field{5}(1) - (U_0/w_0)*u.field{1}(1);
        field(5).varForm.v{2} = @(u, gradU, x, y) phi(x,y)^(2*(1+Theta))*u.field{5}(2) - (U_0/w_0)*u.field{1}(2);
        
        % symmetric
    case 3
        % ( u_r , v_r ) - ( q_f , div(phi^(1+Theta)*v_r )
        field(1).varForm.v{1} = @(u, gradU, x, y) (U_0/w_0)*u.field{1}(1) - (1+Theta)*phi(x,y)^Theta*dxPhi(x,y)*u.field{2}(1);
        field(1).varForm.v{2} = @(u, gradU, x, y) (U_0/w_0)*u.field{1}(2) - (1+Theta)*phi(x,y)^Theta*dyPhi(x,y)*u.field{2}(1);
        field(1).varForm.gradV{1,1} = @(u, gradU, x, y) -phi(x,y)^(1+Theta)*u.field{2}(1);
        field(1).varForm.gradV{2,2} = @(u, gradU, x, y) -phi(x,y)^(1+Theta)*u.field{2}(1);
        
        % -( div(phi^(1+Theta)*u_r), w_f ) - ( (phi/(1-phi))*(q_f-q) , w_f )
        field(2).varForm.v{1} = @(u, gradU, x, y)  -phi(x,y)^(1+Theta)*(gradU.field{1}(1,1) + gradU.field{1}(2,2)) - ...
            (1+Theta)*phi(x,y)^Theta*(dxPhi(x,y)*u.field{1}(1) + dyPhi(x,y)*u.field{1}(2)) - ...
            phi(x,y)*(u.field{2}(1)-u.field{4}(1));
        
        % -( q , divV_m ) + ( 2*(1-phi)*DU_m , DV_m ) - ( (2/3)*(1-phi)*divU_m , divV_m ) - ([(1-phi) ; 0] , v_m)
        field(3).varForm.gradV{1,1} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(1,1) - ((5-2*phi(x,y))/3)*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
        field(3).varForm.gradV{1,2} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
        field(3).varForm.gradV{2,1} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
        field(3).varForm.gradV{2,2} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(2,2) - ((5-2*phi(x,y))/3)*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
        field(3).varForm.v{1} = @(u, gradU, x, y) -(1-phi(x,y));
        
        % -( divU_m, w ) + ( (phi/(1-phi))*(q_f-q) , w )
        field(4).varForm.v{1} = @(u, gradU, x, y) -(gradU.field{3}(1,1) + gradU.field{3}(2,2)) + phi(x,y)*(u.field{2}(1)-u.field{4}(1));
        
        % scaled
    case 4
        % ( u_r , v_r ) - ( phi^(n/2)*q_f , phi^(-1/2)*div(phi^(1+Theta)*V_r) )
        if ~SETUNDERINTEGRATE
            field(1).varForm.v{1} = @(u, gradU, x, y) (U_0/w_0)*u.field{1}(1);
            field(1).varForm.v{2} = @(u, gradU, x, y) (U_0/w_0)*u.field{1}(2);
        end
        
        % ( phi^(-1/2)*div(phi^(1+Theta)*u_r), w_f ) + ( (phi/(1-phi))*(q_f-q) , w_f )
        field(2).varForm.v{1} = @(u, gradU, x, y) u.field{2}(1) - sqrt(phi(x,y))*u.field{4}(1);
        
        % -( q , divV_m ) + ( 2*(1-phi)*DU_m , DV_m ) - ( (2/3)*(1-phi)*divU_m , divV_m ) - ([(1-phi) ; 0] , v_m)
        field(3).varForm.gradV{1,1} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(1,1) - ((5-2*phi(x,y))/3)*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
        field(3).varForm.gradV{1,2} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
        field(3).varForm.gradV{2,1} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
        field(3).varForm.gradV{2,2} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(2,2) - ((5-2*phi(x,y))/3)*(1-phi(x,y))*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
        field(3).varForm.v{1} = @(u, gradU, x, y) -(1-phi(x,y));
        
        % -( divU_m, w ) + ( (phi/(1-phi))*(q_f-q) , w )
        field(4).varForm.v{1} = @(u, gradU, x, y) -(gradU.field{3}(1,1) + gradU.field{3}(2,2)) + sqrt(phi(x,y))*u.field{2}(1)-phi(x,y)*u.field{4}(1);
end