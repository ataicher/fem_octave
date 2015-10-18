C0 = -(w_0/U_0);
    C1 = -(w_0/U_0)*phi_0^(1+Theta);
    C2 = -(w_0/U_0)*phi_0^(2+2*Theta);

    if FORMULATION == 3 || FORMULATION == 4
        field(1).uExact{1} = @(x,y) C1*(4*(x^2-y^2)/(pi*(x^2+y^2)^2) + 1);
        field(1).gradUExact{1,1} = @(x,y) 8*C1/(pi*(x^2+y^2)^3) * (-x^3 + 3*x*y^2);
        field(1).gradUExact{1,2} = @(x,y) 8*C1/(pi*(x^2+y^2)^3) * (3*x^2*y - y^3);
        
        field(1).uExact{2} = @(x,y) C1*8*x*y/(pi*(x^2+y^2)^2);
        field(1).gradUExact{2,1} = @(x,y) 8*C1/(pi*(x^2+y^2)^3) * (3*x^2*y - y^3);
        field(1).gradUExact{2,2} = @(x,y) 8*C1/(pi*(x^2+y^2)^3) * (x^3 - 3*x*y^2);
    else
        field(1).uExact{1} = @(x,y) C2*(4*(x^2-y^2)/(pi*(x^2+y^2)^2) + 1);
        field(1).gradUExact{1,1} = @(x,y) 8*C2/(pi*(x^2+y^2)^3) * (-x^3 + 3*x*y^2);
        field(1).gradUExact{1,2} = @(x,y) 8*C2/(pi*(x^2+y^2)^3) * (3*x^2*y - y^3);
        
        field(1).uExact{2} = @(x,y) C2*8*x*y/(pi*(x^2+y^2)^2);
        field(1).gradUExact{2,1} = @(x,y) 8*C2/(pi*(x^2+y^2)^3) * (3*x^2*y - y^3);
        field(1).gradUExact{2,2} = @(x,y) 8*C2/(pi*(x^2+y^2)^3) * (x^3 - 3*x*y^2);
        
    end
    
    if FORMULATION == 4
        % q_f = sqrt(phi_0)*cos(theta)*(r - (2*B)/r)
        field(2).uExact{1} = @(x,y) (1-phi_0) * sqrt(phi_0) * x*(1 - 4/(pi*(x^2+y^2)) );
        field(2).gradUExact{1,1} = @(x,y) (1-phi_0) * sqrt(phi_0) *(1 - 4/(pi*(x^2+y^2)^2) * (x^2-y^2) );
        field(2).gradUExact{1,2} = @(x,y) (1-phi_0) * sqrt(phi_0) * -4/(pi*(x^2+y^2)^2) * 2*x*y;
    else
        % q_f = cos(theta)*(r - (2*B)/r)
        field(2).uExact{1} = @(x,y) (1-phi_0) * x*(1 - 4/(pi*(x^2+y^2)) );
        field(2).gradUExact{1,1} = @(x,y) (1-phi_0) * (1 - 4/(pi*(x^2+y^2)^2) * (x^2-y^2));
        field(2).gradUExact{1,2} = @(x,y) (1-phi_0) * -4/(pi*(x^2+y^2)) * 2*x*y;
    end
    
    % u_mx = -B*cos(theta)^2
    field(3).uExact{1} = @(x,y) 2*x^2/(pi*(x^2+y^2));
    field(3).gradUExact{1,1} = @(x,y) -4/((x^2+y^2)^2) * x*y^2;
    field(3).gradUExact{1,2} = @(x,y) 4/((x^2+y^2)^2) * -x^2*y;
    % u_my = B*(theta - sin(theta)*cos(theta))
    field(3).uExact{2} = @(x,y) (2/pi)*(atan(y/x) - x*y/(x^2+y^2));
    field(3).gradUExact{2,1} = @(x,y) -4/((x^2+y^2)^2) * y^3;
    field(3).gradUExact{2,2} = @(x,y) -4/((x^2+y^2)^2) * x*y^2;
    
    % q_m = cos(theta)*(r - (2*B)/r)
    field(4).uExact{1} = @(x,y) (1-phi_0) * x * (1 - 4/(pi*(x^2+y^2)) );
    field(4).gradUExact{1,1} = @(x,y) (1-phi_0) * (1 - 4/(pi*(x^2+y^2)) * (x^2-y^2) );
    field(4).gradUExact{1,2} = @(x,y) (1-phi_0) * -4/(pi*(x^2+y^2)) * 2*x*y;
    
    if FORMULATION == 2
        field(5).uExact{1} = @(x,y) C0*(1 - 4/(pi*(x^2+y^2)^2) * (x^2-y^2) );
        field(5).gradUExact{1,1} = @(x,y) 8*C0/(pi*(x^2+y^2)^3) * (-x^3 + 3*x*y^2);
        field(5).gradUExact{1,2} = @(x,y) 8*C0/(pi*(x^2+y^2)^3) * (3*x^2*y - y^3);
        
        field(5).uExact{2} = @(x,y) C0*(1 - 4/(pi*(x^2+y^2)^2) * 2*x*y );
        field(5).gradUExact{2,1} = @(x,y) 8*C0/(pi*(x^2+y^2)^3) * (3*x^2*y - y^3);
        field(5).gradUExact{2,2} = @(x,y) 8*C0/(pi*(x^2+y^2)^3) * (x^3 - 3*x*y^2);
    end
    
    Urx = field(1).uExact{1}; Ury = field(1).uExact{2};
    dxUrx = field(1).gradUExact{1,1}; dyUry = field(1).gradUExact{2,2};
    Qf = field(2).uExact{1};
    Umx = field(3).uExact{1}; Umy = field(3).uExact{2};
    dxUmx = field(3).gradUExact{1,1}; dyUmx = field(3).gradUExact{1,2}; dxUmy = field(3).gradUExact{2,1}; dyUmy = field(3).gradUExact{2,2};
    Qm = field(4).uExact{1};