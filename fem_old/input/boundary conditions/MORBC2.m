% MORBC.m

% boundary condition options
%   1. Exact solution with corner removed
%   2. Exact solution
%   3. zero stress far field condtion with corner removed
%   4. zero stress far field condition
%   5. zero stress far field condition with no fluid flow through left
%   boundary except at the bottom [0 y_f]

switch SETBC
    
    %EXACT SOLUTION WITH CORNER REMOVED
    case 1
        % v_r
        % right edge
        %   Dirichlet - exact solution
        field(1).bndry(1).loc = @(x,y) (x-LX);
        field(1).bndry(1).alpha = [0 0];
        field(1).bndry(1).beta = [1 1];
        field(1).bndry(1).eta = {@(x,y)  Urx(x,y), @(x,y) Ury(x,y)};
        % top edge
        %   Dirichlet - exact solution
        field(1).bndry(2).loc = @(x,y) (y-LY);
        field(1).bndry(2).alpha = [0 0];
        field(1).bndry(2).beta = [1 1];
        field(1).bndry(2).eta = {@(x,y) Ury(x,y) , @(x,y) -Urx(x,y)};
        % bottom edge
        %   Dirichlet - symmetry
        field(1).bndry(3).loc = @(x,y) y;
        field(1).bndry(3).alpha = [0 0];
        field(1).bndry(3).beta = [1 1];
        field(1).bndry(3).eta = {@(x,y) -Ury(x,y) , @(x,y) Urx(x,y)};
        % left edge
        %   Dirichlet - exact solution
        field(1).bndry(4).loc = @(x,y) x;
        field(1).bndry(4).alpha = [0 0];
        field(1).bndry(4).beta = [1 1];
        field(1).bndry(4).eta = {@(x,y) -Urx(x,y) , @(x,y) -Ury(x,y)};
        % corner edges
        %   Dirichlet - exact solution
        field(1).bndry(5).loc = @(x,y) (x-x_0);
        field(1).bndry(5).alpha = [0 0];
        field(1).bndry(5).beta = [1 1];
        field(1).bndry(5).eta = {@(x,y) Urx(x,y), @(x,y) Ury(x,y)};
        field(1).bndry(6).loc = @(x,y) (y-y_0);
        field(1).bndry(6).alpha = [0 0];
        field(1).bndry(6).beta = [1 1];
        field(1).bndry(6).eta = {@(x,y) Ury(x,y), @(x,y) -Urx(x,y)};
        
        % v_m
        % right edge
        %   Dirichlet - exact solution
        field(3).bndry(1).loc = @(x,y) (x-LX);
        field(3).bndry(1).alpha = [0 0];
        field(3).bndry(1).beta = [1 1];
        field(3).bndry(1).eta = {@(x,y) Umx(x,y) , @(x,y) Umy(x,y)};
        % top edge
        %   Dirichlet - exact solution
        field(3).bndry(2).loc = @(x,y) (y-LY);
        field(3).bndry(2).alpha = [0 0];
        field(3).bndry(2).beta = [1 1];
        field(3).bndry(2).eta = {@(x,y) Umy(x,y) , @(x,y) -Umx(x,y)};
        % bottom edge
        %   symmetry
        field(3).bndry(3).loc = @(x,y) y;
        field(3).bndry(3).alpha = [0 0];
        field(3).bndry(3).beta = [1 1];
        field(3).bndry(3).eta = {@(x,y) -Umy(x,y) , @(x,y) Umx(x,y)};
        % left edge
        %   Dirichlet
        field(3).bndry(4).loc = @(x,y) x;
        field(3).bndry(4).alpha = [0 0];
        field(3).bndry(4).beta = [1 1];
        field(3).bndry(4).eta = {@(x,y) 0, @(x,y) -1};
        % corner edges
        %   Dirichlet - exact solution
        field(3).bndry(5).loc = @(x,y) (x-x_0);
        field(3).bndry(5).alpha = [0 0];
        field(3).bndry(5).beta = [1 1];
        field(3).bndry(5).eta = {@(x,y) -Umx(x,y), @(x,y) -Umy(x,y)};
        field(3).bndry(6).loc = @(x,y) (y-y_0);
        field(3).bndry(6).alpha = [0 0];
        field(3).bndry(6).beta = [1 1];
        field(3).bndry(6).eta = {@(x,y) -Umy(x,y), @(x,y) Umx(x,y)};
        
        % EXACT SOLUTION WITH FLUID PRESSURE SET TO INTERPOLANT AT CORNER
    case 2
        % v_r
        % right edge
        %   Dirichlet - exact solution
        field(1).bndry(1).loc = @(x,y) (x-LX);
        field(1).bndry(1).alpha = [0 1];
        field(1).bndry(1).beta = [1 0];
        field(1).bndry(1).eta = {@(x,y)  Urx(x,y), @(x,y) 0};
        % top edge
        %   Dirichlet - exact solution
        field(1).bndry(2).loc = @(x,y) (y-LY);
        field(1).bndry(2).alpha = [0 1];
        field(1).bndry(2).beta = [1 0];
        field(1).bndry(2).eta = {@(x,y) Ury(x,y) , @(x,y) 0};
        % bottom edge
        %   Dirichlet - symmetry
        field(1).bndry(3).loc = @(x,y) y;
        field(1).bndry(3).alpha = [0 1];
        field(1).bndry(3).beta = [1 0];
        field(1).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
        % left edge
        %   Dirichlet - exact solution
        field(1).bndry(4).loc = @(x,y) x;
        field(1).bndry(4).alpha = [0 1];
        field(1).bndry(4).beta = [1 0];
        field(1).bndry(4).eta = {@(x,y) -Urx(x,y) , @(x,y) 0};
        
        % v_m
        % right edge
        %   Dirichlet - exact solution
        field(3).bndry(1).loc = @(x,y) (x-LX);
        field(3).bndry(1).alpha = [0 0];
        field(3).bndry(1).beta = [1 1];
        field(3).bndry(1).eta = {@(x,y) Umx(x,y) , @(x,y) Umy(x,y)};
        % top edge
        %   Dirichlet - exact solution
        field(3).bndry(2).loc = @(x,y) (y-LY);
        field(3).bndry(2).alpha = [0 0];
        field(3).bndry(2).beta = [1 1];
        field(3).bndry(2).eta = {@(x,y) Umy(x,y) , @(x,y) -Umx(x,y)};
        % bottom edge
        %   symmetry
        field(3).bndry(3).loc = @(x,y) y;
        field(3).bndry(3).alpha = [0 1];
        field(3).bndry(3).beta = [1 0];
        field(3).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
        % left edge
        %   Dirichlet
        field(3).bndry(4).loc = @(x,y) x;
        field(3).bndry(4).alpha = [0 0];
        field(3).bndry(4).beta = [1 1];
        field(3).bndry(4).eta = {@(x,y) 0, @(x,y) -1};
        
        % ZERO STRESS FAR FIELD CONDTION WITH CORNER REMOVED
    case 3
        % v_r
        % right edge
        %   Neumann - far field
        field(1).bndry(1).loc = @(x,y) (x-LX);
        field(1).bndry(1).alpha = [1 1];
        field(1).bndry(1).beta = [0 0];
        field(1).bndry(1).eta = {@(x,y) phi(x,y)^(1+Theta)*LX , @(x,y) 0};
        % top edge
        %   Dirichlet - far field
        field(1).bndry(2).loc = @(x,y) (y-LY);
        field(1).bndry(2).alpha = [0 1];
        field(1).bndry(2).beta = [1 0];
        field(1).bndry(2).eta = {@(x,y) 0 , @(x,y) 0};
        % bottom edge
        %   Dirichlet - symmetry
        field(1).bndry(3).loc = @(x,y) y;
        field(1).bndry(3).alpha = [0 1];
        field(1).bndry(3).beta = [1 0];
        field(1).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
        % left edge
        %   Neumann - zero pressure
        field(1).bndry(4).loc = @(x,y) x;
        field(1).bndry(4).alpha = [1 1];
        field(1).bndry(4).beta = [0 0];
        field(1).bndry(4).eta = {@(x,y) 0 , @(x,y) 0};
        % corner edges
        %   Dirichlet - exact solution
        field(1).bndry(5).loc = @(x,y) (x-x_0);
        field(1).bndry(5).alpha = [0 1];
        field(1).bndry(5).beta = [1 0];
        field(1).bndry(5).eta = {@(x,y) -Urx(x,y), @(x,y) 0};
        field(1).bndry(6).loc = @(x,y) (y-y_0);
        field(1).bndry(6).alpha = [0 1];
        field(1).bndry(6).beta = [1 0];
        field(1).bndry(6).eta = {@(x,y) -Ury(x,y), @(x,y) 0};
        
        % v_m
        % right edge
        %   Neumann - zero normal stress
        field(3).bndry(1).loc = @(x,y) (x-LX);
        field(3).bndry(1).alpha = [1 1];
        field(3).bndry(1).beta = [0 0];
        field(3).bndry(1).eta = {@(x,y) x , @(x,y) 0};
        % top edge
        %   Neumann - zero normal stress
        field(3).bndry(2).loc = @(x,y) (y-LY);
        field(3).bndry(2).alpha = [1 1];
        field(3).bndry(2).beta = [0 0];
        field(3).bndry(2).eta = {@(x,y)  x, @(x,y) 0};
        % bottom edge
        %   symmetry
        field(3).bndry(3).loc = @(x,y) y;
        field(3).bndry(3).alpha = [0 1];
        field(3).bndry(3).beta = [1 0];
        field(3).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
        % left edge
        %   Dirichlet
        field(3).bndry(4).loc = @(x,y) x;
        field(3).bndry(4).alpha = [0 0];
        field(3).bndry(4).beta = [1 1];
        field(3).bndry(4).eta = {@(x,y) 0, @(x,y) -1};
        % corner edges
        %   Dirichlet - exact solution
        field(3).bndry(5).loc = @(x,y) (x-x_0);
        field(3).bndry(5).alpha = [0 0];
        field(3).bndry(5).beta = [1 1];
        field(3).bndry(5).eta = {@(x,y) -Umx(x,y), @(x,y) -Umy(x,y)};
        field(3).bndry(6).loc = @(x,y) (y-y_0);
        field(3).bndry(6).alpha = [0 0];
        field(3).bndry(6).beta = [1 1];
        field(3).bndry(6).eta = {@(x,y) -Umy(x,y), @(x,y) Umx(x,y)};
        
        % ZERO STRESS FAR FIELD CONDITION
    case 4
        % v_r
        % right edge
        %   Neumann - far field
        field(1).bndry(1).loc = @(x,y) (x-LX);
        field(1).bndry(1).alpha = [1 1];
        field(1).bndry(1).beta = [0 0];
        field(1).bndry(1).eta = {@(x,y) phi(x,y)^(1+Theta)*LX , @(x,y) 0};
        % top edge
        %   Dirichlet - far field
        field(1).bndry(2).loc = @(x,y) (y-LY);
        field(1).bndry(2).alpha = [0 1];
        field(1).bndry(2).beta = [1 0];
        field(1).bndry(2).eta = {@(x,y) 0 , @(x,y) 0};
        % bottom edge
        %   Dirichlet - symmetry
        field(1).bndry(3).loc = @(x,y) y;
        field(1).bndry(3).alpha = [0 1];
        field(1).bndry(3).beta = [1 0];
        field(1).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
        % left edge
        %   Neumann - zero pressure
        field(1).bndry(4).loc = @(x,y) x;
        field(1).bndry(4).alpha = [1 1];
        field(1).bndry(4).beta = [0 0];
        field(1).bndry(4).eta = {@(x,y) 0 , @(x,y) 0};
        
        % v_m
        % right edge
        %   Neumann - zero normal stress
        field(3).bndry(1).loc = @(x,y) (x-LX);
        field(3).bndry(1).alpha = [1 1];
        field(3).bndry(1).beta = [0 0];
        field(3).bndry(1).eta = {@(x,y) x , @(x,y) 0};
        % top edge
        %   Neumann - zero normal stress
        field(3).bndry(2).loc = @(x,y) (y-LY);
        field(3).bndry(2).alpha = [1 1];
        field(3).bndry(2).beta = [0 0];
        field(3).bndry(2).eta = {@(x,y)  x, @(x,y) 0};
        % bottom edge
        %   symmetry
        field(3).bndry(3).loc = @(x,y) y;
        field(3).bndry(3).alpha = [0 1];
        field(3).bndry(3).beta = [1 0];
        field(3).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
        % left edge
        %   Dirichlet
        field(3).bndry(4).loc = @(x,y) x;
        field(3).bndry(4).alpha = [0 0];
        field(3).bndry(4).beta = [1 1];
        field(3).bndry(4).eta = {@(x,y) 0, @(x,y) -x};
        
        % NO FLUID FLOW THROUGH LEFT BOUNDARY EXCEPT AT BOTTOM [0 y_f]. ZERO
        % STRESS FAR FIELD CONDITION
    case 5
        % v_r
        % right edge
        %   Neumann - far field
        field(1).bndry(1).loc = @(x,y) (x-LX);
        field(1).bndry(1).alpha = [1 1];
        field(1).bndry(1).beta = [0 0];
        field(1).bndry(1).eta = {@(x,y) phi(x,y)^(1+Theta)*LX , @(x,y) 0};
        % top edge
        %   Dirichlet - far field
        field(1).bndry(2).loc = @(x,y) (y-LY);
        field(1).bndry(2).alpha = [0 1];
        field(1).bndry(2).beta = [1 0];
        field(1).bndry(2).eta = {@(x,y) 0 , @(x,y) 0};
        % bottom edge
        %   Dirichlet - symmetry
        field(1).bndry(3).loc = @(x,y) y;
        field(1).bndry(3).alpha = [0 1];
        field(1).bndry(3).beta = [1 0];
        field(1).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
        % left edge
        %  Dirichlet - no flow
        field(1).bndry(4).loc = @(x,y) x + (y<y_f);
        field(1).bndry(4).alpha = [0 1];
        field(1).bndry(4).beta = [1 0];
        field(1).bndry(4).eta = {@(x,y) 0 , @(x,y) 0};
        %   Neumann - zero pressure
        field(1).bndry(5).loc = @(x,y) x + (y>y_f);
        field(1).bndry(5).alpha = [1 1];
        field(1).bndry(5).beta = [0 0];
        field(1).bndry(5).eta = {@(x,y) 0 , @(x,y) 0};
        
        % v_m
        % right edge
        %   Neumann - zero normal stress
        field(3).bndry(1).loc = @(x,y) (x-LX);
        field(3).bndry(1).alpha = [1 1];
        field(3).bndry(1).beta = [0 0];
        field(3).bndry(1).eta = {@(x,y) x , @(x,y) 0};
        % top edge
        %   Neumann - zero normal stress
        field(3).bndry(2).loc = @(x,y) (y-LY);
        field(3).bndry(2).alpha = [1 1];
        field(3).bndry(2).beta = [0 0];
        field(3).bndry(2).eta = {@(x,y)  x, @(x,y) 0};
        % bottom edge
        %   symmetry
        field(3).bndry(3).loc = @(x,y) y;
        field(3).bndry(3).alpha = [0 1];
        field(3).bndry(3).beta = [1 0];
        field(3).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
        % left edge
        %   Dirichlet
        field(3).bndry(4).loc = @(x,y) x;
        field(3).bndry(4).alpha = [0 0];
        field(3).bndry(4).beta = [1 1];
        field(3).bndry(4).eta = {@(x,y) 0, @(x,y) -1};
        % corner edges
        
    case 6
        % v_r
        % right edge
        %   Dirichlet - exact solution
        % field(1).bndry(1).loc = @(x,y) (x-Lx);
        % field(1).bndry(1).alpha = [0 1];
        % field(1).bndry(1).beta = [1 0];
        % field(1).bndry(1).eta = {@(x,y)  Urx(x,y), @(x,y) 0};
        %   Neumann - far field
        field(1).bndry(1).loc = @(x,y) (x-Lx);
        field(1).bndry(1).alpha = [1 1];
        field(1).bndry(1).beta = [0 0];
        field(1).bndry(1).eta = {@(x,y) phi(x,y)^(1+Theta)*Lx , @(x,y) 0};
        % top edge
        %   Dirichlet - exact solution
        % field(1).bndry(2).loc = @(x,y) (y-Ly);
        % field(1).bndry(2).alpha = [0 1];
        % field(1).bndry(2).beta = [1 0];
        % field(1).bndry(2).eta = {@(x,y) Ury(x,y) , @(x,y) 0};
        %   Dirichlet - far field
        field(1).bndry(2).loc = @(x,y) (y-Ly);
        field(1).bndry(2).alpha = [0 1];
        field(1).bndry(2).beta = [1 0];
        field(1).bndry(2).eta = {@(x,y) 0 , @(x,y) 0};
        % bottom edge
        %   Dirichlet - symmetry
        field(1).bndry(3).loc = @(x,y) y;
        field(1).bndry(3).alpha = [0 1];
        field(1).bndry(3).beta = [1 0];
        field(1).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
        % left edge
        %   Dirichlet - exact solution
        % field(1).bndry(4).loc = @(x,y) x;
        % field(1).bndry(4).alpha = [0 1];
        % field(1).bndry(4).beta = [1 0];
        % field(1).bndry(4).eta = {@(x,y) -Urx(x,y) , @(x,y) 0};
        %   Neumann - zero pressure
        field(1).bndry(4).loc = @(x,y) x;
        field(1).bndry(4).alpha = [1 1];
        field(1).bndry(4).beta = [0 0];
        field(1).bndry(4).eta = {@(x,y) 0 , @(x,y) 0};
        % corner edges
        %   Dirichlet - exact solution
        % field(1).bndry(5).loc = @(x,y) (x-x_0);
        % field(1).bndry(5).alpha = [0 1];
        % field(1).bndry(5).beta = [1 0];
        % field(1).bndry(5).eta = {@(x,y) -Urx(x,y), @(x,y) 0};
        % field(1).bndry(6).loc = @(x,y) (y-y_0);
        % field(1).bndry(6).alpha = [0 1];
        % field(1).bndry(6).beta = [1 0];
        % field(1).bndry(6).eta = {@(x,y) -Ury(x,y), @(x,y) 0};
        %   Dirichlet - average flux
        %   v_rx ~= -(w_0/U_0)*2*B*cos(2*theta(x,y))/r2(x,y);
        %   v_ry = -(w_0/U_0)*2*B*sin(2*theta(x,y))/r2(x,y);
        % field(1).bndry(5).loc = @(x,y) (x-x_0);
        % field(1).bndry(5).alpha = [0 1];
        % field(1).bndry(5).beta = [1 0];
        % field(1).bndry(5).eta = {@(x,y) (2*B*w_0/U_0)*1/(x_0^2+y_0^2), @(x,y) 0};
        % field(1).bndry(6).loc = @(x,y) (y-y_0);
        % field(1).bndry(6).alpha = [0 1];
        % field(1).bndry(6).beta = [1 0];
        % field(1).bndry(6).eta = {@(x,y) (2*B*w_0/U_0)*x_0/(y_0*(x_0^2+y_0^2)), @(x,y) 0};
        %
        % v_m
        % right edge
        %   Dirichlet - exact solution
        % field(3).bndry(1).loc = @(x,y) (x-Lx);
        % field(3).bndry(1).alpha = [0 0];
        % field(3).bndry(1).beta = [1 1];
        % field(3).bndry(1).eta = {@(x,y) Umx(x,y) , @(x,y) Umy(x,y)};
        %   Neumann - zero normal stress
        field(3).bndry(1).loc = @(x,y) (x-Lx);
        field(3).bndry(1).alpha = [1 1];
        field(3).bndry(1).beta = [0 0];
        field(3).bndry(1).eta = {@(x,y) x , @(x,y) 0};
        % top edge
        %   Dirichlet - exact solution
        % field(3).bndry(2).loc = @(x,y) (y-Ly);
        % field(3).bndry(2).alpha = [0 0];
        % field(3).bndry(2).beta = [1 1];
        % field(3).bndry(2).eta = {@(x,y) Umy(x,y) , @(x,y) -Umx(x,y)};
        %   Neumann - exact normal stress
        % field(3).bndry(2).loc = @(x,y) (y-Ly);
        % field(3).bndry(2).alpha = [1 1];
        % field(3).bndry(2).beta = [0 0];
        % field(3).bndry(2).eta = {@(x,y) x , @(x,y) 0};
        %   Neumann - zero normal stress
        field(3).bndry(2).loc = @(x,y) (y-Ly);
        field(3).bndry(2).alpha = [1 1];
        field(3).bndry(2).beta = [0 0];
        field(3).bndry(2).eta = {@(x,y)  x, @(x,y) 0};
        % bottom edge
        %   symmetry
        field(3).bndry(3).loc = @(x,y) y;
        field(3).bndry(3).alpha = [0 1];
        field(3).bndry(3).beta = [1 0];
        field(3).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
        % left edge
        %   Dirichlet
        field(3).bndry(4).loc = @(x,y) x;
        field(3).bndry(4).alpha = [0 0];
        field(3).bndry(4).beta = [1 1];
        field(3).bndry(4).eta = {@(x,y) 0, @(x,y) -1};
        % corner edges
        %   Dirichlet - exact solution
        % field(3).bndry(5).loc = @(x,y) (x-x_0);
        % field(3).bndry(5).alpha = [0 0];
        % field(3).bndry(5).beta = [1 1];
        % field(3).bndry(5).eta = {@(x,y) -Umx(x,y), @(x,y) -Umy(x,y)};
        % field(3).bndry(6).loc = @(x,y) (y-y_0);
        % field(3).bndry(6).alpha = [0 0];
        % field(3).bndry(6).beta = [1 1];
        % field(3).bndry(6).eta = {@(x,y) -Umy(x,y), @(x,y) Umx(x,y)};
end