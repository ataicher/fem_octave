
syms theta r psi_m v_m psi_r v_r q q_f R f w_0 U_0 A B

% phi_0 = .024;
% w_0 = phi_0^2*4.9e-6;
% U_0 = 3.18e-9;
% c = w_0/U_0;

% ==alpha = 0;
% alpha = alpha*pi/180;
% A = 2*sin(alpha)^2/(pi-2*alpha-sin(2*alpha));
% B = 2/(pi-2*alpha-sin(2*alpha));

f = -B*theta*cos(theta)
fp = diff(f,theta)
fp2 = diff(f,theta,2)
fp3 = diff(f,theta,3)

psi_m = r*f;

% vr_m = simplify(1/r*diff(psi,theta));
% vtheta_m = simplify(-diff(psi,r));
v_m = simplify([1/r*diff(psi_m,theta) ; -diff(psi_m,r)])

psi_r = -w_0/U_0*(2*B/r + r)*sin(theta);
% psi_r = -w_0/U_0*(2*B/r)*sin(theta);

v_r = simplify([(1/r)*diff(psi_r,theta) ; -diff(psi_r,r)])

R = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];

v_m = simplify(R*v_m)
v_r = simplify(R*v_r)

q_f = (-2*B/r +r)*cos(theta)

q_m = q_f

syms x y gradV_m sigma

v_m = B*[x^2/(x^2+y^2) ; atan(y/x)-x*y/(x^2+y^2)]
q_f = -2*B*x/(x^2+y^2) + x;

gradV_m = simplify([diff(v_m,x) diff(v_m,y)])

q_fI = eye(2)*q_f

sigma = simplify([diff(v_m,x) diff(v_m,y)] - eye(2)*q_f)