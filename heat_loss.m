function [Qloss] = heat_loss(theta, theta_1, P, T, V)

N =2000;
gamma =1.35;
B = .08255;
S = .1143;
rc = 8.5;
Tw = 373;


Vd = S*pi*(B^2)/4;
Vc = Vd/(rc -1);

Pinit = 95 *1000;
Tinit = 294;

Vinit = Vc + Vd;
theta0 = 15;


Sp__ = 2*S*(N/60);
Pm = (Pinit * Vinit^gamma)/(V^gamma);

if theta >= -180 && theta < -theta0
    C1 = 2.28;
    C2 = 0;
elseif theta >= -theta0
    C1 = 2.28;
    C2 = 3.24e-3;
end

W_ = C1 * Sp__ + C2 * Vd * Tinit*(P - Pm)/Pinit/Vinit;

m = .8;
C = 3.26e-3;

h = C * (B^(m-1)) * (P^m) * (T^(.75-1.62*m)) * W_^m;

Ap = (pi*B^2)/4;
Ach = Ap;
Aw = Ap + Ach +(4*V/B);

dqdt = h * Aw *(T- Tw);

Qloss = dqdt * (theta - theta_1)/6/N;
end



