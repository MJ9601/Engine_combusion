clear
clc

B = .08255;
S = .1143;
l =.254;
rc = 8.5;

N = 2000;
a = S/2;
R = l/a;
Vd = 2*a*pi*(B^2)/4;
Vc = Vd/(rc -1);
nw = 3;
aw =5;

% countncountt value
Pinit = 95 *1000;
Tinit = 294;
phi = .8;
Vinit = Vc + Vd;

% constant value
mass_to = Pinit * Vinit /.287/Tinit;
cte = Pinit * Vinit^1.35;
Rc = 8.314;
n = Pinit * Vinit/Rc/ Tinit;
LHV = 43.448;
AF_st_mass = 14.7;

fuel_mass = phi * mass_to / AF_st_mass;  % kg

% Qin = LHV * 1000 * fuel_mass; 
gamma = 1.35;


cv = 0.71;

w1 = 0;
w2 = 0;

% matrcountx
V = zeros(360,1);
P = zeros(360,1);
T = zeros(360,1);


% countncounttcountalcountzatcounton
P(1) = Pinit;
V(1) = Vinit;
T(1) = Tinit;
theta0 = 15;
deltheta = 50;

theta = -180:1:180;
mb_1 =0;
for count=2: size(theta,2)
    V(count) = Vc*(1 + 0.5 *(rc-1)*(R + 1 - cosd(theta(count)) - sqrt(R^2 - (sind(theta(count)))^2)));
    P(count) = P(count-1) * ( (V(count-1)/V(count))^gamma);
    T(count) = T(count-1) * ( (V(count-1)/V(count))^ (gamma-1));


     if theta(count)>= -theta0 && theta(count)<= -theta0 + deltheta
         T_ = T(count);
         
         if theta0 ==0 && deltheta ==0
             mb =fuel_mass;
         else
            mb = (1-exp(-aw*((theta(count)-(-theta0))/deltheta)^(nw+1)))*fuel_mass;
         end
            Qin = LHV * 1000 * (mb - mb_1);
            mb_1 = mb;
        
        T(count) = T(count) + (Qin/mass_to/cv);
        P(count) = P(count) * T(count)/T_;
     end


         
        
        
        
end


V(:) = V(:)*1000;
P(:) = P(:)/101000;
thetax = -180:1:180;


% figure (1)
% plot(thetax, V(1:end))
% figure (2)
% plot(thetax, P(1:end))
% figure (3)
% plot(thetax , T(1:end))
% figure (4)
% plot(V, P)

figure (1)
subplot(3, 1,1);
plot(thetax, V(1:end),'-r')
ylabel('Volume (lit)')
hold on


subplot(3, 1,2);
plot(thetax, P(1:end), '-r')
ylabel('Perssure (atm)')
hold on



subplot(3, 1,3);
plot(thetax , T(1:end),'-r')
xlabel('-180 \leq \theta (degree)\leq 180')
ylabel('Temperature (k)')
hold on


figure (2)
plot(V, P, '-r')
xlabel('Volume (lit)')
ylabel('Perssure (atm)')
hold on

max(T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% considering heat losses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear


B = .08255;
S = .1143;
l =.254;
rc = 8.5;

N = 2000;
a = S/2;
R = l/a;
Vd = 2*a*pi*(B^2)/4;
Vc = Vd/(rc -1);
nw = 3;
aw =5;

% countncountt value
Pinit = 95 *1000;
Tinit = 294;
phi = .8;
Vinit = Vc + Vd;

% constant value
mass_to = Pinit * Vinit /.287/Tinit;
cte = Pinit * Vinit^1.35;
Rc = 8.314;
n = Pinit * Vinit/Rc/ Tinit;
LHV = 43.448;
AF_st_mass = 14.7;
theta0 = 15;
deltheta = 50;

fuel_mass = phi * mass_to / AF_st_mass;  % kg

% Qin = LHV * 1000 * fuel_mass; 
gamma = 1.35;


cv = 0.71;

w1 = 0;
w2 = 0;

% matrcountx
V = zeros(360,1);
P = zeros(360,1);
T = zeros(360,1);
count = 2;

% countncounttcountalcountzatcounton
P(1) = Pinit;
V(1) = Vinit;
T(1) = Tinit;


theta = -180:1:180;
mb_1 =0;
Qloss__ = zeros(1, 360);

for count=2: size(theta,2)
    V(count) = Vc*(1 + 0.5 *(rc-1)*(R + 1 - cosd(theta(count)) - sqrt(R^2 - (sind(theta(count)))^2)));
    P(count) = P(count-1) * ( (V(count-1)/V(count))^gamma);
    T(count) = T(count-1) * ( (V(count-1)/V(count))^ (gamma-1));


     if theta(count)>= -theta0 && theta(count)<= -theta0 + deltheta
         T_ = T(count);
         
         if theta0 ==0 && deltheta ==0
             mb =fuel_mass;
         else
            mb = (1-exp(-aw*((theta(count)-(-theta0))/deltheta)^(nw+1)))*fuel_mass;
         end
            Qin = LHV * 1000 * (mb - mb_1);
            mb_1 = mb;
        
        T(count) = T(count) + (Qin/mass_to/cv);
        P(count) = P(count) * T(count)/T_;
     end
    Qloss = heat_loss(theta(count), theta(count-1), P(count), T(count), V(count));
    Qloss__(count) = Qloss;
    T__ = T(count);
    T(count) = T(count) + (- Qloss/mass_to/cv);
    P(count) = P(count) * T(count)/T__;


         
        
        
        
end


V(:) = V(:)*1000;
P(:) = P(:)/101000;
thetax = -180:1:180;


% figure (1)
% plot(thetax, V(1:end))
% figure (2)
% plot(thetax, P(1:end))
% figure (3)
% plot(thetax , T(1:end))
% figure (4)
% plot(V, P)

figure (1)
subplot(3, 1,1);
plot(thetax, V(1:end),'--k')
ylabel('Volume (lit)')
legend('Without heat loss', 'With heat loss', 'Location','southwest')
legend boxoff

subplot(3, 1,2);
plot(thetax, P(1:end), '--k')
ylabel('Perssure (atm)')
legend('Without heat loss', 'With heat loss', 'Location','northwest')
legend boxoff

subplot(3, 1,3);
plot(thetax , T(1:end),'--k')
xlabel('-180 \leq \theta (degree)\leq 180')
ylabel('Temperature (k)')
legend('Without heat loss', 'With heat loss', 'Location','northwest')
legend boxoff

figure (2)
plot(V, P, '--k')
xlabel('Volume (lit)')
ylabel('Perssure (atm)')
legend('Without heat loss', 'With heat loss', 'Location','northeast')
legend boxoff

figure (3)
plot(theta, Qloss__)
xlabel('-180 \leq \theta (degree)\leq 180')
ylabel('Heat loss')
max(T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%