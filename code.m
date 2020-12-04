clc 
clear


B = 0.08255;
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
Pinit = 95 * 1000;
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

Qin = LHV * 1000 * fuel_mass; 


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



for theta = -180:1:180
    V(count) = Vc*(1 + 0.5 *(rc-1)*(R + 1 - cosd(theta) - sqrt(R^2 - (sind(theta))^2)));
    
    if theta >= -180 && theta <0
        P(count) = cte/V(count)^1.35;
        T(count) = P(count)* V(count)/Rc/n;
        w1 = w1 +((P(count) + P(count-1))*(V(count) - V(count-1))/2);
%         disp(count)
    elseif theta == 0
        T(count) = (Qin/mass_to/cv) + T(180);
        P(count) = n*Rc*T(count)/V(count);
%         disp(count)
    else
        cte = P(182)*(V(182))^1.35;
        P(count) = cte/V(count)^1.35;
        T(count) = P(count)* V(count)/Rc/n;
        w2 = w2 +((P(count) + P(count-1))*(V(count) - V(count-1))/2);
        
    end
    count = count +1;
end
V(:) = V(:)*1000;
P(:) = P(:)/101000;
thetax = -180:1:180;

% 
% figure (1)
% plot(thetax, V(2:end))
% figure (2)
% plot(thetax, P(2:end))
% figure (3)
% plot(thetax , T(2:end))
% figure (4)
% plot(V, P)
% disp((w2 -w1))
% disp('Output power:')
% disp((w2 -w1)/735.499)
% disp('horsepower')

figure (1)
subplot(3, 1,1);
plot(thetax, V(2:end), '-r')
ylabel('Volume (lit)')
hold on

subplot(3, 1,2);
plot(thetax, P(2:end), '-r')
ylabel('Perssure (atm)')
hold on

subplot(3, 1,3);
plot(thetax , T(2:end), '-r')
xlabel('-180 \leq \theta (degree)\leq 180')
ylabel('Temperature (k)')
hold on


figure (2)
plot(V, P, '-r')
xlabel('Volume (lit)')
ylabel('Perssure (atm)')
hold on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% breakcountng code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


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
plot(thetax, V(1:end),'--k')
ylabel('Volume (lit)')
legend('Otto cycle ', 'Wiebe Function', 'Location','southwest')
legend boxoff

subplot(3, 1,2);
plot(thetax, P(1:end), '--k')
ylabel('Perssure (atm)')
legend('Otto cycle ', 'Wiebe Function', 'Location','northwest')
legend boxoff

subplot(3, 1,3);
plot(thetax , T(1:end),'--k')
xlabel('-180 \leq \theta (degree)\leq 180')
ylabel('Temperature (k)')
legend('Otto cycle ', 'Wiebe Function', 'Location','northwest')
legend boxoff

figure (2)
plot(V, P, '--k')
xlabel('Volume (lit)')
ylabel('Perssure (atm)')
legend('Otto cycle ', 'Wiebe Function')
legend boxoff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%