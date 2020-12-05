function [E] = Engine_CFR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Engine Geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%
E.B = 0.08255;                      % Bore(m)
E.S = 0.1143;                       % Stroke(m)
E.l = 0.254;                        % Connecting Rod Length(m)
E.rc = 8.5;                         % Compression ratio

E.N = 2000;                         % Engine Speed (RPM)
E.teta0 = 15;                       % Begining of combustion (BTDC)
E.delteta = 50;                     % Combustion Duration (CA)
E.IVC = 00;                         % Inlet Valve Closing (ABDC)
E.EVO = 00;                         % Exhuast Valve Opening (BBDC)

E.a = E.S/2;                        % Radius of crank(m)
E.R = E.l/E.a;                      % COnnecting Rod/Crank Radius ratio
E.Vd = pi*(E.B^2)*E.S/4;            % Displacement Volume
E.Vc = E.Vd/(E.rc-1);               % Clearance volume(m^3)
E.Vtotal = E.Vc+E.Vd;               % Total Volume of Cylinder (m^3)
E.W = 2 * E.S * E.N /60;            % Piston mean speed (m/sec)
E.h = E.Vc/(pi*(E.B^2)/4);          % Distance between cylinder head and piston @TDC (m)

E.nw = 3;
E.aw = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%% Initial Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
E.Patm = 95;                        % Inlet Pressure (kPa)
E.Tatm = 294;                       % Inlet Temperature (K)
E.Tw = 373;                         % Cylinder-Wall Temperature (K)
E.phi = 0.80;                       % Equivalence ratio

tetatemp = pi*(-180+E.IVC)/180;
V_initial = ...                     % m3
    E.Vc*(1+(0.5*(E.rc-1)*(E.R+1-cos(tetatemp)-sqrt((E.R^2)-(sin(tetatemp)^2)))));
E.mass = ...                        % kg
    (E.Patm*V_initial)/(0.287*E.Tatm);  