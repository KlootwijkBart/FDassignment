% Load data 
load FTISxprt-20200309_flight3.mat % data from test flight


%% Flight data processing
% =========================================================================
%   FLIGHT DATA PROCESSING
% =========================================================================

% Flight Data
FD_t   = flightdata.time.data;                 % time [s]
FD_h   = flightdata.Dadc1_alt.data*0.3048;     % altitude [m]
FD_cas = flightdata.Dadc1_cas.data*0.51444444; % calibrated air speed [m/s]
FD_tas = flightdata.Dadc1_tas.data*0.51444444; % true air speed [m/s]  
FD_TAT = flightdata.Dadc1_tat.data+273.15;     % total air temperature [K]

% Angle of attack
FD_aoa = flightdata.vane_AOA.data;             % angle of attack [deg]

% Euler angles
FD_th  = flightdata.Ahrs1_Pitch.data;          % pitch angle [deg]
FD_phi = flightdata.Ahrs1_Roll.data;           % roll angle [deg]  
   
% Rotation rates
FD_q   = flightdata.Ahrs1_bPitchRate.data;     % pitch rate [deg/s]
FD_p   = flightdata.Ahrs1_bRollRate.data;      % roll rate [deg/s]
FD_r   = flightdata.Ahrs1_bYawRate.data;       % yaw rate [deg/s]

% Control surface deflection
FD_de  = flightdata.delta_e.data;              % elevator deflection [deg]
FD_da  = flightdata.delta_a.data;              % aileron deflection [deg]
FD_dr  = flightdata.delta_r.data;              % rudder deflection [deg]

% Stick input
FD_fe  = flightdata.column_fe.data;            % stick force [N]
FD_Se  = flightdata.column_Se.data;            % stick deflection [deg]

% Fuel used
FD_FUl = flightdata.lh_engine_FU.data*0.45359237;  % fuel used by left engine [kg]
FD_FUr = flightdata.rh_engine_FU.data*0.45359237;  % fuel used by right engine [kg]


%% Compute parameters

% Constant values concerning atmosphere and gravity

rho0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)
gamma  = 1.4;             % ratio of specific heats for air [-]
P0     = rho0*R*Temp0;    % ISA sea-level pressure [Pa]

% Prepare tables
FD_p=[]; FD_rho=[]; FD_eas=[]; FD_mass=[]; 

% Run
for i = 1:length(FD_t)
     
% compute ISA values    
p   = P0*(1+(lambda*FD_h(i)/Temp0))^-(g/(lambda*R)); % static pressure [Pa]
T   = Temp0 + (lambda*FD_h(i));           % equivalent air temperature [K]
rho = p/(R*T);                            % air density [kg/m^3]

% compute equivalent airspeed from true airspeed and densities
eas = FD_tas(i)*(rho/rho0)^0.5;           % [m/s]

% compute instantaneous mass 
m = m0 - FD_FUr(i) - FD_FUl(i);           % ramp mass minus fuel used [kg]

% append values to tables
FD_p    = [FD_p;p];
FD_rho  = [FD_rho;rho];
FD_eas  = [FD_eas;eas];
FD_mass = [FD_mass;m];

end

%% Stationary flight condition

% indices for start and end of different eigenmodes, in the order shown
% below
idxstart = [30311,32380,33441,34661,35131,37261];
idxend   = [32011,32410,33621,35011,35511,39411];

% ask user input 
eigenmode = input(['',...
                   '\n1: Phugoid',...
                   '\n2: Short period',...
                   '\n3: A-periodic roll',...
                   '\n4: Dutch roll',...
                   '\n5: Dutch roll damped',...
                   '\n6: Spiral',...
                   '\n',...
                   '\nWhich eigenmode is to be simulated? ']);

idxstart = idxstart(eigenmode);
idxend   = idxend(eigenmode);

% initial stationary flight conditions @ start of eigenmotion
hp0    = FD_h(idxstart);    % pressure altitude in the stationary flight condition [m]
V0     = FD_eas(idxstart);  % equivalent airspeed in the stationary flight condition [m/sec]
alpha0 = FD_aoa(idxstart)*pi/180;  % angle of attack in the stationary flight condition [rad]
th0    = FD_th(idxstart)*pi/180;   % pitch angle in the stationary flight condition [rad]


%% Cit_par

% Citation 550 - Linear simulation

% Initial aircraft mass
m0     = 60500/9.81;     % ramp mass [kg]

% standard engine fuel flow for BOTH engines
mdot   = 2*0.048;        % [kg/sec] 

% aerodynamic properties
e      = 0.852874;       % Oswald factor [ ]
CD0    = 0.031389;       % Zero lift drag coefficient [ ]
CLa    = 4.97;           % Slope of CL-alpha curve [ ]

% Longitudinal stability
Cma    = -0.5626;        % longitudinal stabilty [ ]
Cmde   = -1.1642;        % elevator effectiveness [ ]

% Aircraft geometry

S      = 30.00;	          % wing area [m^2]
Sh     = 0.2*S;           % stabiliser area [m^2]
Sh_S   = Sh/S;	          % [ ]
lh     = 0.71*5.968;      % tail length [m]
c      = 2.0569;	      % mean aerodynamic cord [m]
lh_c   = lh/c;	          % [ ]
b      = 15.911;	      % wing span [m]
bh     = 5.791;	          % stabilser span [m]
A      = b^2/S;           % wing aspect ratio [ ]
Ah     = bh^2/Sh;         % stabilser aspect ratio [ ]
Vh_V   = 1;		          % [ ]
ih     = -2*pi/180;       % stabiliser angle of incidence [rad]

% Constant values concerning atmosphere and gravity

rho0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)
gamma  = 1.4;             % ratio of specific heats for air [-]
P0     = 101325;          % ISA sea-level pressure [Pa]      

rho    = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
W      = m*g;				                                    % [N]       (aircraft weight)

% Constant values concerning aircraft inertia

muc    = m/(rho*S*c);
mub    = m/(rho*S*b);
KX2    = 0.019;
KZ2    = 0.042;
KXZ    = 0.002;
KY2    = 1.25*1.114;

% Aerodynamic constants

Cmac   = 0;                     % Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa;   		        % Wing normal force slope [ ]
CNha   = 2*pi*Ah/(Ah+2);        % Stabiliser normal force slope [ ]
depsda = 4/(A+2);               % Downwash gradient [ ]

% Lift and drag coefficient

CL = 2*W/(rho*V0^2*S);               % Lift coefficient [ ]
CD = CD0 + (CLa*alpha0)^2/(pi*A*e);  % Drag coefficient [ ]

% Stabiblity derivatives

CX0    = W*sin(th0)/(0.5*rho*V0^2*S);
CXu    = -0.095;
CXa    = -0.47966;
CXadot = +0.08330;
CXq    = -0.28170;
CXde   = -0.03728;

CZ0    = -W*cos(th0)/(0.5*rho*V0^2*S);
CZu    = -0.37616;
CZa    = -5.74340;
CZadot = -0.00350;
CZq    = -5.66290;
CZde   = -0.69612;

Cmu    = +0.06990;
Cmadot = +0.17800;
Cmq    = -8.79415;

CYb    = -0.7500;
CYbdot =  0     ;
CYp    = -0.0304;
CYr    = +0.8495;
CYda   = -0.0400;
CYdr   = +0.2300;

Clb    = -0.10260;
Clp    = -0.71085;
Clr    = +0.23760;
Clda   = -0.23088;
Cldr   = +0.03440;

Cnb    =  +0.1348;
Cnbdot =   0     ;
Cnp    =  -0.0602;
Cnr    =  -0.2061;
Cnda   =  -0.0120;
Cndr   =  -0.0939;

clf()
plot(FD_t,FD_h)


%{ 

Index range of demonstrated eigenmotions: 

Phugoid            30311:32011
Short period       32380:32410
A-periodic roll    33441:33621
Dutch roll         34661:35011
Dutch roll damped  35131:35511
Spiral             37261:39411  

%}


%% Instant
% V0 = Veq_tab(32380);


%% Define State Space Model

% =========================================================================
%   SYMMETRIC EQUATIONS
% =========================================================================

C_1s = [ -2*muc*c/(V0*V0)         0            0             0             ;...
                0        (CZa-2*muc)*(c/V0)    0             0             ;...
                0                 0          -c/V0           0             ;...
                0            Cmadot*c/V0       0   -2*muc*KY2*c*c/(V0*V0)  ];

C_2s = [ CXu/V0   CXa   CZ0       CXa*c/V0      ;...
         CZu/V0   CZa  -CX0  (c/V0)*(CZq+2*muc) ;...
            0      0     0          c/V0        ;...
         Cmu/V0   Cma    0        Cmq*c/V0      ];

C_3s = [ CXde  ;...
         CZde  ;...
           0   ;...
         Cmde  ] ;

% Compute A and B matrices and define state space model    
% State vector: [u, AoA, theta, q]
% Input vector: [delta_e]
As = -1*inv(C_1s)*C_2s;
Bs = -1*inv(C_1s)*C_3s;
Cs = eye(size(As));
Ds = [0];
syss = ss(As,Bs,Cs,Ds);


% =========================================================================
%   ASYMMETRIC EQUATIONS
% =========================================================================

C_1a = [ (CYbdot-2*mub)*b/V0     0              0                     0            ;...
                0            -0.5*b/V0          0                     0            ;...
                0                0      -2*mub*KX2*b^2/V0^2    2*mub*KXZ*b^2/V0^2  ;...
           Cnbdot*b/V0           0       2*mub*KXZ*b^2/V0^2   -2*mub*KZ2*b^2/V0^2  ];

C_2a = [ CYb     CL   CYp*b/(2*V0)  (CYr-4*mub)*b/(2*V0) ;...
          0       0     b/(2*V0)              0          ;...
         Clb      0   Clp*b/(2*V0)       Clr*b/(2*V0)    ;...
         Cnb      0   Cnp*b/(2*V0)       Cnr*b/(2*V0)    ];

C_3a = [ CYda  CYdr  ;...
          0     0    ;...
         Clda  Cldr  ;...
         Cnda  Cndr  ];

% Compute A and B matrices and define state space model    
% State vector: [beta, phi, roll rate, yaw rate]
% Input vector: [delta_a, delta_r]
Aa = -1*inv(C_1a)*C_2a;
Ba = -1*inv(C_1a)*C_3a;
Ca = eye(size(Aa));
Da = [0];
sysa = ss(Aa,Ba,Ca,Da);


%% Plot Results
% =========================================================================
% PLOT
% =========================================================================

%{ 

Index range of demonstrated eigenmotions:

Phugoid            30311:32011
Short period       32380:32410
A-periodic roll    33441:33621
Dutch roll         34661:35011
Dutch roll damped  35131:35511
Spiral             37261:39411  

%}

u_de = FD_de(idxstart:idxend) ; % negative impulse input in elevator deflection
Veq  = FD_eas(idxstart:idxend);
n = length(u_de);
t = linspace(0,(0.1*n),n);
[y,t] = lsim(syss,u_de,t);
clf();
tiledlayout(4,1)

nexttile
plot(t,V0+y(:,1,1)/V0,t,Veq)  % plot speed V 
title('V [m/s]')
legend()

nexttile
plot(t,alpha0+y(:,2,1)) % plot AoA alpha
title('Angle of attack [deg]')

nexttile
plot(t,th0+y(:,3,1))    % plot pitch angle theta
title('Pitch angle [deg]')

nexttile
plot(t,FD_q(1)+y(:,4,1))        % plot pitch rate q
title('Pitch rate [deg/s]')


