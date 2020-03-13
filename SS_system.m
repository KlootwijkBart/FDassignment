% Load data 
load Cit_par_data.mat              % parameters provided
load FTISxprt-20200309_flight3.mat % data from test flight

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

% =========================================================================
% PLOT
% =========================================================================

u = zeros(2000,1); u(1) = -15.9317; % negative impulse input in elevator deflection
t = linspace(0,200,2000)';
[y,t] = lsim(syss,u,t);
clf();
tiledlayout(4,1)

nexttile
plot(t,V0+y(:,1,1)/V0)  % plot speed V 
title('V [m/s]')

nexttile
plot(t,alpha0+y(:,2,1)) % plot AoA alpha
title('Angle of attack [rad]')

nexttile
plot(t,th0+y(:,3,1))    % plot pitch angle theta
title('Pitch angle [rad')

nexttile
plot(t,y(:,4,1))        % plot pitch rate q
title('Pitch rate')


