load Cit_par_data.mat
load FTISxprt-20200309_flight3.mat
% =========================================================================
%   SYMMETRIC EQUATIONS
% =========================================================================

C_1s = [ -2*muc*c/(V0*V0)         0            0             0             ;...
                0        (CZadot-2*muc)*(c/V0)    0             0             ;...
                0                 0          -c/V0           0             ;...
                0            Cmadot*c/V0       0   -2*muc*KY2*c*c/(V0*V0)  ];

C_2s = [ CXu/V0   CXa   CZ0       CXq*c/V0      ;...
         CZu/V0   CZa  -CX0  (c/V0)*(CZq+2*muc) ;...
            0      0     0          c/V0        ;...
         Cmu/V0   Cma    0        Cmq*c/V0      ];

C_3s = [ CXde  ;...
         CZde  ;...
           0   ;...
         Cmde  ] ;

As = -1*inv(C_1s)*C_2s;
Bs = -1*inv(C_1s)*C_3s;



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

     
% Calculate A and B matrices of state-space model     
Aa = -1*inv(C_1a)*C_2a;
Ba = -1*inv(C_1a)*C_3a;


Ca = eye(size(Aa));
Da = [0];

sysa = ss(Aa,Ba,Ca,Da);

[y,t] = step(sysa,10);
plot(t,y(:,:,1))
legend('Sideslip','Bank angle', 'Roll Rate','Yaw Rate')

