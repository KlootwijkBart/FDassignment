load Cit_par_data.mat

% Short period eigenvalues (lambda_1,2)

C_1s_simp = [ (CZadot-2*muc)*c/V0          0       ;...
                Cmadot*c/V0     -2*muc*KY2*c/V0 ];
C_2s_simp = [ CZa (CZq+2*muc)*c/V0   ;...
              Cma     Cmq*c/V0       ];
          
As_simp = -1*inv(C_1s_simp)*C_2s_simp;

%% 
el_a = As_simp(1,1);
el_b = As_simp(1,2);
el_c = As_simp(2,1);
el_d = As_simp(2,2);
Discr = sqrt((-el_a-el_d)^2-4*(el_a*el_d-el_b*el_c));
lambda_1 = ((el_a+el_d)+Discr)/2;
lambda_2 = ((el_a+el_d)-Discr)/2;
%%

% Phugoid eigenvalues (lambda_3,4)

C_1s_simp1 = [ -2*muc*c/(V0*V0)         0            0             0             ;...
                      0                 0            0             0             ;...
                      0                 0          -c/V0           0             ;...
                      0                 0            0             0             ];

C_2s_simp1 = [ CXu/V0   CXa   CZ0          0          ;...
               CZu/V0   CZa    0     (c/V0)*(2*muc)   ;...
                  0      0     0          c/V0        ;...
               Cmu/V0   Cma    0        Cmq*c/V0      ];

As_simp1 = -1*pinv(C_1s_simp1)*C_2s_simp1;   % C_1s_simp is not invertible :(

% Aperiodic roll eigenvalues

lambda_a1 = Clp/(4*mub*KX2*b/V0);

% Dutch roll eigenvalues

C_1a_simp = [ -2*mub*b/V0             0            ;...
                   0        -2*mub*KZ2*b*b/(V0*V0) ];
C_2a_simp = [ CYb -4*mub*b/(2*V0) ;...
              Cnb   Cnr*b/(2*V0)  ];

Aa_simp = -1*inv(C_1a_simp)*C_2a_simp;
el_a1 = Aa_simp(1,1);
el_b1 = Aa_simp(1,2);
el_c1 = Aa_simp(2,1);
el_d1 = Aa_simp(2,2);
Discr1 = sqrt((-el_a1-el_d1)^2-4*(el_a1*el_d1-el_b1*el_c1));

lambda_a2 = ((el_a1+el_d1)+Discr1)/2;
lambda_a3 = ((el_a1+el_d1)-Discr1)/2;

% Spiral motion eigenvalue

C_1a_simp1 = [ 0     0     0 0 ;...
               0 -0.5*b/V0 0 0 ;...
               0     0     0 0 ;...
               0     0     0 0 ];

C_2a_simp1 = [ CYb     CL         0         -4*mub     ;...
               0       0       b/(2*V0)        0       ;...
               Clb     0     Clp*b/(2*V0) Clr*b/(2*V0) ;...
               Cnb     0     Cnp*b/(2*V0) Cnr*b/(2*V0) ];
           
Aa_simp1 = -1*pinv(C_1a_simp1)*C_2a_simp1;

test_matrix = [ -2*muc*c/(V0^2) 0 0;...
                0 0 0;...
                0 -c/V0 0];
test_matrix2 = [ CXu/V0 CZ0 0;...
                CZu/V0 0 2*muc*c/V0;...
                0 0 c/V0];
