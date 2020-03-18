load Cit_par_data.mat

% Short period eigenvalues (lambda_1,2)

C_1s_simp = [ (CZadot-2*muc)*c/V0          0       ;...
                Cmadot*c/V0     -2*muc*KY2*c/V0 ];
C_2s_simp = [ CZa (CZq+2*muc)*c/V0   ;...
              Cma     Cmq*c/V0       ];
          
As_simp = -1*inv(C_1s_simp)*C_2s_simp;

el_a = As_simp(1,1);
el_b = As_simp(1,2);
el_c = As_simp(2,1);
el_d = As_simp(2,2);
Discr = sqrt((-el_a-el_d)^2-4*(el_a*el_d-el_b*el_c));
lambda_1 = ((el_a+el_d)+Discr)/2;
lambda_2 = ((el_a+el_d)-Discr)/2;

% Phugoid eigenvalues (lambda_3,4)

C_1s_simp1 = [ -2*muc*c/(V0*V0)         0            0             0             ;...
                      0                 0            0             0             ;...
                      0                 0          -c/V0           0             ;...
                      0                 0            0             0             ];

C_2s_simp1 = [ CXu/V0   CXa   CZ0          0          ;...
               CZu/V0   CZa    0     (c/V0)*(2*muc)   ;...
                  0      0     0          c/V0        ;...
               Cmu/V0   Cma    0        Cmq*c/V0      ];

As_simp1 = -1*inv(C_1s_simp1)*C_2s_simp1;   % C_1s_simp is not invertible :(
% eig(As_simp1)

% Aperiodic roll eigenvalues