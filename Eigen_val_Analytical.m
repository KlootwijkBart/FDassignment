load FTISxprt-20200309_flight3.mat

%SIM CASES


% Short period eigenvalues #1

A1=4*muc^2*KY2;
B1=-2*muc*(KY2*CZa + Cmadot +Cmq);
C1= CZa*Cmq - 2*muc *Cma;


lamda_sym_spm_1= (V0/c)*(-B1+sqrt(4*A1*C1 - B1^2)*j)/(2*A1);
lamda_sym_spm_2= (V0/c)*(-B1-sqrt(4*A1*C1- B1^2)*j)/(2*A1);

% Phugoid eigenvalues #2
A2 = 2*muc*(CZa*Cmq - 2*muc*Cma);
B2= 2*muc*(CXu*Cma-Cmu*CXa) + Cmq*(CZu*CXa-CXu*CZa);
C2 = CZ0*(Cmu*CZa-CZu*Cma);

lamda_sym_phug_1= (V0/c)*(-B2+sqrt(4*A2*C2 - B2^2)*j)/(2*A2);
lamda_sym_phug_2= (V0/c)*(-B2-sqrt(4*A2*C2 - B2^2)*j)/(2*A2);


%ASYM CASES


% Aperiodic roll

lambda_asym_aperiodic = (V0/b)*Clp/(4*mub*KX2);

%dutch roll

lambda_asym_dutch_1 = (V0/b)*(2*(Cnr+2*KZ2*CYb)+sqrt(64*KZ2*(4*mub*Cnb+CYb*Cnr)-4*(Cnr+2*KZ2*CYb)^2)*j)/(16*mub*KZ2);

lambda_asym_dutch_2 = (V0/b)*(2*(Cnr+2*KZ2*CYb)-sqrt(64*KZ2*(4*mub*Cnb+CYb*Cnr)-4*(Cnr+2*KZ2*CYb)^2)*j)/(16*mub*KZ2);

%Spiral motion

lambda_asym_spiral = (V0/b)*(2*CL*(Clb*Cnr - Cnb*Clr))/(Clp*(CYb*Cnr+4*mub*Cnb)- Cnp*(CYb*Clr+4*mub*Clb));


sym = [lamda_sym_spm_1,lamda_sym_spm_2,lamda_sym_phug_1,lamda_sym_phug_2]
as  = [lambda_asym_aperiodic,lambda_asym_dutch_1,lambda_asym_dutch_2,lambda_asym_spiral]


