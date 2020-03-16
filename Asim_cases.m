load Cit_par_data.mat

% Aperiodic roll

C_1as_simp = [-2*mub*KX2^2*b/V0];

C_2as_simp = [Clp*b/(2*V0)];

As_simp_as = -1*inv(C_1as_simp)*C_2as_simp;
    
lambda_as_1 = Clp/(4*mub*KX2^2);

%dutch roll

C_1as_simp1 = [ -2*mub*b/V0               0;...
                 0               2*mub*KZ2^2*(b/V0)^2];
C_2as_simp1 = [ CYb                  -2*mub*b/V0;...
                Cnb           Cnr*b/(2*V0)];
            

As_simp_as1 = -1*inv(C_1as_simp1)*C_2as_simp1;

lambda_as_2 = (2*(Cnr+2*KZ2^2*CYb)+sqrt(64*KZ2^2*(4*mub*Cnb+CYb*Cnr)-4*(Cnr+2*KZ2^2*CYb)^2))/(16*mub*KZ2^2);      

lambda_as_3 = (2*(Cnr+2*KZ2^2*CYb)-sqrt(64*KZ2^2*(4*mub*Cnb+CYb*Cnr)-4*(Cnr+2*KZ2^2*CYb)^2))/(16*mub*KZ2^2); 

%Spiral motion

C_1as_simp2 = [0 0 0 0;...
               0 -0.5*b/V0 0 0;...
               0 0 0 0;...
               0 0 0 0];

C_2as_simp2 =[CYb CL 0 -2*mub*b/V0;...
              0 0 0 0;...
              Clb 0 0.5*Clp*b/V0 0.5*Clr*b/V0;...
              Cnb 0 0.5*Cnp*b/V0 0.5*Cnr*b/V0];
          
%As_simp_as2 = -1*inv(C_1as_simp2)*C_2as_simp2; singular matrix
          
lambda_as_4 = (2*CL*(Clb*Cnr - Cnb*Clr))/(Clp*(CYb*Cnr+4*mub*Cnb)- Clp*(CYb*Clr+4*mub*Clb));

