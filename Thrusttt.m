P0 = 101325;
T0 = 288.15;
rho0=1.225;
gamma = 1.4;
g = 9.80665;
lambda = -0.0065;
Gasconstant = 287.15;
height = flightdata.Dadc1_alt.data*0.3048;
TAT=flightdata.Dadc1_tat.data;
Cas=flightdata.Dadc1_cas.data*0.514444;
AOA = flightdata.vane_AOA.data;
mach=flightdata.Dadc1_mach.data;
TAT=flightdata.Dadc1_tat.data+273.15;
plst=[];
mlst=[];
ilst=[];
hlst=[];
vlst=[];
AOAlst=[];
machlst=[];
SATlst=[];
Veqlst=[];
Soslst=[];
Tdifflst=[];
height = [13000 12990 12990 12990 13310 13430]*0.3048;
Cas = [248 221 190 159 132 118]*0.5144444;
TAT = [-8 -10.5 -12.5 -14.2 -16.1 -16.8]+273.15;


for i = 1:49900
p = P0*(1+(lambda*height(i)/T0))^-(g/(lambda*Gasconstant));
a = (2/(gamma-1));
b = ((P0)/p);
c = (1 + (((gamma-1)*rho0*Cas(i)^2)/(2*gamma*P0)));
d = (gamma/(gamma-1));
e = (gamma-1)/gamma;
m = ((((((c^d)-1)*b+1)^e)-1)*a)^0.5;
SAT=TAT(i)/(1+((gamma-1)/2)*m^2);
SAT;
Sos=(gamma*Gasconstant*SAT)^0.5;
Vt=m*Sos;
rho=p/(Gasconstant*SAT);
Veq=Vt*(rho/rho0)^0.5;
Tisa=T0 + (lambda*height(i));
Tdiff=Tisa-SAT;

ilst = [ilst,i];
% % plst = [plst,p];
 mlst = [mlst,m];
% % hlst= [hlst,height(i)];
% % vlst= [vlst,Cas(i)];
% % AOAlst=[AOAlst,AOA(i)];
% % machlst=[machlst,mach(i)];
% SATlst=[SATlst,SAT];
% Veqlst=[Veqlst,Veq];
% Soslst=[Soslst,Sos];
Tdifflst= [Tdifflst,Tdiff];
end
disp(Tdifflst)
disp(mlst)
% %  plot(ilst,plst);
  plot(ilst,mlst,ilst,machlst);
% % plot(ilst,hlst)
% % plot(ilst,vlst)
% % plot(ilst,AOAlst)
% %  plot(ilst,machlst)
%  plot(ilst,SATlst)
% plot(ilst,Veqlst)
% % plot(ilst,Soslst)