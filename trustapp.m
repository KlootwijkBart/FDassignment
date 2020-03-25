%%flight data Stationary (excel)
height = [13000 12990 12990 12990 13310 13430 13420 13520 13570 13660 12830 12390 11380 11750 11630]*0.3048;
Cas = [248 221 190 159 132 118 150 140 131 121 160 168 180 153 151]*0.5144444;
TAT = [-8 -10.5 -12.5 -14.2 -16.1 -16.8 -15.2 -16 -16.5 -17.5 -13.5 -12.5 -9.8 -12 -12]+273.15;
angleoa=[1.6 2.4 3.6 5.7 8.7 11 6.4 7.5 8.7 10.2 5.5 4.9 4.1 6.1 6.3]*0.0174533; %[radians]
ffl=[697 585 497 426 432 418 375 373 371 368 385 390 403 394 395]/3600*0.453592;
ffr=[762 640 538 470 478 456 410 408 406 403 420 427 441 431 430]/3600*0.453592;
Fused=[460 497 545 577 622 648 703 719 732 748 770 786 817 835 861]*0.453592;
Deltaelev=[0 0 0 0 0 0 -0.7 -1.1 -1.7 -2.4 -0.3 0.1 0.4 -0.6 -1.4]*0.0174533; %radians
Detlatrim=[0 0 0 0 0 0 2.1 2.1 2.1 2.1 2.1 2.1 2 2 2]*0.0174533; %radians
F_e=[0 0 0 0 0 0 0 -15 -32 -42 18 37 74 0 -35];

%%reference data Stationary (excel)
% 
% height = [5010	5020	5020	5030	5020	5110 6060	6350	6550	6880	6160	5810	5310 5730	5790]*0.3048;
% Cas = [249	221	192	163	130	118 161	150	140	130	173	179	192 161	161]*0.5144444;
% TAT = [12.5	10.5	8.8	7.2	6	5.2 5.5	4.5	3.5	2.5	5.0	6.2	8.2 5.0	5.0]+273.15;
% angleoa=[1.7 2.4	3.6	5.4	8.7	10.6 5.3	6.3	7.3	8.5	4.5	4.1	3.4 5.3	5.3]*0.0174533; %radians
% ffl=[798	673	561	463	443	474 462	458	454	449	465	472	482 471	468]/3600*0.453592;
% ffr=[813	682	579	484	467	499 486	482	477	473	489	496	505 493	490]/3600*0.453592;
% Fused=[360	412	447	478	532	570 664	694	730	755	798	825	846 881	910]*0.453592;
% Deltaelev=[0 0 0 0 0 0 0	-0.4	-0.9	-1.5	0.4	0.6	1 0	-0.5]*0.0174533; %radians
% Detlatrim=[0 0 0 0 0 0 2.8	2.8	2.8	2.8	2.8	2.8	2.8 2.8	2.8]*0.0174533; %radians
% F_e=[0 0 0 0 0 0 0	-23	-29	-46	26	40	83 0	-30];

%%flight data (all data with 10hz)
% height = flightdata.Dadc1_alt.data*0.3048;
% TATref=flightdata.Dadc1_tat.data+273.15;
% SATref=flightdata.Dadc1_sat.data+273.15;
% Cas=flightdata.Dadc1_cas.data*0.514444;
% Tas=flightdata.Dadc1_tas.data*0.514444;
% AOA = flightdata.vane_AOA.data*0.0174533;
% mach=flightdata.Dadc1_mach.data;
% TAT=flightdata.Dadc1_tat.data+273.15;
% Fused1=flightdata.lh_engine_FU.data*0.453592;
% Fused2=flightdata.rh_engine_FU.data*0.453592;

%% constants

P0 = 101325;
T0 = 288.15;
rho0=1.225;
gamma = 1.4;
g = 9.81;
lambda = -0.0065;
Gasconstant = 287.15;
Dynviscosity=0.0000163;

%%plane properties
Total_empty_weight=9165*0.453592;
mass_pass=803;    
% mass_pass=695;          %reference data 
% Fuelt1=4050*0.453592;       %reference data 
Fuelt1=4000*0.453592;   
SurfaceWing=30;
chord=2.0569;
s_fuel_flow=0.048;
s_mass=60500;
Cmd=-1.1642;
Cmtc=-0.0064;
span      = 15.911;	  % wing span [m]
Aspect      = span^2/SurfaceWing;           % wing aspect ratio [ ]

%% list for all values

plst=[];
mlst=[];
ilst=[];
hlst=[];
vlst=[];
AOAlst=[];
machlst=[];
SATlst=[];
SATreflst=[];
Veqlst=[];
Tasreflst=[];
Vtlst=[];
Soslst=[];
Tdifflst=[];
rholst=[];
Cllst=[];
Cl2lst=[];
Cdlst=[];
reynoldslst=[];
Tclst=[];
Tcstlst=[];
Deleveq=[];
Veqredlst=[];
F_nredlst=[];
Liftlst=[];
soslst=[];


%% %calculates all values for each data point
for i = 1:length(height)
    %lift calcuation
    Fuel_left=Fuelt1-Fused(i);
%     Fuel_left=Fuelt1-(Fused1(i)+Fused2(i));   #used for reference data in
%     flightdata (all data)
    Total_weight = Fuel_left+mass_pass+Total_empty_weight;
    Lift=Total_weight*g;
        
    %pressure
    p = P0*(1+(lambda*height(i)/T0))^-(g/(lambda*Gasconstant));
    
    %mach calculation with Pressure and Calibrated airspeed (p,Cas)
    m = (((((((1 + (((gamma-1)*rho0*Cas(i)^2)/(2*gamma*P0)))^(gamma/(gamma-1)))-1)*((P0)/p)+1)^((gamma-1)/gamma))-1)*(2/(gamma-1)))^0.5;
    
    %static air temperature from Total air temperature
    SAT=TAT(i)/(1+((gamma-1)/2)*m^2);
    SATlst=[SATlst,SAT];
    
    %speed of sound from SAT
    Sos=(gamma*Gasconstant*SAT)^0.5;
    
    %True airspeed from mach and Sos
    Vt=m*Sos;
    Vtlst=[Vtlst,Vt];
    
%     Tasreflst=[Tasreflst,Tas(i)];    #true airspeed from reference flight
%     data file

    %Isa standard temperature
    Tisa=T0 + (lambda*height(i));
    %rho from T isa
    rho=p/(Gasconstant*Tisa);
    
    %equivalent airspeed calibrated with rho
    Veq=Vt*(rho/rho0)^0.5;
    
    %difference in SAT and T isa fro thrust calculation
    Tdiff=SAT-Tisa;
    
    %Cl calulation from rho, Veq and surface area
    Cl=Lift/(0.5*rho0*Veq*Veq*SurfaceWing);
    Cl2=Cl^2;
    
    %reynolds number range 
    reynoldslst=[reynoldslst,(chord*rho0*Veq)/Dynviscosity];
    
    %matrix for the thrust exe file
    thrustmat(i,1)=height(i);
    thrustmat(i,2)=m;
    thrustmat(i,3)=Tdiff;
    thrustmat(i,4)=ffl(i);
    thrustmat(i,5)=ffr(i);
    
    %Reduced values for Veq and mach 
    Veqred=Veq*((s_mass/Lift)^0.5);
    Meqred=Veqred/Sos;
    
    %matrix for the thrust exe file but for standard values
    thrustmat(i+length(height),1)=height(i);
    thrustmat(i+length(height),2)=Meqred;
    thrustmat(i+length(height),3)=Tdiff;
    thrustmat(i+length(height),4)=s_fuel_flow;
    thrustmat(i+length(height),5)=s_fuel_flow;
    
    Cllst=[Cllst,Cl];
    Cl2lst=[Cl2lst,Cl2];
    rholst= [rholst,rho];
    Veqlst=[Veqlst,Veq];
    Liftlst=[Liftlst,Lift];
    Veqredlst=[Veqredlst,Veqred];
    mlst=[mlst,m];
    Vtlst=[Vtlst,Vt];
    soslst=[soslst,Sos];
    
end

%% Thrust calcation with exe file (outside the for loop)
save('D:\tu delft\bsc 3rd year\SVV\matlab.dat','thrustmat','-ascii');
system('D:\tu delft\bsc 3rd year\SVV\thrust(1).exe');
tdata=load('D:\tu delft\bsc 3rd year\SVV\thrust.dat');

%% New for loop for last calculations where thrust is needed

for i = 1:length(height)
    T1=tdata(i,1);
    T2=tdata(i,2); 
    T1s=tdata(i+length(height),1);
    T2s=tdata(i+length(height),2);
 
    %Cd from thrust, rho, Veq and surface area
    Cd=(T1+T2)/(0.5*rho0*Veqlst(i)*Veqlst(i)*SurfaceWing);    
    %thrust coefficient NOTE Cd is Tc lift=drag
    Tclst=[Tclst,Cd];
    
    %Standard thrust coefficient
    Tcstlst=[Tcstlst,((T1s+T2s)/(0.5*rho0*Veqlst(i)*Veqlst(i)*SurfaceWing))];
       
    DeltaeEq=Deltaelev(i)-((Cmtc/Cmd)*(((T1s+T2s)/(0.5*rho0*Veqlst(i)*Veqlst(i)*SurfaceWing))-((T1+T2)/(0.5*rho0*Veq*Veq*SurfaceWing))));
    
    Deleveq=[Deleveq,DeltaeEq];  
    F_nredlst=[F_nredlst,F_e(i)*((s_mass/Liftlst(i)))];
    Cdlst=[Cdlst,Cd];
end

%% %Cl zero lift and CD0
CLcoef = polyfit(angleoa(1:6),Cllst(1:6),1);
aoa0lift=-1*CLcoef(2)/CLcoef(1);
CDcoef = polyfit(Cdlst(1:6),Cl2lst(1:6),1);
Cd0=-1*CDcoef(2)/CDcoef(1);
Cd0=0.025;
oswald=(CDcoef(1)*1+CDcoef(2))/((1-Cd0)*pi*Aspect);
number1=CLcoef(1)*0.0174533*5.3+CLcoef(2);

Cdnewlst=[];
for i = 1:6
    Cdnew=Cd0+(Cl2lst(i)/(pi*oswald*Aspect));
    Cdnewlst=[Cdnewlst,Cdnew];
end
ClCdcoef=polyfit(Cdnewlst,Cllst(1:6),2);
xposclcd=(ClCdcoef(3)/(ClCdcoef(1)))^0.5;
yposclcd=ClCdcoef(1)*xposclcd^2+ClCdcoef(2)*xposclcd+ClCdcoef(3);
glideratio=yposclcd/xposclcd;
glidexlst=[0,xposclcd];
glideylst=[0,yposclcd];

%% Plotting all required plots %Cl-Cd

[xlst2new,ylst2new]=lineplot(Cdnewlst,Cllst(1:6),2);
%[xlst2,ylst2]=lineplot(Cdlst(1:6),Cllst(1:6),2);
plot(xlst2new,ylst2new,Cdnewlst,Cllst(1:6),'O');%,xlst2,ylst2);%Cdnewlst,Cllst(1:6),'O',
grid on;
% title('Cd - Cl');
xlabel('Cd');
ylabel('Cl')
% legend('Trendline calculated data','Measured data');
txt="Configuration: Landing gear up, Flaps neutral";
text(Cd0+0.0006,1.15,txt)
txt="Mach range ["+ mlst(6)+" , "+mlst(1)+"]";
text(Cd0+0.0006,1.1,txt)
txt = "Reynolds number range  ["+ min(reynoldslst(1:6))+" , "+max(reynoldslst(1:6))+"]";
text(Cd0+0.0006,1.05,txt)
% txt="Trendline formula calculated data: "+ClCdcoef(1)+" * Cd^2 + "+ClCdcoef(2)+" * Cd + "+ClCdcoef(3);
% text(Cd0+0.0006,1,txt)


%% angle of attack vs cl plot
% [xlst,ylst]=line(angleoa(1:6),Cllst(1:6));
% set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
% plot(angleoa(1:6),Cllst(1:6),'O',xlst,ylst)
% grid on;
% title('Cl - \alpha');
% xlabel('\alpha [rad]')
% ylabel('Cl');
% txt="Configuration: Landing gear up, Flaps up";
% text(aoa0lift-0.03,0.95,txt)
% txt="Mach range ["+ mlst(6)+" , "+mlst(1)+"]";
% text(aoa0lift-0.03,0.9,txt)
% txt = "Reynolds number range  ["+ min(reynoldslst(1:6))+" , "+max(reynoldslst(1:6))+"]";
% text(aoa0lift-0.03,0.85,txt)
% txt="Trendline formula: "+CLcoef(1)+" * \alpha + "+CLcoef(2);
% text(aoa0lift-0.03,0.8,txt)
% txt="\alpha for Cl = 0: "+aoa0lift;
% text(aoa0lift-0.03,0.75,txt)

%% cd vs  Cl^2 plot
% [xlst,ylst]=line(Cdlst(1:6),Cl2lst(1:6));
% plot(Cdlst(1:6),Cl2lst(1:6),'O',xlst,ylst)
% grid on;
% % title('Cd - Cl^2');
% xlabel('Cd');
% ylabel('Cl^2')
% txt="Configuration: Landing gear up, Flaps neutral";
% text(Cd0+0.0006,0.85,txt)
% txt="Mach range ["+ min(mlst(7:15))+" , "+max(mlst(7:15))+"]";
% text(Cd0+0.0006,0.9,txt)
% txt = "Reynolds number range  ["+min(reynoldslst(7:15))+" , "+max(reynoldslst(7:15))+"]";
% text(Cd0+0.0006,0.95,txt)
% txt="Trendline formula: "+CDcoef(1)+" * Cd "+CDcoef(2);
% text(Cd0+0.0006,0.8,txt)


%% Elevator trim deflection vs Equivalent aispeed
%14 to 15 is with mustafa in front
% [xlst,ylst]=line(Veqredlst(14:15),Deleveq(14:15));
% [xlst2,ylst2]=lineplot(Veqredlst(7:13),Deleveq(7:13),2);
% plot(Veqredlst(7:13),Deleveq(7:13),'O',Veqredlst(14:15),Deleveq(14:15),'O',xlst2,ylst2,xlst,ylst)
% set(gca, 'YDir','reverse')
% grid on;
% % title('Elevator trim curve');
% xlabel('V_{eq} (reduced) [m/s]');
% ylabel('Reduced elevator deflection [radians]')
% distancey=(max(Deleveq(7:13))-min(Deleveq(7:13)))/length(Deleveq(7:13));
% txt="Configuration: Landing gear up, Flaps neutral";
% text(min(Veqlst(7:13)),max(Deleveq(7:13)),txt)
% 
% txt="Mach range ["+ min(mlst(7:15))+" , "+max(mlst(7:15))+"]";
% text(min(Veqlst(7:13)),(max(Deleveq(7:13))+(distancey/2)),txt)
% 
% 
% txt = "Reynolds number range  ["+ min(reynoldslst(7:15))+" , "+max(reynoldslst(7:15))+"]";
% text(min(Veqlst(7:13)),max(Deleveq(7:13))+distancey,txt)
% 
% % coef=polyfit(Veqredlst(7:13),Deleveq(7:13),2);
% % txt="Trendline formula: "+coef(1)+" * Veq^2 "+coef(2)+"* Veq + "+coef(3);
% % text(min(Veqlst(7:13)),max(Deleveq(7:13))+distancey*1.5,txt)
% 
% legend('cg after','cg forward');

%% elevator control force vs Equivalent airspeed
% [xlst,ylst]=line(Veqredlst(14:15),F_nredlst(14:15));
% [xlst2,ylst2]=lineplot(Veqredlst(7:13),F_nredlst(7:13),2);
% plot(Veqredlst(7:13),F_nredlst(7:13),'O',Veqredlst(14:15),F_nredlst(14:15),'O',xlst2,ylst2,xlst,ylst)
% set(gca, 'YDir','reverse')
% grid on;
% % title('Elevator control force curve');
% xlabel('V_{eq} (reduced) [m/s]');
% ylabel('F_e (reduced) [N]')
% distancey=(max(Deleveq(7:13))-min(F_nredlst(7:13)))/length(F_nredlst(7:13));
% legend('cg after','cg forward');
% txt="Configuration: Landing gear up, Flaps neutral";
% text(min(Veqlst(7:13))-10,max(F_nredlst(7:13))-distancey*2,txt)
% 
% txt="Mach range ["+ min(mlst(7:15))+" , "+max(mlst(7:15))+"]";
% text(min(Veqlst(7:13))-10,max(F_nredlst(7:13))-distancey,txt)
% 
% txt = "Reynolds number range  ["+ min(reynoldslst(7:15))+" , "+max(reynoldslst(7:15))+"]";
% text(min(Veqlst(7:13))-10,max(F_nredlst(7:13)),txt)

% coef=polyfit(Veqredlst(7:13),F_nredlst(7:13),2);
% txt="Trendline formula: "+coef(1)+" * Veq^2 "+coef(2)+"* Veq + "+coef(3);
% text(min(Veqlst(7:13))-10,max(F_nredlst(7:13))-distancey*3,txt)

%% angle of attack vs elevator trim deflection
% 14 to 15 is with mustafa in front
% [xlst,ylst] = line(angleoa(7:13),Deleveq(7:13));
% coef=polyfit(angleoa(7:13),Deleveq(7:13),1);
% plot(angleoa(7:13),Deleveq(7:13),'O',angleoa(14:15),Deleveq(14:15),'O',xlst,ylst)



%% function to create trendline that end in y=0
function [xlst,ylst] = line(lst1,lst2)
coef = polyfit(lst1,lst2,1);
y0=-1*coef(2)/coef(1);
ylst=[];
    for i = 1:length(lst1)
    y=coef(1)*lst1(i)+coef(2);
    ylst=[ylst,y];
    
    end
ylst=[ylst,y0*coef(1)+coef(2)];
xlst=[lst1,y0];
end

%function for second deg polynomial

function [xlst2,ylst2]= lineplot(lst1,lst2,deg)
coef=polyfit(lst1,lst2,deg);
xspace=linspace(min(lst1)-0.01,max(lst1),100);
xlst2=[];
ylst2=[];
for i = 1:100
xlst2=[xlst2,xspace(i)];
if deg==1
    yspace=coef(1)*xspace(i)+coef(2);
    ylst2=[ylst2,yspace];
end
if deg==2
    yspace=(coef(1)*xspace(i)^2+coef(2)*xspace(i)+coef(3));
    ylst2=[ylst2,yspace];
end
end
end





