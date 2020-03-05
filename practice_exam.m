Practice exam

s = tf('s');
K = 0.8551;
H = 0.5*(2*s+1)/(s^2*(s^2+0.4*s+4));
Hol = K*H;
Hcl = feedback(Hol, 1)
bode(Hcl)

m1 = 12;
m2 = 1000;
k1 = 46000;
k2 = 8000;
b = 2800;

A = [ 0 1 0 0;...
     -(k1/m1)-(k2/m1) -b/m1 k2/m1 b/m1;...
     0 0 0 1;...
     k2/m2 b/m2 -k2/m2 -b/m2];
B = [0;...
    k1/m1;...
    0;...
    0];
C = [0 0 1 0;...
    k2/m2 b/m2 -k2/m2 -b/m2];
D = [0];
sys = ss(A,B,C,D);
    
s = tf('s');
K = 1;
Kr = 0.04;
h1 = 20/((s+1)*(s+4));
Hinner = feedback(K*h1, Kr);
Houter = feedback(Hinner/s, 1);
h = feedback(1, Hinner/s);

t = 0:0.01:60;
y = step(Houter,t);
plot(t,y)

idx = y> max(y)*1.05 | y< max(y)*0.95;
tx = t(idx); tset = tx(end)+0.01;

Hol = Hinner/s;
openloop = K*h/s * 
rltool(Hol)

clear
% parameter values
m = 3.5;
b = 9;
k = 60;
% matrices
A = [ -b/m -k/m; 1 0];
B = [ 1/m; 0];
C = [ 0 1];
D = [ 0 ];
% system
sys = ss(A, B, C, D);


t = 0:0.01:30;

w1 = 3.9357078; % rad/s
w2 = 3.7197762;
u1 = sin(w1*t);
u2 = sin(w2*t);

y1 = lsim(sys,u1,t);
y2 = lsim(sys,u2,t);

plot(t, [y1 y2])


    