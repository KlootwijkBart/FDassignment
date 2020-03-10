szA = size(Aa);
szB = size(Ba);
C = eye(szA);
D = zeros(szB);
sys = ss(Aa,Ba,C,D);
[y,t] = step(sys, 10);

clf();
plot(t,y(:,3,1));