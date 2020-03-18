szA = size(Aa);
szB = size(Ba);
C = eye(szA);
D = zeros(szB);
sys = ss(Aa,Ba,C,D);
[y,t] = step(sys, 10);

clf();
% y: index 1=tijd, index2=state variable, index3=input variable 
plot(t,y(:,:,2));