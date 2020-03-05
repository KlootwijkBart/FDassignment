s = tf('s')

K = 1
Kr = 0.9
h = 20/((s+1)*(s+4))
h2 = 1/s

Hi = feedback(K*h1,Kr)

Hyr = feedback(Hi*h2,1)

a = Hyr.num{1}
b = Hyr.den{1}

Her = feedback(1,Hi*h2)

c = Her.num{1}
d = Her.den{1}

t = 0:0.01:20

y = step(Her,t)

t = lsiminfo(y,t,'SettlingTimeThreshold', 0.05)

x = (20/((s+1)*(s+4)))*(Kr+1/s)
rltool(x)

