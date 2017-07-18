function [Bt] = Bl(I,w1,w2,th,lstrip,x,y,z) 
%calculates sum of B components as function of position and current 

%w1=10e-6; w2=1e-6; lstrip = 50e-6;
a = (w1-z*(w1-w2)/lstrip)*0.5;
b = th/2;

%a is in x-axis (width), b is in y-axis (thickness)
r1 = sqrt((y-b)^2+(x+a)^2);
r2 = sqrt((y+b)^2+(x+a)^2);
r3 = sqrt((y+b)^2+(x-a)^2);
r4 = sqrt((y-b)^2+(x-a)^2);

Phi1 = atan((x+a)/(y-b));
Phi2 = atan((x+a)/(y+b));
Phi3 = atan((x-a)/(y+b));
Phi4 = atan((x-a)/(y-b));

Hx = -I/(8*pi*b*a)*((x+a)*(Phi1-Phi2)-(x-a)*(Phi4-Phi3)+(y+b)*log(r2/r3)-(y-b)*log(r1/r4));
Hy =  I/(8*pi*b*a)*((y+b)*(Phi2-Phi3)-(y-b)*(Phi1-Phi4)+(x+a)*log(r2/r1)-(x-a)*log(r3/r4));


H = sqrt(Hx^2+Hy^2);
Bt = H*4*pi*1e-7;
