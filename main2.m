clear all;
clc;
%%%%%%%%%%%%%%%
%pathogen size
rhiv = 60e-9;
rherp = 75e-9;
%pathogen mass
mhiv = (4/3)*pi*rhiv^3*1000;
mherp = (4/3)*pi*rherp^3*1000;
%set up variables describing channel
%%{
for g = 15:100
w1=g*1e-6; %m
w2=15e-6;
th=0.5e-6;
J0max = 2e10;
%I=0.15; %A
I = J0max*w2*th;
lstrip=200e-6;
xd=0; %particle centered above conductor
y=0.75e-6; %dist. from centre of strip (need 2 calc. w DLVO force equilibrium)
%z-axis is direction along channel
%%}

%set up z-array from 0 to length of strip for later use
res = 0.1; 
z = ((1:res:(lstrip*1e6+1))-1)*1e-6;

%%%%%%%%%%%%%%%%
%calculate magnetic flux density along channel
Bz = @(z) Bl(I,w1,w2,th,lstrip,xd,y,z);

B_dist = arrayfun(Bz, z);
%plot (z,B_dist);

%%%%%%%%%%%%%%%%
%calculate gradient of flux density-squared, along channel
Bz_sqr = @(z) Bl(I,w1,w2,th,lstrip,xd,y,z).^2;
B_sqr_dist = arrayfun(Bz_sqr, z);
%plot (z, B_sqr_dist);
dB2 = diff(B_sqr_dist);
dz = diff(z);
dB2dz = dB2./dz;
z2 = z; z2(end) = [];
%plot(z2,dB2dz)

%%%%%%%%%%%%%%

%%{
%set up variables describing particle
%rp = 0.25e-6; %m
Xi = 3.09;
mag_fract = 1;
rp = 0.25e-6;
%rp = 0.75e-6;
%Xi = 0.18850;
%mag_fract = 0.1168;
V = (4/3)*pi*rp^3;
Vm = V*mag_fract;
%rho = 5530; %kg/m^3
rho = 1100;
Mp = rho*V;
mu0 = 4*pi*1e-7;
dyn_visc=1e-3; %1 for water at 20°, 0.798 for wtr at 30° 
%dyn_visc = 3*10^-3; % approx. 3-4 for blood 
fd =1; %drag coefficient

k_drag = 6*pi*dyn_visc*rp*fd;
k_dragl = 6*pi*dyn_visc*(rp+rhiv)*fd;
%%}

%calculate force array due to gradient
F_mag = dB2dz*(Vm*Xi/(2*mu0));

%plot (z2,F_mag);
%p = polyfit(z2,F_mag,2); %uses 4th order polyfit, if change order also need to change ODE function!
%plot(z2,polyval(p,z2),z2,F_mag);

%%%%%%%%%%%%%%
%calculate traversal time using terminal speed
%case 1: average, using force
%vterm = mean(F_mag)/k_drag;
%vterml = mean(F_mag)/k_dragl;

%for i = 1:200
    %array(i) = F_mag(i)/k_drag;
    %array(i) = auxiliary(F_mag,(i-1)*1e-6,lstrip);
%end
%t=1:200; plot(t,array);

%case 2: average, using gradient:
%up = (rp^2*Xi)*mean(dB2dz)/(9*mu0*dyn_visc);

%ttrav = lstrip/vterm
%ttravl = lstrip/vterml
%deltat = ttravl-ttrav

%taccel = up/(mean(F_mag)/Mp);

%case 3: find actual time taken from v(x) and x
%aux =  @(x) 1/auxiliary(F_mag,x,lstrip);
%aux(3e-6)
%ttrav = integral(aux,0,lstrip);
%{
%case 4: ...different integration approximation...
varray = F_mag./k_drag;
invarray = k_drag./F_mag;
invarrayl = k_dragl./F_mag;
ttrav(g) = trapz(invarray)*(lstrip/length(F_mag));
ttravl(g) = trapz(invarrayl)*(lstrip/length(F_mag));

%case 5: ...test w. rectangular integral approx...
%ttrav2 = sum(invarray)*(lstrip/2000);

%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
%set up ODE for position of particle...
%F_total = F_mag + F_drag; 
%m*acc = Fmag - 6*pi*eta*rp*fd*vel_p
%hence 0 = (m)*d2zdt2 + (6*pi*eta*rp*fd)*dzdt - Fmag(O(x^4))
%odeforce = @(z,t) odefun(k_drag,Mp,z,p2,t);

%[t,z]=ode23t(@odefun,[0,1],[1e-6,5e-6]);
%plot(t,z(:,1))

%%%%%%%%%%%%%%
%diagnose convergence problems with Euler integration

x1=1e-6;
v1=0;
um=1e-6;
dt=1e-5;
n=0;
t =0:dt:0.1;
ad = zeros(length(t));
vd = zeros(length(t));
xd = zeros(length(t));
tic
for t = 0:dt:0.1
    n=n+1;
    if(x1<=lstrip)
    %Ftd(n) = interp1(F_mag,x1/lstrip*200)-k_drag*v1;
    Ftd(n) = mean(F_mag)-k_drag*v1;
    ad(n) = max(0,Ftd(n)/Mp); 
    %ad(n) = a;
    v1 = v1 + ad(n)*dt;
    vd(n) = v1;
    x1 = x1 + v1*dt;
    xd(n)=x1;
    end
    if(x1<=lstrip)
        ad(n)=0; v1 = 0; vd(n)=0; xd(n)=0;
    end
end
toc

t=0:dt:0.1;
subplot(2,2,1)
plot(t,ad,'g')
subplot(2,2,2)
plot(t,vd,'b')
subplot(2,2,3)
plot(t,xd,'r')

%plot((1:35),ttrav)
end
%plot((1:100),ttrav)
%}