clc
clear
close all
p=0.001;%Gaussian white noise intensity 
t_max=2000;
y = wgn(1,t_max,p);%Gaussian white noise
V(1)=0.6;%Initial value of V 
tt=0.1; %the delay time in the deterministic force
v=0.15; %the delay time in the random force
A=0.05; %A and ¦× are the amplitude and frequency of the seasonal oscillations of ZG, respectively.
w=0.0001; %A and ¦× are the amplitude and frequency of the seasonal oscillations of ZG, respectively
 
d0=0.5 % long-term ZG in the absence of vegetation interference and anthropogenic groundwater extraction
dinf=0.4 %the threshold of vegetation tolerance to shallow water table and insufficient aeration of the root zone
dsup=20  %the threshold of the water table below which taproots cannot extract water
 
arf=0.2;  %¦Áregulates the temporal response of the system
a=0.001;  %measures the sensitivity of the Vcc to variations of the ZG;
b=60;     %i.e.,K is a sensitivity coefficient and represents the impact of V and groundwater pumping on ZG;
 
e=0.002;  %%¦Å is the strength of extrinsic fluctuations (i.e. Gaussian white noises)
deta=0.5; %% denotes the strength of the cross-correlation between ¦Å and D 
D=0.003;  %%D are the strength of intrinsic fluctuations (i.e. Gaussian white noises)
%%%
for i=1:2000
heff_V(i)=(1+tt)*(arf*V(i)*(a.*(d0+b.*V(i)+A.*cos(w.*i)-dinf).*(dsup-d0-b.*V(i)-A.*cos(w.*i))));
Geff_V(i)=sqrt(e*(1+v)^2*V(i)^2+2*deta*sqrt(D*e)*(1+v)*V(i)+D);
V(i+1)=V(i)+heff_V(i)+Geff_V(i)*y(i);
if V(i+1)<0
V(i+1)=0;
end
if V(i+1)>1
V(i+1)=1;
end
end
 
figure(1)
plot(1:size(y,2),y,'-k')
xlabel('t')
ylabel('V(t)')
figure(2)
plot(1:size(V,2),V,'-k')
xlabel('t')
ylabel('V(t)')