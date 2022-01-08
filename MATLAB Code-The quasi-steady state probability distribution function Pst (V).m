clc
clear
close all
syms tao arfa a1 K Z0 Zinf Zsup ebsl gui ba Heff Geff Bv FF AA D VV V
D=0.01   %%D are the strength of intrinsic fluctuations (i.e. Gaussian white noises)
ebsl=0.001%%¦Å is the strength of extrinsic fluctuations (i.e. Gaussian white noises)
K=45    %K is a sensitivity coefficient and represents the impact of V and groundwater pumping on ZG;
a1=0.02 % measures the sensitivity of the Vcc to variations of the ZG;
gui=0.15 %vthe delay time in the random force
tao=0.1  %¦Óthe delay time in the deterministic force
arfa=0.15 %¦Áregulates the temporal response of the system
Z0=0.5   %Z0 is the long-term ZG in the absence of vegetation interference and anthropogenic groundwater extraction
Zinf=0.4  %the threshold of vegetation tolerance to shallow water table and insufficient aeration of the root zone
Zsup=20  %the threshold of the water table below which taproots cannot extract water
ba=0.5    %%denotes the strength of the cross-correlation between ¦Å and D 
 
Heff=(1+tao).*(arfa.*VV.*(a1.*(K.*VV+Z0-Zinf).*(Zsup-K.*VV-Z0)-VV))
Geff=sqrt(ebsl.*(1+gui).^2.*V.^2+2.*ba.*sqrt(D.*ebsl).*(1+gui).*V+D)
Bv=ebsl.*(1+gui).^2.*VV.^2+2.*ba.*sqrt(D.*ebsl).*(1+gui).*VV+D
AA=int(Heff./Bv, VV, 0, V)
BB=int(exp(int(Heff./Bv, VV, 0, V))./Geff, V, 0, 1)
N=1./(BB)
PST=(N./Geff).*exp(AA)
VK=[0:0.01:1];
PSTV=subs(PST,V,VK)
plot(VK,PSTV);