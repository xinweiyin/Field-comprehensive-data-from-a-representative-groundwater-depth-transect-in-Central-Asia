clc
clear
close all
PP=[]'; %%%%%%%%%% Model parameter values 
CT=PP(:,1)      %%the salt-stress threshold
b=PP(:,2)       %%the rate of transpiration decline above CT
Czr=PP(:,3)     %%Mean salt concentration in root zone
Czc=PP(:,4)     %%Mean salt concentration in [Zr, Zc] zone
Czg=PP(:,5)     %%Mean salt concentration in [Zc, ZG] zone
ETmax=PP(:,6)   %%Maximum evapotranspiration rate
n=PP(:,7)       %Averaged soil porosity of root zone
Zr=PP(:,8)      %Root zone depth
a=PP(:,9)       %Average rainfall depth
LAMDA=PP(:,10)  %Average rainfall frequency
Data=PP(:,11)   %Threshold of rainfall interception
M=PP(:,12).*10  %Mean salt mass in root zone
s1=PP(:,13)     %sfc--Rainfall exceeding s1 is instantaneously lost via deep percolation at s = s1 
Zg=PP(:,14)      %%groundwater table depth
Zc=PP(:,15)      %%Maximum height of capillary rise
s1x=PP(:,16)     %%Soil moisture at stress point of root zone
s2=PP(:,17)      %%mean relative soil water in stable zone£»
sc2=PP(:,18)     %%virtual wilting point of stable zone
s2x=PP(:,19)     %%Soil moisture at stress point of [Zr, Zc] zone
LAI=PP(:,20)     %Leaf area index (LAI)
sw=PP(:,21)      %degree of soil saturation at wilting point
sx=PP(:,22)      %Mean soil saturation at stress point
%S=PP(:,23)      %Measured mean relative soil water in root zone

lamda=LAMDA.*exp(-Data./a)
sc1=(M./(n.*Zr))./(1./b+CT)./100
st=(M./(n.*Zr))./CT.*10
r=n.*Zr./a
fei=ETmax./(n.*Zr.*s1)
x=1+b.*CT
%%%%%%%%%NT for highly sensitive and moderately tolerant species
V1=lamda./(fei.*x);
V2=r.*s1-r.*sc1;
T1=gamma(V1)-gammainc(V1,V2);
NT=(exp(r.*sc1).*(r.*(s1-sc1)).^(lamda./(fei.*x)))./T1

%%%%%%%%%NT for highly tolerant species
V1=lamda./(fei.*x);
M1=r.*st-r.*sc1;
M2=lamda./fei;
M3=r.*st;
M4=r.*s1;
T0=gamma(V1)-gammainc(V1,M1);
T2=gammainc(M2,M3)-gammainc(M2,M4);
NT=(exp(r.*sc1).*(r.*(st-sc1)).^(lamda./(fei.*x)))./(x.*(T0+T2.*x.*(r.*st).^(lamda./fei)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms s
s=0:0.001:1
p1=NT.*(exp(-r.*s)./(s-sc1)).*((s-sc1)./(st-sc1)).^(lamda./(fei.*x));
p=(sc1<s<=st).*(p1);
%plot(s,p)
hold on
S=s(find(p==max(p))) 
Y=max(p); 
%%%%%%%%%%%%%%%%%%%%%
syms z
Z=Zr./100;  %Zr
B=(integral(@(z)exp(-2.6111.*z),0,+inf))
Z0=0                                     
A=(integral(@(z)exp(-2.6111.*z),0,Z))
Fzr=A/B 
%%%%%%%%%%%%%%%%%%%%%Fzc
syms z
Z=Zg-Zc;                         
B=(integral(@(z)exp(-2.6111.*z),0,+inf)); 
A=(integral(@(z)exp(-2.6111.*z),Zr./100,Z./100));
Fzc=A./B 
%%%%%%%%%%%%%%%%%%%%%Bzr
%BZR=(S-sc1)./(s1x-sc1);
%Bzr1=min(1,BZR);
%Bzr=max(0,Bzr1)
Bzr=0.*(S<sw)+((S-sw)./(sx-sw)).*(S>=sw & S<sx)+1.*(S>=sx)
%%%%%%%%%%%%%%%%%%%%%Bzc
%BZC=(s2-sc2)./(s2x-sc2);
%Bzc1=min(1,BZC)
%Bzc=max(0,Bzc1)
Bzc=0.*(s2<sw)+((s2-sw)./(sx-sw)).*(s2>=sw & s2<sx)+1.*(s2>=sx)
%%%%%%%%%%%%%%%%%%%%%Rzg
syms z
Z=Zg;                                      
B=(integral(@(z)exp(-2.6111.*z),0,+inf));
A=(integral(@(z)exp(-2.6111.*z),Zr./100,Zg./100));
C=Z./(Z-Zc);                              
func=A./B;
Rzg=C.*func     
FF=[Fzr,Fzc]
BB=[Bzr,Bzc]
%%%%%%%%%%%%%%%%%%%%
ETzr=(Czr>CT).*Fzr.*Bzr.*(1-exp(-0.45.*LAI)).*ETmax.*(1-b.*(Czr-CT))+(Czr<=CT).*Fzr.*Bzr.*(1-exp(-0.45.*LAI)).*ETmax  
ETzc=(Czc>CT).*Fzc.*Bzc.*(1-exp(-0.45.*LAI)).*ETmax.*(1-b.*(Czc-CT))+(Czc<=CT).*Fzc.*Bzc.*(1-exp(-0.45.*LAI)).*ETmax  
ETzg1=(Czg>CT).*(1-ETzr-ETzc)*(1-b.*(Czg-CT))+(Czg<=CT).*(1-ETzr-ETzc);      
ETzg2=(Czg>CT).*(Rzg.*(1-exp(-0.45.*LAI)).*ETmax)*(1-b.*(Czg-CT))+(Czg<=CT).*(Rzg.*(1-exp(-0.45.*LAI)).*ETmax);      
ETzg=min(ETzg1,ETzg2);                                                
YY=[ETzr,ETzc,ETzg]                                                          
YYY=YY./(ETzr+ETzc+ETzg)
AAA=YYY'
bar(YYY,'stacked')