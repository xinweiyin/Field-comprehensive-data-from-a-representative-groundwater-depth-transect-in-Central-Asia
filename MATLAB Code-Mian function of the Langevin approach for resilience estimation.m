%%%This code is referenced from paper ¡°Arani, B. M. S., Carpenter, S. R., Lahti, L., van Nes, E. H. & Scheffer, M. Exit time as a measure of ecological resilience. Science 372, 1168 (2021)¡±. The auxiliary codes of this mian code is detailed in Arani et al. (17).
function results = allfigures
addpath('chebfun-master');
results = [];
realmod = [];
asigma = 0.3; 
if plotfigno(1) %Figure 1 extended abstract
adjustplot(gcf, '-d')
if true %can be set to true to generate a new data set
%needs grind for matlab to be installed: download from https://www.sparcs-center.org/grind
% global sigma g_t g_Y
% use may_babak
asigma = 0.3
sigma = asigma;
t = csvread();
x = csvread();
data = table(t, x);
writetable(data, 'may_long.csv')
end
%load fig_main mod
nreplicates=1000 
tab = readtable('may_long.csv') 
L=min(x);
R=max(x);
dt=1;
results = LangevinReconst(tab.x, L, R, 60, 1:5, dt, 'Nadaraya-Watson')
if nreplicates > 0
results = Langevin_bootstrap(tab.x, results, {'error-propagation', 'block'}, nreplicates, struct('MaxStep', 0.001)); 
else
results = results;
end
save('Fig2main_boot1', 'results' , '-v7.3');
else
%load fig_main mod
res = load('Fig2main_boot1');
results = res.results;
results.bootstrap = res.results.bootstrap(2);
mod = langevin_eq(results)        
%clc;
asigma = 0.3 
realmod = openrealmod(asigma, mod.domain)
mod.namex = 'State';        
hax = subplot(6, 1, 1);
mod.plot('potential', 'verticalbar', false);
        
tab_short = readtable('may_long.csv');
hax = subplot(6, 1, 2);
plot(tab_short.x, tab_short.t, 'b-');
realmod = openrealmod(asigma, mod.domain);

hax = subplot(6, 1, 3);
adjustplot('-d')
mod.plot('D1', 'error', true)

hax = subplot(6, 1, 4);
adjustplot('-d')
mod.plot('D2', 'error', true)

hax = subplot(6, 1, 5);
thepdf = mod.pdf;
mod.plot(thepdf);

hax = subplot(6, 1, 6);
mod.plot('mean_exit');
  
end
end
