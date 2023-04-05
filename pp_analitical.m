%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pp_analitical(rho_surf,c0,depth,Temp,eps_H2O)
%% Plot the analyitcal solution of the Poisson-Boltzmann equation
% pp_analitical(rho_surf,c0,depth,Temp,eps_H2O)
% Input
% rho_surf: concentration of the surface charges
% depth: depth of the simulation box
% c0: concentration of the electrolyte in solution
% Temp: temperature
% eps_H2O: relative permittivity


NA=6.022E23;% Avogadro number
q_el=1.60217662e-19; %elementary charge C
eps0=8.85418781762e-12; %permittivity of vacuum F/m
kB=1.380946e-23; %boltzmann constant J/K


%% analytical solution (Poisson-Boltzman for 1:1 mixtures)

sigma=rho_surf*NA*q_el; %C/m2
c0_molec=c0*NA;
ki=((2*c0_molec*q_el^2)/(eps_H2O*eps0*kB*Temp))^0.5;

%%Grahame equation
psi0=(2*kB*Temp/q_el)*asinh(sigma/(8*c0_molec*eps_H2O*eps0*kB*Temp)^0.5);
distance=linspace(0,depth);

y0=(q_el*psi0)/(kB*Temp);
A=exp(y0/2)+1;
B=exp(y0/2)-1;
D=2*kB*Temp/q_el;

psi_completo= @(x) D*log((A+B*exp(-ki*x))./(A-B*exp(-ki*x)));

fun_c=@(potential) c0*exp((-q_el*potential)/(kB*Temp));

plot(fliplr(distance),fun_c(psi_completo(distance)),'linewidth',2)
hold on
fun_c1=@(potential) c0*exp((q_el*potential)/(kB*Temp));

plot(fliplr(distance),fun_c1(psi_completo(distance)),'linewidth',2)
ylabel('Concentration [mol/m^3]')

yyaxis right

plot(fliplr(distance),psi_completo(distance),'linewidth',2)

xlabel('z [m]')

ylabel('\psi [V]')
set(gca,'fontsize',20)
