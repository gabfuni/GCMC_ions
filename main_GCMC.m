%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear 
clc

%% constants
NA=6.022E23;% Avogadro number
q_el=1.60217662e-19; %elementary charge C
eps0=8.85418781762e-12; %permettivity of vacuum F/m
kB=1.380946e-23; %boltzmann constant J/K

eps_H2O=80; %relative permittivity

%% input values
gamma_pm=1; %activity coefficient
Temp=298; %Temperature K
ro_surf=1.044e-07;% surface concentration mol/m2  
rho_salt=10;% concentration of salt in water mol/m3
n=30e-9; %m box dimensions
depth=15E-9; %m box debth
d_ion1=0.36e-9;% diamater of the ion 1 Na+  m
d_ion2=0.33e-9;% diameter of the ion 2Cl-   m
nSteps=50000; % number of steps 


%% 
sigma=ro_surf*NA*q_el; % surface charge C/m2
num_surf=round(ro_surf*n^2*NA); %number of particles on the surface
ki=((2*rho_salt*NA*q_el^2)/(eps_H2O*eps0*kB*Temp))^0.5; %inverse of the Debye length
%%Grahame equation
psi0=(2*kB*Temp/q_el)*asinh(sigma/(8*rho_salt*NA*eps_H2O*eps0*kB*Temp)^0.5);

%% main
for j=1:1 %in series
    tic
    i=1;
    parfor i=1:8 %in parallel
        top=[];
        seed_MC=100*j+i; % seed number
        rng(seed_MC,'twister'); % For reproducibility
        
        %GCMC function
        [top,number_insertion,number_deletion,number_translation] = GCMC_fun3(top,num_surf,rho_salt, n,depth, Temp,eps_H2O,nSteps,gamma_pm,strcat(num2str(i),'_',num2str(j)),d_ion1,d_ion2)

        top_tot{i,j}=top;
        %number of moves performed
        deletion_done_tot{i,j}=number_deletion;
        insertion_done_tot{i,j}=number_insertion;
        translation_done_tot{i,j}=number_translation;
    end
    % duration of the simulation
    GCMC_duration(j)=toc
end
top_tot=top_tot(:)';
pp_flat_plate_v3_GCMC(top_tot,n,depth);  %post processing of the GCMC results
pp_analitical(ro_surf,rho_salt,depth,Temp,eps_H2O);  %analytical solution
save('outputs/output'); % save the simulation files
