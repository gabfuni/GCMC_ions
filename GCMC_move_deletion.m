%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [top, deletion_done]=GCMC_move_deletion(top,n,B, Temp, eps_H2O)
%% insertion move called by GCMC_fun3.m
% [top, deletion_done]=GCMC_move_deletion(top,n,B, Temp, eps_H2O)
% Input
% top: topology resulting before the insertion move
% n: width and length of the simulation box
% Temp: Temperature
% eps_H2O: relative permittivity
% B: exponential factor

% Output
% top: topology resulting after the insertion move
% deletion_done: positive if the deletion move was succeful

deletion_done=0;

kB=1.380946e-23; %boltzmann constant J/K
beta=1/(kB*Temp);

random_delete_ion1=randi([1,length(top.ion1)]);
random_delete_ion2=randi([1,length(top.ion2)]);

% insert molecule and calculate its neighbours
guess.pos_ion1=top.ion1(random_delete_ion1).pos;
guess.d_ion1=top.ion1(random_delete_ion1).d_steric;
guess.charge_ion1=top.ion1(random_delete_ion1).charge;

guess.pos_ion2=top.ion2(random_delete_ion2).pos;
guess.d_ion2=top.ion2(random_delete_ion2).d_steric;
guess.charge_ion2=top.ion2(random_delete_ion2).charge;

[guess, flag, d2] = search_neighbours (top,n,guess,random_delete_ion1,random_delete_ion2);

if flag==1 
    flag=0;
else
    [guess, DU] = calculate_energy (guess,d2,eps_H2O);
    
    acc_probability=1/(exp(B-DU*beta)/((length(top.ion1)+1)*(length(top.ion2)+1)));
    
    R=rand;
    if acc_probability>R 
        top.ion1(random_delete_ion1)=[];
        top.ion2(random_delete_ion2)=[];
        deletion_done=1;
    end
end
