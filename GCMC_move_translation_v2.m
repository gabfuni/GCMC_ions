%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [top, translation_done]=GCMC_move_translation_v2(top,depth,n,Temp,eps_H2O)
%% translation move called by GCMC_fun3.m
% [top, translation_done]=GCMC_move_translation_v2(top,depth,n)
% Input
% top: topology resulting before the insertion move
% n: width and length of the simulation box
% depth: depth of the simulation box
% Temp: Temperature
% eps_H2O: relative permittivity

% Output
% top: topology resulting after the translation move
% translation_done: positive if the translation move was succeful

translation_done=0;
kB=1.380946e-23; %boltzmann constant J/K
beta=1/(kB*Temp);

% insert molecule and find its neighbours

random_move_ion1=randi([1,length(top.ion1)]);
random_move_ion2=randi([1,length(top.ion2)]);

guess.pos_ion1=top.ion1(random_move_ion1).pos;
guess.d_ion1=top.ion1(random_move_ion1).d_steric;
guess.charge_ion1=top.ion1(random_move_ion1).charge;

guess.pos_ion2=top.ion2(random_move_ion2).pos;
guess.d_ion2=top.ion2(random_move_ion2).d_steric;
guess.charge_ion2=top.ion2(random_move_ion2).charge;

[guess] = GCMC_move_traslation_inside(guess,n,depth);

[guess, flag, d2] = search_neighbours (top,n,guess,random_move_ion1,random_move_ion2);

if flag==1 %non devo considerare nulla, le particelle sono troppo vicine
    flag=0;
else
    [guess,DU] = calculate_energy (guess,d2,eps_H2O);

    acc_probability=exp(-DU*beta);%probability of retaining the inserted particle
    
    R=rand;
    if acc_probability>R 
        top.ion1(random_move_ion1).pos=guess.pos_ion1;
        top.ion1(random_move_ion1).charge=guess.charge_ion1;
        top.ion1(random_move_ion1).d_steric=guess.d_ion1;
        top.ion1(random_move_ion1).acc_probability=acc_probability;
        top.ion1(random_move_ion1).DU=DU;
        top.ion2(random_move_ion2).pos=guess.pos_ion2;
        top.ion2(random_move_ion2).charge=guess.charge_ion2;
        top.ion2(random_move_ion2).d_steric=guess.d_ion2;
        top.ion2(random_move_ion2).acc_probability=acc_probability;
        top.ion2(random_move_ion2).DU=DU;
        translation_done=1;
    end
end
