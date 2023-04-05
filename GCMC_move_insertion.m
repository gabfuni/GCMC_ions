%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [top, insertion_done]=GCMC_move_insertion(top,depth,n,d_ion1,charge_ion1,d_ion2,charge_ion2,B, Temp, eps_H2O)
%% insertion move called by GCMC_fun3.m
% [top, insertion_done]=GCMC_move_insertion(top,depth,n,d_ion1,charge_ion1,d_ion2,charge_ion2,B, Temp, eps_H2O)
% Input
% top: topology resulting before the insertion move
% n: width and length of the simulation box
% depth: depth of the simulation box
% Temp: Temperature
% eps_H2O: relative permittivity
% d_ion1: diameter of the ion 1
% d_ion2: diameter of the ion 2
% charge_ion1: charge of ion1
% charge_ion2: charge of ion2
% B: exponential factor

% Output
% top: topology resulting after the insertion move
% insertion_done: positive if the insertion move was succeful


insertion_done=0;
kB=1.380946e-23; %boltzmann constant J/K
beta=1/(kB*Temp);

% insert molecule and calculate its neighbours
guess=initialize_guess(d_ion1,charge_ion1,d_ion2,charge_ion2);
guess.pos_ion1=[n*rand(1,2) (depth-d_ion1/2)*rand(1,1)]; %random counterion in the bulk 
guess.pos_ion2=[n*rand(1,2) (depth-d_ion2/2)*rand(1,1)]; %random counterion in the bulk 

number_move_ion1=NaN; %used for the translation and delation modes, avoid double counting of the guess molecule
number_move_ion2=NaN; %used for the translation and delation modes, avoid double counting of the guess molecule

% searching neighbours of the inserted particle
[guess, flag, d2] = search_neighbours (top,n,guess,number_move_ion1,number_move_ion2);

if flag==1 %skip if particles too close
    flag=0;
else
    %energy calculation
    [guess, DU] = calculate_energy (guess,d2,eps_H2O);
    
    %probability of retaining the inserted particle
    acc_probability=exp(B-DU*beta)/((length(top.ion1)+1)*(length(top.ion2)+1));

    R=rand;
    if acc_probability>R 
        length_previous_ion1=length(top.ion1)+1;
        length_previous_ion2=length(top.ion2)+1;
        top.ion1(length_previous_ion1).pos=guess.pos_ion1;
        top.ion1(length_previous_ion1).charge=guess.charge_ion1;
        top.ion1(length_previous_ion1).d_steric=guess.d_ion1;
        top.ion1(length_previous_ion1).acc_probability=acc_probability;
        top.ion1(length_previous_ion1).DU=DU;
        top.ion2(length_previous_ion2).pos=guess.pos_ion2;
        top.ion2(length_previous_ion2).charge=guess.charge_ion2;
        top.ion2(length_previous_ion2).d_steric=guess.d_ion2;
        top.ion2(length_previous_ion2).acc_probability=acc_probability;
        top.ion2(length_previous_ion2).DU=DU;
        insertion_done=1;
    end
end
