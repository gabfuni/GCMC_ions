%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [guess, flag, d2] = search_neighbours (top,n,guess,number_move_ion1,number_move_ion2)
%% search neighbors close to the inserted particle
%[guess, flag, d2] = search_neighbours (top,n,guess,number_move_ion1,number_move_ion2)
% Input
% top: topology
% n: width and length of the simulation box
% guess: inserted particles
% number_move_ion1: used for the translation and delation modes, it avoids
% double counting of the guess molecule
% number_move_ion2: used for the translation and delation modes, it avoids
% double counting of the guess molecule


% Output
% guess: inserted particles
% flag: 0-1 value
% d2: distance between the inserted ion1 and ion2


flag=0;

[d2] = d_calculation_pbc(guess.pos_ion1, guess.pos_ion2,n);
if d2<(guess.d_ion1+guess.d_ion2)/2
    flag=1;
    return
end
%% surface charge
k=1; %counts the neighbours surf

% energy contribution: surface charges - ions
for j=1:length(top.surf) %search_neighbours

    [guess.n_csurf_ion(k,:),guess.n_charge_csurf_ion(k),guess.distance_csurf_ion(k),flag] =...
        fun_interactions_csurf_ion_pbc(guess.pos_ion1,top.surf(j).pos,top.surf(j).charge,top.surf(j).d_steric,guess.d_ion1,flag,n); 
    if flag==1
        break
    end

    [guess.n_csurf_ion2(k,:),guess.n_charge_csurf_ion2(k),guess.distance_csurf_ion2(k)] =...
        fun_interactions_csurf_ion_pbc(guess.pos_ion2,top.surf(j).pos,top.surf(j).charge,top.surf(j).d_steric,guess.d_ion2,flag,n); 
    if flag==1
        break;
    end
    k=k+1;
end

% energy contribution: ions - ions
if flag==0
    k1=1; %counts the neighbours ion1
    for j=1:length(top.ion1) %search_neighbour
        if j==number_move_ion1
            continue
        end

        [guess.n_ions_ion(k1,:),guess.n_charge_ions_ion(k1),guess.distance_ions_ion(k1),flag] =...
           fun_interactions_ions_ion_pbc(top.ion1(j).pos,guess.pos_ion1,top.ion1(j).charge,top.ion1(j).d_steric,guess.d_ion1,flag,n); 
        if flag==1
            break;
        end

        [guess.n_ions_ion2(k1,:),guess.n_charge_ions_ion2(k1),guess.distance_ions_ion2(k1),flag] =...
            fun_interactions_ions_ion_pbc(top.ion1(j).pos,guess.pos_ion2,top.ion1(j).charge,top.ion1(j).d_steric,guess.d_ion2,flag,n); 
        if flag==1
            break;
        end
        k1=k1+1;
    end
end


if flag==0
    k2=1; %counts the neighbours ion2
    for j=1:length(top.ion2) %search_neighbour
        if j==number_move_ion2
            continue
        end
        [guess.n_ions2_ion(k2,:),guess.n_charge_ions2_ion(k2),guess.distance_ions2_ion(k2),flag] =...
            fun_interactions_ions_ion_pbc(top.ion2(j).pos,guess.pos_ion1,top.ion2(j).charge,top.ion2(j).d_steric,guess.d_ion1,flag,n); 
        if flag==1
            break;
        end
        [guess.n_ions2_ion2(k2,:),guess.n_charge_ions2_ion2(k2),guess.distance_ions2_ion2(k2),flag] =...
            fun_interactions_ions_ion_pbc(top.ion2(j).pos,guess.pos_ion2,top.ion2(j).charge,top.ion2(j).d_steric,guess.d_ion2,flag,n); 
        if flag==1
            break;
        end
        k2=k2+1;
    end
end