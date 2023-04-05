%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [guess, DU] = calculate_energy (guess,d2,eps_H2O)
%% calculate energy variation due to the GCMC move
% [guess,guessU, DU, DU_charge] = calculate_energy (guess,d2,eps_H2O)
% Input:
% guess: inserted particle(s), consider its neighbors
% d2: distance between the inserted ion1 and ion2
% eps_H2O: relative permittivity

% Output:
% guess: inserted particle(s), consider its neighbors
% guessU: energy associated to the particle insertion/deletion or
% traslation
% DU: total energy of the system associated to the GCMC move (in this case
% it is equal to DU_charge)


%variation in the system elecric potential due to the new charge
DU_charge=0; 

% contributoin from surface charges
for i_U=1:length(guess.n_csurf_ion(:,1)) 
    %charge-charge interactions U
    guessU.charge_U_csurf_ion(i_U)=charge_charge_v2(eps_H2O, guess.charge_ion1, guess.n_charge_csurf_ion(i_U), guess.distance_csurf_ion(i_U));
    guessU.charge_U_csurf_ion2(i_U)=charge_charge_v2(eps_H2O, guess.charge_ion2, guess.n_charge_csurf_ion2(i_U), guess.distance_csurf_ion2(i_U));
    DU_charge=DU_charge+guessU.charge_U_csurf_ion(i_U)+guessU.charge_U_csurf_ion2(i_U);
end
% contribution from ions 1
for i_U=1:length(guess.n_ions_ion(:,1))
    guessU.charge_U_ions_ion(i_U)=charge_charge_v2(eps_H2O, guess.charge_ion1, guess.n_charge_ions_ion(i_U), guess.distance_ions_ion(i_U)); 
    guessU.charge_U_ions_ion2(i_U)=charge_charge_v2(eps_H2O, guess.charge_ion2, guess.n_charge_ions_ion2(i_U), guess.distance_ions_ion2(i_U)); 
    DU_charge=DU_charge+guessU.charge_U_ions_ion(i_U)+guessU.charge_U_ions_ion2(i_U);
end
% contribution from ions 2
for i_U=1:length(guess.n_ions2_ion(:,1))
    guessU.charge_U_ions2_ion(i_U)=charge_charge_v2(eps_H2O, guess.charge_ion1, guess.n_charge_ions2_ion(i_U), guess.distance_ions2_ion(i_U));
    guessU.charge_U_ions2_ion2(i_U)=charge_charge_v2(eps_H2O, guess.charge_ion2, guess.n_charge_ions2_ion2(i_U), guess.distance_ions2_ion2(i_U)); 
    DU_charge=DU_charge+guessU.charge_U_ions2_ion(i_U)+guessU.charge_U_ions2_ion2(i_U);
end

% interactions between the two inserted ions
[guessU.charge_U_ion_ion2] = fun_interactions_ion_ion2(d2,guess.charge_ion1,guess.charge_ion2,eps_H2O); 

% total energy
DU=DU_charge+guessU.charge_U_ion_ion2;
end
