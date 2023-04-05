%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [n_ions_ion,n_charge_ions_ion,distance_ions_ion,flag] = fun_interactions_ions_ion_pbc(pos_ion,pos_ion_guess,charge_ion,d_ion,d_ion_guess,flag,n) 
%% collects the distances and features of the esisting ions and the inserted ion
% [n_ions_ion,n_charge_ions_ion,distance_ions_ion,flag] = fun_interactions_ions_ion_pbc(pos_ion,pos_ion_guess,charge_ion,d_ion,d_ion_guess,flag,n)
% Input:
% pos_ion: position of the existing ion
% pos_ions_guess: position of the inserted ion
% d: distance between ion 1 and ion 2
% charge_ion1: charge of ion
% d_ion_diameter of the existing ion
% d_ion_guess: diameter of the inserted ion
% flag: 0-1 value
% n: length or width of the simulation box

% Output:
% n_ions_ion: ions close to the inserted ion
% n_charge_ions_ion: charge of the ions close to the inserted ion
% distance_ions_ion: distance of the ions from the inserted ion
% flag: 0-1 value


[d] = d_calculation_pbc(pos_ion, pos_ion_guess,n);
if ~isnan(d)
    if d<(d_ion_guess+d_ion)/2
        flag=1;
    end
    n_ions_ion=pos_ion;
    n_charge_ions_ion=charge_ion; 
    distance_ions_ion=d;
else
    n_ions_ion=[NaN NaN NaN];
    n_charge_ions_ion=NaN; 
    distance_ions_ion=NaN;
end