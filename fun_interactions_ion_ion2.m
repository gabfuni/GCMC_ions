%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [charge_U_ion_ion2] = fun_interactions_ion_ion2(d,charge_ion,charge_ion2,eps_H2O) 
%% Potential energy between the inserted ion1 and ion2
% [charge_U_ion_ion2] = fun_interactions_ion_ion2(d,charge_ion,charge_ion2,eps_H2O) 
% Input:
% d: distance between ion 1 and ion 2
% charge_ion1: charge of ion1
% charge_ion2: charge of ion2
% eps_H2O: relative permittivity

% Output:
% charge_U_ion_ion2: interaction energy between the inserted ion1 and ion2

if ~isnan(d)
    charge_U_ion_ion2=charge_charge_v2(eps_H2O,charge_ion,charge_ion2,d); %interazione tra testa carica e counterion del surfattante inserito
else
    charge_U_ion_ion2=0;
end