%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [n_csurf_ion,n_charge_csurf_ion,distance_csurf_ion,flag] = fun_interactions_csurf_ion_pbc(pos_ion,pos,charge,d_steric,d_ion,flag,n) 
%% collects the distances and features of the esisting surface charges and the inserted ion
% [n_csurf_ion,n_charge_csurf_ion,distance_csurf_ion,flag] = fun_interactions_csurf_ion_pbc(pos_ion,pos,charge,d_steric,d_ion,flag,n)
% Input:
% pos_ion: position of the ion
% pos: position of the surface charge
% charge: charge of the surface charge
% d_steric: diameter of the surface charge (set to 0)
% d_ion: diameter of the ion
% flag: 0-1 value
% n: length or width of the simulation box

% Output:
% n_csurf_ion: position of the surface charge close to the ion
% n_charge_csurf_ion: charge of the surface charge close to the ion
% distance_csurf_ion: distance between the surface charge and the ion
% flag: 0-1 value

[d] = d_calculation_pbc(pos, pos_ion,n);
if ~isnan(d) && ~isnan(charge) 
    if d<(d_steric+d_ion)/2 %position of the ion not acceptable
        flag=1;
        n_csurf_ion=[NaN NaN NaN];
        n_charge_csurf_ion=NaN; % surface charges close to the ion 
        distance_csurf_ion=NaN;
        return
    end
    n_csurf_ion=pos;
    n_charge_csurf_ion=charge; % surface charges close to the ion 
    distance_csurf_ion=d;
else
    n_csurf_ion=[NaN NaN NaN];
    n_charge_csurf_ion=NaN; % surface charges close to the ion
    distance_csurf_ion=NaN;
end