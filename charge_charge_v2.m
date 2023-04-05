%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Uz] = charge_charge_v2(eps_H2O,Q1, Q2, d)
%% charge-charge interaction
% [U] = charge_charge_v2(eps,Q1, Q2, d)
% Input:
% eps: dielectric constant
% Q1: charge 1
% Q2: charge 2
% d: distance betweeen the two point charges

% Output:
% U: Interaction energy


if isnan(Q1) || isnan(Q2) || isnan(d)
    Uz=0;
else
    eps0=8.85418781762e-12; %F/m
	
    U=@(z) Q1*Q2/(4*pi*eps_H2O*eps0*z);

    Uz=U(d);
end