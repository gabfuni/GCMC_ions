%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function guess=initialize_guess(d_ion1,charge_ion1,d_ion2,charge_ion2)
%% initialize guess particles data
% guess=initialize_guess(d_ion1,charge_ion1,d_ion2,charge_ion2)

% Input
% d_ion1: diameter of the ion 1
% d_ion2: diameter of the ion 2
% charge_ion1: charge of ion1
% charge_ion2: charge of ion2
% B: exponential factor

% Output
% guess: inserted particles

guess.d_ion1=d_ion1;
guess.charge_ion1=charge_ion1;
guess.d_ion2=d_ion2;
guess.charge_ion2=charge_ion2;
end

