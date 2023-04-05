%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [guess_new] = GCMC_move_traslation_inside(guess,n,depth)
%% traslates a molecule, considers pbc and box dimensions
% [guess_new] = GCMC_move_traslation_inside(guess,n,depth)
% Input
% guess: inserted particles
% n: width and length of the simulation box
% depth: depth of the simulation box

% Output
% guess_new: inserted particles updated

rho=0.3e-9; %radius in meters how far is translated the particle

guess_new=guess;

theta=deg2rad(randi([0,360]));
phi=deg2rad(randi([0,180])); 
guess_new.pos_ion1(1)=guess.pos_ion1(1)+rho*cos(theta)*sin(phi);
guess_new.pos_ion1(2)=guess.pos_ion1(2)+rho*sin(theta)*sin(phi);
guess_new.pos_ion1(3)=guess.pos_ion1(3)+rho*cos(phi);

theta=deg2rad(randi([0,360]));
phi=deg2rad(randi([0,180])); 
guess_new.pos_ion2(1)=guess.pos_ion2(1)+rho*cos(theta)*sin(phi);
guess_new.pos_ion2(2)=guess.pos_ion2(2)+rho*sin(theta)*sin(phi);
guess_new.pos_ion2(3)=guess.pos_ion2(3)+rho*cos(phi);

% z coordinate inside the box
while guess_new.pos_ion1(3)>(depth-guess_new.d_ion1/2) || guess_new.pos_ion1(3)<guess_new.d_ion1/2 || guess_new.pos_ion1(2)>n || guess_new.pos_ion1(2)<0|| guess_new.pos_ion1(1)>n || guess_new.pos_ion1(1)<0
    theta=deg2rad(randi([0,360]));
    phi=deg2rad(randi([0,180])); 
    guess_new.pos_ion1(1)=guess.pos_ion1(1)+rho*cos(theta)*sin(phi);
    guess_new.pos_ion1(2)=guess.pos_ion1(2)+rho*sin(theta)*sin(phi);
    guess_new.pos_ion1(3)=guess.pos_ion1(3)+rho*cos(phi);
end
while guess_new.pos_ion2(3)>(depth-guess_new.d_ion2/2) || guess_new.pos_ion2(3)<guess_new.d_ion2/2 || guess_new.pos_ion2(2)>n || guess_new.pos_ion2(2)<0|| guess_new.pos_ion2(1)>n || guess_new.pos_ion2(1)<0
    theta=deg2rad(randi([0,360]));
    phi=deg2rad(randi([0,180])); 
    guess_new.pos_ion2(1)=guess.pos_ion2(1)+rho*cos(theta)*sin(phi);
    guess_new.pos_ion2(2)=guess.pos_ion2(2)+rho*sin(theta)*sin(phi);
    guess_new.pos_ion2(3)=guess.pos_ion2(3)+rho*cos(phi);
end

end