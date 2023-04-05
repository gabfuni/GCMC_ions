function [guess_new] = move_traslation_inside(guess,n,depth)
%% Translate a particle 
% [guess_new] = move_traslation_inside(guess,n,depth)
% Input:
% guess: particle to be translated 
% n: width and length of the simulation box
% depth: depth of the simulation box

% Output:
% guess: translated particle

rho=0.5e-9; %radius in meters for the translation of the particle

teta=deg2rad(randi([0,360]));
phi=deg2rad(randi([0,180])); 
guess_new.pos_ion1(1)=guess.pos_ion1(1)+rho*cos(teta)*sin(phi);
guess_new.pos_ion1(2)=guess.pos_ion1(2)+rho*sin(teta)*sin(phi);
guess_new.pos_ion1(3)=guess.pos_ion1(3)+rho*cos(phi);

teta=deg2rad(randi([0,360]));
phi=deg2rad(randi([0,180])); 
guess_new.pos_ion2(1)=guess.pos_ion2(1)+rho*cos(teta)*sin(phi);
guess_new.pos_ion2(2)=guess.pos_ion2(2)+rho*sin(teta)*sin(phi);
guess_new.pos_ion2(3)=guess.pos_ion2(3)+rho*cos(phi);

% z coordinate inside the box
while guess_new.pos_ion1(3)>depth || guess_new.pos_ion1(3)<0 || guess_new.pos_ion1(2)>n || guess_new.pos_ion1(2)<0|| guess_new.pos_ion1(1)>n || guess_new.pos_ion1(1)<0
    teta=deg2rad(randi([0,360]));
    phi=deg2rad(randi([0,180])); 
    guess_new.pos_ion1(1)=guess.pos_ion1(1)+rho*cos(teta)*sin(phi);
    guess_new.pos_ion1(2)=guess.pos_ion1(2)+rho*sin(teta)*sin(phi);
    guess_new.pos_ion1(3)=guess.pos_ion1(3)+rho*cos(phi);
end
while guess_new.pos_ion2(3)>depth || guess_new.pos_ion2(3)<0 || guess_new.pos_ion2(2)>n || guess_new.pos_ion2(2)<0|| guess_new.pos_ion2(1)>n || guess_new.pos_ion2(1)<0
    teta=deg2rad(randi([0,360]));
    phi=deg2rad(randi([0,180])); 
    guess_new.pos_ion2(1)=guess.pos_ion2(1)+rho*cos(teta)*sin(phi);
    guess_new.pos_ion2(2)=guess.pos_ion2(2)+rho*sin(teta)*sin(phi);
    guess_new.pos_ion2(3)=guess.pos_ion2(3)+rho*cos(phi);
end

end