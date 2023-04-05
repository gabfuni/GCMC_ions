%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [top] = initialize_topology_GCMC_v2(num_surf,depth,n,charge_surf,charge_ion1,charge_ion2,d_ion1,d_ion2)
%% initialize the topology
%[top] = initialize_topology_GCMC_v2(num_surf,depth,n,charge_surf,charge_ion1,charge_ion2,d_ion1,d_ion2)
% num_surf: number of charges at the surface
% n: width and length of the simulation box
% depth: depth of the simulation box
% charge_surf: charge of every charge at the surface (+1)
% charge_ion1: charge of ion1
% charge_ion2: charge of ion2
% d_ion1: diameter of the ion 1
% d_ion2: diameter of the ion 2

% Output:
% top: topology with num_surf of charges at the surface, 1 ion 1 and num_surf+1 ions 2

if num_surf==0 
    disp("error")
    return
else
    top.ion1(1).pos=[n*rand(1,2) depth-d_ion1/2];
    top.ion1(1).charge=charge_ion1;
    top.ion1(1).d_steric=d_ion1;
    for i=1:num_surf
        top.surf(i).pos=[n*rand(1,2) depth];
        top.surf(i).charge=charge_surf;
        top.surf(i).d_steric=0;
        top.ion2(i).pos=[n*rand(1,2) rand*(depth-d_ion2/2)];
        top.ion2(i).charge=charge_ion2;
        top.ion2(i).d_steric=d_ion2;
    end
    top.ion2(i+1).pos=[n*rand(1,2) rand*(depth-d_ion2/2)];
    top.ion2(i+1).charge=charge_ion2;
    top.ion2(i+1).d_steric=d_ion2;
end
