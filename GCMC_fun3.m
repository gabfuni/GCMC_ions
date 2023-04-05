%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [top,number_insertion,number_deletion,number_translation] = GCMC_fun3(top,num_surf,rho_salt, n,depth, Temp,eps_H2O,nSteps,gamma_pm,name_save, d_ion1,d_ion2)
%% solve the GCMC algorithm for a mixture of monovalent ions close to planar surface
% [top,number_insertion,number_deletion,number_translation] = GCMC_fun3(top,num_surf,rho_salt, n,depth, Temp,eps,nSteps,gamma_pm,name_save, d_ion1,d_ion2)
% Input
% top: initial topology (empty vector if the simulation just started)
% num_surf: number of charges at the surface
% rho_salt: concentration of the salt in solution
% n: width and length of the simulation box
% depth: depth of the simulation box
% Temp: Temperature
% eps_H2O: relative permittivity
% nSteps: number of steps of the GCMC algorithm
% gamma_pm: activity coefficient
% name_save: name of the partial outputs to be saved during the simulation
% d_ion1: diameter of the ion 1
% d_ion2: diameter of the ion 2

% Output
% top: topology resulting from the GCMC algorithm
% number_insertion: number of insertions done
% number_deletion: number of deletions done
% number_translation: number of translations done


NA=6.022E23;% Avogadro number
Vol=n^2*depth;%volume of the simulation box
rho_target=rho_salt*NA;%number denisty
B=2*(log(gamma_pm)+log(Vol*rho_target));

q_el=1.60217662e-19; %elementary charge C

charge_surf=+1*q_el; %charge of every charge at the surface (+1)
charge_ion1=charge_surf; %charge of ion1(monovalent) (+1)
charge_ion2=-charge_ion1; %charge of ion2(monovalent) (-1)

%number of moves performed
number_insertion=0;
number_deletion=0;
number_translation=0;

%topology creation
if isempty(top)==1
    [top] = initialize_topology_GCMC_v2(num_surf,depth,n,charge_surf,charge_ion1,charge_ion2,d_ion1,d_ion2);
end

for iii=1:nSteps
    %save partial outputs in /outputs
    if iii==nSteps/20
        save(strcat('outputs/top_',num2str(iii),'_GCMC_fun3_',name_save));
    end
    if iii==10000
        save(strcat('outputs/top_',num2str(iii),'_GCMC_fun3_',name_save));
    end
    if iii==20000
        save(strcat('outputs/top_',num2str(iii),'_GCMC_fun3_',name_save));
    end
    if iii==50000
        save(strcat('outputs/top_',num2str(iii),'_GCMC_fun3_',name_save));
    end
    if iii==75000
        save(strcat('outputs/top_',num2str(iii),'_GCMC_fun3_',name_save));
    end
    if iii==100000
        save(strcat('outputs/top_',num2str(iii),'_GCMC_fun3_',name_save));
    end
    if iii==150000
        save(strcat('outputs/top_',num2str(iii),'_GCMC_fun3_',name_save));
    end
    gcmc_move=rand;
    
    if (gcmc_move > 0.4 && gcmc_move<=0.8)
        %% insertion
        [top, insertion_done]=GCMC_move_insertion(top,depth,n,d_ion1,charge_ion1,d_ion2,charge_ion2,B, Temp, eps_H2O);
        number_insertion=number_insertion+insertion_done;
    elseif gcmc_move<=0.4 && length(top.ion1)>2
        %% deletion
        [top, deletion_done]=GCMC_move_deletion(top,n,B, Temp, eps_H2O);
        number_deletion=number_deletion+deletion_done;
    elseif gcmc_move>0.8 && length(top.ion1)>2
        %% translation
        [top, translation_done]=GCMC_move_translation_v2(top,depth,n, Temp, eps_H2O);
        number_translation=number_translation+translation_done;
    else 
        continue
    end
end
