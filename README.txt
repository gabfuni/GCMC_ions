%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gran Canonical Monte Carlo algorithm to calculate the distribution of ions
close to a charged planar surface

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FOLDERS:
- outputs: contains the .mat files with the topology obtained at the end of the simulation and the partial outputs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTIONS:
- calculate_energy.m: calculates the energy variation due to the GCMC move
- charge_charge_v2.m: calculate the charge-charge interaction energy
- d_calculation_pbc.m: calculate the distance between two particles 
- fun_interactions_csurf_ion_pbc.m: collects the distances and features of the esisting surface charges and the inserted ion
- fun_interactions_ion_ion2.m: potential energy between the inserted ion1 and ion2
- fun_interactions_ions_ion_pbc.m: collects the distances and features of the esisting ions and the inserted ion
- GCMC_fun3.m: solve the GCMC algorithm for a mixture of monovalent ions close to planar surface
- GCMC_move_deletion.m: insertion move called by GCMC_fun3.m
- GCMC_move_insertion.m: insertion move called by GCMC_fun3.m
- GCMC_move_translation_v2.m: translation move called by GCMC_fun3.m
- GCMC_move_traslation_inside.m: traslates a molecule, considers pbc and box dimensions
- initialize_guess.m: initialize guess particles data
- initialize_topology_GCMC_v2.m: initialize the topology
- main_GCMC.m: main MATLAB file to be run
- move_traslation_inside.m: translate a particle
- pp_analitical.m: plot the analyitcal solution of the Poisson-Boltzmann equation
- pp_flat_plate_v3_GCMC.m: plot the solution of the GCMC simulations in two figure
- search_neighbours.m: search neighbors close to the inserted particle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WORKFLOW:
1) open the main program
- Open main_GCMC.m and choose the x-y dimensions of the simualtion box
 (variable n) and the depth, the dielectric constant (eps_H2O, default setting 80),
 the diameter of the ions and the surface concentration of the charges
- Vary the j from 1 to the number of simulations that you wnat to run in series
- Vary the i from 1 to the number of simulations that you wnat to run in parallel
- Run the code :)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%