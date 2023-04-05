%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%

function pp_flat_plate_v3_GCMC(top_tot,n,depth)
%% Plot the solution of the GCMC simulations in two figure
% pp_flat_plate_v3_GCMC(top_tot,n,depth)
% Input
% top_tot: vector with all the topologies simulated
% depth: depth of the simulation box
% n: length or width of the simulation box

% Output
% Figure 1: 3D plot of the topology
% Figure 2: Average distribution of the ions along the z direaction.
% Analytical solution plotted along with the GCMC results


NA=6.022E23;% Avogadro number

% binning
dz=0.1e-9;
binslength=[0:dz:depth];

surf=vertcat(top_tot{1,1}.surf([1:length(top_tot{1,1}.surf)]).pos);
ion1=vertcat(top_tot{1,1}.ion1([1:length(top_tot{1,1}.ion1)]).pos);
ion2=vertcat(top_tot{1,1}.ion2([1:length(top_tot{1,1}.ion2)]).pos);

%plot the position in 3D of the surface charge and of the ions in solution
%of the first calculated topology
figure1=figure;
coltab=[0 0 0.5;0.5 0 0; 0 0.5 0];

scatter3(surf(:,1),surf(:,2),surf(:,3),50,coltab(1,:),'filled')
hold on
scatter3(ion1(:,1),ion1(:,2),ion1(:,3),50,coltab(2,:),'filled')
scatter3(ion2(:,1),ion2(:,2),ion2(:,3),50,coltab(3,:),'filled')

xlabel('m')
ylabel('m')
zlabel('m')
set(gca, 'fontsize',20)

%%plot the distribution of counterions
ions_distribution=zeros(length(top_tot),length(binslength)-1);
ions_distribution2=ions_distribution;
j=1;
for i=1:length(top_tot)
    pos_ions1=vertcat(top_tot{1,i}.ion1([1:length(top_tot{1,i}.ion1)]).pos);
    ions_distribution(i,:)=histcounts(pos_ions1(:,3),[0:dz:depth]);
end

for i=1:length(top_tot)
    pos_ions2=vertcat(top_tot{1,i}.ion2([1:length(top_tot{1,i}.ion2)]).pos);
    ions_distribution2(i,:)=histcounts(pos_ions2(:,3),[0:dz:depth]);
end

z=[0:dz:depth];


min_ions=min(ions_distribution);
max_ions=max(ions_distribution);
mean_ions=mean(ions_distribution);
std_ions=std(ions_distribution);
min_ions2=min(ions_distribution2);
max_ions2=max(ions_distribution2);
mean_ions2=mean(ions_distribution2);
std_ions2=std(ions_distribution2);
z=z(2:end);

ions_distribution_plot=ions_distribution(1,:)/((n)^2*dz*NA);
ions_distribution2_plot=ions_distribution2(1,:)/((n)^2*dz*NA);

figure2=figure;
mean_ions=mean_ions/((n)^2*dz*NA);
std_ions=std_ions/((n)^2*dz*NA);
curve1=mean_ions(:)+std_ions(:);
curve2=mean_ions(:)-std_ions(:);
z2 = [z, fliplr(z)]';

inBetween = [curve1', fliplr(curve2')]';
h2=fill(z2, inBetween, [1 0.7 0.7]);
set(h2,'EdgeColor','none')
set(h2,'facealpha',.8)
hold on

mean_ions2=mean_ions2/((n)^2*dz*NA);
std_ions2=std_ions2/((n)^2*dz*NA);
curve1=mean_ions2(:)+std_ions2(:);
curve2=mean_ions2(:)-std_ions2(:);
z2 = [z, fliplr(z)]';

inBetween = [curve1', fliplr(curve2')]';
h2=fill(z2, inBetween, [0.7 0.7 1]);
set(h2,'EdgeColor','none')
set(h2,'facealpha',.8)
hold on
plot(z,mean_ions(:),'linewidth',1.5,'color',[0.5 0 0])
plot(z,mean_ions2(:),'linewidth',1.5,'color',[0 0 0.5])

xlabel("z [m]")
ylabel("# of charges [mol/m^3]")
set(gca,'fontsize',20)

