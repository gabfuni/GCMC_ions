%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      An overview on modelling approaches for photochemical       %%%%
%%%% and photoelectrochemical solar fuels processes and technologies  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gabriele Falciani, Eliodoro Chiavazzo (eliodoro.chiavazzo@polito.it) %%
%%      Department of Energy, Politecnico di Torino, Turin, Italy       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r] = d_calculation_pbc(top_pos,guess_pos,n)
%% distance between two particles 
% [r] = d_calculation_pbc(top_pos,guess_pos,n)
% Input:
% top_pos: alrady existing particle position
% guess_pos: guess particle position
% n: width and length of the simulation box

% Output:
% r: distance between the particles

if top_pos==guess_pos
    r=NaN;
    return 
end

if isnan(top_pos(1)) || isnan(guess_pos(1))
    r=NaN;
else    
    X = [top_pos(1); guess_pos(1)];
    Y = [top_pos(2); guess_pos(2)];
    
    dX = pdist(X,'euclidean');
    if round(dX/n)==1
        if guess_pos(1)>n/2
            guess_pos(1)=guess_pos(1)-n;
        else
            guess_pos(1)=guess_pos(1)+n;
        end
    end

    dY = pdist(Y,'euclidean');
    if round(dY/n)==1
        if guess_pos(2)>n/2
            guess_pos(2)=guess_pos(2)-n;
        else
            guess_pos(2)=guess_pos(2)+n;
        end
    end
    
    R = [top_pos; guess_pos];
    r = pdist(R,'euclidean');
    
end
end