function [O,Gamma] = O_Gamma(F,H,nz,L,t)
% Author: Oliver Kost, kost@ntis.zcu.cz
%
% Computation of "O" matrix (observability matrix) and the "Gamma" matrix for MDM

nx = size(F{1},1);

Gamma = zeros(sum(nz(t:t+L-1)),(L-1)*nx);
for j=1:L % column
    part_O_Gamma = zeros(sum(nz(t:t+L-1)),nx);
    partF = 1;
    for i=j:L % row
        if i==j
            part_O_Gamma((1:nz(t+i-1))+sum(nz(t:t+i-2)),:) = H{t+i-1};
        elseif i>j
            partF = F{t+i-2} * partF;
            part_O_Gamma((1:nz(t+i-1))+sum(nz(t:t+i-2)),:) = H{t+i-1} * partF;
        end
    end
    if j==1
        O = part_O_Gamma;
    else
        Gamma(:,(1:nx)+nx*(j-2)) = part_O_Gamma;
    end
end

end

