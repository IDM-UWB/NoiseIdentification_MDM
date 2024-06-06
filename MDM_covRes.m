function covRes = MDM_covRes(L,G,z,u,Mat_covRes,version)
% Author: Oliver Kost, kost@ntis.zcu.cz
%
% Nullspace MDM with null(O)*Z
% L: number of measuremnts in Z; User parameter

if ~exist('version','var')
    version = 0; % Known input U
elseif isempty(version)
    version = 0; % Known input U
end

Number = size(z,1);

Awu = Mat_covRes{1};
Avz = Mat_covRes{2};
Xi = Mat_covRes{3};

covRes=cell(Number-L+1,1);
for t=1:Number-L+1
    Z = vertcat(z{t:(t+L-1)}); % Augmented measurement vector
    Res = Avz{t} * Z; 
    if version ~= 1
        U = vertcat(u{t:(t+L-2)}); % Augmented input vector
        Res = Res - Awu{t} * blkdiag(G{t:t+L-2}) * U;
    end
    covRes{t} = Xi{t}*reshape(Res*Res',size(Res,1)^2,1);
end

end





