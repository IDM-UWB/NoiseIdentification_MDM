function covRes = MDM_covRes_LTI(L,G,z,u,Mat_covRes,version)
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

Avz = Mat_covRes{2}{1};
Xi = Mat_covRes{3}{1};
nRes2 = size(Avz,1)^2;

nz = size(z{1},1);
z = horzcat(z{:}); % Measurement vector from cell to matrix
Z=zeros(L*nz,Number-L+1); % Augmented measurement vector
for i=0:L-1
    Z((1+i*nz:(i+1)*nz),1:end)=z(:,(1+i:Number-L+1+i)); 
end

if version ~= 1
    nu = size(G,2);
    u = horzcat(u{:}); % Input vector from cell to matrix
    U=zeros((L-1)*nu,Number-L+1); % Augmented input vector
    for i=0:L-1-1
        U((1+i*nu:(i+1)*nu),1:end)=u(:,(1+i:Number-L+1+i)); 
    end
    Awu = Mat_covRes{1}{1};
    Awu_blkG = Awu * kron(eye(L-1),G);
end

Res = Avz * Z; 
if version ~= 1
    Res = Res - Awu_blkG * U; 
end

covRes=cell(Number-L+1,1); % Covariance of residue vector
for t=1:Number-L+1
    covRes{t} = Xi*reshape(Res(:,t)*Res(:,t)',nRes2,1);
end

end

