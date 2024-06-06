function [A2u,covRes,Mat_covRes,Xi_A2] = MDM_nullO_LTI(L,F,G,E,nz,H,D,z,u,Qb,Rb,version)
% Author: Oliver Kost, kost@ntis.zcu.cz
%
% Requires files: MDM_nullO.m, O_Gamma.m
%
% nullspace MDM with null(O)*Z
% L: number of measuremnts in Z; User parameter
%
% version = 0: null(O)*Z
% version = 1: null([O,Gamma*G])*[Z;U]
% version = 2: (eye-O*pinv(O))*Z
% version = 3: null(O,'r')*Z

if ~exist('version','var')
    version = 0; % Known input U
elseif isempty(version)
    version = 0; % Known input U
end

[nx,nw] = size(E);
nv = size(D,2);
nu = size(G,2);
Number = size(z,1);

F = mat2cell(repmat(F,L-1,1),repmat(nx,1,L-1),nx);
G = mat2cell(repmat(G,L-1,1),repmat(nx,1,L-1),nu);
E = mat2cell(repmat(E,L-1,1),repmat(nx,1,L-1),nw);

H = mat2cell(repmat(H,L,1),repmat(nz,1,L),nx);
D = mat2cell(repmat(D,L,1),repmat(nz,1,L),nv);
nz = repmat(nz,L,1);

[A2u,covRes_1,Mat_covRes,Xi_A2] = MDM_nullO(L,F,G,E,nz,H,D,z,u,Qb,Rb,version);

nz = nz(1);
z = horzcat(z{:}); % Measurement vector from cell to matrix
Z=zeros((L)*nz,Number-L+1); % Augmented measurement vector
for i=0:L-1
    Z((1+i*nz:(i+1)*nz),1:end)=z(:,(1+i:Number-L+1+i)); 
end

if version ~= 1
    u = horzcat(u{:}); % Control vector from cell to matrix
    U=zeros((L-1)*nu,Number-L+1); % Augmented input vector
    for i=0:L-1-1
        U((1+i*nu:(i+1)*nu),1:end)=u(:,(1+i:Number-L+1+i)); 
    end
    Awu = Mat_covRes{1}{1};
    Awu_blkG = Awu * blkdiag(G{1:1+L-2});
end

covRes=cell(Number-L+1,1); % Covariance of residue vector
covRes(1) = covRes_1; % First covariance of residue vector from MDM_nullO

Avz = Mat_covRes{2}{1};
Xi = Mat_covRes{3}{1};
nRes2 = size(Xi,2);

Res = Avz * Z;
if version ~= 1
    Res = Res - Awu_blkG * U; % Residue vector
end
for t=2:Number-L+1  
    covRes{t} = Xi*reshape(Res(:,t)*Res(:,t)',nRes2,1); % Covariance of residue vector
end

end

