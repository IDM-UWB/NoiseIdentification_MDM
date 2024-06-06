function [A2u,covRes,Mat_covRes,Xi_A2] = MDM_nullO(L,F,G,E,nz,H,D,z,u,Qb,Rb,version)
% Author: Oliver Kost, kost@ntis.zcu.cz
%
% Requires files: O_Gamma.m
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

[nx,nw] = size(E{1});
nv = size(D{1},2);
Number = size(nz,1);

nQb = size(Qb,2);
nRb = size(Rb,2);

%%% Replikation matrix for noise covarainces
Psi=zeros(((L-1)*nw+L*nv)^2,nQb+nRb);
for j=1:nQb
    Psi(:,j)=reshape(blkdiag(kron(eye(L-1),Qb{j}),zeros(L*nv)),((L-1)*nw+L*nv)^2,1);
end
for j=1:nRb
    Psi(:,nQb+j)=reshape(blkdiag(zeros((L-1)*nw),kron(eye(L),Rb{j})),((L-1)*nw+L*nv)^2,1);
end
%%% End: Replikation matrix for noise covarainces

O = cell(Number-L+1,1); % Obsevable matrix
Gamma = cell(Number-L+1,1);% Gamma matrix

Awu = cell(Number-L+1,1);
Avz = cell(Number-L+1,1);
A = cell(Number-L+1,1);
A2 = cell(Number-L+1,1);

Xi = cell(Number-L+1,1);
Xi_A2 = cell(Number-L+1,1);
A2u = cell(Number-L+1,1);
for t = 1:Number-L+1
    [O{t},Gamma{t}] = O_Gamma(F,H,nz,L,t); % TODO
    
    if version == 0
        Avz{t} = null(O{t}')';
        if isempty(Avz{t})
            warning('MDM_nullO: Set larger parameter L ')
        end
        
    elseif version == 1
        Avz{t} = null([O{t},Gamma{t} * blkdiag(G{t:t+L-2})]')';
        if isempty(Avz{t})
            warning('MDM_nullO: Set larger parameter L ')
        end
        
    elseif version == 2 % Standard MDM (N=0)
        if rank(O{t})<nx % Observability check
            disp('Set larger parameter L')
            error('MDM_nullO: The state is UNobservable (version=3)')
        end  
        Avz{t} = eye(sum(nz(t:t+L-1)))-O{t}*((O{t}'*O{t})\O{t}'); % Residual maker matrix 
        
    elseif version == 3
        Avz{t} = null(O{t}','r')';
        if isempty(Avz{t})
            warning('MDM_nullO: Set larger parameter L ')
        end           
        
    else 
        error('MDM_nullO: Unspecified version number')
    end
    
    Awu{t} = Avz{t}*Gamma{t};
    
    A = [Awu{t}*blkdiag(E{t:t+L-2}),Avz{t}*blkdiag(D{t:t+L-1})]; % [Aw, Av]
    
    A2{t} = kron(A,A); 
    
    nA2 = size(A2{t},1);
    Xi_part = eye(nA2);
    Xi{t} = Xi_part(1==triu(ones(sqrt(nA2))),:);

    Xi_A2{t} = Xi{t}*A2{t};
    
    A2u{t} = Xi_A2{t}*Psi;
end

covRes = cell(Number-L+1,1); % Covariance of residue vector
for t=1:Number-L+1
    Z = vertcat(z{t:(t+L-1)}); % Augmented measurement vector
    Res = Avz{t} * Z;
    if version ~= 1
        U = vertcat(u{t:(t+L-2)}); % Augmented input vector
        Res = Res - Awu{t} * blkdiag(G{t:t+L-2}) * U;
    end
    covRes{t} = Xi{t}*reshape(Res*Res',size(Res,1)^2,1); 
end

Mat_covRes{1} = Awu;
Mat_covRes{2} = Avz;
Mat_covRes{3} = Xi;
end

