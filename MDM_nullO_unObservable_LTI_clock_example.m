clc
close all
clear all

% MDM, LTI, Unobservable state, Clock-model Example

rng(0) % R2019b

Number = 1e3; % Number of measurements
MC = 1e4; % Number of MC simulations

%%%%%%%%%%%%%%%%%%%% System / Model and Basis matrices %%%%%%%%%%%%%%%%%%%%
if 1
nx = 6;
nw = 6;

nz = 2;
nv = 2;

Ts = 10;
Fs = [1 Ts; 0 1];

F = blkdiag(Fs,Fs,Fs);
u = repmat({0},Number,1);
G = zeros(nx,1);
E = eye(nw); 

H = [1 0 -1 0  0 0;
     1 0  0 0 -1 0];
D = eye(nv);

% nx>rank(obsv(F{1},H{1})) % UNobservable

%%% basis matrices
Qb = cell(0);
Rb = cell(0);

Qb{1} = kron(diag([1 0 0]),[Ts 0; 0 0]);
Qb{2} = kron(diag([1 0 0]),[Ts^3/3 Ts^2/2; Ts^2/2 Ts]);
Qb{3} = kron(diag([0 1 0]),[Ts 0; 0 0]);
Qb{4} = kron(diag([0 1 0]),[Ts^3/3 Ts^2/2; Ts^2/2 Ts]);
Qb{5} = kron(diag([0 0 1]),[Ts 0; 0 0]);
Qb{6} = kron(diag([0 0 1]),[Ts^3/3 Ts^2/2; Ts^2/2 Ts]);

Rb{1} = [1 0;0 0];
Rb{2} = [0 0;0 1];

nQb = size(Qb,2);
nRb = size(Rb,2);

b_true = [6e-19; 5e-21; 2e-18; 3e-20; 7e-19; 4e-21; 8e-18; 1e-17];

Q = zeros(nw);
for idx = 1:nQb
    Q = Q + Qb{idx}*b_true(idx);
end
R = zeros(nv);
for idx = 1:nRb
    R = R + Rb{idx}*b_true(nQb+idx);
end

%%%
end
%%%%%%%%%%%%%%%%% End: System / Model and Basis matrices %%%%%%%%%%%%%%%%%%
 
startOld=now;
for iMC=1:MC
start=now;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data generator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
w = chol(Q)'*randn(nw,Number);
v = chol(R)'*randn(nv,Number); 
x = nan(nx,Number);
z = cell(Number,1);
x(:,1) = ones(nx,1) + randn(nx,1);
for t=1:Number
    if t<Number
        x(:,t+1) = F * x(:,t)  + G * u{t} + E * w(:,t); 
    end
    z{t} = H * x(:,t) + D * v(:,t);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%% End: Data generator %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% MDM matrices and vectors %%%%%%%%%%%%%%%%%%%%%%%%%
if 1
L = 10;
if iMC==1 
    [A2u,covRes,Mat_covRes] = MDM_nullO_LTI(L,F,G,E,nz,H,D,z,u,Qb,Rb);
else
    covRes = MDM_covRes_LTI(L,G,z,u,Mat_covRes);
end
end
%%%%%%%%%%%%%%%%%%%%% End: MDM matrices and vectors %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% MDM - unweighted %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
tic
numMea = Number-L+1;
b_Uw_cov(:,:,iMC) = (numMea*A2u{1}'*A2u{1})\eye(size(A2u{1},2));
b_Uw(:,iMC) = b_Uw_cov(:,:,iMC)*A2u{1}'*sum(horzcat(covRes{:}),2);
time_Uw(iMC) = toc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%% End: MDM - unweighted %%%%%%%%%%%%%%%%%%%%%%%%%%

if mod(iMC,10)==0; disp([' MDM estimate - ', num2str(iMC/(MC)*100),'% ']); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
mean_b_Uw = mean(b_Uw,2);
cov_b_Uw = cov(b_Uw');
est_b_Uw_cov = mean(b_Uw_cov,3);

[b_true,mean_b_Uw]
[diag(cov_b_Uw),sqrt(diag(cov_b_Uw))]

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End: Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
circle = [sin(linspace(0,2*pi,80));cos(linspace(0,2*pi,80))];
figure(1)
for i=1:floor((nQb+nRb)/2)
subplot(1,floor((nQb+nRb)/2),i)
hold all

in = (1:2)+(i-1)*2;

p0 = plot(b_true(in(1)),b_true(in(2)),'mo','MarkerSize',20);

p1 =plot(mean_b_Uw(in(1)),mean_b_Uw(in(2)),'b+','MarkerSize',15);
elips_b_Uw = mean_b_Uw(in) + chol(cov_b_Uw(in,in))'*circle;
p2 = plot(elips_b_Uw(1,:),elips_b_Uw(2,:),'b-');

grid on
%%% Estimates
%plot(b_Uw(in(1),:),b_Uw(in(2),:),'r+','MarkerSize',3);
%%%
end
legend([p0,p1,p2],{'True','MDM - MC sample mean', 'MDM - MC sample covariance'})

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End: Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('Example_UNobservable_clock_null.mat')