clc
close all
clear all

% MDM, LTV, Unobservable state and Unkown control Example

rng(0) % R2019b

Number = 1e3; % Number of measurements
MC = 1e4; % Number of MC simulations

%%%%%%%%%%%%%%%%%%%% System / Model and Basis matrices %%%%%%%%%%%%%%%%%%%%
if 1
nx = 3;
F = cell(Number,1);
G = cell(Number,1);
E = cell(Number,1);
u = cell(Number,1);

H = cell(Number,1);
D = cell(Number,1);
nz = nan(Number,1);

nw = 3;
nv = 3;
b_true = [1;1;-1;2;2;1];
Q = b_true(1)*[1 0 0;0 1 0;0 0 1] + b_true(2)*[0 0 0;0 1 0;0 0 1] + b_true(3)*[0 -1 0;-1 0 -1;0 -1 0];
R = b_true(4)*[1 0 0; 0 0 0;0 0 1] + b_true(5)*[0 0 0; 0 2 0;0 0 0] + b_true(6)*[0 0 1; 0 0 1;1 1 0];

for t=1:Number
    F{t} = [1  2    1;
            0 -1.01 2;
            0  0    1];    
    u{t} = sin(t/Number);
    G{t} = [0; cos(10*t/Number); 1];   
    E{t} = [-3 2 0;2 2 2;5 0 1];

    nz(t) = 3;
    H{t} = [0 1 0;...
            0 0 2;
            0 1 1];
    D{t} = [1 1 0;0 2 1;1 0 -1];
end

% nx>rank(obsv(F{1},H{1})) % UNobservable

%%% basis matrices
Qb = cell(0);
Rb = cell(0);

Qb{1}=[1 0 0;0 1 0;0 0 1];
Qb{2}=[0 0 0;0 1 0;0 0 1];
Qb{3}=[0 -1 0;-1 0 -1;0 -1 0];

Rb{1}=[1 0 0;0 0 0;0 0 1];
Rb{2}=[0 0 0;0 2 0;0 0 0];
Rb{3}=[0 0 1;0 0 1;1 1 0];

nQb = size(Qb,2);
nRb = size(Rb,2);
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
        x(:,t+1) = F{t} * x(:,t)  + G{t} * u{t} + E{t} * w(:,t); 
    end
    z{t} = H{t} * x(:,t) + D{t} * v(:,t);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%% End: Data generator %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% MDM matrices and vectors %%%%%%%%%%%%%%%%%%%%%%%%%
if 1
L = 2;
version = 1; % Unknown control
if iMC==1 
    [A2u,covRes,Mat_covRes] = MDM_nullO(L,F,G,E,nz,H,D,z,[],Qb,Rb,version);
else
    covRes = MDM_covRes(L,[],z,[],Mat_covRes,version);
end
end
%%%%%%%%%%%%%%%%%%%%% End: MDM matrices and vectors %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% MDM - unweighted %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
tic
A2u_Uw = vertcat(A2u{:});
b_Uw_cov(:,:,iMC) = (A2u_Uw'*A2u_Uw)\eye(size(A2u_Uw,2));
b_Uw(:,iMC) = b_Uw_cov(:,:,iMC)*A2u_Uw'*vertcat(covRes{:});
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
est_elips_b_Uw = mean_b_Uw(in) + chol(est_b_Uw_cov(in,in))'*circle;
p3 = plot(est_elips_b_Uw(1,:),est_elips_b_Uw(2,:),'-.','color',[.1 .1 .5]);

grid on
%%% Estimates
%plot(b_Uw(in(1),:),b_Uw(in(2),:),'r+','MarkerSize',3);
%%%
end
legend([p0,p1,p2,p3],{'True','MDM - MC sample mean', 'MDM - MC sample covariance', 'MDM - Est covariance'})

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End: Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('Example_UNobservable_UnknownU_null.mat')