clc
close all
clear all

% MDM, LTV, Example TAC and SYSID 

rng(0) % R2019b

Number = 1e3; % Number of measurements
MC = 1e4; % Number of MC simulations

%%%%%%%%%%%%%%%%%%%% System / Model and Basis matrices %%%%%%%%%%%%%%%%%%%%
if 1
nx = 1;
F = cell(Number,1);
G = cell(Number,1);
E = cell(Number,1);
u = cell(Number,1);

H = cell(Number,1);
D = cell(Number,1);
nz = nan(Number,1);

nw = 1;
nv = 1;
Q = 2;
R = 1;
for t=1:Number
    F{t} = 0.8-0.1*sin(7*pi*t/Number);    
    G{t} = 1;   
    u{t} = sin(t/Number);   
    E{t} = 1; 

    nz(t) = 1;
    H{t} = 1+0.99*sin(100*pi*t/Number);
    D{t} = 1;   
end

%%% basis matrices
Qb = cell(0);
Rb = cell(0);

Qb{1}=1;
Rb{1}=1;

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
x(:,1) = 1 + randn(nx,1);
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
if iMC==1 
    [A2u,covRes,Mat_covRes,Xi_A2] = MDM_nullO(L,F,G,E,nz,H,D,z,u,Qb,Rb);
else
    covRes = MDM_covRes(L,G,z,u,Mat_covRes);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%% MDM - weighted(Uw) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
%%%%%% Matrix for weighting %%%%%
if iMC==1
    for timeShift = 0:L-1 
        [EwvLS4_fun{timeShift+1},QRu_sim] = EwvLS4_compute_Fast(L,timeShift,nw,nv);
    end
end
%%% End: Matrix for weighting %%%
tic
%%% 4. moments - Uw estimate is used for weight matrix %%%
nMNumber = size(1:Number-L+1,2);
EwvLS4_all = 0;
for timeShift = 0:L-1
    EwvLS4_Uw = double(subs(EwvLS4_fun{timeShift+1},QRu_sim,b_Uw(:,iMC)));
    EwvLS4_part = kron([zeros(timeShift,nMNumber);eye(nMNumber-timeShift,nMNumber)], EwvLS4_Uw);
    EwvLS4_all = EwvLS4_all + EwvLS4_part;
    if timeShift>0
        EwvLS4_all = EwvLS4_all + EwvLS4_part';
    end 
end
%%% 4. moments - Uw estimate is used for weight matrix %%%

blkdiag_Xi_A2 = blkdiag(Xi_A2{:});
cov_covRes = blkdiag_Xi_A2 * EwvLS4_all * blkdiag_Xi_A2';
inv_cov_covRes = cov_covRes\eye(size(cov_covRes,1));

A2u_We = vertcat(A2u{:});

b_We_cov(:,:,iMC) = (A2u_We'*inv_cov_covRes*A2u_We)\eye(size(A2u_We,2));
b_We(:,iMC) = b_We_cov(:,:,iMC)*A2u_We'*inv_cov_covRes*vertcat(covRes{:});
time_We(iMC) = toc;
end
%%%%%%%%%%%%%%%%%%%%%%%%% End: MDM - weighted(Uw) %%%%%%%%%%%%%%%%%%%%%%%%%

if mod(iMC,10)==0; disp([' MDM estimate - ', num2str(iMC/(MC)*100),'% ']); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
mean_b_Uw = mean(b_Uw,2);
cov_b_Uw = cov(b_Uw');
est_b_Uw_cov = mean(b_Uw_cov,3);

mean_b_We = mean(b_We,2);
cov_b_We = cov(b_We');
est_b_We_cov = mean(b_We_cov,3);

[[Q;R],mean_b_Uw,mean_b_We]

[diag(cov_b_Uw),diag(cov_b_We)]

[diag(est_b_Uw_cov),diag(est_b_We_cov)]

[norm(cov_b_Uw-est_b_Uw_cov)/norm(est_b_Uw_cov),...
 norm(cov_b_We-est_b_We_cov)/norm(est_b_We_cov)]

[mean(time_Uw), mean(time_We);...
 mean(time_Uw)/mean(time_Uw), mean(time_We)/mean(time_Uw)]
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End: Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
figure(1)
hold all
grid on
xlim([1.75 2.25])
ylim([0.85 1.15])
xlabel('Q','Interpreter','tex')
ylabel('R','Interpreter','tex')

circle = [sin(linspace(0,2*pi,80));cos(linspace(0,2*pi,80))];

p0 = plot(Q,R,'ko','MarkerSize',40,'LineWidth',0.5);

p1 = plot(mean_b_Uw(1),mean_b_Uw(2),'b+','MarkerSize',30);
elips_b_Uw = mean_b_Uw + chol(cov_b_Uw)'*circle;
p2 = plot(elips_b_Uw(1,:),elips_b_Uw(2,:),'b-','LineWidth',1.2);
est_elips_b_Uw = mean_b_Uw + chol(est_b_Uw_cov)'*circle;
p3 = plot(est_elips_b_Uw(1,:),est_elips_b_Uw(2,:),'b--','LineWidth',1.2);

p4 = plot(mean_b_We(1),mean_b_We(2),'x','color',[0 .5 .2],'MarkerSize',30);
elips_b_We = mean_b_We + chol(cov_b_We)'*circle;
p5 = plot(elips_b_We(1,:),elips_b_We(2,:),'-','color',[0 .5 .2],'LineWidth',1.2);
est_elips_b_We = mean_b_We + chol(est_b_We_cov)'*circle;
p6 = plot(est_elips_b_We(1,:),est_elips_b_We(2,:),'--','color',[1 0 0],'LineWidth',1.2);

pn = plot(0,nan,'w.');

legend([p0,p1,p4,pn,p2,p5,pn,p3,p6],{'True','Uw - MC sample mean', 'We - MC sample mean',...
                                     '    ','Uw - MC sample std' , 'We - MC sample std',...
                                     '    ','Uw - average est std', 'We - average est std'},'NumColumns',3)
%%% Estimates
%plot(b_Uw(1,:),b_Uw(2,:),'r+','MarkerSize',3);
%plot(b_We(1,:),b_We(2,:),'ro','MarkerSize',3);
%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End: Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('Example_TAC_null.mat')