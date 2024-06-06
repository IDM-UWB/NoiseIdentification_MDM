function [EwvLS4_fun,QRu_sim,EwvLS4] = EwvLS4_compute_Fast(L,timeShift,nw,nv,Q,R)
% Author: Oliver Kost, kost@ntis.zcu.cz
%
% Requires files: EwvLS4_compute.m
%
% Computation of the structure and value of 4-th moment of the augmented noise vector for different time shift  
% 4-th moment: kron([W(t);V(t)],[W(t);V(t)]) * kron([W(t+timeShift);V(t+timeShift)],[W(t+timeShift);V(t+timeShift)])

    if exist('EwvLS4_compute_Data.mat')==2
        y=load('EwvLS4_compute_Data.mat'); 
        if ~isempty(find(vecnorm([L,timeShift,nw,nv]-y.Pall,1,2)==0))
            EwvLS4_fun = y.D(L,timeShift+1,nw+1,nv+1).EwvLS4_fun;
            QRu_sim = y.D(L,timeShift+1,nw+1,nv+1).QRu_sim;
        else
            [EwvLS4_fun,QRu_sim] = EwvLS4_compute(L,timeShift,nw,nv);

            Pall=[y.Pall;[L,timeShift,nw,nv]];
            D=y.D;
            D(L,timeShift+1,nw+1,nv+1).EwvLS4_fun = EwvLS4_fun;
            D(L,timeShift+1,nw+1,nv+1).QRu_sim = QRu_sim;
            save('EwvLS4_compute_Data.mat','D','Pall')   
            disp(['Saved - ',num2str([L,timeShift,nw,nv])])
        end
    else
        [EwvLS4_fun,QRu_sim] = EwvLS4_compute(L,timeShift,nw,nv);

        Pall=[L,timeShift,nw,nv];
        D(L,timeShift+1,nw+1,nv+1).EwvLS4_fun=EwvLS4_fun;
        D(L,timeShift+1,nw+1,nv+1).QRu_sim=QRu_sim;
        save('EwvLS4_compute_Data.mat','D','Pall')
        disp(['Saved - ',num2str([L,timeShift,nw,nv])])    
    end

    if exist('Q','var') && exist('R','var')
        QRu = [Q(tril(ones(nw))==1);R(tril(ones(nv))==1)];
        EwvLS4 = double(subs(EwvLS4_fun,QRu_sim,QRu));
    else
        EwvLS4 = [];
    end
end




    
