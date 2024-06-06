function [EwvLS4_fun,QRu,EwvLS4] = EwvLS4_compute(L,timeShift,nw,nv,Qt,Rt)
% Author: Oliver Kost, kost@ntis.zcu.cz
%
% Computation of the structure and value of 4-th moment of the augmented noise vector for different time shift  
% 4-th moment: kron([W(t);V(t)],[W(t);V(t)]) * kron([W(t+timeShift);V(t+timeShift)],[W(t+timeShift);V(t+timeShift)])
%
% timeShift = 0:L-1 % All nonzero cross-covariances

Qu = sym('Q',[nw*(nw+1)/2 1],'real'); 
Ru = sym('R',[nv*(nv+1)/2 1],'real');
QRu = [Qu;Ru];

if L>timeShift % vectors overlap -> fourth moment is non-zero

QRLSu = [kron(ones(L-1+timeShift,1),Qu);kron(ones(L+timeShift,1),Ru)];

Q=sym(zeros(nw)); Q(tril(ones(nw))==1) = Qu; Q = Q+Q'-diag(diag(Q));
R=sym(zeros(nv)); R(tril(ones(nv))==1) = Ru; R = R+R'-diag(diag(R));
QRL = blkdiag(kron(eye(L-1),Q),kron(eye(L),R));
QRLS = blkdiag(kron(eye(L-1+timeShift),Q),kron(eye(L+timeShift),R));

w = sym('w%d_%d',[nw L-1+timeShift],'real');
v = sym('v%d_%d',[nv L+timeShift],'real');
wv = [w(:);v(:)];
wvL = [reshape(w(:,1:L-1),(L-1)*nw,1);reshape(v(:,1:L),L*nv,1)]; 
wvL2 = wvL*wvL';
wvS = [reshape(w(:,(1:L-1)+timeShift),(L-1)*nw,1);reshape(v(:,(1:L)+timeShift),L*nv,1)]; 
wvS2 = wvS*wvS';
wvLS4 = wvS2(:)*wvL2(:)';

%%% 4-th central moment replacement
wvLSu4 = [];
for iLS=1:L+timeShift
    if iLS<L+timeShift
        w2 = w(:,iLS)*w(:,iLS)';
    else
        w2 = [];
    end
    v2 = v(:,iLS)*v(:,iLS)';
    wvLSu4 = [wvLSu4; unique(w2(:)*w2(:)');...
                      unique(v2(:)*v2(:)')];
end
for k=1:size(wvLSu4(:),1)    
    %%% Power of the elements w1_1^2*w2_1^2 -> powEl=[1 1 2 2]
    powEl=[]; 
    for i=1:size(wv,1)
        [~, wvLSu4_part] = coeffs(wvLSu4(k),wv(i)); % Eliminate all elements except the i-th
        
        if wvLSu4_part~=1 % i-th element does not occur
            for l=1:4 % up to 4-th power/moment
                wvLSu4_part = wvLSu4_part/wv(i);
                if wvLSu4_part==1
                    powEl=[powEl,ones(1,l)*i]; % w2_1^2 -> powEl=[powEl, 2 2]
                    break
                end
            end
            if size(powEl,2)==4 
                break
            end
        end
    end
    %%% End: Power of the elements w1_1^2*w2_1^2 -> powEl=[1 1 2 2]

    %%% Isserlis' theorem
    isser = 0; 
    for i=2:4 
        r=2:4; r(i-1)=[];
        isser = isser + QRLS(powEl(1),powEl(i))*QRLS(powEl(r(1)),powEl(r(2)));
    end
    idx_wvLS4 = find(wvLS4==wvLSu4(k)); % indexes of k-th 4-th-moment
    wvLS4(idx_wvLS4) = isser;
    %%% End: Isserlis' theorem

end
%%% End: 4-th central moment replacement

%%% covariance replacement
w2 = [];
v2 = [];
for iLS=1:L+timeShift
    if iLS<L+timeShift
        w2 = [ w2 ; unique(w(:,iLS)*w(:,iLS)','stable')];
    end
    v2 = [ v2 ; unique(v(:,iLS)*v(:,iLS)','stable')];  
end
wvLSu2 = [w2;v2]; 
wvLS4=subs(wvLS4,wvLSu2,QRLSu);
%%% End : covariance replacement

%%% Zeroing of first moments
wvLS4 = subs(wvLS4,wv,zeros((L-1+timeShift)*nw+(L+timeShift)*nv,1));
%%% End: Zeroing of first moments

EwvLS4_fun = wvLS4-QRL(:)*QRL(:)';

if exist('Qt','var') && exist('Rt','var')
    QRtu = [Qt((tril(ones(nw))==1));Rt((tril(ones(nv))==1))];
    EwvLS4 = double(subs(EwvLS4_fun,QRu,QRtu));
else
    EwvLS4 = [];
end

else % 'if L>timeShift' no vector overlap -> zero fourth moment
    EwvLS4 = zeros(((L-1)*nw+L*nv)^2);
    EwvLS4_fun = EwvLS4;
end

end

