function [L,X,U,w,Z] = nestsample(counts,t,K)

% clear all
% close all

C_min = min(counts);
C_max = max(counts);
%% Objects
N = 100;
A = zeros(0,K);
while size(A,1) < N
    a = rand(1,K-1);
    if sum(a) < 1
        a(K) = 1-sum(a);
        A = [A;a];
        clear a
    end
end

F = zeros(0,K);
while size(F,1) < N
    f = rand(1,K)*(C_max-C_min) + C_min;
    if all(diff(f)>0)
        F = [F;f];
    end
end
    


%% Iterate
N = 100; % no of objects
nt = 100; % no of iterations

theta(1:N,:) = [A(1:N,:) F(1:N,:)];
Z0 = 0;X0 = 1;L0 = 0; % Initialize
% log lkhs
LLk = loglike(counts,t,theta);

for i = 1:nt
    [Y,I] = min(LLk);
    L(i) = Y;
    X(i) = exp(-i/N);
    U(i,:) = theta(I,:);
    
    % sample new object from prior constrained by llk > Y
    new_Lk = -inf;
    while new_Lk < Y
        a=[];
        while length(a) < K
            a = rand(1,K-1);
            if sum(a) < 1
                a(K) = 1-sum(a);
            end
        end
        f = [1 0];
        while ~all(diff(f)>0)
            f = rand(1,K)*(C_max-C_min) + C_min;
        end
        new_Lk = loglike(counts,t,[a f]);
    end
    LLk(I) = new_Lk;
    theta(I,:) = [a f];
end

figure(K)
plot(X,L,'*')
C = log(1e-300);
Z = trapz(X,exp(L-C));

for i = 1:length(X)
    if i ~=1
        w(i) = trapz([X(i) X(i-1)],[exp(L(i)-C) exp(L(i-1)-C)]);
    else
        w(i) = trapz([X(i) X0],[exp(L(i)-C) exp(L0-C)]);
    end
end
w = w/Z;
end

%%
function LLk = loglike(counts,t,theta)

K = size(theta,2)/2;
N = size(theta,1);
LLk = [];
for n = 1:N
    Lk=0;
    for i = 1:length(counts)
        temp = 0;
        for j = 1:K
            temp = temp + (theta(n,j)*exp(-1*theta(n,K+j)*t)*(theta(n,K+j)*t)^counts(i))/factorial(counts(i));
        end
        Lk = Lk + log(temp);
    end
    LLk = [LLk Lk];
end
end


