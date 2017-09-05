% clear all
close all

%% test Data
mu1 = 50;
mu2 = 100;
t=1;
counts = [poissrnd(mu1,50,1); poissrnd(mu2,100,1)];

%%
nK = 3;
for k = 1:nK
    [L(k,:),X(k,:),U,w,Z(k)] = nestsample(counts,t,k);
%     display(['Mean values are' num2str(sum(w.U))])
end

mean(U)

figure(nK+1)
s=subplot(1,1,1);
plot(1:nK,Z/Z(1),'*');
set(s,'XTick',[1 2 3])
xlabel('number of states')
ylabel('Bayes Factor')
set(get(s,'Title'),'string',['Model: k = ' num2str(nK)])