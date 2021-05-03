%%Diving tracks into categories 
%%Matrix A(:,1)=Index of Original Tracks .  A(:,2)=probability to belong to
%%state 1,A(:,3)= to state 2 ....
%%X_Raw

period=1.0;
state=1;
PB=0;%0.0118;%0.0079;    %In Frames 0.0119  %0.045 LG=0.009 in frames
j=state;%4

n=optimalSize;
B=results.state(n).posteriorProb;

clear X_length
j=1;
for i=1:length(X_raw)
    a=length(X_raw{i});
    if length(X_raw{i})>=7
        X_length(j)=a;
        j=j+1;
    end
end


A=[trackInfo.splitIndex,B];
Dat=zeros(max(A(:,1)),n+1);
for i=1:max(A(:,1))
    l=find(A(:,1)==i);
    
    prob=sum(A(l,2:end),1);
    prob=prob./sum(prob);
    prob_st(:,i)=prob;
    Dat(i,:)=[X_length(i),prob];  %Probability and length 
end


    
clear histv CDF y1 tt 

[y1 histv]=histwv(Dat(:,1), Dat(:,j), 7, max(Dat(:,1)), max(Dat(:,1))) ;
tt=7:1:length(y1)+6; %Shortest track
for i=1:length(y1)
    CDF(i)=sum(y1(1:i));
end

[CDF,m1,n1] = unique(CDF,'first');


tt=tt(m1);
 P2=exp(-PB*tt); %H2B Photobleaching parameter


CDF=CDF/CDF(end);

CDF=1-CDF;

CDF=CDF./P2;
CDF=CDF/CDF(1);

tt=tt*period;
plot(tt,CDF,'o')
hold on
final=length(tt)-1;

set(gca,'xscale','log')
set(gca,'yscale','log')

init=1;
%%%%Fit 1 two exponential            
%% Double exponential function of the form (A*(f*(exp(-gamma*x)+(1-f)*(exp(-eta*x)))) with parameters f=p(1) gamma=p(2) eta=p(3) A=p(4)
ft= @(p,xdata)p(4)*(p(1)*(exp(-p(2)*xdata))+(1-p(1))*(exp(-p(3)*xdata)));

%%Fit it to the difference between predictive values and actual values
fitft = @(p)ft(p,tt(init:final))-CDF(init:final);

%%Now the fit for the double exponential
lb=[0,0,0,0];
ub=[1,Inf,Inf,Inf];
p0 = 0.1*ones(1,4); 
problem1 = createOptimProblem('lsqnonlin','x0',p0,'objective',fitft,...
    'lb',lb,'ub',ub);
[Parameter1,errormulti1] = run(MultiStart,problem1,50);

options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
[~,~,Res,~,~,~,J] = lsqnonlin(fitft,p0,lb,ub,options);
CI_Parameter1=nlparci(Parameter1,Res,'jacobian',J);


figure;
scatter(tt,CDF,'.')
hold on
plot(tt,ft(Parameter1,tt),'g')
set(gca,'xscale','log')
set(gca,'yscale','log')
title('Fit 1 (Double Exponential)')
xlabel('Time (Secs)')
ylabel('Survival Distribution')
hold off


%%%%Fit 2 Power-Law
        

%% Power Law Fit of the form A*x^(-b) with parameter A=p(1) b=p(2)
ft2= @(p,xdata)(p(1).*xdata.^(-p(2)));

%%Fit it to the difference between predictive values and actual values
fitft2 = @(p)ft2(p,tt(init:final))-CDF(init:final);

%%Now the fit for the double exponential
lb=[0,0];
ub=[Inf,Inf];
p0 = 0.1*ones(1,2); 
problem2 = createOptimProblem('lsqnonlin','x0',p0,'objective',fitft2,...
    'lb',lb,'ub',ub);
[Parameter2,errormulti2] = run(MultiStart,problem2,50);

[~,~,Res,~,~,~,J] = lsqnonlin(fitft2,Parameter2,lb,ub);
CI_Parameter2=nlparci(Parameter2,Res,'jacobian',J);

figure;
scatter(tt,CDF,'.')
hold on
plot(tt,ft2(Parameter2,tt),'g')
set(gca,'xscale','log')
set(gca,'yscale','log')
title('Fit 2 (Power-Law)')
xlabel('Time (Secs)')
ylabel('CDF')
hold off
% max(CDF(t(init:final)).'-curve2(t(init:final)))



t=tt;
final=length(t)-2;

res1=(((ft(Parameter1,t(init:final))-CDF(init:final)).^2));
res2=(((ft2(Parameter2,t(init:final))-CDF(init:final)).^2));

RSS1=sum(((ft(Parameter1,t(init:final))-CDF(init:final)).^2));   %Double Exponential
RSS2=sum(((ft2(Parameter2,t(init:final))-CDF(init:final)).^2));  %Power Law


BIC11=length(t(init:final))*log(RSS1/length(t(init:final)))+3*log(length(t(init:final)));
BIC12=length(t(init:final))*log(RSS2/length(t(init:final)))+2*log(length(t(init:final)));


DeltaBIC2=-(BIC12-BIC11)

DeltaBIC2/(log(length(t(init:final))))




