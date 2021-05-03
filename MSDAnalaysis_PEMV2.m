n=optimalSize; %number of states

prompt = {'Enter state number:','Enter probability threshold:'};
dlgtitle = 'Parameter Input';
dims = [1 35];
definput = {'1','0.6'};
answer = inputdlg(prompt,dlgtitle,dims,definput)

A=results.state(n).posteriorProb;
th=str2num(answer{2});
i=str2num(answer{1});

% L1=find(A(:,1)>=th);
% L2=find(A(:,2)>=th);
% L3=find(A(:,3)>=th);
% L4=find(A(:,4)>=th);
% L5=find(A(:,5)>=th);
% L6=find(A(:,6)>=th);
% L7=find(A(:,7)>=th);
% L8=find(A(:,8)>=th);
L=find(A(:,i)>=th);

tracks=cell(length(L),1);
zz=[];
for i=1:(length(L))%length(L)
%     B=results.X{L(i)};
    B=X{L(i)};
    x=B(:,1)*0.117;%0.117 for 140X obj and 0.164 for 100x silicon
    y=B(:,2)*0.117;%0.117
    T=1:length(x);
    T=T*0.012;
    tracks{i}=[T.' x y];
    
    %Added to calculate distribution
    %D1=sqrt(x.^2+y.^2);
    %zz=[zz;abs(diff(D1))];
end

ma = msdanalyzer(2, 'µm', 's');
ma = ma.addAll(tracks);

% figure;
% [hps, ha] = ma.plotTracks;
% ma.labelPlotTracks(ha);

ma = ma.computeMSD;
%figure
%hmsd = ma.plotMeanMSD(gca, true);


%%%%
%%Error Bars
figure;
mmsd = ma.getMeanMSD;
t = mmsd(:,1);
x = mmsd(:,2);
dx = mmsd(:,3) ./ sqrt(mmsd(:,4));
    errorbar(t, x, dx, 'k')


    
    