clear all
clc

%files = dir('*.mat');
files = uigetfile('*.mat', 'MultiSelect', 'on');
k=1;
for i=1:length(files)
%     eval(['load ' files(i).name]);

    eval(['load ' files{i}]);
    
    A=Results.Tracking.Tracks;
   
    for j=1:A(end,4)
        L=find(A(:,4)==j);
        x=A(L,1);%*0.117;
        y=A(L,2);%*0.117;
        X{k}=[x y];
        k=k+1;
    end
end
save('TracksPIV.mat','X');


%Pixel Size 0.117