%% the script to show the algorithm of DTW
clear;clc;close all;

% two series of points as the two curves
c1=[0,0,1,1,2,4,3,1,0];
c2=[0,1,1,1,4,4,2,3,1,1,0];

% build the matrix C_ij
l1=length(c1);l2=length(c2);
Cmat=inf*ones(l1,l2);
Cmat(1,1)=(c1(1)-c2(1))^2;

% squared euclidean to be the metric
dp=@(x,y)((x-y).^2);

% do not apply sliding window in the illustration
w=max(l1,l2);

for i=1:l1
    for j=max(i-w,1):min(l2,i+w)
        Cmat(i,j)=dp(c1(i),c2(j));
        if i>1 && j>1
            Cmat(i,j)=Cmat(i,j)+min([Cmat(i-1,j),Cmat(i,j-1),Cmat(i-1,j-1)]);
        elseif i>1 && j==1
            Cmat(i,j)=Cmat(i,j)+Cmat(i-1,j);
        elseif i==1 && j>1
            Cmat(i,j)=Cmat(i,j)+Cmat(i,j-1);
        end
    end
end
% in the end, we should have the last points connected
dist=Cmat(l1,l2);

% save as a csv file
csvwrite('Cmat.csv',[[0,c2];[c1',Cmat]])

idxmat=[l1,l2];
not11=1;
while not11
    i=idxmat(1,1);j=idxmat(1,2);
    if i>1 && j>1
        [~,rt]=min([Cmat(i-1,j-1),Cmat(i-1,j),Cmat(i,j-1)]);
        temp=[i-1,j-1;i-1,j;i,j-1];
        idxmat=[temp(rt,:);idxmat];
    elseif i==1 && j>1
        idxmat=[i,j-1;idxmat];
    elseif i>1 && j==1
        idxmat=[i-1,j;idxmat];
    else
        not11=0;
    end
end
idx1=idxmat(:,1)';
idx2=idxmat(:,2)';

% show the original and aligned curves
figure(1)
% original
% subplot(2,1,1)
set(gcf, 'position', [100 100 1250 500]);
plot(1:l1,c1,'b-');
hold on
plot(1:l2,c2,'r-');
hold off
legend('Curve 1','Curve 2','location','best');
set(gca,'fontsize',20)

% aligned
% subplot(2,1,2)
figure(2)
set(gcf, 'position', [100 100 1250 500]);
plot(1:length(idx1),c1(idx1),'b-');
hold on
plot(1:length(idx2),c2(idx2),'r-');
hold off
legend('Curve 1','Curve 2','location','best');
set(gca,'fontsize',20)

% % as a comparison, the dtw function in matlab
% figure(3)
% dtw(c1,c2);




