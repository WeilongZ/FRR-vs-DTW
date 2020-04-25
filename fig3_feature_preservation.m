%% Show that DTW will change both curves
clear;clc;close all;

% generate two functions
t=linspace(0,1,501);
dt=mean(diff(t));
f1=5*sin(t.^3*pi*4);
f2=sin(t*pi*4);
figure()
set(gcf, 'position', [100 100 800 600]);
plot(t,f1,'b-',t,f2,'r-','LineWidth',2);
% title('Original f_1 and f_2');
legend('f_1','f_2','location','best');
% xlabel('t');
% ylabel('f');
set(gca,'FontSize',18);
set(gca,'Ytick',[-6 -3 0 3 6]);ylim([-6 6]);
% ----------------- time warpping: match q2 to q1 -----------------

% FRR method: --------------------
% obtain SRVF
df1=gradient(f1,dt);df2=gradient(f2,dt);
q1=df1./sqrt(abs(df1)+eps);q2=df2./sqrt(abs(df2)+eps);

% gam12=DPTW(q2,q1,t,3);
gam12=DynamicProgrammingQ(q1,q2,0,0);
gam12=(gam12-gam12(1))/(gam12(end)-gam12(1));
q1tw=interp1(t,q1,gam12).*sqrt(gradient(gam12,dt));
f1tw=interp1(t,f1,gam12);
figure()
set(gcf, 'position', [100 100 800 600]);
plot(t,f1tw,'b-',t,f2,'r-','LineWidth',2);
% title('warpped f_1 and f_2: FRR');
legend('Aligned f_1','f_2','location','best');
% xlabel('t');
% ylabel('f');
set(gca,'FontSize',18);
set(gca,'Ytick',[-6 -3 0 3 6]);ylim([-6 6]);

% DTW method: -------------------
[dist_dtw,idx1,idx2]=DTW_2(f1,f2,1000,1);

figure()
set(gcf, 'position', [100 100 800 600]);
plot(f1(idx1),'LineWidth',2);
hold on
plot(f2(idx2),'LineWidth',2);
% title({'Warpped f_1 and Warpped f_2: DTW',['DTW distance = ',num2str(dist_dtw)]});
legend('Aligned f_1','Aligned f_2','location','best');
% ylabel('f');
set(gca,'FontSize',18);
set(gca,'Xtick',[0 300 600 900]);xlim([0 900]);
set(gca,'Ytick',[-6 -3 0 3 6]);ylim([-6 6]);


% ----------------- time warpping: match q1 to q2 -----------------

% FRR method: --------------------

% gam12=DPTW(q2,q1,t,3);
gam21=DynamicProgrammingQ(q2,q1,0,0);
gam21=(gam21-gam21(1))/(gam21(end)-gam21(1));
q2tw=interp1(t,q2,gam21).*sqrt(gradient(gam21,dt));
f2tw=interp1(t,f2,gam21);
figure()
set(gcf, 'position', [100 100 800 600]);
plot(t,f1,'b-',t,f2tw,'r-','LineWidth',2);
% title('warpped f_1 and f_2: FRR');
legend('f_1','Aligned f_2','location','best');
% xlabel('t');
% ylabel('f');
set(gca,'FontSize',18);
set(gca,'Ytick',[-6 -3 0 3 6]);ylim([-6 6]);

% % DTW method: -------------------
% [dist_dtw,idx2,idx1]=DTW_2(f2,f1,1000,1);
% 
% figure()
% plot(f1(idx1));
% hold on
% plot(f2(idx2));
% title({'warpped f_1 and f_2: DTW',['DTW distance = ',num2str(dist_dtw)]});
% legend('f_1','f_2');
% ylabel('f');
% set(gca,'FontSize',14);




