clear all;close all;
n = 10; fig = 1;lambda = 0;
t = 0:0.01:1;
A = (rand(1,n)*4)+1;
alpha = rand(1,n)*2+1;

F = zeros(n,size(t,2));
for i = 1:n
    F(i,:) = A(i)*sin(2*pi*t.^alpha(i));
end

figure(fig)
fig = fig + 1;
plot(t,F,'b-','LineWidth',2);
set(gca,'fontsize',18);
set(gca,'Ytick',[-6 -3 0 3 6]);ylim([-6 6]);
set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1]);xlim([0 1]);

f = F';
binsize = mean(diff(t));
[M, N] = size(f);
f0 = f;

for i = 1:N
    q(:,i) = gradient(f(:,i), binsize)./sqrt(abs(gradient(f(:,i), binsize))+eps);
end

%%% set initial using the original f space
disp(sprintf('\n Initializing...\n'));
mnq = mean(q,2);
dqq = sqrt(sum((q - mnq*ones(1,N)).^2,1));
[ignore, min_ind] = min(dqq);
mq = q(:,min_ind); 
mf = f(:,min_ind);

for k = 1:N
    q_c = q(:,k,1)'; mq_c = mq';
    gam0 = DynamicProgrammingQ(q_c/norm(q_c),mq_c/norm(mq_c),lambda,0);
    gam(k,:) = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale
end
gamI = SqrtMeanInverse(gam);
gamI_dev = gradient(gamI, 1/(M-1));
%mq = interp1(t, mq, (t(end)-t(1)).*gamI + t(1))'.*sqrt(gamI_dev');
mf = interp1(t, mf, (t(end)-t(1)).*gamI + t(1))';
mq = gradient(mf, binsize)./sqrt(abs(gradient(mf, binsize))+eps);


%%% compute mean
disp(sprintf(' Computing Karcher mean of %d functions in SRVF space...\n',N));
ds = inf; 
MaxItr = 15;
for r = 1:MaxItr
    disp(sprintf('updating step: r=%d', r)); 
    if r == MaxItr
        disp(sprintf('maximal number of iterations is reached. \n'));
    end   
    
    %%%% Matching Step %%%%
    clear gam gam_dev;
    % use DP to find the optimal warping for each function w.r.t. the mean
    for k = 1:N
        q_c = q(:,k,1)'; mq_c = mq(:,r)';
        gam0 = DynamicProgrammingQ(q_c/norm(q_c),mq_c/norm(mq_c),lambda,0);
        gam(k,:) = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale
        gam_dev(k,:) = gradient(gam(k,:), 1/(M-1));
%        q(:,k,r+1) = interp1(t, q(:,k,1), (t(end)-t(1)).*gam(k,:) + t(1))'.*sqrt(gam_dev(k,:)');
        f(:,k,r+1) = interp1(t, f(:,k,1), (t(end)-t(1)).*gam(k,:) + t(1))';
        q(:,k,r+1) = gradient(f(:,k,r+1), binsize)./sqrt(abs(gradient(f(:,k,r+1), binsize))+eps);
    end
    
    ds(r+1) = sum(trapz(t, (mq(:,r)*ones(1,N)-q(:,:,r+1)).^2)) + ...
                lambda*sum(trapz(t, (1-sqrt(gam_dev')).^2));
    
    %%%% Minimization Step %%%%
    % compute the mean of the matched function
    mq(:,r+1) = mean(q(:,:,r+1),2);
    
    qun(r) = norm(mq(:,r+1)-mq(:,r))/norm(mq(:,r));
    if qun(r) < 1e-2 && r >=10
        break;
    end
end

if lambda == 0
    disp(sprintf('additional run for adjustment')); 
    r = r+1;
    for k = 1:N
        q_c = q(:,k,1)'; mq_c = mq(:,r)';
        gam0 = DynamicProgrammingQ(q_c/norm(q_c),mq_c/norm(mq_c),lambda,0);
        gam(k,:) = (gam0-gam0(1))/(gam0(end)-gam0(1));  % slight change on scale
        gam_dev(k,:) = gradient(gam(k,:), 1/(M-1));      
    end   
    gamI = SqrtMeanInverse(gam);
    gamI_dev = gradient(gamI, 1/(M-1));
    mq(:,r+1) = interp1(t, mq(:,r), (t(end)-t(1)).*gamI + t(1))'.*sqrt(gamI_dev');
    for k = 1:N
        q(:,k,r+1) = interp1(t, q(:,k,r), (t(end)-t(1)).*gamI + t(1))'.*sqrt(gamI_dev');
        f(:,k,r+1) = interp1(t, f(:,k,r), (t(end)-t(1)).*gamI + t(1))';  
    %    q(:,k,r+1) = gradient(f(:,k,r+1), binsize)./sqrt(abs(gradient(f(:,k,r+1), binsize))+eps);
        gam(k,:) = interp1(t, gam(k,:), (t(end)-t(1)).*gamI + t(1));
    end     
end

% aligned data
fn = f(:,:,r+1);

figure(fig)
fig = fig + 1;
plot((0:M-1)/(M-1), gam,'k-', 'linewidth', 2);
axis square;
set(gca,'fontsize',18);
set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1]);xlim([0 1]);
set(gca,'Ytick',[0 0.2 0.4 0.6 0.8 1]);ylim([0 1]);
% title('Warping functions');

figure(fig)
fig = fig + 1;
plot(t, fn,'r-', 'LineWidth',2);
set(gca,'fontsize',18);
% title('Warped data');
set(gca,'Ytick',[-6 -3 0 3 6]);ylim([-6 6]);
set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1]);xlim([0 1]);

mean_f0 = mean(f0, 2); 
std_f0 = std(f0, 0, 2);
figure(fig)
fig = fig + 1;clf;
plot(t, mean_f0, 'b-', 'linewidth', 2); hold on;
plot(t, mean_f0+std_f0, 'r-', 'linewidth', 1);
plot(t, mean_f0-std_f0, 'g-', 'linewidth', 1);
set(gca,'fontsize',18);
title('Original data: Mean \pm STD');
set(gca,'Ytick',[-6 -3 0 3 6]);ylim([-6 6]);
set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1]);xlim([0 1]);


mean_fn = mean(fn, 2); 
std_fn = std(fn, 0, 2);
figure(fig)
fig = fig + 1;clf;
plot(t, mean_fn, 'b-', 'linewidth', 2); hold on;
plot(t, mean_fn+std_fn, 'r-', 'linewidth', 1);
plot(t, mean_fn-std_fn, 'g-', 'linewidth', 1);
set(gca,'fontsize',18);
title('Warped data: Mean \pm STD');
set(gca,'Ytick',[-6 -3 0 3 6]);ylim([-6 6]);
set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1]);xlim([0 1]);

COV_F = cov(F);
[V_F D_F] = svd(COV_F);
COV_fn = cov(fn');
[V_fn D_fn] = svd(COV_fn);

% Karcher mean of warping
mu_gam = SqrtMeanInverse(gam);
mu_gam = invertGamma(mu_gam);
% SRVF of warping function
for i = 1:n
    psi(i,:) = sqrt(diff(gam(i,:))/(1/(M-1)));
end
psi_mu = sqrt(diff(mu_gam)/(1/(M-1)));
% compute tangent vector
for i = 1:n
    theta = acos(sum(psi_mu.*psi(i,:))*(1/(M-1)));
    v(i,:) = (theta/sin(theta))*(psi(i,:) - cos(theta)*psi_mu);
end

COV_gam = cov(v);
[V_gam D_gam] = svd(COV_gam);
FP_gam_V = V_gam(:,1);
tmp = cos(1)*psi_mu + (sin(1))*FP_gam_V';
tmp_gam = [0 cumsum(tmp.*tmp)]/M;
FP_gam = (tmp_gam-min(tmp_gam))/(max(tmp_gam)-min(tmp_gam));
% FP_gam = invertGamma(tmp_gam);

figure(fig)
fig = fig + 1;
plot(t,V_F(:,1),'linewidth', 2);
set(gca,'fontsize',18);
title(['First CP for original: ' num2str(D_F(1,1)/sum(sum(D_F)))]);

figure(fig)
fig = fig + 1;
plot(t,V_fn(:,1),'linewidth', 2);
set(gca,'fontsize',18);
title(['First CP for aligned: ' num2str(D_fn(1,1)/sum(sum(D_fn)))]);

figure(fig)
fig = fig + 1;
plot((0:M-1)/(M-1),FP_gam,'linewidth', 2);
set(gca,'fontsize',18);
title(['First CP for warping functions: ' num2str(D_gam(1,1)/sum(sum(D_gam)))]);

mean_f0 = mean(f0, 2); 
fp_f0 = V_F(:,1);
figure(fig)
fig = fig + 1;clf;
p1 = plot(t, mean_f0, 'b-', 'linewidth', 2); hold on;
p2 = plot(t, mean_f0+sqrt(D_F(1,1))*fp_f0, 'r-', 'linewidth', 2);
p3 = plot(t, mean_f0-sqrt(D_F(1,1))*fp_f0, 'r-', 'linewidth', 2);
p4 = plot(t, mean_f0+2*sqrt(D_F(1,1))*fp_f0, 'g-', 'linewidth', 2);
p5 = plot(t, mean_f0-2*sqrt(D_F(1,1))*fp_f0, 'g-', 'linewidth', 2);
set(gca,'fontsize',18);
% title('Original data: Mean \pm FP');
set(gca,'Ytick',[-6 -3 0 3 6]);ylim([-7 7]);
set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1]);xlim([0 1]);
l1=legend([p1 p3 p5],'Mean','Mean$\pm\sqrt{\lambda_1}f_1$','Mean$\pm2\sqrt{\lambda_1}f_1$');
set(l1,'Interpreter','latex');
mean_fn = mean(fn, 2); 
fp_fn = V_fn(:,1);
figure(fig)
fig = fig + 1;clf;
p1 = plot(t, mean_fn, 'b-', 'linewidth', 2); hold on;
p2 = plot(t, mean_fn+sqrt(D_fn(1,1))*fp_fn, 'r-', 'linewidth', 2);
p3 = plot(t, mean_fn-sqrt(D_fn(1,1))*fp_fn, 'r-', 'linewidth', 2);
p4 = plot(t, mean_fn+2*sqrt(D_fn(1,1))*fp_fn, 'g-', 'linewidth', 2);
p5 = plot(t, mean_fn-2*sqrt(D_fn(1,1))*fp_fn, 'g-', 'linewidth', 2);
set(gca,'fontsize',18);
% title('Warped data: Mean \pm FP');
set(gca,'Ytick',[-6 -3 0 3 6]);ylim([-7 7]);
set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1]);xlim([0 1]);
% legend('Mean','Mean + FP','Mean - FP','Mean + 2*FP','Mean - 2*FP');
l1=legend([p1 p3 p5],'Mean','Mean$\pm\sqrt{\lambda_1^{(a)}}f_1^{(a)}$','Mean$\pm2\sqrt{\lambda_1^{(a)}}f_1^{(a)}$');
set(l1,'Interpreter','latex');

FP_gam_V = V_gam(:,1);
tmp_fp11 = psi_mu + sqrt(D_gam(1,1))*FP_gam_V';
tmp_fp10 = psi_mu - sqrt(D_gam(1,1))*FP_gam_V';
tmp_fp21 = psi_mu + 2*sqrt(D_gam(1,1))*FP_gam_V';
tmp_fp20 = psi_mu - 2*sqrt(D_gam(1,1))*FP_gam_V';

tmp_11 = cos(1)*psi_mu + (sin(1))*tmp_fp11;
tmp_10 = cos(1)*psi_mu + (sin(1))*tmp_fp10;
tmp_21 = cos(1)*psi_mu + (sin(1))*tmp_fp21;
tmp_20 = cos(1)*psi_mu + (sin(1))*tmp_fp20;

tmp_gam_11 = [0 cumsum(tmp_11.*tmp_11)]/M;
tmp_gam_10 = [0 cumsum(tmp_10.*tmp_10)]/M;
tmp_gam_21 = [0 cumsum(tmp_21.*tmp_21)]/M;
tmp_gam_20 = [0 cumsum(tmp_20.*tmp_20)]/M;

FP_gam_11 = (tmp_gam_11-min(tmp_gam_11))/(max(tmp_gam_11)-min(tmp_gam_11));
FP_gam_10 = (tmp_gam_10-min(tmp_gam_10))/(max(tmp_gam_10)-min(tmp_gam_10));
FP_gam_21 = (tmp_gam_21-min(tmp_gam_21))/(max(tmp_gam_21)-min(tmp_gam_21));
FP_gam_20 = (tmp_gam_20-min(tmp_gam_20))/(max(tmp_gam_20)-min(tmp_gam_20));

figure(fig)
fig = fig + 1;clf;
p1 = plot(t, mu_gam, 'b-', 'linewidth', 2); hold on;
p2 = plot(t, FP_gam_11, 'r-', 'linewidth', 2);
p3 = plot(t, FP_gam_10, 'r-', 'linewidth', 2);
p4 = plot(t, FP_gam_21, 'g-', 'linewidth', 2);
p5 = plot(t, FP_gam_20, 'g-', 'linewidth', 2);
axis square;
set(gca,'fontsize',18);
% title('Warped data: Mean \pm FP');
set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1]);xlim([0 1]);
l1=legend([p1 p3 p5],'Mean','Mean$\pm\sqrt{\lambda_1^{(p)}}f_1^{(p)}$','Mean$\pm2\sqrt{\lambda_1^{(p)}}f_1^{(p)}$');
set(l1,'Interpreter','latex');

