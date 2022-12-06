% first assign your TR to dt, for example
dt = 0.1; % in seconds -- this should result in functions with the same temporal
% resolution as the main dataset
% then several excerts from spm_get_bf.m
% canonical hemodynamic response function
%---------------------------------------------------------------
p = [6 16 1 1 6 0 32];
bf = spm_hrf(dt,p);

[bf6 p] = spm_hrf(dt, [6 16 1 1 6 0 32]);
[bf5 p] = spm_hrf(dt, [5 16 1 1 6 0 32]);
[bf4 p] = spm_hrf(dt, [4 16 1 1 6 0 32]);
[bf3 p] = spm_hrf(dt, [3 16 1 1 6 0 32]);
[bf2 p] = spm_hrf(dt, [2 16 1 1 6 0 32]);
[bf1 p] = spm_hrf(dt, [1 16 1 1 6 0 32]);

min=-0.1;
max=1;

BF1 = rescale(bf1,min,max);
BF2 = rescale(bf2,min,max);
BF3 = rescale(bf3,min,max);
BF4 = rescale(bf4,min,max);
BF5 = rescale(bf5,min,max);
BF6 = rescale(bf6,min,max);


tiledlayout(1,6) % Requires R2019b or later
nexttile
plot(BF6, 'k', 'LineWidth',2.0);
nexttile
plot(BF6, 'k'); hold on 
plot(BF5, 'b');
nexttile
plot(BF6, 'k'); hold on 
plot(BF4, 'b');
nexttile
plot(BF6, 'k'); hold on 
plot(BF3, 'b');
nexttile
plot(BF6, 'k'); hold on 
plot(BF2, 'b');
nexttile
plot(BF6, 'k'); hold on 
plot(BF1, 'b');

plot(BF1, 'LineWidth',2.0); hold on
plot(BF2, 'LineWidth',2.0); hold on
plot(BF3, 'LineWidth',2.0); hold on 
plot(BF4, 'LineWidth',2.0); hold on 
plot(BF5, 'LineWidth',2.0); hold on
plot(BF6, 'LineWidth',2.0); hold on
xticks([0 50 100 150 200 250 300 350]); hold on
xticklabels({'0','5','10','15','20','25','30', '35'});
xlabel('time (s)'); hold on
ylabel('arbitrary units'); hold on
exportgraphics(gcf,'mani_HRFs.eps','ContentType','vector');
% add time derivative
%---------------------------------------------------------------
dp = 1;
p(6) = p(6) + dp;
D = (bf(:,1) - spm_hrf(dt,p))/dp;
bf = [bf D(:)];
p(6) = p(6) - dp;
% add dispersion derivative
%--------------------------------------------------------
dp = 0.01;
p(3) = p(3) + dp;
D = (bf(:,1) - spm_hrf(dt,p))/dp;
bf = [bf D(:)];

plot(mean(PARA(1,:))*bf(:,1))

LC_hrf = [mean(PARA(1,:))*bf(:,1)+mean(PARA(2,:))*bf(:,2)+mean(PARA(3,:))*bf(:,3)]; 

LC_hrf=rescale(LC_hrf)
spm_hrf = rescale(bf(:,1))

plot(spm_hrf); hold on
plot(LC_hrf)
