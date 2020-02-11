%========================================================================%
%                        Title: Code for Econ 771 PS3                    %
%                        Author: Alex Marsh                              %
%                        Date: 03/04/2020                                %
%========================================================================%

%%
% Set Parameters and Objects
rng(123); %set rng seed

M     = 5000;         %set number of simulations
n     = 30;           %set sample size
beta  = [1 2 3]';     %true beta values
alpha = 0.05;         %set test size
K     = size(beta,1); %store size of beta vector

%%
% Exercise 1: Part 1a)

beta_hats1 = zeros(3,M); %preallocate beta_hat matrix
se_hats1   = zeros(3,M); %preallocate standard error matrix

%--------------------------------------------------------------------%
%                     Monte Carlo Simulation Loop                    %
%--------------------------------------------------------------------%

for i = 1:M
    X               = ones(n,3);              %intialize X
    X(:,2:3)        = randn(n,2);             %gen X2 and X3
    Y               = X*beta+randn(n,1);      %gen Y
    beta_hats1(:,i) = (X'*X)\(X'*Y);          %estimate beta_hat
    es              = Y-X*beta_hats1(:,i);    %form residuals
    sigma2          = ((es'*es)/(n-K));       %estimate var of errors
    vcov_mat        = sigma2*((X'*X)\eye(K)); %calculate vcov matrix
    se_hats1(:,i)   = sqrt(diag(vcov_mat));   %store standard errors
end

%%
% Exercise 1: Part 1b)

t_stats1 = (beta_hats1-beta)./se_hats1; %calculate T-stats

%generate figure for beta_1
figure();
histogram(t_stats1(1,:),'Normalization','pdf');
title('Histogram of Standardized $\hat{\beta}_1$','interpreter','latex')
hold on

x = min(t_stats1(1,:)):0.1:max(t_stats1(1,:));
y = normpdf(x);
plot(x,y);

hold off

%generate figure for beta_2
figure();
histogram(t_stats1(2,:),'Normalization','pdf');
title('Histogram of Standardized $\hat{\beta}_2$','interpreter','latex')
hold on

x = min(t_stats1(2,:)):0.1:max(t_stats1(2,:));
y = normpdf(x);
plot(x,y);

hold off

%generate figure for beta_3
figure();
histogram(t_stats1(3,:),'Normalization','pdf');
title('Histogram of Standardized $\hat{\beta}_3$','interpreter','latex')
hold on
x = min(t_stats1(3,:)):0.1:max(t_stats1(3,:));
y = normpdf(x);
plot(x,y);
hold off

%%
% Exercise 1: Part 1c)
p_vals_c1    = 2*(1-tcdf(abs(t_stats1),n-K)); %calculate pvalues
decisions_c1 = p_vals_c1<=alpha;              %acccept or reject
sizes1       = mean(decisions_c1,2)           %calculate size

%%
% Exercise 1: Part 1d)
t_stats1d    = beta_hats1./se_hats1;          %calculate new tstats
p_vals_d1    = 2*(1-tcdf(abs(t_stats1d),n-K));%calculate pvalues
decisions_d1 = p_vals_d1<=alpha;              %accept or reject
power1       = mean(decisions_d1,2)           %calculate power


%%
% Exercise 1: Part 2a)
rng(123); %reset rng seed to get the same data

beta_hats2 = zeros(3,M); %preallocate beta_hat matrix 
se_hats2   = zeros(3,M); %preallocate standard error matrix

%--------------------------------------------------------------------%
%                     Monte Carlo Simulation Loop                    %
%--------------------------------------------------------------------%

for i = 1:M
    X               = ones(n,3);                  %initialize X
    X(:,2:3)        = randn(n,2);                 %genereate X2 & X2
    e_td            = randn(n,1).*(X(:,2)+X(:,3));%form errors
    Y               = X*beta+e_td;                %form Y
    beta_hats2(:,i) = (X'*X)\(X'*Y);              %estimate beta_hat
    es              = Y-X*beta_hats2(:,i);        %form residuals
    sigma2          = ((es'*es)/(n-K));           %estimate var errors
    vcov_mat        = sigma2*((X'*X)\eye(K));     %calc vcov matrix
    se_hats2(:,i)   = sqrt(diag(vcov_mat));       %store std errors
end

%%
% Exercise 1: Part 2b)

t_stats2 = (beta_hats2-beta)./se_hats2; %calculate t-stats

%generate figure for beta_1
figure();
histogram(t_stats2(1,:),'Normalization','pdf');
title('Histogram of Standardized $\hat{\beta}_1$','interpreter','latex')
hold on

x = min(t_stats2(1,:)):0.1:max(t_stats2(1,:));
y = normpdf(x);
plot(x,y);

hold off


%generate figure for beta_1
figure();
histogram(t_stats2(2,:),'Normalization','pdf');
title('Histogram of Standardized $\hat{\beta}_2$','interpreter','latex')
hold on

x = min(t_stats2(2,:)):0.1:max(t_stats2(2,:));
y = normpdf(x);
plot(x,y);

hold off


%generate figure for beta_1
figure();
histogram(t_stats2(3,:),'Normalization','pdf');
title('Histogram of Standardized $\hat{\beta}_3$','interpreter','latex')
hold on

x = min(t_stats2(3,:)):0.1:max(t_stats2(3,:));
y = normpdf(x);
plot(x,y);

hold off

%%
% Exercise 1: Part 2c)
p_vals_c2    = 2*(1-tcdf(abs(t_stats2),n-K)); %calc pvalues
decisions_c2 = p_vals_c2<=alpha;              %accecpt or reject
sizes2       = mean(decisions_c2,2)           %calc sizes

%%
% Part 2d)
t_stats2d    = beta_hats2./se_hats2;           %calc new t-stats
p_vals_d2    = 2*(1-tcdf(abs(t_stats2d),n-K)); %calc pvalues
decisions_d2 = p_vals_d2<=alpha;               %accept or reject
power2       = mean(decisions_d2,2)            %calc power

