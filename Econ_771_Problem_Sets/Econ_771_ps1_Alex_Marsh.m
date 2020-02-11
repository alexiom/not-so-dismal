%===================================================%
%           Title: Code for Econ 771 PS1            %
%           Author: Alex Marsh                      %
%           Date: 01/20/2020                        %
%===================================================%

%%
% import and format data
fpath = strcat('/Users/alexander/Documents',...
    '/School/Classes_Current',...
    '/ECON_771_Econometrics/Problem_Sets',...
    '/ps1.csv');
Data = readtable(fpath,'ReadVariableNames',true);
Y  = Data.Y;
X1 = Data.X1;
X2 = Data.X2;
N  = length(Y);

X = [ones(N,1) X1 X2];
clear fpath Data;

%%
% make initial plot
sz = 25;
subplot(1,2,1)
scatter(X1,Y,sz)
title('X1 versus Y')
xlabel('X1')
ylabel('Y')

subplot(1,2,2)
scatter(X2,Y,sz)
title('X2 versus Y')
xlabel('X2')
ylabel('Y')

%%
% estimate b_hat
beta_hat = (X'*X)\(X'*Y)

%%
% calculate R^2
es = Y - X*beta_hat;
M0 = eye(N) - ones(N)/N;
R2 = 1-(es'*es)/(Y'*M0*Y)

%%
% add y_hat to plot
Y_hat = X*beta_hat;

subplot(1,2,1)
scatter(X1,Y,sz,'filled');
title('X1 versus Y')
xlabel('X1')
ylabel('Y')

hold on 
scatter(X1,Y_hat,sz,'filled');
legend('Y','Y\_hat','Location','best')

subplot(1,2,2)
scatter(X2,Y,sz,'filled');
title('X2 versus Y')
xlabel('X2')
ylabel('Y')

hold on 
scatter(X2,Y_hat,sz,'filled');
legend('Y','Y\_hat','Location','best')
hold off
