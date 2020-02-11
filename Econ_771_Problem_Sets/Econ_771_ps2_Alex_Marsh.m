%========================================================================%
%                        Title: Code for Econ 771 PS2                    %
%                        Author: Alex Marsh                              %
%                        Date: 01/28/2020                                %
%========================================================================%

%% Import and format data
% 
fpath = strcat('/Users/alexander/Documents',...
    '/School/Classes_Current',...
    '/ECON_771_Econometrics/Problem_Sets',...
    '/nerlove.xls');

Data  = readtable(fpath, "UseExcel", false);
TC    = Data.TC;    %store Total Cost
Q     = Data.Q;     %store Output
PL    = Data.PL;    %store Price of Labor
PF    = Data.PF;    %store Price of Fuel
PK    = Data.PK;    %store Price of Capital
N     = length(TC); %store Sample Size

TC_ln = log(TC); %form log(TC)
Q_ln  = log(Q);  %form log(Q)
PL_ln = log(PL); %form log(PL)
PF_ln = log(PF); %form log(PF)
PK_ln = log(PK); %form log(PK)

clear fpath Data;

%% Question b) estimate (1.7.4)
% 
X = [ones(N,1) Q_ln PL_ln PK_ln PF_ln]; %form X matrix
K = size(X,2);                          %store number of coefs
Y = TC_ln;                              %form Y vector

beta_hat_b = (X'*X)\(X'*Y);             %estimate beta
es_b       = Y-X*beta_hat_b;            %store residuals
sigma2_b   = (es_b'*es_b)/(N-K);        %estimate sigma^2
vcov_mat   = sigma2_b*((X'*X)\eye(K));  %form vcov matrix
se_hat_b   = sqrt(diag(vcov_mat));      %store standard errors
beta_hat_b
se_hat_b
%%
% The estimates are identical to those in the text. 

%% Question c) estimate (1.7.6)
% 

TCPF_ln = TC_ln - PF_ln; %form log(TC/PF)
PKPF_ln = PK_ln - PF_ln; %form log(PK/PF)
PLPF_ln = PL_ln - PF_ln; %form log(PL/PF)

X = [ones(N,1) Q_ln PLPF_ln PKPF_ln];  %form new X matrix
K = size(X,2);                         %store number of coefs
Y = TCPF_ln;                           %store Y vector

beta_hat_c = (X'*X)\(X'*Y);            %estimate beta
es_c       = Y-X*beta_hat_c;           %form residuals
sigma2_c   = (es_c'*es_c)/(N-K);       %estimate sigma^2
vcov_mat   = sigma2_c*((X'*X)\eye(K)); %estimate vcov matrix
se_hat_c   = sqrt(diag(vcov_mat));     %store standard errors
beta_hat_c
se_hat_c

%%
% The results are similar to those in the paper and identical to the ones
% in the text. The results in the paper can be found on page 176 in Table
% 3, especifically regression I.

%% Question d) 
% 
Xd1 = X(1:29,:);         %store X matrix for group 1
Yd1 = Y(1:29);           %store Y vector for group 1
Xd2 = X((1:29)+29*1,:);  %store X matrix for group 2
Yd2 = Y((1:29)+29*1);    %store Y vector for group 2
Xd3 = X((1:29)+29*2,:);  %store X matrix for group 3
Yd3 = Y((1:29)+29*2);    %store Y vector for group 3
Xd4 = X((1:29)+29*3,:);  %store X matrix for group 4
Yd4 = Y((1:29)+29*3);    %store Y vector for group 4
Xd5 = X((1:29)+29*4,:);  %store X matrix for group 5
Yd5 = Y((1:29)+29*4);    %store Y vector for group 5

beta_hat_d1 = (Xd1'*Xd1)\(Xd1'*Yd1); %estimate beta for group 1
beta_hat_d2 = (Xd2'*Xd2)\(Xd2'*Yd2); %estimate beta for group 2
beta_hat_d3 = (Xd3'*Xd3)\(Xd3'*Yd3); %estimate beta for group 3
beta_hat_d4 = (Xd4'*Xd4)\(Xd4'*Yd4); %estimate beta for group 4
beta_hat_d5 = (Xd5'*Xd5)\(Xd5'*Yd5); %estimate beta for group 5

%store all group beta_hats in one vector
beta_hat_d = [beta_hat_d1; beta_hat_d2; ...
    beta_hat_d3; beta_hat_d4; beta_hat_d5];

es_d1 = Yd1-Xd1*beta_hat_d1; %form residuals for group 1
es_d2 = Yd2-Xd2*beta_hat_d2; %form residuals for group 2
es_d3 = Yd3-Xd3*beta_hat_d3; %form residuals for group 3
es_d4 = Yd4-Xd4*beta_hat_d4; %form residuals for group 4
es_d5 = Yd5-Xd5*beta_hat_d5; %form residuals for group 5

sigma2_d1 = (es_d1'*es_d1)/(29-K); %estimate sigma^2 for group 1
sigma2_d2 = (es_d2'*es_d2)/(29-K); %estimate sigma^2 for group 2
sigma2_d3 = (es_d3'*es_d3)/(29-K); %estimate sigma^2 for group 3
sigma2_d4 = (es_d4'*es_d4)/(29-K); %estimate sigma^2 for group 4
sigma2_d5 = (es_d5'*es_d5)/(29-K); %estimate sigma^2 for group 5


%store estimates of returns to scale in one vector
RTS = 1./[beta_hat_d1(2) beta_hat_d2(2) beta_hat_d3(2)...
    beta_hat_d4(2) beta_hat_d5(2)]

%store estimates of sigma^2 in one vector
sigma2_d  = [sigma2_d1 sigma2_d2 sigma2_d3...
    sigma2_d4 sigma2_d5]


%%
% Both the returns to scale and the variance of the errors decrease as
% output increases.

%% Part e) 
% 
X_s = [X(1:29,:) zeros(29,4*4);
    zeros(29,4) X((1:29)+1*29,:) zeros(29,4*3);
    zeros(29,4*2) X((1:29)+2*29,:) zeros(29,4*2);
    zeros(29,4*3) X((1:29)+3*29,:) zeros(29,4*1);
    zeros(29,4*4) X((1:29)+4*29,:)]; % form stacked X matrix 

beta_hat_e = (X_s'*X_s)\(X_s'*Y);    %estimate beta 
es_e       = Y-X_s*beta_hat_e;       %form residuals

%calculate total SSR for part d
SSR_d = es_d1'*es_d1+es_d2'*es_d2+es_d3'*es_d3+...
    es_d4'*es_d4+es_d5'*es_d5;
SSR_e = es_e'*es_e; %calculate SSR for part e

[beta_hat_e beta_hat_d]
[SSR_d SSR_e]

%%
% It can be seen that beta_hat_d and beta_hat_e are the same as well as
% SSR_d and SSR_e. This is because Seemingly Unrelated Regressions and OLS
% equation-by-equation are equivalent (at least for point estimates).


%% Part f)
% 

SSR_R = es_c'*es_c;  %calcuclate SSR_R
SSR_U = es_e'*es_e;  %calcculate SSR_U
r     = 16;          %store number of restrictions
K     = size(X_s,2); %store number of parameters

F_stat_f = ((SSR_R-SSR_U)/r)/(SSR_U/(N-K)) %calculate F-stat
p_value_f = 1-fcdf(F_stat_f,r,N-K)         %calculate p-value

%%
% As can be seen above, there are 16 restrictions (intuitively, 16 equals
% signs). The F-statistic is $6.0678$ which can be rejected at every common
% significance level as the p-value is $1.0107e-09$; that is, we can
% reject the hypothesis that _all_ the coefficients are the same across
% groups. 
%% Part g)
%

X_s2 = [X(:,1).*(((1:145)<=29))' X(:,2).*(((1:145)<=29))' ...
    X(:,1).*((1:145)<=29*2&(1:145)>29*1)' ...
    X(:,2).*((1:145)<=29*2&(1:145)>29*1)'...
    X(:,1).*((1:145)<=29*3&(1:145)>29*2)' ...
    X(:,2).*((1:145)<=29*3&(1:145)>29*2)' ...
    X(:,1).*((1:145)<=29*4&(1:145)>29*3)'...
    X(:,2).*((1:145)<=29*4&(1:145)>29*3)'...
    X(:,1).*((1:145)<=29*5&(1:145)>29*4)' ...
    X(:,2).*((1:145)<=29*5&(1:145)>29*4)'...
    X(:,3:4)];

beta_hat_g = (X_s2'*X_s2)\(X_s2'*Y) %estimate Model 3
es_g = Y-X_s2*beta_hat_g;           %form residuals

SSR_R = es_g'*es_g; %calculate SSR_R
r = 8;              %store number of restrictions

F_stat_g = ((SSR_R-SSR_U)/r)/(SSR_U/(N-K)) %calculate F-stat
p_value_g = 1-fcdf(F_stat_g,r,N-K)         %calculate p-value

%%
% We cannot reject the null that the elasticities are equal across groups.
% This together with the rejected null in part f suggests that $\beta_1$
% (the intercept) and $\beta_2$ (which is $\frac{1}{r}$) do vary by group while the
% price elasticities do not vary by group.

%% Part h)
% 

X = [ones(N,1) Q_ln Q_ln.^2 PLPF_ln PKPF_ln]; %form X matrix
K = size(X,2);

beta_hat_h = (X'*X)\(X'*Y)       %estimate beta using OLS
es_h       = Y-X*beta_hat_h;     %form residuals

sz = 25;
figure();
scatter(Q_ln,es_h,sz,'filled');
title('OLS residuals versus log(Q)')
xlabel('log(Q)')
ylabel('residuals')

W = 1./((0.0565+2.1377./Q));          %form W vector
W = diag(W);                          %diagonalize W

beta_hat_h_w = (X'*W*X)\(X'*W*Y)      %estimate beta using WLS
es_h_w     = Y-X*beta_hat_h_w;        %form residuals

figure();
scatter(Q_ln,es_h_w,sz,'filled');
title('WLS residuals versus log(Q)')
xlabel('log(Q)')
ylabel('residuals')
