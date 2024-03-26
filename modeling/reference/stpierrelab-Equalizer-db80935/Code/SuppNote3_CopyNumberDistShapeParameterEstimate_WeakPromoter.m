% This script generates data fpr Supplementary Note 3. 
% Estimate the standard deviation of the log normal plasmid copy number distribution using open loop data
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

% Read single-cell data
Array1=load('..\ExperimentData\SuppNote3_PGK_SingleCellExpression.mat');

OpenLoop = Array1.single_cell_A9_20190419;

OpenLoop=OpenLoop(OpenLoop>0);
pHat1 = gamfit(OpenLoop);

x = 0:10:max(OpenLoop);
k = pHat1(1)
theta = pHat1(2);
y = gamcdf(x,k,theta);
% plot distribution of the openloop expression level
figure(1)
histogram(OpenLoop,'Normalization','cdf')
xlabel('Expression level (a.u.)')
ylabel('CDF')
title('Distribution of expression level of the PGK openloop circuit')
hold on
plot(x,y,'LineWidth',3)
hold off
legend('Unregulated PGK circuit expression level CDF','Fitted gamma CDF')

x2 = 1:100;
mu = 2.63;
sigma = pHat1(2);
y2 = lognpdf(x2,mu,sigma);

X = sprintf('The shape parameter of plasmid copy number distribution is %d ',pHat1(1));
disp(X)
