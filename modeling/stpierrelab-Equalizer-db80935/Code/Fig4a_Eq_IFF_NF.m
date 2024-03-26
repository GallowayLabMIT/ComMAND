% This script generates Figure 4 A. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

color = {'b' 'k' 'k' };
mark = {'s','s','s'};
construct = {'CMV Openloop plasmid' 'Equalizer plasmid' 'NF cell line' 'PGK Openloop plasmid' 'UBC Openloop plasmid' 'CMV cell line'};

parentFdr = '..\ExperimentData\';

% Load data
load('..\ExperimentData\Fig4a_Eq_IFF_NF.mat')

% Find the Eq data point with the lowest CV
[~,Eq_min_index]=min(nanmean(Eq_cov,1));
Eq_exp_minCV = mean(Eq_avg(:,Eq_min_index));
Eq_COV = Eq_cov(:,Eq_min_index);

[~,	NF_min_index]=min(nanmean(NF_cov,1));
NF_COV = NF_cov(:,NF_min_index);

IFF_COV = IFF_cov;

% Normalize the expression to the Eq mean expression with the lowest CV
Eq_avg = Eq_avg./Eq_exp_minCV;
IFF_avg = IFF_avg./Eq_exp_minCV;
NF_avg = NF_avg./Eq_exp_minCV;

Eq_xerror = nanstd(Eq_avg)/sqrt(size(Eq_avg,1));
Eq_yerror = nanstd(Eq_cov)/sqrt(size(Eq_cov,1));

NF_xerror = nanstd(NF_avg)/sqrt(size(NF_avg,1));
NF_yerror = nanstd(NF_cov)/sqrt(size(NF_cov,1));

IFF_xerror = nanstd(IFF_avg)/sqrt(size(IFF_avg,1));
IFF_yerror = nanstd(IFF_cov)/sqrt(size(IFF_cov,1));

% Sort data point based on expression level for Eq and NF
Eq_avg = mean(Eq_avg,1);
Eq_cov = mean(Eq_cov,1);
NF_avg = mean(NF_avg,1);
NF_cov = mean(NF_cov,1);
IFF_avg = mean(IFF_avg,1);
IFF_cov = mean(IFF_cov,1);
[Eq_avg_sorted, Eq_order]  = sort(Eq_avg);
Eq_cov_sorted = Eq_cov(Eq_order);
Eq_xerror = Eq_xerror(Eq_order);
Eq_yerror = Eq_yerror(Eq_order);

[NF_avg_sorted, NF_order]  = sort(NF_avg);
NF_cov_sorted = NF_cov(NF_order);
NF_xerror = NF_xerror(NF_order);
NF_yerror = NF_yerror(NF_order);

figure(1)

e=errorbar(Eq_avg_sorted,Eq_cov_sorted,Eq_yerror,Eq_yerror,Eq_xerror,Eq_xerror,'--s','CapSize',18,'MarkerSize',15,'MarkerEdgeColor', 'b','MarkerFaceColor','b','LineWidth',4);
e.Color = 'k';
hold on

e=errorbar(NF_avg_sorted,NF_cov_sorted,NF_yerror,NF_yerror,NF_xerror,NF_xerror,'--s','CapSize',18,'MarkerSize',15,'MarkerEdgeColor', 'g','MarkerFaceColor','g','LineWidth',4);
e.Color = 'k';

e=errorbar(IFF_avg,IFF_cov,IFF_yerror,IFF_yerror,IFF_xerror,IFF_xerror,'--s','CapSize',18,'MarkerSize',15,'MarkerEdgeColor', 'm','MarkerFaceColor','m','LineWidth',4);
e.Color = 'k';
hold off

set(gca, 'XScale', 'log')
ylim([70 180])
xlim([0.5 100])
hold off
ax = gca;
ax.YAxis.Color = 'k';
ax.YAxis.LineWidth = 2;
ax.YAxis.FontSize = 20;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 20;
ax.TickLength = [0.025, 0.025];
pbaspect([1 1 1])
set(gca,'box','off')
hold off

% stats
[H1, pValue1, W1] = swtest(Eq_COV, 0.05)
[H2, pValue2, W2] = swtest(NF_COV, 0.05)
[H3, pValue3, W3] = swtest(IFF_COV, 0.05)
[h1,p_1] = vartest2(Eq_COV,NF_COV)
[h2,p_2] = vartest2(Eq_COV,IFF_COV)
p1=ranksum(Eq_COV,NF_COV,'Tail','left','method','exact')
p2=ranksum(Eq_COV,IFF_COV,'Tail','left','method','exact')


% two-way ANOVA followed by post-hoc Bonferonni
figure(2)
x = [Eq_COV, NF_COV,IFF_COV];
[p,tbl,stats] = anova2(x,3,'off')
c = multcompare(stats,'CType','bonferroni')

ax = gca;
ax.YAxis.Color = 'k';
ax.YAxis.LineWidth = 2;
ax.YAxis.FontSize = 20;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 20;
ax.TickLength = [0.02, 0.02];
pbaspect([1 1 1])
hold off