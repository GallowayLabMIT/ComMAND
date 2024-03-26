% This script generates Supplementary Figure 15. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

figure(1)
X = categorical({'Equalizer-L', 'Equalizer-L\_multipromoter'});
X = reordercats(X,{'Equalizer-L', 'Equalizer-L\_multipromoter'});
Y = [0.066509232	0.153819652].*100;
bar(X,Y, 0.2)
hold on
xt = get(gca, 'XTick');
ax = gca;
ax.YAxis.Color = 'k';
ax.YAxis.LineWidth = 2;
ax.YAxis.FontSize = 18;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 18;
ax.Title.FontSize = 20;
ylim([0 40])
pbaspect([1 1.4 1])
ylabel('Intrinsic Noise(%CV)')
title('Set 1: original parameters')
box off
hold off

figure(2)
X = categorical({'Equalizer-L', 'Equalizer-L\_multipromoter'});
X = reordercats(X,{'Equalizer-L', 'Equalizer-L\_multipromoter'});
Y = [0.152627241	0.302705985].*100;
bar(X,Y, 0.2)
hold on
xt = get(gca, 'XTick');
ax = gca;
ax.YAxis.Color = 'k';
ax.YAxis.LineWidth = 2;
ax.YAxis.FontSize = 18;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 18;
ax.Title.FontSize = 20;
ylim([0 40])
pbaspect([1 1.4 1])
ylabel('Intrinsic Noise(%CV)')
title('Set 2: amplified translational parameters')
box off
hold off

figure(3)
X = categorical({'Equalizer-L', 'Equalizer-L\_multipromoter'});
X = reordercats(X,{'Equalizer-L', 'Equalizer-L\_multipromoter'});
Y = [0.131650382	0.315014774].*100;
bar(X,Y, 0.2)
hold on
xt = get(gca, 'XTick');
ax = gca;
ax.YAxis.Color = 'k';
ax.YAxis.LineWidth = 2;
ax.YAxis.FontSize = 18;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 18;
ax.Title.FontSize = 20;
ylim([0 40])
pbaspect([1 1.4 1])
ylabel('Intrinsic Noise(%CV)')
title('Set 3: amplified transcriptional parameters')
box off
hold off