% This script generates Figure 4 B/C. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

miRNAdissociation = 0.303;
RISC = 1.7e5;

celltot = 10000; % Assume a total of 10000 cells
copynumber_cell = round(gamrnd(0.5716,120,celltot,1)); % copynumber in each cell,std of 1.522 fit from PGK A9 well on 20190419
copynumber_cell = copynumber_cell(copynumber_cell>0);
copynumber_cell_sorted = sort(copynumber_cell);
expression_vector=zeros(celltot,1);

sbioloadproject ../Models/Equalizer_model m1
csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1e7;

extracell_inducer = [0 0.1 0.5 1 5 10 50 100];  % in the unit of ng/ml
dox_influx = 0.156.*extracell_inducer;

copynumber = unique(copynumber_cell_sorted);
POI = zeros(length(copynumber),length(extracell_inducer));

Eq_CV = zeros(length(extracell_inducer),1);
Eq_mean = zeros(length(extracell_inducer),1);
count = 0;
for i = 1:length(copynumber)
    for j = 1:length(extracell_inducer)
        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});
        
        set(m1.Reaction(11).KineticLaw.Parameters,'Value',dox_influx(j));
        set(m1.Reaction(13).KineticLaw.Parameters,'Value',miRNAdissociation);
        set(m1.species(10),'InitialAmount',dox_influx(j)/3.33e-4);
        set(m1.species(13),'InitialAmount',RISC);
        
        simdata = sbiosimulate(m1,csObj,vObj1);
        [~, stateData] = selectbyname(simdata, 'Cell.POI');
        POI(i,j) = stateData(end);
        count = count + 1;
    end
end

for i = 1:length(extracell_inducer)
    explvl = POI(:,i)';
    for j = 1:length(copynumber_cell)
        if copynumber_cell(j)==0
            expression_vector(j)=0;
        else
            expression_vector(j) = explvl(copynumber==copynumber_cell(j));
        end
    end
    Eq_CV(i) = std(expression_vector(:))/mean(expression_vector)*100;
    Eq_mean(i) = mean(expression_vector);
end


sbioloadproject ../Models/NF_model m1
csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1e7;

extracell_inducer_NF = [0 1 5 10 50 100 500 1000];%in the unit of ng/ml
dox_influx = 0.156.*extracell_inducer_NF;

copynumber = unique(copynumber_cell_sorted);
POI = zeros(length(copynumber),length(extracell_inducer_NF));

NF_CV = zeros(length(extracell_inducer_NF),1);
NF_mean = zeros(length(extracell_inducer_NF),1);
count = 0;

for i = 1:length(copynumber)
    for j = 1:length(extracell_inducer_NF)
        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});
        
        set(m1.Reaction(6).KineticLaw.Parameters,'Value',dox_influx(j));
        set(m1.species(5),'InitialAmount',dox_influx(j)/3.33e-4);
        
        simdata = sbiosimulate(m1,csObj,vObj1);
        [~, stateData] = selectbyname(simdata, 'Cell.POI');
        POI(i,j) = stateData(end);
        count = count + 1;
    end
end

for i = 1:length(extracell_inducer_NF)
    explvl = POI(:,i)';
    for j = 1:length(copynumber_cell)
        if copynumber_cell(j)==0
            expression_vector(j)=0;
        else
            expression_vector(j) = explvl(copynumber==copynumber_cell(j));
        end
    end
    NF_CV(i) = std(expression_vector(:))/mean(expression_vector)*100;
    NF_mean(i) = mean(expression_vector);
end


sbioloadproject ../Models/IFF_model m1
csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1e7;
RNA_deg = 2.88e-4;

POI = zeros(length(copynumber),1);

count = 0;

for i = 1:length(copynumber)

    namevObj1 = strcat('v1_',num2str(count));
    vObj1 = addvariant(m1,namevObj1);
    addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});
    set(m1.Reaction(10).KineticLaw.Parameters,'Value',miRNAdissociation);
%     set(m1.Reaction(13).KineticLaw.Parameters,'Value',0);
    set(m1.species(10),'InitialAmount',RISC);
    simdata = sbiosimulate(m1,csObj,vObj1);
    [~, stateData] = selectbyname(simdata, 'Cell.POI');
    POI(i) = stateData(end);
    count = count + 1;
  
end


explvl = POI';
for j = 1:length(copynumber_cell)
    if copynumber_cell(j)==0
        expression_vector(j)=0;
    else
        expression_vector(j) = explvl(copynumber==copynumber_cell(j));
    end
end
IFF_CV = std(expression_vector)/mean(expression_vector)*100;
IFF_mean = mean(expression_vector);

% Load experiment data
Array1=load('../ExperimentData/Fig4b_Eq_CV.mat');
Eq_CV_exp = Array1.cov;
Eq_yerror = std(Eq_CV_exp,1)/sqrt(size(Eq_CV_exp,1));
Eq_CV_exp = mean(Eq_CV_exp,1);
Array1=load('../ExperimentData/Fig4b_NF_CV.mat');
NF_CV_exp = Array1.cov;
NF_CV_std = std(NF_CV_exp,1);
NF_CV_std = NF_CV_std(1,(1:8));
NF_yerror = NF_CV_std/sqrt(size(NF_CV_exp,1));
NF_CV_exp = mean(NF_CV_exp,1);
NF_CV_exp = NF_CV_exp(1,(1:8));
Array1=load('../ExperimentData/Fig4c_IFF_CV.mat');
IFF_CV_exp = Array1.cov;
IFF_yerror = std(IFF_CV_exp,1)/sqrt(size(IFF_CV_exp,1));
IFF_CV_exp = mean(IFF_CV_exp,1);

figure(14)
extracell_inducer = [0.01 0.1 0.5 1 5 10 50 100];       % in the unit of ng/ml
extracell_inducer_NF = [0.01 1 5 10 50 100 500 1000];   % in the unit of ng/ml
plot(extracell_inducer,(Eq_CV.^2+69.6333^2).^0.5,'--bs','LineWidth',4,'MarkerSize',15,'MarkerEdgeColor', 'b')
hold on 
plot(extracell_inducer_NF,(NF_CV.^2+71.4333^2).^0.5,'--rs','LineWidth',4,'MarkerSize',15,'MarkerEdgeColor', 'r')
e=errorbar(extracell_inducer_NF,NF_CV_exp,NF_yerror,'--s','CapSize',18,'MarkerSize',15,'MarkerEdgeColor', 'r','MarkerFaceColor','r','LineWidth',4);
e.Color = 'k';
e=errorbar(extracell_inducer,Eq_CV_exp,Eq_yerror,'--s','CapSize',18,'MarkerSize',15,'MarkerEdgeColor', 'b','MarkerFaceColor','b','LineWidth',4);
e.Color = 'k';
xt = get(gca, 'XTick');
ax = gca;
ax.YAxis.Color = 'k';
ax.YAxis.LineWidth = 2;
ax.YAxis.FontSize = 20;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 20;
ax.TickLength = [0.025, 0.025];
ax.XAxis.FontSize = 18;
ax.Title.FontSize = 20;
ylim([70 180])
xlim([5e-3 1000])
xlabel('[doxycycline] (ng/mL)')
ylabel('Cell-to-cell variation (%CV)')
title('Total cell-to-cell variation')
legend_list = {'Equalizer\_simulation','Equalizer\_experiment'};
legend(legend_list,'Location','northeast','FontSize', 16); 
legend box off 
set(gca, 'XScale', 'log')
pbaspect([1 1 1])
box off
hold off


figure(15)
X = categorical({'IFF experiment', 'IFF simulation'});
X = reordercats(X,{'IFF experiment', 'IFF simulation'});
Y = [IFF_CV_exp,(IFF_CV^2+71.4333^2)^0.5];
bar(X,Y, 0.2)
hold on
er = errorbar(X,Y,[IFF_yerror 0],[IFF_yerror 0]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 3;  
er.CapSize = 10;  
xt = get(gca, 'XTick');
ax = gca;
ax.YAxis.Color = 'k';
ax.YAxis.LineWidth = 2;
ax.YAxis.FontSize = 18;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 18;
ax.Title.FontSize = 20;
ylim([70 180])
ylabel('Cell-to-cell variation (%CV)')
title('Total cell-to-cell variation')
box off
hold off
save('../SimulationData/Fig4bc_Predicted_Total_CV_miRNA_0.303_theta_120_NF.mat','Eq_CV','extracell_inducer','NF_CV','IFF_CV','extracell_inducer_NF')
