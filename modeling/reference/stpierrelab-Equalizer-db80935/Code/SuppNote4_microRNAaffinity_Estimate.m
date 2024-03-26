% This script estimates RISC concentration and microRNA affinity. 
% Data used in Supplementary Note 4. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

% Initial guesses for RISC concentration and miRNA
RISC = 1.7E5;                   % unit: molecule/cell
miRNA0 = [0.0184 0.184 1.84];   % unit: 1/(molecule*second)
% pairs = combvec(RISC0,miRNA0);

x_fmin = zeros(length(miRNA0),1);
f_eval = zeros(length(miRNA0),1);

k = 0.5716;         % fitted from CopyNumberDistribution_StandardDeviationEstimate.m
theta = 120;        % fitted from CopyNumberDistribution_MeanEsitmate.m
celltot = 10000;    % Assume a total of 10000 cells
copynumber_cell = round(gamrnd(k,theta,celltot,1)); % copynumber in each cell
copynumber_cell = copynumber_cell(copynumber_cell>0);
copynumber_cell_sorted = sort(copynumber_cell);
copynumber = unique(copynumber_cell_sorted);

for i = 1:length(miRNA0)
    pair0 = miRNA0(i);
    options = optimset('Display','iter','PlotFcns',@optimplotfval);
    cal_diff = @(x)expressionlevel(x,celltot,copynumber,copynumber_cell_sorted,RISC);
    [x_fmin(i),f_eval(i)] = fminsearch(cal_diff,pair0,options);
end

save('../SimulationData/SupplementaryNote4_microRNAaffinity_estimate_theta120.mat','miRNA0','x_fmin','f_eval')


function f=expressionlevel(pair,celltot,copynumber,copynumber_cell_sorted,RISC)
% values cannot be negative

miRNAaffinity = max(0,pair)

RNA_deg = 2.88e-4;

sbioloadproject ../Models/Equalizer_model m1
csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1e7;

extracell_inducer = [0 1 5 10 50 100];%in the unit of ng/m
dox_influx = 0.156.*extracell_inducer;

%vector to hold COV at each induction level
EQ_Meanexp = zeros(1, length(extracell_inducer));
expression_vector=zeros(celltot,1);

POI = zeros(length(copynumber),length(extracell_inducer));
count = 0;

for i = 1:length(copynumber)
    for j = 1:length(extracell_inducer)
        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});
        
        set(m1.Reaction(11).KineticLaw.Parameters,'Value',dox_influx(j));
        set(m1.species(10),'InitialAmount',dox_influx(j)/3.33e-4);
    
        set(m1.species(13),'InitialAmount',RISC); %set RISC concentration
        set(m1.Reaction(13).KineticLaw.Parameters,'Value',miRNAaffinity); %set miRNA affinity
        
        set(m1.Reaction(12).KineticLaw.Parameters,'Value',RNA_deg); %set mRNA degradation rate
        set(m1.Reaction(20).KineticLaw.Parameters,'Value',RNA_deg); %set mRNA-microRNA degradation rate
        
        simdata = sbiosimulate(m1,csObj,vObj1);
        [~, stateData] = selectbyname(simdata, 'Cell.POI');
        POI(i,j) = stateData(end);
        count = count + 1;
    end
end

for i = 1:length(extracell_inducer)
    explvl = POI(:,i)';
    for j = 1:length(copynumber_cell_sorted)
        if copynumber_cell_sorted(j)==0
            expression_vector(j)=0;
        else
            expression_vector(j) = explvl(copynumber==copynumber_cell_sorted(j));
        end
    end
    EQ_Meanexp(i) = mean(expression_vector);
end

sbioloadproject ../Models/NF_model m1
csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1e7;

%vector to hold COV at each induction level
NF_Meanexp = zeros(1, length(extracell_inducer));
expression_vector=zeros(celltot,1);

POI = zeros(length(copynumber),length(extracell_inducer));
count = 0;

for i = 1:length(copynumber)
    for j = 1:length(extracell_inducer)
        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});
        
        set(m1.Reaction(6).KineticLaw.Parameters,'Value',dox_influx(j));
        set(m1.species(5),'InitialAmount',dox_influx(j)/3.33e-4);
        
        set(m1.Reaction(7).KineticLaw.Parameters,'Value',RNA_deg); %set mRNA degradation rate
        
        simdata = sbiosimulate(m1,csObj,vObj1);
        [~, stateData] = selectbyname(simdata, 'Cell.POI');
        POI(i,j) = stateData(end);
        count = count + 1;
    end
end

for i = 1:length(extracell_inducer)
    explvl = POI(:,i)';
    for j = 1:length(copynumber_cell_sorted)
        if copynumber_cell_sorted(j)==0
            expression_vector(j)=0;
        else
            expression_vector(j) = explvl(copynumber==copynumber_cell_sorted(j));
        end
    end
    NF_Meanexp(i) = mean(expression_vector);
end

%Load experimental Equalizer and NF mean expression at the same [dox]
% load('../ExperimentData/SuppNote4_Equalizer_expression.mat','avg')
% EQ_expmean = mean(avg);
% load('../ExperimentData/SuppNote4_NF_expression.mat','avg')
% NF_expmean = mean(avg);
load('../ExperimentData/SuppNote4_Equalizer_expression_Fig4a.mat','avg')
EQ_expmean = mean(avg);
load('../ExperimentData/SuppNote4_NF_expression_Fig4a.mat.mat','avg')
NF_expmean = mean(avg);

exp_ratio = NF_expmean./EQ_expmean;
simulation_ratio = NF_Meanexp./EQ_Meanexp;

f = sumsqr(exp_ratio-simulation_ratio);
end
