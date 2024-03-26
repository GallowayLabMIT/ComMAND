% This script generates Supplementary Figure 11. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

miRNAdissociation = [0.303];
RISC = 1.7e5;
count=0;

extracellular_inducer = [0 0.1 0.5 1 5 10 50 100]; % in the unit of ng/m
legend_list = {'Equalizer','IFF', 'NF', 'Unregulated'};

copynumber_plot = 1:5:1000;
    
POI_Eq = zeros(length(copynumber_plot),length(miRNAdissociation));

%Numerical differentiation to get local log sensitivity
dPOIdCN_Eq = zeros(size(POI_Eq));

dox_influx = 0.156.*extracellular_inducer(4);

for i = 1:length(copynumber_plot)
    for j= 1:length(miRNAdissociation)
        sbioloadproject ../Models/Equalizer_model m1
        csObj = m1.addconfigset('newStopTimeConfigSet');
        csObj.StopTime = 1e7;
        
        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber_plot(i)});

        vObj = [vObj1];

        set(m1.Reaction(11).KineticLaw.Parameters,'Value',dox_influx);
        set(m1.Reaction(13).KineticLaw.Parameters,'Value',miRNAdissociation(j));
        set(m1.species(10),'InitialAmount',dox_influx/3.33e-4);
        set(m1.species(13),'InitialAmount',RISC);

        simdata = sbiosimulate(m1,csObj,vObj);
        [~, stateData] = selectbyname(simdata, 'Cell.POI');
        POI_Eq(i,j) = stateData(end);
    end
end

for i = 1:length(copynumber_plot)
    for j = 1:length(miRNAdissociation)
        if i == 1
            dPOIdCN_Eq(i,j) = (-POI_Eq(i+2,j)+4*POI_Eq(i+1,j)-3*POI_Eq(i,j))/2*(i/POI_Eq(i,j));
        elseif i == length(copynumber_plot)
            dPOIdCN_Eq(i,j) = (3*POI_Eq(i,j)-4*POI_Eq(i-1,j)+POI_Eq(i-2,j))/2*(i/POI_Eq(i,j));
        else 
            dPOIdCN_Eq(i,j) = (POI_Eq(i+1,j)-POI_Eq(i-1,j))/2*(i/POI_Eq(i,j));
       end
    end
end

figure(5)
dosage_compensation_Eq = 1./dPOIdCN_Eq(2:end);
plot(copynumber_plot(2:end),dosage_compensation_Eq,'b-','LineWidth',3)
hold on

POI_IFF = zeros(length(copynumber_plot),length(miRNAdissociation));

for i = 1:length(copynumber_plot)
    for j= 1:length(miRNAdissociation)

        sbioloadproject ../Models/IFF_model  m1
        csObj = m1.addconfigset('newStopTimeConfigSet');
        csObj.StopTime = 1e7;

        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber_plot(i)});

        vObj = [vObj1];

        set(m1.Reaction(10).KineticLaw.Parameters,'Value',miRNAdissociation(j));
        set(m1.species(10),'InitialAmount',RISC);

        simdata = sbiosimulate(m1,csObj,vObj);
        [time, stateData] = selectbyname(simdata, 'Cell.POI');
        POI_IFF(i,j) = stateData(end);
        count = count + 1;
    end
end


%Numerical differentiation to get local log sensitivity
dPOIdCN_IFF = zeros(size(POI_IFF));

for i = 1:length(copynumber_plot)
    for j = 1:length(miRNAdissociation)
        if i == 1
            dPOIdCN_IFF(i,j) = (-POI_IFF(i+2,j)+4*POI_IFF(i+1,j)-3*POI_IFF(i))/2*(i/POI_IFF(i,j));
        elseif i == length(copynumber_plot)
            dPOIdCN_IFF(i,j) = (3*POI_IFF(i,j)-4*POI_IFF(i-1,j)+POI_IFF(i-2,j))/2*(i/POI_IFF(i,j));
        else
            dPOIdCN_IFF(i,j) = (POI_IFF(i+1,j)-POI_IFF(i-1,j))/2*(i/POI_IFF(i,j));
        end
    end
end
dosage_compensation_IFF = 1./dPOIdCN_IFF(2:end);
plot(copynumber_plot(2:end),dosage_compensation_IFF,'g-','LineWidth',3)

POI = zeros(length(copynumber_plot),1);
dPOIdCN_NF = zeros(size(POI));

extracellular_inducer = extracellular_inducer(6);
dox_influx = 0.156.*extracellular_inducer;

for i = 1:length(copynumber_plot)

    sbioloadproject ../Models/NF_model  m1
    csObj = m1.addconfigset('newStopTimeConfigSet');
    csObj.StopTime = 1e7;

    namevObj1 = strcat('v1_',num2str(count));
    vObj1 = addvariant(m1,namevObj1);
    addcontent(vObj1,{'species','gene','InitialAmount',copynumber_plot(i)});

    set(m1.Reaction(6).KineticLaw.Parameters,'Value',dox_influx);
    set(m1.species(5),'InitialAmount',dox_influx/3.33e-4);

    simdata = sbiosimulate(m1,csObj,vObj1);
    [time, stateData] = selectbyname(simdata, 'Cell.POI');
    POI(i) = stateData(end);
    count = count + 1;

end


%Numerical differentiation to get local log sensitivity
for j = 1:length(copynumber_plot)
    if j == 1 
        dPOIdCN_NF(j) = (-POI(j+2)+4*POI(j+1)-3*POI(j))/2*(j/POI(j));
    elseif j == length(copynumber_plot)
        dPOIdCN_NF(j) = (3*POI(j)-4*POI(j-1)+POI(j-2))/2*(j/POI(j));
    else 
        dPOIdCN_NF(j) = (POI(j+1)-POI(j-1))/2*(j/POI(j));
    end
end

dosage_compensation_NF = 1./dPOIdCN_NF(2:end);
plot(copynumber_plot(2:end),dosage_compensation_NF,'c-','LineWidth',3)

POI = zeros(length(copynumber_plot),1);

for i = 1:length(copynumber_plot)

    sbioloadproject ../Models/Constitutive_model  m1
    csObj = m1.addconfigset('newStopTimeConfigSet');
    csObj.StopTime = 1e7;

    namevObj1 = strcat('v1_',num2str(count));
    vObj1 = addvariant(m1,namevObj1);

    vObj = [vObj1];

    set(m1.species(1),'InitialAmount',copynumber_plot(i));
    
    simdata = sbiosimulate(m1,csObj,vObj);
    [time, stateData] = selectbyname(simdata, 'Cell.POI');
    POI(i) = stateData(end);
    OL_POI = stateData;
    count = count + 1;
end

%Numerical differentiation to get local log sensitivity
dPOIdCN = zeros(size(POI));

for j = 1:length(copynumber_plot)
    if j == 1
        dPOIdCN(j) = (-POI(j+2)+4*POI(j+1)-3*POI(j))/2*(j/POI(j));
    elseif j == length(copynumber_plot)
        dPOIdCN(j) = (3*POI(j)-4*POI(j-1)+POI(j-2))/2*(j/POI(j));
    else 
        dPOIdCN(j) = (POI(j+1)-POI(j-1))/2*(j/POI(j));
    end
end
    
dosage_compensation_openloop = 1./dPOIdCN(2:end);
plot(copynumber_plot(2:end),dosage_compensation_openloop,'k-','LineWidth',3)

ax = gca;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 18;
ax = gca;
xlabel('Plasmid copy number')
ylabel('Predicted gene dosage compensation')

ylim([0 30])
xlim([1 500])
yticks([1 10 20 30 40])
xticks([1 100 200 300 400 500])
box off

legend(legend_list,'Location','northeast','FontSize', 20); 
legend box off 
hold off
save('../SimulationData/SuppFig11_logsensitivity_miRNA0.303_dox1.mat','copynumber_plot','miRNAdissociation','copynumber','dosage_compensation_Eq','dosage_compensation_NF','dosage_compensation_IFF','dosage_compensation_openloop')
