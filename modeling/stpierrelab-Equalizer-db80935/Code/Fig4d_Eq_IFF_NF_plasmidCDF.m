% This script generates Figure 4 D. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

extracell_inducer = [1];
dox_influx = 0.156.*extracell_inducer;

copynumber = 0:1:1000;
POI_Eq = zeros(length(copynumber),length(extracell_inducer));
RISC_Eq = zeros(length(copynumber),length(extracell_inducer));
TetR_Eq = zeros(length(copynumber),length(extracell_inducer));
POI_Eq_2 = zeros(length(copynumber),length(extracell_inducer));
RISC_Eq_2 = zeros(length(copynumber),length(extracell_inducer));
TetR_Eq_2 = zeros(length(copynumber),length(extracell_inducer));
count = 0;

for i = 1:length(copynumber)
    for k = 1:length(extracell_inducer)
        sbioloadproject ../Models/Equalizer_model m1
        csObj = m1.addconfigset('newStopTimeConfigSet');
        csObj.StopTime = 1e7;
        
        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});

        vObj = [vObj1];

        set(m1.Reaction(11).KineticLaw.Parameters,'Value',dox_influx(k));
        set(m1.species(10),'InitialAmount',dox_influx(k)/3.33e-4);

        simdata = sbiosimulate(m1,csObj,vObj);
        [~, stateData] = selectbyname(simdata, 'Cell.POI');
        POI_Eq(i,k) = stateData(end);
        [~, stateData] = selectbyname(simdata, 'Cell.RISC');
        RISC_Eq(i,k) = stateData(end);
    end
end


sbioloadproject ../Models/NF_model m1

csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1e7;

extracell_inducer = 10;
dox_influx = 0.156.*extracell_inducer;


POI_NF = zeros(length(copynumber),length(extracell_inducer));
TetR_NF = zeros(length(copynumber),length(extracell_inducer));
count = 0;

for i = 1:length(copynumber)
    for j = 1:length(extracell_inducer)
        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});
        
        namevObj2 = strcat('v2_',num2str(count));
        vObj2 = addvariant(m1,namevObj2);
        addcontent(vObj2,{'species','Inducer','InitialAmount',dox_influx(j)/3.33e-4});
        
        vObj = [vObj1 vObj2];
        
        set(m1.Reaction(6).KineticLaw.Parameters,'Value',dox_influx(j));

        simdata = sbiosimulate(m1,csObj,vObj);
        [~, stateData] = selectbyname(simdata, 'Cell.POI');
        POI_NF(i,j) = stateData(end);
        [~, stateData] = selectbyname(simdata, 'Cell.TetR');
        TetR_NF(i,j) = stateData(end);
        count = count + 1;
    end
end

dPOIdCN = zeros(size(POI_NF));
for i = 1:length(extracell_inducer)
    for j = 1:length(copynumber)
        if j == 1
            dPOIdCN(j,i) = (-POI_NF(j+2,i)+4*POI_NF(j+1,i)-3*POI_NF(j,i))/2*(j/POI_NF(j,i));
        elseif j == length(copynumber)
            dPOIdCN(j,i) = (3*POI_NF(j,i)-4*POI_NF(j-1,i)+POI_NF(j-2,i))/2*(j/POI_NF(j,i));
        else
            dPOIdCN(j,i) = (POI_NF(j+1,i)-POI_NF(j-1,i))/2*(j/POI_NF(j,i));
        end
    end
end

sbioloadproject ../Models/IFF_model m1
csObj = m1.addconfigset('newStopTimeConfigSet');

csObj.StopTime = 1e7;

extracell_inducer = 10;
dox_influx = 0.156.*extracell_inducer;

miRNAaffinity = [0.303];
RISC = 1.7e5;

POI_IFF = zeros(length(copynumber),length(miRNAaffinity));
RISC_IFF = zeros(length(copynumber),length(miRNAaffinity));
count = 0;

for i = 1:length(copynumber)
    for k = 1:length(miRNAaffinity)

        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});
        set(m1.Reaction(10).KineticLaw.Parameters,'Value',miRNAaffinity(k));
        set(m1.species(10),'InitialAmount',RISC);
        simdata = sbiosimulate(m1,csObj,vObj1);
        [~, stateData] = selectbyname(simdata, 'Cell.POI');
        POI_IFF(i,k) = stateData(end);
        [~, stateData] = selectbyname(simdata, 'Cell.RISC');
        RISC_IFF(i,k) = stateData(end);
        count = count + 1;
    end
end



sbioloadproject ../Models/Constitutive_model m1

csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1e7;

POI_openloop = zeros(1,length(copynumber));
count = 0;
for i = 1:length(copynumber)
    namevObj1 = strcat('v1_',num2str(count));
    vObj1 = addvariant(m1,namevObj1);
    addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});

    vObj = [vObj1];

    simdata = sbiosimulate(m1,csObj,vObj);
    [time, stateData] = selectbyname(simdata, 'Cell.POI');
    POI_openloop(i) = stateData(end);

    count = count + 1;
end
leakage=0.25;
transcription = @(x) (183^4./(183^4+x.^4)+leakage)./(1+leakage);

figure(47)
yyaxis left
loglog(copynumber,POI_openloop,'Color','k','LineWidth',3)
hold on

loglog(copynumber,POI_Eq,'b--','LineWidth',3)  
loglog(copynumber,POI_NF,'r-','LineWidth',3)
loglog(copynumber,POI_IFF,'g-','LineWidth',3)  
xt = get(gca, 'XTick');

k = 0.5716;
theta = 120;
y2 = 1-gamcdf(copynumber,k,theta);
y2 = y2./(y2(2))*100; % normalized so that look at cells have p>=1

% plot distribution of the openloop expression level
yyaxis right
plot(copynumber,y2,'m--','LineWidth',3)
ylim([1e-2 100])
ax = gca;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 20;
ax.TickLength = [0.025, 0.025];
box off
pbaspect([1 1 1])
hold off

save('../SimulationData/Fig4d_Predicted_expression.mat','copynumber','POI_openloop','POI_Eq','POI_NF','POI_IFF')
