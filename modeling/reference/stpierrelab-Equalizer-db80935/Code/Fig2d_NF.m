% This script generates Figure 2 D. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

sbioloadproject ../Models/NF_model m1

csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1e7;

extracell_inducer = 1;
dox_influx = 0.156.*extracell_inducer;

copynumber = 0:1:1000;
POI = zeros(length(copynumber),length(extracell_inducer));
TetR = zeros(length(copynumber),length(extracell_inducer));
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
        POI(i,j) = stateData(end);
        [~, stateData] = selectbyname(simdata, 'Cell.TetR');
        TetR(i,j) = stateData(end);
        count = count + 1;
    end
end

dPOIdCN = zeros(size(POI));
for i = 1:length(extracell_inducer)
    for j = 1:length(copynumber)
        if j == 1
            dPOIdCN(j,i) = (-POI(j+2,i)+4*POI(j+1,i)-3*POI(j,i))/2*(j/POI(j,i));
        elseif j == length(copynumber)
            dPOIdCN(j,i) = (3*POI(j,i)-4*POI(j-1,i)+POI(j-2,i))/2*(j/POI(j,i));
        else 
            dPOIdCN(j,i) = (POI(j+1,i)-POI(j-1,i))/2*(j/POI(j,i));
        end
    end
end

sbioloadproject ../Models/NF_model_ideal m1 m1

csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1e7;

OutDir = 'C:\Research\Membrane Localization Improvement Project\Model\Simbiology Models\MatFiles\';

POI_ideal = zeros(length(copynumber),1);
TetR_ideal = zeros(length(copynumber),1);
dox_influx_ideal = 0.156.*extracell_inducer;
count = 0;
for i = 1:length(copynumber)
    namevObj1 = strcat('v1_',num2str(count));
    vObj1 = addvariant(m1,namevObj1);
    addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});

    namevObj2 = strcat('v2_',num2str(count));
    vObj2 = addvariant(m1,namevObj2);
    addcontent(vObj2,{'species','Inducer','InitialAmount',dox_influx_ideal/3.33e-4});

    vObj = [vObj1 vObj2];

    set(m1.Reaction(6).KineticLaw.Parameters,'Value',dox_influx_ideal);

    simdata = sbiosimulate(m1,csObj,vObj);
    [time, stateData] = selectbyname(simdata, 'Cell.POI');
    POI_ideal(i) = stateData(end);
    [~, stateData] = selectbyname(simdata, 'Cell.TetR');
    TetR_ideal(i) = stateData(end);
    count = count + 1;
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
transcription = @(x,leakage) (183^4./(183^4+x.^4)+leakage)./(1+leakage);

figure(1)
yyaxis left
loglog(copynumber,POI_openloop(:),'k-','LineWidth',3)
hold on
loglog(copynumber,POI_ideal(:),'b-','LineWidth',3)
loglog(copynumber,POI(:),'r-','LineWidth',3)
xt = get(gca, 'XTick');
ylim([1e2 1e6])
xlim([1 1e3])
yyaxis right
transcription_ideal = transcription(TetR_ideal,0);
transcription_leaky = transcription(TetR,0.25);
loglog(copynumber,transcription_ideal,'b--','LineWidth',3)
loglog(copynumber,transcription_leaky,'r--','LineWidth',3)
ylim([0.004 1])
set(gca, 'YScale', 'log')
ax = gca;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 20;
ax.TickLength = [0.03, 0.03];
box off
pbaspect([1 .8 1])
hold off
save('../SimulationData/Fig2d_NF.mat','copynumber','POI_ideal','POI','transcription_ideal','transcription_leaky')
