% This script generates Figure 4 G. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

sbioloadproject ../Models/NF_model m1

csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1e7;

extracell_inducer = 10;
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

extracell_inducer = [1];
dox_influx = 0.156.*extracell_inducer;

miRNA = [0.303];
RISC = 1.7e5;

copynumber = 0:1:1000;
POI_Eq = zeros(length(copynumber),length(miRNA),length(extracell_inducer));
RISC_Eq = zeros(length(copynumber),length(miRNA),length(extracell_inducer));
TetR_Eq = zeros(length(copynumber),length(miRNA),length(extracell_inducer));
count = 0;

for i = 1:length(copynumber)
    for j= 1:length(miRNA)
        for k = 1:length(extracell_inducer)
            sbioloadproject ../Models/Equalizer_model m1
            csObj = m1.addconfigset('newStopTimeConfigSet');
            csObj.StopTime = 5e7;
            
            namevObj1 = strcat('v1_',num2str(count));
            vObj1 = addvariant(m1,namevObj1);
            addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});

            vObj = [vObj1];

            set(m1.Reaction(11).KineticLaw.Parameters,'Value',dox_influx(k));
            set(m1.Reaction(13).KineticLaw.Parameters,'Value',miRNA(j));
            set(m1.species(10),'InitialAmount',dox_influx(k)/3.33e-4);
            set(m1.species(13),'InitialAmount',RISC);

            simdata = sbiosimulate(m1,csObj,vObj);
            [~, stateData] = selectbyname(simdata, 'Cell.POI');
            POI_Eq(i,j,k) = stateData(end);
            [~, stateData] = selectbyname(simdata, 'Cell.RISC');
            RISC_Eq(i,j,k) = stateData(end);
            [~, stateData] = selectbyname(simdata, 'Cell.TetR');
            TetR_Eq(i,j,k) = stateData(end);
        end
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
transcription = @(x,leakage) (183^4./(183^4+x.^4)+leakage)./(1+leakage);

figure(9)
Eq_rate = transcription(TetR_Eq,0.25)*100;
loglog(copynumber,Eq_rate,'b-','LineWidth',3)
hold on
NF_rate = transcription(TetR,0.25)*100;
loglog(copynumber,NF_rate,'r-','LineWidth',3)
IFF_rate = 100;
yline(IFF_rate,'g-','LineWidth',3)
Openloop_rate = 100;
yline(Openloop_rate,'k--','LineWidth',3)
xlim([1 1e3])
ylim([0.5 100])
set(gca, 'YScale', 'log')
ax = gca;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 20;
ax.TickLength = [0.03, 0.03];
box off
pbaspect([1 1 1])
hold off
save('../SimulationData/Fig4g_Transcription.mat','copynumber','Eq_rate','IFF_rate','NF_rate','Openloop_rate')
