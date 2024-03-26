% This script generates Figure 4 E. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

extracell_inducer = [1];
dox_influx = 0.156.*extracell_inducer;

copynumber = 0:1:1000;
POI_Eq = zeros(length(copynumber),length(extracell_inducer));
RISC_Eq = zeros(length(copynumber),length(extracell_inducer));
TetR_Eq = zeros(length(copynumber),length(extracell_inducer));
mRNA_Eq = zeros(length(copynumber),length(extracell_inducer));
mRNAmicroRNA_Eq = zeros(length(copynumber),length(extracell_inducer));
microRNARISC_Eq = zeros(length(copynumber),length(extracell_inducer));

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
        [~, stateData] = selectbyname(simdata, 'Cell.TetR');
        TetR_Eq(i,k) = stateData(end);
        [~, stateData] = selectbyname(simdata, 'Cell.mRNA');
        mRNA_Eq(i,k) = stateData(end);
        [~, stateData] = selectbyname(simdata, 'Cell.mRNAmicroRNA');
        mRNAmicroRNA_Eq(i,k) = stateData(end);
        [~, stateData] = selectbyname(simdata, 'Cell.microRNA-RISC');
        microRNARISC_Eq(i,k) = stateData(end);
    end
end

sbioloadproject ../Models/NF_model m1

csObj = m1.addconfigset('newStopTimeConfigSet');
csObj.StopTime = 1e7;

extracell_inducer = 10;
% extracell_inducer = 0.5;
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

miRNAdissociation = [0.303];
RISC = 1.7e5;

POI_IFF = zeros(length(copynumber),length(miRNAdissociation));
RISC_IFF = zeros(length(copynumber),length(miRNAdissociation));
mRNA_IFF = zeros(length(copynumber),length(miRNAdissociation));
mRNAmicroRNA_IFF = zeros(length(copynumber),length(miRNAdissociation));
microRNARISC_IFF = zeros(length(copynumber),length(miRNAdissociation));
count = 0;

for i = 1:length(copynumber)
    for k = 1:length(miRNAdissociation)

        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});
        set(m1.Reaction(10).KineticLaw.Parameters,'Value',miRNAdissociation(k));
        set(m1.species(10),'InitialAmount',RISC);
        simdata = sbiosimulate(m1,csObj,vObj1);
        [~, stateData] = selectbyname(simdata, 'Cell.POI');
        POI_IFF(i,k) = stateData(end);
        [~, stateData] = selectbyname(simdata, 'Cell.RISC');
        RISC_IFF(i,k) = stateData(end);
        [~, stateData] = selectbyname(simdata, 'Cell.mRNA');
        mRNA_IFF(i,j,k) = stateData(end);
        [~, stateData] = selectbyname(simdata, 'Cell.mRNAmicroRNA');
        mRNAmicroRNA_IFF(i,j,k) = stateData(end);
        [~, stateData] = selectbyname(simdata, 'Cell.microRNA-RISC');
        microRNARISC_IFF(i,j,k) = stateData(end);
        count = count + 1;
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

transcription = @(x,leakage) (183^4./(183^4+x.^4)+leakage)./(1+leakage);

figure(47)
Eq_rate = transcription(TetR_Eq,0.25).*(mRNA_Eq.*3.33e-4./(mRNAmicroRNA_Eq.*0.002)./(3.33e-4/2.88e-4)).*100;
loglog(copynumber,Eq_rate,'b-','LineWidth',3)
hold on
NF_rate = transcription(TetR_NF,0.25).*100;
loglog(copynumber,transcription(TetR_NF,0.25).*100,'r-','LineWidth',3)
IFF_rate = (mRNA_IFF.*3.33e-4./(mRNAmicroRNA_IFF.*0.002)./(3.33e-4/2.88e-4)).*100;
loglog(copynumber,IFF_rate,'g-','LineWidth',3)
Openloop_rate = 100;
yline(100,'k--','LineWidth',3)
ylim([.5 100])
xlim([1 1000])
ax = gca;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 20;
ax.TickLength = [0.025, 0.025];
box off
pbaspect([1 1 1])
hold off

save('../SimulationData/Fig4e_ExpressionRate.mat','copynumber','Eq_rate','IFF_rate','NF_rate','Openloop_rate')