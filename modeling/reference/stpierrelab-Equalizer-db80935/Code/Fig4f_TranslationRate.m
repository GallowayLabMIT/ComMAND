% This script generates Figure 4 F. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

extracell_inducer = [1];
dox_influx = 0.156.*extracell_inducer;

miRNA = [0.303];
RISC = 1.7e5;

copynumber = 0:1:1000;
POI_Eq = zeros(length(copynumber),length(miRNA),length(extracell_inducer));
RISC_Eq = zeros(length(copynumber),length(miRNA),length(extracell_inducer));
TetR_Eq = zeros(length(copynumber),length(miRNA),length(extracell_inducer));
mRNA_Eq = zeros(length(copynumber),length(miRNA),length(extracell_inducer));
mRNAmicroRNA_Eq = zeros(length(copynumber),length(miRNA),length(extracell_inducer));
microRNARISC_Eq = zeros(length(copynumber),length(miRNA),length(extracell_inducer));

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
            [~, stateData] = selectbyname(simdata, 'Cell.mRNA');
            mRNA_Eq(i,j,k) = stateData(end);
            [~, stateData] = selectbyname(simdata, 'Cell.mRNAmicroRNA');
            mRNAmicroRNA_Eq(i,j,k) = stateData(end);
            [~, stateData] = selectbyname(simdata, 'Cell.microRNA-RISC');
            microRNARISC_Eq(i,j,k) = stateData(end);
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

copynumber = 0:1:1000;
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
    end
end


figure(21)
Eq_rate = (mRNA_Eq(:).*3.33e-4./(mRNAmicroRNA_Eq(:).*0.002)./(3.33e-4/2.88e-4))*100;
loglog(copynumber,Eq_rate,'b-','LineWidth',3)
hold on
IFF_rate = (mRNA_IFF(:).*3.33e-4./(mRNAmicroRNA_IFF(:).*0.002)./(3.33e-4/2.88e-4))*100;
loglog(copynumber,IFF_rate,'g-','LineWidth',3)
NF_rate = 100;
yline(NF_rate,'r-','LineWidth',3)
Openloop_rate = 100;
yline(Openloop_rate,'k--','LineWidth',3)
set(gca,'Yscale','log')
xlim([1 1e3])
ylim([.5 100])
ax = gca;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 20;
ax.TickLength = [0.03, 0.03];
box off
pbaspect([1 1 1])
hold off
save('../SimulationData/Fig4f_TranslationRate.mat','copynumber','Eq_rate','IFF_rate','NF_rate','Openloop_rate')