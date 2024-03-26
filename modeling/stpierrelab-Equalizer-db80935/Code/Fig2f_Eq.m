% This script generates Figure 2 F. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

extracell_inducer = [1];
dox_influx = 0.156.*extracell_inducer;

miRNAdissociation = [0.3];
RISC = 1.7e5;

copynumber = 0:1:1000;
POI_Eq = zeros(length(copynumber),length(miRNAdissociation),length(extracell_inducer));
RISC_Eq = zeros(length(copynumber),length(miRNAdissociation),length(extracell_inducer));
TetR_Eq = zeros(length(copynumber),length(miRNAdissociation),length(extracell_inducer));
count = 0;

for i = 1:length(copynumber)
    for j= 1:length(miRNAdissociation)
        for k = 1:length(extracell_inducer)
            sbioloadproject ../Models/Equalizer_model m1
            csObj = m1.addconfigset('newStopTimeConfigSet');
            csObj.StopTime = 1e7;
            
            namevObj1 = strcat('v1_',num2str(count));
            vObj1 = addvariant(m1,namevObj1);
            addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});

            vObj = [vObj1];

            set(m1.Reaction(11).KineticLaw.Parameters,'Value',dox_influx(k));
            set(m1.Reaction(13).KineticLaw.Parameters,'Value',miRNAdissociation(j));
            set(m1.species(10),'InitialAmount',dox_influx(k)/3.33e-4);
            set(m1.species(13),'InitialAmount',RISC);

            simdata = sbiosimulate(m1,csObj,vObj);
            [~, stateData] = selectbyname(simdata, 'Cell.POI');
            POI_Eq(i,j,k) = stateData(end);
            [~, stateData] = selectbyname(simdata, 'Cell.RISC');
            RISC_Eq(i,j,k) = stateData(end);
            [~, stateData] = selectbyname(simdata, 'Cell.RISC');
            RISC_Eq(i,j,k) = stateData(end);
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
leakage=0.25;
transcription = @(x) (183^4./(183^4+x.^4)+leakage)./(1+leakage);


for i = 1:length(miRNAdissociation)
    figure(i)
    loglog(copynumber,POI_openloop(:),'Color','k','LineWidth',3)
    hold on
    % loglog(copynumber,POI_Eq(:,1,1),'g-','LineWidth',3)
    loglog(copynumber,POI_Eq(:,i,1),'b--','LineWidth',3)  
    xt = get(gca, 'XTick');
    % ylim([0 1e5])
    xlim([1 1e3])
    ax = gca;
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 2;
    ax.XAxis.FontSize = 20;
    ax.TickLength = [0.03, 0.03];
    box off
    pbaspect([1 .8 1])
    hold off
end

save('../SimulationData/Fig2f_Eq.mat','miRNAdissociation','copynumber','POI_openloop','POI_Eq')
