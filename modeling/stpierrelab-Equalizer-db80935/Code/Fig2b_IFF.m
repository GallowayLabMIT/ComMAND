% This script generates Figure 2 B. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

sbioloadproject ../Models/IFF_model m1
csObj = m1.addconfigset('newStopTimeConfigSet');

csObj.StopTime = 1e7;

miRNAdissociation = [0.3];
RISC = 1.7e5;

copynumber = 0:1:1000;
POI_IFF = zeros(length(copynumber),length(miRNAdissociation));
RISC_IFF = zeros(length(copynumber),length(miRNAdissociation));
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

for i = 1:length(miRNAdissociation)
    figure(i)
    yyaxis left
    loglog(copynumber,POI_openloop(:),'Color','k','LineWidth',3)
    hold on
    loglog(copynumber,POI_IFF(:,i),'g','LineWidth',3)
    xt = get(gca, 'XTick');
    % ylim([0 1e5])
    xlim([1 1e3])
    yyaxis right
    loglog(copynumber,RISC_IFF(:,i),'b','LineWidth',3)

    ylim([-1e4 18e4])
    ax = gca;
    ax.XAxis.Color = 'k';
    ax.XAxis.LineWidth = 2;
    ax.XAxis.FontSize = 20;
    ax.TickLength = [0.03, 0.03];
    box off
    pbaspect([1 .8 1])
    hold off
end
save('../SimulationData/Fig2b_IFF.mat','copynumber','POI_openloop','RISC_IFF')
