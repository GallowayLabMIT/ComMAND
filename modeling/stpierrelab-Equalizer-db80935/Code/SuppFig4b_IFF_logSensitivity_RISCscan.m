% This script generates Supplementary Figure 4 B. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

RISC_list = 1.7e+05*[0.2,0.6,1,1.4,1.8];
count=0;
maxcopy=1000;
copynumber_plot = 1:5:maxcopy;

miRNAdissociation = 0.3;

legend_list = {'-80% [RISC]', '-40% [RISC]', 'approximated [RISC]', '+40% [RISC]',  '+80% [RISC]'};
%Numerical differentiation to get local log sensitivity
dPOIdCN = zeros(length(RISC_list),length(copynumber_plot));
for m = 1:length(RISC_list)
    RISC = RISC_list(m);

    POI = zeros(length(copynumber_plot),1);
    TetR_Eq = zeros(length(copynumber_plot),1);

    for i = 1:length(copynumber_plot)

        sbioloadproject ../Models/IFF_model m1
        csObj = m1.addconfigset('newStopTimeConfigSet');
        csObj.StopTime = 1e8;

        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber_plot(i)});

        vObj = [vObj1];

        set(m1.Reaction(10).KineticLaw.Parameters,'Value',miRNAdissociation);
        set(m1.species(10),'InitialAmount',RISC);

        simdata = sbiosimulate(m1,csObj,vObj);
        [time, stateData] = selectbyname(simdata, 'Cell.POI');
        POI(i) = stateData(end);
        IFF_POI = stateData;
        count = count + 1;

    end


    for j = 1:length(copynumber_plot)
       if j == 1
       dPOIdCN(m,j) = (-POI(j+2)+4*POI(j+1)-3*POI(j))/2*(j/POI(j));
       elseif j == length(copynumber_plot)
       dPOIdCN(m,j) = (3*POI(j)-4*POI(j-1)+POI(j-2))/2*(j/POI(j));
       else 
       dPOIdCN(m,j) = (POI(j+1)-POI(j-1))/2*(j/POI(j));
       end
    end

    figure(1)
    plot(copynumber_plot(2:end),1./dPOIdCN(m,2:end),'LineWidth',3)
    hold on 

end

xt = get(gca, 'XTick');
ax = gca;
ax.YAxis.Color = 'k';
ax.YAxis.LineWidth = 2;
ax.YAxis.FontSize = 18;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 18;
xlabel('DNA copy number')
ylabel('Predicted gene dosage compensation')
legend(legend_list,'Location','northeast','FontSize', 20); 
legend box off 
yticks([1 10 20 30 40 50 60 70 80])
xticks([1 100 200 300 400 500])
ylim([0 80])
xlim([1 500])

box off
hold off
predicted_compensation = 1./dPOIdCN;
save('../SimulationData/SupplementaryFig4_IFF_RISCscan.mat','copynumber_plot','RISC_list','predicted_compensation')


