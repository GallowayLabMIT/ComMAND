% This script generates Supplementary Figure 4 E. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

miRNA = 0.3*[0.2,0.6,1,1.4,1.8];
maxcopy=1000;
copynumber_plot = 1:5:maxcopy;

best_inducer =  1;%in the unit of ng/m
legend_list = {'-80% microRNA affinity', '-40% microRNA affinity', 'approximated microRNA affinity', '+40% microRNA affinity',  '+80% microRNA affinity'};
count = 1;
POI = zeros(length(miRNA),length(copynumber_plot));
min_index_Eq = zeros(length(miRNA),1);
%Numerical differentiation to get local log sensitivity
dPOIdCN = zeros(size(POI));
for m = 1:length(miRNA)
    miRNAdissociation = miRNA(m);

    dox_influx = 0.156.*best_inducer;



    for i = 1:length(copynumber_plot)

        sbioloadproject ../Models/Equalizer_model m1
        csObj = m1.addconfigset('newStopTimeConfigSet');
        csObj.StopTime = 1e7;

        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber_plot(i)});

        vObj = [vObj1];

        set(m1.Reaction(11).KineticLaw.Parameters,'Value',dox_influx);
        set(m1.Reaction(13).KineticLaw.Parameters,'Value',miRNAdissociation);
        set(m1.species(10),'InitialAmount',dox_influx/3.33e-4);

        simdata = sbiosimulate(m1,csObj,vObj);
        [~, stateData] = selectbyname(simdata, 'Cell.POI');
        POI(m,i) = stateData(end);
        count = count + 1;

    end     

    for j = 1:length(copynumber_plot)
       if j == 1
       dPOIdCN(m,j) = (-POI(m,j+2)+4*POI(m,j+1)-3*POI(m,j))/2*(j/POI(m,j));
       elseif j == length(copynumber_plot)
       dPOIdCN(m,j) = (3*POI(m,j)-4*POI(m,j-1)+POI(m,j-2))/2*(j/POI(m,j));
       else 
       dPOIdCN(m,j) = (POI(m,j+1)-POI(m,j-1))/2*(j/POI(m,j));
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
save('../SimulationData/SupplementaryFig4_Eq_microRNAaffinityScan.mat','copynumber_plot','miRNA','predicted_compensation')

