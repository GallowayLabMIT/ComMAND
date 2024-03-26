% This script generates Supplementary Figure 4 A. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

count = 0;

extracellular_inducer = [0 1 5 10 50];  % in the unit of ng/m
legend_list = {'dox = 0','dox = 1','dox = 5','dox = 10','dox = 50','dox = 100','dox = 500','dox = 1000'};
copynumber_plot = 1:5:1000;
    
POI = zeros(length(copynumber_plot),length(extracellular_inducer));
dPOIdCN = zeros(size(POI));
miRNAdissociation = 0.3;
for m = 1:length(extracellular_inducer)

    dox_influx = 0.156.*extracellular_inducer(m);
   
    for i = 1:length(copynumber_plot)

        sbioloadproject ../Models/Equalizer_model m1
        csObj = m1.addconfigset('newStopTimeConfigSet');
      
        csObj.StopTime = 1e7; 
        
        namevObj1 = strcat('v1_',num2str(count));
        vObj1 = addvariant(m1,namevObj1);
        addcontent(vObj1,{'species','gene','InitialAmount',copynumber_plot(i)});

        set(m1.Reaction(11).KineticLaw.Parameters,'Value',dox_influx);
        set(m1.species(10),'InitialAmount',dox_influx/3.33e-4);
        set(m1.Reaction(13).KineticLaw.Parameters,'Value',miRNAdissociation);
        
        simdata = sbiosimulate(m1,csObj,vObj1);
        [time, stateData] = selectbyname(simdata, 'Cell.POI');
        POI(i,m) = stateData(end);
        count = count + 1;
    end
    % Numerical differentiation to get local log sensitivity

    for j = 1:length(copynumber_plot)
       if j == 1 
       dPOIdCN(j,m) = (-POI(j+2,m)+4*POI(j+1,m)-3*POI(j,m))/2*(j/POI(j,m));
       elseif j == length(copynumber_plot)
       dPOIdCN(j,m) = (3*POI(j,m)-4*POI(j-1,m)+POI(j-2,m))/2*(j/POI(j,m));
       else 
       dPOIdCN(j,m) = (POI(j+1,m)-POI(j-1,m))/2*(j/POI(j,m));
%        dPOIdCN(j,m) = (-POI(j+2,m)+8*POI(j+1,m)-8*POI(j-1,m)+POI(j-2,m))/12*(j/POI(j,m));
       end
    end
end

figure(1)
for i = 1:length(extracellular_inducer)
    plot(copynumber_plot(2:end),1./dPOIdCN(2:end,i),'LineWidth',3)
    hold on
end

xt = get(gca, 'XTick');
ax = gca;
ax.XAxis.Color = 'k';
ax.XAxis.LineWidth = 2;
ax.XAxis.FontSize = 18;
xlabel('DNA copy number')
ylabel('Predicted gene dosage compensation')
legend(legend_list,'Location','northeast','FontSize', 16); 
yticks([1 10 20 30 40 50 60 70 80])
xticks([1 100 200 300 400 500])
ylim([0 80])
xlim([1 500])
box off

predicated_compensation = 1./dPOIdCN;
save('../SimulationData/SupplementaryFig4_Eq_logsensitivity_DoxScan.mat','predicated_compensation','copynumber_plot','extracellular_inducer')
