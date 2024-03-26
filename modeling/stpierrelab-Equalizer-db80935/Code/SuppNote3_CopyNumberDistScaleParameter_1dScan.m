% This script generates data fpr Supplementary Note 3. 
% 
% Jin Yang, harvey1@mit.edu
% St-Pierre Lab (stpierrelab.com) Apr. 2021

%%Generate the population with log normal copy number distribution
leakage = 0.25;
% shape parameter
k=0.5716;
theta = 60:1:200;
MSE = zeros(length(leakage),length(theta));

%Load experimental data
load('../ExperimentData/SuppNote3_NF_expression.mat')
expmean = mean(avg);
% normalize the expression level
expmean = expmean/max(expmean);

for m = 1:length(leakage)
    for n = 1:length(theta)
        celltot = 10000; %Assume a total of 1000 cells
        copynumber_cell = round(gamrnd(k,theta(n),celltot,1)); %copynumber in each cell,std of 1.522 fit from PGK A9 well on 20190419
        copynumber_cell = copynumber_cell(copynumber_cell>0);
        copynumber_cell_sorted = sort(copynumber_cell);
        expression_vector=zeros(celltot,1);
        sbioloadproject ../Models/NF_model m1
        csObj = m1.addconfigset('newStopTimeConfigSet');
        csObj.StopTime = 1e6;
        RNA_deg = 2.88e-4;
        extracell_inducer = [0 1 5 10 50 100 500 1000];%in the unit of ng/m
        dox_influx = 0.156.*extracell_inducer;

        %vector to hold COV at each induction level
        Meanexp = zeros(1, length(extracell_inducer));

        copynumber = unique(copynumber_cell_sorted);
        POI = zeros(length(copynumber),length(extracell_inducer));
        count = 0;
        for i = 1:length(copynumber)
            for j = 1:length(extracell_inducer)
                namevObj1 = strcat('v1_',num2str(count));
                vObj1 = addvariant(m1,namevObj1);
                addcontent(vObj1,{'species','gene','InitialAmount',copynumber(i)});

                set(m1.Reaction(6).KineticLaw.Parameters,'Value',dox_influx(j));
                set(m1.species(5),'InitialAmount',dox_influx(j)/3.33e-4);
                set(m1.Parameters(5),'Value',leakage(m))
                
                simdata = sbiosimulate(m1,csObj,vObj1);
                [~, stateData] = selectbyname(simdata, 'Cell.POI');
                POI(i,j) = stateData(end);
                count = count + 1;
            end
        end

        for i = 1:length(extracell_inducer)
            explvl = POI(:,i)';
            for j = 1:length(copynumber_cell)
                if copynumber_cell(j)==0
                    expression_vector(j)=0;
                else
                    expression_vector(j) = explvl(copynumber==copynumber_cell(j));
                end
            end
            Meanexp(i) = mean(expression_vector);
        end

        % normalize the simulated expression level
        Meanexp = Meanexp/max(Meanexp);

        MSE(m,n) = sumsqr(Meanexp-expmean)/length(extracell_inducer);
    end
end
plot(theta,-MSE)
save('../SimulationData/SupplementaryNote3_theta_1dScan.mat','MSE','leakage','theta')