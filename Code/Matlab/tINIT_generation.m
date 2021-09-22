function [modelUI,modelI]=tINIT_generation(A549_Healthy_mean,A549_Infected_mean)

%Load model and convert model
load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);

%% healthy
data_healthy.genes=ihuman.genes;
data_healthy.tissues={'UI'};
data_healthy.levels=A549_Healthy_mean; 
data_healthy.threshold=prctile(A549_Healthy_mean ,75) 
tissueUI = 'UI';
taskStructureH = parseTaskList('metabolicTasks_Essential.xlsx');
params.TimeLimit=5000;  
tINIT_model_H_A549 = getINITModel2(ihuman, tissueUI, [], [], data_healthy, [], true, [], true, true, taskStructureH, params, []);


%% infected
refModelI=model_infected(ihuman);
refModelI=changeObjective(refModelI,{'Cov2VBOF'},1);
data_infected.genes=refModelI.genes;
data_infected.tissues={'I'};
data_infected.levels=A549_Infected_mean;
data_infected.threshold=prctile(A549_Infected_mean ,75) 
tissueI='I';
taskStructureI = parseTaskList('metabolicTasks_Essential_infected.xlsx');
tINIT_model_I_A549 = getINITModel2(refModelI, tissueI, [], [], data_infected, [], true, [], true, true, taskStructureI, params, []);


modelUI=tINIT_model_H_A549
modelI=tINIT_model_I_A549
checkTasks(modelUI, [], true, false, false, essentialTasksUI);

modelI=removeMets(modelI,'sarscov2s');
modelI=model_infected(modelI);
modelI=addBoundaryMets(modelI);
checkTasks(modelI, [], true, false, false, essentialTasksI);

%Simplify the model to add the exchange reaction and make it worksing. 
modelUI=simplifyModel(modelUI);
modelI=simplifyModel(modelI);
end

