function [modelUI,modelI]=INIT_generation_good(A549_Healthy_mean,A549_Infected_mean)

%Load model and convert model
load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);

%%healthy
data_healthy.genes=ihuman.genes;
data_healthy.tissues={'UI'}; 
data_healthy.levels=A549_Healthy_mean; 
data_healthy.threshold=prctile(A549_Healthy_mean ,75) 
essentialTasksH=parseTaskList('metabolicTasks_Essential.xlsx');
tissueH = 'UI';
params.TimeLimit=5000;

%Run the getINITModel2 command. The execution takes around 2 hours.
tINIT_model_H_A549 = getINITModel2(ihuman, tissueH, [], [], data_healthy, [], true, [], true, true, essentialTasksH, params, []);


%% infected

model_I=model_infected(ihuman);
model_I=changeObjective(model_I,{'Cov2VBOF'},1);
data_infected.genes=model_I.genes;
data_infected.tissues={'I'};
data_infected.levels=A549_Infected_mean;
data_infected.threshold=prctile(A549_Infected_mean ,75) 
essentialTasksI=parseTaskList('metabolicTasks_Essential_infected.xlsx');
tissueI='I';
params.TimeLimit=5000;  % additional optimization parameters for the INIT algorithm. 5000 s is optimal. 
tINIT_model_I_A549 = getINITModel2(model_I, tissueI, [], [], data_infected, [], true, [], true, true, [], params, []);


%%
modelUI=tINIT_model_H_A549
modelI=tINIT_model_I_A549
checkTasks(modelUI, [], true, false, false, essentialTasksUI);
modelI=addBoundaryMets(modelI);
checkTasks(modelI, [], true, false, false, essentialTasksI);

%Simplify the model to add the exchange reaction and make it worksing. 
modelUI=simplifyModel(modelUI);
modelI=simplifyModel(modelI);
end

