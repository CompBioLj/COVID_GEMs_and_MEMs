function[Gimme_model_H, Gimme_model_I]=generate_GIMME(model_Healthy,NHBE_Healthy_mean,essentialTasks_H,Cobra_infected_model,NHBE_Infected_mean,essentialTasks_I)
%% Generation of healthy model
load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
%load global model
exprData.gene=model_Healthy.genes;
exprData.value=NHBE_Healthy_mean;
rxn_expression= mapExpressionToReactions(model_Healthy,exprData);
th_lb=prctile(rxn_expression,75);
for k=1:length(rxn_expression)
    if isnan(rxn_expression(k))==1
        rxn_expression(k)=-1;
    end
end
options.solver='GIMME';
options.expressionRxns=rxn_expression;
options.threshold=th_lb;
Gimme_NHBE_H=createTissueSpecificModel(model_Healthy, options);

Gimme_NHBE_H=ravenCobraWrapper(Gimme_NHBE_H)%put in Raven format
%checkTasks(Gimme_NHBE_H, [], true, false, false, essentialTasks_H);
refModel=ihuman;
%Remove exchange reactions and reactions already included in the INIT
%model
refModelNoExc = removeReactions(refModel,union(Gimme_NHBE_H.rxns,getExchangeRxns(refModel)),true,true);
paramsFT = [];  
[outModel,addedRxnMat] = fitTasks(Gimme_NHBE_H,refModelNoExc,[],true,[],essentialTasks_H,paramsFT);
addedRxnsForTasks = refModelNoExc.rxns(any(addedRxnMat,2));

% The model can now perform all the tasks defined in the task list.
Gimme_NHBE_H = outModel;
Gimme_model_H=Gimme_NHBE_H;

%% Generation of Infected model
exprData.gene=Cobra_infected_model.genes;
exprData.value=NHBE_Infected_mean;
rxn_expression_I= mapExpressionToReactions(Cobra_infected_model,exprData);
th_lb=prctile(rxn_expression_I,75);
for k=1:length(rxn_expression_I)
    if isnan(rxn_expression_I(k))==1
        rxn_expression_I(k)=-1;
    end
end
options.solver='GIMME';
options.expressionRxns=rxn_expression_I;
options.threshold=th_lb;
Gimme_NHBE_I=createTissueSpecificModel(Cobra_infected_model, options);
%Gimme_NHBE_I.id='Gimme_NHBE_I'

Gimme_NHBE_I=ravenCobraWrapper(Gimme_NHBE_I)%put Raven structure

checkTasks(Gimme_NHBE_I, [], true, false, false, essentialTasks_I);


refModel=Cobra_infected_model;
refModel=ravenCobraWrapper(Cobra_infected_model);
refModel=changeObjective(refModel,{'Cov2VBOF'},1);
essentialTasks_imat=essentialTasks_I;
initModel=Gimme_NHBE_I
%initModel.rxns{ismember(initModel.rxns(:,1),'EX_sarscov2s[s]')}='EX_sarscov2s'
%Remove exchange reactions and reactions already included in the INIT
%model
refModelNoExc = removeReactions(refModel,union(initModel.rxns,getExchangeRxns(refModel)),true,true);
paramsFT = [];  
[outModel,addedRxnMat] = fitTasks(initModel,refModelNoExc,[],true,[],essentialTasks_imat,paramsFT);

addedRxnsForTasks = refModelNoExc.rxns(any(addedRxnMat,2));
% The model can now perform all the tasks defined in the task list.
Gimme_NHBE_I = outModel;
Gimme_model_I=Gimme_NHBE_I;

% %infected
modelI=Gimme_NHBE_I

modelI=simplifyModel(modelI);
modelI.b=repelem(0,length(modelI.b));
modelI=changeObjective(modelI,'Cov2VBOF',1)
checkObjective(modelI)

% %healthy
modelUI=Gimme_NHBE_H
checkObjective(modelUI)
modelUI=simplifyModel(modelUI);

% %% convert to Cobra structure 
Gimme_model_H=ravenCobraWrapper(modelUI)
Gimme_model_I=ravenCobraWrapper(modelI)

%% Save model
Gimme_model_H.modelID='Gimme_NHBE_H'
Gimme_model_I.modelID='Gimme_NHBE_I'
end

