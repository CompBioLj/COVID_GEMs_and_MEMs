function[IMAT_model_H, iMAT_model_I]=Generation_iMAT_alt(model_Healthy,NHBE_Healthy_mean,essentialTasks_H,Cobra_infected_model,NHBE_Infected_mean,essentialTasks_I)
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

options.solver='iMAT';
options.expressionRxns=rxn_expression;
options.threshold_lb=th_lb;
options.threshold_ub=th_lb;
iMAT_NHBE_H=createTissueSpecificModel(model_Healthy, options);

iMAT_NHBE_H=ravenCobraWrapper(iMAT_NHBE_H)%put in Raven format

checkTasks(iMAT_NHBE_H, [], true, false, false, essentialTasks_H);
refModel=ihuman;
%Remove exchange reactions and reactions already included in the iMAT
%model
refModelNoExc = removeReactions(refModel,union(iMAT_NHBE_H.rxns,getExchangeRxns(refModel)),true,true);
paramsFT = [];  
[outModel,addedRxnMat] = fitTasks(iMAT_NHBE_H,refModelNoExc,[],true,[],essentialTasks_H,paramsFT);

addedRxnsForTasks = refModelNoExc.rxns(any(addedRxnMat,2));

% The model can now perform all the tasks defined in the task list.
iMAT_NHBE_H = outModel;

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
options.solver='iMAT';
options.expressionRxns=rxn_expression_I;
options.threshold_lb=th_lb;
options.threshold_ub=th_lb;
iMAT_NHBE_I=createTissueSpecificModel(Cobra_infected_model, options);

iMAT_NHBE_I=ravenCobraWrapper(iMAT_NHBE_I)%put Raven structure

checkTasks(iMAT_NHBE_I, [], true, false, false, essentialTasks_I);

refModel=Cobra_infected_model;
refModel=ravenCobraWrapper(Cobra_infected_model);
refModel=changeObjective(refModel,{'Cov2VBOF'},1);
essentialTasks_imat=essentialTasks_I;
initModel=iMAT_NHBE_I
%Remove exchange reactions and reactions already included in the iMAT
%model
refModelNoExc = removeReactions(refModel,union(initModel.rxns,getExchangeRxns(refModel)),true,true);
paramsFT = [];  
[outModel,addedRxnMat] = fitTasks(initModel,refModelNoExc,[],true,[],essentialTasks_imat,paramsFT);
addedRxnsForTasks = refModelNoExc.rxns(any(addedRxnMat,2));
% The model can now perform all the tasks defined in the task list.
iMAT_NHBE_I = outModel;

%%  Constrain media by HAM media
%infected
modelI=iMAT_NHBE_I
checkTasks(modelI, [], true, false, false, essentialTasks_I);
modelI=simplifyModel(modelI);
modelI.b=repelem(0,length(modelI.b))';
modelI=changeObjective(modelI,'Cov2VBOF',1)
checkObjective(modelI)
%healthy
modelUI=iMAT_NHBE_H
modelUI=simplifyModel(modelUI);

%% convert to Cobra structure 
IMAT_model_H=ravenCobraWrapper(modelUI)
iMAT_model_I=ravenCobraWrapper(modelI)

%% Save model
IMAT_model_H.modelID='IMAT_model_H'
iMAT_model_I.modelID='iMAT_model_I'
end