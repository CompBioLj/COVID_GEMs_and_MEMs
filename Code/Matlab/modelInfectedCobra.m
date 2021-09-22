function Cobra_infected_model=modelInfectedCobra(model_Healthy,cobra_vbof)
Cobra_infected_model=model_Healthy
%cytosole
met='sarscov2c[c]';
metList=table2array(cobra_vbof(:,1));
stoicList=table2array(cobra_vbof(:,2));
Cobra_infected_model.metNames(ismember(Cobra_infected_model.mets,met))={'VBOF'};
Cobra_infected_model=addMetabolite(Cobra_infected_model,met);
Cobra_infected_model=addReaction(Cobra_infected_model,'Cov2VBOF','metaboliteList',metList,'stoichCoeffList',stoicList,'reversible',false);
%extracellular
met='sarscov2s[s]';
Cobra_infected_model=addMetabolite(Cobra_infected_model,met);
Cobra_infected_model.metNames(ismember(Cobra_infected_model.mets,met))={'VBOF'};
%Add transport reaction from cytosol to extracellular
Cobra_infected_model=addReaction(Cobra_infected_model,'VBOFt','reactionFormula','sarscov2c[c] -> sarscov2s[s]');
%Add transport reaction from extracellular to boundary
Cobra_infected_model=addExchangeRxn(Cobra_infected_model,'sarscov2s[s]');
Cobra_infected_model=changeObjective(Cobra_infected_model,{'Cov2VBOF'},1);
end