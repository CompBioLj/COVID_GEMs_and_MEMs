function infected_model=model_infected(model)


infected_model=model
vbof=readtable('VBOF.csv');
metList=table2array(vbof(:,1));
stoicList=table2array(vbof(:,2));

%add Sarscov2 biomass objectif function in cytolsol
met='sarscov2c';

infected_model=addMetabolite(infected_model,met);
position=ismember(infected_model.mets,met)
infected_model.metNames(position)={'VBOF'}
infected_model.metComps(position)=4
infected_model=addReaction(infected_model,'Cov2VBOF','metaboliteList',metList,'stoichCoeffList',stoicList,'reversible',false);

%add Sarscov2 in extracellular space and exchange
met2='sarscov2s';
infected_model=addMetabolite(infected_model,met2);
position2=ismember(infected_model.mets,met2)
infected_model.metNames(position2)={'VBOF'};
infected_model.metComps(position2)=1; 
infected_model=addReaction(infected_model,'VBOFt','reactionFormula','sarscov2c -> sarscov2s');
infected_model=addExchangeRxn(infected_model,'sarscov2s');

end
