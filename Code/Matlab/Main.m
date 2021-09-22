clear all
%initialise Cobra toolbox
initCobraToolbox(false)
changeCobraSolver('gurobi','all');
setRavenSolver('gurobi');
%Load model and convert model
load('Human-GEM.mat');
ihuman = addBoundaryMets(ihuman);
essentialTasks_H=parseTaskList('metabolicTasks_Essential.xlsx');
essentialTasks_I=parseTaskList('metabolicTasks_Essential_infected.xlsx');
model_Healthy=ravenCobraWrapper(ihuman)
cobra_vbof=readtable('VBOF_cobra.csv');
Cobra_infected_model=modelInfectedCobra(model_Healthy,cobra_vbof)



Database=input('which data do you want to use :\n- NHBE \n- A549 \n- Calu-3 \n- Lung\n- 293T\n- All\n Write your choice: ','s');
while(ismember(['NHBE','A549','Calu-3','Lung','293T','All'],Database)==0)
    Database=input('Input not valid \n Please chose between this :\n- NHBE \n- A549 \n- Calu-3 \n- 293T\n- Lung\n- 293T\n- All\n  Write your choice:','s')
end
 model_choice=input('which model do you want to use :\n- Gimme \n- iMAT \n- tINIT \n- INIT\n- All\n' ,'s');
        while(ismember(['Gimme','IMAT','tINIT','INIT','All'],model_choice)==0)
            model_choice=input('Input not valid \n Please chose between this :n- Gimme \n- iMAT \n- tINIT \n- INIT\n','s')
        end
switch (Database)
    case 'NHBE'
         load('NHBE_Healthy_mean.mat');
        NHBE_Healthy_mean=cell2mat(NHBE_Healthy_mean);
        NHBE_Healthy_mean_INIT=NHBE_Healthy_mean;
        for k=1:length(NHBE_Healthy_mean_INIT)
            if NHBE_Healthy_mean_INIT(k)==-1
                NHBE_Healthy_mean_INIT(k)=NaN;
            end
        end       
        load('NHBE_Infected_mean.mat');
        NHBE_Infected_mean=cell2mat(NHBE_Infected_mean);
        NHBE_Infected_mean_INIT=NHBE_Infected_mean;
        for k=1:length(NHBE_Infected_mean_INIT)
            if NHBE_Infected_mean_INIT(k)==-1
                NHBE_Infected_mean_INIT(k)=NaN
            end
        end
        
        model_choice=input('which model do you want to use :\n- Gimme \n- iMAT \n- tINIT \n- INIT\n- All\n' ,'s');
        while(ismember(['Gimme','IMAT','tINIT','INIT','All'],model_choice)==0)
            model_choice=input('Input not valid \n Please chose between this :n- Gimme \n- iMAT \n- tINIT \n- INIT\n','s')
        end
        file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
        file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');
        switch (model_choice)
            case 'Gimme'
                [Gimme_model_H_NHBE_alt, Gimme_model_I_NHBE_alt]=Generation_Gimme_alt(model_Healthy,NHBE_Healthy_mean,essentialTasks_H,Cobra_infected_model,NHBE_Infected_mean,essentialTasks_I);
                Gimme_model_H_NHBE_alt.rules(204)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_NHBE_alt.rxns{204})));
                writeCbModel(Gimme_model_H_NHBE_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_NHBE_alt, 'format','mat', 'fileName', file_name_Infected)
                
            case 'IMAT'
                [iMAT_model_H_NHBE_alt, iMAT_model_I_NHBE_alt]=Generation_iMAT_alt(model_Healthy,NHBE_Healthy_mean,essentialTasks_H,Cobra_infected_model,NHBE_Infected_mean,essentialTasks_I);
                iMAT_model_I_NHBE_alt.rules(2872)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),iMAT_model_I_NHBE_alt.rxns{2872})));   
                writeCbModel(iMAT_model_H_NHBE_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(iMAT_model_I_NHBE_alt, 'format','mat', 'fileName', file_name_Infected)
            
            case 'tINIT'
                [tINIT_model_H_NHBE,tINIT_model_I_NHBE]=tINIT_generation(NHBE_Healthy_mean_INIT,NHBE_Infected_mean_INIT)
                tINIT_model_H_NHBE=ravenCobraWrapper(tINIT_model_H_NHBE);
                tINIT_model_I_NHBE=rmfield(tINIT_model_I_NHBE,'rules');
                tINIT_model_I_NHBE=ravenCobraWrapper(tINIT_model_I_NHBE);
                writeCbModel(tINIT_model_H_NHBE, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_NHBE, 'format','mat', 'fileName', file_name_Infected)
                
            case 'INIT'
                [INIT_model_H_NHBE,INIT_model_I_NHBE]=INIT_generation_good(NHBE_Healthy_mean_INIT,NHBE_Infected_mean_INIT)
                INIT_model_H_NHBE=ravenCobraWrapper(INIT_model_H_NHBE);
                INIT_model_I_NHBE=rmfield(INIT_model_I_NHBE,'rules');
                INIT_model_I_NHBE=ravenCobraWrapper(INIT_model_I_NHBE);
                writeCbModel(INIT_model_H_NHBE, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(INIT_model_I_NHBE, 'format','mat', 'fileName', file_name_Infected)
                
            case'All'
                %Gimme
                model_choice='Gimme'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');    
                [Gimme_model_H_NHBE_alt, Gimme_model_I_NHBE_alt]=Generation_Gimme_alt(model_Healthy,NHBE_Healthy_mean,essentialTasks_H,Cobra_infected_model,NHBE_Infected_mean,essentialTasks_I);
                Gimme_model_H_NHBE_alt.rules(204)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_NHBE_alt.rxns{204})));
                writeCbModel(Gimme_model_H_NHBE_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_NHBE_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %IMAT
                model_choice='IMAT'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');    
                [iMAT_model_H_NHBE_alt, iMAT_model_I_NHBE_alt]=Generation_iMAT_alt(model_Healthy,NHBE_Healthy_mean,essentialTasks_H,Cobra_infected_model,NHBE_Infected_mean,essentialTasks_I);
                iMAT_model_I_NHBE_alt.rules(2872)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),iMAT_model_I_NHBE_alt.rxns{2872})));   
                writeCbModel(iMAT_model_H_NHBE_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(iMAT_model_I_NHBE_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %tINIT
                model_choice='tINIT'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');    
                [tINIT_model_H_NHBE,tINIT_model_I_NHBE]=tINIT_generation(NHBE_Healthy_mean_INIT,NHBE_Infected_mean_INIT)
                tINIT_model_H_NHBE=ravenCobraWrapper(tINIT_model_H_NHBE);
                tINIT_model_I_NHBE=rmfield(tINIT_model_I_NHBE,'rules');
                tINIT_model_I_NHBE=ravenCobraWrapper(tINIT_model_I_NHBE);
                writeCbModel(tINIT_model_H_NHBE, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_NHBE, 'format','mat', 'fileName', file_name_Infected)
                
                
                %INIT
                model_choice='INIT'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');    
                [INIT_model_H_NHBE,INIT_model_I_NHBE]=INIT_generation_good(NHBE_Healthy_mean_INIT,NHBE_Infected_mean_INIT)
                INIT_model_H_NHBE=ravenCobraWrapper(INIT_model_H_NHBE);
                INIT_model_I_NHBE=rmfield(INIT_model_I_NHBE,'rules');
                INIT_model_I_NHBE=ravenCobraWrapper(INIT_model_I_NHBE);
                writeCbModel(INIT_model_H_NHBE, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(INIT_model_I_NHBE, 'format','mat', 'fileName', file_name_Infected)     
        end
        
    case 'A549'
        load('A549_Healthy_mean.mat');
        A549_Healthy_mean=cell2mat(A549_Healthy_mean);
        A549_Healthy_mean_INIT=A549_Healthy_mean;
        for p=1:length(A549_Healthy_mean)
            if isnan(A549_Healthy_mean(p))==1
                A549_Healthy_mean(p)=-1
            end
        end
        
        load('A549_Infected_mean.mat');
        A549_Infected_mean=cell2mat(A549_Infected_mean);
        A549_Infected_mean_INIT=A549_Infected_mean;
        for p=1:length(A549_Infected_mean)
            if isnan(A549_Infected_mean(p))==1
                A549_Infected_mean(p)=-1
            end
        end
        
        model_choice=input('which model do you want to use :\n- Gimme \n- iMAT \n- tINIT \n- INIT\n- All\n' ,'s');
        while(ismember(['Gimme','IMAT','tINIT','INIT','All'],model_choice)==0)
            model_choice=input('Input not valid \n Please chose between this :n- Gimme \n- iMAT \n- tINIT \n- INIT\n','s')
        end
        file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
        file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');        
        switch (model_choice)
            case'Gimme'
                [Gimme_model_H_A549_alt, Gimme_model_I_A549_alt]=Generation_Gimme_alt(model_Healthy,A549_Healthy_mean,essentialTasks_H,Cobra_infected_model,A549_Infected_mean,essentialTasks_I);               
                Gimme_model_H_A549_alt.rules(206)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_A549_alt.rxns{206})));
                writeCbModel(Gimme_model_H_A549_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_A549_alt, 'format','mat', 'fileName', file_name_Infected)
                
            case'IMAT'
                [IMAT_model_H_A549_alt, iMAT_model_I_A549_alt]=Generation_iMAT_alt(model_Healthy,A549_Healthy_mean,essentialTasks_H,Cobra_infected_model,A549_Infected_mean,essentialTasks_I);
                iMAT_model_I_A549_alt.rules(2889)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_A549_alt.rxns{2889})));
                writeCbModel(IMAT_model_H_A549_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(iMAT_model_I_A549_alt, 'format','mat', 'fileName', file_name_Infected)
                
            case'tINIT'
                [tINIT_model_H_A549,tINIT_model_I_A549]=tINIT_generation(A549_Healthy_mean_INIT,A549_Infected_mean_INIT)
                tINIT_model_I_A549=rmfield(tINIT_model_I_A549,'rules')
                tINIT_model_I_A549=ravenCobraWrapper(tINIT_model_I_A549)
                writeCbModel(tINIT_model_H_A549, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_A549, 'format','mat', 'fileName', file_name_Infected)
            case 'INIT'
                [INIT_model_H_A549,INIT_model_I_A549]=INIT_generation_good(A549_Healthy_mean_INIT,A549_Infected_mean_INIT)
                INIT_model_H_A549=ravenCobraWrapper(INIT_model_H_A549)
                writeCbModel(INIT_model_H_A549, 'format','mat', 'fileName', 'INIT_model_H_A549_alt.mat')
                INIT_model_I_A549=rmfield(INIT_model_I_A549,'rules')
                INIT_model_I_A549=ravenCobraWrapper(INIT_model_I_A549)
                writeCbModel(INIT_model_I_A549, 'format','mat', 'fileName', 'INIT_model_I_A549_alt.mat')
                
            case'All'
                %Gimme
                model_choice='Gimme'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');     
                [Gimme_model_H_A549_alt, Gimme_model_I_A549_alt]=Generation_Gimme_alt(model_Healthy,A549_Healthy_mean,essentialTasks_H,Cobra_infected_model,A549_Infected_mean,essentialTasks_I);               
                Gimme_model_H_A549_alt.rules(206)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_A549_alt.rxns{206})));
                writeCbModel(Gimme_model_H_A549_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_A549_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %iMAT
                model_choice='IMAT'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');     
                [IMAT_model_H_A549_alt, iMAT_model_I_A549_alt]=Generation_iMAT_alt(model_Healthy,A549_Healthy_mean,essentialTasks_H,Cobra_infected_model,A549_Infected_mean,essentialTasks_I);
                iMAT_model_I_A549_alt.rules(2889)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_A549_alt.rxns{2889})));
                writeCbModel(IMAT_model_H_A549_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(iMAT_model_I_A549_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %tINIT
                model_choice='tINIT'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');    
                [tINIT_model_H_A549,tINIT_model_I_A549]=tINIT_generation(A549_Healthy_mean_INIT,A549_Infected_mean_INIT)
                tINIT_model_I_A549=rmfield(tINIT_model_I_A549,'rules')
                tINIT_model_I_A549=ravenCobraWrapper(tINIT_model_I_A549)
                writeCbModel(tINIT_model_H_A549, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_A549, 'format','mat', 'fileName', file_name_Infected)
                
                %INIT
                model_choice='INIT'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');    
                [INIT_model_H_A549,INIT_model_I_A549]=INIT_generation_good(A549_Healthy_mean_INIT,A549_Infected_mean_INIT)
                INIT_model_H_A549=ravenCobraWrapper(INIT_model_H_A549)
                writeCbModel(INIT_model_H_A549, 'format','mat', 'fileName', 'INIT_model_H_A549_alt.mat')
                INIT_model_I_A549=rmfield(INIT_model_I_A549,'rules')
                INIT_model_I_A549=ravenCobraWrapper(INIT_model_I_A549)
                writeCbModel(INIT_model_I_A549, 'format','mat', 'fileName', 'INIT_model_I_A549_alt.mat')   
        end
        
    case 'CALU'
        
        load('CALU_Healthy_mean.mat');
        CALU_Healthy_mean=cell2mat(CALU_Healthy_mean);
        CALU_Healthy_mean_INIT=CALU_Healthy_mean;
        for p=1:length(CALU_Healthy_mean)
            if isnan(CALU_Healthy_mean(p))==1
                CALU_Healthy_mean(p)=-1
            end
        end
        load('CALU_Infected_mean.mat');
        CALU_Infected_mean=cell2mat(CALU_Infected_mean);
        CALU_Infected_mean_INIT=CALU_Infected_mean;
        for p=1:length(CALU_Infected_mean)
            if isnan(CALU_Infected_mean(p))==1
                CALU_Infected_mean(p)=-1
            end
        end
        model_choice=input('which model do you want to use :\n- Gimme \n- iMAT \n- tINIT \n- INIT\n- All\n' ,'s');
        while(ismember(['Gimme','IMAT','tINIT','INIT','All'],model_choice)==0)
            model_choice=input('Input not valid \n Please chose between this :n- Gimme \n- iMAT \n- tINIT \n- INIT\n','s')
        end
        file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
        file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');        
        switch (model_choice)
            case'Gimme'
                [Gimme_model_H_CALU_alt, Gimme_model_I_CALU_alt]=Generation_Gimme_alt(model_Healthy,CALU_Healthy_mean,essentialTasks_H,Cobra_infected_model,CALU_Infected_mean,essentialTasks_I);
                Gimme_model_H_CALU_alt.rules(206)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_CALU_alt.rxns{206})));
                writeCbModel(Gimme_model_H_CALU_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_CALU_alt, 'format','mat', 'fileName', file_name_Infected)
           
            case'IMAT'
                [IMAT_model_H_CALU_alt, iMAT_model_I_CALU_alt]=Generation_iMAT_alt(model_Healthy,CALU_Healthy_mean,essentialTasks_H,Cobra_infected_model,CALU_Infected_mean,essentialTasks_I);
                writeCbModel(IMAT_model_H_CALU_alt, 'format','mat', 'fileName', file_name_healthy)
                iMAT_model_I_CALU_alt.rules(2822)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_CALU_alt.rxns{2822})));
                writeCbModel(iMAT_model_I_CALU_alt, 'format','mat', 'fileName', file_name_Infected)

            case'tINIT'
                [tINIT_model_H_CALU,tINIT_model_I_CALU]=tINIT_generation(CALU_Healthy_mean_INIT,CALU_Infected_mean_INIT)
                tINIT_model_H_CALU=ravenCobraWrapper(tINIT_model_H_CALU)
                tINIT_model_I_CALU=rmfield(tINIT_model_I_CALU,'rules')
                tINIT_model_I_CALU=ravenCobraWrapper(tINIT_model_I_CALU)
                writeCbModel(tINIT_model_H_CALU, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_CALU, 'format','mat', 'fileName', file_name_Infected)
                
            case'INIT'
                [INIT_model_H_Lung,INIT_model_I_Lung]=INIT_generation_good(Lung_Healthy_mean_INIT,Lung_Infected_mean_INIT)
                INIT_model_H_CALU=ravenCobraWrapper(INIT_model_H_CALU)
                INIT_model_I_CALU=rmfield(INIT_model_I_CALU,'rules')
                INIT_model_I_CALU=ravenCobraWrapper(INIT_model_I_CALU)
                writeCbModel(INIT_model_H_CALU, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(INIT_model_I_CALU, 'format','mat', 'fileName', file_name_Infected)
                
            case'All'
                %Gimme
                model_choice='Gimme'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');    
                [Gimme_model_H_CALU_alt, Gimme_model_I_CALU_alt]=Generation_Gimme_alt(model_Healthy,CALU_Healthy_mean,essentialTasks_H,Cobra_infected_model,CALU_Infected_mean,essentialTasks_I);
                Gimme_model_H_CALU_alt.rules(206)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_CALU_alt.rxns{206})));
                writeCbModel(Gimme_model_H_CALU_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_CALU_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %IMAT
                model_choice='IMAT'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');    
                [IMAT_model_H_CALU_alt, iMAT_model_I_CALU_alt]=Generation_iMAT_alt(model_Healthy,CALU_Healthy_mean,essentialTasks_H,Cobra_infected_model,CALU_Infected_mean,essentialTasks_I);
                writeCbModel(IMAT_model_H_CALU_alt, 'format','mat', 'fileName', file_name_healthy)
                iMAT_model_I_CALU_alt.rules(2822)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_CALU_alt.rxns{2822})));
                writeCbModel(iMAT_model_I_CALU_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %tINIT
                model_choice='tINIT'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');    
                [tINIT_model_H_CALU,tINIT_model_I_CALU]=tINIT_generation(CALU_Healthy_mean_INIT,CALU_Infected_mean_INIT)
                tINIT_model_H_CALU=ravenCobraWrapper(tINIT_model_H_CALU)
                tINIT_model_I_CALU=rmfield(tINIT_model_I_CALU,'rules')
                tINIT_model_I_CALU=ravenCobraWrapper(tINIT_model_I_CALU)
                writeCbModel(tINIT_model_H_CALU, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_CALU, 'format','mat', 'fileName', file_name_Infected)
                
                %INIT
                model_choice='INIT'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');    
                [INIT_model_H_Lung,INIT_model_I_Lung]=INIT_generation_good(Lung_Healthy_mean_INIT,Lung_Infected_mean_INIT)
                INIT_model_H_CALU=ravenCobraWrapper(INIT_model_H_CALU)
                INIT_model_I_CALU=rmfield(INIT_model_I_CALU,'rules')
                INIT_model_I_CALU=ravenCobraWrapper(INIT_model_I_CALU)
                writeCbModel(INIT_model_H_CALU, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(INIT_model_I_CALU, 'format','mat', 'fileName', file_name_Infected)
                

        end
    case'Lung'
        load('Lung_Healthy_mean.mat');
        Lung_Healthy_mean=cell2mat(Lung_Healthy_mean);
        Lung_Healthy_mean_INIT=Lung_Healthy_mean;
        for p=1:length(Lung_Healthy_mean)
            if isnan(Lung_Healthy_mean(p))==1
                Lung_Healthy_mean(p)=-1
            end
        end
        load('Lung_Infected_mean.mat');
        Lung_Infected_mean=cell2mat(Lung_Infected_mean);
        Lung_Infected_mean_INIT=Lung_Infected_mean;
        for p=1:length(Lung_Infected_mean)
            if isnan(Lung_Infected_mean(p))==1
                Lung_Infected_mean(p)=-1
            end
        end
         model_choice=input('which model do you want to use :\n- Gimme \n- iMAT \n- tINIT \n- INIT\n- All\n' ,'s');
        while(ismember(['Gimme','IMAT','tINIT','INIT','All'],model_choice)==0)
            model_choice=input('Input not valid \n Please chose between this :n- Gimme \n- iMAT \n- tINIT \n- INIT\n','s')
        end
        file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
        file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');  
        switch(model_choice)
            case'Gimme'
                [Gimme_model_H_Lung_alt, Gimme_model_I_Lung_alt]=Generation_Gimme_alt(model_Healthy,Lung_Healthy_mean,essentialTasks_H,Cobra_infected_model,Lung_Infected_mean,essentialTasks_I); 
                Gimme_model_H_Lung_alt.rules(214)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_Lung_alt.rxns{214})));
                writeCbModel(Gimme_model_H_Lung_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_Lung_alt, 'format','mat', 'fileName', file_name_Infected)

            case'IMAT'
                [IMAT_model_H_Lung_alt, iMAT_model_I_Lung_alt]=Generation_iMAT_alt(model_Healthy,Lung_Healthy_mean,essentialTasks_H,Cobra_infected_model,Lung_Infected_mean,essentialTasks_I);
                writeCbModel(IMAT_model_H_Lung_alt, 'format','mat', 'fileName', file_name_healthy)
                iMAT_model_I_Lung_alt.rules(3121)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_Lung_alt.rxns{3121})));
                writeCbModel(iMAT_model_I_Lung_alt, 'format','mat', 'fileName', file_name_Infected)

            case'tINIT'
                [tINIT_model_H_Lung,tINIT_model_I_Lung]=tINIT_generation(Lung_Healthy_mean_INIT,Lung_Infected_mean_INIT)
                tINIT_model_H_Lung=ravenCobraWrapper(tINIT_model_H_Lung)
                tINIT_model_I_Lung=rmfield(tINIT_model_I_Lung,'rules')
                tINIT_model_I_Lung=ravenCobraWrapper(tINIT_model_I_Lung)
                writeCbModel(tINIT_model_H_Lung, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_Lung, 'format','mat', 'fileName', file_name_Infected)

            case'INIT'
                [INIT_model_H_Lung,INIT_model_I_Lung]=INIT_generation_good(Lung_Healthy_mean_INIT,Lung_Infected_mean_INIT)
                INIT_model_H_Lung=ravenCobraWrapper(INIT_model_H_Lung)
                writeCbModel(INIT_model_H_Lung, 'format','mat', 'fileName', 'INIT_model_H_Lung_alt.mat')
                INIT_model_I_Lung=rmfield(INIT_model_I_Lung,'rules')
                INIT_model_I_Lung=ravenCobraWrapper(INIT_model_I_Lung)
                writeCbModel(INIT_model_I_Lung, 'format','mat', 'fileName', 'INIT_model_I_Lung_alt.mat')
                
            case'All'
                %Gimme
                model_choice='Gimme'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');    
                [Gimme_model_H_Lung_alt, Gimme_model_I_Lung_alt]=Generation_Gimme_alt(model_Healthy,Lung_Healthy_mean,essentialTasks_H,Cobra_infected_model,Lung_Infected_mean,essentialTasks_I); 
                Gimme_model_H_Lung_alt.rules(214)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_Lung_alt.rxns{214})));
                writeCbModel(Gimme_model_H_Lung_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_Lung_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %IMAT
                model_choice='IMAT'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');    
                [IMAT_model_H_Lung_alt, iMAT_model_I_Lung_alt]=Generation_iMAT_alt(model_Healthy,Lung_Healthy_mean,essentialTasks_H,Cobra_infected_model,Lung_Infected_mean,essentialTasks_I);
                writeCbModel(IMAT_model_H_Lung_alt, 'format','mat', 'fileName', file_name_healthy)
                iMAT_model_I_Lung_alt.rules(3121)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_Lung_alt.rxns{3121})));
                writeCbModel(iMAT_model_I_Lung_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %tINIT
                model_choice='tINIT'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');    
                [tINIT_model_H_Lung,tINIT_model_I_Lung]=tINIT_generation(Lung_Healthy_mean_INIT,Lung_Infected_mean_INIT)
                tINIT_model_H_Lung=ravenCobraWrapper(tINIT_model_H_Lung)
                tINIT_model_I_Lung=rmfield(tINIT_model_I_Lung,'rules')
                tINIT_model_I_Lung=ravenCobraWrapper(tINIT_model_I_Lung)
                writeCbModel(tINIT_model_H_Lung, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_Lung, 'format','mat', 'fileName', file_name_Infected)
                
                %INIT
                model_choice='INIT'
                file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');    
                [INIT_model_H_Lung,INIT_model_I_Lung]=INIT_generation_good(Lung_Healthy_mean_INIT,Lung_Infected_mean_INIT)
                INIT_model_H_Lung=ravenCobraWrapper(INIT_model_H_Lung)
                writeCbModel(INIT_model_H_Lung, 'format','mat', 'fileName', 'INIT_model_H_Lung_alt.mat')
                INIT_model_I_Lung=rmfield(INIT_model_I_Lung,'rules')
                INIT_model_I_Lung=ravenCobraWrapper(INIT_model_I_Lung)
                writeCbModel(INIT_model_I_Lung, 'format','mat', 'fileName', 'INIT_model_I_Lung_alt.mat')
                

        end
    case'293T'
        load('T_Healthy_mean.mat');
        T_Healthy_mean=cell2mat(T_Healthy_mean);
        T_Healthy_mean_INIT=T_Healthy_mean;
        for p=1:length(T_Healthy_mean)
            if isnan(T_Healthy_mean(p))==1
                T_Healthy_mean(p)=-1
            end
        end
        load('T_Infected_mean.mat');
        T_Infected_mean=cell2mat(T_Infected_mean);
        T_Infected_mean_INIT=T_Infected_mean;
        for p=1:length(T_Infected_mean)
            if isnan(T_Infected_mean(p))==1
                T_Infected_mean(p)=-1
            end
        end
        
        model_choice=input('which model do you want to use :\n- Gimme \n- iMAT \n- tINIT \n- INIT\n- All\n' ,'s');
        while(ismember(['Gimme','IMAT','tINIT','INIT','All'],model_choice)==0)
            model_choice=input('Input not valid \n Please chose between this :n- Gimme \n- iMAT \n- tINIT \n- INIT\n','s')
        end
        file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
        file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');  
        switch(model_choice)
            case'Gimme'
                [Gimme_model_H_293T_alt, Gimme_model_I_293T_alt]=Generation_Gimme_alt(model_Healthy,T_Healthy_mean,essentialTasks_H,Cobra_infected_model,T_Infected_mean,essentialTasks_I);
                Gimme_model_H_293T_alt.rules(221)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_293T_alt.rxns{221})));
                writeCbModel(Gimme_model_H_293T_alt, 'format','mat', 'fileName', 'Gimme_model_H_293T_alt.mat')
                writeCbModel(Gimme_model_I_293T_alt, 'format','mat', 'fileName', 'Gimme_model_I_293T_alt.mat')

            case'IMAT'
                [IMAT_model_H_293T_alt, iMAT_model_I_293T_alt]=Generation_iMAT_alt(model_Healthy,T_Healthy_mean,essentialTasks_H,Cobra_infected_model,T_Infected_mean,essentialTasks_I);
                iMAT_model_I_293T_alt.rules(2854)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_293T_alt.rxns{2854})));
                writeCbModel(iMAT_model_I_293T_alt, 'format','mat', 'fileName', 'iMAT_model_I_293T_alt.mat')
                writeCbModel(IMAT_model_H_293T_alt, 'format','mat', 'fileName', 'iMAT_model_H_293T_alt.mat')

            case'tINIT'
                [tINIT_model_H_293T,tINIT_model_I_293T]=tINIT_generation(T_Healthy_mean_INIT,T_Infected_mean_INIT)
                tINIT_model_H_293T=ravenCobraWrapper(tINIT_model_H_293T)
                tINIT_model_I_293T=rmfield(tINIT_model_I_293T,'rules')
                tINIT_model_I_293T=ravenCobraWrapper(tINIT_model_I_293T)
                writeCbModel(tINIT_model_H_293T, 'format','mat', 'fileName', 'tINIT_model_H_293T_alt.mat')
                writeCbModel(tINIT_model_I_293T, 'format','mat', 'fileName', 'tINIT_model_I_293T_alt.mat')

            case'INIT'
                [INIT_model_H_293T,INIT_model_I_293T]=INIT_generation_good(T_Healthy_mean_INIT,T_Infected_mean_INIT)
                INIT_model_H_293T=ravenCobraWrapper(INIT_model_H_293T)
                INIT_model_I_293T=rmfield(INIT_model_I_293T,'rules')
                INIT_model_I_293T=ravenCobraWrapper(INIT_model_I_293T)
                writeCbModel(INIT_model_H_293T, 'format','mat', 'fileName', 'INIT_model_H_293T_alt.mat')
                writeCbModel(INIT_model_I_293T, 'format','mat', 'fileName', 'INIT_model_I_293T_alt.mat')
                
            case'All'
                %Gimme
                [Gimme_model_H_293T_alt, Gimme_model_I_293T_alt]=Generation_Gimme_alt(model_Healthy,T_Healthy_mean,essentialTasks_H,Cobra_infected_model,T_Infected_mean,essentialTasks_I);
                Gimme_model_H_293T_alt.rules(221)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_293T_alt.rxns{221})));
                writeCbModel(Gimme_model_H_293T_alt, 'format','mat', 'fileName', 'Gimme_model_H_293T_alt.mat')
                writeCbModel(Gimme_model_I_293T_alt, 'format','mat', 'fileName', 'Gimme_model_I_293T_alt.mat')
                
                %IMAT
                [IMAT_model_H_293T_alt, iMAT_model_I_293T_alt]=Generation_iMAT_alt(model_Healthy,T_Healthy_mean,essentialTasks_H,Cobra_infected_model,T_Infected_mean,essentialTasks_I);
                iMAT_model_I_293T_alt.rules(2854)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_293T_alt.rxns{2854})));
                writeCbModel(iMAT_model_I_293T_alt, 'format','mat', 'fileName', 'iMAT_model_I_293T_alt.mat')
                writeCbModel(IMAT_model_H_293T_alt, 'format','mat', 'fileName', 'iMAT_model_H_293T_alt.mat')
                
                %tINIT
                [tINIT_model_H_293T,tINIT_model_I_293T]=tINIT_generation(T_Healthy_mean_INIT,T_Infected_mean_INIT)
                tINIT_model_H_293T=ravenCobraWrapper(tINIT_model_H_293T)
                tINIT_model_I_293T=rmfield(tINIT_model_I_293T,'rules')
                tINIT_model_I_293T=ravenCobraWrapper(tINIT_model_I_293T)
                writeCbModel(tINIT_model_H_293T, 'format','mat', 'fileName', 'tINIT_model_H_293T_alt.mat')
                writeCbModel(tINIT_model_I_293T, 'format','mat', 'fileName', 'tINIT_model_I_293T_alt.mat')
                
                %INIT
                [INIT_model_H_293T,INIT_model_I_293T]=INIT_generation_good(T_Healthy_mean_INIT,T_Infected_mean_INIT)
                INIT_model_H_293T=ravenCobraWrapper(INIT_model_H_293T)
                INIT_model_I_293T=rmfield(INIT_model_I_293T,'rules')
                INIT_model_I_293T=ravenCobraWrapper(INIT_model_I_293T)
                writeCbModel(INIT_model_H_293T, 'format','mat', 'fileName', 'INIT_model_H_293T_alt.mat')
                writeCbModel(INIT_model_I_293T, 'format','mat', 'fileName', 'INIT_model_I_293T_alt.mat')
     
        end
    case'All'
        load('T_Healthy_mean.mat');
        T_Healthy_mean=cell2mat(T_Healthy_mean);
        T_Healthy_mean_INIT=T_Healthy_mean;
        for p=1:length(T_Healthy_mean)
            if isnan(T_Healthy_mean(p))==1
                T_Healthy_mean(p)=-1
            end
        end
        load('T_Infected_mean.mat');
        T_Infected_mean=cell2mat(T_Infected_mean);
        T_Infected_mean_INIT=T_Infected_mean;
        for p=1:length(T_Infected_mean)
            if isnan(T_Infected_mean(p))==1
                T_Infected_mean(p)=-1
            end
        end
        
        load('A549_Healthy_mean.mat');
        A549_Healthy_mean=cell2mat(A549_Healthy_mean);
        A549_Healthy_mean_INIT=A549_Healthy_mean;
        for p=1:length(A549_Healthy_mean)
            if isnan(A549_Healthy_mean(p))==1
                A549_Healthy_mean(p)=-1
            end
        end
        load('A549_Infected_mean.mat');
        A549_Infected_mean=cell2mat(A549_Infected_mean);
        A549_Infected_mean_INIT=A549_Infected_mean;
        for p=1:length(A549_Infected_mean)
            if isnan(A549_Infected_mean(p))==1
                A549_Infected_mean(p)=-1
            end
        end
        
        load('CALU_Healthy_mean.mat');
        CALU_Healthy_mean=cell2mat(CALU_Healthy_mean);
        CALU_Healthy_mean_INIT=CALU_Healthy_mean;
        for p=1:length(CALU_Healthy_mean)
            if isnan(CALU_Healthy_mean(p))==1
                CALU_Healthy_mean(p)=-1
            end
        end
        load('CALU_Infected_mean.mat');
        CALU_Infected_mean=cell2mat(CALU_Infected_mean);
        CALU_Infected_mean_INIT=CALU_Infected_mean;
        for p=1:length(CALU_Infected_mean)
            if isnan(CALU_Infected_mean(p))==1
                CALU_Infected_mean(p)=-1
            end
        end
        
        load('Lung_Healthy_mean.mat');
        Lung_Healthy_mean=cell2mat(Lung_Healthy_mean);
        Lung_Healthy_mean_INIT=Lung_Healthy_mean;
        for p=1:length(Lung_Healthy_mean)
            if isnan(Lung_Healthy_mean(p))==1
                Lung_Healthy_mean(p)=-1
            end
        end
        load('Lung_Infected_mean.mat');
        Lung_Infected_mean=cell2mat(Lung_Infected_mean);
        Lung_Infected_mean_INIT=Lung_Infected_mean;
        for p=1:length(Lung_Infected_mean)
            if isnan(Lung_Infected_mean(p))==1
                Lung_Infected_mean(p)=-1
            end
        end
        
        load('NHBE_Healthy_mean.mat');
        NHBE_Healthy_mean=cell2mat(NHBE_Healthy_mean);
        NHBE_Healthy_mean_INIT=NHBE_Healthy_mean;
        for k=1:length(NHBE_Healthy_mean_INIT)
            if NHBE_Healthy_mean_INIT(k)==-1
                NHBE_Healthy_mean_INIT(k)=NaN;
            end
        end
                        
        load('NHBE_Infected_mean.mat');
        NHBE_Infected_mean=cell2mat(NHBE_Infected_mean);
        NHBE_Infected_mean_INIT=NHBE_Infected_mean;
        for k=1:length(NHBE_Infected_mean_INIT)
            if NHBE_Infected_mean_INIT(k)==-1
                NHBE_Infected_mean_INIT(k)=NaN
            end
        end
        model_choice=input('which model do you want to use :\n- Gimme \n- iMAT \n- tINIT \n- INIT\n- All\n' ,'s');
        while(ismember(['Gimme','IMAT','tINIT','INIT','All'],model_choice)==0)
            model_choice=input('Input not valid \n Please chose between this :n- Gimme \n- iMAT \n- tINIT \n- INIT\n','s')
        end
        file_name_healthy=strcat(model_choice,'_model_H_',Database,'_alt.mat');
        file_name_Infected=strcat(model_choice,'_model_I_',Database,'_alt.mat');  
        switch(model_choice)
            case'Gimme'
                %NHBE
                Database2='NHBE'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');    
                [Gimme_model_H_NHBE_alt, Gimme_model_I_NHBE_alt]=Generation_Gimme_alt(model_Healthy,NHBE_Healthy_mean,essentialTasks_H,Cobra_infected_model,NHBE_Infected_mean,essentialTasks_I);
                Gimme_model_H_NHBE_alt.rules(204)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_NHBE_alt.rxns{204})));
                writeCbModel(Gimme_model_H_NHBE_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_NHBE_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %A549
                Database2='A549'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');    
                [Gimme_model_H_A549_alt, Gimme_model_I_A549_alt]=Generation_Gimme_alt(model_Healthy,A549_Healthy_mean,essentialTasks_H,Cobra_infected_model,A549_Infected_mean,essentialTasks_I);               
                Gimme_model_H_A549_alt.rules(206)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_A549_alt.rxns{206})));
                writeCbModel(Gimme_model_H_A549_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_A549_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %CAlU
                Database2='CALU'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');
                [Gimme_model_H_CALU_alt, Gimme_model_I_CALU_alt]=Generation_Gimme_alt(model_Healthy,CALU_Healthy_mean,essentialTasks_H,Cobra_infected_model,CALU_Infected_mean,essentialTasks_I);
                Gimme_model_H_CALU_alt.rules(206)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_CALU_alt.rxns{206})));
                writeCbModel(Gimme_model_H_CALU_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_CALU_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %Lung
                Database2='Lung'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');
                [Gimme_model_H_Lung_alt, Gimme_model_I_Lung_alt]=Generation_Gimme_alt(model_Healthy,Lung_Healthy_mean,essentialTasks_H,Cobra_infected_model,Lung_Infected_mean,essentialTasks_I); 
                Gimme_model_H_Lung_alt.rules(214)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_Lung_alt.rxns{214})));
                writeCbModel(Gimme_model_H_Lung_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_Lung_alt, 'format','mat', 'fileName', file_name_Infected)

                %293T
                Database2='293T'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');
                [Gimme_model_H_293T_alt, Gimme_model_I_293T_alt]=Generation_Gimme_alt(model_Healthy,T_Healthy_mean,essentialTasks_H,Cobra_infected_model,T_Infected_mean,essentialTasks_I);
                Gimme_model_H_293T_alt.rules(221)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_293T_alt.rxns{221})));
                writeCbModel(Gimme_model_H_293T_alt, 'format','mat', 'fileName', 'Gimme_model_H_293T_alt.mat')
                writeCbModel(Gimme_model_I_293T_alt, 'format','mat', 'fileName', 'Gimme_model_I_293T_alt.mat')
                
            case'IMAT'
                %NHBE
                Database2='NHBE'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [iMAT_model_H_NHBE_alt, iMAT_model_I_NHBE_alt]=Generation_iMAT_alt(model_Healthy,NHBE_Healthy_mean,essentialTasks_H,Cobra_infected_model,NHBE_Infected_mean,essentialTasks_I);
                iMAT_model_I_NHBE_alt.rules(2872)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),iMAT_model_I_NHBE_alt.rxns{2872})));   
                writeCbModel(iMAT_model_H_NHBE_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(iMAT_model_I_NHBE_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %A549
                Database2='A549'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [IMAT_model_H_A549_alt, iMAT_model_I_A549_alt]=Generation_iMAT_alt(model_Healthy,A549_Healthy_mean,essentialTasks_H,Cobra_infected_model,A549_Infected_mean,essentialTasks_I);
                iMAT_model_I_A549_alt.rules(2889)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_A549_alt.rxns{2889})));
                writeCbModel(IMAT_model_H_A549_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(iMAT_model_I_A549_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %CALU
                Database2='CALU'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [IMAT_model_H_CALU_alt, iMAT_model_I_CALU_alt]=Generation_iMAT_alt(model_Healthy,CALU_Healthy_mean,essentialTasks_H,Cobra_infected_model,CALU_Infected_mean,essentialTasks_I);
                writeCbModel(IMAT_model_H_CALU_alt, 'format','mat', 'fileName', file_name_healthy)
                iMAT_model_I_CALU_alt.rules(2822)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_CALU_alt.rxns{2822})));
                writeCbModel(iMAT_model_I_CALU_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %Lung
                Database2='Lung'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [IMAT_model_H_Lung_alt, iMAT_model_I_Lung_alt]=Generation_iMAT_alt(model_Healthy,Lung_Healthy_mean,essentialTasks_H,Cobra_infected_model,Lung_Infected_mean,essentialTasks_I);
                writeCbModel(IMAT_model_H_Lung_alt, 'format','mat', 'fileName', file_name_healthy)
                iMAT_model_I_Lung_alt.rules(3121)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_Lung_alt.rxns{3121})));
                writeCbModel(iMAT_model_I_Lung_alt, 'format','mat', 'fileName', file_name_Infected)

                %293T
                Database2='293T'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [IMAT_model_H_293T_alt, iMAT_model_I_293T_alt]=Generation_iMAT_alt(model_Healthy,T_Healthy_mean,essentialTasks_H,Cobra_infected_model,T_Infected_mean,essentialTasks_I);
                iMAT_model_I_293T_alt.rules(2854)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_293T_alt.rxns{2854})));
                writeCbModel(iMAT_model_I_293T_alt, 'format','mat', 'fileName', 'iMAT_model_I_293T_alt.mat')
                writeCbModel(IMAT_model_H_293T_alt, 'format','mat', 'fileName', 'iMAT_model_H_293T_alt.mat')
                
            case'tINIT'
                %NHBE
                Database2='NHBE'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [tINIT_model_H_NHBE,tINIT_model_I_NHBE]=tINIT_generation(NHBE_Healthy_mean_INIT,NHBE_Infected_mean_INIT)
                tINIT_model_H_NHBE=ravenCobraWrapper(tINIT_model_H_NHBE);
                tINIT_model_I_NHBE=rmfield(tINIT_model_I_NHBE,'rules');
                tINIT_model_I_NHBE=ravenCobraWrapper(tINIT_model_I_NHBE);
                writeCbModel(tINIT_model_H_NHBE, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_NHBE, 'format','mat', 'fileName', file_name_Infected)
                
                %A549
                Database2='A549'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [tINIT_model_H_A549,tINIT_model_I_A549]=tINIT_generation(A549_Healthy_mean_INIT,A549_Infected_mean_INIT)
                tINIT_model_I_A549=rmfield(tINIT_model_I_A549,'rules')
                tINIT_model_I_A549=ravenCobraWrapper(tINIT_model_I_A549)
                writeCbModel(tINIT_model_H_A549, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_A549, 'format','mat', 'fileName', file_name_Infected)
                
                %CALU
                Database2='CALU'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [tINIT_model_H_CALU,tINIT_model_I_CALU]=tINIT_generation(CALU_Healthy_mean_INIT,CALU_Infected_mean_INIT)
                tINIT_model_H_CALU=ravenCobraWrapper(tINIT_model_H_CALU)
                tINIT_model_I_CALU=rmfield(tINIT_model_I_CALU,'rules')
                tINIT_model_I_CALU=ravenCobraWrapper(tINIT_model_I_CALU)
                writeCbModel(tINIT_model_H_CALU, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_CALU, 'format','mat', 'fileName', file_name_Infected)
                
                %Lung
                Database2='Lung'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [tINIT_model_H_Lung,tINIT_model_I_Lung]=tINIT_generation(Lung_Healthy_mean_INIT,Lung_Infected_mean_INIT)
                tINIT_model_H_Lung=ravenCobraWrapper(tINIT_model_H_Lung)
                tINIT_model_I_Lung=rmfield(tINIT_model_I_Lung,'rules')
                tINIT_model_I_Lung=ravenCobraWrapper(tINIT_model_I_Lung)
                writeCbModel(tINIT_model_H_Lung, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_Lung, 'format','mat', 'fileName', file_name_Infected)

                %293T
                [tINIT_model_H_293T,tINIT_model_I_293T]=tINIT_generation(T_Healthy_mean_INIT,T_Infected_mean_INIT)
                tINIT_model_H_293T=ravenCobraWrapper(tINIT_model_H_293T)
                tINIT_model_I_293T=rmfield(tINIT_model_I_293T,'rules')
                tINIT_model_I_293T=ravenCobraWrapper(tINIT_model_I_293T)
                writeCbModel(tINIT_model_H_293T, 'format','mat', 'fileName', 'tINIT_model_H_293T_alt.mat')
                writeCbModel(tINIT_model_I_293T, 'format','mat', 'fileName', 'tINIT_model_I_293T_alt.mat')
                
                
            case'INIT'
                %NHBE
                Database2='NHBE'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [INIT_model_H_NHBE,INIT_model_I_NHBE]=INIT_generation_good(NHBE_Healthy_mean_INIT,NHBE_Infected_mean_INIT)
                INIT_model_H_NHBE=ravenCobraWrapper(INIT_model_H_NHBE);
                INIT_model_I_NHBE=rmfield(INIT_model_I_NHBE,'rules');
                INIT_model_I_NHBE=ravenCobraWrapper(INIT_model_I_NHBE);
                writeCbModel(INIT_model_H_NHBE, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(INIT_model_I_NHBE, 'format','mat', 'fileName', file_name_Infected)
                
                %A549
                Database2='A549'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [INIT_model_H_A549,INIT_model_I_A549]=INIT_generation_good(A549_Healthy_mean_INIT,A549_Infected_mean_INIT)
                INIT_model_H_A549=ravenCobraWrapper(INIT_model_H_A549)
                writeCbModel(INIT_model_H_A549, 'format','mat', 'fileName', 'INIT_model_H_A549_alt.mat')
                INIT_model_I_A549=rmfield(INIT_model_I_A549,'rules')
                INIT_model_I_A549=ravenCobraWrapper(INIT_model_I_A549)
                writeCbModel(INIT_model_I_A549, 'format','mat', 'fileName', 'INIT_model_I_A549_alt.mat')
                
                %CALU
                Database2='CALI'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [INIT_model_H_Lung,INIT_model_I_Lung]=INIT_generation_good(Lung_Healthy_mean_INIT,Lung_Infected_mean_INIT)
                INIT_model_H_CALU=ravenCobraWrapper(INIT_model_H_CALU)
                INIT_model_I_CALU=rmfield(INIT_model_I_CALU,'rules')
                INIT_model_I_CALU=ravenCobraWrapper(INIT_model_I_CALU)
                writeCbModel(INIT_model_H_CALU, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(INIT_model_I_CALU, 'format','mat', 'fileName', file_name_Infected)
                
                %Lung
                [INIT_model_H_Lung,INIT_model_I_Lung]=INIT_generation_good(Lung_Healthy_mean_INIT,Lung_Infected_mean_INIT)
                INIT_model_H_Lung=ravenCobraWrapper(INIT_model_H_Lung)
                writeCbModel(INIT_model_H_Lung, 'format','mat', 'fileName', 'INIT_model_H_Lung_alt.mat')
                INIT_model_I_Lung=rmfield(INIT_model_I_Lung,'rules')
                INIT_model_I_Lung=ravenCobraWrapper(INIT_model_I_Lung)
                writeCbModel(INIT_model_I_Lung, 'format','mat', 'fileName', 'INIT_model_I_Lung_alt.mat')
                
               %293T
               [INIT_model_H_293T,INIT_model_I_293T]=INIT_generation_good(T_Healthy_mean_INIT,T_Infected_mean_INIT)
               INIT_model_H_293T=ravenCobraWrapper(INIT_model_H_293T)
               INIT_model_I_293T=rmfield(INIT_model_I_293T,'rules')
               INIT_model_I_293T=ravenCobraWrapper(INIT_model_I_293T)
               writeCbModel(INIT_model_H_293T, 'format','mat', 'fileName', 'INIT_model_H_293T_alt.mat')
               writeCbModel(INIT_model_I_293T, 'format','mat', 'fileName', 'INIT_model_I_293T_alt.mat')
                
            case'ALL'
                %'Gimme'
                %NHBE
                Database2='NHBE'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');    
                [Gimme_model_H_NHBE_alt, Gimme_model_I_NHBE_alt]=Generation_Gimme_alt(model_Healthy,NHBE_Healthy_mean,essentialTasks_H,Cobra_infected_model,NHBE_Infected_mean,essentialTasks_I);
                Gimme_model_H_NHBE_alt.rules(204)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_NHBE_alt.rxns{204})));
                writeCbModel(Gimme_model_H_NHBE_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_NHBE_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %A549
                Database2='A549'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');    
                [Gimme_model_H_A549_alt, Gimme_model_I_A549_alt]=Generation_Gimme_alt(model_Healthy,A549_Healthy_mean,essentialTasks_H,Cobra_infected_model,A549_Infected_mean,essentialTasks_I);               
                Gimme_model_H_A549_alt.rules(206)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_A549_alt.rxns{206})));
                writeCbModel(Gimme_model_H_A549_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_A549_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %CAlU
                Database2='CALU'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');
                [Gimme_model_H_CALU_alt, Gimme_model_I_CALU_alt]=Generation_Gimme_alt(model_Healthy,CALU_Healthy_mean,essentialTasks_H,Cobra_infected_model,CALU_Infected_mean,essentialTasks_I);
                Gimme_model_H_CALU_alt.rules(206)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_CALU_alt.rxns{206})));
                writeCbModel(Gimme_model_H_CALU_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_CALU_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %Lung
                Database2='Lung'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');
                [Gimme_model_H_Lung_alt, Gimme_model_I_Lung_alt]=Generation_Gimme_alt(model_Healthy,Lung_Healthy_mean,essentialTasks_H,Cobra_infected_model,Lung_Infected_mean,essentialTasks_I); 
                Gimme_model_H_Lung_alt.rules(214)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_Lung_alt.rxns{214})));
                writeCbModel(Gimme_model_H_Lung_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(Gimme_model_I_Lung_alt, 'format','mat', 'fileName', file_name_Infected)

                %293T
                Database2='293T'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');
                [Gimme_model_H_293T_alt, Gimme_model_I_293T_alt]=Generation_Gimme_alt(model_Healthy,T_Healthy_mean,essentialTasks_H,Cobra_infected_model,T_Infected_mean,essentialTasks_I);
                Gimme_model_H_293T_alt.rules(221)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),Gimme_model_H_293T_alt.rxns{221})));
                writeCbModel(Gimme_model_H_293T_alt, 'format','mat', 'fileName', 'Gimme_model_H_293T_alt.mat')
                writeCbModel(Gimme_model_I_293T_alt, 'format','mat', 'fileName', 'Gimme_model_I_293T_alt.mat')
                
            %'IMAT'
                %NHBE
                Database2='NHBE'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [iMAT_model_H_NHBE_alt, iMAT_model_I_NHBE_alt]=Generation_iMAT_alt(model_Healthy,NHBE_Healthy_mean,essentialTasks_H,Cobra_infected_model,NHBE_Infected_mean,essentialTasks_I);
                iMAT_model_I_NHBE_alt.rules(2872)=model_Healthy.rules(find(ismember(model_Healthy.rxns(:,1),iMAT_model_I_NHBE_alt.rxns{2872})));   
                writeCbModel(iMAT_model_H_NHBE_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(iMAT_model_I_NHBE_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %A549
                Database2='A549'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [IMAT_model_H_A549_alt, iMAT_model_I_A549_alt]=Generation_iMAT_alt(model_Healthy,A549_Healthy_mean,essentialTasks_H,Cobra_infected_model,A549_Infected_mean,essentialTasks_I);
                iMAT_model_I_A549_alt.rules(2889)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_A549_alt.rxns{2889})));
                writeCbModel(IMAT_model_H_A549_alt, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(iMAT_model_I_A549_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %CALU
                Database2='CALU'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [IMAT_model_H_CALU_alt, iMAT_model_I_CALU_alt]=Generation_iMAT_alt(model_Healthy,CALU_Healthy_mean,essentialTasks_H,Cobra_infected_model,CALU_Infected_mean,essentialTasks_I);
                writeCbModel(IMAT_model_H_CALU_alt, 'format','mat', 'fileName', file_name_healthy)
                iMAT_model_I_CALU_alt.rules(2822)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_CALU_alt.rxns{2822})));
                writeCbModel(iMAT_model_I_CALU_alt, 'format','mat', 'fileName', file_name_Infected)
                
                %Lung
                Database2='Lung'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [IMAT_model_H_Lung_alt, iMAT_model_I_Lung_alt]=Generation_iMAT_alt(model_Healthy,Lung_Healthy_mean,essentialTasks_H,Cobra_infected_model,Lung_Infected_mean,essentialTasks_I);
                writeCbModel(IMAT_model_H_Lung_alt, 'format','mat', 'fileName', file_name_healthy)
                iMAT_model_I_Lung_alt.rules(3121)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_Lung_alt.rxns{3121})));
                writeCbModel(iMAT_model_I_Lung_alt, 'format','mat', 'fileName', file_name_Infected)

                %293T
                Database2='293T'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [IMAT_model_H_293T_alt, iMAT_model_I_293T_alt]=Generation_iMAT_alt(model_Healthy,T_Healthy_mean,essentialTasks_H,Cobra_infected_model,T_Infected_mean,essentialTasks_I);
                iMAT_model_I_293T_alt.rules(2854)=Cobra_infected_model.rules(find(ismember(Cobra_infected_model.rxns(:,1),iMAT_model_I_293T_alt.rxns{2854})));
                writeCbModel(iMAT_model_I_293T_alt, 'format','mat', 'fileName', 'iMAT_model_I_293T_alt.mat')
                writeCbModel(IMAT_model_H_293T_alt, 'format','mat', 'fileName', 'iMAT_model_H_293T_alt.mat')
                
            %'tINIT'
                %NHBE
                Database2='NHBE'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [tINIT_model_H_NHBE,tINIT_model_I_NHBE]=tINIT_generation(NHBE_Healthy_mean_INIT,NHBE_Infected_mean_INIT)
                tINIT_model_H_NHBE=ravenCobraWrapper(tINIT_model_H_NHBE);
                tINIT_model_I_NHBE=rmfield(tINIT_model_I_NHBE,'rules');
                tINIT_model_I_NHBE=ravenCobraWrapper(tINIT_model_I_NHBE);
                writeCbModel(tINIT_model_H_NHBE, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_NHBE, 'format','mat', 'fileName', file_name_Infected)
                
                %A549
                Database2='A549'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [tINIT_model_H_A549,tINIT_model_I_A549]=tINIT_generation(A549_Healthy_mean_INIT,A549_Infected_mean_INIT)
                tINIT_model_I_A549=rmfield(tINIT_model_I_A549,'rules')
                tINIT_model_I_A549=ravenCobraWrapper(tINIT_model_I_A549)
                writeCbModel(tINIT_model_H_A549, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_A549, 'format','mat', 'fileName', file_name_Infected)
                
                %CALU
                Database2='CALU'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [tINIT_model_H_CALU,tINIT_model_I_CALU]=tINIT_generation(CALU_Healthy_mean_INIT,CALU_Infected_mean_INIT)
                tINIT_model_H_CALU=ravenCobraWrapper(tINIT_model_H_CALU)
                tINIT_model_I_CALU=rmfield(tINIT_model_I_CALU,'rules')
                tINIT_model_I_CALU=ravenCobraWrapper(tINIT_model_I_CALU)
                writeCbModel(tINIT_model_H_CALU, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_CALU, 'format','mat', 'fileName', file_name_Infected)
                
                %Lung
                Database2='Lung'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [tINIT_model_H_Lung,tINIT_model_I_Lung]=tINIT_generation(Lung_Healthy_mean_INIT,Lung_Infected_mean_INIT)
                tINIT_model_H_Lung=ravenCobraWrapper(tINIT_model_H_Lung)
                tINIT_model_I_Lung=rmfield(tINIT_model_I_Lung,'rules')
                tINIT_model_I_Lung=ravenCobraWrapper(tINIT_model_I_Lung)
                writeCbModel(tINIT_model_H_Lung, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(tINIT_model_I_Lung, 'format','mat', 'fileName', file_name_Infected)

                %293T
                [tINIT_model_H_293T,tINIT_model_I_293T]=tINIT_generation(T_Healthy_mean_INIT,T_Infected_mean_INIT)
                tINIT_model_H_293T=ravenCobraWrapper(tINIT_model_H_293T)
                tINIT_model_I_293T=rmfield(tINIT_model_I_293T,'rules')
                tINIT_model_I_293T=ravenCobraWrapper(tINIT_model_I_293T)
                writeCbModel(tINIT_model_H_293T, 'format','mat', 'fileName', 'tINIT_model_H_293T_alt.mat')
                writeCbModel(tINIT_model_I_293T, 'format','mat', 'fileName', 'tINIT_model_I_293T_alt.mat')
                
                
            %'INIT'
                %NHBE
                Database2='NHBE'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [INIT_model_H_NHBE,INIT_model_I_NHBE]=INIT_generation_good(NHBE_Healthy_mean_INIT,NHBE_Infected_mean_INIT)
                INIT_model_H_NHBE=ravenCobraWrapper(INIT_model_H_NHBE);
                INIT_model_I_NHBE=rmfield(INIT_model_I_NHBE,'rules');
                INIT_model_I_NHBE=ravenCobraWrapper(INIT_model_I_NHBE);
                writeCbModel(INIT_model_H_NHBE, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(INIT_model_I_NHBE, 'format','mat', 'fileName', file_name_Infected)
                
                %A549
                Database2='A549'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [INIT_model_H_A549,INIT_model_I_A549]=INIT_generation_good(A549_Healthy_mean_INIT,A549_Infected_mean_INIT)
                INIT_model_H_A549=ravenCobraWrapper(INIT_model_H_A549)
                writeCbModel(INIT_model_H_A549, 'format','mat', 'fileName', 'INIT_model_H_A549_alt.mat')
                INIT_model_I_A549=rmfield(INIT_model_I_A549,'rules')
                INIT_model_I_A549=ravenCobraWrapper(INIT_model_I_A549)
                writeCbModel(INIT_model_I_A549, 'format','mat', 'fileName', 'INIT_model_I_A549_alt.mat')
                
                %CALU
                Database2='CALI'
                file_name_healthy=strcat(model_choice,'_model_H_',Database2,'_alt.mat');
                file_name_Infected=strcat(model_choice,'_model_I_',Database2,'_alt.mat');  
                [INIT_model_H_Lung,INIT_model_I_Lung]=INIT_generation_good(Lung_Healthy_mean_INIT,Lung_Infected_mean_INIT)
                INIT_model_H_CALU=ravenCobraWrapper(INIT_model_H_CALU)
                INIT_model_I_CALU=rmfield(INIT_model_I_CALU,'rules')
                INIT_model_I_CALU=ravenCobraWrapper(INIT_model_I_CALU)
                writeCbModel(INIT_model_H_CALU, 'format','mat', 'fileName', file_name_healthy)
                writeCbModel(INIT_model_I_CALU, 'format','mat', 'fileName', file_name_Infected)
                
                %Lung
                [INIT_model_H_Lung,INIT_model_I_Lung]=INIT_generation_good(Lung_Healthy_mean_INIT,Lung_Infected_mean_INIT)
                INIT_model_H_Lung=ravenCobraWrapper(INIT_model_H_Lung)
                writeCbModel(INIT_model_H_Lung, 'format','mat', 'fileName', 'INIT_model_H_Lung_alt.mat')
                INIT_model_I_Lung=rmfield(INIT_model_I_Lung,'rules')
                INIT_model_I_Lung=ravenCobraWrapper(INIT_model_I_Lung)
                writeCbModel(INIT_model_I_Lung, 'format','mat', 'fileName', 'INIT_model_I_Lung_alt.mat')
                
               %293T
               [INIT_model_H_293T,INIT_model_I_293T]=INIT_generation_good(T_Healthy_mean_INIT,T_Infected_mean_INIT)
               INIT_model_H_293T=ravenCobraWrapper(INIT_model_H_293T)
               INIT_model_I_293T=rmfield(INIT_model_I_293T,'rules')
               INIT_model_I_293T=ravenCobraWrapper(INIT_model_I_293T)
               writeCbModel(INIT_model_H_293T, 'format','mat', 'fileName', 'INIT_model_H_293T_alt.mat')
               writeCbModel(INIT_model_I_293T, 'format','mat', 'fileName', 'INIT_model_I_293T_alt.mat')
                         
        end
        
        
        


        
        
        
        
        
        
        
        
        
        
        