clc,clearvars,close all
%initCobraToolbox

%% Model
%Set solver, in this case, gurobi
changeCobraSolver('gurobi');

%Load model of Bifidobaterium infantis ATCC 15697
load('iLR578.mat')

%Add experimentally observed ratios for ac:lac (uncomment for simulations)
% model = addRatioReaction(model, {'EX_ac(e)','EX_lac_L(e)'}, [2 2.5]);
% model = addRatioReaction(model, {'EX_ac(e)','EX_lac_L(e)'}, [2 3.0]);
model = addRatioReaction(model, {'EX_ac(e)','EX_lac_L(e)'}, [2 3.5]);

%% Computation of maximum yields
% Define objective function for the analysis: 1) AGORA, 2) Schopping et al.
biomassEQ=1;                          
switch biomassEQ
    case 1
        bioRxn='Biomass_growth_AGORA';
    case 2
        bioRxn='Biomass_growth_Schopping';
end
model=changeObjective(model,bioRxn);

%PM (g/mol)
PM_Csource=[342.3, 488.39, 707.6, 633.6];
PM_mets=[76.09, 59.044, 46.07, 46.03, 90.08, 118.09]';

%Parameters and data
EX_rxns = {'EX_lcts(e)','EX_Lacto_N_Neotetraose(e)','EX_3fucosyllactose(e)',...
    'EX_6sialyllactose(e)'};
rS=[1,0.48,0.7,0.54];                                   % lactose-equivalent substrate conversion
mu_exp=[0.171, 0.098, 0.031, 0.043];                    % exp specific growth rates (1/h)
Y_exp=[0.0096, 0.0104, 0.0105, 0.0083];                 % exp yields (gdcw/mmol)
rS_exp=mu_exp./Y_exp;                                   % exp rates (mmol/gdcw/h)

%Initialize variables
yields=zeros(4,3);
yields(:,1)=1e3*Y_exp./PM_Csource;                              % gdcw/g
maxMu=zeros(4,3);
maxMu(:,1)=mu_exp;
met_yields=zeros(6,4);
max_yields=zeros(6,4);
min_yields=zeros(6,4);

% Optimality parameter for FVA
alpha = 0.99;

% Threshold for GIMME
prob = 0.3;

% Additionally observed ratios
model = addRatioReaction(model, {'EX_ac(e)','EX_etoh(e)'}, [1 100]);
model = addRatioReaction(model, {'EX_ac(e)','EX_for(e)'}, [1 100]);

%Mets names
Met_names={'EX_12ppd_S(e)';'EX_ac(e)';'EX_etoh(e)';'EX_for(e)';...
    'EX_lac_L(e)';'EX_succ(e)'};
model=changeRxnBounds(model,Met_names,0,'l');

% Main loop for calculations before expression data integration
for ix = 1:4
    modelTemp=changeRxnBounds(model,EX_rxns{ix},-rS(ix));
    fba_solution=optimizeCbModel(modelTemp);
    yields(ix,2)=fba_solution.f/(rS_exp(ix)/1e3)/PM_Csource(ix);
    maxMu(ix,2)=fba_solution.f;
    met_yields(:,ix)=PM_mets.*(fba_solution.x(ismember(model.rxns,Met_names))./...
        fba_solution.f)/1e3;
    modelTemp=changeRxnBounds(modelTemp,bioRxn,alpha*fba_solution.f,'b');
    max_sol_temp=zeros(6,1);
    min_sol_temp=zeros(6,1);
    for jx = 1:6
        modelTemp=changeObjective(modelTemp,Met_names{jx});
        max_solution=optimizeCbModel(modelTemp,'max');
        max_sol_temp(jx)=max_solution.f;
        min_solution=optimizeCbModel(modelTemp,'min');
        min_sol_temp(jx)=min_solution.f;
    end
    max_yields(:,ix)=PM_mets.*(max_sol_temp/...
        fba_solution.f)/1e3;
    min_yields(:,ix)=PM_mets.*(min_sol_temp/...
        fba_solution.f)/1e3;
end

% Save values
LCT_yields=met_yields(:,1);
LNNT_yields=met_yields(:,2);
FL3_yields=met_yields(:,3);
SL6_yields=met_yields(:,4);
Table2=table(Met_names,LCT_yields,LNNT_yields,FL3_yields,SL6_yields);

%% Expression data integration using GIMME
%Make a matrix from transcriptomic information in excel
[n,s] = xlsread('./data/Normalized_RNA_counts_final.xlsx','normalized_counts_final');
%Take the names of the genes from de table
features = s(2:end,1);
%Make a list with highly expressed genes (log10(rpkm)>3,5) w/o considering
%ribosomal proteins and tRNA 
he={'Blon_0235','Blon_0307','Blon_0308',...
    'Blon_0309','Blon_0310','Blon_0382',...
    'Blon_0612','Blon_0679','Blon_0694',...
    'Blon_0840','Blon_0900','Blon_1089',...
    'Blon_1095','Blon_1722','Blon_1738',...
    'Blon_1836','Blon_1922','Blon_1923',...
    'Blon_2152','Blon_2211','Blon_2215','Blon_2257'};
%Make a list with the rows in which these genes are
pos_he = find(ismember(features, he));
%Make a matrix with the highly expressed genes in LAC, HMO early, HMO mid,
%HMO late, LNnT, 3FL, 6SL (columns in that order)
norm_counts_he = n(pos_he, :);

%Use GIMME to minimize the use of low expression reactions
%First, intersect the list of genes with the genes mapped onto the model
[inter, posa, posb] = intersect(features, model.genes);
%Change the lower bound of lactose to make it the only carbon source
model1=changeRxnBounds(model,'EX_lcts(e)',-rS(1),'l');
%Optimize de model with this carbon source
maxLCTS=optimizeCbModel(model1,'max');
minLCTS=optimizeCbModel(model1,'min');
%Make a list with the normalized gene counts in model (as indexed when you
%intersect) from the transcriptomic information
norm_counts_lcts = n(posa, 1);
%Arrange an array with the sorted normalized gene counts and their corresponding
%index from the transcriptomic information
[sorted, ind_sorted] = sort(norm_counts_lcts);
%Set the threshold in the quantile given by prob
threshold1 = quantile(norm_counts_lcts, prob);
%Make a list with the normalized gene counts from the transcriptomic data
%for this condition
expression = n(:,1);
%Make a vector for the genes in the model and filled it with -1 
data = -1*ones(size(model1.genes));

for i = 1:length(model1.genes)
   if ismember(model1.genes(i), features)
        data(i) = expression((ismember(features, model1.genes{i})));
   end
end
expressionDataLAC = struct('gene', {features}, 'value', expression);
[expressionRxns,parsedGPR] = mapExpressionToReactions(model, expressionDataLAC);
options_GIMME = struct('solver', 'GIMME', 'threshold', threshold1, 'expressionRxns', expressionRxns);

tissueModel_GIMME_LCTS = createTissueSpecificModel(model1, options_GIMME);
fba_GIMME_LCTS = optimizeCbModel(tissueModel_GIMME_LCTS);
yields(1,3)=fba_GIMME_LCTS.f/(rS_exp(1)/1e3)/PM_Csource(1);
maxMu(1,3)=fba_GIMME_LCTS.f;

save('GIMME_iLR578_LCTS.mat','tissueModel_GIMME_LCTS','expressionRxns','fba_GIMME_LCTS')

rxnID=tissueModel_GIMME_LCTS.rxns;
fluxValues=fba_GIMME_LCTS.x;
LCTS=table(rxnID,fluxValues);
writetable(LCTS,'LCTS.csv','Delimiter',',');

%% LNNT
model2=changeRxnBounds(model,'EX_Lacto_N_Neotetraose(e)',-rS(2),'l');
maxLNnT=optimizeCbModel(model2,'max');
minLNnT=optimizeCbModel(model2,'min');

norm_counts_LNNT = n(posa, 5);
% [sorted, ind_sorted] = sort(norm_counts_LNNT);

threshold3 = quantile(norm_counts_LNNT, prob);

expression = n(:,5);
data = -1*ones(size(model2.genes));

for i = 1:length(model2.genes)
   if ismember(model2.genes(i), features)
        data(i) = expression((ismember(features, model2.genes{i})));
   end
end
expressionDataLNNT = struct('gene', {features}, 'value', expression);
[expressionRxns, parsedGPR] = mapExpressionToReactions(model, expressionDataLNNT);
options_GIMME = struct('solver', 'GIMME', 'threshold', threshold3, 'expressionRxns', expressionRxns);

tissueModel_GIMME_LNNT = createTissueSpecificModel(model2, options_GIMME);
fba_GIMME_LNNT = optimizeCbModel(tissueModel_GIMME_LNNT);
yields(2,3)=fba_GIMME_LNNT.f/(rS_exp(2)/1e3)/PM_Csource(2);
maxMu(2,3)=fba_GIMME_LNNT.f;

save('GIMME_iLR578_LNNT.mat','tissueModel_GIMME_LNNT','expressionRxns','fba_GIMME_LNNT')

rxnID=tissueModel_GIMME_LNNT.rxns;
fluxValues=fba_GIMME_LNNT.x;
LNNT=table(rxnID,fluxValues);
writetable(LNNT,'LNNT.csv','Delimiter',',');


%% 3FL
model3=changeRxnBounds(model,'EX_3fucosyllactose(e)',-rS(3),'l');
max3FL=optimizeCbModel(model3,'max');
min3FL=optimizeCbModel(model3,'min');

norm_counts_3FL = n(posa, 6);
% [sorted, ind_sorted] = sort(norm_counts_3FL);

threshold5 = quantile(norm_counts_3FL, prob);

expression = n(:,6);
data = -1*ones(size(model3.genes));

for i = 1:length(model3.genes)
   if ismember(model3.genes(i), features)
        data(i) = expression((ismember(features, model3.genes{i})));
   end
end
expressionData3FL = struct('gene', {features}, 'value', expression);
[expressionRxns, parsedGPR] = mapExpressionToReactions(model, expressionData3FL);
options_GIMME = struct('solver', 'GIMME', 'threshold', threshold5, 'expressionRxns', expressionRxns);

tissueModel_GIMME_3fl = createTissueSpecificModel(model3, options_GIMME);
fba_GIMME_3fl = optimizeCbModel(tissueModel_GIMME_3fl);
yields(3,3)=fba_GIMME_3fl.f/(rS_exp(3)/1e3)/PM_Csource(3);
maxMu(3,3)=fba_GIMME_3fl.f;

save('GIMME_iLR578_3FL.mat','tissueModel_GIMME_3fl','expressionRxns','fba_GIMME_3fl')

rxnID=tissueModel_GIMME_3fl.rxns;
fluxValues=fba_GIMME_3fl.x;
FL3=table(rxnID,fluxValues);
writetable(FL3,'3FL.csv','Delimiter',',');

%% 6SL
model4=changeRxnBounds(model,'EX_6sialyllactose(e)',-rS(4),'l');
max6SL=optimizeCbModel(model4,'max');
min6SL=optimizeCbModel(model4,'min');

norm_counts_6SL = n(posa, 7);
% [sorted, ind_sorted] = sort(norm_counts_6SL);

threshold6 = quantile(norm_counts_6SL, prob);

expression = n(:,7);
data = -1*ones(size(model4.genes));

for i = 1:length(model4.genes)
   if ismember(model4.genes(i), features)
        data(i) = expression(find(ismember(features, model4.genes{i})));
   end
end
expressionData6SL = struct('gene', {features}, 'value', expression);
[expressionRxns, parsedGPR] = mapExpressionToReactions(model, expressionData6SL);
options_GIMME = struct('solver', 'GIMME', 'threshold', threshold6, 'expressionRxns', expressionRxns);

tissueModel_GIMME_6sl = createTissueSpecificModel(model4, options_GIMME);
fba_GIMME_6sl = optimizeCbModel(tissueModel_GIMME_6sl);
yields(4,3)=fba_GIMME_6sl.f/(rS_exp(4)/1e3)/PM_Csource(4);
maxMu(4,3)=fba_GIMME_6sl.f;

save('GIMME_iLR578_6SL.mat','tissueModel_GIMME_6sl','expressionRxns','fba_GIMME_6sl')

rxnID=tissueModel_GIMME_6sl.rxns;
fluxValues=fba_GIMME_6sl.x;
SL6=table(rxnID,fluxValues);
writetable(SL6,'6SL.csv','Delimiter',',');

%% Display final results
% Yields
C_source={'Lactose';'LNNT';'3fl';'6sl'};
Exp_yield=yields(:,1)*1e3;
Sim_yield_no_expression=yields(:,2)*1e3;
Sim_yield_with_expression=yields(:,3)*1e3;
Table1=table(C_source,Exp_yield,Sim_yield_no_expression,Sim_yield_with_expression)
Table2

% Max. specific growth rates
Exp_mu=maxMu(:,1);
Sim_mu_no_expression=maxMu(:,2);
Sim_mu_with_expression=maxMu(:,3);
Table_mu=table(C_source,Exp_mu,Sim_mu_no_expression,Sim_mu_with_expression)

%% Flux variability analysis
% Rename variables for simple presentation
min_yield_Lactose=min_yields(:,1);
min_yield_LNnT=min_yields(:,2);
min_yield_3FL=min_yields(:,3);
min_yield_6SL=min_yields(:,4);
Table_min_yields=table(Met_names,min_yield_Lactose,min_yield_LNnT,min_yield_3FL,min_yield_6SL)

max_yield_Lactose=max_yields(:,1);
max_yield_LNnT=max_yields(:,2);
max_yield_3FL=max_yields(:,3);
max_yield_6SL=max_yields(:,4);
Table_max_yields=table(Met_names,max_yield_Lactose,max_yield_LNnT,max_yield_3FL,max_yield_6SL)
