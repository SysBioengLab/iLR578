%% Model
%Set solver, in this case, gurobi
changeCobraSolver('gurobi');

%Load model of Bifidobaterium infantis ATCC 15697
load('iLR578.mat')

%Add ratios with acetate
model = addRatioReaction(model, {'EX_ac(e)','EX_lac_L(e)'}, [2 3]);
% model = addRatioReaction(model, {'EX_ac(e)','EX_succ(e)'}, [1 100]);

% Define objective function for the analysis: 1) AGORA, 2) Schopping et al.
biomassEQ=1;                          
switch biomassEQ
    case 1
        bioRxn='Biomass_growth_AGORA';
    case 2
        bioRxn='Biomass_growth_Schopping';
end
model=changeObjective(model,bioRxn);

%% Single gene knockout
model1=changeRxnBounds(model,'EX_lcts(e)',-1,'l');
[grRatio1, grRateKO1, grRateWT1, hasEffect1, delRxns1, fluxSolution1] = singleGeneDeletion(model1, 'FBA');

model2=changeRxnBounds(model,'EX_Lacto_N_Neotetraose(e)',-0.46,'l');
[grRatio2, grRateKO2, grRateWT2, hasEffect2, delRxns2, fluxSolution2] = singleGeneDeletion(model2, 'FBA');

model3=changeRxnBounds(model,'EX_3fucosyllactose(e)',-0.67,'l');
[grRatio3, grRateKO3, grRateWT3, hasEffect3, delRxns3, fluxSolution3] = singleGeneDeletion(model3, 'FBA');

model4=changeRxnBounds(model,'EX_6sialyllactose(e)',-0.52,'l');
[grRatio4, grRateKO4, grRateWT4, hasEffect4, delRxns4, fluxSolution4] = singleGeneDeletion(model4, 'FBA');

save singleKO_results.mat 