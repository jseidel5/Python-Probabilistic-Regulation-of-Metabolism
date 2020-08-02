import cobra
import libsbml
import pandas as pd
from scipy import stats
import numpy as np
from cobra.flux_analysis import flux_variability_analysis
from optlang.symbolics import add
import os.path
import re


def read_data(model_path,expression_path,regnet_path,regnet_sheet = 0,expr_sheet = 0):
    ## reads metabolic model/network
    model = cobra.io.read_sbml_model(model_path)
    ## each column in excel file is expression for 1 gene; gene name is headline
    raw_expression = pd.read_excel(expression_path,sheet_name = expr_sheet) # sheet = 1 -> 2nd sheet
    ## regnet is excel file with 2 columns with headers 'regulator' and 'target'
    regnet = pd.read_excel(regnet_path,sheet_name = regnet_sheet)
    return model,raw_expression,regnet
    
def impute_and_quantile_expression(raw_expression):
    #TODO need to impute
    raw_expression_copy = raw_expression.copy()
    norm_expression = quantileNormalize(raw_expression.iloc[:,1:])
    raw_expression_copy.update(norm_expression) 
    return raw_expression_copy
    
def binarize_expression(norm_expression,binarize_thresh,path):
    path2 = '.\\' + path + '\\bin_expr.pkl'
    if os.path.isfile(path2) == False:
        norm_expression_copy = norm_expression.iloc[:,1:].copy()
        expression_ids = norm_expression.copy()
        quantile = 0
        quantile_list = []
        for column in norm_expression_copy:
            quantile = np.quantile(norm_expression_copy[column],q = binarize_thresh)
            quantile_list.append(quantile)
            norm_expression_copy[column] = norm_expression_copy[column].map(lambda x: 1 if x >= quantile else 0 )
        expression_ids.update(norm_expression_copy)
        expression_ids.to_pickle(path2)
    else:
        expression_ids = pd.read_pickle(path2)
    return expression_ids

#########################

## taken from https://github.com/ShawnLYU/Quantile_Normalize

def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df

###########################

class Regulator:
    def __init__(self,name):
        self.name = name
        self.targets = []
        self.norm_expression = []
        self.bin_expression = [] ## binarized expression
        self.rxn_constrain_prob = {}
        #self.in_model = [] ## if Regulator is in model get -> KO
        
    def get_targets(self,regnet):
        reduced_df = regnet.loc[regnet['regulator'] == self.name] 
        self.targets = reduced_df['target'].to_list()
        return self.targets

    def make_target_objects(self):
        target_objects = {}
        target_ref = []
        for target in self.targets:
            target_objects[target] = Target(target,self)
        return target_objects
    
    def get_norm_expression(self,norm_expression):
        # could probably be nicer; creates nested list [[]]???
        nested_expression = norm_expression.loc[norm_expression['id'] == self.name].values.tolist()
        ## debug
        #print('debug:',nested_expression)
        ##
        self.norm_expression = nested_expression[0][1:] # remove fist element -> id
    
    def get_bin_expression(self,bin_expression):
        #same as get normexpression
        nested_expression = bin_expression.loc[bin_expression['id'] == self.name].values.tolist()
        self.bin_expression = nested_expression[0][1:]
    
    def test_in_model(self,model):
        try:
            model.genes.get_by_id(self.name)
        except:
            self.in_model = False
        else:
            self.in_model = True
    
    def get_rxn_constrain_prob(self,target_objects):
        for k,v in target_objects.items(): # loop through genes associated with regulator
            if v.is_significant == True:
                for rxn in v.reactions: # loop through rxn associated with gene
                    if rxn.id in self.rxn_constrain_prob: # if its already in dic only change if interaction_prob is lower
                        if v.interaction_prob < self.rxn_constrain_prob[rxn.id]:
                            self.rxn_constrain_prob[rxn.id] = v.interaction_prob
                    else:
                        self.rxn_constrain_prob[rxn.id] = v.interaction_prob
    
    def knockout(self,model,biomass_rxn,path):
        # without 'with:' so exchanges can be set outside scope; before running PROM_wrapper()
        standard_fba = model.optimize()
        obj_rxn = model.reactions.get_by_id(biomass_rxn)
        path2 = '.\\' + path + '\\fva.pkl'
        if os.path.isfile(path2) == False:
            #print('Analysing Flux Variability')
            fva_result = flux_variability_analysis(model,loopless=True)
            fva_result.applymap(lambda x:x if abs(x) >= 0.00001 else 0)
            fva_result.to_pickle(path2)
        else:
            #print('Excel file with FVA results found')
            fva_result = pd.read_pickle(path2)
        with model:
            objective_variables = []
            for rxn_id,rxn_prob in self.rxn_constrain_prob.items():
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.fva_min = fva_result.loc[rxn.id,'minimum']
                rxn.fva_max = fva_result.loc[rxn.id,'maximum']
                    
                ## only constrain rxn if associated genes are constrained 
                ## change gpr = rxn.... and regex if GPR is saved somewhere else and or genes have other name format
                
                #gpr = rxn._gene_reaction_rule # for iSynCJ816
                gpr = rxn.gene_reaction_rule # for iNJ661
                
                for target in self.targets:
                    if target in gpr:
                        gpr = gpr.replace(target,'True')
                gpr = re.sub(r'Rv[a-zA-Z0-9]*','False',gpr) # for iNJ661
                #gpr = re.sub(r'SGL_RS[a-zA-Z0-9]*','False',gpr)  # for iSynCJ816 
                #print('here:',gpr)
                
                ## Penalty Variable Weights
                ## ALPHA
                if rxn.flux < -10e-6:
                    rxn.alpha = model.problem.Variable(name = rxn.id + '_' + 'alpha',lb = 0, type = 'continuous')
                    weight = max(abs(rxn.fva_min),0.001)
                    rxn.alpha_weight = min(
                        -1*abs(standard_fba.objective_value)/abs(weight),
                        -1*abs(standard_fba.objective_value)/abs(rxn.fva_min)
                    )
                    objective_variables.append(rxn.alpha * rxn.alpha_weight)
                    lower_bound = min(rxn.fva_min * rxn_prob,-10e-6)
                    alpha_constraint = model.problem.Constraint(
                        rxn.flux_expression - rxn.alpha,
                        lb = lower_bound,
                        name = rxn.id + '_' + 'alpha_cons'
                    )
                    model.add_cons_vars(alpha_constraint)
                ## BETA
                elif rxn.flux > 10e-6:
                    rxn.beta = model.problem.Variable(name = rxn.id + '_' + 'beta',lb = 0, type = 'continuous')
                    weight = max(abs(rxn.fva_max),0.001)
                    rxn.beta_weight = min(-1*abs(standard_fba.objective_value)/abs(weight),0)
                    objective_variables.append(rxn.beta * rxn.beta_weight)    
                    upper_bound = max(rxn.fva_max * rxn_prob,10e-6)
                    beta_constraint = model.problem.Constraint(
                        rxn.flux_expression - rxn.beta,
                        ub = upper_bound,
                        name = rxn.id + '_' + 'beta_cons'
                    )
                    model.add_cons_vars(beta_constraint)
            ##debug 
            #print('obj_vars',objective_variables)
            ## Objective
            prom_objective = model.problem.Objective(
                obj_rxn.flux_expression + add(objective_variables) ,
                direction = 'max'
            )  
            model.objective = prom_objective
            #print('prom_objective',model.objective)
            ko_solution = model.optimize()
            #for rxn_id,rxn_prob in self.rxn_constrain_prob.items():
            #    print('{} flux:'.format(rxn_id),standard_fba[rxn_id])
        return ko_solution,objective_variables,standard_fba
                
class Target(Regulator):
    def __init__(self,name,regulator):
        self.regulator = regulator
        self.name = name
        self.reactions = []
        self.norm_expression = []
        self.bin_expression = []
        self.reactions = []
        
        #self.in_model 
        #self.is_significant
        #self.interaction_prob
    
    def get_reactions(self,model):
        gene = model.genes.get_by_id(self.name)
        self.reactions = [rxn for rxn in gene.reactions]
    
    def do_ks_test(self):
        targ_where_reg_on = []   # normalized expr values for target where binarized expression of regulator is 1/0
        targ_where_reg_off = []
        for index,expr in enumerate(self.regulator.bin_expression):
            if expr == 1:
                targ_where_reg_on.append(self.norm_expression[index])
            else:
                targ_where_reg_off.append(self.norm_expression[index])
        if stats.ks_2samp(targ_where_reg_on,targ_where_reg_off).pvalue <= 0.05:
            self.is_significant = True
        else:
            self.is_significant = False
    
    def calc_interaction_prob(self):
        if self.is_significant == False:
            self.interaction_prob = 1
        else:
            k_bin_targ_on_where_reg_off = 0
            k_bin_reg_off = 0
            for index,expr in enumerate(self.regulator.bin_expression):
                if expr == 0:
                    k_bin_reg_off += 1
                    if self.bin_expression[index] == 1:
                        k_bin_targ_on_where_reg_off += 1
            self.interaction_prob = k_bin_targ_on_where_reg_off/k_bin_reg_off
        
def PROM_wrapper(reg2KO,orig_model,raw_expression,regnet,biomass_rxn,detail_print=False,path,binarize_thresh = 0.33):
    model = orig_model.copy()
    regulator = Regulator(name = reg2KO)
    regulator.get_targets(regnet)
    target_objects = regulator.make_target_objects()
    #norm_expression = impute_and_quantile_expression(raw_expression) 
    norm_expression = raw_expression
    regulator.get_norm_expression(norm_expression)
    regulator.test_in_model(model)
    #print('raw expr:',raw_expression)
    #print('norm expr:',norm_expression)
    bin_expression = binarize_expression(norm_expression,binarize_thresh,path)
    #print(quantile_list)
    #print('bin expr:',bin_expression)
    regulator.get_bin_expression(bin_expression)
    
    #print('reg norm expr:',regulator.norm_expression)
    #print('reg bin expr:',regulator.bin_expression)
    
    for k,v in target_objects.items():
        v.test_in_model(model)
        v.get_norm_expression(norm_expression)
        v.get_reactions(model)
        v.get_bin_expression(bin_expression)
        v.do_ks_test()
        v.calc_interaction_prob()
        if detail_print == True:
            print(v.name,'regulated by:',v.regulator.name)
            print(' |-in model:',v.in_model)
            print(' |-significant interaction:',v.is_significant)
            print(' |-PROM prob:',v.interaction_prob)
            #print('norm expr:',v.norm_expression)
            #print('bin expr:',v.bin_expression)
            print(' |-reactions:',[rxn.id for rxn in v.reactions])
            print(' |___________________________________')
    regulator.get_rxn_constrain_prob(target_objects)
    prom_solution,penalty_variables,standard_fba = regulator.knockout(model,biomass_rxn,path)
    if detail_print == True:
        print('rxn_constrain_probs:',regulator.rxn_constrain_prob)
        print('penalty variables:',penalty_variables)
    print('------------------------------------')
    print('Original objective value:',standard_fba.objective_value)
    print('PROM objective values',prom_solution.objective_value)
    print('------------------------------------')
