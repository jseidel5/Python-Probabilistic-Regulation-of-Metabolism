import cobra
import libsbml
import pandas as pd
from scipy import stats
import numpy as np

def read_data(model_path,expression_path,regnet_path,regnet_sheet = 0,expr_sheet = 0):
    ## reads metabolic model
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
    
def binarize_expression(norm_expression):
    pass

#def calc_PROM_probability(reg_expr,targ_expr):
#    pass

def PROM_KO(reg2KO,model):
    pass

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
        #self.in_model = [] ## if Regulator is in model get -> KO
        
    def get_targets(self,regnet):
        reduced_df = regnet.loc[regnet['regulator'] == self.name] 
        self.targets = reduced_df['target'].to_list()
        return self.targets

    def make_target_objects(self):
        target_objects = {}
        for target in self.targets:
            target_objects[target] = Target(target,self.name)
        return target_objects
    
    def get_norm_expression(self,norm_expression):
        # could probably be nicer; creates nested list [[]]???
        nested_expression = norm_expression.loc[norm_expression['id'] == self.name].values.tolist()
        self.norm_expression = nested_expression[0][1:] # remove fist element -> id
        return self.norm_expression
    
    def get_bin_expression(self,bin_expression):
        pass
    
    def test_in_model(self,model):
        try:
            model.genes.get_by_id(self.name)
        except:
            self.in_model = False
        else:
            self.in_model = True
    
class Target(Regulator):
    def __init__(self,name,regulator):
        self.regulator = regulator
        self.name = name
        self.reactions = []
        self.norm_expression = []
        self.reactions = []
        #self.in_model = None
        #self.is_significant
        #self.interaction_prob
    
    def get_reactions(self,model):
        gene = model.genes.get_by_id(self.name)
        self.reactions = [rxn.id for rxn in gene.reactions]
    
    def calc_interaction_prob():
        pass
    
def PROM_wrapper(reg2KO,model,raw_expression,regnet,binarize_thresh = 0.33):
    regulator = Regulator(name = reg2KO)
    regulator.get_targets(regnet)
        #print('re_norm_expr',regulator.norm_expression)
        #print(regulator.targets)
    target_objects = regulator.make_target_objects()
        #print(raw_expression.iloc[:,1:])
    norm_expression = impute_and_quantile_expression(raw_expression) 
        #print('norm_expr_2',norm_expression.loc[norm_expression['id'] =='Rv0117'])
    regulator.get_norm_expression(norm_expression)
        #print(regulator.norm_expression)
    regulator.test_in_model(model)
    #print(regulator.in_model)
    for k,v in target_objects.items():
        v.test_in_model(model)
        v.get_norm_expression(norm_expression)
        v.get_reactions(model)
        print(v.name,'regulated by:',v.regulator)
        print(v.in_model)
        #print(v.norm_expression)
        #print(v.reactions)
        
 
 
 
model,raw_expression,regnet = read_data(
    expression_path=r'C:\Users\seidel\Documents\Integration_RMN_Bachelor\results\01.07.2020 PROM 3.0\promdata\expressiondata_iNJ661.xlsx',
    expr_sheet=0,
    model_path=r'C:\Users\seidel\Documents\Integration_RMN_Bachelor\results\01.07.2020 PROM 3.0\promdata\iNJ661.xml',
    regnet_path=r'C:\Users\seidel\Documents\Integration_RMN_Bachelor\results\01.07.2020 PROM 3.0\promdata\regulator_target_iNJ661.xlsx'
)       
        
        
PROM_wrapper('Rv0117',model,raw_expression,regnet)
