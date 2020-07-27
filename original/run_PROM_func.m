function [f,v,f0,v0,status1,lostxns,model_parameters] =  PROM(model,expression,expressionid,regulator,targets,litevidence,prob_prior)
%% [f,v,status1,lostxns] = promv2(model,expression,expressionid,regulator,targets,litevidence,prob_prior);
%% [f,v,f0,v0,f_ko,v_ko,status1,lostxns,weights11_all,A_all,dxdt0_all,lb11_all,ub11_all,ctype1_all] =  run_PROM_func(model,expression,expressionid,regulator,targets,litevidence,prob_prior)
% The PROM algorithm predicts the growth phenotype and the flux response
% after transcriptional perturbation, given a metabolic and regulatory
% network.
% INPUTS
% Model is the metabolic model obtained from COBRA toolbox through readcbmodel command
%
% Gene expression data - rows - genes,columns - conditions; no need to normalize or impute
%
% Expressionid - an array of identifiers for each row/gene should be included
%h
% regulatory network - format - cell array of regulators and matching target genes
% example   Regulator = {'RegA'; 'RegB' ; 'RegC'};  Targets =
% {'GeneA';'GeneB';'GeneC'}
% note that the names or identifiers used in the regulatory data should
% match the names/ids given for gene expression data
%
%OPTIONAL - litevidence & prob_prior -> these should have the same length
%as the regulator/target arrays;
% high confidence interactions (not necessarily based on literature) should be flagged as one in litevidence
%array and the rest should be set as 0.
% Prob_prior should be set values between 0 and 1 for those interactions with litevidence. ( other values in
%the array would be ignored)
%
% OUTPUT - the algorithm gives the growth rate (f) and flux response (v) after knock out of all
% the regulators in the regulatory model; status is the glpk solver status 
% the status should be 5 for glpk; if its not then check solver error log
% lostxns gives the interactions that could not be quantified based on the
% threshold set for binarization
% the program would shoot a warning if the threshold chosen is bad. The
% default value (0.2 - 0.4) should work for most cases 
%
% EXAMPLE
% load mtbpromdata
% [f_ko,v_ko] = prom(model,expression,expressionid,regulator,targets);
% this is the tuberculosis model used in the PNAS paper

% DETAILED OUTPUT (Midani)
% f0               model cost based on wild-type model 
% f                model cost based on penalized-PROM model
% f_ko             model cost based on non-penalized-PROM model 
% v0               fluxes based on wild-type model
% v                fluxes based on penalized-PROM model
% v_ko             fluxes based on non-penalized-PROM model  
% status1          linear programming solver status
% lostxns          interactions that could not be quantified based on binarization threshold
% weights11_all    objective function updated after each TF knockout
% A_all            stoichiometry updated after each TF knockout (should not change)
% dxdt0_all        FBA RHS updated after each TF knckout
% lb11_all         FBA lower bounds (for penalized modle) after each TF knockout
% ub11_all         FBA lower bounds (for penalized modle) after each TF knockout
% ctype1_all       FBA sense of each constraint (i.e. minimize vs maximize); should be maximize for growth, and minimize for penalties. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 5, litevidence = [];prob_prior = [];
elseif (nargin < 5) || (nargin == 6),
    error('Incorrect number of input arguments to %s',mfilename)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('initializing data')
regulated = targets; % targets of TrmB
[tfnames,b,m]= unique(regulator);  % names of the unique targets; b maps the first call of each unique regulator, while m maps each unique regulator to each of its calls. 
weights      = model.c; 
stoic        = model.S; 
S            = model.S; 
ctype        = repmat('S',size(model.b));
lbff         = model.lb; 
ubff         = model.ub;
dxdt         = model.b; 
param.msglev = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finding rxn/gene position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[u,v] = find(model.rxnGeneMat); 
% u - reaction; v - genes
rxnpos = u;    % size equal to v
genelist = v;
clear u v

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% need to find the reactions that each gene has influence on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%d%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnumsinexpsn = expressionid; %genes with expression data
litevidence = logical(litevidence);   % genes with literature evidence 
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%d%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Baseline Wild-type FBA of model: Must match _fba_only results in _fba_only/_main_output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%d%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scou = 1;
lbg  = model.lb; 
ubg  = model.ub;
a1   = [S, zeros(size(S,1),length(ubg)) , zeros(size(S,1),length(ubg)) ]; %size 643x2340
a2   = sparse([eye(length(ubg)), eye(length(ubg)),zeros(length(ubg))]);   %size 780x2340
a3   = sparse([eye(length(ubg)), zeros(length(ubg)),-eye(length(ubg))]);  %size 780x2340
A    = [a1;a2;a3]; %size 2203x2340 
weights11 = [weights;zeros(2*length(lbg),1)]; %size 2340 (optimal,lb,ub)
weights00 = [weights;zeros(2*length(lbg),1)]; %size 2340 (optimal,lb,ub)
lb11      = [-1e8*ones(length(lbg),1);zeros(length(lbg),1);zeros(length(lbg),1)]; %size 2340 lb of optimal 
ub11      = [1e8*ones(length(lbg),1);zeros(length(lbg),1);zeros(length(lbg),1)];  %size 2340 ub of optimal
dxdt0     = [zeros(size(S,1),1);lbg;ubg]; %size 2203 RHS: Sv=dxdt0 with original model.lb and model.ub as RHS elements
ctype1    = [repmat('S',size(S,1),1);repmat('L',size(lbg,1),1);repmat('U',size(lbg,1),1)];

[v0,f0] = glpk(-weights11,A,dxdt0,lb11,ub11,ctype1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Probabilities using a Global Threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% impute and quantile expression data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
lost_xn = false(size(regulated));  %vector of zeros

disp('finding probabilities')
cou = 1;cou1 = 1;cou3 = 1;
data = expression;                  %expression data e.g. 2400x361
data = knnimpute(data);             %expression data imputed for missing values
data = quantilenorm(data);          %expresssion data quantile normalized
data1 = data;                       %imputed-normalized data copy 
datathresh = quantile(data(:),0.33);%id threshold for al 2400 genes using 0.33th-percentile. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% binarize expression data based on thereshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if datathresh < 0,
    data(data>=datathresh) = 1;
    data(data < datathresh) = 0;
else
    data(data < datathresh) = 0;
    data(data>=datathresh) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% modulate gene expression based on regulator status (on/off)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for  i = 1:length(regulated) % for each target 
    k = find(ismember(bnumsinexpsn,regulated(i)));   %is epxression data available for target
    l = find(ismember(bnumsinexpsn,regulator(i)));   %is expression data available for its regulator (trmB)
    if ~isempty(k) & ~isempty(l)    %if both are available
        te = data1(k,:);            %te is expression vector of target
        te1 = data1(l,:);           %te1 is expression vector of regulator
        
        tec = data(k,:);            %tec is binarized expression vector of target
        tec1 = data(l,:);           %tec1 is binarized expression vector of regulator 
        
        cou1 = cou1 + 1;
        %do they (expression vector target when regulator is on AND 
        %expression vector of target when regulator is off)-l
        %follow same distirbution?
        try kstest2(te(tec1 == 1),te(tec1== 0));   
            
            if  (kstest2(te(tec1 == 1),te(tec1== 0)) == 1), % if different distribution, then 
                % prob is count of expression for target when regulator is
                % off over total count of times regulator is off
                prob1 = sum(tec(tec1 == 0))/length(tec(tec1==0)); 
                probtfgene(i) = prob1; %save probability of all targets
                cou = cou + 1;
                %   this formula also gives the same answer  - (sum(~tec1(tec == 1))/length(tec1(tec==1))) * (sum(tec)/sum(~tec1))
                
            else
                
                probtfgene(i) = 1;  % no effect
                
            end
            
        catch ERRLG    % cant be estimated from microarray ; if it has strong evidence, consider setting this to zero 
            probtfgene(i) = 1;
            lost_xn(i) = 1;            
            
        end
    else
        probtfgene(i) = 1; %If no evidence exists for effect of regulator on target, assume target is always on
                    lost_xn(i) = 1;
    end
end
probtfgene = probtfgene(:);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check if too many interactions were unclear/missing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if there is problem with binarizartion.. usually there ll be some
% interactions that'd be missed, but if most of it (like >75% ) are...
if (sum(lost_xn) > 0.75*length(probtfgene))   
   % missed then its time to change the threshold
   datathreshflag = 1;
   disp('change binarization threshold')
else
    datathreshflag = 0;
end

if ~isempty(litevidence)
    % you could set those interactions that you think have strong
    % literature evidence to have predefined probabilities
    probtfgene(litevidence) = prob_prior(litevidence);  
end                                                             

toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Baseline Wild-type FBA of model: Should match _fba_only results in _fba_only/_main_output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lbf   = lbff; %original model.lb
ubf   = ubff; %original model.ub
[v,f] = glpk(-weights,stoic,dxdt,lbf,ubf,ctype); % use model.lb as constraints with zeros RHS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  run PROM for each regulator knockout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('running PROM')
count    = 1; 
thresh   = 10^(-6); 
mthresh  = 10^(-3);
allgenes = [model.genes;tfnames]; %size 488+69 (last 69 are repeats)
[ir,posgenelist] = ismember(regulated,model.genes); %[membership, index in model.genes for each target]
hw = waitbar(0);
kappa = 1;
vm    = zeros(size(v));
bnumstobekoed   = tfnames;   
% bnumstobekoed -  gives the geneids of the genes to be knocked out 
% by default it knocksout all the tfs in the model one by one
% bnumstobekoed = [tfs_that_are_off;tf_that u want to knockout];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  in our case, it is ony one TF (trmB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ci = 1:length(bnumstobekoed)
    
    lbg = lbf; ubg = ubf; % original model bounds
    lb11 = [-1e8*ones(length(lbg),1);zeros(length(lbg),1);zeros(length(lbg),1)]; %uniform model bounds
    ub11 = [1e8*ones(length(lbg),1);zeros(length(lbg),1);zeros(length(lbg),1)];  %uniform model bounds
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% check if it's a metabolic or regulatory gene or both
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if any(strcmpi(model.genes,bnumstobekoed(ci))) %is trmB in the model gene list
        %if so, find the reaction (numbers) associated with this TF (tmrB)
        temppos = rxnpos(genelist == find(strcmp(model.genes,bnumstobekoed(ci))));  
    	for jj = 1:length(temppos)      % for each reaction
            if model.rev(temppos(jj))   % if reaction is reversible
                lbg(temppos) = -thresh; % specify its bounds to be very cloes to zero 
                ubg(temppos) = thresh;
            else
                lbg(temppos) = -thresh; %why specify lb if rxn is irreversible (numerical issues?)
            end
        end
        
    end
    
    %solve for fluxes using original model bounds
    [v1,fk(ci)]  = glpk(-weights,S,dxdt,lbg,ubg,ctype); 
    
    % check if gene to be koed is a TF with targets in model.
    % gene can be knocked out but might not have any targets, 
    % hence this is a necessary check
    if any(ismember(tfnames,bnumstobekoed(ci))), 
        
              tfstate = logical(zeros(size(tfnames)));             % all tf state are initialized to false 
        tfstate(find(ismember(tfnames,bnumstobekoed(ci)))) = 1;    % modify only current TF to TRUE
                    k = find(ismember(regulator,tfnames(tfstate)));% find the TF involvement in regulator-target pairs
             tempgene = regulated(k);                              % find its targets
        tempgeneprobs = probtfgene(k);                             % find probabilty of expression of targets if this TF is off
          tempgenepos = posgenelist(k);                            % find index in model.gene for each target
           temprxnpos = rxnpos(ismember(genelist,tempgenepos));    % find reactions associated with these targets
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %% Constrain flux through the reactions associated with these genes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          x = true(size(model.genes));        % initialize all genes in model to true
        [isInModel,geneInd] = ismember(tempgene,model.genes); % for each target, find its index in model.genes
                 x(geneInd) = false;                          % then, knock it out by setting it to false
               constrainRxn = false(length(temprxnpos),1);    % initialize list of TF-associated reaction to false 
               
        % Figure out if any of the reaction states are changed
        for j = 1:length(temprxnpos)               % for each TF-associated reaction
            if (~eval(model.rules{temprxnpos(j)})) % if none of its associated genes are on (based on user-defined rules), then it will be constrained 
                constrainRxn(j) = true;
            end
        end
        
        tempgeneprobs(tempgenepos == 0)  = '';     %if a target is not member of the gene model, then it won't/can't be modified
          tempgenepos(tempgenepos == 0)  = '';
        % temprxnpos has the rxns that are going to  be affected by a TF
        % (i.e. TF-associated reactions)
        % krxnpos are the rxns that will be affected by a target gene
        % (i.e. gene-asscoiated reactions)
        % we loop around all the genes.. 
        
        
        for l = 1:length(tempgenepos)   % for each target(ed) gene
            if ~isnan(tempgeneprobs(l)) % if it has a numerical probability associated with regulator off status
                % find its associated reactions (i.e. krxpos)
                krxnpos = ismember(temprxnpos,rxnpos(ismember(genelist,tempgenepos(l))));
                for m = 1:length(temprxnpos) % for all TF-associated reactions
                    if krxnpos(m)            % if the reaction is associated with the current looping target gene 
                        tempsubsysmid = 0;   % midpoint is zero
                        
                        if constrainRxn(m)                   % if reaction is to be constrained
                            if (tempgeneprobs(l) < 1)        % if its 1 no use in changing the bounds - might as well save time
                                if (tempgeneprobs(l) ~= 0)   % if its zero no point in estimating vm again - saves time.. but cant include in the above statement because you have to change the bounds
                                    if ~vm(temprxnpos(m))    % done to save time - if estimated already use it
                                       
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
                                        %% Determine bounds for this TF-associated reaction using two FBA runs (one max LP, one min LP) 
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                                        weights1 = weights; 
                                        lbv = lbf;                     %original model.lb
                                        ubv = ubf;                     %original model.ub
                                        grwthpos = find(weights == 1); %find biomass reaction
                                        lbv(grwthpos) = v(grwthpos);   %always increase biomass reaction 
                                        weights1(temprxnpos(m)) = -1;  %minimize current constraining reaction
                                        [v11,fva1]  = glpk(-weights1,S,dxdt,lbv,ubv,ctype); %solve v11 and fva1
                                        weights1(temprxnpos(m)) = 1;   %maximize current constraining reaction
                                        [v12,fva2]  = glpk(-weights1,S,dxdt,lbv,ubv,ctype); %solve v12 and fva2
                                        
                                        if  v(temprxnpos(m)) < 0
                                            %if original flux is negative, minimze it further by selecting the flux from either
                                            %solutions (where you minimize/maximize reaction while maximizing biomass),
                                            %or original wildtype flux. 
                                            vm(temprxnpos(m)) = min([v11(temprxnpos(m)),v12(temprxnpos(m)),v(temprxnpos(m))]);
                                        elseif v(temprxnpos(m)) > 0
                                            %if original flux is positive, maximize it further by selecting the flux from either
                                            %solutions (where you minimize/maximize reaction while maximizing biomass),
                                            %or original wildtype flux. 
                                            vm(temprxnpos(m)) = max([v11(temprxnpos(m)),v12(temprxnpos(m)),v(temprxnpos(m))]);
                                        else
                                            %if original flux is zero, maximize it further by selecting the absolute flux from either
                                            %solutions (where you minimize/maximize reaction while maximizing biomass),
                                            %or original wildtype flux.  
                                            vm(temprxnpos(m)) = max([abs(v11(temprxnpos(m))),abs(v12(temprxnpos(m))),abs(v(temprxnpos(m)))]);
                                        end
                                        
                                        
                                    end
                                end
 
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
                                %% Update knockout model bounds and objective function with penalties for exceeding flux bounds
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                
                                % vm(temprxnpos(m)) is optimized flux through reaction "m" 
                                % so, xx = flux * (1-e(-tempsubsysmid) + flux * FluxProb * e(-tempsubsysmid)
                                xx = ( vm(temprxnpos(m))*(1 - exp(- tempsubsysmid) ) + ( vm(temprxnpos(m))*tempgeneprobs(l)*exp(-tempsubsysmid)));
                                    % if flux is negative
                                    %   lbf/ubf here is original model.lb/ub
                                    %   lbg/ubg here is initialized to be very close to zero, but updated when scrutinizing each reaction
                                    %   ub11 here is set to maximum possible bound 
                                if  v(temprxnpos(m)) < 0
                                    tem = max([lbf(temprxnpos(m)),xx,lbg(temprxnpos(m))]);  % make sure we arent violating the original bounds; also get the higher value if there were multiple modifications for the rxn
                                    lbg(temprxnpos(m)) = min([tem,-thresh]);                % prevents the solver from crashing (can't be too close to zero)
                                    ub11(1*length(ubg) + temprxnpos(m)) = 1000;             % modify RHS of the upper bound for the lower bound to be 1000 for the particular rxn. 
                                    weights11(1*length(ubg) + temprxnpos(m)) = ((-1)*kappa/abs(vm(temprxnpos(m))))*abs(f0);   % v0 f0 are the wild type values; kappa is set to 1 by default. 
                                    vv = max([abs(vm(temprxnpos(m))),mthresh]);
                                    % update objective function: add a penalty for exceeding a lower bound
                                    weights11(1*length(ubg) + temprxnpos(m)) = min([((-1)*kappa/abs(vv))*abs(f0), weights11(1*length(ubg) + temprxnpos(m)) ]);
                                elseif v(temprxnpos(m)) > 0
                                    % Recall: ubf is original; ubg is modified for knockouts as you loop for each reaction
                                    tem = min([xx,ubf(temprxnpos(m)),ubg(temprxnpos(m))]);  % choose tem as minimum of xx (which is optimal estimate) or upper bound of original or knockout model
                                    ubg(temprxnpos(m)) = max(tem,thresh);                   % choose whichever is higher the above or threshold; so if is negligible, move it up a bit. 
                                    ub11(2*length(ubg) + temprxnpos(m)) = 1000;             % modify RHS of upper bound of the upper bound to be 1000 for the particular rxn. 
                                    vv = max([abs(vm(temprxnpos(m))),mthresh]);             % choose absolute of optimal estimate or a predfined threshold; so if negligible, move it up a bit. 
                                    % update objective function: add a penalty for exceeding an upper bound
                                    weights11(2*length(ubg) + temprxnpos(m)) = min([((-1)*kappa/abs(vv))*abs(f0), weights11(2*length(ubg) + temprxnpos(m)) ]); 
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
    %% Run FBA on knockout-modified model   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    % fluxes are at steady state, hence dxdt0 are zero, but flux bounds (lbg, ubg) are not zero. 
    dxdt0 = [zeros(size(S,1),1);lbg;ubg]; %RHS constraints of FBA model
    % weights11 include your objective function constraint and constraints by the PROM knockout simulation
    [v00(ci,:),f00(ci),status1(ci)] = glpk(-weights11,A,dxdt0,lb11,ub11,ctype1);

    coun  = 1;        
    lbh   = lbg; % final lb of knockout model for current TF (ci)
    ubh   = ubg; % final lb of knockout model for current TF (ci)
    
    % track PROM FBA parameters for current TF (ci)
    weights11_all(:,ci) = weights11; 
    A_all{ci}           = A;
    dxdt0_all(:,ci)     = dxdt0;
    lb11_all(:,ci)      = lb11;
    ub11_all(:,ci)      = ub11;
    ctype1_all{ci}      = ctype1;
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
    %% Check for solver failures  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    while ((status1(ci) == 105) || (v00(ci,661) < 0)) % if solver failed, or objective function flux is negative 
        lbh(lbh ~= lbf) = lbh(lbh ~= lbf) - 1E-3;
        ubh(ubh ~= ubf) = ubh(ubh ~= ubf) + 1E-3;
        dxdt0 = [zeros(size(S,1),1);lbh;ubh];
        [v00(ci,:),f00(ci),status1(ci)] = glpk(-weights11,A,dxdt0,lb11,ub11,ctype1);
        coun = coun + 1;
        if (coun > 100),
            % if its impossible to estimate, then
            % check the unweighted g.r and the one with max weight prom -
            % if very less difference use it - warn the user about
            % the problem at the iteration number - ci;
            [v3,f3,status3] = glpk(-weights,S,dxdt,lbg,ubg,ctype); % original model
            [v30,f30,status30] = glpk(-weights00,A,dxdt0,lb11,ub11,ctype1); %original model without modifying weights for TF-associated reaction. 
            if abs((f3-f30)/abs(f3)) < 0.1
                f00(ci) = f3;
                
            else
                disp(' problem in'); disp(ci);break;  % if that doesnt work,  display a warning
            end
            
            disp('check'); disp(ci); break;
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
    %% Run FBA on non-penalized model   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    % knockout bounds
    lbg_st(ci,:) = lbg;
    ubg_st(ci,:) = ubg;
    
    %original bounds
    lb_st(ci,:) = lbff;
    ub_st(ci,:) = ubff;
    
    %original PROM but with knockout bounds
    [v2(ci,:),f1(ci),status] = glpk(-weights,S,dxdt,lbg,ubg,ctype);

    % close status bar
    ktime = toc;
    waitpar = [num2str(ceil(ci/length(tfnames)*100)),'% complete. Time taken:',num2str(ceil(ktime)),' secs'];
    waitbar(ci/length(tfnames),hw,waitpar);
    
    
    ff00(scou,ci) = v00(ci,661);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if datathreshflag
            if all(lost_xn(k))   % if none of the probabilities of a gene can be estimated, then ->
                v00(ci,:) = NaN;
                f1(ci) = NaN;
                v2(ci,:) = NaN;
                f00(ci) = NaN;
                %break;
            end
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear tempgenepos tempgeneprobs temprxnpos k
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
%% Final output values  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
%original PROM with knockout bounds
f_ko = -f1'; 
v_ko = v2;

%original PROM with penalized knockout bounds 
f00_ko(:,scou) = v00(:,661);
v00_ko = v00;

f = f00_ko; v = v00_ko;

lostxns(:,scou) = lost_xn;

model_parameters = {weights11_all,dxdt0_all,lb11_all,ub11_all,ctype1_all};
end
