%% 




parentpath = pwd;
addpath([parentpath '/eigenspace_analysis/']);
addpath([parentpath '/regularization/']);
addpath([parentpath '/plotsFn/']);
addpath([parentpath '/uncontraint_LS/']);
addpath([parentpath '/constraint_LS/']);

%% save data to your local DIR:  
% making DIR to your root DIR
if ispc
    SAVE_DIR = [getenv('USERPROFILE'), '\DataAnalysis\dartr_Fredholm\output\'];
else
    SAVE_DIR = [getenv('HOME'),'/DataAnalysis/dartr_Fredholm/output/'];     
end
if ~exist(SAVE_DIR,'dir'), mkdir(SAVE_DIR); end    
addpath(SAVE_DIR);


