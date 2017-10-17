% example simulation script

% all simulations coded by Adrian Ponce-Alvarez (adapted by Katharina Glomb)
% ref: Deco G, Ponce-Alvarez A, Hagmann P, Romani G, Mantini D, Corbetta M (2014) 
% How local excitation-inhibition ratio impacts the whole brain dynamics. 
% JNeurosci 34:7886–7898.

clearvars
% load your SC/EC matrix
load('SC_matrix',SC)

% global coupling parameters; a reasonable range is 0.5 to about 5 or 6,
% but depends on your SC
G_array = 0.5:0.1:5;

% simulation time in TR=2s-frames; should match your data
T_sim = 661;
S = 24; % should match your number of subjects

% get values for feedback inhibition to initialize; save and run only once,
% will save time with simulations if you have to repeat them
J_fname = 'weight_file';
if ~exist(J_fname,'file')
    [JI,~,~] = Get_balanced_weights(SC,G);
else
    load(J_fname,JI)
end

% can be parallelized with parfor
for s=1:S
    fname_suffix = ['_' num2str(s)];
    [all_TC,meanrates,stdrates] = DMF_sim(T_sim,SC,G,JI);
    save(['simulation',fname_suffix],'all_TC','meanrates','stdrates','G')
end

