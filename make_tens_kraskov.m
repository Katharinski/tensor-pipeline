% script for making tensors using Kraskov et al.'s MI(1) method
% do this for the 24 youngest of Petra's datasets, and also the simulated 
% data for SC and EC (for all G)

clearvars
w = 120; % window width for tensor, in s
% load correct time courses
% empirical: all subjects, task ID 49, simulated: all G's for one
% simulation, SC task IDs 1-24 (3 are missing), EC task IDs 25-48
%sge_task_id = str2num(getenv('SGE_TASK_ID'));
run = 49;%sge_task_id;
%run=25;
% runs that went wrong previously - all surrogates
if run<=48
    sim = 1;
else
    sim = 0;
end

if ~sim
    bds = load_petradata();
    [T,N,x] = size(bds); % emp data: x=#subj/sessions; sim data: x=#WE
    % global regression
    display('Warning! Global regression ON!\n')
    regr_signal = zeros(size(bds));
    for s=1:x
        regr_signal(:,:,s) = Global_regression(bds(:,:,s)')';
    end
    bds = regr_signal;
    %%%%%
    save_fname = 'tens_kraskov_petra_regr';
else 
    if (run<=24) 
        fname = ['sims_SC/gold10_S_',num2str(run),'.mat'];
        save_fname = ['tens_kraskov_SC',num2str(run)];
    elseif (run>24 && run<=48)
        fname = ['sims_EC/gold10_S_',num2str(run),'.mat'];
        save_fname = ['tens_kraskov_EC',num2str(run)];
    end
    if exist('fname','var')
        load(fname,'all_TC');
        bds = all_TC;
        clear all_TC
        [T,N,x] = size(bds); 
    end
end
if exist('save_fname','var')
    % real data
%     alltens = zeros(N,N,round((T-w/2)),x);
%     for i=1:x
%         fprintf('Real data, %i\n',i)
%         tens = MI_kraskov(bds(:,:,i),120);
%         alltens(:,:,:,i) = tens;
%     end
%     save(save_fname,'alltens');
    
    % surrogate data
    alltens_surr0 = zeros(N,N,round((T-w/2)),x);
    data_surr0 = zeros(size(bds));
    alltens_surr1 = zeros(N,N,round((T-w/2)),x);
    data_surr1 = zeros(size(bds));
    
    for s=1:x
        fprintf('Surrogate data, %i\n',s)
        surr_bds0 = surrogates_cov(bds(:,:,s),0);
        data_surr0(:,:,s) = surr_bds0;
        tens = MI_kraskov(surr_bds0,w);
        alltens_surr0(:,:,:,s) = tens;
        surr_bds1 = surrogates_cov(bds(:,:,s),1);
        data_surr1(:,:,s) = surr_bds1;
        tens = MI_kraskov(surr_bds1,w);
        alltens_surr1(:,:,:,s) = tens;
    end
    save([save_fname,'_surr'],'data_surr0','alltens_surr0','data_surr1',...
        'alltens_surr1','w');
end




