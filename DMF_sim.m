% function to perform simulations with the dynamic mean field model 
% ref: Deco G, Ponce-Alvarez A, Hagmann P, Romani G, Mantini D, Corbetta M (2014) 
% How local excitation-inhibition ratio impacts the whole brain dynamics. 
% JNeurosci 34:7886–7898.
% code by Adrian Ponce, adpated by Katharina Glomb

function [all_TC,meanrates,stdrates] = DMF_sim(T_sim,SC,G,JI)

N = size(SC,1);

% simulation parameters for optimization
dtt  = 1e-3;  % Sampling rate of simulated neuronal activity (s) for BOLD
dt = 0.1; % time step (ms) for simulation

% some parameters different for final simulation with optimal params
tmaxJ = 10000; % simulation time for balance
tspanJ = length(0:dt:tmaxJ);
rounds_balance = 500;

taun = 100; % tau_exc (NMDA)
taug = 10; % tau_inh (GABA)
gamma = 0.641;
sig = 0.01;
JN = 0.15; % J_NMDA
J = ones(N,1);
I0 = 0.382;
Jexte = 1.; % W_E
Jexti = 0.7; % W_I

curr = zeros(tmaxJ,N);
delta = 0.02*ones(N,1);

% part that gets split
rng shuffle % make sure initial conditions are different
tmax = T_sim*2000 + 20000; % time for final sim.; fs=0.5Hz; BOLD() cuts off ~20s
tspan = length(0:dt:tmax);

% since code was written for Hagmann SC, matrix has to be rescaled
sumHagSC = 15.3014;
% CHag = load('Hagmann_matrix');
% sumHagSC = sum(CHag(:));
fac = sumHagSC/sum(SC(:));
SC = SC*fac;

% create arrays
Jfinal = zeros(N,length(G)); % record adjusted local feedback inhibition (J)
all_TC = zeros(T_sim,N,length(G)); % for recording TCs of each G
meanrates = zeros(N,length(G)); % record maximum firing rate
stdrates = zeros(N,length(G)); % record std of firing rates

index = 1; % count G-values (WE)
for g=G
    fprintf('G = %g\n',g);
    w = 1.4*ones(N,1);
    if exist('JI','var')
        J = JI(:,G==g);
    end
    for b=1:rounds_balance
        fprintf('Simulation for J-adjustment, round %i.\n',b);
        sn = 0.001*ones(N,1);
        sg = 0.001*ones(N,1);
        nn = 1;
        j = 0;
        for s=2:1:tspanJ
            xn = I0*Jexte+w*JN.*sn+g*JN*SC*sn-J.*sg;
            xg = I0*Jexti+JN*sn-sg;
            rn = phie(xn);
            rg = phii(xg);
            sn = sn+dt*(-sn/taun+(1-sn)*gamma.*rn./1000.)+sqrt(dt)*sig*randn(N,1);
            sn(sn>1) = 1;
            sn(sn<0) = 0;
            sg = sg+dt*(-sg/taug+rg./1000.)+sqrt(dt)*sig*randn(N,1);
            sg(sg>1) = 1;
            sg(sg<0) = 0;
            j = j+1;
            if j==10
                % use curr instead of the usual neuro_act for
                % reasons unknown to me (but correctness confirmed
                % by Gustavo)
                curr(nn,:) = xn'-125/310;
                nn = nn+1;
                j=0;
            end
        end
        
        currm = mean(curr,1);
        
        % adjust J separately for each region
        flag = 0;
        for n=1:1:N
            if abs(currm(n)+0.026)>0.005
                if currm(n)<-0.026
                    J(n) = J(n)-delta(n);
                    delta(n) = delta(n)-0.001;
                    if delta(n)<0.001
                        delta(n) = 0.001;
                    end
                else
                    J(n) = J(n)+delta(n);
                end
            else
                flag = flag+1;
            end
        end
        
        fprintf('%i out of %i regions are adjusted correctly.\n',flag,N);
        
        if flag==N
            break;
        end
    end
    
    
    Jfinal(:,index) = J;
    
    % Final Simulation
    fprintf('Starting final simulation. This will take a while...\n');
    neuro_act = zeros(tmax,N);
    rates = zeros(tmax,N);
    
    sn = 0.001*ones(N,1);
    sg = 0.001*ones(N,1);
    
    % Warm-Up (reach steady state):
    %---------
    
    for s=2:1:5000
        xn = I0*Jexte+w*JN.*sn+g*JN*SC*sn-J.*sg;
        xg = I0*Jexti+JN*sn-sg;
        rn = phie(xn);
        rg = phii(xg);
        sn = sn+dt*(-sn/taun+(1-sn)*gamma.*rn./1000.)+sqrt(dt)*sig*randn(N,1);
        sn(sn>1) = 1;
        sn(sn<0) = 0;
        sg = sg+dt*(-sg/taug+rg./1000.)+sqrt(dt)*sig*randn(N,1);
        sg(sg>1) = 1;
        sg(sg<0) = 0;
    end
    
    % Actual simulation
    %---------
    nn=1;
    j=0;
    
    for s = 2:1:tspan
        xn = I0*Jexte+w*JN.*sn+g*JN*SC*sn-J.*sg;
        xg = I0*Jexti+JN*sn-sg;
        rn = phie(xn);
        rg = phii(xg);
        sn = sn+dt*(-sn/taun+(1-sn)*gamma.*rn./1000.)+sqrt(dt)*sig*randn(N,1);
        sn(sn>1) = 1;
        sn(sn<0) = 0;
        sg = sg+dt*(-sg/taug+rg./1000.)+sqrt(dt)*sig*randn(N,1);
        sg(sg>1) = 1;
        sg(sg<0) = 0;
        j = j+1;
        if j==10
            neuro_act(nn,:) = xn';
            %neuro_act(nn,:) = rn';
            rates(nn,:) = rn';
            nn = nn+1;
            j = 0;
        end
    end
    
    nn = nn-1;
    meanrates(:,index) = mean(rates);
    stdrates(:,index) = std(rates);
    
    %%%% BOLD empirical
    fprintf('Simulation done, computing BOLD signal.\n');
    % Friston BALLOON MODEL
    T = nn*dtt; % Total time in seconds
    
    B = BOLD(T,neuro_act(1:nn,1)'); % B=BOLD activity, bf=Foutrier transform, f=frequency range)
    BOLD_act = zeros(length(B),N);
    BOLD_act(:,1) = B;
    
    for nnew=2:N
        %  nnew
        B = BOLD(T,neuro_act(1:nn,nnew));
        BOLD_act(:,nnew) = B;
    end
    
    bds = BOLD_act(1:2000:end,:);
    all_TC(:,:,index) = bds(1:T_sim,:); % for recording all TCs for all Gs
    index = index+1;
end
end
