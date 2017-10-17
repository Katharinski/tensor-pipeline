% function to match features and templates by computing the confusion
% matrix between each pair and the resulting kappa index
% input: features #brain regions x #features
% output: matrix of kappa indices 

% add additional input argument "option" and set to 2 in order to group
% some templates and reduce their number 

function kappa_matrix = match_templates(features,templates,levels)
%load('rsn_weights','rsn_weights')
%load('Human_66','Order') 
%template = rsn_weights(:,Order)'; % results in nregions x ntemplatesindex = 0:levels-1;
index = 0:levels-1;
% quantize templates
qtemp = zeros(size(templates));
for i=1:size(templates,2)
    limits = linspace(min(templates(:,i)),max(templates(:,i)),levels+1);   
    partition = limits(2:end-1);
    qtemp(:,i) = quantiz(templates(:,i),partition);
end
% quantize features
qfeat = zeros(size(features));
for i=1:size(features,2)
    limits = linspace(min(features(:,i)),max(features(:,i)),levels+1);
    partition = limits(2:end-1);
    qfeat(:,i) = quantiz(features(:,i),partition);
end
% match each feature with all templates
kappa_matrix = zeros(size(templates,2),size(features,2));
for i=1:size(features,2)
    for j=1:size(templates,2)
        % build confusion matrix (see Wikipedia entry)
        conf_matrix = zeros(levels);
        for k1=index
            for k2=index
                conf_matrix(k1+1,k2+1) = nnz((qfeat(:,i)==k1) & (qtemp(:,j)==k2));
            end
        end
        totsum = sum(conf_matrix(:));
        Pr_a = trace(conf_matrix)/totsum;
        Pr_e =  sum((sum(conf_matrix)/totsum).*(sum(conf_matrix,2)/totsum)');
        kappa_matrix(j,i) = (Pr_a-Pr_e)/(1-Pr_e);
    end
end
end
