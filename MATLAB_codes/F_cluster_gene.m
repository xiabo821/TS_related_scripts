 function [cluster_exp,cluster_exp_UMI,gene_cluster_idx] = F_cluster_gene(input_matrix,cell_cluster,k,varargin)
%% [cluster_exp,cluster_exp_UMI,gene_cluster_idx] = F_cluster_gene(input_matrix,cell_cluster,k,varargin)
%  F_cluster_gene firstly merge the cells of same type and then cluster the
%  genes by kmeans clustering. 
%  
%  required functions:
%       plot_areaerrorbar.m
%       cbrewer.m
%  
%  Input: 
%       input matrix: gene(row)-by-cell(column)
%       cell_cluster: overall the cluster should be C1-Cn
%       k           : k for kmeans clustering 
%  varargin:
%       'include_unexpressed_gene': whether include unexpressed_gene.
%                           if 'yes', unexpressed_gene will be an aditional
%                           cluster 0(C0), generating k+1 clusters: C0 and C1-Cn.
%                           default is 'no';
%       'expression_cutoff_n': define expressed/unexpressed genes by minimumly
%                           expressed number of cells. By doing so, the 
%                           function firstly counts how many cell expressed
%                           this gene. If the number of expressed cells of 
%                           this gene is >= n, this gene is defined as expressed. 
%                           n should be integers, from 3 to inf.
%                           Default number is 3.
%       'umi_ratio_cutoff': define expressed/unexpressed genes by minimum
%                           expression level in the highest expressed cell 
%                           cluster. By doing so, the function first
%                           calculate the mean-UMI for each cluster, and
%                           then find out the max mean-UMI for each gene. If 
%                           the max mean-UMI is larger than the assigned ratio,
%                           it will be defined as expressed genes.
%                           Usually the ratio is assigned from 0.1 to 1, the cutoff
%                           indicating at least 1/10 UMIs in total of 10 cells.
%                           Default number of the ratio is 0; 
%       'ylimit'           : y-axis limit for the plots.
%  Output:
%       cluster_exp : expression matrix of gene(row)-by-cluster(column,
%                     C1-Cn). The number represents the average of UMIs of
%                     each cell cluster after median-normalization (column).
%       gene_cluster_idx : cluster index of all genes.


%% defaul settings
default_include_unexpressed_gene = 'no';
valid_include_unexpressed_gene = {'no','yes'};
check_include_unexpressed_gene = @(x) any(validatestring(x,valid_include_unexpressed_gene));

default_expression_cutoff_n  = 3;
default_expression_cutoff_ratio  = 0;
default_ylimit = [-4 4];

p = inputParser;
addOptional(p,'include_unexpressed_gene',default_include_unexpressed_gene,check_include_unexpressed_gene)
addOptional(p,'ylimit',default_ylimit)
addParameter(p,'expression_cutoff_n',default_expression_cutoff_n,@isnumeric)
addParameter(p,'umi_ratio_cutoff',default_expression_cutoff_ratio,@isnumeric)
parse(p,varargin{:})

default_include_unexpressed_gene = p.Results.include_unexpressed_gene;
disp([' Include unexpressed gene as an extra cluster?  ',default_include_unexpressed_gene])
default_ylimit = p.Results.ylimit;
default_expression_cutoff_n = p.Results.expression_cutoff_n;
default_expression_cutoff_ratio = p.Results.umi_ratio_cutoff;
if strcmp('yes',cellstr(default_include_unexpressed_gene))
    fprintf(' Expressed gene: at least express in %d cells.\n',default_expression_cutoff_n)
end


%% merge cells and setup inputs
stage = cell_cluster;
A   =   F_normalize(input_matrix,'median');
A_1 =       A;
A_1(A_1~=0)=1;
%Define expression level of gene expression across main cell types.
clear cluster_expression
for i= 1:length(unique(stage))
    cluster_expression(:,i) = mean(A(:,find(stage==i))')';
    cluster_expression_UMI(:,i) = mean(input_matrix(:,find(stage==i))')';
end

if strcmp('no',cellstr(default_include_unexpressed_gene))
    B = cluster_expression;
elseif strcmp('yes',cellstr(default_include_unexpressed_gene))
    %index of expressed/unexpressed protein-coding genes
    exp_idx1   = find(sum(A_1')>=default_expression_cutoff_n)';
    exp_idx2   = find(max(cluster_expression_UMI')>=default_expression_cutoff_ratio)';
    exp_idx    = intersect(exp_idx1,exp_idx2);
    unexp_idx  = setdiff(1:length(A(:,1)), exp_idx);
    fprintf('Detected %d unexpressed genes as Cluster 0.\n',length(unexp_idx))
    
end


%% kmeans clustersing of genes

% zscore transformation of expressed genes to find gene clusters.
%B_zscore  = zscore( (median(sum(cluster_expression))*bsxfun(@rdivide,cluster_expression,sum(cluster_expression)))' )';
B_zscore = zscore(cluster_expression')';
exp       = B_zscore(exp_idx,:);
% kmeans clustering of expressed prot_coding genes
rng(1)  % for reproducibility
km_exp_clusters = kmeans(exp,k); %kmeans into k clusters


%% visulization


%plot all genes (including unexpressed genes)
cluster_index = zeros(length(A(:,1)),1);
cluster_index(exp_idx) = km_exp_clusters;

fprintf('Detected %d clusters; Number of each cluster is: \n',length(unique(cluster_index)))
    j = length(unique(cluster_index));
if length(unexp_idx)==0
    for i=1:j
    fprintf('          Cluster %d:    %d \n',i,length(find(cluster_index==i)))
    end
elseif length(unexp_idx)>0
    fprintf('          Cluster 0:    %d \n',length(unexp_idx))
    for i=1:(j-1)
    fprintf('          Cluster %d:    %d \n',i,length(find(cluster_index==i)))
    end
end

%setup colormaps for clusting
cd ./cbrewer
RdBu_cm = flipud(cbrewer('div','RdBu',50));
Set3_cm = flipud(cbrewer('qual','Set3',length(unique(cluster_index))));
cd ..



%heatmap
[i,xi]=sort(cluster_index);
xticklength=length(B_zscore(1,:));

figure;
subplot('position',[0.25 0.1 0.7 0.8]);
imagesc(B_zscore(xi,:));colorbar;colormap(RdBu_cm)
xlim([0.5, xticklength+0.5]);
XTick = 1:xticklength;
xlabel('Cell clusters');
title('Gene expression by gene groups')
%draw heatmap of gene clusters
f=subplot('position',[0.05 0.1 0.03 0.8]);
imagesc(i);xticks([]);yticks([]);colormap(f,Set3_cm);
ylabel('Gene clusters');


% plot expressed gene modules
a=floor(sqrt(length(unique(cluster_index))))+1;
j = length(unique(cluster_index));

figure;
if length(unexp_idx)==0
for i = 1:j
    subplot(a,a,i);
    y = B_zscore(find(cluster_index==i),:);
    x =  1:xticklength;
    shadedErrorBar(x, mean(y,1),std(y), 'lineprops',{'b-o','markerfacecolor','g'});
    xlim([0.5, xticklength+0.5]);
    ylim(default_ylimit);
    h_gca = gca;
    h_gca.XTick = 1:xticklength;
    title(sprintf('Cluster %d (n=%d)',i,length(find(cluster_index==i))));
end

elseif length(unexp_idx)>0
for i = 1:j
    subplot(a,a,i);
    y = B_zscore(find(cluster_index==(i-1)),:);
    x =  1:xticklength;
    shadedErrorBar(x, mean(y,1),std(y), ...
        'lineprops',{'c-o','MarkerSize',3,...
        'Color',[214,96,77]/255,...
        'MarkerFaceColor',[178,24,43]/255});
    ylim(default_ylimit);    
    xlim([0.5, xticklength+0.5]);
    h_gca = gca;
    h_gca.XTick = 1:xticklength;
    title(sprintf('Cluster %d (n=%d)',i-1,length(find(cluster_index==(i-1)))));
end
end
% figure;hold on;
% shadedErrorBar(x, y, {@mean, @(x) 1*std(x)}, {'-b', 'LineWidth', 2}, 0);


%% output results
cluster_exp      = cluster_expression;
cluster_exp_UMI  = cluster_expression_UMI;
gene_cluster_idx = cluster_index;
