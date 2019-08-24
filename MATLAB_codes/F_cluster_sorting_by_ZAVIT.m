function [cluster_idx] = F_cluster_sorting_by_ZAVIT(input_matrix,gene_cluster,angle,varargin)
%% [cluster_idx] = F_cluster_sorting_by_ZAVIT(input_matrix,gene_cluster,angle,varargin)
%  F_cluster_sorting_by_ZAVIT uses F_cluster_gene results to sort the gene
%  clusters by zavit, then sort the cluster and plot.
%
%  required functions:
%       plot_areaerrorbar.m
%       cbrewer.m
%       zavit_angle.m
%  Input: 
%       input_matrix: cluster-merged expression matrix from F_cluster_gene
%                     outputs. The matrix is the average UMI matrix of each
%                     cluster. The F_cluster_sorting_by_ZAVIT funtion will
%                     firstly calculate the zscore of genes across all
%                     stages before ordering the gene clusters.
%       gene_cluster: gene cluster indeces from F_cluster_gene outputs. The
%                     cluster could be 1-n (without unexpressed C0 cluster),
%                     or 0-n integeters (including unexpressed C0 cluster).
%       angle       : angle for zavit ordering, ranges from 0-360
%  varargin:
%       'ylimit'           : y-axis limit for the plots.
%  Output:
%       cluster_idx : redefined cluster index number. Same clusters with
%                     F_cluster_gene output 'cluster_idx', but just cluster
%                     number is different
%
%  Bo Xia
%


%% defaul settings and inputs

default_ylimit = [-4 4];

p = inputParser;
addOptional(p,'ylimit',default_ylimit)
parse(p,varargin{:})

default_ylimit = p.Results.ylimit;

% inputs
exp = zscore(input_matrix')';
cluster =  gene_cluster;
angle   =  angle;

%% kmeans clustering of expressed prot_coding genes

%if clusters contain unexpressed genes as the C0 cluster
if min(cluster)==0
    cluster = cluster +1;
end

%Sort the the kmeans clusters by averaging expression of each cluster
for i = 1 : length(unique(cluster))
    exp_norm_clust (i,:) = mean (exp(find(cluster==i),:),1);
end

%Cluster order by ZAVIT_angle=angle for visulization
cluster_order = zavit_angle(exp_norm_clust,angle);
figure;imagesc (exp_norm_clust(cluster_order,:));
xlabel('cell cluster (merged)')
ylabel('gene clusters')
title ('Raw ZAVIT sorting')

%reorder the kmeans clusters
zavit_cluster = cluster;
for i = 1 : length(unique(cluster))
    zavit_cluster( find(cluster == cluster_order(i)) ) = i;
end

%return the cluster number. If C0 exists, remain the C0 cluster 0.
if min(gene_cluster)==0
    zavit_cluster = zavit_cluster - 1;
elseif min(gene_cluster)>0
    zavit_cluster = zavit_cluster;
end

%% visualization after zavit reordering

%colormap
cd ./cbrewer
RdBu_cm = flipud(cbrewer('div','RdBu',50));
Blues_cm = cbrewer('seq','Blues',50);
cd ..

xticklength=length(exp(1,:));
[i,xi]=sort(zavit_cluster);

figure;
%draw heatmap as the genes sorted;
subplot('position',[0.25 0.1 0.65 0.85]);
imagesc(exp(xi,:));colorbar;colormap(RdBu_cm)
h_gca = gca;
xlim([0.5, xticklength+0.5]);
h_gca.XTick = 1:xticklength;
xlabel('Cell clusters');
title('Gene expression by gene groups')
%draw heatmap of gene clusters
f=subplot('position',[0.1 0.1 0.05 0.85]);
imagesc(i);xticks([]);yticks([]);colormap(f,Blues_cm);
ylabel('Gene clusters');

%plot all genes cluster by cluter

% visualization all reordered genes
%a=ceil(sqrt(length(unique(zavit_cluster))));
j = length(unique(zavit_cluster));
n = unique(zavit_cluster);
figure;
for i = 1:j
    subplot(j,1,i);
    y = exp(find(zavit_cluster==n(i)),:);
    x =  1:xticklength;
    shadedErrorBar(x, mean(y,1),std(y), ...
        'lineprops',{'c-o','MarkerSize',3,...
        'Color',[214,96,77]/255,...
        'MarkerFaceColor',[178,24,43]/255});
    ylim(default_ylimit);    
    xlim([0.5, xticklength+0.5]);
    h_gca = gca;
    h_gca.XTick = 1:xticklength;
    title(sprintf('Cluster %d (n=%d)',n(i),length(find(zavit_cluster==n(i))))); 
end

%suptitle('Clusters of expressed genes')


%% output

    cluster_idx = zavit_cluster;

