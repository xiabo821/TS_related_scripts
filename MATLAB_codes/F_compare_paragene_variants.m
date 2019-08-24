function [ranksumtest_p,varargout] = F_compare_paragene_variants(variants_number,variants_genelength,variants_genename,gene_cluster,input_genename,color_map,varargin)
%% [ranksumtest_p,varargout] = F_compare_paragene_variants(variants_number,variants_genelength,variants_genename,...
%                        gene_cluster,input_genename,color_map,varargin);
%  F_compare_paragene_variants is wrote for 'transcriptional-scanning' manuscript.
%  Overall, this function use [variants] as input, [gene_cluster] as
%  grouping/staging method, [reference_group] as the reference for
%  comparison. Then this function plots the results as defined below and
%  generates ranksum test p_values with Bonferroni correction.
%   To plot results, please add 'figure' function before 'F_compare_paragene_variants'
%   otherwise the funtion will overwrite previous results.
%   This funciton is very similar with the F_compare_variants.m, with
%   difference in ploting styles.
%
%  required functions:
%
%  Input:
%       variants_number   : Variants number of each gene
%       variants_genelength : gene length of corresponding variants
%       variants_genename : gene names of corresponding variants
%       input_gene_cluster: gene cluster indeces. The cluster should be 0-n
%                          integeters, where 0 stands for the unexpressed
%                          C0 cluster. C0 cluster will be used as control
%                          for comparison and test.
%       input_genename    : gene names of input genes corresponding to
%                           cluster indeces. 
%                           Note: the variants_genename and input_genename
%                           should have the same name type.
%       color_map         : color map used to plot the variants by clusters
%  varargin:
%       'ylimit'          : y-axis limit for the plots.
%       'visualization'   : Show plots or not. default 'yes'.
%                           Optional, 'no'
%  Output:
%       ranksum_p     : ranksum test p_values compared with unexpressed
%                       genes. The p_value undregoes Bonferroni adjustment
%                       by multiplying the number of tests in each category.
%  Bo Xia



%% Parse defaul settings and setup inputs
defaultVisualization = 'yes';
validVisualization = {'no','yes'};
checkVisualization = @(x) any(validatestring(x,validVisualization));

default_ylimit = [0 300];

p = inputParser;
addOptional(p,'visualization',defaultVisualization,checkVisualization)
addOptional(p,'ylimit',default_ylimit)
parse(p,varargin{:})

default_ylimit = p.Results.ylimit;
defaultVisualization = p.Results.visualization;

    vnum = variants_number;
    vglength = variants_genelength;
    vgname = variants_genename;
    gclu   = gene_cluster;
    gname  = input_genename;
    cm     = color_map;

%% Compute variants per kb

% Compute the variants
[~,i,j]=intersect(gname,vgname,'stable');
if length(j)==0
    fprintf('Error! Gene name mismatches!')
    return
end
variants=vnum(j,:)./vglength(j,:)*1000;
gclu = gclu(i);
%Compute the unexpressed group for statistical test
C0 = variants(find(gclu==0));
C0 = C0( ~any( isnan( C0 ) | isinf( C0 ), 2 ) ,: );  C0(isoutlier(C0))=[];

%compute the stats for each expression categories
for b=1:length(unique(gclu))
    Clu = variants(find(gclu==(b-1)));
    Clu = Clu( ~any( isnan( Clu ) | isinf( Clu ), 2 ) ,: );  Clu(isoutlier(Clu))=[];
    ranksumtest_p(b) = ranksum(C0,Clu) * (length(unique(gclu))-1) ;   %multiply by number of tests, Bonferroni adjustment
    Global_mean(b) = mean(Clu);  %mean of variants
    Global_ci99(b,1:2) = bootci(10000,{@mean,Clu},'alpha',0.01)-mean(Clu); %calculate 99% CI by bootstrap with 10000times.
end

ranksumtest_p(find(ranksumtest_p>1)) = 1;
    %output variable
    asym_score_marker=num2cell(ranksumtest_p);
%     asym_score_marker(find(ranksumtest_p<0.01)) = {'*'};
%     asym_score_marker(find(ranksumtest_p>=0.01)) = {'n.s.'};
%     asym_score_marker(1) = {''};


    

%% Visualization
n = length(unique(gclu));

if strcmp('yes',cellstr(defaultVisualization))
    [x,y]=sort(gclu);%Sort the stage number so that to do grouped boxplot.
    for i = 1:n
        bar(i,Global_mean(i),0.7,'FaceColor',cm(i,:));hold on;
    end

    %figure;  -->to allow subplots
%    boxplot(1:n,variants(y),x,'Symbol','','Notch','on','Colors',cm);hold on;
%    violinplot(variants(y),x);hold on;
    errorbar(1:n,Global_mean(:),Global_ci99(:,1),Global_ci99(:,2),...
        '--','Color','k','Marker','.','MarkerSize',10,'LineWidth',1,'CapSize',4);hold on;
    ylim(default_ylimit);  
    xlim([0.3 n+0.7]);
    xticks(1:n);
    xticklabels({'C0','C1','C2','C3','C4','C5','C6','C7','C8'});
    text(1:n,repelem(max(default_ylimit)*0.95,n),asym_score_marker);hold off;

    set(gca,'color','none');
    box off
    
end

%% outputs
ranksumtest_p= ranksumtest_p';
varargout{1} = Global_mean';
varargout{2} = Global_ci99';