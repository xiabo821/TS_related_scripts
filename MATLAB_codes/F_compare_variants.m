function OUTPUT_STRUCTURE = F_compare_variants(variants_number,variants_genelength,variants_genename,gene_cluster,input_genename,color_map,varargin)
%% OUTPUT_STRUCTURE = F_compare_variants(variants_number,variants_genelength,variants_genename,...
%                                        gene_cluster,input_genename,color_map,varargin);
%  F_compare_variants is wrote for 'transcriptional-scanning' manuscript.
%  Overall, this function use [variants] as input, [gene_cluster] as
%  grouping/staging method, [reference_group] as the reference for
%  comparison. Then this function plots the results as defined below and
%  generates ranksum test p_values with Bonferroni correction.
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
%       'ylimit'          : y-axis limit for the plots. Default is the 99%
%                           of the input data list.
%       'visualization'   : Show boxplots results of variants or not.
%                           Default 'yes'. Optional, 'no'
%       'plot_distribution': Show plots of kernel distribution of vatiants
%                            or not. default 'no'.
%                            Optional, 'yes'
%  Output:
%       OUTPUT_STRUCTURE, containing the below three variables:
%       ranksum_p     : ranksum test p_values compared with unexpressed
%                       genes. The p_value undregoes Bonferroni adjustment
%                       by multiplying the number of tests in each category.
%       mean_value    : average value of the mutation rates
%       ci99          : 99% confidence intervals calculated by bootstrap
%                       method.
%  Bo Xia



%% Parse defaul settings and setup inputs
defaultVisualization = 'yes';
validVisualization = {'no','yes'};
checkVisualization = @(x) any(validatestring(x,validVisualization));

defaultPlotdistribution = 'no';
validPlotdistribution = {'no','yes'};
checkPlotdistribution = @(x) any(validatestring(x,validPlotdistribution));

default_ylimit = [0 300.1];

p = inputParser;
addOptional(p,'visualization',defaultVisualization,checkVisualization)
addOptional(p,'plot_distribution',defaultPlotdistribution,checkPlotdistribution)
addOptional(p,'ylimit',default_ylimit)
parse(p,varargin{:})

default_ylimit = p.Results.ylimit;
defaultVisualization = p.Results.visualization;
defaultPlotdistribution = p.Results.plot_distribution;

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


%For checking individual gene mutation rates.
check_individual_gene = 0;
if check_individual_gene
    checklist = {'SRY'}
    human_gene_info = struct2table(tdfread('C:\Users\Bo\Box Sync\Pap33_Spermatogenesis\Analysis\2017-11-03_human_TESE\Original file\human_genenames.tsv'));
%     keys = strmatch('HLA',cellstr(human_gene_info.name))
%     cellstr(human_gene_info.name(keys,:))
    [~,~,keys] = intersect(checklist,cellstr(human_gene_info.name),'stable'); 
    [~,~,keys] = intersect(cellstr(human_gene_info.ensembl_id(keys,:)),vgname(j),'stable'); 
    variants(keys)
    figure;boxplot(variants(keys))
    mean(variants(keys))
    median(variants(keys))
end

%compute the stats for each expression categories
%plot distribution of variants
for b=1:length(unique(gclu))
    Clu = variants(find(gclu==(b-1)));
    Clu = Clu( ~any( isnan( Clu ) | isinf( Clu ), 2 ) ,: );  Clu(isoutlier(Clu))=[];
    ranksumtest_p(b) = ranksum(C0,Clu) * (length(unique(gclu))-1) ;   %multiply by number of tests, Bonferroni adjustment
    Global_mean(b) = mean(Clu);  %mean of variants
    Global_median(b) = median(Clu);  %mean of variants
    Global_ci99(b,1:2) = bootci(10000,{@mean,Clu},'alpha',0.01)-mean(Clu); %calculate 99% CI by bootstrap with 10000times.
      %plot kernel density of each category
%       [f,xi] = ksdensity(Clu); 
%       plot(xi,f,'Color',cm(b,:),'LineWidth',1.5);hold on
end
  hold off;
  legend({'C0','C1','C2','C3','C4','C5','C6','C7','C8'});
ranksumtest_p(find(ranksumtest_p>1)) = 1;
    %output variable
     asym_score_marker=num2cell(ranksumtest_p);
%     asym_score_marker(find(ranksumtest_p<0.05)) = {'*'};
%     asym_score_marker(find(ranksumtest_p>=0.05)) = {'n.s.'};
%     asym_score_marker(1) = {''};

    

%% Visualization
n = length(unique(gclu));

if strcmp('yes',cellstr(defaultVisualization))
    [x,y]=sort(gclu);%Sort the stage number so that to do grouped boxplot.
    
    figure;  % -->to allow subplots
    boxplot(variants(y),x,'Symbol','','Notch','on','Colors',cm);hold on;
%    violinplot(variants(y),x);hold on;

    errorbar(1:n,Global_mean(:),Global_ci99(:,1),Global_ci99(:,2),...
        '-','Color',[0.1137 0.5686 0.7529],'Marker','.','MarkerSize',10,'LineWidth',1,'CapSize',4);hold on;
    if default_ylimit(2) ==300.1
        default_ylimit = [0 prctile(variants, 99)];
    end
    ylim(default_ylimit); 
    xlim([0.3 n+0.7]);
    xticks(1:n);
    xticklabels({'C0','C1','C2','C3','C4','C5','C6','C7','C8'});
    text(1:n,repelem(max(default_ylimit)*0.95,n),asym_score_marker);hold on;
    text(1,max(default_ylimit)*0.85,strcat('n= ',num2str(histcounts(gclu)))); 
    hold off;
    set(gca,'color','none');
    box off
end



if strcmp('yes',cellstr(defaultPlotdistribution))
    figure;
    for b=1:length(unique(gclu))
        Clu = variants(find(gclu==(b-1)));
        Clu = Clu( ~any( isnan( Clu ) | isinf( Clu ), 2 ) ,: );  Clu(isoutlier(Clu))=[];
        %plot kernel density of each category
        [f,xi] = ksdensity(Clu);
        plot(xi,f,'Color',cm(b,:),'LineWidth',1.5);hold on
        line([Global_median(b) Global_median(b)],[0 0.04],'LineStyle','--','Color',cm(b,:),'LineWidth',0.5);hold on
        line([Global_mean(b) Global_mean(b)],[0 0.04],'Color',cm(b,:),'LineWidth',0.5);hold on
    end
    legend({'C0','C1','C2','C3','C4','C5','C6','C7','C8'});
    for b=1:length(unique(gclu))
        line([Global_median(b) Global_median(b)],[0 0.04],'LineStyle','--','Color',cm(b,:),'LineWidth',0.5);hold on
        line([Global_mean(b) Global_mean(b)],[0 0.04],'Color',cm(b,:),'LineWidth',0.5);hold on
    end
    set(gca,'color','none');
        hold off; box off
end
%% outputs
OUTPUT_STRUCTURE.ranksumtest_p= ranksumtest_p';
OUTPUT_STRUCTURE.mean = Global_mean';
OUTPUT_STRUCTURE.ci99 = Global_ci99';