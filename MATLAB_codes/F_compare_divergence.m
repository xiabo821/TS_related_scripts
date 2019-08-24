function [ranksumtest_p,varargout] = F_compare_divergence(variants_number,variants_genename,...
                                                        gene_cluster,input_genename,...
                                                        varargin)
%% [ranksumtest_p,varargout] = F_compare_divergence(variants_number,variants_genename,...
%                                                   gene_cluster,input_genename,...
%                                                   varargin);
%  F_compare_divergence is wrote for 'transcriptional-scanning' manuscript.
%  Overall, this function use [variants] as input, [gene_cluster] as
%  grouping/staging method, [reference_group] as the reference for
%  comparison. Then this function plots (errorbar) the results as defined below and
%  generates ranksum test p_values with Bonferroni correction.
%   To plot results, please add 'figure' function before 'F_compare_divergence'
%   otherwise the funtion will overwrite previous results.
%
%  required functions:
%
%  Input:
%       variants_number   : divergence/dN/dS number or anything else of
%                           each gene on the rows. This can be a single column 
%                           of numbers or multiple columns. Each column will
%                           produce a errorbar plot curve 
%       variants_genename : gene names of corresponding variants
%       input_gene_cluster: gene cluster indeces. The cluster should be 0-n
%                          integeters, where 0 stands for the unexpressed
%                          C0 cluster. C0 cluster will be used as control
%                          for comparison and test.
%       input_genename    : gene names of input genes corresponding to
%                           cluster indeces. 
%                           Note: the variants_genename and input_genename
%                           should have the same name type.
%  varargin:
%       'ylimit'          : y-axis limit for the plots.
%       'color_map'       : color map used to plot the variants by
%                           clusters. RGB colormap numbers, [0 1], with 3 columns.
%       'legend_name'     : cell array of columns names for variants_number
%                           input matrix.
%       'show_significance': default 'yes'. If 'no', will show the
%                            significance as * for ranksum test p_value <0.001
%                            If containing multiple rows, the *s from top
%                            to bottom lays correspond to input variance
%                            column 1 to the end.
%                           
%  Output:
%       ranksum_p     : ranksum test p_values compared with unexpressed
%                       genes. The p_value undregoes Bonferroni adjustment
%                       by multiplying the number of tests in each category.
%  Bo Xia
%


%% Parse defaul settings and setup inputs

default_ylimit = [0 300];
default_colormap = [];
default_legendname = '';
defaultSignificance = 'yes';
validSignificance = {'no','yes'};
checkSignificance = @(x) any(validatestring(x,validSignificance));

p = inputParser;
addOptional(p,'ylimit',default_ylimit,@isnumeric)
addOptional(p,'color_map',default_colormap)
addOptional(p,'legend_name',default_legendname)
addOptional(p,'show_significance',defaultSignificance,checkSignificance)

parse(p,varargin{:})

%parameters
default_ylimit      = p.Results.ylimit;
default_colormap    = p.Results.color_map;
default_legendname  = p.Results.legend_name;
%inputs
    vnum   = variants_number;
        v_type_num = length(vnum(1,:));
    vgname = variants_genename;
    gclu   = gene_cluster;
    gname  = input_genename;

%% Compute variants per kb

% Compute the variants
[~,i,j]=intersect(gname,vgname,'stable');
if isempty(j)
    fprintf('Error! Gene name mismatches!')
    return
end

variants = vnum(j,:);
gclu = gclu(i);
%Compute the unexpressed group for statistical test
C0 = variants(find(gclu==0));
C0 = C0( ~any( isnan( C0 ) | isinf( C0 ), 2 ) ,: );  C0(isoutlier(C0))=[];

            
%compute the stats for each expression categories
for a=1:v_type_num
    for b=1:length(unique(gclu))
        Clu = variants(find(gclu==(b-1)),a);
        Clu = Clu( ~any( isnan( Clu ) | isinf( Clu ), 2 ) ,: );  Clu(isoutlier(Clu))=[];
        ranksumtest_p(b,a) = ranksum(C0,Clu) * (length(unique(gclu))-1) ;   %multiply by number of tests, Bonferroni adjustment
        Global_mean(b,a) = mean(Clu);  %mean of variants
        Global_ci99(b,1:2,a) = bootci(10000,{@mean,Clu},'alpha',0.01)-mean(Clu); %calculate 99% CI by bootstrap with 10000times.
    end
end

ranksumtest_p(find(ranksumtest_p>1)) = 1;
    %output variable
    asym_score_marker=num2cell(ranksumtest_p);
    asym_score_marker(find(ranksumtest_p<0.001)) = {'*'};
    asym_score_marker(find(ranksumtest_p>=0.001)) = {'n.s.'};
    
    
    
    


%% Visualization
n = length(unique(gclu));
    
    %figure;  -->to allow subplots
    for b=1:v_type_num
        f=errorbar(1:n,Global_mean(:,b),Global_ci99(:,1,b),Global_ci99(:,2,b),'LineWidth',1);hold on
        xlim([0.5 n+0.5]); 
        if ~isempty(default_colormap)
           f.Color = default_colormap(b,:);
        end
    end
    
    if ~isempty(default_legendname)
            legend(default_legendname(1:v_type_num));
        elseif ~isempty(default_ylimit)
            ylim(default_ylimit);
        end
    if strcmp('yes',cellstr(defaultSignificance))
    h = max(ylim);
    for i=1:v_type_num
        text(1:n,repelem( h*(1-i*0.02),n),asym_score_marker(:,i));
    end
    end
    xticks(1:n);xticklabels({'C0','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10'});
    set(gca,'color','none');
    box off;
    hold off;
    


%% outputs
ranksumtest_p;
varargout{1} = Global_mean;
varargout{2} = Global_ci99;