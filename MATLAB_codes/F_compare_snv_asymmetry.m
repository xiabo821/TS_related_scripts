function OUTPUT_STRUCTURE = F_compare_snv_asymmetry(...
                                    fn_snv_type_number,fn_base_frequency,...
                                    input_gene_cluster,input_genename,...
                                    varargin)
%% OUTPUT_STRUCTURE = F_compare_snv_asymmetry(...
%                                             fn_snv_number,...
%                                             fn_base_frequency,...
%                                             input_gene_cluster,...
%                                             input_genename,...
%                                             varargin);
%  F_compare_snv_asymmetry is wrote for 'transcriptional-scanning' manuscript.
%  Overall, this function use [SNV frequency] as input, [gene_cluster] as
%  grouping/staging method, [reference_group] as the reference for
%  comparison. Then this function plots the results of SNV frequency
%  between coding strand and template strand, as well as asymmetry scores,
%  and generates ranksum test p_values with Bonferroni correction.
%
%  required functions: No more!
%           barwitherr.m    Copyright@ Martina F. Callaghan
%
%  Input:
%       fn_snv_type_number: file name of .tsv file (genes_by_SNVtypenumbers).   
%                           First column should be gene names. Other rows
%                           contains SNV number of each mutation type, including: 
%                           A-T,A-G,A-C,T-A,T-G,T-C,G-A,G-T,G-C,C-A,C-T,C-G
%       fn_base_frequency : file name of .tsv file which contains A/T/G/C 
%                           numbers of a given gene/sequence. First column 
%                           should be genenames. Note: base_frequency gene name  
%                           and snv_number gene names should be the same type. 
%       input_gene_cluster: gene cluster indeces. The cluster should be 0-n
%                           integeters, where 0 stands for the unexpressed
%                           C0 cluster. 0 (C0) is required because C0 cluster will 
%                           be used as control for comparison and statistical tests.
%       input_genename    : gene names of input genes corresponding to
%                           cluster indeces. 
%                           Note 1: the snv_number gene name and input_genename
%                           should have the same name type. 
%                           Note 2: input_gene_cluster and input_genename
%                           should have exactly one-to-one matching to each other.
%  varargin:
%       'snv_asymmetry_order'   : SNV asymmetry calculation.
%                           default:{'AT','TA','AG','TC','TG','AC','CT','GA','GT','CA','CG','GC'}
%                           stands for: A>T/T>A; A>G/T>C; T>G/A>C; C>T/G>A; G>T/C>A; C>G/G>C;
%       'Mode'                  : Analysis mode, meaning calculating
%                                 mutation rates by single genes (default,'single_gene');
%        [This version          or by metagenes (option 'meta_INTEGER',INTEGER stands for an integer,
%        onlysupports            e.g.'meta_20', means merging mutation stats and reference 
%        'single_gene'               frequencies of meta genes for every random 20 genes);
%        mode.]                 or by adding mutation/base_frequency of each gene cluster 
%                                    all togeter (option, 'meta_all').           
%                                  
%       'Subtraction'           : Option for involving substraction for the
%                                 fn_snv_type_number input.
%                                 Default is 'no'. 
%                                 Note: if option 'yes', do need to assign
%                                 fn_snv_type_number_ref and
%                                 fn_base_frequency_ref, which should have
%                                 exactly the same data organization with the
%                                 fn_snv_type_number file and
%                                 fn_base_frequency  file, respectively.
%       'fn_snv_type_number_ref': The optional input file used as reference
%                                 to subtract from the the fn_snv_type_number
%                                 input. Note this file should to have
%                                 exactly the same data organization with the
%                                 fn_snv_type_number file.
%       'fn_base_frequency_ref' : The optional input file used as reference
%                                 to subtract from the fn_base_frequency
%                                 input. Note this file should to have
%                                 exactly the same data organization with the
%                                 fn_base_frequency file.
%       'FaceColor'             : color map used to bar-plot the asymmetry scores.
%       'Title_prefix'          : Set title prefix of all the plots. 
%  
%  Output: OUTPUT_STRUCTURE, containing the below sub_variables:
%       SNVtype_by_column  : The column names for each SNVtype corresponding
%                            to the values listed together. Same with the input SNV order. 
%       gene_cluster_name  : The gene clusters used in the analysis. 
%                            Same with the input gene cluster names.
%       CSvsTS_test_p      : paired sample t-test p_values for each SNVtype and each
%                            gene cluster. CS_TS_test_p < 0.0001 is also marked
%                            as '*' in the strand_SNV_level barwitherr plots.
%       expVSunexp_asym_score_test_p: ranksum test p_values of asymmetrys between
%                            expressed gene clusters vs unexpressed gene
%                            cluster. Asymmetry_test_p < 0.001 is also marked
%                            as '*' in the asymmetry_score barwitherr plots.
%       strand_asym_stat   : A x-y-z matrix. Average mutation rates (Z1) and  
%                            confidence intervals (Z2 & Z3) for specific gene 
%                            clusters (x-axis) distinguished by SNV type (y-axis). 
%       asym_score_stat    : A x-y-z matrix. Average asymmetry scores (Z1) and 
%                            their confidence intervals (Z2 & Z3) for specific 
%                            gene clusters (x-axis) distinguished by SNV type pairs (y-axis). 
%
%  Bo Xia

%% Parse defaul settings and setup inputs
fprintf('****************************************\n')
time1 = datetime;
fprintf('%s     Analysis started.\n',datestr(time1))

%default asymmetry calculations: 
%                       A>T/T>A;A>G/T>C;T>G/A>C;C>T/G>A;G>T/C>A;C>G/G>C;
default_snv_asymmetry = {'AT','TA','AG','TC','TG','AC','CT','GA','GT','CA','CG','GC'};
default_face_color = [0.5529    0.8275    0.7804;
    1.0000    1.0000    0.7020;
    0.7451    0.7294    0.8549;
    0.9843    0.5020    0.4471;
    0.5020    0.6941    0.8275;
    0.9922    0.7059    0.3843;
    0.7020    0.8706    0.4118;
    0.9882    0.8039    0.8980;
    0.8510    0.8510    0.8510;
    0.7373    0.5020    0.7412;
    0.8000    0.9216    0.7725;
    1.0000    0.9294    0.4353];

default_Mode         = 'single_gene';
defaultSubtraction = 'no';
validSubtraction = {'no','yes'};
checkSubtraction = @(x) any(validatestring(x,validSubtraction));

defaultfn_snv_type_number_ref = {''};
defaultfn_base_frequency_ref = {''};
default_title_prefix         = '';

p = inputParser;
addOptional(p,'snv_asymmetry_order',default_snv_asymmetry)
addOptional(p,'Mode',default_Mode)
addOptional(p,'Subtraction',defaultSubtraction,checkSubtraction)
addOptional(p,'fn_snv_type_number_ref',defaultfn_snv_type_number_ref)
addOptional(p,'fn_base_frequency_ref',defaultfn_base_frequency_ref)
addOptional(p,'FaceColor',default_face_color)
addOptional(p,'Title_prefix',default_title_prefix)
parse(p,varargin{:})

%default numbers
default_snv_asymmetry = p.Results.snv_asymmetry_order;
default_Mode = p.Results.Mode;
defaultSubtraction = p.Results.Subtraction;
defaultfn_snv_type_number_ref = p.Results.fn_snv_type_number_ref;
defaultfn_base_frequency_ref = p.Results.fn_base_frequency_ref;
default_face_color = p.Results.FaceColor;
default_title_prefix    = p.Results.Title_prefix;

%setup input files
    %[SNV_number,SNV_genenames,SNV_type]
    [vnum,vgnames,vtype_names] = tsv_import(fn_snv_type_number);
    %[base_frequency,genenames,base_type]
    [vglength,vglenames,base_names] = tsv_import(fn_base_frequency);
    
    %input sample check if involves subtraction of input SNV_type numbers.
    if strcmp('yes',cellstr(defaultSubtraction))
        if strcmp('',cellstr(defaultfn_snv_type_number_ref)) | strcmp('',cellstr(defaultfn_base_frequency_ref)) 
            fprintf('Need to assign the "Subtraction" reference inputs.\n')
            return
        else
            [vnum_ref,vgnames_ref,vtype_names_ref] = tsv_import(defaultfn_snv_type_number_ref);
            [vglength_ref,vglenames_ref,base_names_ref] = tsv_import(defaultfn_base_frequency_ref);
            [~,a,b]=intersect(vgnames,vgnames_ref,'stable');
            [~,c,d]=intersect(vtype_names,vtype_names_ref,'stable');
            [~,e,f]=intersect(vglenames,vglenames_ref,'stable');
            [~,g,h]=intersect(base_names,base_names_ref,'stable');
            if sum(a-b)==0 & sum(c-d)==0 
                vnum     = vnum - vnum_ref;
            else
                fprintf('Inputs mismatch error! \n fn_snv_type_number and fn_snv_type_number_ref should have identical data organization.\n')
                return
            end
            if sum(e-f)==0 & sum(g-h)==0 
                vglength = vglength - vglength_ref;
            else
                fprintf('Inputs mismatch error! \n fn_base_frequency and fn_base_frequency_ref should have identical data organization.\n')
                return
            end
        end
    end

    %other default values
    gclu   = input_gene_cluster;
    gname  = cellstr(input_genename);
    asym_order  = default_snv_asymmetry;
    face_cm     = default_face_color;

%% set up SNV frequency variants

[valid_genenames,i,j]=intersect(vgnames,vglenames,'stable');
vnum     = vnum(i,:);
vglength = vglength(j,:);
if length(j)==0
    fprintf('Error! Mismatched gene name of input SNV/base_frequency files!')
    return
end
%calculate SNV frequency per kilo-bases for each gene.
vnum_perkb = [];
vnum_perkb(:,strcmp('AT',asym_order)) = vnum(:,strcmp('AT',vtype_names))...
                                        ./   vglength(:,strcmp('A',base_names)) *1000;
vnum_perkb(:,strcmp('AG',asym_order)) = vnum(:,strcmp('AG',vtype_names))...
                                        ./   vglength(:,strcmp('A',base_names)) *1000;
vnum_perkb(:,strcmp('AC',asym_order)) = vnum(:,strcmp('AC',vtype_names))...
                                        ./   vglength(:,strcmp('A',base_names)) *1000;
vnum_perkb(:,strcmp('TA',asym_order)) = vnum(:,strcmp('TA',vtype_names))...
                                        ./   vglength(:,strcmp('T',base_names)) *1000;
vnum_perkb(:,strcmp('TG',asym_order)) = vnum(:,strcmp('TG',vtype_names))...
                                        ./   vglength(:,strcmp('T',base_names)) *1000;
vnum_perkb(:,strcmp('TC',asym_order)) = vnum(:,strcmp('TC',vtype_names))...
                                        ./   vglength(:,strcmp('T',base_names)) *1000;
vnum_perkb(:,strcmp('GA',asym_order)) = vnum(:,strcmp('GA',vtype_names))...
                                        ./   vglength(:,strcmp('G',base_names)) *1000;
vnum_perkb(:,strcmp('GT',asym_order)) = vnum(:,strcmp('GT',vtype_names))...
                                        ./   vglength(:,strcmp('G',base_names)) *1000;
vnum_perkb(:,strcmp('GC',asym_order)) = vnum(:,strcmp('GC',vtype_names))...
                                        ./   vglength(:,strcmp('G',base_names)) *1000;
vnum_perkb(:,strcmp('CA',asym_order)) = vnum(:,strcmp('CA',vtype_names))...
                                        ./   vglength(:,strcmp('C',base_names)) *1000;
vnum_perkb(:,strcmp('CT',asym_order)) = vnum(:,strcmp('CT',vtype_names))...
                                        ./   vglength(:,strcmp('C',base_names)) *1000;
vnum_perkb(:,strcmp('CG',asym_order)) = vnum(:,strcmp('CG',vtype_names))...
                                        ./   vglength(:,strcmp('C',base_names)) *1000;

%remove NaN and Inf numbers, which indicates no ref sequence stats from the data
valid_genenames = valid_genenames( ~any( isnan( vnum_perkb ) | isinf( vnum_perkb ), 2 ) );
vnum_perkb = vnum_perkb( ~any( isnan( vnum_perkb ) | isinf( vnum_perkb ), 2 ) ,: );



%% Compute SNV per kb per gene strand and paired_sample t-test between coding and template strands
%setup the input gene list
[valid_geneid,a,b]=intersect(gname,valid_genenames,'stable');
gclu= gclu(a);
vnum_perkb = vnum_perkb(b,:);

%check gene names
fprintf('  %d out of %d input_genes are included for downstream analysis!\n',...
                    length(valid_geneid), length(gname)  );
if length(a)==0
    fprintf('Error! Input gene name mismatches!')
    return
end
if min(gclu)~=0
    fprintf('Error! input_gene_cluster does not have C0!')
    return
end



%Calculate strand asymmetry test p_values (use paired-sample t-test with Bonferroni adjustment)
for i = 1:(length(asym_order)./2)
    mut_type = char(asym_order((i*2-1):i*2));
    fprintf('Started processing %s>%s/%s>%s mutation pair!\n',...
             mut_type(1,1),mut_type(1,2),mut_type(2,1),mut_type(2,2));
    asym_type(i) = cellstr(sprintf('%s>%s/%s>%s',mut_type(1,1),mut_type(1,2),mut_type(2,1),mut_type(2,2)));
    asym_type_2(i) = strcat(asym_order(2*i-1),{'_'},asym_order(2*i));
    for j = unique(gclu)'
        gc = find(gclu==j);
        %strand_asym_ranksum_p(j+1,i) = ranksum(vnum_perkb(gc,2*i-1),vnum_perkb(gc,2*i)) ...
                                       %* length(unique(gclu));
        cs = vnum_perkb(gc,2*i-1);
        ts = vnum_perkb(gc,2*i);
        [~,cs_out,~,~,~] = filloutliers(cs,'previous');
        [~,ts_out,~,~,~] = filloutliers(ts,'previous');
        quolif_gene = setdiff(1:length(cs),union(cs_out,ts_out));
        %remove outlier for doing paired-sample ttest.
        [~,strand_asym_p_value(j+1,i)] = ttest(cs(quolif_gene),ts(quolif_gene),'Alpha',0.01);
        %[~,strand_asym_p_value(j+1,i)] = ttest(vnum_perkb(gc,2*i-1),vnum_perkb(gc,2*i),'Alpha',0.01);
        strand_asym_p_value(j+1,i) = strand_asym_p_value(j+1,i) * length(unique(gclu));
        %multiply by number of tests, Bonferroni adjustment
        gc_name(j+1) = strcat({'C'},num2str(j));    
    end
end
    %export the paired-sample t-test pvalues of strand SNV_levels
    strand_asym_p_value(find(strand_asym_p_value>1)) = 1;
    SNV_strand_asym_p_value = array2table(strand_asym_p_value,...
                                  'VariableNames',asym_type_2,...
                                  'RowNames',gc_name);
    SNV_strand_asym_marker=num2cell(strand_asym_p_value);
    SNV_strand_asym_marker(find(strand_asym_p_value<=0.000001)) = {'**'};
    SNV_strand_asym_marker(find(strand_asym_p_value<=0.01 & strand_asym_p_value>0.000001)) = {'*'};
    SNV_strand_asym_marker(find(strand_asym_p_value>0.01)) = {'n.s.'};

    
%% Calculate strand asymmetry stats and VISUALIZATION
%calculate strand asymmetry stats 
for i = 1: length(asym_order)
    for j = unique(gclu)'
        vgc = vnum_perkb(find(gclu==j),i);
        strand_asym_stat(j+1,i,1) = mean(vgc);
        ci = bootci(10000,{@mean,vgc},'alpha',0.01)-mean(vgc);  %mean of strand_SNV
        strand_asym_stat(j+1,i,2) = ci(1);  %lower CI-99% of strand_SNV
        strand_asym_stat(j+1,i,3) = ci(2);  %upper CI-99% of strand_SNV
    end
end

%plot with barwitherr.m
l = length(asym_order)/2; %plot numbers
n = length(unique(gclu)); %gene cluster numbers
Z =  [1 2;3 4;5 6;7 8;9 10;11 12];
figure;
for a=1:l
    subplot(1,l,a);
    barplot = barwitherr(strand_asym_stat(:,Z(a,:),2:3),strand_asym_stat(:,Z(a,:),1),0.9);
    barplot(1).FaceColor = [244 165 130]/256;
    barplot(2).FaceColor = [146 197 222]/256;
    xlim([0.5 n+0.5]); hold on;
    line_y = mean(strand_asym_stat(1,Z(a,:),1));
    ylim([0 line_y*1.25]);
    line([0.5 n+0.5],[line_y line_y],'Color',[227 26 28]/256,'LineStyle','--');
    text(1:n,repelem(line_y*1.09,n),SNV_strand_asym_marker(:,a));hold off;
    title(asym_type(a));box off;
    xticklabels(gc_name);
    set(gca,'color','none');
end
suptitle([default_title_prefix ' [Germline mutation rates (per kb)]'])

time2 = datetime;
fprintf('%s     Finished the first part.\n',datestr(time2))


%% asymmetry score calculations, ranksum test, and VISUALIZATION
%Calculate asymmetry scores
for i = 1:(length(asym_order)./2)
    asymmetry_score(:,i) = log2((vnum_perkb(:,2*i-1)+1) ./ (vnum_perkb(:,2*i)+1) );
end

%Calculate strand asymmetry ranksum test p_values 
for i = 1: (length(asym_order)./2)
    vgc0 = asymmetry_score(find(gclu==0),i);
        for j = unique(gclu)'
            vgc = asymmetry_score(find(gclu==j),i);
            % ranksum test of asymmetry scores versus unexpressed gene cluster.
            asym_score_ranksum_p(j+1,i) = ranksum(vgc0,vgc)...
                                          * (length(asym_order)./2 - 1);
                                          %multiply by number of tests, Bonferroni adjustment
            asym_score_stat(j+1,i,1) = mean(vgc);
            ci = bootci(10000,{@mean,vgc},'alpha',0.01)-mean(vgc);  %mean of asymetry_score
            asym_score_stat(j+1,i,2) = ci(1);  %lower CI-99% of asymetry_score
            asym_score_stat(j+1,i,3) = ci(2);  %upper CI-99% of asymetry_score
        end
end

    %export the ranksum test pvalues comparing asymmetry scores bwteen
    %expressed and unexpressed gene clusters.
    asym_score_ranksum_p(find(asym_score_ranksum_p>1)) = 1;
    %output variable
    asym_score_ranksum_test_p = array2table(asym_score_ranksum_p,...
                                  'VariableNames',asym_type_2,...
                                  'RowNames',gc_name);
    asym_score_marker=num2cell(asym_score_ranksum_p);
    asym_score_marker(find(asym_score_ranksum_p<0.000001)) = {'**'};
    asym_score_marker(find(asym_score_ranksum_p<=0.01 & asym_score_ranksum_p>0.000001)) = {'*'};
    asym_score_marker(find(asym_score_ranksum_p>0.01)) = {'n.s.'};
    asym_score_marker(1,:) = {''};

%plot with barwitherr.m
l = length(asym_order)/2; %plot numbers
n = length(unique(gclu)); %gene cluster numbers
ylimit(1) = min(min(asym_score_stat(:,:,1)+asym_score_stat(:,:,2)))*1.15;
ylimit(2) = max(max(asym_score_stat(:,:,1)+asym_score_stat(:,:,3)))*1.1;
% face_cm = [[0.50,0.50,0.50;0.40,0.76,0.647;0.988,0.553,0.384;...
%             0.5529,0.62745,0.796;0.9058,0.5412,0.7647;0.651,0.847,0.329]];
figure;
for a=1:l
    subplot(1,l,a);    
    for i = 1:n
        bar(i,asym_score_stat(i,a,1),0.8,'FaceColor',face_cm(i,:));hold on;
    end
    hold on;
    errorbar(1:n,asym_score_stat(:,a,1),asym_score_stat(:,a,2),asym_score_stat(:,a,3),'.k');
     hold on;
    xlim([0.5 n+0.5]); 
    ylim(ylimit);
    text(1:n,repelem(max(max(asym_score_stat(:,:,1)))*1.05,n),asym_score_marker(:,a));hold off;
    title(asym_type(a));box off;
    xticks(1:n);   xticklabels(gc_name);
    set(gca,'color','none');
end
suptitle([default_title_prefix ' [Asymmetry scores]'])

time3=datetime;
fprintf('%s     Finished the second part.\n',datestr(time3))
fprintf('Check out the plots!\n ')
fprintf('The analysis took: ')
time3-time1
fprintf('****************************************\n')


%% outputs
OUTPUT_STRUCTURE.SNVtype_by_column           = default_snv_asymmetry;
OUTPUT_STRUCTURE.gene_cluster_name           = gc_name;
OUTPUT_STRUCTURE.CSvsTS_test_p               = SNV_strand_asym_p_value;
OUTPUT_STRUCTURE.expVSunexp_asym_score_test_p= asym_score_ranksum_test_p;
OUTPUT_STRUCTURE.strand_asym_stat            = strand_asym_stat;
OUTPUT_STRUCTURE.asym_score_stat             = asym_score_stat;



%% functions
function [matrix,row_names,column_names] = tsv_import(filename)
%%import filename.tsv as three matlab variables:
% filename should be like: 'filename.tsv'
% 1. matrix(gene*cell)
% 2. row_names
% 3. column_names   
%   First (n) columns can be gene info and should have one column titled as 
%   'gene_id', for which be exported as the 'row_names';
%   The rest of the columns (from n+1 to end) should contain SNV number or
%   reference types, for which will be exported as the 'column_names'.

%%-->import all data

imported_data= importdata(filename);
matrix       = imported_data.data;
row_names    = imported_data.textdata(2:end,strcmp('gene_id',imported_data.textdata(1,:)));
column_names = imported_data.textdata(1,strcmp('',imported_data.textdata(2,:)));