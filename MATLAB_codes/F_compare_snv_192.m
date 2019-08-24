function OUTPUT_STRUCTURE = F_compare_snv_192(fn_snv_type_number,...
                                           fn_base_frequency,...
                                           input_gene_cluster,...
                                           input_genename,...
                                           varargin)
%% F_compare_snv_192
%  OUTPUT_STRUCTURE = F_compare_snv_192(fn_snv_type_number,...
%                                    fn_base_frequency,...
%                                    input_gene_cluster,...
%                                    input_genename,...
%                                    varargin)
%
%  F_compare_snv_192 is wrote for 'transcriptional-scanning' manuscript.
%  Overall, this function use [SNV frequency] with immediate sequence 
%  context as input, [gene_cluster] as grouping/staging method. 
%  Then this function plots the results of SNV frequency by gene groups
%  between coding strand and template strand, as well as asymmetry scores.
%
%  required functions: 
%           seqrcomplement: From the MATLAB Bioinfomatic toolbox for nucleotide
%                           sequence analysis. Note: it works with char input type.
%                            
%  Input:
%       fn_snv_type_number: file name of .tsv file (genes_by_SNVtypenumbers).   
%                           First (n) columns can be gene info and should 
%                           include one column with column title as 'gene_id';
%                           The rest of the columns (from n+1 to end) should 
%                           contain SNV number of 192 mutation types with column titles, e.g.: 
%                           ATCG, indicating a T-to-G mutation where T is in the ATC context.
%       fn_base_frequency : file name of .tsv file which contains 64 types of XAY/XTY/XGY/XCY 
%                           tri-bucleotide frequencies in the reference genome of a given gene/sequence.  
%                           First (n) columns can be gene info and should 
%                           include one column with title as 'gene_id';
%                           The rest of the columns (from n+1 to end) should 
%                           contain 64 types of tri-bucleotide frequencies   
%                           with column titles as XAY/XTY/XGY/XCY, e.g.: 
%                           ATC, indicating a ATC tri-nucleotide reference where T is in the ATC context.
%       input_gene_cluster: gene cluster indeces. The cluster should be 0-n
%                           integeters, where 0 stands for the unexpressed
%                           C0 cluster. C0 cluster will be used as control
%                           for comparison and test.
%       input_genename    : gene names of input genes corresponding to
%                           cluster indeces. 
%                           Note 1: the snv_number gene name and input_genename
%                           should have the same name type. 
%                           Note 2: input_gene_cluster and input_genename
%                           should have exactly one-to-one matching to each other.
%       NOTE:               Keep the input of 'gene_id' in "fn_snv_type_number"
%                           and "fn_base_frequency" and "input_genename"
%                           all consistent, i.e. Ensembl ID.
%  varargin:
%       'snv_asymmetry_order': SNV asymmetry calculation.
%                           default:{'AT','TA','AG','TC','TG','AC','CT','GA','GT','CA','CG','GC'}
%                           stands for: A>T/T>A; A>G/T>C; T>G/A>C; C>T/G>A; G>T/C>A; C>G/G>C;
%       'sequence_context_order': default reference context order.
%                           default:{'A','G','T','C'}
%                           With "default_snv_asymmetry" above, the 16 sequence context is
%                           displayed according to A-T-G-C order. e.g. for A>T/T>A:
%                           AAAT,AATT,AAGT,AACT,TAAT,TATT,TAGT,TACT,
%                           GAAT,GATT,GAGT,GACT,CAAT,CATT,CAGT,CACT
%       'Addup_SNV_file'  : String arry of file names of Addup_SNV.tsv file 
%                        (genes_by_SNVtypenumbers). It could have as many
%                        as possible Addup_SNV_file file names.
%                        Option for adding up additional SNV.tsv file to
%                        the input SNV.tsv files. The rows and columns
%                        should stand for exactly the same with the input
%                        file, fn_snv_type_number, meaning:
%                        first column should be gene names; other columns
%                        contains SNV number of each mutation type, including: 
%                        A>T,A>G,A>C,T>A,T>G,T>C,G>A,G>T,G>C,C>A,C>T,C>G according to 5'&3' bases.
%       'FaceColor'       : Color map used to bar-plot the asymmetry scores by 'input_gene_cluster'.
%       'MutRate_ylimit'  : Set y-limit of the mutation rates plots. 
%       'MutRate_ylimit2' : Set y-limit of the mutation rates plots.
%                           Specifically for CpG mutations, which are ~10 folds higher than others. 
%       'AsySco_ylimit'   : Set y-limit of the asymmetry score plots. 
%       'Title_prefix'    : Set title prefix of all the plots. 
%  Output:
%    OUTPUT_STRUCTURE is a output structure file which contains the below arrays/matrices:
%       coding_mut_level  : Mutation level on the coding strand by each
%                           mutation type. e.g. for ATG>C mutation level on CS:
%                           calculated as total SNV count per 1k ref sequence.
%       coding_ref_freq   : Total count of reference tribucleotides on
%                           coding strand of the genes in the corresponding gene cluster.
%       coding_mut_rate   : Mutation rate on the coding strand, normalized
%                           by per 1k reference allele.
%       template_mut_level: Mutation level on the template strand by each
%                           mutation type. e.g. for ATG>C mutation level on TS:
%                           the ATG>C mutation on TS means CAT>G on CS, 
%                           calculated as total SNV count per 1k ref sequence.
%       template_ref_freq : Total count of reference tribucleotides on
%                           template strand of the genes in the corresponding gene cluster.
%       template_mut_rate : Mutation rate on the template strand, normalized
%                           by per 1k reference allele.
%       asymmetry_score   : log2(coding_mut_rate/template_mut_rate).
%       mutation_type     : Corresponds to row names of the above 3 outputs:
%                           First row is the reference; second row is tri-base
%                           mutation type in which the middle base change to
%                           the mutation base; thrrd row is the single-base mut type;
%       gene_cluster_annot: Corresponds to column names of the top 3 outputs.
%
%  Bo Xia 05/31/2019

%% Parse defaul settings and setup inputs
fprintf('****************************************\n')
time1 = datetime;
fprintf('%s   192-mutation type analysis started.\n',datestr(time1))

%default asymmetry calculations: 
%                       A>T/T>A;A>G/T>C;T>G/A>C;C>T/G>A;G>T/C>A;C>G/G>C;
%default_snv_asymmetry = {'AT','TA','AG','TC','AC','TG','CT','GA','GT','CA','CG','GC'}; %Don't change!
default_snv_asymmetry = {'AT','TA','AG','TC','TG','AC','CT','GA','GT','CA','CG','GC'}; %Don't change!
%default_context_order = {'A','T','G','C'}; %Don't change!
default_context_order = {'A','G','T','C'}; %Don't change!
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
default_MutRate_ylimit  = [0 1000.1]; 
default_MutRate_ylimit2 = [0 1000.1]; 
default_AS_ylimit       = [-0.5 1]; 
default_title_prefix    = '';
default_addup_SNV_file  = [];

p = inputParser;
addOptional(p,'snv_asymmetry_order',default_snv_asymmetry)
addOptional(p,'sequence_context_order',default_context_order)
addOptional(p,'FaceColor',default_face_color)
addOptional(p,'MutRate_ylimit',default_MutRate_ylimit)
addOptional(p,'MutRate_ylimit2',default_MutRate_ylimit2)
addOptional(p,'AsySco_ylimit',default_AS_ylimit)
addOptional(p,'Title_prefix',default_title_prefix)
addOptional(p,'Addup_SNV_file',default_addup_SNV_file)
parse(p,varargin{:})

%default numbers
default_snv_asymmetry   = p.Results.snv_asymmetry_order;
default_context_order   = p.Results.sequence_context_order;
default_face_color      = p.Results.FaceColor;
default_MutRate_ylimit  = p.Results.MutRate_ylimit;
default_MutRate_ylimit2 = p.Results.MutRate_ylimit2;
default_AS_ylimit       = p.Results.AsySco_ylimit;
default_title_prefix    = p.Results.Title_prefix;
default_addup_SNV_file    = p.Results.Addup_SNV_file;

%setup input files
    %[SNV_number,SNV_genenames,SNV_type]
    [vnum,vgnames,vtype_names] = tsv_import(fn_snv_type_number);
    %[base_frequency,genenames,base_type]
    [base_fq,bfgnames,base_names] = tsv_import(fn_base_frequency);

    
%input sample check if need to Add-up of input SNV_type numbers.
for i = 1:length(default_addup_SNV_file)
    addup = default_addup_SNV_file(i);
    fprintf('	Add up additional SNV count file %d.\n',i)
    fprintf('       File %d name is: \n       %s\n',i,addup)
    %Addup_SNV_file1
    [vnum_a1,vgnames_a1,vtype_names_a1] = tsv_import(addup);
    [~,a,b]=intersect(vgnames,vgnames_a1,'stable');
    [~,c,d]=intersect(vtype_names,vtype_names_a1,'stable');
    if sum(a-b)==0 & sum(c-d)==0
        vnum     = vnum + vnum_a1;
        fprintf('       File %d add up successfully.\n',i)
    else
        fprintf('   !!Add up file %d has mismatch error! \n',i)
        fprintf('   Addup_SNV_file and input SNV count file hould match data organization.\n')
    end
end

mutated_genes = find(sum(vnum,2)>0);
vnum    = vnum(mutated_genes,:);
vgnames = vgnames(mutated_genes);

    %other default values
    gclu            = input_gene_cluster;       %input gene cluster
    gname           = cellstr(input_genename);  %gname and gclu should be one-to-one.
    unique_gclu     = unique(input_gene_cluster);%unique gene cluster annotation
    asym_order      = default_snv_asymmetry;    %asymmetry pairs
    context_order   = default_context_order;    %base order
    face_cm         = default_face_color;       %face color of asymmetry score plots
    mutrate_ylimit  = default_MutRate_ylimit;        %ylimit of asymmetry score plots
    as_ylimit       = default_AS_ylimit;        %ylimit of asymmetry score plots

%Report input settings
    fprintf('Input gene_clusters includes:\n')
    for i = 1:length(unique_gclu)
        fprintf('    Cluster %s\n',num2str(unique_gclu(i)))
    end
    fprintf('Default mutation order as: \n')
    for mut = 1:(length(asym_order)./2) %for every X->Y mutation type, e.g. A>G/T>C
        mut_type = char(asym_order((mut*2-1):mut*2));
        fprintf('    %s>%s/%s>%s\n',...
                 mut_type(1,1),mut_type(1,2),mut_type(2,1),mut_type(2,2));
    end
    
    
    
%% set up SNV frequency variants
%1st filter: select genes with valid base frequency and valid SNV data
[tem_genenames,i,j]=intersect(vgnames,bfgnames,'stable');
vnum     = vnum(i,:);
base_fq = base_fq(j,:);
if length(j)==0
    fprintf('Error! Invalid gene names between input SNV/base_frequency files!\n')
    fprintf('Double check and keep the "gene_id" in "fn_snv_type_number" and "fn_base_frequency" consistent!\n')
    fprintf('Function stopped!\n')
    return
end
%2nd filter: select genes overlapping with input gene cluster
[valid_genenames,i,j]=intersect(gname,tem_genenames,'stable');
gclu    = gclu(i);
vnum    = vnum(j,:);
base_fq = base_fq(j,:);
if length(j)==0
    fprintf('Error! Invalid gene names between "input_genename" and SNV/base_frequency "gene_id"!\n')
    fprintf('Double check and keep the "gene_id" and "input_genename" all consistent!\n')
    fprintf('Function stopped!\n')
    return
end

%% Compute SNV per kb per mutation type with tri-bucleotide ref, with refering to strand
%Set output files:
gene_cluster_annot = {};
mutation_type      = []; %output will be two rows of all mutation type, 1st row as ref and 2nd row as ref>alt;
coding_mut_level   = []; 
coding_ref_freq    = []; 
coding_mut_rate    = []; 
template_mut_level = []; 
template_ref_freq  = []; 
template_mut_rate  = []; 
asymmetry_score    = []; 
%Large nested for loop for doing this
for mut = 1:(length(asym_order)./2) %for every X->Y mutation type, e.g. A>G/T>C
    mut_type = char(asym_order((mut*2-1):mut*2));
    fprintf('Started processing %s>%s/%s>%s mutation pair!\n',...
             mut_type(1,1),mut_type(1,2),mut_type(2,1),mut_type(2,2));
    for three = 1:length(context_order) %for every base at 3' of X
        for five = 1:length(context_order) %for every base at 5' of X
            input_ref = [char(context_order(five)),mut_type(1,1),char(context_order(three))];
            %Append the mutation type
            mutation_type = [mutation_type [cellstr(input_ref);...
                                            cellstr(sprintf('%s>%s',input_ref,mut_type(1,2)));...
                                            cellstr(sprintf('%s>%s/%s>%s',mut_type(1,1),mut_type(1,2),mut_type(2,1),mut_type(2,2)))]]; 
            %Calculate the mutation levels of all gene cluster
            for c = 1:length(unique_gclu)
                temp1 = strcmp([input_ref,mut_type(1,2)],vtype_names);
                temp2 = strcmp(input_ref,base_names);
                coding_level(c,1) =  sum(vnum(find(gclu == unique_gclu(c)),temp1));
                coding_ref(c,1)   =  sum(base_fq(find(gclu == unique_gclu(c)),temp2));
                coding_rate(c,1)  =  coding_level(c,1) / coding_ref(c,1) *1000;
                temp3 = strcmp([seqrcomplement(input_ref),seqrcomplement(mut_type(1,2))],vtype_names);
                temp4 = strcmp(seqrcomplement(input_ref),base_names);
                template_level(c,1) =  sum(vnum(find(gclu == unique_gclu(c)),temp3));
                template_ref(c,1)   =  sum(base_fq(find(gclu == unique_gclu(c)),temp4));
                template_rate(c,1)  =  template_level(c,1) / template_ref(c,1)  *1000;
                assymscore(c,1) =  log2(coding_rate(c,1)./template_rate(c,1));
                gene_cluster_annot(c,1) = cellstr(sprintf('GeneClass%s',num2str(unique_gclu(c))));
            end
            %Append the mutation level and assymmetry score
            coding_mut_level   = [coding_mut_level coding_level]; 
            coding_ref_freq    = [coding_ref_freq coding_ref]; 
            coding_mut_rate    = [coding_mut_rate coding_rate]; 
            template_mut_level = [template_mut_level template_level]; 
            template_ref_freq  = [template_ref_freq template_ref]; 
            template_mut_rate  = [template_mut_rate template_rate]; 
            asymmetry_score    = [asymmetry_score assymscore]; 
        end
    end
end


%% VISUALIZATION

%subplot by gene cluster
n = length(unique_gclu);
m = length(mutation_type(2,:));
CS_TS_facecolor = [240 83 59;54 144 192]/256;
backcolor = [229 245 249;254 235 232;230 230 230;239 237 245;229 245 224;254 232 200]/256;
%% Figure 1,plot mutation level by type
figure;
for a=1:n
    subplot(n,1,a);
    data = [coding_mut_level(a,:);template_mut_level(a,:)]';
    data = data./sum(sum(data));
    lim = max(max(data)) * 1.1;
    %Step1: set background color
    for t = 1:6
        area([t*16-15.5 t*16+0.5],[lim lim],...
            'FaceColor',backcolor(t,:),'LineStyle','none'); hold on;
        line([t*16+0.5 t*16+0.5],[0 lim],'Color',[0.4 0.4 0.4]);
        line([t*16-7.5 t*16-7.5],[0 lim],'Color',[0.6 0.6 0.6],'LineStyle','--');
        hold on;
        if a==1
        text(t*16-14, lim*0.95, mutation_type(3,t*16),'FontSize',10);hold on;
        end
    end
    %Step2: plot informative values
    b=bar(data);hold on;
    b(1).FaceColor = CS_TS_facecolor(1,:);
    b(2).FaceColor = CS_TS_facecolor(2,:);
    %Step3: set axis properties
    ax = gca;
    ylim([0 lim]);
    ylabel(gene_cluster_annot(a)); 
    xlim([0.2 m+0.8]); xticks(1:m); 
    set(gca,'color','none'); box off;
    set(ax,'Xticklabel',[],'TickLength',[0.005, 0.01],'TickDir', 'out','fontsize',5);
    ax.LabelFontSizeMultiplier = 2 ;
    if a==n
        xticklabels(mutation_type(2,:)); xtickangle(90);
        xlabel('Mutation type with reference to sequence context');
    end
    hold off;
end
suptitle([default_title_prefix ' [Germline mutation levels (total number)]'])

%% Figure 2, plot mutation rate by type. 
%Because of the high mutation rates of  CpG>TpG, the Y-axis of the plots will be scalled.
if default_MutRate_ylimit(2)~=1000.1
    lim = mutrate_ylimit(2);
else
    lim = ceil(max(max([coding_mut_rate;template_mut_rate]))./10)*10;
end
figure;
for a=1:n
    left = 0.1;   width = 0.85;
    hei1 = 0.75*0.9/n;  hei2 = 0.22*0.9/n;
    bot1 = 1 - 0.9*a/n;  bot2 = 1 - 0.9*a/n + hei1;
    data = [coding_mut_rate(a,:);template_mut_rate(a,:)]';
    subplot('position',[left bot1 width hei1]); 
    %Step1: set background color
    for t = 1:6
        area([t*16-15.5 t*16+0.5],[lim lim],...
            'FaceColor',backcolor(t,:),'LineStyle','none'); hold on;
        line([t*16+0.5 t*16+0.5],[0 lim],'Color',[0.4 0.4 0.4]);
        line([t*16-7.5 t*16-7.5],[0 lim],'Color',[0.6 0.6 0.6],'LineStyle','--');
        hold on;
        if a==1
        text(t*16-14, lim*0.95, mutation_type(3,t*16),'FontSize',10);hold on;
        end
    end
    %Step2: plot informative values
    b=bar(data);hold on;
    b(1).FaceColor = CS_TS_facecolor(1,:);
    b(2).FaceColor = CS_TS_facecolor(2,:);
    %Step3: set axis properties
    ax = gca;
    ylim([0 lim]);
    ylabel(gene_cluster_annot(a)); 
    xlim([0.2 m+0.8]); xticks(1:m); 
    set(gca,'color','none'); box off;
    set(ax,'Xticklabel',[],'TickLength',[0.005, 0.01],'TickDir', 'out','fontsize',5);
    ax.LabelFontSizeMultiplier = 2 ;
    if a==n
        xticklabels(mutation_type(2,:)); xtickangle(90);
        xlabel('Mutation type with reference to sequence context');
    end
    hold off;
    %
    subplot('position',[left bot2 width hei2]); 
    %Step1: set background color
    lim2 = default_MutRate_ylimit2;
    if lim2(2) == 1000.1
        lim2(2) = ceil(max(max([coding_mut_rate;template_mut_rate]))./10)*10;
        lim2(1) = 0.5*lim2(2);
    end
    for t = 1:6
        area([t*16-15.5 t*16+0.5],[lim2(2) lim2(2)],...
            'FaceColor',backcolor(t,:),'LineStyle','none'); hold on;
        line([t*16+0.5 t*16+0.5],[0 lim2(2)],'Color',[0.4 0.4 0.4]);
        line([t*16-7.5 t*16-7.5],[0 lim2(2)],'Color',[0.6 0.6 0.6],'LineStyle','--');
        hold on;
    end
    %Step2: plot informative values
    b=bar(data);hold on;
    b(1).FaceColor = CS_TS_facecolor(1,:);
    b(2).FaceColor = CS_TS_facecolor(2,:);
    %Step3: set axis properties
    ax = gca;
    ylim(lim2);
    xlim([0.2 m+0.8]); 
    set(gca,'color','none'); box off;
    set(ax,'XColor','none','xtick',[],'Xticklabel',[],'TickLength',[0.005, 0.01],'TickDir', 'out','fontsize',5);
    
end
suptitle([default_title_prefix ' [Germline mutation rates (per kb)]'])

%% Figure 3: plot asymmetry scores
figure;
for a=1:n
    subplot(n,1,a); 
    lim = as_ylimit(2);
    %Step1: set background color
    for t = 1:6
        area([t*16-15.5 t*16+0.5],[lim lim],...
            'FaceColor',backcolor(t,:),'LineStyle','none'); hold on;
        area([t*16-15.5 t*16+0.5],[as_ylimit(1) as_ylimit(1)],...
            'FaceColor',backcolor(t,:),'LineStyle','none'); hold on;
        line([t*16+0.5 t*16+0.5],[as_ylimit(1) lim],'Color',[0.4 0.4 0.4]);
        line([t*16-7.5 t*16-7.5],[as_ylimit(1) lim],'Color',[0.6 0.6 0.6],'LineStyle','--');
        hold on;
        if a==1
        text(t*16-14, lim*0.95, mutation_type(3,t*16),'FontSize',10);hold on;
        end
    end
    %Step2: plot asymmetry scores
    line([0.2 m+0.8],[0 0],'Color',[0 0 0]); hold on;
    asy  = asymmetry_score(a,:)';
    asy2 = asy;
    cs_0 = find(isinf(asy) & (asy<0));
    ts_0 = find(isinf(asy) & (asy>0));
    asy2(cs_0) = 0;   %Coding strand mutation rate = 0;
    asy2(ts_0) = 0;   %Template strand mutation rate = 0;
    bar(asy2,0.7,'FaceColor',face_cm(a,:));hold on;
    %Step2.5: Label missing values in the AS plots
    for z = 1:length(cs_0)
        text(cs_0(z), 0, 'CS=0','FontSize',4.5,'Rotation',90);hold on;
    end
    for z = 1:length(ts_0)
        text(ts_0(z), 0, 'TS=0','FontSize',4.5,'Rotation',90);hold on;
    end
    %Step3: set axis properties
    ax = gca;
    ylim(as_ylimit);
    ylabel(gene_cluster_annot(a)); 
    xlim([0.2 m+0.8]); xticks(1:m); 
    set(gca,'color','none'); box off;
    set(ax,'Xticklabel',[],'TickLength',[0.005, 0.01],'TickDir', 'out','fontsize',5);
    ax.LabelFontSizeMultiplier = 2 ;
    if a==n
        xticklabels(mutation_type(2,:)); xtickangle(90);
        xlabel('Mutation type with reference to sequence context');
    end
    hold off;
end
suptitle([default_title_prefix ' [Asymmetry scores]'])

time2=datetime;
fprintf('%s     Finished the second part.\n',datestr(time2))
fprintf('Check out the plots!\n ')
fprintf('The analysis took: ')
time2-time1
fprintf('****************************************\n')


%% outputs as a MATLAB structure
%Set output files :
OUTPUT_STRUCTURE.gene_cluster_annot = gene_cluster_annot;
OUTPUT_STRUCTURE.mutation_type      = mutation_type; %three rows: ref, tri-base mutation type, single-base mut types;
OUTPUT_STRUCTURE.coding_mut_level   = coding_mut_level; 
OUTPUT_STRUCTURE.coding_ref_freq    = coding_ref_freq; 
OUTPUT_STRUCTURE.coding_mut_rate    = coding_mut_rate; 
OUTPUT_STRUCTURE.template_mut_level = template_mut_level; 
OUTPUT_STRUCTURE.template_ref_freq  = template_ref_freq; 
OUTPUT_STRUCTURE.template_mut_rate  = template_mut_rate; 
OUTPUT_STRUCTURE.asymmetry_score    = asymmetry_score;


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

