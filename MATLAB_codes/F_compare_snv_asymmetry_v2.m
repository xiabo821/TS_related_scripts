function OUTPUT_STRUCTURE = F_compare_snv_asymmetry_v2(...
                                    fn_snv_type_number,fn_base_frequency,...
                                    input_gene_cluster,input_genename,...
                                    varargin)
%% OUTPUT_STRUCTURE = F_compare_snv_asymmetry_v2(...
%                                             fn_snv_number,...
%                                             fn_base_frequency,...
%                                             input_gene_cluster,...
%                                             input_genename,...
%                                             varargin);
%  F_compare_snv_asymmetry_v2 is wrote for 'transcriptional-scanning' manuscript.
%  Overall, this function use [SNV frequency] as input, [gene_cluster] as
%  grouping/staging method, [reference_group] as the reference for
%  comparison. Then this function plots the results of SNV frequency
%  between coding strand and template strand, as well as asymmetry scores,
%  and generates ranksum test p_values with Bonferroni correction.
%  
%  ############################################################################
%  NOTE1: THE MAJOR DIFFERENCE OF this V2 (compared with v1, F_compare_snv_asymmetry) 
%         is that V2 can select the Mode of calculating the mutation. 
%         For details, check out varargin: 'Mode'.
%  NOTE2: v2 also allows option to plot the results by adding genes up in
%         each gene cluster, generating a metagene of all genes.
%  NOTE3: F_compare_snv_asymmetry_v2 also allows analysing gene clusters
%         which do not have an unexpressed-gene cluster(C0). In this case,
%         the function with just do the plotting without statistical comparison.
%  NOTE4: F_compare_snv_asymmetry_v2 also allows option to plot the results
%         by box plot all together, without considering the mutation types.
%         To do this, it requires the F_compare_variants.m function file 
%         in the same folder. 
%  2019-07-11 by Bo Xia
%  ############################################################################
%
%  Required functions: 
%           barwitherr.m            Copyright@ Martina F. Callaghan
%           F_compare_variants.m    Copyright@ Bo Xia
%  Input:
%    fn_snv_type_number: file name of .tsv file (genes_by_SNVtypenumbers).   
%                        First column should be gene names. Other rows
%                        contains SNV number of each mutation type, including: 
%                        A>T,A>G,A>C,T>A,T>G,T>C,G>A,G>T,G>C,C>A,C>T,C>G
%    fn_base_frequency : file name of .tsv file which contains A/T/G/C 
%                        numbers of a given gene/sequence. First column 
%                        should be genenames. Note: base_frequency gene name  
%                        and snv_number gene names should be the same type. 
%    input_gene_cluster: gene cluster indeces. The cluster should be 0-n
%                        integeters, where 0 stands for the unexpressed
%                        C0 cluster. If there is a C0 cluster, then this gene 
%                        cluster will be used as control for comparison
%                        and statistical test. If not, then no comparison
%                        will be performed or indicated in the plots.
%    input_genename    : gene names of input genes corresponding to
%                        cluster indeces. 
%                        Note 1: the snv_number gene name and input_genename
%                        should have the same name type, e.g. Ensembl IDs. 
%                        Note 2: input_gene_cluster and input_genename
%                        should have exactly one-to-one matching to each other.
%  varargin:
%    'snv_asymmetry_order': SNV asymmetry calculation.
%          default:{'AT','TA','AG','TC','TG','AC','CT','GA','GT','CA','CG','GC'}
%          stands for: A>T/T>A; A>G/T>C; T>G/A>C; C>T/G>A; G>T/C>A; C>G/G>C;
%    'Mode'            : Option for analysis mode, meaning calculating mutation
%                        rates by single genes (default,'single_gene'); Or,
%                        by metagenes (option 'meta_INTEGER',INTEGER stands 
%                        for an integer, e.g.'meta_20', means merging mutation 
%                        stats and reference frequencies of meta genes for 
%                        every random 20 genes); while 'meta_1' will be the same 
%                        with 'single_gene' mode since it works on
%                        individual genes.
%                        Default: 'single_gene', same with 'meta_1';
%                        Option: 'meta_10' or any other 'meta_INTEGER'.
%    'Include_meta_all': Option for including a 'meta_all' results, meaning
%                        adding up mutation/base_frequency of each gene cluster 
%                        all togeter as one meta gene. Then plot the
%                        mutation rate and asymmetry score results without
%                        error bars/confidence intervals.
%                        Default: 'no';
%                        Option: 'yes'; No other terms allowed.
%    'F_compare_variants': Option for including a 'F_compare_variants' results. 
%                        This function allows boxplotting the mutation rates
%                        across gene clusters. Check F_compare_variants.m
%                        for details. 
%                        Default: 'no';
%                        Option: 'yes'; Will include the analysis by
%                        comparing the variants level in total.
%                        Option: 'only'; Will only compare and plot the
%                        variants level in total without generating further
%                        downstream analysis into asymmetries, nor will
%                        generate output stats of the data.
%    'Norm_length'     : Option for setting the normalization length.
%                        Default: 1000, meaning normalize to 1000 bases.
%    'Addup_SNV_file'  : String arry of file names of Addup_SNV.tsv file 
%                        (genes_by_SNVtypenumbers). It could have as many
%                        as possible Addup_SNV_file file names.
%                        Option for adding up additional SNV.tsv file to
%                        the input SNV.tsv files. The rows and columns
%                        should stand for exactly the same with the input
%                        file, fn_snv_type_number, meaning:
%                        first column should be gene names; other columns
%                        contains SNV number of each mutation type, including: 
%                        A>T,A>G,A>C,T>A,T>G,T>C,G>A,G>T,G>C,C>A,C>T,C>G.
%    'Subtraction'     : Option for involving substraction for the
%                        fn_snv_type_number input.
%                        Default is 'no'. 
%                        Note: if option 'yes', do need to assign
%                              fn_snv_type_number_ref and
%                              fn_base_frequency_ref, which should have
%                              exactly the same data organization with the
%                              fn_snv_type_number file and
%                              fn_base_frequency  file, respectively.
%    'fn_snv_type_number_ref': The optional input file used as reference
%                        to subtract from the the fn_snv_type_number input. 
%                        Note this file should to have exactly the same data 
%                        organization with the fn_snv_type_number file.
%    'fn_base_frequency_ref' : The optional input file used as reference
%                        to subtract from the fn_base_frequency input.
%                        Note this file should to have exactly the same data
%                        organization with the fn_base_frequency file.
%    'FaceColor'       : color map used to bar-plot the asymmetry scores.
%    'Ylim_fold_factor': Multiply the ylim of mutation rate plots so that to
%                        adjust the ylim for visualization. Default is 1.
%    'Title_prefix'    : Set title prefix of all the plots. 
%  
%  Output: OUTPUT_STRUCTURE, containing the below sub_variables:
%    SNVtype_by_column : The column names for each SNVtype corresponding
%                        to the values listed together. 
%                        Same with the input SNV order. 
%    gene_cluster_name : The gene clusters used in the analysis. 
%                        Same with the input gene cluster names.
%    CSvsTS_test_p     : paired sample t-test p_values for each SNVtype and each
%                        gene cluster. CS_TS_test_p < 0.0001 is also marked
%                        as '*' in the strand_SNV_level barwitherr plots.
%    expVSunexp_asym_score_test_p: ranksum test p_values of asymmetrys between
%                        expressed gene clusters vs unexpressed gene
%                        cluster. Asymmetry_test_p < 0.001 is also marked
%                        as '*' in the asymmetry_score barwitherr plots.
%    strand_asym_stat  : A x-y-z matrix. Average mutation rates (Z1) and  
%                        confidence intervals (Z2 & Z3) for specific gene 
%                        clusters (x-axis) distinguished by SNV type (y-axis). 
%    asym_score_stat   : A x-y-z matrix. Average asymmetry scores (Z1) and 
%                        their confidence intervals (Z2 & Z3) for specific 
%                        gene clusters (x-axis) distinguished by SNV type 
%                        pairs (y-axis). 
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
default_norm_length  = 1000;
default_F_compare_variants = 'no';
defaultSubtraction = 'no';
default_Include_meta_all = 'no';
validOption = {'no','yes','only'};
checkOption = @(x) any(validatestring(x,validOption));

defaultfn_snv_type_number_ref = '';
defaultfn_base_frequency_ref  = '';
default_title_prefix          = '';
defaulty_lim_fold_factor       =1;
default_addup_SNV_file  = [];

p = inputParser;
addOptional(p,'snv_asymmetry_order',default_snv_asymmetry)
addOptional(p,'Mode',default_Mode)
addOptional(p,'Norm_length',default_norm_length,@isnumeric)
addOptional(p,'F_compare_variants',default_F_compare_variants,checkOption)
addOptional(p,'Subtraction',defaultSubtraction,checkOption)
addOptional(p,'Include_meta_all',default_Include_meta_all,checkOption)
addOptional(p,'fn_snv_type_number_ref',defaultfn_snv_type_number_ref)
addOptional(p,'fn_base_frequency_ref',defaultfn_base_frequency_ref)
addOptional(p,'FaceColor',default_face_color)
addOptional(p,'Ylim_fold_factor',defaulty_lim_fold_factor,@isnumeric)
addOptional(p,'Title_prefix',default_title_prefix)
addOptional(p,'Addup_SNV_file',default_addup_SNV_file)
parse(p,varargin{:})

%default numbers
default_snv_asymmetry = p.Results.snv_asymmetry_order;
default_Mode = p.Results.Mode;
default_norm_length = p.Results.Norm_length;
default_Include_meta_all = p.Results.Include_meta_all;
default_F_compare_variants = p.Results.F_compare_variants;
defaultSubtraction = p.Results.Subtraction;
defaultfn_snv_type_number_ref = p.Results.fn_snv_type_number_ref;
defaultfn_base_frequency_ref = p.Results.fn_base_frequency_ref;
default_face_color = p.Results.FaceColor;
defaulty_lim_fold_factor = p.Results.Ylim_fold_factor;
default_title_prefix    = p.Results.Title_prefix;
default_addup_SNV_file    = p.Results.Addup_SNV_file;


%% setup input files
%[SNV_number,SNV_genenames,SNV_type]
[vnum,vgnames,vtype_names] = tsv_import(fn_snv_type_number);
%[base_frequency,genenames,base_type]
[vglength,vglenames,base_names] = tsv_import(fn_base_frequency);

%input sample check if involves subtraction of input SNV_type numbers.
if strcmp('yes',cellstr(defaultSubtraction))
    if (strcmp('',cellstr(defaultfn_snv_type_number_ref)) || ...
        strcmp('',cellstr(defaultfn_base_frequency_ref)))
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
            fprintf('!!Inputs mismatch error! \n')
            fprintf('fn_snv_type_number and fn_snv_type_number_ref should match data organization.\n')
            return
        end
        if sum(e-f)==0 & sum(g-h)==0
            vglength = vglength - vglength_ref;
        else
            fprintf('!!Inputs mismatch error! \n')
            fprintf('fn_base_frequency and fn_base_frequency_ref should match data organization.\n')
            return
        end
    end
end

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

%Filter genes by selecting the genes with mutations.
mutated_genes = find(sum(vnum,2)>0);
vnum    = vnum(mutated_genes,:);
vgnames = vgnames(mutated_genes);

%other default values
gclu   = input_gene_cluster;
gname  = cellstr(input_genename);
asym_order  = default_snv_asymmetry;
face_cm     = default_face_color;



%% Set up SNV frequency variants and input gene list

[v_genenames,i,j]=intersect(vgnames,vglenames,'stable');
vgnames  = vgnames(i);
vnum     = vnum(i,:);
vglenames= vglenames(j);
vglength = vglength(j,:);
if length(j)==0
    fprintf('Error! Mismatched gene names in the input SNV/base_frequency files!')
    return
end
%setup the input gene list
[valid_geneid,a,b]=intersect(gname,v_genenames,'stable');
%check gene names
if length(a)==0
    fprintf('  !Error! Mismatched gene names in <input_genename> \n');
    fprintf('              and in the <SNV/base_frequency files>!\n');
    return
else
    fprintf('  %d out of %d input_genes are included for downstream analysis!\n',...
               length(valid_geneid), length(gname)  );
end
valid_geneid = valid_geneid;
gclu         = gclu(a);
gname        = gname(a);
vnum         = vnum(b,:);
vglength     = vglength(b,:);
v_genenames  = v_genenames(b);
if min(gclu)~=0
    fprintf('  !Attention! input_gene_cluster does not have C0!\n');
    fprintf('      Indicating that no unexpressed gene cluster!\n');
    fprintf('      There will not be statistical comparison between gene clusters!\n');
end



%% Calculate SNV frequency depending on Mode. 

if strcmp('single_gene',cellstr(default_Mode))
    %calculate SNV frequency per kilo-bases for each gene.
    fprintf('NOTE: you are performing mutation asymmetry analysis based on individual genes.\n')
    default_Mode = 'meta_1'; meta_n = 1;
    vnum_perkb = [];
    for i = 1:length(asym_order)
        m = char(asym_order(i));
        ref = cellstr(m(1));
        vnum_perkb(:,strcmp(asym_order(i),asym_order)) = vnum(:,strcmp(asym_order(i),vtype_names))...
                                    ./   vglength(:,strcmp(ref,base_names)) * default_norm_length;
    end
    %remove NaN and Inf numbers, which indicates no ref sequence stats from the data
    valid =  ~any( isnan( vnum_perkb ) | isinf( vnum_perkb ), 2 );
    vnum_perkb   = vnum_perkb(valid,:);
    valid_geneid = valid_geneid(valid);
    vnum         = vnum(valid,:);
    vglength     = vglength(valid,:);
    v_genenames  = v_genenames(valid);
    gclu         = gclu(valid);
    gname        = gname(valid);
    
elseif strncmp('meta',cellstr(default_Mode),4)
    % default_Mode = 'meta_10'
    meta_n = split(default_Mode,'_');
    meta_n = str2num(char(meta_n(2)));
    if meta_n>0
        fprintf('NOTE: Performing mutation asymmetry analysis based on meta_genes.\n')
        fprintf('    Meta_genes: average stats of every %d gene(s).\n',meta_n)
        if meta_n == 1
        fprintf('    Same with <single_gene> mode, or default mode.\n')
        end
    else
        fprintf('ERROR: Meta_gene-based mutation asymmetry analysis has problems!')
        fprintf('     Check you input of <Mode>. Should be like, e.g.: <meta_10> \n')
        return
    end
        %reassign gene cluster: cluster by cluster in a for loop:
    gclu_type = unique(gclu)'; 
    glu_meta      = [];
    vnum_meta     = []; 
    vglength_meta = [];
    for j = 1:length(gclu_type)
        gc = find(gclu==gclu_type(j));
        valid_geneid_temp = valid_geneid(gc);
        gclu_temp         = gclu(gc);
        vnum_temp         = vnum(gc,:);
        vglength_temp     = vglength(gc,:);
        
        %for each cluster, the genes are compiled for every meta_n genes
        meta_gene_n = ceil(length(gc)./meta_n);
        gclu_temp2 =[]; vnum_temp2 =[]; vglength_temp2 =[];
        for k = 1:(meta_gene_n-1)
            %calculate mean on columns, but not on rows
            gclu_temp2(k,:)     = mean(gclu_temp((k-1)*meta_n+1:k*meta_n , :),1); 
            vnum_temp2(k,:)     = mean(vnum_temp((k-1)*meta_n+1:k*meta_n , :),1);
            vglength_temp2(k,:) = mean(vglength_temp((k-1)*meta_n+1:k*meta_n , :),1);
        end
            k=meta_gene_n;
            gclu_temp2(k,:)     = mean(gclu_temp((k-1)*meta_n+1:end , :),1);
            vnum_temp2(k,:)     = mean(vnum_temp((k-1)*meta_n+1:end , :),1);
            vglength_temp2(k,:) = mean(vglength_temp((k-1)*meta_n+1:end , :),1);
        glu_meta      = [glu_meta; gclu_temp2];
        vnum_meta     = [vnum_meta; vnum_temp2];
        vglength_meta = [vglength_meta; vglength_temp2];
    end
    %calculate SNV frequency per kilo-bases for each meta gene.
    vnum_perkb = [];
    gclu     = glu_meta;
    vnum     = vnum_meta;
    vglength = vglength_meta;
    for i = 1:length(asym_order)
        m   = char(asym_order(i));
        ref = cellstr(m(1));
        vnum_perkb(:,strcmp(asym_order(i),asym_order)) = ...
                                vnum(:,strcmp(asym_order(i),vtype_names))...
                           ./   vglength(:,strcmp(ref,base_names)) ...
                           *    default_norm_length;
    end
end




%% F_compare_variants option: boxplotting results by gene cluster

base = {'A','C','G','T'};
if strcmp('yes',cellstr(default_F_compare_variants))
    fprintf('NOTE: you choose to INCLUDE comparing the mutation level across clusters.\n')
    [~,~,b]=intersect(base,base_names);
    variants_genelength = sum(vglength(:,b),2);
    
    variants_level = F_compare_variants(sum(vnum,2),variants_genelength,v_genenames,...
                                        gclu,gname,face_cm);
    title([default_title_prefix ' [mutation rates (per kb)]'])
    
    OUTPUT_STRUCTURE.variants_level_structure           = variants_level;
end

if strcmp('only',cellstr(default_F_compare_variants))
    fprintf('NOTE: you choose to ONLY compare the mutation level across clusters.\n')
    fprintf('      The function will stop after performing this analysis.\n')
    [~,~,b]=intersect(base,base_names);
    variants_genelength = sum(vglength(:,b),2);
    variants_level = F_compare_variants(sum(vnum,2),variants_genelength,v_genenames,...
                                        gclu,gname,face_cm);
    title([default_title_prefix ' [mutation rates (per kb)]']);
    OUTPUT_STRUCTURE.variants_level_structure           = variants_level;
    return
end


%% Compute SNV per kb per gene strand and paired_sample t-test between coding and template strands

%Calculate strand asymmetry test p_values (use paired-sample t-test with Bonferroni adjustment)
for i = 1:(length(asym_order)./2)
    mut_type = char(asym_order((i*2-1):i*2));
    fprintf('Started processing %s>%s/%s>%s mutation pair!\n',...
             mut_type(1,1),mut_type(1,2),mut_type(2,1),mut_type(2,2));
    asym_type(i) = cellstr(sprintf('%s>%s/%s>%s',mut_type(1,1),mut_type(1,2),...
                                                 mut_type(2,1),mut_type(2,2)));
    asym_type_2(i) = strcat(asym_order(2*i-1),{'_'},asym_order(2*i));
    gclu_type = unique(gclu)'; 
    for j = 1:length(gclu_type)
        gc = find(gclu==gclu_type(j));
        %strand_asym_ranksum_p(j+1,i) = ranksum(vnum_perkb(gc,2*i-1),vnum_perkb(gc,2*i)) ...
                                       %* length(unique(gclu));
        cs = vnum_perkb(gc,2*i-1);
        ts = vnum_perkb(gc,2*i);
        [~,cs_out,~,~,~] = filloutliers(cs,'previous');
        [~,ts_out,~,~,~] = filloutliers(ts,'previous');
        quolif_gene = setdiff(1:length(cs),union(cs_out,ts_out));
        %remove outlier for doing paired-sample ttest.
        [~,strand_asym_p_value(j,i)] = ttest(cs(quolif_gene),ts(quolif_gene),'Alpha',0.01);
        %[~,strand_asym_p_value(j,i)] = ttest(vnum_perkb(gc,2*i-1),vnum_perkb(gc,2*i),'Alpha',0.01);
        strand_asym_p_value(j,i) = strand_asym_p_value(j,i) * length(unique(gclu));
        %multiply by number of tests, Bonferroni adjustment
        gc_name(j) = strcat({'C'},num2str(gclu_type(j)));    
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
    gclu_type = unique(gclu)'; 
    for j = 1:length(gclu_type)
        vgc = vnum_perkb(find(gclu==gclu_type(j)),i);
        vgc = vgc(~any( isnan(vgc) | isinf(vgc), 2 ));
%         strand_asym_stat(j,i,1) = mean(vgc);
%         ci = bootci(10000,{@mean,vgc},'alpha',0.01)-mean(vgc);  %mean of strand_SNV
        ci = bootci(10000,{@mean,vgc},'alpha',0.01);
        strand_asym_stat(j,i,1) = mean(ci);
        ci = ci - mean(ci);  %mean of strand_SNV
        strand_asym_stat(j,i,2) = ci(1);  %lower CI-99% of strand_SNV
        strand_asym_stat(j,i,3) = ci(2);  %upper CI-99% of strand_SNV
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
    lim_y = line_y*1.25 * defaulty_lim_fold_factor;
    ylim([0 lim_y]);
    line([0.5 n+0.5],[line_y line_y],'Color',[227 26 28]/256,'LineStyle','--');
    text(1:n,repelem(lim_y,n),SNV_strand_asym_marker(:,a));hold off;
    title(asym_type(a));box off;
    xticklabels(gc_name);
    set(gca,'color','none');
end
suptitle([default_title_prefix '(meta: ' num2str(meta_n) ' gene)' ' [Germline mutation rates (per kb)]'])

time2 = datetime;
fprintf('%s     Finished the first part.\n',datestr(time2))


%% asymmetry score calculations, ranksum test, and VISUALIZATION
%Calculate asymmetry scores
for i = 1:(length(asym_order)./2)
    %asymmetry_score(:,i) = log2((vnum_perkb(:,2*i-1)+1) ./ (vnum_perkb(:,2*i)+1) );
    asymmetry_score(:,i) = log2((vnum_perkb(:,2*i-1)) ./ (vnum_perkb(:,2*i)) );
end

%Calculate strand asymmetry ranksum test p_values 
asym_score_ranksum_test_p = [];
for i = 1: (length(asym_order)./2)
    vgc0 = asymmetry_score(find(gclu==0),i);
    gclu_type = unique(gclu)';
        for j = 1:length(gclu_type)
            vgc = asymmetry_score(find(gclu==gclu_type(j)),i);
            vgc = vgc(~any( isnan(vgc) | isinf(vgc) | vgc==0, 2 ));
            % ranksum test of asymmetry scores versus unexpressed gene cluster, if exist.
            if min(gclu_type) == 0
            asym_score_ranksum_p(j,i) = ranksum(vgc0,vgc)...
                                          * (length(asym_order)./2 - 1);
                                          %multiply by number of tests, Bonferroni adjustment
            end
%             asym_score_stat(j,i,1) = mean(vgc);
%             ci = bootci(10000,{@mean,vgc},'alpha',0.01)-mean(vgc);  %mean of asymetry_score
            ci = bootci(10000,{@mean,vgc},'alpha',0.01);  %mean of asymetry_score
            asym_score_stat(j,i,1) = mean(ci);
            ci = ci - mean(ci);  %mean of asymetry_score
            asym_score_stat(j,i,2) = ci(1);  %lower CI-99% of asymetry_score
            asym_score_stat(j,i,3) = ci(2);  %upper CI-99% of asymetry_score
        end
end

    %export the ranksum test pvalues comparing asymmetry scores bwteen
    %expressed and unexpressed gene clusters, if unexpressed gene cluster exists.
    if min(gclu_type) == 0
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
    end
    
     
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
    if min(gclu_type) == 0
    text(1:n,repelem(max(max(asym_score_stat(:,:,1)))*1.05,n),asym_score_marker(:,a));hold off;
    end
    title(asym_type(a));box off;
    xticks(1:n);   xticklabels(gc_name);
    set(gca,'color','none');
end
suptitle([default_title_prefix '(meta: ' num2str(meta_n) ' gene)' ' [Asymmetry scores]'])

time3=datetime;
fprintf('%s     Finished the second part.\n',datestr(time3))



%% Check whether to plot meta_all results.
if strcmp('yes',cellstr(default_Include_meta_all))
    fprintf('NOTE: Performing mutation asymmetry analysis of meta_all mode.\n')
    fprintf('      Meaning averaging stats of all genes in a given cluster.\n')
    %reassign gene cluster: cluster by cluster in a for loop:
    gclu_type = unique(gclu)';
    glu_meta_all      = [];
    vnum_meta_all     = [];
    vglength_meta_all = [];
    for j = 1:length(gclu_type)
        %for each cluster, average the vnum and vglength
        %calculate mean on columns, but not on rows
        gc = find(gclu==gclu_type(j));
        gclu_temp         = mean(gclu(gc),1);
        vnum_temp         = mean(vnum(gc,:),1);
        vglength_temp     = mean(vglength(gc,:),1);
        
        glu_meta_all      = [glu_meta_all; gclu_temp];
        vnum_meta_all     = [vnum_meta_all; vnum_temp];
        vglength_meta_all = [vglength_meta_all; vglength_temp];
    end
    %calculate SNV rate per gene cluster and plot.
    vnum_perkb_all = [];
    for i = 1:length(asym_order)
        m   = char(asym_order(i));
        ref = cellstr(m(1));
        vnum_perkb_all(:,strcmp(asym_order(i),asym_order)) = ...
                      vnum_meta_all(:,strcmp(asym_order(i),vtype_names))...
                 ./   vglength_meta_all(:,strcmp(ref,base_names)) ...
                 *    default_norm_length;
    end
    %plot with bar
    l = length(asym_order)/2; %plot numbers
    n = length(unique(gclu)); %gene cluster numbers
    Z =  [1 2;3 4;5 6;7 8;9 10;11 12];
    figure;
    for a=1:l
        subplot(1,l,a);
        barplot = bar(vnum_perkb_all(:,Z(a,:)),0.9);
        barplot(1).FaceColor = [244 165 130]/256;
        barplot(2).FaceColor = [146 197 222]/256;
        xlim([0.5 n+0.5]); hold on;
        line_y = mean(vnum_perkb_all(1,Z(a,:)));
        lim_y = line_y*1.25 * defaulty_lim_fold_factor;
        ylim([0 lim_y]);
        line([0.5 n+0.5],[line_y line_y],'Color',[227 26 28]/256,'LineStyle','--');hold off;
        title(asym_type(a));box off;
        xticklabels(gc_name);
        set(gca,'color','none');
    end
    suptitle([default_title_prefix '(meta: all genes)' ' [Germline mutation rates (per kb)]'])
    
    % calculate CS_vs_TS asymmetry scores per gene cluster and plot.
    for i = 1:(length(asym_order)./2)
        %asymmetry_score_all(:,i) = log2((vnum_perkb_all(:,2*i-1)+1) ./ (vnum_perkb_all(:,2*i)+1) );
        asymmetry_score_all(:,i) = log2((vnum_perkb_all(:,2*i-1)) ./ (vnum_perkb_all(:,2*i)) );
    end
    %plot with bar
    l = length(asym_order)/2; %plot numbers
    n = length(unique(gclu)); %gene cluster numbers
    ylimit(1) = min(min(asymmetry_score_all(:,:)))*1.1;
    ylimit(2) = max(max(asymmetry_score_all(:,:)))*1.1;
    % face_cm = [[0.50,0.50,0.50;0.40,0.76,0.647;0.988,0.553,0.384;...
    %             0.5529,0.62745,0.796;0.9058,0.5412,0.7647;0.651,0.847,0.329]];
    figure;
    for a=1:l
        subplot(1,l,a);
        for i = 1:n
            bar(i,asymmetry_score_all(i,a),0.8,'FaceColor',face_cm(i,:));hold on;
        end
        hold on;
        xlim([0.5 n+0.5]);
        ylim(ylimit);
        title(asym_type(a));box off;
        xticks(1:n);   xticklabels(gc_name);
        set(gca,'color','none');
    end
    suptitle([default_title_prefix '(meta: all genes)' ' [Asymmetry scores]'])
    
    fprintf('%s     Finished the meta_all analysis part.\n',datestr(time3))
end

%% outputs
fprintf('Generating the output stats:\n')

OUTPUT_STRUCTURE.SNVtype_by_column           = default_snv_asymmetry;
OUTPUT_STRUCTURE.SNV_asym_type               = asym_type;
OUTPUT_STRUCTURE.gene_cluster_name           = gc_name;
OUTPUT_STRUCTURE.MutRate_CSvsTS_stat         = strand_asym_stat;
OUTPUT_STRUCTURE.MutRate_CSvsTS_test_p       = SNV_strand_asym_p_value;
OUTPUT_STRUCTURE.AsymScore_expVSunexp_stat   = asym_score_stat;
OUTPUT_STRUCTURE.AsymScore_expVSunexp_test_p = asym_score_ranksum_test_p;

if strcmp('yes',cellstr(default_Include_meta_all))
OUTPUT_STRUCTURE.MutRate_CSvsTS_metaall_stat = vnum_perkb_all;
OUTPUT_STRUCTURE.AsymScore_metaall           = asymmetry_score_all;
end

time4=datetime;
fprintf('Check out the plots!\n ')
fprintf('The analysis took: ')
time4-time1
fprintf('Done!\n')
fprintf('****************************************\n')

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