
%gene sets clusters data
%includes the expression level of human genes from the former step.
load('YOUR_PATH\h_m_stage_expression.mat');
% need 'h_gene_cluster','human_genename','m_gene_cluster','mouse_genename'
    human_gene_info = struct2table(tdfread('YOUR_PATH\human_genes.tsv'));
    [~,~,xi] = intersect(human_genename,cellstr(human_gene_info.name),'stable');
    h_gname  = cellstr(human_gene_info.ensembl_id(xi,:));

    mouse_gene_info = struct2table(tdfread('YOUR_PATH\mouse_genes.tsv'));
    [~,~,xi] = intersect(mouse_genename,cellstr(mouse_gene_info.name),'stable');
    m_gname  = cellstr(mouse_gene_info.ensembl_id(xi,:));



cd 'YOUR_PATH\cbrewer'
%cbrewer()
    expgene_cm = [29,180,170]/256;
    genecluster_edgecm = [0.15 0.15 0.15;cbrewer('qual','Dark2',5)];
    genecluster_facecm = [0.5 0.5 0.5;cbrewer('qual','Set2',5)];
    genecluster_facecm_ctl = [0.7 0.7 0.7;cbrewer('qual','Pastel2',5)];
    RdBu_cm = flipud(cbrewer('div','RdBu',50));
    Blues_cm = cbrewer('seq','Blues',50);
    Blues_cm = flipud(cbrewer('seq','Blues',50));
    RdBu_cm = RdBu_cm(5:45,:);
    PuBuGn_cm = cbrewer('seq','PuBuGn',50);
    PuBuGn_cm = cbrewer('seq','YlGnBu',50);
    exunex_cm = [0.5 0.5 0.5;0.25 0.68 0.46];
    paired_cm = cbrewer('qual','Paired',6);
    paired_cm = flipud(paired_cm(3:6,:));
    soma_cm = [0.15 0.15 0.15;0.3 0.3 0.3; 0.14 0.54 0.27; 0.17 0.44 0.71];
    GTEX_cm = cbrewer('qual','Dark2',53);
cd 'YOUR_PATH'



    
%% 201907 Human mutation asymmetry revisit. Done! 
% human_dbSNP_variants_in_total
%human gene body
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
     human_genebody = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_gene_cluster,...
                                     h_gname,...
                                     'F_compare_variants','yes',...
                                     'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Human-genebody'); 
%human gene body(1000GenomeSNV)
    fn1='YOUR_PATH\human_1000G_genebody_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
     human_genebody_1000G = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_gene_cluster,...
                                     h_gname,...
                                     'F_compare_variants','only',...
                                     'FaceColor',genecluster_edgecm,...
                                     'Title_prefix','Human-genebody-1000G'); 
%human intron
    fn1='YOUR_PATH\human_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.tsv'
     human_intron = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_gene_cluster,...
                                     h_gname,...
                                      'F_compare_variants','yes',...
                                    'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Human-intron'); 
%human upstream
    fn1='YOUR_PATH\human_upstream_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.upstream.singleFreq.tsv'
     human_upstream = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_gene_cluster,...
                                     h_gname,...
                                      'F_compare_variants','yes',...
                                     'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Human-upstream'); 
%human downstream
    fn1='YOUR_PATH\human_downstream_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.downstream.singleFreq.tsv'
     human_downstream = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_gene_cluster,...
                                     h_gname,...
                                      'F_compare_variants','yes',...
                                     'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Human-downstream'); 

                                 
% errorbar plot upstream-intron-downstream asymmetry scores together
ds = human_downstream.AsymScore_expVSunexp_stat;
us = human_upstream.AsymScore_expVSunexp_stat;
intron = human_intron.AsymScore_expVSunexp_stat;
gc_name = human_intron.gene_cluster_name;
asym_type = mouse_intron.SNV_asym_type;

figure;
for a=1:6
    subplot(1,6,a);    
    line([0.5 6.5],[0 0],'Color',[0 0 0]); hold on;
    %downstream
    errorbar(1:6,ds(:,a,1),ds(:,a,2),ds(:,a,3),'-d',...
        'Color',[0.98 .57 .44],'MarkerSize',5,'LineWidth',1); hold on;
    %upstream
    errorbar(1:6,us(:,a,1),us(:,a,2),us(:,a,3),'-s',...
        'Color',[0.13 .44 .71],'MarkerSize',5,'LineWidth',1); hold on;
    %genebody
    errorbar(1:6,intron(:,a,1),intron(:,a,2),intron(:,a,3),'-o',...
        'Color',[0.79 .093 .11],'MarkerSize',5,'LineWidth',1); hold on;
    title(asym_type(a));box off;
    xlim([0.5 6.5]); 
    ylim([-0.25 0.6]);
    xticks(1:6);   xticklabels(gc_name);
    set(gca,'color','none');
end

%Done





%% Human mutation asymmetry control/sensitivity analysis: SI Figure3. DONE! 

    cd 'YOUR_PATH\cbrewer'
    genecluster_edgecm = [0.15 0.15 0.15;cbrewer('qual','Dark2',5)];
    genecluster_facecm = [0.5 0.5 0.5;cbrewer('qual','Set2',5)];
    RdBu_cm = flipud(cbrewer('div','RdBu',50));
    RdBu_cm = RdBu_cm(5:45,:);
    exunex_cm = [0.5 0.5 0.5;0.25 0.68 0.46];
    cd 'YOUR_PATH\Analysis'

if replicates_sensitivity
% cluster genes with human_TESE1 or human_TESE2 germ cells
%     %cluster genes by clusters of cells
%     load('YOUR_PATH\human_tese1and2_replicates_germcells.mat');
%     A = h_germcell_sm;     
%     % define gene clusters by replicates
%     rep1 = find(h_germcell_repidx<3);
%     rep2 = find(h_germcell_repidx>2); 
%     donor1_cellcluster = h_germcell_repidx(rep1);
%     donor2_cellcluster = h_germcell_repidx(rep2);
%     histcounts(donor1_cellcluster)
%     histcounts(donor2_cellcluster)
%     rng(1)
%     [~,~,h_cluster_index1]=F_cluster_gene(   A(:,rep1),h_germ_IDX(rep1),5,...
%                                              'expression_cutoff_n',5,...
%                                              'umi_ratio_cutoff',0,...
%                                              'include_unexpressed_gene','yes',...
%                                              'ylimit',[-2 3]);
%     rep1_cluster = h_cluster_index1;
%     rep1_cluster(find(h_cluster_index1==1))=2;
%     rep1_cluster(find(h_cluster_index1==2))=5;
%     rep1_cluster(find(h_cluster_index1==3))=1;
%     rep1_cluster(find(h_cluster_index1==4))=4;
%     rep1_cluster(find(h_cluster_index1==5))=3;
%     
%     rng(1)
%     [~,~,h_cluster_index2]=F_cluster_gene(   A(:,rep2),h_germ_IDX(rep2),5,...
%                                              'expression_cutoff_n',5,...
%                                              'umi_ratio_cutoff',0,...
%                                              'include_unexpressed_gene','yes',...
%                                              'ylimit',[-2 3]);
%     rep2_cluster = h_cluster_index2;
%     rep2_cluster(find(h_cluster_index2==1))=2;
%     rep2_cluster(find(h_cluster_index2==2))=1;
%     rep2_cluster(find(h_cluster_index2==3))=5;
%     rep2_cluster(find(h_cluster_index2==4))=4;
%     rep2_cluster(find(h_cluster_index2==5))=3;
    
    % gene clusters come from merged TESE1_TESE2 gene clutersing results
    %  rep1_cluster; rep2_cluster;
        %import gene signatures
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
        F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     rep1_cluster,...
                                     h_gname,...
                                     'F_compare_variants','only',...
                                     'FaceColor',genecluster_edgecm,...
                                     'Title_prefix','Human-genebody TESE1 clustering'); 

        F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     rep2_cluster,...
                                     h_gname,...
                                     'F_compare_variants','only',...
                                     'FaceColor',genecluster_edgecm,...
                                     'Title_prefix','Human-genebody TESE2 clustering'); 

    %  rep1_cluster and rep2_cluster consensus clusters;    
    conse_genecluster = [];
    for i=1:6
        temp=find(rep1_cluster==(i-1) & rep2_cluster==(i-1));
        conse_genecluster = [conse_genecluster;temp];
    end
    gene_cluster = rep1_cluster(conse_genecluster);
    [~,~,xi] = intersect(human_genename(conse_genecluster),cellstr(human_gene_info.name),'stable');
    h_gname_rep  = cellstr(human_gene_info.ensembl_id(xi,:));
        F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     gene_cluster,...
                                     h_gname_rep,...
                                     'F_compare_variants','only',...
                                     'FaceColor',genecluster_edgecm,...
                                     'Title_prefix','Human-genebody TESE1&2 consensus clustering'); 
        
    
end  %Done




if gene_clustering_sensitivity
% cluster genes by different parameters
load('YOUR_PATH\human_tese1and2_replicates_germcells.mat');
    A = h_germcell_sm;

    %cluster genes with k=1:10; 
    %expression_cutoff: n=5,r=0; n=5,r=0.1; n=10,r=0; n=10,r=0.1; 
    cutoff=[5 0;5 0.1; 10 0;10 0.1];
    for j=1:4
    for k=1:10
    rng(1)
    [~,~,h_gene_cluster_st]=F_cluster_gene(A,h_germ_IDX,k,...
        'expression_cutoff_n',cutoff(j,1),...
        'umi_ratio_cutoff',cutoff(j,2),...
        'include_unexpressed_gene','yes',...
        'ylimit',[-2 3]);    
         N_unexp(j) = length(find(h_gene_cluster_st==0));

    %Calculate mutation rates and ratios to the unexpressed gene group.  
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
       mr = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_gene_cluster_st,...
                                     h_gname,...
                                     'F_compare_variants','only',...
                                     'FaceColor',genecluster_edgecm,...
                                     'Title_prefix','Human-genebody TESE1 clustering'); 

    temp = mr.variants_level_structure.mean;    
        mean_var_ratios{j,k} = temp(2:end)./temp(1);
    end
    end
    
    figure;
    for j=1:4
    for k=1:10 
    subplot(4,10,10*(j-1)+k)
    mean_value = mean(mean_var_ratios{j,k});
    temp = 1.5-rand(k,1);
    scatter(temp,mean_var_ratios{j,k},10,[0.5 0.5 0.5],'filled');hold on;
    errorbar(1,mean_value,std(mean_var_ratios{j,k}),...
            'Color','k','LineWidth',1);hold on;
    line(0.5:1.5,[mean_value mean_value],'Color','k','LineWidth',1.5);hold off;
    xlim([0.3 1.7]);    ylim([0.55 1]); 
    xticks(1);xticklabels(strcat('k=',num2str(k)))
        set(gca,'color','none');box off;
        if k>1
        ax1 = gca;  
        ax1.YAxis.Visible = 'off';
        elseif k==1
        ylabel(strcat('cutoff: n=',num2str(cutoff(j,1)),'ratio=',...
               num2str(cutoff(j,2)),'N(unexp)=',num2str(N_unexp(j)))  )
        end
    end
    end

end  %Done
    
  




%% Add soma-expressed genes. Fig.2 & SI Fig.3

% load somatic cells and germ cells
    load('YOUR_PATH\human_tese1and2_replicates_soma.mat');
    load('YOUR_PATH\h_m_stage_expression.mat');
    
    soma_mat = [soma_sertoli  soma_leydig  soma_myoid  soma_macrophage soma_endothelia];
    soma_cluster(1:51) = 1;
    soma_cluster(52:143) = 2;
    soma_cluster(144:303) = 3;
    soma_cluster(304:327) = 4;
    soma_cluster(328:379) = 5;
    soma_name = {'Sertoli','Leydig','Myoid','Macrophage','Endothelia'};
    B= soma_mat;     B(find(B>0))=1;
    
    human_gene_info = struct2table(tdfread('YOUR_PATH\human_genes.tsv'));
    [~,~,xj] = intersect(human_genename,cellstr(human_gene_info.name),'stable');
    h_gname  = cellstr(human_gene_info.ensembl_id(xj,:));
    
    soma_together = 1;
    if soma_together
        %binarize the gene clustering as germexp-somaexp, germexp-somaunexp,
        %                                germunexp-somaexp, germunexp-somaunexp
        soma_gene_cluster(1:length(h_gene_cluster),1) = 0;%initialize the gene indeces.
        %    soma_gene_cluster(find(max_soma>=0.1 & sum(B')>10))=1;
        min_cell_cutoff = 5;   %round(size(B,2)*0.05)
        soma_gene_cluster(find(sum(B')>min_cell_cutoff))=1;
        germ_gene_cluster = h_gene_cluster;
        germ_gene_cluster(find(h_gene_cluster>0))=1;
        input_gene_cluster = soma_gene_cluster;
        input_gene_cluster(find(h_gene_cluster==0 & soma_gene_cluster==0)) = 0; %germunexp-somaunexp
        input_gene_cluster(find(h_gene_cluster==0 & soma_gene_cluster >0)) = 1; %germunexp-somaexp
        input_gene_cluster(find(h_gene_cluster >0 & soma_gene_cluster==0)) = 2; %germexp-somaunexp
        input_gene_cluster(find(h_gene_cluster >0 & soma_gene_cluster >0)) = 3; %germexp-somaexp
        histcounts(input_gene_cluster)
        
        %Plot the germ-soma SNV plots
        fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
        fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
        figure;
        soma_together = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_cluster,...
                                     h_gname,...
                                     'F_compare_variants','only',...
                                     'FaceColor',soma_cm,...
                                     'Title_prefix','Genebody'); 
        ylim([0 270]);
    end
    

%For individual somatic cells.
figure;
only_germ_expressed = 0;
%binarize the gene clustering as germunexp/germexp:(somaexp, somaunexp)
for c=1:5
        C = B(:,find(soma_cluster==c));
        min_cell_cutoff = round(length(find(soma_cluster==c))*0.05)
        soma_gene_cluster(1:length(h_gene_cluster),1) = 0;%initialize the gene indeces.
        
        soma_gene_cluster(find(sum(C')>min_cell_cutoff))=1;
        if only_germ_expressed
        %h_gene_cluster should be >0
            input_gene_cluster = soma_gene_cluster(find(h_gene_cluster>0));
            input_genename = h_gname(find(h_gene_cluster>0));
        else
            input_gene_cluster = soma_gene_cluster;
            input_gene_cluster(find(h_gene_cluster==0 & soma_gene_cluster==0)) = 0; %germunexp-somaunexp
            input_gene_cluster(find(h_gene_cluster==0 & soma_gene_cluster >0)) = 1; %germunexp-somaexp
            input_gene_cluster(find(h_gene_cluster >0 & soma_gene_cluster==0)) = 2; %germexp-somaunexp
            input_gene_cluster(find(h_gene_cluster >0 & soma_gene_cluster >0)) = 3; %germexp-somaexp
            gene_count(c,1:4) = histcounts(input_gene_cluster)
            if gene_count(c,2)<2
                input_gene_cluster(find(input_gene_cluster>0))=input_gene_cluster(find(input_gene_cluster>0))-1;
                histcounts(input_gene_cluster)
            end
        input_genename = h_gname;
        end
    
    % Plot human genebody variants
   subplot(1,5,c)
   genebody=1;
if genebody
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
    %figure;
    soma(c).var = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_cluster,...
                                     input_genename,...
                                     'F_compare_variants','only',...
                                     'FaceColor',soma_cm,...
                                     'Title_prefix','Genebody'); 
    ylim([0 270]);
end
%     hold on;
%     text(1:2,repelem(20,2),num2cell(histcounts(input_gene_cluster)));hold off;
%     xticklabels({'somaunexp','somaexp'});
%     xtickangle(45)
%     suptitle(soma_name(c))
end



%###################
%Calculate the asymmetry scores of unexpressed genes across germ cells and somatic cells 
for c=1:5
        C = B(:,find(soma_cluster==c));
        min_cell_cutoff = round(length(find(soma_cluster==c))*0.05)
        soma_gene_cluster(1:length(h_gene_cluster),1) = 0;%initialize the gene indeces.
        soma_gene_cluster(find(sum(C')>min_cell_cutoff))=1;
   % Calculate the asymmetry scores using human intronic variants
   intron=1;
if intron
    fn1='YOUR_PATH\human_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.tsv'
     soma(c).var = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     soma_gene_cluster,...
                                     h_gname,...
                                      'F_compare_variants','no',...
                                    'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Human-intron-soma'); 
end
end

%Plot the asymmetry scores of unexp genes according to the germ cell and somatic cells
%This plot requires calculating the asymmetry scores of expressed and unexpressed genes
%in the testicular somatic cells and other somatic cells included for the study.
unexp_sample={'Germ','Sert','Ley','Myo','Mac','Endo','Guo','Hema','Mbr','Cort','kid'};
asym_type = soma(c).var.SNV_asym_type;
unexp_asyscore(1).stat = human_intron.AsymScore_expVSunexp_stat;
for c=1:5
    unexp_asyscore(c+1).stat = soma(c).var.AsymScore_expVSunexp_stat;
end

figure;
l = 6; %plot numbers
n = 6; %Number of samples with analysed unexpressed genes
for a=1:l
    subplot(1,l,a);    
    for i = 1:n
        bar(i,unexp_asyscore(i).stat(1,a,1),0.7,'FaceColor',[0.7 0.7 0.7]);hold on;
        errorbar(i,unexp_asyscore(i).stat(1,a,1),unexp_asyscore(i).stat(1,a,2),unexp_asyscore(i).stat(1,a,3),'.k');
    end
    hold off;
    xlim([0.5 n+0.5]); 
    ylim([-0.05 0.35]);
    title(asym_type(a));box off;
    xticks(1:n);   xticklabels(unexp_sample(1:n)); xtickangle(45);
    set(gca,'color','none');
end
suptitle(['Asymmetry scores of unexpressed genes across cell types'])





%% Add MSCI analysis
    % Load data: need 'h_gene_cluster','human_genename','m_gene_cluster','mouse_genename'
%    load('YOUR_PATH\h_m_stage_expression.mat');
    human_gene_info = struct2table(tdfread('YOUR_PATH\human_genes.tsv'));
    sexchr_genes = [strmatch('Y',human_gene_info.chromosome,'exact'); ...
                    strmatch('X',human_gene_info.chromosome,'exact')];
    
    [~,xi,xj] = intersect(human_genename,cellstr(human_gene_info.name(sexchr_genes,:)),'stable');
    h_X_geneid  = cellstr(human_gene_info.ensembl_id(sexchr_genes(xj),:));
    h_X_gname  = human_genename(xi,:);
    sexchr_gene_cluster = h_gene_cluster(xi); %gene cluster of X-chr genes
    histcounts(sexchr_gene_cluster)
    
    % plot X chr genes
        % heatmap
    exp = zscore(h_cluster_expression(xi,:)')';
    xticklength=length(exp(1,:));
    [i,j]=sort(sexchr_gene_cluster);
    figure;
    %draw heatmap as the genes sorted;
    subplot('position',[0.25 0.1 0.65 0.85]);
    imagesc(exp(j,:));colorbar;colormap(RdBu_cm)
    h_gca = gca;
    xlim([0.5, xticklength+0.5]);
    h_gca.XTick = 1:xticklength;
    xlabel('Cell clusters');
      h_gca.YTick = 1:length(h_X_gname);
      yticklabels(h_X_gname(j))
    title('Human sex-chr genes expression by clusters')
    %draw heatmap of gene clusters
    f=subplot('position',[0.1 0.1 0.05 0.85]);
    imagesc(i);xticks([]);yticks([]);colormap(f,genecluster_facecm);
    ylabel('Gene clusters');
    
    %by visualizing the expression, redefine the gene cluster as:
    %unexpressed, pre-MSCI-exp, post-MSCI-exp
    gene_cluster = sexchr_gene_cluster;
%     %sexchr_gene_cluster([find(gene_cluster==1);find(gene_cluster==2)])=1;
%     sexchr_gene_cluster([find(gene_cluster==3);find(gene_cluster==4);find(gene_cluster==5)])=3;
%     histcounts(sexchr_gene_cluster)
     sexchr_gene_cluster([find(gene_cluster==1);find(gene_cluster==2)])=1;
     sexchr_gene_cluster([find(gene_cluster==3);find(gene_cluster==4);find(gene_cluster==5)])=2;
     histcounts(sexchr_gene_cluster)
    
    
     
     
    % Plot germline SNV rates 
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
    %figure;
    MSCI_var = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     sexchr_gene_cluster,...
                                     h_X_geneid,...
                                     'F_compare_variants','only',...
                                     'FaceColor',genecluster_edgecm,...
                                     'Title_prefix','MSCI'); 
                                 
    ylim([0 270]);
     h_X_gname  = human_genename(xi,:);
    sexchr_gene_cluster = h_gene_cluster(xi); %gene cluster of X-chr genes
   





%% Add GTEx gene expression and mutation analaysis control. Fig.2 $SI Fig.1. %Done.

if GTEx_gene_number

% Plot GTEx gene number
    load('YOUR_PATH\GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.mat')
    %GTEx data is in TPM counts

    %imput gene signatures
    human_gene_info = struct2table(tdfread('YOUR_PATH\human_genes.tsv'));
    i=strmatch('protein',cellstr(human_gene_info.type));
    [gene,~,j] = intersect(cellstr(human_gene_info.name(i,:)),...
                        cellstr(GTExAnalysis20160115v7RNASeQCv1_gene(:,2)),'stable'); 
    GTEx_samplename = cellstr(GTExAnalysis20160115v7RNASeQCv1_Sample);
    GTEx_sampleexp  = GTExAnalysis20160115v7RNASeQCv1(j,:);
    
    %Check gene expression in GTEx.
    gix = strmatch('APOBEC',gene)
    [~,j]=setdiff(GTExAnalysis20160115v7RNASeQCv1_Sample,'Testis');
    [~,~,gix] = intersect(human_genename(find(h_gene_cluster ==0)),gene,'stable');
    B = GTEx_sampleexp(gix,j);
    length(find(max(B')>=1))
    figure;
    imagesc(log2(1+B),[0 5]);
    yticks(1:length(gix));yticklabels(gene(gix));  %yticklabels('')
    xticks(1:52);xticklabels(GTEx_samplename(j));xtickangle(45);
    colormap(flipud(Blues_cm(5:50,:)));
    xlabel('Human GTEx samples')
    title('GTEx gene expression(log2TPM) of spermatogenesis unexpressed genes')
    
    [~,~,k]=intersect(human_genename(find(h_gene_cluster ==0)),cellstr(human_gene_info.name),'stable');
    unexp_gene_id = cellstr(human_gene_info.ensembl_id(k,:));
    test=split(GTExAnalysis20160115v7RNASeQCv1_gene(:,1),'.');
    [~,~,k]= setdiff(unexp_gene_id,test(:,1),'stable');


minexp = [0.01 0.05 0.1 0.2 0.5 1 2 5 10];
clear gene_number
figure;
for i = 1:53
    for j=1:length(minexp)
        gene_number(j,i) = length(find(GTEx_sampleexp(:,i)>minexp(j)));
    end
    plot(1:length(minexp), flipud(gene_number(:,i))','LineWidth',1,'Color',[0.5 0.5 0.5] );hold on;
end
plot(1:length(minexp), flipud(gene_number(:,strmatch('Testi',GTEx_samplename)))',...
        '-o','LineWidth',3,'Color',[0.8 0 0]);hold on;
[x,xi]=sort(gene_number(1,:),'descend');
text(repmat(length(minexp),10,1),x(1:10),GTEx_samplename(xi(1:10)));
hold off;
xlim([0.5 length(minexp)+0.5]);xticks(1:length(minexp));
xticklabels(flip(minexp));
xlabel('Minimum expression level (mean TPM) cutoffs')
ylabel('Deteted genes expressed higher than cutoffs')
title('Human genes expressed across cell types in GTEx')
set(gca,'color','none')


    
    
end %Done


if GTEX_exp_grouping
%%  % use GTEx expression data only for gene clustering
    % load germ cells expression matrix and GTEx expression
    load('YOUR_PATH\h_m_stage_expression.mat');
    load('YOUR_PATH\GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.mat')
    %GTEx data is in TPM counts

    %imput gene signatures
    human_gene_info = struct2table(tdfread('YOUR_PATH\human_genes.tsv'));
    [input_GTEXgenes,~,i] = intersect(human_genename,...
                        cellstr(GTExAnalysis20160115v7RNASeQCv1_gene(:,2)),'stable'); 
    [~,j]=setdiff(GTExAnalysis20160115v7RNASeQCv1_Sample,'Testis');
    k=strmatch('Testis',cellstr(GTExAnalysis20160115v7RNASeQCv1_Sample));
    GTEx_samplename = cellstr(setdiff(GTExAnalysis20160115v7RNASeQCv1_Sample,'Testis'));
    GTEX_other_exp = GTExAnalysis20160115v7RNASeQCv1(i,j);
    GTEX_testis_exp = GTExAnalysis20160115v7RNASeQCv1(i,k);
    [~,~,xj] = intersect(input_GTEXgenes,cellstr(human_gene_info.name),'stable');
    input_genename  = cellstr(human_gene_info.ensembl_id(xj,:));

    
    length(find(GTEX_testis_exp<1))
    % Generate gene index
    max_gtex = max(GTEX_other_exp')';
    clear gene_cluster
    gene_cluster(1:length(input_genename)) = zeros;
    cutoff = 0.1;
    gene_cluster(find(GTEX_testis_exp<cutoff & max_gtex<cutoff)) = 0;
    gene_cluster(find(GTEX_testis_exp<cutoff & max_gtex>cutoff)) = 1;
    gene_cluster(find(GTEX_testis_exp>cutoff & max_gtex<cutoff)) = 2;
    gene_cluster(find(GTEX_testis_exp>cutoff & max_gtex>cutoff)) = 3;
    histcounts(gene_cluster)
    
    
    % plot heatmap of these genes expression pattern in GTEx
    exp = log10(GTEX_exp+1);
    xticklength=length(exp(1,:));
    [i,xi]=sort(gene_cluster);
    figure;
    %draw heatmap as the genes sorted;
    subplot('position',[0.25 0.1 0.65 0.85]);
    imagesc(exp(xi,:));colorbar;colormap(RdBu_cm)
    h_gca = gca;
    xlim([0.5, xticklength+0.5]);
    h_gca.XTick = 1:xticklength;
    xticklabels(GTEx_samplename);xtickangle(45);
    xlabel('Cell clusters');
    title('scRNA-seq unexp-gene expression from GTEx')
    %draw heatmap of gene clusters
    f=subplot('position',[0.1 0.1 0.05 0.85]);
    imagesc(i');xticks([]);yticks([]);colormap(f,genecluster_facecm);
    ylabel('Gene clusters');
    

%% test for each GTEx tissue. Fig.2. %Done
    load('YOUR_PATH\h_m_stage_expression.mat');
    load('YOUR_PATH\GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.mat')
    %GTEx data is in TPM counts

    %imput gene signatures
    human_gene_info = struct2table(tdfread('YOUR_PATH\human_genes.tsv'));
    [input_GTEXgenes,~,i] = intersect(human_genename,...
                        cellstr(GTExAnalysis20160115v7RNASeQCv1_gene(:,2)),'stable'); 
    GTEX_exp = GTExAnalysis20160115v7RNASeQCv1(i,:);
    [~,~,xj] = intersect(input_GTEXgenes,cellstr(human_gene_info.name),'stable');
    input_genename  = cellstr(human_gene_info.ensembl_id(xj,:));
    
    %setup the binary geneclusters as expressed or not expressed
%     GTEX_exp = GTEX_exp(find(sum(GTEX_exp')>0),:);
%     input_genename = input_genename(find(sum(GTEX_exp')>0),:);
    
        GTEX_exp_unexp_cluster = GTEX_exp;
        GTEX_exp_unexp_cluster(GTEX_exp_unexp_cluster>=0.1) = 1;
        GTEX_exp_unexp_cluster(GTEX_exp_unexp_cluster<0.1) = 0;
    
    %collect stats of GTEX exp versus unexp datasets
    
    figure;
    for i=1:53
    subplot(1,10,i-50);
    i
    gene_cluster = GTEX_exp_unexp_cluster(:,i);
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv';
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv';
    temp1 = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     gene_cluster,...
                                     input_genename,...
                                     'F_compare_variants','only',...
                                     'FaceColor',exunex_cm,...
                                     'Title_prefix',''); 
    ylim([0 290]);xlim([0.5 2.5]);yticklabels('')
    mean_GTEX_var(i,:) = temp1.variants_level_structure.mean;
%      fn1='YOUR_PATH\human_intron_12_result.12.stranded.tsv'
%      fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.tsv'
%     temp2 = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
%                                      gene_cluster,...
%                                      input_genename,... 
%                                      'F_compare_variants','only',...
%                                      'FaceColor',exunex_cm,...
%                                      'Title_prefix','GTEx-intron'); 
%     mean_GTEX_var_intron(i,:) = temp2.variants_level_structure.mean;
    end
    
    figure;
    for i=1:53
        scatter(mean_GTEX_var(i,2),mean_GTEX_var(i,1),50,i,'filled');hold on;
    end
    hold off; colormap(GTEX_cm)
    xlim([120 123]);    ylim([120 138]);  
    legend(GTExAnalysis20160115v7RNASeQCv1_Sample)
    title('GTEx-exp/unexp germline variants')
    ylabel('Germline mutations per kb in Unexp-genes')
    xlabel('Germline mutations per kb in Exp-genes')
    text(122,136,'<-Testis');
    set(gca,'Color','none')
    
    
    zscore_ratio=zscore(mean_GTEX_var(:,1)./mean_GTEX_var(:,2));
    figure; 
    violinplot(mean_GTEX_var(:,1)./mean_GTEX_var(:,2));
    %violinplot(zscore_ratio);
    %label Testis
    labels = sprintfc('%s', GTExAnalysis20160115v7RNASeQCv1_Sample(49));
    text(1.01,mean_GTEX_var(49,1)./mean_GTEX_var(49,2),labels);
    text(1.1,mean_GTEX_var(49,1)./mean_GTEX_var(49,2),num2str(zscore_ratio(49)));
    %label Ovary
    labels = sprintfc('%s', GTExAnalysis20160115v7RNASeQCv1_Sample(40));
    text(1.01,mean_GTEX_var(40,1)./mean_GTEX_var(40,2),labels);
    %label EBV-transformed lymphocytes
    labels = sprintfc('%s', GTExAnalysis20160115v7RNASeQCv1_Sample(22));
    text(1.01,mean_GTEX_var(22,1)./mean_GTEX_var(22,2),labels);
    %label whole blood
    labels = sprintfc('%s', GTExAnalysis20160115v7RNASeQCv1_Sample(53));
    text(1.01,mean_GTEX_var(53,1)./mean_GTEX_var(53,2),labels);
    ylabel('Germline mutations ratios(Unexp versus Exp-genes)')

    
    
    
end %Done






%% Multiple gene families. Fig.2 & SI Fig.3. %Done
    cd 'YOUR_PATH\cbrewer'
    genecluster_edgecm = [0.15 0.15 0.15;cbrewer('qual','Dark2',5)];
    genecluster_facecm = [0.5 0.5 0.5;cbrewer('qual','Set2',5)];
    RdBu_cm = flipud(cbrewer('div','RdBu',50));
    RdBu_cm = RdBu_cm(5:45,:);
    exunex_cm = [0.5 0.5 0.5;0.25 0.68 0.46];
    ten_type = cbrewer('qual','Dark2',10);
    cd 'YOUR_PATH'
    
%     %cluster genes by clusters of cells
%     load('YOUR_PATH\human_tese_replicates_germcells_ordered.mat');
%     A = h_germcell_sm_ordered;
%     %cluster genes with k=5
%     [h_cluster_expression,h_cluster_exp_meanUMI,h_gene_cluster]=F_cluster_gene(A,h_germ_IDX_ordered,5,...
%         'expression_cutoff_n',10,...
%         'umi_ratio_cutoff',0.1,...
%         'include_unexpressed_gene','yes',...
%         'ylimit',[-2 3]);
    
    %load gene_info file
    human_gene_info = struct2table(tdfread('YOUR_PATH\human_genes.tsv'));
    [~,~,xi] = intersect(human_genename,cellstr(human_gene_info.name),'stable');
    h_gname  = cellstr(human_gene_info.ensembl_id(xi,:));
    load('YOUR_PATH\h_m_stage_expression.mat');
    families = {'Ankyrin repeat domain containing (ANKRD)',...
                'Armadillo-like helical domain containing (ARMH)',...
                'Basic helix-loop-helix proteins (BHLH)',...
                'BTB domain containing (BTBD)',...
                'Cadherins',...
                'CD molecules (CD)',...
                'EF-hand domain containing',...
                'Fibronectin type III domain containing',...
                'GPCR, Class A rhodopsin-like(excluding ORs)',...
                'GPCR, Class A rhodopsin-like(Olfactory receptors)',...
                'Heat shock proteins',...
                'Helicases',...
                'Histones',...
                'Homeoboxes',...
                'Immunoglobulin superfamily domain containing',...
                'Ion channels',...
                'PDZ domain containing (PDZ)',... %
                'PHD finger proteins',... %
                'Pleckstrin homology domain containing (PLEKH)',... %
                'Ras small GTPases superfamily',...
                'Ring finger proteins',...
                'RNA binding motif containing (RBM)',...
                'Solute carriers (SLC)',...
                'WD repeat domain containing (WDR)',...
                'Zinc fingers C2H2-type',...
                'Zinc fingers - other',...
                'T cell receptor gene'};
    % Below gene families qualifies the criteria(>10 genes expressed/unexpressed)
        families = {'GPCR, Class A rhodopsin-like(Olfactory receptors)',...
                'Basic helix-loop-helix proteins (BHLH)',...
                'GPCR, Class A rhodopsin-like(excluding ORs)',...
                'CD molecules (CD)',...
                'EF-hand domain containing',...
                'Homeoboxes',...
                'Immunoglobulin superfamily domain containing',...
                'Ion channels',...
                'Solute carriers (SLC)',...
                'Zinc fingers - other'};
     % Reorder by alphabet of gene family names         
        families = {
                'Basic helix-loop-helix proteins (BHLH)',...
                'CD molecules (CD)',...
                'EF-hand domain containing',...
                'GPCR, Class A rhodopsin-like(excluding ORs)',...
                'Homeoboxes',...
                'Immunoglobulin superfamily domain containing',...
                'Ion channels',...
                'Solute carriers (SLC)',...
                'Zinc fingers - other'};
            
clear h_ranksumtest_p var_stats_mean var_stats_ci99
for  f = 1:27
    f
    %set up input gene family members
    filename = strcat('YOUR_PATH\Bioinform_Datasets\HUHO_Gene_N_C\',string(families(f)),'.txt');
    A = readtable(filename); 
    [~,~,idx]=intersect(cellstr(A.ApprovedSymbol),human_genename,'stable');
    input_genename = h_gname(idx);
    input_gene_cluster = h_gene_cluster(idx);
    input_gene_cluster(find(input_gene_cluster>0)) = 1;
    C_counts = histcounts(input_gene_cluster)
    
    if min(C_counts)>=10 & length(C_counts)==2
    % plot heatmap of OR genes expression pattern
%     exp = F_normalize(h_cluster_exp_meanUMI(idx,:),'tpm_log10');
%     xticklength=length(exp(1,:));
%     [i,xi]=sort(input_gene_cluster);
%     figure;
%     %draw heatmap as the genes sorted;
%     subplot('position',[0.25 0.1 0.65 0.85]);
%     imagesc(exp(xi,:));colorbar;colormap(RdBu_cm)
%     h_gca = gca;
%     xlim([0.5, xticklength+0.5]);
%     h_gca.XTick = 1:xticklength;
%      yticks(1:length(input_genename));
%      yticklabels(input_genename);
%     xlabel('Cell clusters');
%     title('Human gene expression by clusters')
%     %draw heatmap of gene clusters
%     f=subplot('position',[0.1 0.1 0.05 0.85]);
%     imagesc(i);xticks([]);yticks([]);colormap(f,genecluster_facecm);
%     ylabel('Gene clusters');
    
    
    genebody = 1;
    intron  = 1;
    genebody_asym = 0;
    if genebody
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
     F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_cluster,...
                                     input_genename,...
                                     'F_compare_variants','only',...
                                     'FaceColor',exunex_cm,...
                                     'Title_prefix','Genebody'); 
    end
    titlename= string(families(f))
    g_n = string(sprintf('unexp(C0,n=%d)-exp(C1,n=%d)',C_counts(1),C_counts(2)))
    suptitle({titlename;g_n})

    if intron
        % intron variation
     fn1='YOUR_PATH\human_intron_12_result.12.stranded.tsv'
     fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.tsv'
      F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_cluster,...
                                     input_genename,...
                                      'F_compare_variants','only',...
                                    'FaceColor',exunex_cm,...
                                     'Title_prefix','Intron'); 
        
        %input gene signatures
%         [~,i,j]=intersect(input_genename,human_gene_name_info(:,1),'stable');
%         gene_cluster = input_gene_cluster(i); %gene cluster
%         genename = human_gene_name_info(j,2);
%         exunex_cm; %gene cluster colormap
%         
%         % plot the overal variants level
%         [h_ranksumtest_p(:,2,f),var_stats_mean(2,1:n,f),var_stats_ci99(3:4,1:n,f)] = F_compare_paragene_variants(variants_number(a,1),...
%             variants_genelength(b,1),...
%             variants_genename(a),...
%             gene_cluster,genename,...
%             exunex_cm,'ylimit',[0 350]);
%         xticklabels({'C0','C1','C2','C3','C4','C5'});
    end 
    titlename= string(families(f))
    g_n = string(sprintf('unexp(C0,n=%d)-exp(C1,n=%d)',C_counts(1),C_counts(2)))
    suptitle({titlename;g_n})

    end

end

% plot the publishable results-genebody mutations
%Save the stats into two structures.
GF_gb = struct;
GF_intron = struct;

figure;
for  f =1:10
    %set up input gene family members
    filename = strcat('YOUR_PATH\HUHO_Gene_N_C\',string(families(f)),'.txt');
    A = readtable(filename); 
    [~,~,idx]=intersect(cellstr(A.ApprovedSymbol),human_genename,'stable');
    input_genename = h_gname(idx);
    input_gene_cluster = h_gene_cluster(idx);
    input_gene_cluster(find(input_gene_cluster>0)) = 1;
    C_counts = histcounts(input_gene_cluster)

    if min(C_counts)>=10 & length(C_counts)==2    
    
    subplot(1,10,f) 
        % plot the overal variants level
    genebody = 1;
    intron  = 0;

    if genebody
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
    GF_gb(f).var = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_cluster,...
                                     input_genename,...
                                     'F_compare_variants','only',...
                                     'FaceColor',exunex_cm,...
                                     'Title_prefix','Genebody'); 
    ylim([0 400]);
    end

    if intron
        % intron variation
     fn1='YOUR_PATH\human_intron_12_result.12.stranded.tsv'
     fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.tsv'
      GF_intron(f).var = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_cluster,...
                                     input_genename,...
                                      'F_compare_variants','only',...
                                    'FaceColor',exunex_cm,...
                                     'Title_prefix','Intron'); 
    ylim([0 400]);
    end 
    
    title(string(families(f)))
        
        if f>1
        ax1 = gca;  
        ax1.YAxis.Visible = 'off';
        elseif f==1
        ylabel('Germline mutations per kb')
        end
    
    end

end
suptitle('Genebody SNV level across gene families')
suptitle('Intronic SNV level across gene families')

input = GF_intron;
input = GF_gb;
figure;
for i=1:10
    var_stats_mean = input(i).var.variants_level_structure.mean;
    var_stats_ci99 = input(i).var.variants_level_structure.ci99;
    subplot(1,10,i);
        for j = 1:2
        bar(j,var_stats_mean(j),0.7,'FaceColor',exunex_cm(j,:));hold on;
        end
    errorbar(1:2,var_stats_mean(:),var_stats_ci99(1,:),var_stats_ci99(2,:),...
        '--','Color','k','Marker','.','MarkerSize',10,'LineWidth',1,'CapSize',4);hold off;
    ylim([0 300]);  xlim([0.3 2.7]);
    xticks(1:2); xticklabels({'C0','C1'});
    set(gca,'color','none');
    box off
        if i>1
        ax1 = gca;  
        ax1.YAxis.Visible = 'off';
        elseif i==1
        ylabel('Germline mutations per kb')
        end
end
%NOTE: p-value needs to be corrected for multiple test through Bonferroni methods.
suptitle('Genebody SNV level across gene families')
suptitle('Intronic SNV level across gene families')
     
     



    
%% Remove bidirectional gene pairs. SI Figure4. %DONE 
    %select bidirectional gene pairs
    human_gene_feature = struct2table(tdfread('YOUR_PATH\human_ensembl90_gene_features.txt','tab'));
    mouse_gene_feature = struct2table(tdfread('YOUR_PATH\mouse_ensembl90_gene_features.txt','tab'));   
    human_gene_feature.Gene_stable_ID=cellstr(human_gene_feature.Gene_stable_ID);
    human_gene_feature.Gene_name=cellstr(human_gene_feature.Gene_name);
    human_gene_feature.Chromosome_name=cellstr(human_gene_feature.Chromosome_name);
    mouse_gene_feature.Gene_stable_ID=cellstr(mouse_gene_feature.Gene_stable_ID);
    mouse_gene_feature.Gene_name=cellstr(mouse_gene_feature.Gene_name);
    mouse_gene_feature.Chromosome_name=cellstr(mouse_gene_feature.Chromosome_name);
    save('YOUR_PATH\gene_feature.mat',...
         'mouse_gene_feature', 'human_gene_feature');
    
    %load gene feature and setup initiates
    load('YOUR_PATH\gene_feature.mat');
    human = 1;
    mouse = 1;
    
if human
    gene_feature = human_gene_feature;
    gap_length_cutoff   = 5000; %starting point of gene pairs smaller than 5k
    chr = unique(gene_feature.Chromosome_name);
    %set up the empty output gene lists
    confounded_genes        = [];%including bidirecntional gene pairs
    bidirectional_genepairs = [];%
    for i=1:length(chr)
        chr(i)
        chr_gene_feature = gene_feature(strmatch(chr(i),gene_feature.Chromosome_name,'exact'),:);
        gene_start = chr_gene_feature.Gene_start;
        gene_start(find(chr_gene_feature.Strand==-1)) = chr_gene_feature.Gene_end(find(chr_gene_feature.Strand==-1));
        [gene_start_sorted,j]=sort(gene_start);
        chr_gene_feature = chr_gene_feature(j,:);
        if length(chr_gene_feature.Strand)>=2
        gap_length = gene_start_sorted(2:end) - gene_start_sorted(1:end-1);
        pair_direction = chr_gene_feature.Strand(2:end) - chr_gene_feature.Strand(1:end-1);
        idx1 = find(pair_direction==2 & gap_length<gap_length_cutoff);
        idx2 = idx1+1; %bona fide head-to-head bidirectional gene pairs
        idx  = find(pair_direction==-2);
        idx3 = [];
        for k=1:length(idx)
            idx_temp = [];
            %select confounded - strand genes   NOT RIGHT!
            if (chr_gene_feature.Gene_start(idx(k)) < chr_gene_feature.Gene_start(idx(k)+1) &...
                chr_gene_feature.Gene_end(idx(k)) > chr_gene_feature.Gene_end(idx(k)+1))
                idx_temp =  idx(k) + 1; 
            %select confounded + strand genes
            elseif (chr_gene_feature.Gene_start(idx(k)) > chr_gene_feature.Gene_start(idx(k)+1) &...
                    chr_gene_feature.Gene_end(idx(k)) <= chr_gene_feature.Gene_end(idx(k)+1))
                idx_temp =  idx(k); 
            %select confounded + strand genes
            elseif (chr_gene_feature.Gene_start(idx(k)) > chr_gene_feature.Gene_start(idx(k)+1) &...
                        chr_gene_feature.Gene_end(idx(k)) > chr_gene_feature.Gene_end(idx(k)+1))
                idx_temp =  [idx(k);idx(k)+1]; 
            end
            idx3 = [idx3;idx_temp];
        end
        
        gene_list1 = chr_gene_feature.Gene_name(idx1);
        gene_list2 = chr_gene_feature.Gene_name(idx2);
        gene_list3 = chr_gene_feature.Gene_name(idx3);
        bidirectional_genepairs = [bidirectional_genepairs;[gene_list1 gene_list2]];
        confounded_genes = [confounded_genes;gene_list1;gene_list2;gene_list3];
        end
    end
    clear chr chr_gene_feature gap_length gene_start gene_start_sorted pair_direction
    clear gene_list1 gene_list2 gene_list3 i idx idx1 idx2 idx3 idx_temp j k 
    
    human_confounded_genes = unique(confounded_genes);
    human_bidirectional_genes = bidirectional_genepairs;
    
    
    load('YOUR_PATH\gene_feature.mat')
    load('YOUR_PATH\h_m_stage_expression.mat');
    human_gene_info = struct2table(tdfread('YOUR_PATH\human_genes.tsv'));
    [~,~,xi] = intersect(human_genename,cellstr(human_gene_info.name),'stable');
    [~,~,xj] = intersect(human_confounded_genes,cellstr(human_gene_info.name),'stable');
    [h_gname_wo_confound_genes,xk]  = setdiff(cellstr(human_gene_info.ensembl_id(xi,:)),cellstr(human_gene_info.ensembl_id(xj,:)),'stable');
    h_g_cluster_wo_confound_genes = h_gene_cluster(xk);
    
    
%human intron
    fn1='YOUR_PATH\human_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.tsv'
     human_intron2 = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_g_cluster_wo_confound_genes,...
                                     h_gname_wo_confound_genes,...
                                    'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Human-intron w/o_confounded_genes'); 
 
%human upstream
    fn1='YOUR_PATH\human_upstream_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.upstream.singleFreq.tsv'
     human_upstream2 = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_g_cluster_wo_confound_genes,...
                                     h_gname_wo_confound_genes,...
                                     'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Human-upstream w/o_confounded_genes'); 
%human downstream
    fn1='YOUR_PATH\human_downstream_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.downstream.singleFreq.tsv'
     human_downstream2 = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_g_cluster_wo_confound_genes,...
                                     h_gname_wo_confound_genes,...
                                     'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Human-downstream w/o_confounded_genes'); 

                                 
% errorbar plot upstream-intron-downstream asymmetry scores together
ds = human_downstream2.AsymScore_expVSunexp_stat;
us = human_upstream2.AsymScore_expVSunexp_stat;
intron = human_intron2.AsymScore_expVSunexp_stat;
gc_name = human_intron2.gene_cluster_name;
asym_type = human_intron2.SNV_asym_type;

figure;
for a=1:6
    subplot(1,6,a);    
    line([0.5 6.5],[0 0],'Color',[0 0 0]); hold on;
    %downstream
    errorbar(1:6,ds(:,a,1),ds(:,a,2),ds(:,a,3),'-d',...
        'Color',[0.98 .57 .44],'MarkerSize',5,'LineWidth',1); hold on;
    %upstream
    errorbar(1:6,us(:,a,1),us(:,a,2),us(:,a,3),'-s',...
        'Color',[0.13 .44 .71],'MarkerSize',5,'LineWidth',1); hold on;
    %genebody
    errorbar(1:6,intron(:,a,1),intron(:,a,2),intron(:,a,3),'-o',...
        'Color',[0.79 .093 .11],'MarkerSize',5,'LineWidth',1); hold on;
    title(asym_type(a));box off;
    xlim([0.5 6.5]); 
    ylim([-0.25 0.6]);
    xticks(1:6);   xticklabels(gc_name);
    set(gca,'color','none');
end
end
  

if mouse
    gene_feature = mouse_gene_feature;
    gap_length_cutoff   = 5000; %starting point of gene pairs smaller than 5k
    chr = unique(gene_feature.Chromosome_name);
    %set up the empty output gene lists
    confounded_genes        = [];%including bidirecntional gene pairs
    bidirectional_genepairs = [];%
    for i=1:length(chr)
        chr(i)
        chr_gene_feature = gene_feature(strmatch(chr(i),gene_feature.Chromosome_name,'exact'),:);
        gene_start = chr_gene_feature.Gene_start;
        gene_start(find(chr_gene_feature.Strand==-1)) = chr_gene_feature.Gene_end(find(chr_gene_feature.Strand==-1));
        [gene_start_sorted,j]=sort(gene_start);
        chr_gene_feature = chr_gene_feature(j,:);
        if length(chr_gene_feature.Strand)>=2
        gap_length = gene_start_sorted(2:end) - gene_start_sorted(1:end-1);
        pair_direction = chr_gene_feature.Strand(2:end) - chr_gene_feature.Strand(1:end-1);
        idx1 = find(pair_direction==2 & gap_length<gap_length_cutoff);
        idx2 = idx1+1; %bona fide head-to-head bidirectional gene pairs
        idx  = find(pair_direction==-2);
        idx3 = [];
        for k=1:length(idx)
            idx_temp = [];
            %select confounded - strand genes   NOT RIGHT!
            if (chr_gene_feature.Gene_start(idx(k)) < chr_gene_feature.Gene_start(idx(k)+1) &...
                chr_gene_feature.Gene_end(idx(k)) > chr_gene_feature.Gene_end(idx(k)+1))
                idx_temp =  idx(k) + 1; 
            %select confounded + strand genes
            elseif (chr_gene_feature.Gene_start(idx(k)) > chr_gene_feature.Gene_start(idx(k)+1) &...
                    chr_gene_feature.Gene_end(idx(k)) <= chr_gene_feature.Gene_end(idx(k)+1))
                idx_temp =  idx(k); 
            %select confounded + strand genes
            elseif (chr_gene_feature.Gene_start(idx(k)) > chr_gene_feature.Gene_start(idx(k)+1) &...
                        chr_gene_feature.Gene_end(idx(k)) > chr_gene_feature.Gene_end(idx(k)+1))
                idx_temp =  [idx(k);idx(k)+1]; 
            end
            idx3 = [idx3;idx_temp];
        end
        
        gene_list1 = chr_gene_feature.Gene_name(idx1);
        gene_list2 = chr_gene_feature.Gene_name(idx2);
        gene_list3 = chr_gene_feature.Gene_name(idx3);
        bidirectional_genepairs = [bidirectional_genepairs;[gene_list1 gene_list2]];
        confounded_genes = [confounded_genes;gene_list1;gene_list2;gene_list3];
        end
    end
    clear chr chr_gene_feature gap_length gene_start gene_start_sorted pair_direction
    clear gene_list1 gene_list2 gene_list3 i idx idx1 idx2 idx3 idx_temp j k 
    
    mouse_confounded_genes = unique(confounded_genes);
    mouse_bidirectional_genes = bidirectional_genepairs;
    
    load('YOUR_PATH\h_m_stage_expression.mat');
    mouse_gene_info = struct2table(tdfread('YOUR_PATH\mouse_genes.tsv'));
    [~,~,xi] = intersect(mouse_genename,cellstr(mouse_gene_info.name),'stable');
    [~,~,xj] = intersect(mouse_confounded_genes,cellstr(mouse_gene_info.name),'stable');
    [m_gname2,xk]  = setdiff(cellstr(mouse_gene_info.ensembl_id(xi,:)),cellstr(mouse_gene_info.ensembl_id(xj,:)),'stable');
    m_g_cluster= m_gene_cluster(xk);
    
    cd 'YOUR_PATH\cbrewer'
    genecluster_facecm = [0.5 0.5 0.5;cbrewer('qual','Set2',5)];
    genecluster_facecm_ctl = [0.7 0.7 0.7;cbrewer('qual','Pastel2',5)];
    cd 'YOUR_PATH'
    
%mouse intron
    fn1='YOUR_PATH\mouse_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\mouse_ensembl90.intron.singleFreq.tsv'
     mouse_intron2 = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     m_g_cluster,...
                                     m_gname2,...
                                    'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Mouse-intron'); 
                                 
%mouse upstream region
    fn1='YOUR_PATH\mouse_upstream_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\mouse_ensembl90.upstream.singleFreq.tsv'
     mouse_upstream2 = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     m_g_cluster,...
                                     m_gname2,...
                                    'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Mouse-upstream'); 

%mouse downstream region
    fn1='YOUR_PATH\mouse_downstream_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\mouse_ensembl90.downstream.singleFreq.tsv'
     mouse_downstream2 = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     m_g_cluster,...
                                     m_gname2,...
                                    'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Mouse-downstream'); 

    end
  
    save('YOUR_PATH\gene_feature.mat',...
         'human_bidirectional_genes','human_confounded_genes',...
         'mouse_bidirectional_genes','mouse_confounded_genes','-append');  
    

if human_TCR_genes
    load('YOUR_PATH\human_tese1and2_replicates_germcells.mat');
    %human TCR related genes. (Nature Reviews Molecular Cell Biology volume 9, pages 958–970 (2008))
    marker={'EP300','ERCC1','ERCC2','ERCC3','ERCC4','ERCC5','ERCC6','ERCC8',...
        'HMGN1','LIG1','LIG3','PCNA','POLD1','POLD2','POLD3','POLD4','POLK',...
        'RFC1','RPA1','RPA2','RPA3','TCEA1','USP7','UVSSA','XAB2','XRCC1','XPA'};
    [~,~,goi]=intersect(marker,human_genename,'stable');
    figure;
    %TCR gene_by_gene
    left = 0.15; bottom = 0.45; width = 0.75; height = 0.5;
    s=subplot('position',[left bottom width height]); 
    imagesc(log2(1+h_germcell_sm(goi,h_germcell_pseudotimeorder)),[0 5]);colorbar;%colormap(flipud(corrcoef_map));
    set(gca, 'ytick', 1:length(goi));colormap(s,Blues_cm(5:50,:));
    set(gca, 'yticklabel', human_genename(goi));set(gca,'xtick','');
    ylabel('Human genes involved in transcription-coupled repair');
    %Germ cell IDX
    left = 0.15; bottom = 0.40; width = 0.75; height = 0.03;
    s=subplot('position',[left bottom width height]); 
    imagesc(h_germ_IDX(h_germcell_pseudotimeorder)');colormap(s,h_germtype_cm);
    set(gca,'xtick','');set(gca,'ytick','');
    %Monocle2 pseudotime
    left = 0.15; bottom = 0.36; width = 0.75; height = 0.03;
    s=subplot('position',[left bottom width height]); 
    imagesc(h_germcell_pseudotime');colormap(s,PuBuGn_cm(11:50,:));
    set(gca,'xtick','');set(gca,'ytick','');
    %violinplot of total TCR gene expression 
    sum_TCR=log2(1+sum(h_germcell_sm(goi,:)));
    left = 0.15; bottom = 0.04; width = 0.75; height = 0.3;
    s=subplot('position',[left bottom width height]); 
    violinplot(sum_TCR,h_germ_IDX)
%         for i=1:14
%         plot(find(h_germ_IDX==i),sum_TCR(find(h_germ_IDX==i)),...
%             'o','MarkerEdgeColor',h_germtype_cm(i,:));hold on;
%         end
%         xlim([0 length(sum_TCR)]);
        box off; hold off;set(gca,'color','none');ylim([0 8]);
        ylabel('log2(TCR genes)');set(gca,'xtick','');
        
    
    
end  %Done



     

%% Mouse mutation asymmetry analysis. SI Figure4. %DONE 
if mouse_variation_intotal

%mouse intron
    fn1='YOUR_PATH\mouse_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\mouse_ensembl90.intron.singleFreq.tsv'
     mouse_intron = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     m_gene_cluster,...
                                     m_gname,...
                                    'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Mouse-intron'); 
                                 
%mouse upstream region
    fn1='YOUR_PATH\mouse_upstream_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\mouse_ensembl90.upstream.singleFreq.tsv'
     mouse_upstream = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     m_gene_cluster,...
                                     m_gname,...
                                    'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Mouse-upstream'); 

%mouse downstream region
    fn1='YOUR_PATH\mouse_downstream_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\mouse_ensembl90.downstream.singleFreq.tsv'
     mouse_downstream = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     m_gene_cluster,...
                                     m_gname,...
                                    'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Mouse-downstream'); 
        
                                 
% errorbar plot upstream-intron-downstream asymmetry scores together
ds = mouse_downstream.AsymScore_expVSunexp_stat;
us = mouse_upstream.AsymScore_expVSunexp_stat;
intron = mouse_intron.AsymScore_expVSunexp_stat;
gc_name = mouse_intron.gene_cluster_name
asym_type = mouse_intron.SNV_asym_type  

figure;
for a=1:6
    subplot(1,6,a);    
    line([0.5 6.5],[0 0],'Color',[0 0 0]); hold on;
    %downstream
    errorbar(1:6,ds(:,a,1),ds(:,a,2),ds(:,a,3),'-d',...
        'Color',[0.98 .57 .44],'MarkerSize',5,'LineWidth',1); hold on;
    %upstream
    errorbar(1:6,us(:,a,1),us(:,a,2),us(:,a,3),'-s',...
        'Color',[0.13 .44 .71],'MarkerSize',5,'LineWidth',1); hold on;
    %genebody
    errorbar(1:6,intron(:,a,1),intron(:,a,2),intron(:,a,3),'-o',...
        'Color',[0.79 .093 .11],'MarkerSize',5,'LineWidth',1); hold on;
    title(asym_type(a));box off;
    xlim([0.5 6.5]); 
    ylim([-0.15 0.35]);
    xticks(1:6);   xticklabels(gc_name);
    set(gca,'color','none');
end


end  

     


     
     
%% Controls: Watson/Crick strand random gene clustering. SI Figure5. %DONE!
if asymmetry_score_controls
%do not specify strand of genes
    fn1='YOUR_PATH\human_intron_12_result.12.unstrand.tsv'
    fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.unstranded.tsv'
     human_intron_unstrand = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_gene_cluster,...
                                     h_gname,...
                                      'F_compare_variants','no',...
                                    'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Human-intron-unstrand '); 
                                 
%random assign gene cluster
    figure;
    s=subplot(1,2,1);
    explode = [1 0 0 0 0 0];
    pie(histcounts(h_gene_cluster),explode);colormap(s,genecluster_facecm)
        legend({'C0','C1','C2','C3','C4','C5'});
    s=subplot(1,2,2);
    h_gene_stage_shuffle=h_gene_cluster(randperm(length(h_gene_cluster)));
    pie(histcounts(h_gene_stage_shuffle));colormap(s,genecluster_facecm_ctl)
        legend({'C0','C1','C2','C3','C4','C5'});

    fn1='YOUR_PATH\human_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.tsv'
     human_intron_shuffle = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_gene_stage_shuffle,...
                                     h_gname,...
                                      'F_compare_variants','yes',...
                                    'FaceColor',genecluster_facecm_ctl,...
                                     'Title_prefix','Human-intron-shuffle-genecluster '); 
end                                 
    





%% human_variation_by_expressionlevel Done!

if human_variation_by_expressionlevel
    %% LOAD data
    %######NEED TO ADD CODES OF HOW TO GENERATE THE FILES BELOW
    % human variants
load('YOUR_PATH\h_m_stage_expression.mat');

    cd 'YOUR_PATH\cbrewer'
    %cbrewer()
    OrRd=cbrewer('seq','OrRd',12);
    explevel_edgecm = [0.15 0.15 0.15;OrRd(4:11,:)];
    explevel_facecm = [0.5 0.5 0.5;OrRd(4:11,:)];
    genecluster_facecm = [0.5 0.5 0.5;cbrewer('qual','Set2',5)];
    cd 'YOUR_PATH'
    
%% set up the expression level dependency
    
    h_mean_exp = log2(0.0001+mean(h_cluster_exp_meanUMI')');
    h_expressed_mean = h_mean_exp(find(h_gene_cluster>0));
        h_expressed_mean = h_mean_exp(find(h_gene_cluster>0));
        [f,xi] = ksdensity(h_expressed_mean);
        figure;plot(xi,f,'LineWidth',1.5);set(gca,'Color','none');box off;
        ylabel('Kernel gene density');xlabel('Expression level(log2)');
        xlim([-12 10]);xticks([-10 -8 -6 -4 -2 0 2 4 6 8]);
        
    h_expressed_idx  = find(h_gene_cluster>0);
    h_explevel_cluster(1:length(h_gene_cluster),1) = 0;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=-8) )) = 1;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=-6 & h_expressed_mean>-8) )) = 2;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=-4 & h_expressed_mean>-6) )) = 3;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=-2 & h_expressed_mean>-4) )) = 4;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=0 & h_expressed_mean>-2) )) = 5;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=2 & h_expressed_mean>0) )) = 6;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=4 & h_expressed_mean>2) )) = 7;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean>4) )) = 8;
    histcounts(h_explevel_cluster)
    
    %define gene cluster composation in each expression level cluster
    for i=1:8
        temp = h_gene_cluster(find(h_explevel_cluster == i));
        percentage(i,:) = histcounts(temp)./length(temp);
    end
    figure;bar(percentage,'stacked');colormap(genecluster_facecm(2:6,:));
    set(gca,'Color','none');xlim([0 9]);
    


if define_detailed_explevel_stage
    h_explevel_cluster2(1:length(h_gene_cluster),1) = 0;
    h_expressed_idx  = find(h_gene_cluster>0);
    h_explevel_cluster2( h_expressed_idx( find(h_expressed_mean<=-8) )) = 1;
    h_explevel_cluster2( h_expressed_idx( find(h_expressed_mean>4) )) = 102;

    mid_exp_range_idx = h_expressed_idx( find(h_expressed_mean<=-8 & h_expressed_mean>4 )); 
    mid_exp_range = prctile(h_expressed_mean,1):prctile(h_expressed_mean,99);
    for i = 1:100
        h_explevel_cluster2( h_expressed_idx( find( ...
            h_expressed_mean<=(-8+0.12*i) & ...
            h_expressed_mean>(-8+0.12*(i-1)) ))...
            ) = i+1;
    end    
end


%% mutation assymetry by expression level
    load('YOUR_PATH\h_m_stage_expression.mat');
    human_gene_info = struct2table(tdfread('YOUR_PATH\human_genes.tsv'));
    [~,~,xi] = intersect(human_genename,cellstr(human_gene_info.name),'stable');
    h_gname  = cellstr(human_gene_info.ensembl_id(xi,:));

%human intron: explevel cluster plot overall
    fn1='YOUR_PATH\human_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.tsv'
     human_intron_explevel = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_explevel_cluster,...
                                     h_gname,...
                                     'F_compare_variants','only',...
                                     'FaceColor',explevel_facecm,...
                                     'Title_prefix','Human-intron_explevel'); 


%human intron: explevel cluster asymmetry scores
    fn1='YOUR_PATH\human_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.tsv'
     human_intron_explevel = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_explevel_cluster,...
                                     h_gname,...
                                     'F_compare_variants','no',...
                                     'Ylim_fold_factor',1.55,...
                                     'FaceColor',explevel_facecm,...
                                     'Title_prefix','Human-intron'); 

%human intron: 100 explevel cluster for modeling
    fn1='YOUR_PATH\human_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.tsv'
     human_intron_explevel2 = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_explevel_cluster2,...
                                     h_gname,...
                                     'F_compare_variants','only',...
                                     'FaceColor',explevel_facecm,...
                                     'Title_prefix','Human-intron_explevel2');   

    mean_var = human_intron_explevel2.variants_level_structure.mean;
    ci_99    = human_intron_explevel2.variants_level_structure.ci99;
    figure;
    line([0 103],[mean_var(1) mean_var(1)],'Color','k','LineWidth',2);hold on;
    plot(smooth(mean_var(2:end),3),'k','LineWidth',2);hold on;
    plot(smooth((ci_99(1,2:end)+mean_var(2:end)'),3),'b');hold on;
    plot(smooth((ci_99(2,2:end)+mean_var(2:end)'),3),'b');hold off;
    ylim([50 150]);  xlim([-5 108]);
    set(gca,'Color','none');box off;
    xlabel('germline expression level')
    
    
    
end






%% Cross-species gene evolution analysis! %Done
try_Ape_to_human_Id_dNdS =1;
Rodents_to_mouse_Id_dNdS =1;

if try_Ape_to_human_Id_dNdS
    load('YOUR_PATH\Apes_to_human_Id_dNdS.mat');
    load('YOUR_PATH\h_m_stage_expression.mat');
        clear Variance dS dN dN_vs_dS
    %%Calculate input Variance,dS,dN_vs_dS;
    taxon_name={'Chimpazee','Gorilla','Bonobo','Orangutan','Gibbon'};
    Variance =100-Apes_to_human_Id_dNdS(:,[1 4 7 10 13]);
    dS       =Apes_to_human_Id_dNdS(:,[3 6 9 12 15]);
    dN       =Apes_to_human_Id_dNdS(:,[2 5 8 11 14]);
    dN_vs_dS =Apes_to_human_Id_dNdS(:,[2 5 8 11 14])./Apes_to_human_Id_dNdS(:,[3 6 9 12 15]);
    %input gene name and gene cluster
    human_genename;
    gene_cluster = h_explevel_cluster;   %h_gene_cluster; h_explevel_cluster

    
    figure;
    subplot(1,4,1)
    [ranksumtest_p] = F_compare_divergence(dS,Apes_to_human_Id_dNdS_genenames,...
                                                     gene_cluster,human_genename,...
                                                     'ylimit',[0 0.08]);
                                             ylabel('dS');ylim([0 0.08]);

    subplot(1,4,2)
    [ranksumtest_p] = F_compare_divergence(Variance,Apes_to_human_Id_dNdS_genenames,...
                                                     gene_cluster,human_genename,...
                                                     'ylimit',[0 8]);
                                             ylabel('Divergence');ylim([0 8]);
    subplot(1,4,3)
    [ranksumtest_p] = F_compare_divergence(dN,Apes_to_human_Id_dNdS_genenames,...
                                                     gene_cluster,human_genename,...
                                                     'ylimit',[0 0.03]);
                                             ylabel('dN');ylim([0 0.03]);
    subplot(1,4,4)
    [ranksumtest_p] = F_compare_divergence(dN_vs_dS,Apes_to_human_Id_dNdS_genenames,...
                                                     gene_cluster,human_genename,...
                                                     'ylimit',[0 0.6],...
                                                 'legend_name',taxon_name);
                                             ylabel('dN/dS');ylim([0 0.6]);
    suptitle('Apes-to-Human divergence data')
    
    
    %dN/dS ratio percentage
    gene_cluster = h_explevel_cluster;   %h_gene_cluster; h_explevel_cluster
    figure;
    for a=1:length(dN_vs_dS(1,:))
        for b=1:length(unique(gene_cluster))
            NS=dN_vs_dS(find(gene_cluster==(b-1)),a);
            NS = NS( ~any( isnan( NS ) | isinf( NS ), 2 ) ,: );%NS(isoutlier(NS))=[];
            NS_ratio(b,1)=length(find(NS<=0.5))/length(NS);  %strong purifying selection
            NS_ratio(b,2)=(length(find(NS>0.5))-length(find(NS>1)))/length(NS);  %strong purifying selection
            NS_ratio(b,3)=length(find(NS>1))/length(NS);  %positive selection 
        end
        subplot(1,length(dN_vs_dS(1,:)),a);
        bar(NS_ratio,'stacked');title(taxon_name(a));
        xlim([0.2 length(unique(gene_cluster))+0.8]);
    end
    

    
    %Fast evolving protein coding gene list
    [evo_gene,i,j]=intersect(Apes_to_human_Id_dNdS_genenames,human_genename,'stable');
    length(i)
    fast_evo_gene=evo_gene(find(nanmean(dN_vs_dS(i,:)')>1 & nanmean(dN_vs_dS(i,:)')<inf));

    for clu = 1:9
    cluster_number = clu-1;  %For unexpressed genes
    positive_hit=intersect(fast_evo_gene,human_genename(j(find(h_explevel_cluster(j)==cluster_number))),'stable');
    posi_in_sample_x = length(positive_hit)
    sample_N = length(human_genename(j(find(h_explevel_cluster(j)==cluster_number))))
    posi_in_population_K = length(fast_evo_gene)
    population_M = length(evo_gene)
    hygecdf_pvalue(1,clu) = hygecdf(posi_in_sample_x,population_M,posi_in_population_K,sample_N,'upper')
    %hygecdf(x,M,K,N)
    end
    
    
    
    %Large divergent protein coding genes
    [evo_gene,i,j] = intersect(Apes_to_human_Id_dNdS_genenames,human_genename,'stable');
    length(i)
    [~,x] = sort(nanmean(Variance(i,:)'),'descend','MissingPlacement','last');
    fast_evo_gene = evo_gene(x(1:round(length(evo_gene)*0.1)));

    for clu =1:9
    cluster_number = clu-1;  %For unexpressed genes
    positive_hit=intersect(fast_evo_gene,human_genename(j(find(h_explevel_cluster(j)==cluster_number))),'stable');
    posi_in_sample_x = length(positive_hit)
    sample_N = length(human_genename(j(find(h_explevel_cluster(j)==cluster_number))))
    posi_in_population_K = length(fast_evo_gene)
    population_M = length(evo_gene)
    hygecdf_pvalue(2,clu) = hygecdf(posi_in_sample_x,population_M,posi_in_population_K,sample_N,'upper')
    %hygecdf(x,M,K,N)
    end
    
    
end


if Rodents_to_mouse_Id_dNdS
    load('YOUR_PATH\Rodents_to_mouse_Id_dNdS.mat')
    load('YOUR_PATH\h_m_stage_expression.mat');
    clear Variance dS dN dN_vs_dS
     %%Calculate input Variance,dS,dN_vs_dS;
    taxon_name={'AlgerianMouse','RyukyuMouse','ShrewMouse','Rat','NA-deerMouse'};
    Variance=100-Rodents_to_mouse_IddNdS(:,[1 4 7 10 13]);
    dS      =Rodents_to_mouse_IddNdS(:,[3 6 9 12 15]);
    dN      =Rodents_to_mouse_IddNdS(:,[2 5 8 11 14]);
    dN_vs_dS=Rodents_to_mouse_IddNdS(:,[2 5 8 11 14])./Rodents_to_mouse_IddNdS(:,[3 6 9 12 15]);
    %input gene name and gene cluster
    mouse_genename;
    m_gene_cluster;   %m_gene_cluster; m_gene_exlevel_stage;
        
    
    
   figure;
    subplot(1,4,1)
    [ranksumtest_p] = F_compare_divergence(dS,Rodents_to_mouse_IddNdS_names,...
                                                     m_gene_cluster,mouse_genename,...
                                                    'legend_name',taxon_name);
                                             ylabel('dS');ylim([0 0.4]);
    subplot(1,4,2)
    [ranksumtest_p] = F_compare_divergence(Variance,Rodents_to_mouse_IddNdS_names,...
                                                     m_gene_cluster,mouse_genename,...
                                                    'legend_name',taxon_name);
                                             ylabel('Divergence');ylim([0 25]);
    subplot(1,4,3)
    [ranksumtest_p] = F_compare_divergence(dN,Rodents_to_mouse_IddNdS_names,...
                                                     m_gene_cluster,mouse_genename,...
                                                    'legend_name',taxon_name);
                                             ylabel('dN');ylim([0 0.1]);
    subplot(1,4,4)
    [ranksumtest_p] = F_compare_divergence(dN_vs_dS,Rodents_to_mouse_IddNdS_names,...
                                                     m_gene_cluster,mouse_genename,...
                                                 'legend_name',taxon_name);
                                             ylabel('dN/dS');ylim([0 0.4]);
    suptitle('Rodents-to-Mouse divergence data')

     
        
end

human_genename(find(h_explevel_cluster==8))

if human_mouse_evolution
    %define evolution score
    %zscore of divergence of human-chimpazee, mouse-rat pairs.
    load('YOUR_PATH\human_mouse_orthologue_pairs_20171018.mat');
    load('YOUR_PATH\h_m_stage_expression.mat');
    %setup the conserved gene list
    [~,i_H,~]=intersect(human_mouse_orthologue_pairs_20171018(:,1),human_genename,'stable');
    [~,i_M,~]=intersect(human_mouse_orthologue_pairs_20171018(:,2),mouse_genename,'stable');
    H_M_cons=intersect(i_H,i_M);
    [~,~,i_H]=intersect(human_mouse_orthologue_pairs_20171018(H_M_cons,1),human_genename,'stable');
    [~,~,i_M]=intersect(human_mouse_orthologue_pairs_20171018(H_M_cons,2),mouse_genename,'stable');
    gcl(:,1) = h_gene_cluster(i_H);
    gcl(:,2) = m_gene_cluster(i_M);
    
    load('YOUR_PATH\Apes_to_human_Id_dNdS.mat');
    load('YOUR_PATH\Rodents_to_mouse_Id_dNdS.mat')
    
    taxon_name={'AlgerianMouse','RyukyuMouse','ShrewMouse','Rat','NA-deerMouse'};
    Variance=100-Rodents_to_mouse_IddNdS(:,[1 4 7 10 13]);
    dS      =Rodents_to_mouse_IddNdS(:,[3 6 9 12 15]);
    idx = ~isnan(Variance(:,4));  %refer to rat
    div_m=Variance(idx,4)./mean(Variance(idx,4));
    div_m=zscore(Variance(idx,4));
    idx = ~isnan(dS(:,4));  %refer to rat
    dS_m=dS(idx,4)./mean(dS(idx,4));
    m_name = Rodents_to_mouse_IddNdS_names(idx);

    intersect(mouse_genename(i_M),Rodents_to_mouse_IddNdS_names,'stable');
    
    taxon_name={'Chimpazee','Gorilla','Bonobo','Orangutan','Gibbon'};
    Variance =100-Apes_to_human_Id_dNdS(:,[1 4 7 10 13]);
    dS       =Apes_to_human_Id_dNdS(:,[3 6 9 12 15]);
    idx  = ~isnan(Variance(:,1));  %refer to Chimpanzee
    div_h=Variance(idx,1)./mean(Variance(idx,1));
    div_h=zscore(Variance(idx,1));
    h_name = Apes_to_human_Id_dNdS_genenames(idx);
    idx  = ~isnan(dS(:,1));  %refer to Chimpanzee
    dS_h =dS(idx,1)./mean(dS(idx,1));
    dS_h =zscore(dS(idx,1));
    h_name = Apes_to_human_Id_dNdS_genenames(idx);
    
    
    gcl_idx{1} = intersect(find(gcl(:,1)==0),find(gcl(:,2)>0 ));%human_unexpressed, mouse_expressed
    gcl_idx{2} = intersect(find(gcl(:,1)>0),find(gcl(:,2)==0 ));%human_expressed, mouse_expressed
    gcl_idx{3} = intersect(find(gcl(:,1)==0),find(gcl(:,2)==0 ));%human_unexpressed, mouse_unexpressed
    gcl_idx{4} = intersect(find(gcl(:,1)==1),find(gcl(:,2)==1 ));%human_expressed_c1, mouse_expressed_c1
    gcl_idx{5} = intersect(find(gcl(:,1)==2),find(gcl(:,2)==2 ));%human_expressed_c2, mouse_expressed_c2
    gcl_idx{6} = intersect(find(gcl(:,1)==3),find(gcl(:,2)==3 ));%human_expressed_c3, mouse_expressed_c3
    gcl_idx{7} = intersect(find(gcl(:,1)==4),find(gcl(:,2)==4 ));%human_expressed_c4, mouse_expressed_c4
    gcl_idx{8} = intersect(find(gcl(:,1)==5),find(gcl(:,2)==5 ));%human_expressed_c5, mouse_expressed_c5

    
    
    
    figure;
    for i=1:8
    subplot(1,8,i);   
    idx = gcl_idx{i};
    [~,~,xi]=intersect(human_genename(i_H(idx)),h_name,'stable');
            asym_score_stat(1,1,i) = mean(div_h(xi));
            asym_score_stat(1,2:3,i) = bootci(10000,{@mean,div_h(xi)},'alpha',0.01)-mean(div_h(xi));;
    [~,~,xj]=intersect(mouse_genename(i_M(idx)),m_name,'stable');
            asym_score_stat(2,1,i) = mean(div_m(xi));
            asym_score_stat(2,2:3,i) = bootci(10000,{@mean,div_m(xi)},'alpha',0.01)-mean(div_m(xi));;

            asym_score_ranksum_p(i) = ranksum(div_h(xi),div_m(xi));
            bar(1,asym_score_stat(1,1,i),0.8,'FaceColor',[252,141,89]/255);hold on;
            bar(2,asym_score_stat(2,1,i),0.8,'FaceColor',[145,207,96]/255);hold on;
            errorbar(1:2,asym_score_stat(1:2,1,i),asym_score_stat(1:2,2,i),asym_score_stat(1:2,3,i),'.k');hold on;
        xlim([0.5 2.5]); 
            text(1.5,repelem(max(asym_score_stat(1:2,1,1))*1.1,1),strcat('p=',num2str(asym_score_ranksum_p(i))));hold off;

    end
    
    
    
    
     taxon_name={'Chimpazee','Gorilla','Bonobo','Orangutan','Gibbon'};
    Variance =100-Apes_to_human_Id_dNdS(:,[1 4 7 10 13]);
    dS       =Apes_to_human_Id_dNdS(:,[3 6 9 12 15]);
   figure;
    for i=1:5
        subplot(1,5,i)
    idx  = ~isnan(Variance(:,i));  %refer to Chimpanzee
    div_h=Variance(idx,i)./mean(Variance(idx,i));
    h_name = Apes_to_human_Id_dNdS_genenames(idx);

    [ranksumtest_p] = F_compare_divergence(div_h,h_name,...
                                           gcl(:,1),human_genename(i_H),...
                                                 'legend_name','chimp');
                                             ylabel('dS-normalized');
    end
    
end





%% highly expressed genes GO term. 
load('YOUR_PATH\h_m_stage_expression.mat');
    h_mean_exp = log2(0.0001+mean(h_cluster_exp_meanUMI')');
    h_expressed_mean = h_mean_exp(find(h_gene_cluster>0));
        h_expressed_mean = h_mean_exp(find(h_gene_cluster>0));
        
    h_expressed_idx  = find(h_gene_cluster>0);
    h_explevel_cluster(1:length(h_gene_cluster),1) = 0;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=-8) )) = 1;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=-6 & h_expressed_mean>-8) )) = 2;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=-4 & h_expressed_mean>-6) )) = 3;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=-2 & h_expressed_mean>-4) )) = 4;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=0 & h_expressed_mean>-2) )) = 5;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=2 & h_expressed_mean>0) )) = 6;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=4 & h_expressed_mean>2) )) = 7;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean>4) )) = 8;
    histcounts(h_explevel_cluster)
    
    %save the h_explevel_cluster results for further usage
    save('YOUR_PATH\h_m_stage_expression.mat',...
        'h_explevel_cluster', '-append');

    
    unexpressed_genes = human_genename(find(h_explevel_cluster==0));
    high_expressed_genes = human_genename(find(h_explevel_cluster==8));

    unexp_gene = readtable('\GOPROCESS_unexp.xlsx');
    Go_unexp = {'G protein-coupled receptor signaling pathway',...
                'detection of stimulus',...
                'keratinization',...
                'immune response',...
                'sensory perception',...
                'defense response',...
                'humoral immune response',...
                'innate immune response',...
                'adaptive immune response',...
                'response to biotic stimulus'};
    [~,~,xi] = intersect(Go_unexp,cellstr(unexp_gene.Description),'stable');
    figure;barh(flipud(-log10(unexp_gene.FDRQ_value(xi))));
    yticks(1:10);yticklabels(flip(Go_unexp));
    ylim([0 11]); 
    set(gca,'Color','none');box off;
    xlabel('-log10(FDRQ_value)')


    highexp_gene = readtable('\GOPROCESS_highexp.xlsx');
    Go_highexp = {'reproductive process',...
                'gamete generation',...
                'protein localization to organelle',...
                'RNA catabolic process',...
                'spermatogenesis',...
                'peptide metabolic process',...
                'developmental process involved in reproduction',...
                'nucleic acid metabolic process'};
    [~,~,xi] = intersect(Go_highexp,cellstr(highexp_gene.Description),'stable');
    figure;barh(flipud(-log10(highexp_gene.FDRQ_value(xi))));
    ylim([0 9]); yticks(1:8);yticklabels(flip(Go_highexp));
    set(gca,'Color','none');box off;
    xlabel('-log10(FDRQ_value)')
    

    
    

%% 2019-05 192-type of mutations, considering sequence context (Human & Mouse)
%## LOAD data
    load('YOUR_PATH\h_m_stage_expression.mat');
    human_gene_info = struct2table(tdfread('YOUR_PATH\human_genes.tsv'));
    [~,~,xi] = intersect(human_genename,cellstr(human_gene_info.name),'stable');
    h_gname  = cellstr(human_gene_info.ensembl_id(xi,:));

    mouse_gene_info = struct2table(tdfread('YOUR_PATH\mouse_genes.tsv'));
    [~,~,xi] = intersect(mouse_genename,cellstr(mouse_gene_info.name),'stable');
    m_gname  = cellstr(mouse_gene_info.ensembl_id(xi,:));

    cd 'YOUR_PATH\cbrewer'
    %cbrewer()
    expgene_cm = [29,180,170]/256;
    genecluster_facecm = [0.5 0.5 0.5;cbrewer('qual','Set2',5)];
    genecluster_facecm_ctl = [0.7 0.7 0.7;cbrewer('qual','Pastel2',5)];
    cd 'YOUR_PATH'

if human_genebody
    mutnumb='YOUR_PATH\human_gene_body_result.192.stranded.tsv'
    trifreq='YOUR_PATH\human_ensembl90.genebody.triFreq.tsv'
    %plot by gene cluster
    human_gene_body = F_compare_snv_192(  mutnumb,   trifreq,...
                               h_gene_cluster,...
                               h_gname,...
                               'FaceColor',genecluster_facecm,... 
                               'MutRate_ylimit',[0 151],... 
                               'AsySco_ylimit',[-0.35 1]); 

    %plot pyrimidine rules
    F_triref_asymmetry_score_scatter(human_gene_body,'scatter','yes',...
            'paired_comparison','yes','ylimit',[-0.35 1],'sample_name','Human genebody')
    
    %plot by expression level
    human_gene_body_explevel = F_compare_snv_192(  mutnumb,   trifreq,...
                               h_explevel_cluster,...
                               h_gname,...
                               'FaceColor',explevel_facecm,... 
                               'MutRate_ylimit',[0 151],... 
                               'AsySco_ylimit',[-0.35 1]); 
    %Plot the AS_by_explevel of all types.
    figure;
    for j=1:6
        subplot(1,6,j);
    for i = 1:16
        plot(1:9,test_explevel.asymmetry_score(:,j*16-16+i));hold on;
    end
    xlim([0.5 9.5]);ylim([-.35 1]);
    end    
end

if human_intron
    mutnumb='YOUR_PATH\human_intron_192_result.192.stranded.tsv'
    trifreq='YOUR_PATH\human_ensembl90.intron.triFreq.tsv'
    %plot by gene cluster
    human_intron = F_compare_snv_192(  mutnumb,   trifreq,...
                               h_gene_cluster,...
                               h_gname,...
                               'FaceColor',genecluster_facecm,... 
                               'MutRate_ylimit',[0 151],... 
                               'AsySco_ylimit',[-0.35 1]); 
    %plot pyrimidine rules
    F_triref_asymmetry_score_scatter(human_intron,'scatter','no',...
            'paired_comparison','yes','ylimit',[-0.35 1],...
            'plot_gene_cluster',[2:6],'sample_name','Human intron')

    %Human intron plot on all human expressed genes.
    gc=h_gene_cluster;
    gn = h_gname;    
    gc(gc>0)=1;
    expgene_cm = [29,180,170]/256;
    human_intron_allexpgenes = F_compare_snv_192(  mutnumb,   trifreq,...
                               gc(find(gc==1)),...
                               gn(find(gc==1)),...
                               'FaceColor', expgene_cm,... 
                               'MutRate_ylimit',[0 151],... 
                               'AsySco_ylimit',[-0.35 1],... 
                               'Title_prefix','Human-intron-allexpgenes'); 
                           
    %plot by expression level
    human_intron_explevel = F_compare_snv_192(  mutnumb,   trifreq,...
                               h_explevel_cluster,...
                               h_gname,...
                               'FaceColor',explevel_facecm,... 
                               'MutRate_ylimit',[0 151],... 
                               'AsySco_ylimit',[-0.35 1]); 
    %Plot the AS_by_explevel of all types.
    figure;
    for j=1:6
        subplot(1,6,j);
    for i = 1:16
        plot(1:9,human_intron_explevel.asymmetry_score(:,j*16-16+i));hold on;
    end
    xlim([0.5 9.5]);ylim([-.35 1]);
    end    
end



%Mouse genebody
if mouse_genebody
    mutnumb='YOUR_PATH\mouse_gene_body_result.192.stranded.tsv'
    trifreq='YOUR_PATH\mouse_ensembl90.genebody.triFreq.tsv'
    mouse_192 = F_compare_snv_192(  mutnumb,   trifreq,...
                               m_gene_cluster,...
                               m_gname,...
                               'FaceColor',genecluster_facecm,... 
                               'MutRate_ylimit',[0 30],... 
                               'AsySco_ylimit',[-0.25 0.5]); 
    %plot pyrimidine rules
    F_triref_asymmetry_score_scatter(mouse_192,'scatter','yes',...
            'paired_comparison','yes','ylimit',[-0.25 0.5],'sample_name','Mouse genebody')
end

if mouse_intron
    mutnumb='YOUR_PATH\mouse_intron_192_result.192.stranded.tsv'
    trifreq='YOUR_PATH\mouse_ensembl90.intron.triFreq.tsv'
    mouse_intron = F_compare_snv_192(  mutnumb,   trifreq,...
                               m_gene_cluster,...
                               m_gname,...
                               'FaceColor',genecluster_facecm,... 
                               'MutRate_ylimit',[0 30],... 
                               'AsySco_ylimit',[-0.25 0.5]); 
    %plot pyrimidine rules
    F_triref_asymmetry_score_scatter(mouse_192,'scatter','yes',...
            'paired_comparison','yes','ylimit',[-0.25 0.5],'sample_name','Mouse intron')

    %plot on all mouse expressed genes.
    gc=m_gene_cluster;
    gn = m_gname;    
    gc(gc>0)=1;
    mouse_intron_allexpgenes = F_compare_snv_192(  mutnumb,   trifreq,...
                               gc(find(gc==1)),...
                               gn(find(gc==1)),...
                               'FaceColor', expgene_cm,... 
                               'MutRate_ylimit',[0 30],... 
                               'AsySco_ylimit',[-0.2 0.35],... 
                               'Title_prefix','Mouse-intron-allexpgenes'); 

end        
    



    
%% Human denovo mutations analysis (All)
% These analysis based on the the control samples of An et al, Science, 2018
% and the Iceland trio WGS dataset (Jonsson et al, Nature, 2017)

%## LOAD gene expression data
    load('YOUR_PATH\h_m_stage_expression.mat');
    human_gene_info = struct2table(tdfread('YOUR_PATH\human_genes.tsv'));
    [~,~,xi] = intersect(human_genename,cellstr(human_gene_info.name),'stable');
    h_gname  = cellstr(human_gene_info.ensembl_id(xi,:));

    cd 'YOUR_PATH\cbrewer'
    %cbrewer()
    genecluster_edgecm = [0.15 0.15 0.15;cbrewer('qual','Dark2',5)];
    OrRd=cbrewer('seq','OrRd',12);
    explevel_edgecm = [0.15 0.15 0.15;OrRd(4:11,:)];
    expgene_cm = [29,180,170]/256;
    genecluster_facecm = [0.5 0.5 0.5;cbrewer('qual','Set2',5)];
    genecluster_facecm_ctl = [0.7 0.7 0.7;cbrewer('qual','Pastel2',5)];
    cd 'YOUR_PATH'

%NOTE: the denovo mutations only include autosomes. 
prot_idx = strmatch('protein_coding',cellstr(human_gene_info.type));
autosome={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'};
auto = ismember(cellstr(human_gene_info.chromosome(prot_idx,:)),autosome);
auto_prot_genes = cellstr(human_gene_info.ensembl_id(prot_idx(auto),:));
[~,~,j]=intersect(auto_prot_genes,h_gname,'stable');

if DNM_in_total
% An_ctl_&_Jonsson_all DNM rates across gene clusters
fn1="YOUR_PATH\Jonsson_2017_Nature_SNV_all_result.12.stranded.tsv"
fn2="YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv"
fn3="YOUR_PATH\An_2018_SNV_control_result.12.stranded.tsv"
    %by gene clusters
    An_Jonsson_all = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_gene_cluster(j),...
                                     h_gname(j),...
                                     'Addup_SNV_file',fn3,...
                                     'Mode','single_gene',...
                                     'F_compare_variants','only',...
                                     'Include_meta_all','no',...
                                     'Ylim_fold_factor',1.0,...
                                     'FaceColor',genecluster_edgecm,... 
                                     'Title_prefix','An-Jonsson-all'); 
    %by expression level cluster
    An_Jonsson_all = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     h_explevel_cluster(j),...
                                     h_gname(j),...
                                     'Addup_SNV_file',fn3,...
                                     'Mode','single_gene',...
                                     'F_compare_variants','yes',...
                                     'Include_meta_all','no',...
                                     'Ylim_fold_factor',1.0,...
                                     'FaceColor',explevel_edgecm,... 
                                     'Title_prefix','An-Jonsson-all'); 
end
       

if DNM_192
% An_ctl_&_Jonsson_all DNM rates asymmetries referring to adjacent bases in
%   the EXPRESSED AUTOSOME genes.
    gc = h_gene_cluster(j); %select only the autosome genes.
    gn = h_gname(j);  
    gc(gc>0)=1; %Set the autosome expressed genes as cluster 1.
    An_ctl ="YOUR_PATH\An_2018_SNV_control_result.192.stranded.tsv"
    Jonsson="YOUR_PATH\Jonsson_2017_Nature_SNV_all_result.192.stranded.tsv"
    trifreq="YOUR_PATH\human_ensembl90.genebody.triFreq.tsv"
    %plot by gene cluster
    An_Jonsson_all_192 = F_compare_snv_192(  An_ctl,   trifreq,...
                               gc(find(gc==1)),...
                               gn(find(gc==1)),...
                               'Addup_SNV_file',Jonsson,...
                               'FaceColor',expgene_cm,... 
                               'MutRate_ylimit',[0 0.11],... 
                               'MutRate_ylimit2',[0.3 1],... 
                               'AsySco_ylimit',[-0.5 1],... 
                                'Title_prefix','An Ctl & Jonsson'); 
%plot pyrimidine rules
%     F_triref_asymmetry_score_scatter(An_Jonsson_all_192,'scatter','no',...
%             'paired_comparison','yes','ylimit',[-0.5 1],...
%             'plot_gene_cluster',[1],'sample_name','An Ctl & Jonsson DNM')
        
%For comparison, calculate the dbSNP mutation rates asymmetries referring  
%   to adjacent bases in the expressed AUTOSOME genes.
    mutnumb="YOUR_PATH\human_intron_192_result.192.stranded.tsv"
    trifreq="YOUR_PATH\human_ensembl90.intron.triFreq.tsv"
    human_intron_autosome_192 = F_compare_snv_192(  mutnumb,   trifreq,...
                               gc(find(gc==1)),...
                               gn(find(gc==1)),...
                               'FaceColor', expgene_cm,... 
                               'MutRate_ylimit',[0 151],... 
                               'AsySco_ylimit',[-0.35 1],... 
                               'Title_prefix','Human-intron-allexpAUTOSOMEgenes'); 
end


if DNM_and_SNP_regression
    %The analysis depends on the results of in DNM_192 section
    %Excluding the C>T mutation rates at the CpG context, which is too high.
    %DNM = An_Jonsson_all_DNM; %AS: 0.794 %CS: 0.982 %TS: 0.987
    %DNM = An_ctl_DNM; %AS: 0.647 %CS: 0.97 %TS: 0.977
    %DNM = Jonsson_DNM; %AS: 0.548 %CS: 0.954 %TS: 0.944
    
    DNM = An_Jonsson_all_192; 
    %R2: AS: 0.749; CS: 0.975; TS: 0.982
    %R2: CS_w/t_C>T@CpG: 0.986; TS_w/t_C>T@CpG: 0.983
    %F-test p-value: AS: 6.2e-30; CS: 4.24e-74; TS: 1.05e-80
    %F-test p-value: CS_w/t_C>T@CpG: 1.77e-88; TS_w/t_C>T@CpG: 1.5e-85
                              

muttype_face = [189,215,231;252,174,145;204,204,204;203,201,226;186,228,179;254,217,142]/256;
muttype_edge = [33,113,181;203,24,29;82,82,82;106,81,163;35,139,69;204,76,2]/256;

    % Asymmetry score regression
    x = DNM.asymmetry_score;
    y = human_intron_autosome_192.asymmetry_score;
    corrcoef(x,y)
    dbSNP_to_DNM=fitlm(x,y)
    figure; plot(dbSNP_to_DNM); hold on;
    for i=1:6
        scatter(x(i*16-15:i*16),y(i*16-15:i*16),'MarkerEdgeColor',muttype_edge(i,:),...
                                                'MarkerFaceColor',muttype_face(i,:));hold on
    end
    hold off; box off;
    xlim([-0.5 1]);ylim([-0.5 1]);
    axis equal; 
    ax = gca;
    set(ax,'Color','none')
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    xlabel('Asymmetry scores (DNM)');
    ylabel('Asymmetry scores (dbSNP SNV)');
    title('DNM-to-dbSNP SNV Asymmetry scores regression');

    % Coding strand mutation rates regression
    mr_range = [1:96];
    mr_range = [1:52,57:96];
    
    x = DNM.coding_mut_rate(mr_range)* ...
        sum(human_intron_autosome_192.coding_mut_level(mr_range)) ...
        / sum(DNM.coding_mut_level(mr_range));
    y = human_intron_autosome_192.coding_mut_rate(mr_range);
    corrcoef(x,y)
    dbSNP_to_DNM=fitlm(x,y)
    figure; plot(dbSNP_to_DNM); hold on;
    for i=1:6
        if i<4
        scatter(x(i*16-15:i*16),y(i*16-15:i*16),'MarkerEdgeColor',muttype_edge(i,:),...
                                                'MarkerFaceColor',muttype_face(i,:));hold on;
        elseif i==4
        scatter(x(49:60),y(49:60),'MarkerEdgeColor',muttype_edge(i,:),...
                                  'MarkerFaceColor',muttype_face(i,:));hold on;
        elseif i>4
        scatter(x(i*16-15-4:i*16-4),y(i*16-15-4:i*16-4),'MarkerEdgeColor',muttype_edge(i,:),...
                                                'MarkerFaceColor',muttype_face(i,:));hold on;
        end
    end
    hold off; box off;
    xlim([0 150]);ylim([0 150]);
    ax = gca;
    set(ax,'Color','none')
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    xlabel('CS mutation rates (DNM, normalized by total #SNV)');
    ylabel('CS mutation rates (dbSNP SNV)');
    title('DNM-to-dbSNP SNV CS mutation rates regression');
    
    % Template strand mutation rates regression
    mr_range = [1:96];
    mr_range = [1:52,57:96];
    x = DNM.template_mut_rate(mr_range)* ...
        sum(human_intron_autosome_192.template_mut_level(mr_range)) ...
        / sum(DNM.template_mut_level(mr_range));
    y = human_intron_autosome_192.template_mut_rate(mr_range);
    corrcoef(x,y)
    dbSNP_to_DNM=fitlm(x,y)
    figure; plot(dbSNP_to_DNM); hold on;
    for i=1:6
        if i<4
        scatter(x(i*16-15:i*16),y(i*16-15:i*16),'MarkerEdgeColor',muttype_edge(i,:),...
                                                'MarkerFaceColor',muttype_face(i,:));hold on;
        elseif i==4
        scatter(x(49:60),y(49:60),'MarkerEdgeColor',muttype_edge(i,:),...
                                  'MarkerFaceColor',muttype_face(i,:));hold on;
        elseif i>4
        scatter(x(i*16-15-4:i*16-4),y(i*16-15-4:i*16-4),'MarkerEdgeColor',muttype_edge(i,:),...
                                                'MarkerFaceColor',muttype_face(i,:));hold on;
        end
    end
    hold off; box off;
    xlim([0 105]);ylim([0 105]);
    ax = gca;
    set(ax,'Color','none')
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    xlabel('TS mutation rates (DNM, normalized by total #SNV)');
    ylabel('TS mutation rates (dbSNP SNV)');
    title('DNM-to-dbSNP SNV TS mutation rates regression');
    
    
end







