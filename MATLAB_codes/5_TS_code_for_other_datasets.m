    load('YOUR_PATH\h_m_stage_expression.mat');
    human_gene_info = struct2table(tdfread('YOUR_PATH\human_genes.tsv'));
    [~,~,xj] = intersect(human_genename,cellstr(human_gene_info.name),'stable');
    h_gname  = cellstr(human_gene_info.ensembl_id(xj,:));

%% 1. Human spermatogenesis from Guo et al, CellRes, 2018

 fn='YOUR_PATH\Other scRNA-seq data\Guo et al, CellRes, 2018, human spermatogenesis\GSE112013_Combined_UMI_table.txt'
    [Guo_sperm,Guo_cn,Guo_gn] = F_tsv_import(fn);   

if try_Guo
    rng(1) % for reproducibility
    [~,Guo_PCA_score,~,~,Guo_tSNE_score] = F_kmeans(Guo_sperm,...
                                       'k_for_kmeans',13,...
                                       'SelectVariableGene','yes',...
                                       'number_of_pc',50,...
                                       'plot_PCA','yes',...
                                       'plot_tSNE','yes',...
                                       'visual_mode','off'...
                                       );
    A = log10(1+F_normalize(Guo_sperm,'median'));
    A_zscore = zscore(A,0,2);
    marker=cellstr({'WT1','SOX9','DDX4','CSF1','IGF1','FGFR3','MYH11','ACTA2','SYCP3','CD68','CD163','ACRV1','PECAM1','VWF','PRM1'});
    figure;
    for i=1:15
        goi = strmatch(marker(i),Guo_gn,'exact');
        if (length(goi)~=1)
            continue
        else
        subplot(5,3,i);
         scatter(Guo_tSNE_score(:,1),Guo_tSNE_score(:,2),1,A_zscore(goi,:),'filled');
        %xlabel('tSNE1'); ylabel('tSNE2');
        colormap(PuBuGn_cm(10:end,:));title(Guo_gn(goi));
        set(gca,'color','none');xticklabels({});yticklabels({})
        end
    end
end

    %Use the cell cluster from original paper.
    fn = "YOUR_PATH\Other scRNA-seq data\Guo et al, CellRes, 2018, human spermatogenesis\Guo_adult_cells_cluster.tsv"
    [guo_cluster,guo_row_names,guo_column_names] = F_tsv_import(fn);
        
    [~,i,j] = intersect(Guo_cn,guo_row_names,'stable');
    if sum(i-j)==0
    guo_cluster = guo_cluster(:,3);
    c = guo_cluster;
    guo_cluster(find(c==11))= 9;
    guo_cluster(find(c==9)) = 10;
    guo_cluster(find(c==10))= 11;
    end
    guo_cluster_name = {'SSCs','Sg','Sc-1','Sc-2','RS','ES','Sperm-1','Sperm-2',...
                        'PMC','Macrophage','Endothelia','Sertoli','Ledig'};
%Alread saved: Guo_sperm,Guo_cn,Guo_gn,guo_cluster,guo_row_names,guo_column_names, guo_cluster_name
load('YOUR_PATH\Other scRNA-seq data\Guo et al, CellRes, 2018, human spermatogenesis\Guo2018_human_spermatogenesis.mat')
    





%#######Compare with our study    
    %Plot corrcoef of gene clusters
    load('YOUR_PATH\h_m_stage_expression.mat');
    load('YOUR_PATH\human_TESE_Bioreps_20181031.mat');
    [~,x,y] = intersect(human_genename,Guo_gn,'stable');
    A   =   F_normalize(human_tese(x,:),'median');
    for i= 1:19
        h_tese_cluster_expression(:,i) = mean(A(:,find(h_germ_som==i))')';
        h_tese_cluster_expression_UMI(:,i) = mean(human_tese(x,find(h_germ_som==i))')';
    end
    A   =   F_normalize(Guo_sperm(y,:),'median');
    for i= 1:13
        guo_cluster_expression(:,i) = mean(A(:,find(guo_cluster==i))')';
        guo_cluster_expression_UMI(:,i) = mean(Guo_sperm(y,find(guo_cluster==i))')';
    end

    H_guo = [h_tese_cluster_expression_UMI guo_cluster_expression_UMI];
    H_guo_corr=corrcoef(F_normalize(H_guo,'FT'));
    figure;
    imagesc(flipud(H_guo_corr(20:32,1:19)),[0.5 0.95]);
    colormap(RdBu_cm); axis equal; set(gca,'color','none');ylim([0.5 13.5])
    set(gca,'xtick',1:19); set(gca,'ytick',1:13); 
    xticklabels({'Sg-1','Sg-2','Sc-1','Sc-2','Sc-3','Sc-4','RS-1','RS-2','RS-3','RS-4',...
           'ES-1','ES-2','ES-3','ES-4','PMC','Macrophage','Endothelia','Sertoli','Ledig','Other'});
    xtickangle(90)
    yticklabels(flip(guo_cluster_name));
    xlabel('Human spermatogenic cell clusters');
    ylabel('Guo et al, 2018 cell clusters');
    title('Human-Guo2018 Correlation-UMI');

    %Plot mutation rates across gene clusters.
    %10X Genomics only select the expressed genes, thus we trace back the
    %full gene list of 19831 genes, as used by our study. The missing genes
    %are pre-filled with zero.
    [~,x,y] = intersect(human_genename,Guo_gn,'stable');
    guo_germ(1:19831,1:length(find(guo_cluster<=8))) = 0; 
    guo_germ(x,:)= Guo_sperm(y,find(guo_cluster<=8));
    guo_germ_cluster = guo_cluster(find(guo_cluster<=8));
    %cluster genes with k=5 which is optimal
    rng(1)
    [guo_germ_cluexp,guo_germ_cluexpUMI,guo_germ_cluindex]=F_cluster_gene(guo_germ,guo_germ_cluster,5,...
                                                                    'expression_cutoff_n',5,...
                                                                    'umi_ratio_cutoff',0,...
                                                                    'include_unexpressed_gene','yes',...
                                                                    'ylimit',[-2 3]);
    % reorder the gene cluster by visulization of expression stage
    guo_gene_cluster = guo_germ_cluindex;
    guo_gene_cluster(find(guo_germ_cluindex==1))=1;
    guo_gene_cluster(find(guo_germ_cluindex==2))=2;
    guo_gene_cluster(find(guo_germ_cluindex==3))=4;
    guo_gene_cluster(find(guo_germ_cluindex==4))=3;
    guo_gene_cluster(find(guo_germ_cluindex==5))=5;
    %human gene body plot using the Guo2018 gene clusters
    human_gene_info = struct2table(tdfread('YOUR_PATH\human_genes.tsv'));
    [~,~,xi] = intersect(human_genename,cellstr(human_gene_info.name),'stable');
    h_gname  = cellstr(human_gene_info.ensembl_id(xi,:));
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
     human_genebody_guo = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     guo_gene_cluster,...
                                     h_gname,...
                                     'F_compare_variants','yes',...
                                     'FaceColor',genecluster_edgecm,...
                                     'Title_prefix','Human-genebody-by-Guo2018'); 





%plot germ cell gene number
    [~,i,j] = intersect(human_genename,Guo_gn,'stable');
    A = Guo_sperm(j,:);  A(A>0)=1; 
    figure; violinplot(sum(A),guo_cluster)
    box off; set(gca,'color','none');
    ylabel('Gene number'); xlim([0.1 13.9]);
    xticklabels(guo_cluster_name)
    xtickangle(45)
    title('Human-Guo2018 gene number by clusters');
    
    %Count germ cell gene numbre
    [~,i,j] = intersect(human_genename,Guo_gn,'stable');
    A = Guo_sperm(j,:);  A(A>0)=1; 
    c=guo_cluster; 
    n=length(unique(c)); 
    clear gene_num gene_num2 temp
    for i =1:n
        B = A(:,find(c==i));
        C = sum(B);
        gene_num(i,1) = size(B,2);
        gene_num(i,2) = round(mean(C));
        gene_num(i,3) = round(prctile(C,5));
        gene_num(i,4) = round(prctile(C,25));
        gene_num(i,5) = round(median(C));
        gene_num(i,6) = round(prctile(C,75));
        gene_num(i,7) = round(prctile(C,95));
        if size(B,2)<=200
        gene_num(i,8) = length(find(sum(B,2)>=2));
        else
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:gene_num(i,1),200,'Replace',false)),2)>=2));
        end
        gene_num(i,8) = round(mean(temp));
        end
        %bootstrap 100 times.
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:gene_num(i,1),1000,'Replace',true)),2)>=5));
        end
        gene_num(i,9) = round(mean(temp));
    end
        B = A(:,find(c<9));
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:size(B,2),1000,'Replace',true)),2)>=5));
        end
        gene_num2(1) = round(mean(temp))
        gene_num2(2) = round(std(temp));
    
        
        
% binarize the soma gene clustering as germexp(17941 genes):(somaexp, somaunexp)
c = 1;
A = Guo_sperm(:,find(guo_cluster>8)); A(A>0)=1; 
gn = Guo_gn;
input_gene = human_genename(find(h_gene_cluster>0));
input_gene_id = h_gname(find(h_gene_cluster>0));
input_gene_clu(1:length(input_gene_id)) = 0;
[~,~,i] = intersect(gn(find(sum(A')>=5)),input_gene,'stable');
input_gene_clu(i) = 1;
    %figure;
    subplot(1,5,c)
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
    F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_clu,...
                                     input_gene_id,...
                                     'F_compare_variants','only',...
                                     'FaceColor',exunex_cm,...
                                     'Title_prefix','Genebody'); 
    ylim([0 250]);  title('Guo2018');


% binarize the soma gene and plot germunexp/germexp:(somaexp, somaunexp)
c = 1;
A = Guo_sperm(:,find(guo_cluster>8)); A(A>0)=1; 
gn = Guo_gn;
input_gene = human_genename;
input_gene_id = h_gname;
soma_gene_cluster(1:length(input_gene_id)) = 0;%initialize the gene indeces.
min_cell_cutoff = 5
[~,~,i] = intersect(gn(find(sum(A')>=min_cell_cutoff)),input_gene,'stable');
soma_gene_cluster(i) = 1;

input_gene_cluster = soma_gene_cluster;
input_gene_cluster(find(h_gene_cluster==0 & soma_gene_cluster==0)) = 0; %germunexp-somaunexp
input_gene_cluster(find(h_gene_cluster==0 & soma_gene_cluster >0)) = 1; %germunexp-somaexp
input_gene_cluster(find(h_gene_cluster >0 & soma_gene_cluster==0)) = 2; %germexp-somaunexp
input_gene_cluster(find(h_gene_cluster >0 & soma_gene_cluster >0)) = 3; %germexp-somaexp
gene_count2(c,1:4) = histcounts(input_gene_cluster)
if gene_count2(c,2)<2
    input_gene_cluster(find(input_gene_cluster>0))=input_gene_cluster(find(input_gene_cluster>0))-1;
    histcounts(input_gene_cluster)
end
    
   % Plot human genebody variants
   %figure;
    subplot(1,5,c)
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
    soma_external(c).var = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_cluster,...
                                     input_gene_id,...
                                     'F_compare_variants','only',...
                                     'FaceColor',soma_cm,...
                                     'Title_prefix','Guo2018'); ylim([0 300]);
    %Plot asymmetry scores of unexpressed genes
    fn1='YOUR_PATH\human_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.tsv'
    soma_external(c).var  = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     soma_gene_cluster,...
                                     input_gene_id,...
                                     'F_compare_variants','no',...
                                     'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Guo2018');     

                                 
                                 

%% 2. Pelline et al, NatComm, 2019, human hematopoiesis
 fn='YOUR_PATH\Other scRNA-seq data\Pelline et al, NatComm, 2019, human hematopoiesis\GSM3305359_HSC.raw_counts.tsv'
    [PellineHSC,PellineHSC_cn,PellineHSC_gn] = F_tsv_import(fn);   
 fn='YOUR_PATH\Other scRNA-seq data\Pelline et al, NatComm, 2019, human hematopoiesis\GSM3305360_MPP.raw_counts.tsv'
    [PellineMPP,PellineMPP_cn,PellineMPP_gn] = F_tsv_import(fn);   
 fn='YOUR_PATH\Other scRNA-seq data\Pelline et al, NatComm, 2019, human hematopoiesis\GSM3305361_MLP.raw_counts.tsv'
    [PellineMLP,PellineMLP_cn,PellineMLP_gn] = F_tsv_import(fn);   
 fn='YOUR_PATH\Other scRNA-seq data\Pelline et al, NatComm, 2019, human hematopoiesis\GSM3305362_PreBNK.raw_counts.tsv'
    [PellinePreBNK,PellinePreBNK_cn,PellinePreBNK_gn] = F_tsv_import(fn);   
 fn='YOUR_PATH\Other scRNA-seq data\Pelline et al, NatComm, 2019, human hematopoiesis\GSM3305363_MEP.raw_counts.tsv'
    [PellineMEP,PellineMEP_cn,PellineMEP_gn] = F_tsv_import(fn);  
 fn='YOUR_PATH\Other scRNA-seq data\Pelline et al, NatComm, 2019, human hematopoiesis\GSM3305364_CMP.raw_counts.tsv'
    [PellineCMP,PellineCMP_cn,PellineCMP_gn] = F_tsv_import(fn);  
 fn='YOUR_PATH\Other scRNA-seq data\Pelline et al, NatComm, 2019, human hematopoiesis\GSM3305365_GMP.raw_counts.tsv'
    [PellineGMP,PellineGMP_cn,PellineGMP_gn] = F_tsv_import(fn);  
 fn='YOUR_PATH\Other scRNA-seq data\Pelline et al, NatComm, 2019, human hematopoiesis\GSM3305366_LinNegCD34PosCD164Pos.raw_counts.tsv'
    [CD34PosCD164Pos,CD34PosCD164Pos_cn,CD34PosCD164Pos_gn] = F_tsv_import(fn);   
 fn='YOUR_PATH\Other scRNA-seq data\Pelline et al, NatComm, 2019, human hematopoiesis\GSM3305367_LinNegCD34NegCD164high.raw_counts.tsv'
    [CD34NegCD164high,CD34NegCD164high_cn,CD34NegCD164high_gn] = F_tsv_import(fn);  
 fn='YOUR_PATH\Other scRNA-seq data\Pelline et al, NatComm, 2019, human hematopoiesis\GSM3305368_LinNegCD34lowCD164high.raw_counts.tsv'
    [CD34lowCD164high,CD34lowCD164high_cn,CD34lowCD164high_gn] = F_tsv_import(fn);  
 fn='YOUR_PATH\Other scRNA-seq data\Pelline et al, NatComm, 2019, human hematopoiesis\GSM3305369_LinNegCD34NegCD164low.raw_counts.tsv'
    [CD34NegCD164low,CD34NegCD164low_cn,CD34NegCD164low_gn] = F_tsv_import(fn);  
    
    [~,a,b] = intersect(human_genename,PellineHSC_gn,'stable');
    [~,c,d] = intersect(human_genename,CD34NegCD164low_gn,'stable');
    [h_gm_idx,e,f] = intersect(a,c,'stable');
    Pelline_gn = human_genename(h_gm_idx);
Pelline = [PellineHSC((b(e)),:) PellineMPP((b(e)),:) PellineMLP((b(e)),:) PellinePreBNK((b(e)),:) ...
           PellineMEP((b(e)),:) PellineCMP((b(e)),:) PellineGMP((b(e)),:) ...
           CD34PosCD164Pos(d(f),:) CD34NegCD164high(d(f),:) CD34lowCD164high(d(f),:) CD34NegCD164low(d(f),:)];
Pelline_cellclu = [1*ones(size(PellineHSC,2),1)' ...
                   2*ones(size(PellineMPP,2),1)' ...
                   3*ones(size(PellineMLP,2),1)' ...
                   4*ones(size(PellinePreBNK,2),1)' ...
                   5*ones(size(PellineMEP,2),1)' ...
                   6*ones(size(PellineCMP,2),1)' ...
                   7*ones(size(PellineGMP,2),1)' ...
                   8*ones(size(CD34PosCD164Pos,2),1)' ...
                   9*ones(size(CD34NegCD164high,2),1)' ...
                   10*ones(size(CD34lowCD164high,2),1)' ...
                   11*ones(size(CD34NegCD164low,2),1)'];
%Alread saved: Pelline,Pelline_cellclu,Pelline_gn
load('YOUR_PATH\Other scRNA-seq data\Pelline et al, NatComm, 2019, human hematopoiesis\Pelline_human_hematopoiesis.mat')


    A = Pelline; B=A; A(A>0)=1; 
    length(find(sum(A,2)>(0.005*size(A,2))))
    length(find(sum(A(:,datasample(1:21412,2552,'Replace',false)),2)>5))
    c = Pelline_cellclu;  
    figure; violinplot(sum(A),c,'ShowData',false)
    box off; set(gca,'color','none');
    ylabel('Gene number'); xlim([0.1 11.9])
    xticklabels({'HSC','MPP','MLP','preB/NK','MEP','CMP','GMP',...
                 'CD34PosCD164Pos','CD34NegCD164high','CD34lowCD164high','CD34NegCD164low'});
    xtickangle(45)
    title('Human hematopoiesis (Pelline et al, NatComm, 2019)')

% Count gene numbsers
    A = Pelline; A(A>0)=1; 
    c=Pelline_cellclu;
    n=length(unique(c)); 
    clear gene_num gene_num2
    for i =1:n
        B = A(:,find(c==i));
        C = sum(B);
        gene_num(i,1) = size(B,2);
        gene_num(i,2) = round(mean(C));
        gene_num(i,3) = round(prctile(C,5));
        gene_num(i,4) = round(prctile(C,25));
        gene_num(i,5) = round(median(C));
        gene_num(i,6) = round(prctile(C,75));
        gene_num(i,7) = round(prctile(C,95));
        if size(B,2)<=200
        gene_num(i,8) = length(find(sum(B,2)>=2));
        else
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:gene_num(i,1),200,'Replace',false)),2)>=2));
        end
        gene_num(i,8) = round(mean(temp));
        end
       %bootstrap 100 times.
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:gene_num(i,1),1000,'Replace',true)),2)>=5));
        end
        gene_num(i,9) = round(mean(temp));
    end
        B = A;
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:size(B,2),1000,'Replace',true)),2)>=5));
        end
        gene_num2(1) = round(mean(temp))

% binarize the gene clustering as germexp(17941 genes):(somaexp, somaunexp)
A = Pelline; A(A>0)=1; 
gn = Pelline_gn;
input_gene = human_genename(find(h_gene_cluster>0));
input_gene_id = h_gname(find(h_gene_cluster>0));
input_gene_clu(1:length(input_gene_id)) = 0;
[~,~,i] = intersect(gn(find(sum(A')>=5)),input_gene,'stable');
input_gene_clu(i) = 1;
    %figure;
    subplot(1,5,2)
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
    F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_clu,...
                                     input_gene_id,...
                                     'F_compare_variants','only',...
                                     'FaceColor',exunex_cm,...
                                     'Title_prefix','Genebody'); 
    ylim([0 250]);
    title('Pelline2019')                 
     

    
    
% binarize the soma gene and plot germunexp/germexp:(somaexp, somaunexp)
c = 2;
A = Pelline; A(A>0)=1; 
gn = Pelline_gn;
input_gene = human_genename;
input_gene_id = h_gname;
soma_gene_cluster(1:length(input_gene_id)) = 0;%initialize the gene indeces.
min_cell_cutoff = 5
[~,~,i] = intersect(gn(find(sum(A')>=min_cell_cutoff)),input_gene,'stable');
soma_gene_cluster(i) = 1;

input_gene_cluster = soma_gene_cluster;
input_gene_cluster(find(h_gene_cluster==0 & soma_gene_cluster==0)) = 0; %germunexp-somaunexp
input_gene_cluster(find(h_gene_cluster==0 & soma_gene_cluster >0)) = 1; %germunexp-somaexp
input_gene_cluster(find(h_gene_cluster >0 & soma_gene_cluster==0)) = 2; %germexp-somaunexp
input_gene_cluster(find(h_gene_cluster >0 & soma_gene_cluster >0)) = 3; %germexp-somaexp
gene_count2(c,1:4) = histcounts(input_gene_cluster)
if gene_count2(c,2)<2
    input_gene_cluster(find(input_gene_cluster>0))=input_gene_cluster(find(input_gene_cluster>0))-1;
    histcounts(input_gene_cluster)
end
    
   % Plot human genebody variants
   %figure;
    subplot(1,5,c)
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
    soma_external(c).var = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_cluster,...
                                     input_gene_id,...
                                     'F_compare_variants','only',...
                                     'FaceColor',soma_cm,...
                                     'Title_prefix','Pelline2019'); ylim([0 300]);
    %Plot asymmetry scores of unexpressed genes
    fn1='YOUR_PATH\human_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.tsv'
    soma_external(c).var  = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     soma_gene_cluster,...
                                     input_gene_id,...
                                     'F_compare_variants','no',...
                                     'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Pelline2019');  
      
                                 


    
%% 3. Le Manno et al, Cell, 2016, Human developing midbrain week7-week10
    %Setup input files
    fn='YOUR_PATH\Other scRNA-seq data\Le Manno et al, Cell, 2016, human developing midbrain\GSE76381_EmbryoMoleculeCounts.cef_clean.txt'
    midbrain=importdata(fn);
    midbrain_cn = split(midbrain.textdata(1));
    midbrain_cn = midbrain_cn(2:end);
    midbrain_gn = midbrain.textdata(2:end);
    midbrain=midbrain.data(:,2:end);
    [~,a,b] = intersect(human_genename,midbrain_gn,'stable');
    midbrain = midbrain(b,:);  
    midbrain_gn = midbrain_gn(b,:);

    midbrain_uniq_cn = unique(midbrain_cn);
    midbrain_clu(1:length(midbrain_cn)) = 0;
    for i=1:length(midbrain_uniq_cn)
        idx = strmatch(midbrain_uniq_cn(i),midbrain_cn,'exact');
        midbrain_clu(idx) = i;
    end
    histcounts(midbrain_clu)
    
%Alread saved: midbrain,midbrain_cn,midbrain_gn,midbrain_uniq_cn,midbrain_clu
load('YOUR_PATH\Other scRNA-seq data\Le Manno et al, Cell, 2016, human developing midbrain\LeManno_human_developing_midbrain.mat')
    
    %Plot protein coding gene number across cell types.
    A = midbrain;  
    A(A>0)=1; length(find(sum(A,2)>(0.005*size(A,2))))
    c = midbrain_clu;  
    figure; violinplot(sum(A),c,'ShowData',false)
    box off; set(gca,'color','none');
    ylabel('Gene number'); xlim([0.1 26.9]);ylim([0 6000])
    xticklabels(midbrain_uniq_cn);
    xtickangle(45)
    title('Human developing midbrain(Le Manno et al, Cell, 2016)')
        
% Count gene numbsers
    A = midbrain; A(A>0)=1; 
    c=midbrain_clu;
    n=length(unique(c)); 
    clear gene_num gene_num2
    for i =1:n
        B = A(:,find(c==i));
        C = sum(B);
        gene_num(i,1) = size(B,2);
        gene_num(i,2) = round(mean(C));
        gene_num(i,3) = round(prctile(C,5));
        gene_num(i,4) = round(prctile(C,25));
        gene_num(i,5) = round(median(C));
        gene_num(i,6) = round(prctile(C,75));
        gene_num(i,7) = round(prctile(C,95));
        if size(B,2)<=200
        gene_num(i,8) = length(find(sum(B,2)>=2));
        else
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:gene_num(i,1),200,'Replace',false)),2)>=2));
        end
        gene_num(i,8) = round(mean(temp));
        end
        %bootstrap 100 times.
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:gene_num(i,1),1000,'Replace',true)),2)>=5));
        end
        gene_num(i,9) = round(mean(temp));
    end
        B = A(:,find(c>1));
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:size(B,2),1000,'Replace',true)),2)>=5));
        end
        gene_num2(1) = round(mean(temp))   
    
% binarize the gene clustering as germexp(17941 genes):(somaexp, somaunexp)
A = midbrain; A(A>0)=1; 
gn = midbrain_gn;
input_gene = human_genename(find(h_gene_cluster>0));
input_gene_id = h_gname(find(h_gene_cluster>0));
input_gene_clu(1:length(input_gene_id)) = 0;
[~,~,i] = intersect(gn(find(sum(A')>=5)),input_gene,'stable');
input_gene_clu(i) = 1;
    %figure;
    subplot(1,5,3)
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
    F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_clu,...
                                     input_gene_id,...
                                     'F_compare_variants','only',...
                                     'FaceColor',exunex_cm,...
                                     'Title_prefix','Genebody'); 
    ylim([0 250]);
    title('Le Manno2016')
                

    
    
% binarize the soma gene and plot germunexp/germexp:(somaexp, somaunexp)
c = 3;
A = midbrain; A(A>0)=1; 
gn = midbrain_gn;
input_gene = human_genename;
input_gene_id = h_gname;
soma_gene_cluster(1:length(input_gene_id)) = 0;%initialize the gene indeces.
min_cell_cutoff = 5
[~,~,i] = intersect(gn(find(sum(A')>=min_cell_cutoff)),input_gene,'stable');
soma_gene_cluster(i) = 1;

input_gene_cluster = soma_gene_cluster;
input_gene_cluster(find(h_gene_cluster==0 & soma_gene_cluster==0)) = 0; %germunexp-somaunexp
input_gene_cluster(find(h_gene_cluster==0 & soma_gene_cluster >0)) = 1; %germunexp-somaexp
input_gene_cluster(find(h_gene_cluster >0 & soma_gene_cluster==0)) = 2; %germexp-somaunexp
input_gene_cluster(find(h_gene_cluster >0 & soma_gene_cluster >0)) = 3; %germexp-somaexp
gene_count2(c,1:4) = histcounts(input_gene_cluster)
if gene_count2(c,2)<2
    input_gene_cluster(find(input_gene_cluster>0))=input_gene_cluster(find(input_gene_cluster>0))-1;
    histcounts(input_gene_cluster)
end
    
   % Plot human genebody variants
   %figure;
    subplot(1,5,c)
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
    soma_external(c).var = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_cluster,...
                                     input_gene_id,...
                                     'F_compare_variants','only',...
                                     'FaceColor',soma_cm,...
                                     'Title_prefix','Le Manno2016'); ylim([0 300]);
    %Plot asymmetry scores of unexpressed genes
    fn1='YOUR_PATH\human_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.tsv'
    soma_external(c).var  = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     soma_gene_cluster,...
                                     input_gene_id,...
                                     'F_compare_variants','no',...
                                     'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Le Manno2016');  
                                 

                                 
                                 
%% 4. Nowakowski et al, Science, 2017, human cortex
    %Setup input files
fn='YOUR_PATH\Other scRNA-seq data\Nowakowski et al, Science, 2017, human cortex\exprMatrix.tsv'
[Nowakowski,Nowakowski_cn,Nowakowski_gn] = F_tsv_import(fn);     
fn='YOUR_PATH\Other scRNA-seq data\Nowakowski et al, Science, 2017, human cortex\aap8809_Nowakowski_SM-Tables-S1-S11.xlsx'
Nowakowski_sample = readtable(fn,'Sheet','Table3 - CellMetadata');
    [~,a,b] = intersect(human_genename,Nowakowski_gn,'stable');
    [~,c,d] = intersect(cellstr(Nowakowski_sample.Cell),Nowakowski_cn,'stable');
Nowakowski = Nowakowski(b,:);
Nowakowski_gn = Nowakowski_gn(b);
Nowakowski_cn = Nowakowski_cn(d);
Nowakowski_clu_cn = cellstr(Nowakowski_sample.WGCNAcluster(c,:));
%Concatenate similar cell clusters
Nowakowski_clu_cn2 = regexprep(Nowakowski_clu_cn,'[\d"]','');

    Nowakowski_clu_uniq_cn2 = unique(Nowakowski_clu_cn2);
    Nowakowski_clu2(1:length(Nowakowski_clu_cn2)) = 0;
    for i=1:length(Nowakowski_clu_uniq_cn2)
        idx = strmatch(Nowakowski_clu_uniq_cn2(i),Nowakowski_clu_cn2,'exact');
        Nowakowski_clu2(idx) = i;
    end
    histcounts(Nowakowski_clu2)
    
%Alread saved: Nowakowski,Nowakowski_gn,Nowakowski_cn,Nowakowski_clu_cn,...
%Nowakowski_clu_uniq_cn,Nowakowski_clu,Nowakowski_clu_cn2,Nowakowski_clu_uniq_cn2,Nowakowski_clu2
load('YOUR_PATH\Other scRNA-seq data\Nowakowski et al, Science, 2017, human cortex\Nowakowski_human_cortex.mat')
    
    %Plot protein coding gene number across cell types.
    A = Nowakowski;  
    A(A>0)=1; length(find(sum(A,2)>(0.005*size(A,2))))
    c = Nowakowski_clu2;  
    figure; violinplot(sum(A),c,'ShowData',false)
    box off; set(gca,'color','none');
    ylabel('Gene number'); xlim([0.1 27.9]);ylim([0 8000])
    xticklabels(Nowakowski_clu_uniq_cn2);
    xtickangle(45)
    title('Human  cortex (Nowakowski et al, Science, 2017)')
        
% Count gene numbsers
    A = Nowakowski; A(A>0)=1; 
    c=Nowakowski_clu2;
    n=length(unique(c)); 
    clear gene_num gene_num2
    for i =1:n
        B = A(:,find(c==i));
        C = sum(B);
        gene_num(i,1) = size(B,2);
        gene_num(i,2) = round(mean(C));
        gene_num(i,3) = round(prctile(C,5));
        gene_num(i,4) = round(prctile(C,25));
        gene_num(i,5) = round(median(C));
        gene_num(i,6) = round(prctile(C,75));
        gene_num(i,7) = round(prctile(C,95));
        if size(B,2)<=200
        gene_num(i,8) = length(find(sum(B,2)>=2));
        else
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:gene_num(i,1),200,'Replace',false)),2)>=2));
        end
        gene_num(i,8) = round(mean(temp));
        end
       %bootstrap 100 times.
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:gene_num(i,1),1000,'Replace',true)),2)>=5));
        end
        gene_num(i,9) = round(mean(temp));
    end
        B = A(:,find(c>1));
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:size(B,2),1000,'Replace',true)),2)>=5));
        end
        gene_num2(1) = round(mean(temp))   
    
% binarize the gene clustering as germexp(17941 genes):(somaexp, somaunexp)
A = Nowakowski; A(A>0)=1; 
gn = Nowakowski_gn;
input_gene = human_genename(find(h_gene_cluster>0));
input_gene_id = h_gname(find(h_gene_cluster>0));
input_gene_clu(1:length(input_gene_id)) = 0;
[~,~,i] = intersect(gn(find(sum(A')>=5)),input_gene,'stable');
input_gene_clu(i) = 1;
    %figure;
    subplot(1,5,4)
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
    F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_clu,...
                                     input_gene_id,...
                                     'F_compare_variants','only',...
                                     'FaceColor',exunex_cm,...
                                     'Title_prefix','Genebody'); 
    ylim([0 250]);
    title('Nowakowski2017')
        
        

    
% binarize the soma gene and plot germunexp/germexp:(somaexp, somaunexp)
c = 4;
A = Nowakowski; A(A>0)=1; 
gn = Nowakowski_gn;
input_gene = human_genename;
input_gene_id = h_gname;
soma_gene_cluster(1:length(input_gene_id)) = 0;%initialize the gene indeces.
min_cell_cutoff = 5
[~,~,i] = intersect(gn(find(sum(A')>=min_cell_cutoff)),input_gene,'stable');
soma_gene_cluster(i) = 1;

input_gene_cluster = soma_gene_cluster;
input_gene_cluster(find(h_gene_cluster==0 & soma_gene_cluster==0)) = 0; %germunexp-somaunexp
input_gene_cluster(find(h_gene_cluster==0 & soma_gene_cluster >0)) = 1; %germunexp-somaexp
input_gene_cluster(find(h_gene_cluster >0 & soma_gene_cluster==0)) = 2; %germexp-somaunexp
input_gene_cluster(find(h_gene_cluster >0 & soma_gene_cluster >0)) = 3; %germexp-somaexp
gene_count2(c,1:4) = histcounts(input_gene_cluster)
if gene_count2(c,2)<2
    input_gene_cluster(find(input_gene_cluster>0))=input_gene_cluster(find(input_gene_cluster>0))-1;
    histcounts(input_gene_cluster)
end
    
   % Plot human genebody variants
   %figure;
    subplot(1,5,c)
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
    soma_external(c).var = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_cluster,...
                                     input_gene_id,...
                                     'F_compare_variants','only',...
                                     'FaceColor',soma_cm,...
                                     'Title_prefix','Nowakowski2017'); ylim([0 300]);
    %Plot asymmetry scores of unexpressed genes
    fn1='YOUR_PATH\human_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.tsv'
    soma_external(c).var  = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     soma_gene_cluster,...
                                     input_gene_id,...
                                     'F_compare_variants','no',...
                                     'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Nowakowski2017');  
    
        
%% 5. Hochane et al, PlosBio, 2019, human fetal kidney
fn='YOUR_PATH\Other scRNA-seq data\Hochane et al, PlosBio, 2019, human fetal kidney\w16.tsv'
Hochane_emw16_kidney = dlmread(fn,'\t');
fn='YOUR_PATH\Other scRNA-seq data\Hochane et al, PlosBio, 2019, human fetal kidney\w16_genes.tsv'
Hochane_w16_gn_raw = importdata(fn);
for i=1:33694
    temp = split(char(Hochane_w16_gn_raw(i,:)));
Hochane_w16_gn(i,1:2) = temp;
end
fn='YOUR_PATH\Other scRNA-seq data\Hochane et al, PlosBio, 2019, human fetal kidney\w16_barcodes.tsv'
Hochane_w16_barcode = importdata(fn);
fn='YOUR_PATH\Other scRNA-seq data\Hochane et al, PlosBio, 2019, human fetal kidney\w16_barcodes_celltypes.csv'
Hochane_w16_sample = readtable(fn);

    [~,a,b] = intersect(human_genename,Hochane_w16_gn(:,2),'stable');
    [~,c,d] = intersect(Hochane_w16_sample.cell_barcode,Hochane_w16_barcode,'stable');
Hochane_emw16_kidney = Hochane_emw16_kidney(b,d);
Hochane_w16_gn = Hochane_w16_gn(b);
Hochane_w16_cn = Hochane_w16_sample.cell_type(c);

    Hochane_w16_uniq_cn = unique(Hochane_w16_cn);
    Hochane_w16_clu(1:length(Hochane_w16_cn)) = 0;
    for i=1:length(Hochane_w16_uniq_cn)
        idx = strmatch(Hochane_w16_uniq_cn(i),Hochane_w16_cn,'exact');
        Hochane_w16_clu(idx) = i;
    end
    histcounts(Hochane_w16_clu)

%Alread saved:
%Hochane_emw16_kidney,Hochane_w16_gn,Hochane_w16_cn,Hochane_w16_uniq_cn,Hochane_w16_clu
load('YOUR_PATH\Other scRNA-seq data\Hochane et al, PlosBio, 2019, human fetal kidney\Hochane_human_fetal_kidney.mat')


    A = Hochane_emw16_kidney; A(A>0)=1; 
    length(find(sum(A,2)>(0.005*size(A,2))))
    length(find(sum(A(:,datasample(1:8500,500,'Replace',false)),2)>5))
    c = Hochane_w16_clu;  
    figure; violinplot(sum(A),c,'ShowData',false)
    box off; set(gca,'color','none');
    ylabel('Gene number'); xlim([0.1 22.9])
    xticklabels(Hochane_w16_uniq_cn);
    xtickangle(45)
    title('Human embryonic kidney (Hochane et al,PlosBio,2019) ')

% Count gene numbsers
    A = Hochane_emw16_kidney; A(A>0)=1; 
    c=Hochane_w16_clu;    
    n=length(unique(c)); 
    clear gene_num gene_num2 temp
    for i =1:n
        B = A(:,find(c==i));
        C = sum(B);
        gene_num(i,1) = size(B,2);
        gene_num(i,2) = round(mean(C));
        gene_num(i,3) = round(prctile(C,5));
        gene_num(i,4) = round(prctile(C,25));
        gene_num(i,5) = round(median(C));
        gene_num(i,6) = round(prctile(C,75));
        gene_num(i,7) = round(prctile(C,95));
        if size(B,2)<=200
        gene_num(i,8) = length(find(sum(B,2)>=2));
        else
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:gene_num(i,1),200,'Replace',false)),2)>=2));
        end
        gene_num(i,8) = round(mean(temp));
        end
       %bootstrap 100 times.
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:gene_num(i,1),1000,'Replace',true)),2)>=5));
        end
        gene_num(i,9) = round(mean(temp));
    end
        B = A;
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:size(B,2),1000,'Replace',true)),2)>=5));
        end
        gene_num2(1) = round(mean(temp))

% binarize the gene clustering as germexp(17941 genes):(somaexp, somaunexp)
A = Hochane_emw16_kidney; A(A>0)=1; 
gn = Hochane_w16_gn;
input_gene = human_genename(find(h_gene_cluster>0));
input_gene_id = h_gname(find(h_gene_cluster>0));
input_gene_clu(1:length(input_gene_id)) = 0;
[~,~,i] = intersect(gn(find(sum(A')>=5)),input_gene_id,'stable');
input_gene_clu(i) = 1;
    %figure;
    subplot(1,5,5)
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
    F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_clu,...
                                     input_gene_id,...
                                     'F_compare_variants','only',...
                                     'FaceColor',exunex_cm,...
                                     'Title_prefix','Genebody'); 
    ylim([0 250]);
    title('Hochane2018')           
        
           
 
    
% binarize the soma gene and plot germunexp/germexp:(somaexp, somaunexp)
c = 5;
A = Hochane_emw16_kidney; A(A>0)=1; 
gn = Hochane_w16_gn;
input_gene = human_genename;
input_gene_id = h_gname;
soma_gene_cluster(1:length(input_gene_id)) = 0;%initialize the gene indeces.
min_cell_cutoff = 5
[~,~,i] = intersect(gn(find(sum(A')>=min_cell_cutoff)),input_gene_id,'stable');
soma_gene_cluster(i) = 1;

input_gene_cluster = soma_gene_cluster;
input_gene_cluster(find(h_gene_cluster==0 & soma_gene_cluster==0)) = 0; %germunexp-somaunexp
input_gene_cluster(find(h_gene_cluster==0 & soma_gene_cluster >0)) = 1; %germunexp-somaexp
input_gene_cluster(find(h_gene_cluster >0 & soma_gene_cluster==0)) = 2; %germexp-somaunexp
input_gene_cluster(find(h_gene_cluster >0 & soma_gene_cluster >0)) = 3; %germexp-somaexp
gene_count2(c,1:4) = histcounts(input_gene_cluster)
if gene_count2(c,2)<2
    input_gene_cluster(find(input_gene_cluster>0))=input_gene_cluster(find(input_gene_cluster>0))-1;
    histcounts(input_gene_cluster)
end
    
   % Plot human genebody variants
   %figure;
    subplot(1,5,c)
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
    soma_external(c).var = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_cluster,...
                                     input_gene_id,...
                                     'F_compare_variants','only',...
                                     'FaceColor',soma_cm,...
                                     'Title_prefix','Hochane2018'); ylim([0 300]);
    %Plot asymmetry scores of unexpressed genes
    fn1='YOUR_PATH\human_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.intron.singleFreq.tsv'
    soma_external(c).var  = F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     soma_gene_cluster,...
                                     input_gene_id,...
                                     'F_compare_variants','no',...
                                     'FaceColor',genecluster_facecm,...
                                     'Title_prefix','Hochane2018');      

    
%% [Not included!SMART-seq] Zhong et al, Nature, 2018, human developing prefrontal cortex
 fn='YOUR_PATH\Other scRNA-seq data\Zhong et al, Nature, 2018, human developing prefrontal cortex\GSE104276_readme_sample_barcode.xlsx'
Zhong_PFC_sample = readtable(fn,'Sheet','SampleInfo');
 fn='YOUR_PATH\Other scRNA-seq data\Zhong et al, Nature, 2018, human developing prefrontal cortex\GSE104276_all_pfc_2394_UMI_count_NOERCC.txt'
[Zhong_PFC,Zhong_PFC_cn,Zhong_PFC_gn] = F_tsv_import(fn);   
    [~,a,b] = intersect(human_genename,Zhong_PFC_gn,'stable');
    [~,c,d] = intersect(cellstr(Zhong_PFC_sample.Var1),Zhong_PFC_cn,'stable');
Zhong_PFC = Zhong_PFC(b,d);
Zhong_PFC_gn = Zhong_PFC_gn(b);
Zhong_PFC_cn = Zhong_PFC_cn(d);

    Zhong_PFC_uniq_cn = unique(cellstr(Zhong_PFC_sample.cell_types(c,:)));
    Zhong_PFC_clu(1:length(Zhong_PFC_cn)) = 0;
    for i=1:length(Zhong_PFC_uniq_cn)
        idx = strmatch(Zhong_PFC_uniq_cn(i),cellstr(Zhong_PFC_sample.cell_types(c,:)),'exact');
        Zhong_PFC_clu(idx) = i;
    end
    histcounts(Zhong_PFC_clu)
    
%Alread saved: Zhong_PFC,Zhong_PFC_gn,Zhong_PFC_cn,Zhong_PFC_uniq_cn,Zhong_PFC_clu,Zhong_PFC_sample
load('YOUR_PATH\Other scRNA-seq data\Zhong et al, Nature, 2018, human developing prefrontal cortex\Zhong_human_developing_prefrontal_cortex.mat')

    %Plot protein coding gene number across cell types.
    A = Zhong_PFC; 
    A(A>0)=1; length(find(sum(A,2)>(0.005*size(A,2))))
    c = Zhong_PFC_clu;  
    figure; violinplot(sum(A),c,'ShowData',true)
    box off; set(gca,'color','none');
    ylabel('Gene number'); xlim([0.1 6.9]);%ylim([0 6000])
    xticklabels(Zhong_PFC_uniq_cn);
    xtickangle(45)
    title('Human developing prefrontal cortex(Zhong etal,Nature,2018)')
    
% Count gene numbsers
    A = Zhong_PFC; A(A>0)=1; 
    c=Zhong_PFC_clu;
    n=length(unique(c)); 
    clear gene_num gene_num2
    for i =1:n
        B = A(:,find(c==i));
        C = sum(B);
        gene_num(i,1) = size(B,2);
        gene_num(i,2) = round(mean(C));
        gene_num(i,3) = round(prctile(C,5));
        gene_num(i,4) = round(prctile(C,25));
        gene_num(i,5) = round(median(C));
        gene_num(i,6) = round(prctile(C,75));
        gene_num(i,7) = round(prctile(C,95));
        if size(B,2)<=100
        gene_num(i,8) = length(find(sum(B,2)>=5));
        gene_num(i,9) = 0;
        else
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:gene_num(i,1),100,'Replace',false)),2)>5));
        end
        gene_num(i,8) = round(mean(temp));
        gene_num(i,9) = round(std(temp));
        end
    end
        B = A;
        for j=1:100
            temp(j) = length(find(sum(B(:,datasample(1:size(B,2),1000,'Replace',false)),2)>5));
        end
        gene_num2(1) = round(mean(temp))    
    
% binarize the gene clustering as germexp(17941 genes):(somaexp, somaunexp)
A = Zhong_PFC; A(A>0)=1; 
gn = Zhong_PFC_gn;
input_gene = human_genename(find(h_gene_cluster>0));
input_gene_id = h_gname(find(h_gene_cluster>0));
input_gene_clu(1:length(input_gene_id)) = 0;
[~,~,i] = intersect(gn(find(sum(A')>=5)),input_gene,'stable');
input_gene_clu(i) = 1;
    figure;
    fn1='YOUR_PATH\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\human_ensembl90.genebody.singleFreq.tsv'
    F_compare_snv_asymmetry_v2(  fn1,   fn2,...
                                     input_gene_clu,...
                                     input_gene_id,...
                                     'F_compare_variants','only',...
                                     'FaceColor',exunex_cm,...
                                     'Title_prefix','Genebody'); 
    ylim([0 250]);
    title('Zhong,2018,Genebody')   
 
                                 