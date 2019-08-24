

%%GENERAL processing the data
if Q1 % inport the raw UMI matrices of human_TESE replicates
    % Human TESE1 first replicate (2017-11-03)
    fn='YOUR_PATH\rawdata\TESE1_1_gene_expression.tsv'
    [A1,cn1,g1] = F_inDrop_filter(fn,[0,1000,0.2,0.2]);   
    cn1=strcat('TESE1_1_',cn1);
    % Human TESE1 second replicate (2018-08-17)
    fn='YOUR_PATH\rawdata\TESE1_2_gene_expression.tsv'
    [A2,cn2,g2] = F_inDrop_filter(fn,[0,1000,0.2,0.2]);
    cn2=strcat('TESE1_2_',cn2);
    
    % Human TESE2 first replicate (2018-10-30)
    fn='YOUR_PATH\rawdata\TESE2_1_gene_expression.tsv'
    [A3,cn3,g3] = F_inDrop_filter(fn,[0,1000,0.2,0.2]);   
    % Human TESE2 second replicate (2018-10-30)
    fn='YOUR_PATH\rawdata\TESE2_2_gene_expression.tsv'
    [A4,cn4,g4] = F_inDrop_filter(fn,[0,1000,0.2,0.2]);
    
    %import gene signatures
    fn='YOUR_PATH\rawdata\human_genes.tsv'
    a          = tdfread(fn,'\t');
    human_genename = cellstr(a.name(strmatch('protein_coding',cellstr(a.type),'exact'),:));
    [~,~,prot_name_idx] = intersect(human_genename,g3,'stable');
    
    %set up output matrix, cell names and protein-coding gene names
    human_tese=[A1(prot_name_idx,:) A2(prot_name_idx,:) A3(prot_name_idx,:) A4(prot_name_idx,:)];
    human_cellname = [cn1 cn2 cn3 cn4];
    %specify replicates
    h1=1:length(A1(1,:));
    h2=(1+length(A1(1,:))) : (length(A1(1,:))+length(A2(1,:)));
    h3=(1+length(A1(1,:))+length(A2(1,:))) : (length(A1(1,:))+length(A2(1,:))+length(A3(1,:)));
    h4=(1+length(A1(1,:))+length(A2(1,:))+length(A3(1,:))) : ...
       (length(A1(1,:))+length(A2(1,:))+length(A3(1,:))+length(A4(1,:)));
    human_cellidx(h1)=1;
    human_cellidx(h2)=2;
    human_cellidx(h3)=3;
    human_cellidx(h4)=4;
    
    %knnsmooth
    k = 3;
    num_of_pc = 10;
    A1_sm = knn_smooth(A1(prot_name_idx,:),k,num_of_pc);
    A2_sm = knn_smooth(A2(prot_name_idx,:),k,num_of_pc);
    A3_sm = knn_smooth(A3(prot_name_idx,:),k,num_of_pc);
    A4_sm = knn_smooth(A4(prot_name_idx,:),k,num_of_pc);
    human_tese_sm = [A1_sm A2_sm A3_sm A4_sm];
    
    %save the results for further usage
    save('YOUR_PATH\human_TESE.mat',...
        'human_tese', 'human_tese_sm', 'human_genename','human_cellname','human_cellidx');


end

    %setup colormap
    cd 'YOUR_PATH\cbrewer'
    %cbrewer()
    type_cm = cbrewer('qual','Paired',20);
    type_cm_10 = cbrewer('qual','Paired',20);
    reps_cm = cbrewer('qual','Paired',12);
   % reps_cm = [reps_cm(1:2,:);reps_cm(5:6,:)];
    h_germtype_cm = cbrewer('qual','Paired',14);
    RdBu_cm = flipud(cbrewer('div','RdBu',50));
    Blues_cm = cbrewer('seq','Blues',50);
    genecluster_cm = [0.5 0.5 0.5;cbrewer('qual','Set2',5)];
    genecluster_facecm = [0.5 0.5 0.5;cbrewer('qual','Set2',5)];
    BuPu_cm = cbrewer('seq','BuPu',50);
    PuBuGn_cm = cbrewer('seq','PuBuGn',50);
        Reds_cm = cbrewer('seq','Reds',20);
    cd 'YOUR_PATH'
    


if Clustering_and_selecting_cells %---Relates to Fig.1,S1
    %load data from initial processing
    load('YOUR_PATH\human_TESE.mat');
	
	
	%Calculate sample stats
    A = human_tese;
    A(A>0)=1;
    for i = 4
        B=sum(A(:,find(human_cellidx==i)));
        B=sum(A);
        mean(B)
        median(B)
        prctile(B,5)
        prctile(B,25)
        prctile(B,75)
        prctile(B,95)
    gene_names = human_genename;
    B=A(:,find(human_cellidx==i));
    mito = [strmatch('MT-',gene_names);strmatch('mt-',gene_names)];
    mito_ratio=mean((sum(B(mito,:))./sum(B)))
    %ribo_ratio. Work for human and mouse
    ribo = [strmatch('RPL',gene_names);strmatch('RPS',gene_names);...
        strmatch('Rpl',gene_names);strmatch('Rps',gene_names)]; 
    ribo_ratio=mean((sum(B(ribo,:))./sum(B)))
    end
    %sample replicates
    corrcoef(mean(A(:,find(human_cellidx==1))'),mean(A(:,find(human_cellidx==2))'))
    corrcoef(mean(A(:,find(human_cellidx==3))'),mean(A(:,find(human_cellidx==4))'))
    corrcoef(mean(A(:,find(human_cellidx<=2))'),mean(A(:,find(human_cellidx>2))'))
	
	
	
	
	
    gname = human_genename;
    A = log10(1+F_normalize(human_tese_sm,'median'));
    A_zscore = zscore(A,0,2);
                                   
    %kmeans clustering with Elbow-method optimal k
    rng(1) % for reproducibility
    [h_kmeans_IDX,h_PCA_score,~,~,h_tSNE_score] = F_kmeans(human_tese,...
                                       'knnsmooth','yes',...
                                       'k_for_knnsmooth',3,...
                                       'FTtransPCA','yes',...
                                       'SelectVariableGene','yes',...
                                       'number_of_pc',50,...
                                       'plot_PCA','yes',...
                                       'plot_tSNE','no',...
                                       'visual_mode','on'...
                                       );  


    %kmeans clustering with Elbow-method optimal k
    rng(1) % for reproducibility
    [h_kmeans_IDX,h_PCA_score,~,~,h_tSNE_score] = F_kmeans(human_tese,...
                                       'SelectVariableGene','yes',...
                                       'number_of_pc',50,...
                                       'plot_PCA','yes',...
                                       'plot_tSNE','no',...
                                       'visual_mode','on'...
                                       );

    %plot replicates
    figure;gscatter(h_tSNE_score(:,1),h_tSNE_score(:,2),human_cellidx',reps_cm,'.',10);
    figure;gscatter(h_PCA_score(:,1),h_PCA_score(:,2),human_cellidx',reps_cm,'.',6);
           set(gca,'color','none');xticks([]);yticks([]);
           title('human donor reps');legend({'D1-T1','D1-T2','D2-T1','D2-T2'})
    figure;gscatter(h_PCA_score(:,1),h_PCA_score(:,2),human_cellidx','rgmb','.',5);
    figure;scatter3(h_PCA_score(:,1),h_PCA_score(:,2),h_PCA_score(:,3),10,human_cellidx,'filled');

%     gif('test.gif','DelayTime',0.02,'frame',gcf)
%     for k = 1:72
%         view([5*k 10])
%         gif
%     end
    
    %Color with genes of interest
%      marker=cellstr({'FGFR3','DMRT1','CTCFL','SYCP3','SPO11','PIWIL1','ACRV1','TNP1','PRM1'});
%     marker=cellstr({'ARX','CYP11A1','IGF1','WT1','SOX9','FGFR3','DMRT1','STRA8','GFRA1'});
%     marker=cellstr({'ID2','IL7R','CD52','THY1','RORA','APOE','DAB2','CD74','ADGRE1'});
%     marker=cellstr({'VWF','TIE1','TEK','APOE','CD74','ID2','JUND','STRA8','CTCFL'});
%     marker=cellstr({'WT1','SOX9','ARX','IGF1','MYH11','ACTA2','CD68','CD163','FGFR3'});
     marker=cellstr({'WT1','SOX9','DDX4','CSF1','IGF1','FGFR3','MYH11','ACTA2','SYCP3','CD68','CD163','ACRV1','PECAM1','VWF','PRM1'});
    figure;
    for i=1:15
        goi = strmatch(marker(i),gname,'exact');
        subplot(5,3,i);
         scatter(h_PCA_score(:,1),h_PCA_score(:,2),1,A_zscore(goi,:),'filled');
        %xlabel('tSNE1'); ylabel('tSNE2');
        colormap(PuBuGn_cm(10:end,:));title(gname(goi));
        set(gca,'color','none');xticklabels({});yticklabels({})
    end
    

    %% Remove somatic cells based on marker genes
    human_germcell_colormap = [    
    0.3281    0.1523    0.5313 %Sg-1
    0.5000    0.4492    0.6719 %Sg-2
    0.6094    0.2188    0.3789 %Sc-1
    0.7734    0.1758    0.2930 %Sc-2
    0.9297    0.3906    0.2695 %Sc-3
    0.9844    0.6563    0.3672 %Sc-4
    0.9922    0.8750    0.5430 %RS-1
    0.9922    0.9375    0.6523 %RS-2
    0.9336    0.9688    0.6484 %RS-3
    0.7305    0.8867    0.6250 %RS-4
    0.4453    0.7734    0.6406 %ES-1
    0.3242    0.6719    0.6758 %ES-2
    0.2109    0.5508    0.7266 %ES-3
    0.0137    0.3828    0.6738 %ES-4
];

    %by visulizating the scatter plots (PCA&tSNE) above, remove somatic cells and unknown cells.
    
    h_germ  = [find(h_kmeans_IDX==16);find(h_kmeans_IDX==24);find(h_kmeans_IDX==13);find(h_kmeans_IDX==14);...
             find(h_kmeans_IDX==20);find(h_kmeans_IDX==5);find(h_kmeans_IDX==15);find(h_kmeans_IDX==21);...
             find(h_kmeans_IDX==9);find(h_kmeans_IDX==6);find(h_kmeans_IDX==12);find(h_kmeans_IDX==1);...
             find(h_kmeans_IDX==19);find(h_kmeans_IDX==17)];
    h_soma  = [find(h_kmeans_IDX==2);find(h_kmeans_IDX==4);find(h_kmeans_IDX==7);...
               find(h_kmeans_IDX==8);find(h_kmeans_IDX==10);find(h_kmeans_IDX==11);find(h_kmeans_IDX==18)];
    %redefine germ-somatic cell clusters
    h_germ_som(1:length(h_kmeans_IDX)) = 0;
    h_germ_som(find(h_kmeans_IDX==16)) = 1;
    h_germ_som(find(h_kmeans_IDX==24)) = 2; 
    h_germ_som(find(h_kmeans_IDX==13))  = 3;
    h_germ_som(find(h_kmeans_IDX==14)) = 4;
    h_germ_som(find(h_kmeans_IDX==20)) = 5;
    h_germ_som(find(h_kmeans_IDX==5))  = 6;
    h_germ_som(find(h_kmeans_IDX==15)) = 7;
    h_germ_som(find(h_kmeans_IDX==21))  = 8;
    h_germ_som(find(h_kmeans_IDX==9)) = 9;
    h_germ_som(find(h_kmeans_IDX==6)) = 10;
    h_germ_som(find(h_kmeans_IDX==12)) = 11;
    h_germ_som(find(h_kmeans_IDX==1))  = 12;
    h_germ_som(find(h_kmeans_IDX==19))  = 13;
    h_germ_som(find(h_kmeans_IDX==17))  = 14;
    h_germ_som(h_soma)                  = 15;  %merge somatic cells as a single cluster
    h_germ_IDX=h_germ_som(h_germ)';
    
    cm=[0.8 0.8 0.8;human_germcell_colormap;0.55 0.55 0.55];
    figure;
    for i=1:16
    scatter(h_PCA_score(find(h_germ_som==(i-1)),1),h_PCA_score(find(h_germ_som==(i-1)),2),...
            10,'filled','MarkerFaceColor',cm(i,:),'MarkerFaceAlpha',0.8);hold on;
    end
    xlabel('PC1'); ylabel('PC2');legend on;
    set(gca,'color','none');    box off;
    title('human testicular cells');
    
    %Play with germ cells
    h_germcell = human_tese(:,h_germ);
    h_germcellname = human_cellname(h_germ);
    h_germcell_repidx  = human_cellidx(h_germ);
    %knn_smooth with k=3, and number_of_PC as 20
    h_germcell_sm = knn_smooth(h_germcell,3,20);
    
    %generate tSNE embedding of germcells
    rng(1) % for reproducibility
    [~,h_germ_PCAscore,~,~,h_germ_tSNE_score] = F_kmeans(h_germcell_sm,...
        'SelectVariableGene','yes','number_of_pc',50,'plot_PCA','yes','plot_tSNE','no','visual_mode','on');
    figure;gscatter(h_germ_PCAscore(:,1),h_germ_PCAscore(:,2),h_kmeans_IDX(h_germ),h_germtype_cm,'.',5);
    figure;gscatter(h_germ_PCAscore(:,1),h_germ_PCAscore(:,2),h_germ_IDX,human_germcell_colormap,'.',5);
    %Plot germ cell markers
    B = log10(1+F_normalize(h_germcell,'median'));
    B_zscore = zscore(B,0,2);
    marker=cellstr({'FGFR3','DMRT1','CTCFL','SYCP3','SPO11','PIWIL1','ACRV1','TNP1','PRM1'});
    figure;
    for i=1:9
        goi = strmatch(marker(i),gname,'exact');
        subplot(3,3,i);scatter(h_germ_tSNE_score(:,1),h_germ_tSNE_score(:,2),3,B_zscore(goi,:),'filled');
        %xlabel('tSNE1'); ylabel('tSNE2');
        colormap(PuBuGn_cm);title(gname(goi));
    end
    
    
    %soma from TESE1 & TESE2
    rng(1) % for reproducibility, use k=5 based prior knowledge and early exploratory of the data.
    [h_soma_IDX,~,~,~,h_soma_tsne] = F_kmeans(human_tese(:,h_soma),'k_for_kmeans',5,...
                                       'SelectVariableGene','yes',...
                                       'number_of_pc',50,...
                                       'plot_PCA','no',...
                                       'plot_tSNE','yes',...
                                       'visual_mode','on'...
                                       );
    figure;gscatter(h_soma_tsne(:,1),h_soma_tsne(:,2),human_cellidx(h_soma),type_cm_10,'.',10);
    marker=cellstr({'WT1','SOX9','DDX4','CSF1','IGF1','FGFR3','MYH11','ACTA2','SYCP3','CD68','CD163','ACRV1','PECAM1','VWF','PRM1'});
    marker=cellstr({'WT1','SOX9','DDX4','CSF1','IGF1','FGFR3','MYH11','ACTA2','SYCP3','CD68','CXCL16','CXCR6','PECAM1','VWF','PRM1'});
    figure;
    for i=1:15
        goi = strmatch(marker(i),gname,'exact');
        %goi = gname(strmatch('LAMP2',gname))
        subplot(5,3,i);
         scatter(h_soma_tsne(:,1),h_soma_tsne(:,2),2,A_zscore(goi,h_soma),'filled');
        %xlabel('tSNE1'); ylabel('tSNE2');
        colormap(PuBuGn_cm(10:end,:));title(gname(goi));
        ylim([-32 30]); xlim([-11 35]); 
        set(gca,'color','none');xticklabels({});yticklabels({})
    end
    %by visulizating the scatter plots (PCA&tSNE) above, remove somatic cells and unknown cells.
    human_somatic_cells = human_tese(:,h_soma);
    soma_sertoli = human_somatic_cells(:,find(h_soma_IDX==4));
    soma_leydig = human_somatic_cells(:,find(h_soma_IDX==5));
    soma_macrophage = human_somatic_cells(:,find(h_soma_IDX==2));
    soma_myoid = human_somatic_cells(:,find(h_soma_IDX==1));
    soma_endothelia = human_somatic_cells(:,find(h_soma_IDX==3));

    
    %save the results for further analysis
    save('YOUR_PATH\human_tese1and2_replicates_germcells.mat',...
        'h_germ','h_germcell', 'h_germcell_sm','h_germcellname','human_genename','h_germ_IDX',...
        'h_germcell_repidx','h_germtype_cm');
    save('YOUR_PATH\human_tese1and2_replicates_soma.mat',...
        'h_soma','h_soma_tsne','human_somatic_cells','h_soma_IDX','human_genename',...
        'soma_sertoli','soma_leydig', 'soma_macrophage','soma_myoid','soma_endothelia');


    %% save germ cell info into .csv files, for Monocle2 pipeline
    csvwrite('YOUR_PATH\files_for_monocle\TESE1AND2\human_germcell_exp_unsm.csv', h_germcell);
    csvwrite('YOUR_PATH\files_for_monocle\TESE1AND2\human_germcell_exp_k3sm.csv', h_germcell_sm);
    cell2csv('YOUR_PATH\files_for_monocle\TESE1AND2\human_germcell_genename.csv', ['GENENAME';human_genename]);
    cell2csv('YOUR_PATH\files_for_monocle\TESE1AND2\human_germcell_cellname.csv', ['CELLNAME';human_cellname(h_germ)']);
    %Monocle pseudotime colormap
    cell2csv('YOUR_PATH\files_for_monocle\TESE1AND2\pseudotime_colormap.csv', cellstr(rgb2hex(PuBuGn_cm)));
    

    
    %Monocle2 results--Use k=3 smooth 
    fn='YOUR_PATH\Monocle2_results\human_tese1and2\humanMerge1234_smoothK-3_pseudotime-Order.csv';
    human_k3smooth_pseudotimeOrder = readtable(fn);
    
%     human_order = human_unsmooth_pseudotimeOrder.ORDER;
    h_germcell_pseudotimeorder = human_k3smooth_pseudotimeOrder.ORDER;
    h_germcell_pseudotime = human_k3smooth_pseudotimeOrder.Pseudotime;
    
    for i=1:14
        mean_time(i)=mean(h_germcell_pseudotime(find(h_germ_IDX==i)));
    end
    
    load('YOUR_PATH\human_tese1and2_replicates_germcells.mat');
    h_germcell_sm_trans = F_normalize(F_normalize(h_germcell_sm,'median'),'FT');
    
    figure;
    scatter(h_germ_PCAscore(:,1),h_germ_PCAscore(:,2),10,h_germcell_pseudotime,'filled');colormap(PuBuGn_cm(11:50,:))
        
    %SAVE the MOnocle2 pseudotime order with the germcells files
    save('YOUR_PATH\human_tese1and2_replicates_germcells.mat',...
         'h_germcell_pseudotime', 'h_germcell_pseudotimeorder','-append');
 
    %% save germ cell info into .csv files, for RNA_velocity pipeline
    %2D_space_coordination
    CELLNAME=h_germcellname';
    V1=h_germ_PCAscore(:,1);
    V2=h_germ_PCAscore(:,2);
    CLUSTER=h_germ_IDX;
    twoD_space_coordination = table(CELLNAME,V1,V2,CLUSTER);
    writetable(twoD_space_coordination,'YOUR_PATH\files_for_velocity\human\TESE1AND2\2D_space_coordination.csv');
    %gene_name
    cell2csv('YOUR_PATH\files_for_velocity\human\TESE1AND2\gene_name.csv', ['GENENAME';human_genename]);
    %cell_name
    human_rep_1_cellname = erase(CELLNAME(strmatch('TESE1_1',CELLNAME)),'TESE1_1_');
    cell2csv('YOUR_PATH\files_for_velocity\human\TESE1AND2\human_rep_1_cellname.csv', ['CELLNAME';human_rep_1_cellname]);
    human_rep_2_cellname = erase(CELLNAME(strmatch('TESE1_2',CELLNAME)),'TESE1_2_');
    cell2csv('YOUR_PATH\files_for_velocity\human\TESE1AND2\human_rep_2_cellname.csv', ['CELLNAME';human_rep_2_cellname]);
    human_rep_3_cellname = erase(CELLNAME(strmatch('TESE2_1',CELLNAME)),'TESE2_1_');
    cell2csv('YOUR_PATH\files_for_velocity\human\TESE1AND2\human_rep_3_cellname.csv', ['CELLNAME';human_rep_3_cellname]);
    human_rep_4_cellname = erase(CELLNAME(strmatch('TESE2_2',CELLNAME)),'TESE2_2_');
    cell2csv('YOUR_PATH\files_for_velocity\human\TESE1AND2\human_rep_4_cellname.csv', ['CELLNAME';human_rep_4_cellname]);
    %UMI_matrix
    csvwrite('YOUR_PATH\files_for_velocity\human\TESE1AND2\human_UMI_matrix_unsm.csv', h_germcell);
    %germcell_clusters
    %human_germcell cluster: pseudotime colormap
    cell2csv('YOUR_PATH\files_for_velocity\human\TESE1AND2\human_germcell_colormap.csv', cellstr(rgb2hex(human_germcell_colormap)));
    
end



if clusering_genes  %---Relates to Fig.1,S1
    %% Grouping germ cell expressed genes: use both TESE donors samples as input
    %cluster genes by clusters of cells
    load('YOUR_PATH\human_tese1and2_replicates_germcells.mat');
    
    A = h_germcell_sm;
    %cluster genes with k=5 which is optimal
    rng(1)
    [h_cluster_expression,h_cluster_exp_meanUMI,h_cluster_index]=F_cluster_gene(A,h_germ_IDX,5,...
                                                                    'expression_cutoff_n',5,...
                                                                    'umi_ratio_cutoff',0,...
                                                                    'include_unexpressed_gene','yes',...
                                                                    'ylimit',[-2 3]);
    % reorder the gene cluster by visulization of expression stage
    h_gene_cluster = h_cluster_index;
    h_gene_cluster(find(h_cluster_index==1))=1;
    h_gene_cluster(find(h_cluster_index==2))=4;
    h_gene_cluster(find(h_cluster_index==3))=5;
    h_gene_cluster(find(h_cluster_index==4))=2;
    h_gene_cluster(find(h_cluster_index==5))=3;
     
    % define gene clusters by replicates
    rep1 = find(h_germcell_repidx<3);
    rep2 = find(h_germcell_repidx>2); 
    donor1_cellcluster = h_germcell_repidx(rep1);
    donor2_cellcluster = h_germcell_repidx(rep2);
    histcounts(donor1_cellcluster)
    histcounts(donor2_cellcluster)
    rng(1)
    [~,~,h_cluster_index1]=F_cluster_gene(   A(:,rep1),h_germ_IDX(rep1),5,...
                                             'expression_cutoff_n',5,...
                                             'umi_ratio_cutoff',0,...
                                             'include_unexpressed_gene','yes',...
                                             'ylimit',[-2 3]);
    rep1_cluster = h_cluster_index1;
    rep1_cluster(find(h_cluster_index1==1))=2;
    rep1_cluster(find(h_cluster_index1==2))=5;
    rep1_cluster(find(h_cluster_index1==3))=1;
    rep1_cluster(find(h_cluster_index1==4))=4;
    rep1_cluster(find(h_cluster_index1==5))=3;
    
    rng(1)
    [~,~,h_cluster_index2]=F_cluster_gene(   A(:,rep2),h_germ_IDX(rep2),5,...
                                             'expression_cutoff_n',5,...
                                             'umi_ratio_cutoff',0,...
                                             'include_unexpressed_gene','yes',...
                                             'ylimit',[-2 3]);
    rep2_cluster = h_cluster_index2;
    rep2_cluster(find(h_cluster_index2==1))=2;
    rep2_cluster(find(h_cluster_index2==2))=1;
    rep2_cluster(find(h_cluster_index2==3))=5;
    rep2_cluster(find(h_cluster_index2==4))=4;
    rep2_cluster(find(h_cluster_index2==5))=3;

    % plot human_TESE1 and human_TESE2 results
    for i=1:6
        for j=1:6
            genegroup(j,i)=length(find(rep1_cluster==(i-1) & rep2_cluster==(j-1)));
        end
    end
    figure;imagesc(genegroup);colormap(Reds_cm)
    yticklabels({'C0','C1','C2','C3','C4','C5'});ylabel('human TESE2 gene clusters');
    xticklabels({'C0','C1','C2','C3','C4','C5'});xlabel('human TESE1 gene clusters');
    title('human TESE1 & TESE2 gene cluster overlapping')
    % reorder the cluster and plot
    % angle = 20 works.
    % h_zavit_cluster = F_cluster_sorting_by_ZAVIT(h_cluster_expression,h_cluster_index,120,'ylimit',[-2 3]);

    %% Visulization2

    figure;
    explode = [1 0 0 0 0 0];
    pie(histcounts(h_gene_cluster),explode);colormap(genecluster_facecm)
        legend({'C0','C1','C2','C3','C4','C5'});

    % heatmap
    exp = zscore(h_cluster_expression')';
    xticklength=length(exp(1,:));
    [i,xi]=sort(h_gene_cluster);
    figure;
    %draw heatmap as the genes sorted;
    subplot('position',[0.25 0.1 0.65 0.85]);
    imagesc(exp(xi,:),[-2.5 2.5]);colorbar;colormap(RdBu_cm)
    h_gca = gca;
    xlim([0.5, xticklength+0.5]);
    h_gca.XTick = 1:xticklength;
    xlabel('Cell clusters');
    title('Human gene expression by clusters')
    %draw heatmap of gene clusters
    f=subplot('position',[0.1 0.1 0.05 0.85]);
    imagesc(i);xticks([]);yticks([]);colormap(f,genecluster_cm);
    ylabel('Gene clusters');
    
    % gene cluster plots
    j = length(unique(h_gene_cluster));
    n = unique(h_gene_cluster);
    figure;
    for i = 1:j
        subplot(j,1,i);
        y = exp(find(h_gene_cluster==n(i)),:);
        x =  1:xticklength;
        shadedErrorBar(x, mean(y,1),std(y), ...
            'lineprops',{'c-o','MarkerSize',3,...
            'Color',genecluster_facecm(i,:),...
            'MarkerFaceColor',[178,24,43]/255});
        ylim([-2 3.3]);
        xlim([0.5, xticklength+0.5]);
        h_gca = gca;
        h_gca.XTick = 1:xticklength;
        title(sprintf('Cluster %d (n=%d)',n(i),length(find(h_gene_cluster==n(i)))));
    end
    
    %Save the gene expression of stages and 
    save('YOUR_PATH\h_m_stage_expression.mat',...
        'h_cluster_expression', 'h_cluster_exp_meanUMI','h_gene_cluster', 'human_genename','-append');
    
    
    
end



