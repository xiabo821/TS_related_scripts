%%mouse data analaysis for revision

%%GENERAL processing the data
if Q1 % inport the UMI matrices of mouse testes replicates
    % First replicate 
    fn='YOUR_PATH\rawdata\Mouse_4m_1_gene_expression.tsv';
    [A1,cn1,g1] = F_inDrop_filter(fn,[0,1000,0.2,0.2]);
    cn1=strcat('Mouse_4m_1_',cn1);
    % Second replicate 
    fn='YOUR_PATH\rawdata\Mouse_4m_2_gene_expression.tsv';
    [A2,cn2,g2] = F_inDrop_filter(fn,[0,1000,0.2,0.2]);
    cn2=strcat('Mouse_4m_2_',cn2);

    %import gene signatures
    fn='YOUR_PATH\rawdata\mouse_genes.tsv';
    a          = tdfread(fn,'\t');
    mouse_gene_info_Ensembl91 = [cellstr(a.name) cellstr(a.ensembl_id) ...
                                 cellstr(a.chromosome) cellstr(a.source) cellstr(a.type)];
    mouse_genename = cellstr(a.name(strmatch('protein_coding',cellstr(a.type),'exact'),:));
    [~,~,prot_name_idx] = intersect(mouse_genename,g1,'stable');
    
    %set up output matrix, cell names and protein-coding gene names
    mouse_4m = [A1(prot_name_idx,:) A2(prot_name_idx,:)];
    mouse_cellname = [cn1 cn2];
    %specify replicates
    h1=1:length(A1(1,:));
    h2=(1+length(A1(1,:))) : (length(A1(1,:))+length(A2(1,:)));
    mouse_cellidx(h1)=1;
    mouse_cellidx(h2)=2;
    
    %knnsmooth
    k = 3;
    num_of_pc = 10;
    A1_sm = knn_smooth(A1(prot_name_idx,:),k,num_of_pc);
    A2_sm = knn_smooth(A2(prot_name_idx,:),k,num_of_pc);
    mouse_4m_sm = [A1_sm A2_sm];
    
    %save the results for further analysis
    save('YOUR_PATH\mouse_4m_replicates.mat',...
        'mouse_4m', 'mouse_4m_sm', 'mouse_genename','mouse_cellname','mouse_cellidx',...
        'mouse_gene_info_Ensembl91');

end

    %setup colormap
    cd 'YOUR_PATH\cbrewer'
    %cbrewer()
    type_cm = cbrewer('qual','Paired',19);
    m_germtype_cm = cbrewer('qual','Paired',13);
    h_germtype_cm = cbrewer('qual','Paired',14);
    h_germtype_cm = cbrewer('qual','Accent',14);
    RdBu_cm = flipud(cbrewer('div','RdBu',50));
    genecluster_cm = [0.32 0.32 0.32;cbrewer('qual','Set2',5)];
    Blues_cm = cbrewer('seq','Blues',50);
    BuPu_cm = cbrewer('seq','BuPu',50);
    PuBuGn_cm = cbrewer('seq','PuBuGn',50);
    cd 'YOUR_PATH'
    
        figure;
    scatter(m_germ_tSNE_score(:,1),m_germ_tSNE_score(:,2),10,m_germ_IDX,'filled');
    colormap(m_germtype_cm)


if mouse_4m  %---Relates to Fig.S2A,C,D
    %load data from initial processing
    load('YOUR_PATH\mouse_4m_replicates.mat');
    gname = mouse_genename;
    A = log10(1+F_normalize(mouse_4m,'median'));
    A_zscore = zscore(A,0,2);


    %kmeans clustering with Elbow-method optimal k
    rng(1) % for reproducibility
    [m_kmeans_IDX,m_PCA_score,~,~,m_tSNE_score] = F_kmeans(mouse_4m,...
                                       'FTtransPCA','no',...
                                       'number_of_pc',50,...
                                       'plot_PCA','yes',...
                                       'plot_tSNE','yes'...
                                       ); 
     %color replicates                          
    figure;gscatter(m_tSNE_score(:,1),m_tSNE_score(:,2),mouse_cellidx','rb','.',10);
    figure;gscatter(m_PCA_score(:,1),m_PCA_score(:,2),mouse_cellidx','rb','.',10);

    %Color with genes of interest
   marker=cellstr({'Arx','Igf1','Dmrt1','Stra8','Sycp3','Piwil1','Acrv1','Tnp1','Prm1'});
    marker=cellstr({'Arx','Cyp11a1','Igf1','Wt1','Sox9','Myh11','Acta2','Fgfr3','Dmrt1'});
    figure;
    for i=1:9
        goi = strmatch(marker(i),gname,'exact');
        subplot(3,3,i);scatter(m_tSNE_score(:,1),m_tSNE_score(:,2),5,A_zscore(goi,:),'filled');
       % xlabel('tSNE1'); ylabel('tSNE2');
        colormap(RdBu_cm);title(gname(goi));
    end   
    
    
    %Remove somatic cells based on marker genes
    %by visulizating the scatter plots (PCA&tSNE) above, remove somatic cells.
    m_remove= find(m_kmeans_IDX==8 | m_kmeans_IDX==9 | m_kmeans_IDX==15);
    m_germ  = [find(m_kmeans_IDX==12);find(m_kmeans_IDX==16);find(m_kmeans_IDX==14);find(m_kmeans_IDX==10);...
             find(m_kmeans_IDX==5);find(m_kmeans_IDX==1);find(m_kmeans_IDX==11);find(m_kmeans_IDX==6);...
             find(m_kmeans_IDX==13);find(m_kmeans_IDX==4);find(m_kmeans_IDX==2);find(m_kmeans_IDX==7);find(m_kmeans_IDX==3)];
    %redefine germ-somatic cell clusters
    m_germ_som(1:length(m_kmeans_IDX))  = 0;
    m_germ_som(find(m_kmeans_IDX==12))  = 1;
    m_germ_som(find(m_kmeans_IDX==16))  = 2;
    m_germ_som(find(m_kmeans_IDX==14))  = 3;
    m_germ_som(find(m_kmeans_IDX==10))  = 4;
    m_germ_som(find(m_kmeans_IDX==5))   = 5;
    m_germ_som(find(m_kmeans_IDX==1))   = 6;
    m_germ_som(find(m_kmeans_IDX==11))  = 7;
    m_germ_som(find(m_kmeans_IDX==6))   = 8;
    m_germ_som(find(m_kmeans_IDX==13))  = 9;
    m_germ_som(find(m_kmeans_IDX==4))   = 10;
    m_germ_som(find(m_kmeans_IDX==2))   = 11;
    m_germ_som(find(m_kmeans_IDX==7))   = 12;
    m_germ_som(find(m_kmeans_IDX==3))   = 13;
    m_germ_IDX=m_germ_som(m_germ)';

    figure;gscatter(m_tSNE_score(:,1),m_tSNE_score(:,2),m_germ_som,type_cm,'.',10);
    xlabel('tSNE1'); ylabel('tSNE2');
    legend({'SomaticCell','GermCell'});
    title('Seperate Somatic vs Germ cells');

        
    %Play with germ cells
    m_germcell = mouse_4m(:,m_germ);
    m_germcellname = mouse_cellname(m_germ);
    m_germcell_repidx  = mouse_cellidx(m_germ);
    %knn_smooth with k=3, and number_of_PC as 20
    m_germcell_sm = knn_smooth(m_germcell,3,20);


    mouse_germcell_colormap = [    
    0.3281    0.1523    0.5313
    0.5000    0.4492    0.6719 
    0.6094    0.2188    0.3789
    0.7734    0.1758    0.2930
    0.9297    0.3906    0.2695 
    0.9844    0.6563    0.3672 
    0.9922    0.8750    0.5430 
    0.9922    0.9375    0.6523
    0.9336    0.9688    0.6484
    0.7305    0.8867    0.6250 
    0.4453    0.7734    0.6406
    0.3242    0.6719    0.6758
    0.1123    0.4668    0.7002];

    rng(1) % for reproducibility
    [~,m_germ_PCAscore,~,~,~] = F_kmeans(m_germcell_sm,...
        'SelectVariableGene','yes','number_of_pc',50,...
        'plot_PCA','yes','plot_tSNE','no','visual_mode','on');
    %plot replicates
   marker=cellstr({'Arx','Igf1','Dmrt1','Stra8','Sycp3','Piwil1','Acrv1','Tnp1','Prm1'});
    figure;
    for i=1:9
        goi = strmatch(marker(i),gname,'exact');
        subplot(3,3,i);scatter(m_germ_PCAscore(:,1),m_germ_PCAscore(:,2),5,A_zscore(goi,m_germ),'filled');
        colormap(RdBu_cm);title(gname(goi));
    end
    figure;gscatter(m_germ_tSNE_score(:,1),m_germ_tSNE_score(:,2),m_germ_IDX,mouse_germcell_colormap,'.',6);
    figure;gscatter(m_germ_PCAscore(:,1),m_germ_PCAscore(:,2),m_germ_IDX,mouse_germcell_colormap,'.',6);
    
    %% save germ cell info into .csv files, for Monocle2 pipeline
    csvwrite('YOUR_PATH\files_for_monocle\mouse_reps\mouse_germcell_exp_unsm.csv', m_germcell);
    csvwrite('YOUR_PATH\files_for_monocle\mouse_reps\mouse_germcell_exp_k3sm.csv', m_germcell_sm);
    cell2csv('YOUR_PATH\files_for_monocle\mouse_reps\mouse_germcell_genename.csv', ['GENENAME';mouse_genename]);
    cell2csv('YOUR_PATH\files_for_monocle\mouse_reps\mouse_germcell_cellname.csv', ['CELLNAME';m_germcellname']);
    
    %save the results for further analysis
    save('YOUR_PATH\mouse_4m_replicates_germcells.mat',...
        'm_germ', 'm_germcell', 'm_germcell_sm', 'mouse_genename','m_germ_IDX','m_germcellname',...
        'm_germ_PCAscore','mouse_germcell_colormap');

    
    %Monocle2 results
%     fn='C:\Users\Bo\Box Sync\Pap33_Spermatogenesis\Analysis\Science_1st_review\Monocle2_results\mouse_4month_rep1_&_rep2\mouse_smoothK-0_pseudotime-Order.csv';
%     mouse_unsmooth_pseudotimeOrder = readtable(fn);
    fn='YOUR_PATH\Monocle2_results\mouse_4month_rep1_&_rep2\mouse_smoothK-3_pseudotime-Order.csv';
    mouse_k3smooth_pseudotimeOrder = readtable(fn);
    
%     mouse_order = mouse_unsmooth_pseudotimeOrder.ORDER;
    mouse_order = mouse_k3smooth_pseudotimeOrder.ORDER;
    m_pseudotime = mouse_k3smooth_pseudotimeOrder.Pseudotime;

    load('YOUR_PATH\mouse_4m_replicates_germcells.mat');
    m_germcell_sm_trans = F_normalize(F_normalize(m_germcell_sm,'median'),'Anscombe');
%     germcell_sm_trans = F_normalize(F_normalize(germcell_sm,'median'),'FT');
%     germcell_sm_trans = F_normalize(germcell_sm,'FT');
    %check cell ordering
    [~,j] = sort(mouse_order);
    figure;scatter(m_germ_tSNE_score(:,1),m_germ_tSNE_score(:,2),10,j,'filled'); colormap(jet)
    figure;imagesc(corrcoef(m_germcell_sm_trans(:,j)));colormap(RdBu_cm)
    figure;imagesc(corrcoef(m_germcell_sm_trans));colormap(RdBu_cm)
    
    figure;
    subplot('position',[0.07 0.15 0.05 0.8]);
        imagesc(m_germ_IDX(j));colormap(m_germtype_cm);xticks({})
    s=subplot('position',[0.15 0.15 0.8 0.8]);
        imagesc(corrcoef(m_germcell_sm_trans(:,j)));colormap(s,RdBu_cm);xticks({});yticks({}) 
    title('Corrcoef with k3smooth,FT-transformed matrix')
    s=subplot('position',[0.15 0.07 0.8 0.05]);
        imagesc(m_pseudotime(j)');colormap(s,PuBuGn_cm);yticks([]);
        
    %SAVE the ordered germ cells from Monocle2-ordered cells
    m_germcell_ordered = m_germcell(:,j);
    m_germcell_sm_ordered = m_germcell_sm(:,j);
    m_germcellname_ordered = m_germcellname(j);
    m_germ_tSNE_score_ordered  = m_germ_tSNE_score(j,:);
    m_germ_IDX_ordered = m_germ_IDX(j);
    save('YOUR_PATH\mouse_4m_replicates_germcells_ordered.mat',...
        'm_germcell_ordered', 'm_germcell_sm_ordered', 'mouse_genename','m_pseudotime',...
        'm_germcellname_ordered','m_germ_IDX_ordered','m_germ_tSNE_score_ordered','m_germtype_cm');
    %test if the Monocle2-ordered cells are good
    figure;
    scatter(m_germ_tSNE_score(:,1),m_germ_tSNE_score(:,2),10,m_germ_IDX,'filled');
    colormap(m_germtype_cm)
    figure;
    scatter(m_germ_tSNE_score_ordered(:,1),m_germ_tSNE_score_ordered(:,2),10,m_germ_IDX(j),'filled');
    colormap(m_germtype_cm)
    figure;
    gscatter(m_germ_tSNE_score_ordered(:,1),m_germ_tSNE_score_ordered(:,2),m_germ_IDX_ordered,m_germtype_cm,'.',10);
            
        
    
end 


if prepare_velocity_files

    %2D_space_coordination
    CELLNAME=m_germcellname';
    V1=m_germ_PCAscore(:,1);
    V2=m_germ_PCAscore(:,2);
    CLUSTER=m_germ_IDX;
    twoD_space_coordination = table(CELLNAME,V1, V2,CLUSTER);
    writetable(twoD_space_coordination,'YOUR_PATH\files_for_velocity\mouse\2D_space_coordination.csv');
    %gene_name
    cell2csv('YOUR_PATH\files_for_velocity\mouse\gene_name.csv', ['GENENAME';mouse_genename]);
    %cell_name
    mouse_rep_1_cellname = erase(CELLNAME(strmatch('Mouse_4m_1_',CELLNAME)),'Mouse_4m_1_');
    cell2csv('YOUR_PATH\files_for_velocity\mouse\mouse_rep_1_cellname.csv', ['CELLNAME';mouse_rep_1_cellname]);
    mouse_rep_2_cellname = erase(CELLNAME(strmatch('Mouse_4m_2_',CELLNAME)),'Mouse_4m_2_');
    cell2csv('YOUR_PATH\files_for_velocity\mouse\mouse_rep_2_cellname.csv', ['CELLNAME';mouse_rep_2_cellname]);
    %UMI_matrix
    csvwrite('YOUR_PATH\files_for_velocity\mouse\mouse_UMI_matrix_unsm.csv', m_germcell);
    %germcell_clusters
    %mouse_germcell cluster: pseudotime colormap
    cell2csv('YOUR_PATH\files_for_velocity\mouse\mouse_germcell_colormap.csv', cellstr(rgb2hex(mouse_germcell_colormap)));
    
end



if clusering_genes  %---Relates to Fig.S2C,D
    %% Grouping germ cell expressed genes
    %cluster genes by clusters of cells
    load('YOUR_PATH\mouse_4m_replicates_germcells_ordered.mat');
    
    A = m_germcell_sm_ordered;
    A = m_germcell_sm;
    
    %cluster genes with k=5
    rng(1)
    [m_cluster_expression,m_cluster_exp_meanUMI,m_cluster_index]=F_cluster_gene(A,m_germ_IDX,5,...
        'expression_cutoff_n',3,...
        'include_unexpressed_gene','yes',...
        'ylimit',[-2 3]);

    
        %% Visulization
    % reorder the gene cluster by visulization of expression stage
    m_gene_cluster = m_cluster_index;
    m_gene_cluster(find(m_cluster_index==5)) = 1;
    m_gene_cluster(find(m_cluster_index==4)) = 2;
    m_gene_cluster(find(m_cluster_index==3)) = 3;
    m_gene_cluster(find(m_cluster_index==2)) = 4;
    m_gene_cluster(find(m_cluster_index==1)) = 5;
    
    figure;
    explode = [1 0 0 0 0 0];
    pie(histcounts(m_gene_cluster),explode);colormap(genecluster_cm)
        legend({'C0','C1','C2','C3','C4','C5'});
    
    % heatmap
    exp = zscore(m_cluster_expression')';
    xticklength=length(exp(1,:));
    [i,xi]=sort(m_gene_cluster);
    figure;
    %draw heatmap as the genes sorted;
    subplot('position',[0.25 0.1 0.65 0.85]);
    imagesc(exp(xi,:),[-2.5 2.5]);colorbar;colormap(RdBu_cm)
    h_gca = gca;
    xlim([0.5, xticklength+0.5]);
    h_gca.XTick = 1:xticklength;
    xlabel('Cell clusters');
    title('Mouse gene expression by clusters')
    %draw heatmap of gene clusters
    f=subplot('position',[0.1 0.1 0.05 0.85]);
    imagesc(i);xticks([]);yticks([]);colormap(f,genecluster_cm);
    ylabel('Gene clusters');
    
    % gene cluster plots
    j = length(unique(m_gene_cluster));
    n = unique(m_gene_cluster);
    figure;
    for i = 1:j
        subplot(j,1,i);
        y = exp(find(m_gene_cluster==n(i)),:);
        x =  1:xticklength;
        shadedErrorBar(x, mean(y,1),std(y), ...
            'lineprops',{'c-o','MarkerSize',3,...
            'Color',genecluster_facecm(i,:),...
            'MarkerFaceColor',[178,24,43]/255});
        ylim([-2 3]);
        xlim([0.5, xticklength+0.5]);
        h_gca = gca;
        h_gca.XTick = 1:xticklength;
        title(sprintf('Cluster %d (n=%d)',n(i),length(find(m_gene_cluster==n(i)))));
    end
    
    
    %Save the gene expression of stages and 
        save('YOUR_PATH\h_m_stage_expression.mat',...
        'm_cluster_expression','m_cluster_exp_meanUMI', 'm_gene_cluster', 'mouse_genename','-append');

end



