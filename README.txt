####################
The collected files and scripts all relate to the manuscript as detailed below:


Title: Widespread transcriptional scanning in testes modulates gene evolution rates
Authors: Bo Xia1,4 Yun Yan1, Maayan Baron1, Florian Wagner1, Dalia Barkley1, Marta Chiodin1, Sang Y. Kim5, David L. Keefe3, Joseph P. Alukal3, Jef D. Boeke2,4 and Itai Yanai1,2
Affiliations:
1 Institute for Computational Medicine, NYU Langone Health, New York, NY 10016, USA
2 Department of Biochemistry and Molecular Pharmacology, NYU Langone Health, New York, NY 10016, USA
3 Department of Obstetrics and Gynecology, NYU Langone Health, New York, NY 10016, USA
4 Institute for Systems Genetics, NYU Langone Health, New York, NY 10016, USA
5 Department of Pathology, NYU Langone Health, New York, NY 10016, USA

Cite: Xia, B., Yan, Y., Baron, M., Wagner, F., Barkley, D., Chiodin, M., Kim, S.Y., Keefe, D.L., Alukal, J.P., Boeke, J.D. and Yanai, I., 2020. Widespread transcriptional scanning in the testis modulates gene evolution rates. Cell, 180(2), pp.248-262.




####################
CONTENTS:


######
Folders:
\Data              : files containing files related the scRNA-seq and germline SNV analysis.
\MATLAB_codes	   : cross species data between human and apes
\SNV		   : files containing the germline variants extracted from the public database.


###
\DATA\
	\HUGO_Gene_N_C	   : large gene family member lists extracted from HUGO Gene Nomenclature Committee
	\Monocle2_results  : Monocle2_pseudotime inference results
Other data:
	Apes_to_human_Id_dNdS.mat
	GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.mat
	human_mouse_orthologue_pairs_20171018.mat
	human_ensembl90_gene_features.txt
	human_genes.tsv
	mouse_ensembl90_gene_features.txt
	mouse_genes.tsv


###
\MATLAB_codes\
Main scripts:
	1_TS_code_for_human_spermatogenesis.m
	2_TS_code_for_mouse_spermatogenesis.m
	3_TS_code_for_human_mouse_spermatogenesis_comparison.m
	4_TS_code_for_evolution_collected.m
	5_TS_code_for_other_datasets.m
Related functions:
	F_calculate_snv_asymmetry_score.m
	F_cluster_gene.m
	F_cluster_sorting_by_ZAVIT.m
	F_compare_divergence.m
	F_compare_paragene_variants.m
	F_compare_snv_192.m
	F_compare_snv_asymmetry.m
	F_compare_snv_asymmetry_v2.m
	F_compare_variants.m
	F_inDrop_filter.m
	F_kmeans.m
	F_normalize.m
	F_triref_asymmetry_score_scatter.m
	F_tsv_import.m
	knn_smooth.m
Third party functions:
	\cbrewer
	cbrewer.m
	barwitherr.m
	cell2csv.m
	kmeans_opt.m
	plot_areaerrorbar.m
	rgb2hex.m
	scatter_kde.m
	shadedErrorBar.m
	Violin.m
	violinplot.m
	zavit_angle.m


###
\SNV\
	\Processing_SNV_scripts  : Contains scripts for processing raw germline mutations into clean SNVs.
	\Counting_SNV_scripts    : Contains scripts for counting the SNVs with/without considering adjacent bases.
	\human_gene_segmentation : Contains bed files of human gene segments (gene body, upstream, downstream, intron).
	\mouse_gene_segmentation : Contains bed files of mouse gene segments (gene body, upstream, downstream, intron).
	\de_novo_SNVs		 : Contains SNV counts using de novo mutation data from trio-WGS results.
	human_1000G_SNV_stat.txt
	human_SNV_stat.txt
	mouse_SNV_stat.txt
	human_1000G_genebody_12_result.12.stranded.tsv
	human_downstream_12_result.12.stranded.tsv
	human_ensembl90.downstream.singleFreq.tsv
	human_ensembl90.genebody.singleFreq.tsv
	human_ensembl90.intron.singleFreq.tsv
	human_ensembl90.intron.singleFreq.unstranded.tsv
	human_ensembl90.intron.triFreq.tsv
	human_ensembl90.upstream.singleFreq.tsv
	human_gene_body_result.12.stranded.tsv
	human_intron_12_result.12.stranded.tsv
	human_intron_12_result.12.unstrand.tsv
	human_intron_192_result.192.stranded.tsv
	human_intron_192_result.192.unstrand.tsv
	human_upstream_12_result.12.stranded.tsv
	mouse_downstream_12_result.12.stranded.tsv
	mouse_downstream_12_result.12.stranded.tsv
	mouse_ensembl90.intron.singleFreq.tsv
	mouse_ensembl90.intron.singleFreq.unstranded.tsv
	mouse_ensembl90.intron.triFreq.tsv
	mouse_ensembl90.upstream.singleFreq.tsv
	mouse_intron_12_result.12.stranded.tsv
	mouse_intron_12_result.12.unstrand.tsv
	mouse_intron_192_result.192.stranded.tsv
	mouse_intron_192_result.192.unstrand.tsv
	mouse_upstream_12_result.12.stranded.tsv
	

####################
Contact:
Bo Xia:     Bo.Xia@nyulangone.org; xiabo90@gmail.com
Itai Yanai: Itai.Yanai@nyulangone.org
2019-08-21



