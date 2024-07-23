	P=5E-8
	KB=10000
	r2=0.001
	plink=/home/software/plink/plink
	ref=/home/data/ref/EUR
	home=/home/data/summarydata
	cd /home/data/summarydata/index_snp
	for exp in AD_Rheumatoid_Arthritis_European_2022_NG AD_Rheumatoid_Arthritis_pos_European_2022_NG; do
	${plink} --bfile ${ref}/EUR_1000Genome_phase3_all --clump-snp-field SNP --clump-field P --clump ${home}/${exp}.txt.gz --clump-p1 ${P} --clump-r2 ${r2} --clump-kb ${KB}
	mv plink.clumped ${exp}_${P}_${r2}_${KB}_clump.txt
	rm plink.log
	rm plink.nosex
	done


	P=5E-8
	KB=10000
	r2=0.001
	plink=/home/software/plink/plink
	ref=/home/data/ref/EUR
	home=/home/data/summarydata
	cd /home/data/summarydata/index_snp
	for exp in Blood_Pheweb_EUR_Malignant_lymphoma Blood_finnUKBB_C3_CLL_EXALLC_GRCh37_qc Blood_finnUKBB_C3_DLBCL_EXALLC_GRCh37_qc Blood_finnUKBB_CD2_FOLLICULAR_LYMPHOMA_EXALLC_GRCh37_qc Blood_finnUKBB_CD2_HODGKIN_LYMPHOMA_EXALLC_GRCh37_qc Blood_finnUKBB_CD2_TNK_LYMPHOMA_EXALLC_GRCh37_qc Blood_finnUKBB_CD2_NONFOLLICULAR_LYMPHOMA_EXALLC_GRCh37_qc Blood_finnUKBB_CD2_NONHODGKIN_NAS_EXALLC_GRCh37_qc; do
	${plink} --bfile ${ref}/EUR_1000Genome_phase3_all --clump-snp-field SNP --clump-field P --clump ${home}/${exp}.txt.gz --clump-p1 ${P} --clump-r2 ${r2} --clump-kb ${KB}
	mv plink.clumped ${exp}_${P}_${r2}_${KB}_clump.txt
	rm plink.log
	rm plink.nosex
	done





	for exp in AD_Rheumatoid_Arthritis_European_2014_nature; do
	for P in 5E-8; do
	for KB in 10000; do
	for r2 in 0.001; do
	plink=/f/software/plink.exe
	ref=/f/reference/EUR
	home=/f/Summarydata/autoimmune/20220720qc
	cd /f/Summarydata/autoimmune/index_snp
	${plink} --bfile ${ref}/EUR_1000Genome_phase3_all --clump-snp-field SNP --clump-field P --clump ${home}/${exp}.txt --clump-p1 ${P} --clump-r2 ${r2} --clump-kb ${KB}
	mv plink.clumped ${exp}_${P}_${r2}_${KB}_clump.txt
	rm plink.log
	rm plink.nosex
	done
	done
	done
	done

	for exp in AD_Pernicious_Anemia; do
	for P in 5E-8; do
	for KB in 10000; do
	for r2 in 0.001; do
	plink=/f/software/plink.exe
	ref=/f/reference/EUR
	home=/f/Summarydata/C_data_autoimmune/20220720qc
	cd /f/Summarydata/C_data_autoimmune/index_snp
	${plink} --bfile ${ref}/EUR_1000Genome_phase3_all --clump-snp-field SNP --clump-field P --clump ${home}/${exp}.txt --clump-p1 ${P} --clump-r2 ${r2} --clump-kb ${KB}
	mv plink.clumped ${exp}_${P}_${r2}_${KB}_clump.txt
	rm plink.log
	rm plink.nosex
	done
	done
	done
	done

	for exp in AD_Asthma_adultonset_2019AJHG_maf_GRCh37 AD_Asthma_childonset_2019AJHG_maf_GRCh37 AD_Rheumatoid_Arthritis_European_2022_NG AD_Rheumatoid_Arthritis_pos_European_2022_NG; do
	for P in 5E-8; do
	for KB in 10000; do
	for r2 in 0.001; do
	plink=/f/software/plink.exe
	ref=/f/reference/EUR
	home=/f/Summarydata/C_data_autoimmune/20220720qc
	cd /f/Summarydata/C_data_autoimmune/index_snp
	${plink} --bfile ${ref}/EUR_1000Genome_phase3_all --clump-snp-field SNP --clump-field P --clump ${home}/${exp}.txt.gz --clump-p1 ${P} --clump-r2 ${r2} --clump-kb ${KB}
	mv plink.clumped ${exp}_${P}_${r2}_${KB}_clump.txt
	rm plink.log
	rm plink.nosex
	done
	done
	done
	done


	for exp in Blood_Pheweb_EUR_Malignant_lymphoma; do
	for P in 5E-8; do
	for KB in 10000; do
	for r2 in 0.001; do
	plink=/f/software/plink.exe
	ref=/f/reference/EUR
	home=/f/Summarydata/C_data_blood
	cd /f/Summarydata/C_data_blood/index_snp
	${plink} --bfile ${ref}/EUR_1000Genome_phase3_all --clump-snp-field SNP --clump-field P --clump ${home}/${exp}.txt.gz --clump-p1 ${P} --clump-r2 ${r2} --clump-kb ${KB}
	mv plink.clumped ${exp}_${P}_${r2}_${KB}_clump.txt
	rm plink.log
	rm plink.nosex
	done
	done
	done
	done

