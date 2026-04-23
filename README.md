# meta_viena
Pipelines to obtain reads, proteins and metapeptides 

## Filter out human contamination from metagenomic reads

1. First, check your Illumina reads using **FASTQC**, determine your read length and how many low quality bases are both on the front and tail of the reads.

2. Run **fastp**, we will assumme that both the front and tail have 5 low quality bases, and that our reads are 150 bp, using 16 threads. 

```
fastp --in1 [R1] --in2 [R2] --out1 [trimmed_R1] --out2 [trimmed_R2] -z 9 -V -f 5 -F 5 -D --dup_calc_accuracy 6 -l 140 -c -w 16
```

3. Map the trimmed reads to the GRCh38 human genome reference (GCF_000001405.40) using **bwa-mem2**

```
bwa-mem2 mem -o [mapping.sam] -t 16 [GRCh38 reference] [trimmed_R1] [trimmed_R2]
```

4. Filter out the human contamination using **exclude_human_illumina.pl** and **SAMTools**

```
perl /media/databases/exclude_human_illumina.pl [mapping.sam]
```

5. Your reads are now ready for further analysis and/or assembly. Delete [mapping.sam]

## Assemble your reads using PLASS and filter out partial proteins

1. Assuming your reads are already trimmed and filtered out of human contamination, run **PLASS**

```
plass assemble [R1] [R2] [assembly.faa] [temp_dir] --remove-tmp-files
```

2. Filter out partial proteins using **plass_orf_filter.pl**

```
perl plass_orf_filter.pl [assembly.faa] [PREFIX] [minimun protein length] [maximum protein length]
```

3. After the script is done, you will have 3 output files in your current directory. [PREFIX]_complete_orfs.faa is the one appropiate for further usage.
	- [PREFIX]_complete_orfs.faa
	- [PREFIX]_semi-partial_orfs.faa
	- [PREFIX]_partial_orfs.faa

## Assemble your reads using metaSPAdes and filter out short contigs

1. Run **metaSPAdes** using your filtered reads

```
metaspades.py -o [metaspades_assembly] -1 [R1] -2 [R2] -t 16 -m [RAM in Gb]
```

2. Filter out short contigs using **filter_by_length.pl** (should be >1000 bp)

```
perl filter_by_length.pl metaspades_assembly/contigs.fasta [minimun contig length] [maximun contig length] > [filtered_assembly.fa]
```

## Assemble your reads using PenguiN and filter out short contigs

1. Run **PenguiN** using your filtered reads

```
penguin guided_nuclassemble [R1] [R2] [output] [temp] --remove-tmp-files
```

2. Filter out short contigs using **filter_by_length.pl** (should be >1000 bp)

```
perl filter_by_length.pl penguin_contigs.fasta [minimun contig length] [maximun contig length] > [filtered_assembly.fa]
```

## Get ORFs from the metaspades and penguin assemblies, and translate them into proteins using Prodigal

1. Run Prodigal

```
prodigal -a [proteins_output.faa] -c -f gff -i [assembly_input.fasta] -m -o [prodigal_output.gff] -p meta
```

## Derreplicate your proteins into a non redundant set

1. Concatenate all proteins

```
cat plass_proteins.faa metaspades_proteins.faa penguin_proteins.faa > all_proteins.faa
```

2. Run MMSeqs2 in linclust mode to obtain a NR protein set

```
mmseqs createdb all_proteins.faa all_DB
mmseqs linclust all_DB all_DB_clu tmp --alignment-mode 3 --min-seq-id 0.99 -c 0.99 --cov-mode 0 --cluster-mode 2 --threads 20 --realign --remove-tmp-files
mmseqs createsubdb all_DB_clu all_DB all_DB_clu_rep
mmseqs convert2fasta all_DB_clu_rep all_nr_proteins.faa
```

## Make a metapeptides NR database

1. Recover all ORFs contained in the forward reads using orfM

```
orfm [R1] > orfs_1.fa
````

2. Recover all ORFs contained in the reverse reads using orfM

```
orfm [R2] > orfs_2.fa
````

3. Rename and concatenate orfM ORFs

```
perl change_headers_orfm.pl orfs_1.fa orfs_2.fa > orfM_orfs.faa
```

4. Run sixgill with your reads to get additional metapeptides

```
sixgill_build --minlength 10 --minlongesttryppenlen 7 --minqualscore 30 --minreadcount 2 --outfasta [sixgill_peptides.faa] --nogzipout --minorflength 40 --out [sixgill_peptides_DB]
```

5. Concatenate all metapeptides and make a csv

```
cat orfM_orfs.faa sixgill_peptides.faa > all_metapeptides.faa
echo 'seq_name,peptide' > csv_header
grep ">" all_metapeptides.faa > all_headers
sed -i 's/>//g' all_headers
grep -v ">" all_metapeptides > all_seqs
paste -d "," all_headers all_seqs > all_pre_table
cat csv_header all_pre_table > all_metapeptides.csv
rm csv_header all_headers all_seqs all_pre_table
```

6. Randomize and derreplicate metapeptides

```
python randomize_peptides_PA_v3.py -i all_metapeptides.csv -o all_random_metapeptides.csv --log randomization_log --seed 139808 --chunksize 1000000
```

7. Make final NR metapeptides fasta files for TDA

```
perl make_TDA_fasta_from_random_out.pl all_random_metapeptides.csv [TDA_fasta_prefix] [Targets_prefix] [Decoys_prefix] 30000000
``` 

## Quantify genes/proteins/peptides in copies per cell and get related taxonomies

1. Run metaSPAdes using your reads

```
metaspades.py -o [metaspades_assembly] -1 [R1] -2 [R2] -t 16 -m [RAM in Gb]
```

2. Filter out short contigs using filter_by_length.pl

```
perl filter_by_length.pl metaspades_assembly/contigs.fasta [minimun contig length] [maximun contig length] > [filtered_assembly.fa]
```

3. Bin the contigs into MAGs using any tool(s) you prefer. Personal reccomendations: MetaBAT2, VAMB, SemiBin2, Metadecoder, Binny and COMEBin

4. Run **CAT** with metaSPAdes assembled contigs to get contig taxonomic classifications

```
CAT_pack contigs -c metaspades_contigs.fa -d cat_database/db/ -t cat_database/tax/ --no_stars -n 16 --sensitive --block_size 2 --tmpdir temp -o contigs_CAT
```

5. Run **BAT** with the MAG collection (the MAGs should share the extension e.g. .fa and be in the same directory) to get MAG taxonomic classifications

```
CAT_pack bins -b [MAG dir] -d cat_database/db/ -t cat_database/tax/ -s [.fa] -o MAGs_BAT --no_stars -n 16 --sensitive --block_size 6 --tmpdir temp/
```

6. Use **ARGs_OAP** to build a database from peptides multifasta (DDA or DIA peptides) 

```
args_oap make_db -i peptides.faa
```

7. Run **ARGs_OAP** to get reads that map to the proteins to quantify. 

- Read files must be in a single directory and be called sample_name_1.fastq.gz and sample_name_2.fastq.gz
- Create a structure file for your peptides (tab-separated), where the first column has the header Peptide and all peptide names, and any ammount of extra columns (at least one) where cluster types are specicified
> e.g.

| Peptide | Cluster |
| ----------- | ----------- |
| pep_A | cluster_A |
| pep_B | cluster_A |
| pep_C | cluster_B |
| pep_D | cluster_B |

```
args_oap stage_one -i [reads_directory] -o [args_oap_out] -t 16 -f fastq --database peptides.faa
args_oap stage_two -i [args_oap_out] -t 16 --database peptides.faa --structure1 peptides_structure.txt
```

8. Run **fasta_to_fastq.pl** to recover the reads mapping to the peptides from the args_oap output

```
perl fasta_to_fastq.pl args_oap_out/extracted.filtered.fa [RAT_1.fastq] [RAT_2.fastq] [RAT_single.fastq]
```

9. Get taxonomic classifications from the paired end reads using **RAT**

```
CAT_pack reads -c metaspades_contigs.fa -b [MAG dir] -s [.fa] -t cat_database/tax/ -m mcr -o [RAT_paired] -1 [RAT_1.fastq] -2 [RAT_2.fastq] -d cat_database/db/ --no_stars -n 16 --sensitive --block_size 2 --tmpdir temp --c2c contigs_CAT.contig2classification.txt --b2c MAGs_BAT.bin2classification.txt
```

10. Use the modified CAT_pack script included in this repo to run RAT using single end reads

```
CAT_pack reads -c metaspades_contigs.fa -b [MAG dir] -s [.fa] -t cat_database/tax/ -m mcr -o [RAT_single] -1 [RAT_single.fastq] -d cat_database/db/ --no_stars -n 16 --sensitive --block_size 2 --tmpdir temp/ --c2c contigs_CAT.contig2classification.txt --b2c MAGs_BAT.bin2classification.txt
```

11. Run **SingleM** in microbial_fraction mode using all the metagenomic reads (NOT THE RAT READS) to get the estimated prokaryote genome size and number of prokaryotic bases

```
singlem pipe -1 [sample_name_1.fastq.gz] -2 [sample_name_2.fastq.gz] -p [taxonomic_profile] --otu-table [OTU_table] --threads 16
singlem microbial_fraction -1 [sample_name_1.fastq.gz] -2 [sample_name_2.fastq.gz] -p [taxonomic profile] > sample_name_smf
```

12. Create a **MMSeqs2** database for your peptides

```
mmseqs createdb peptides.faa peptides_DB
```

13. Run **MMSeqs2** using the RAT reads as queries against the proteins

```
mmseqs easy-search [RAT_1.fastq] peptides_DB [RAT_1.m8] tmp --alignment-mode 3 -s 7 --format-output "query,target,pident,qcov,tcov,evalue,bits,qlen,tlen,alnlen" --remove-tmp-files
mmseqs easy-search [RAT_2.fastq] peptides_DB [RAT_2.m8] tmp --alignment-mode 3 -s 7 --format-output "query,target,pident,qcov,tcov,evalue,bits,qlen,tlen,alnlen" --remove-tmp-files
mmseqs easy-search [RAT_single.fastq] peptides_DB [RAT_single.m8] tmp --alignment-mode 3 -s 7 --format-output "query,target,pident,qcov,tcov,evalue,bits,qlen,tlen,alnlen" --remove-tmp-files
```

14. Concatenate files

```
cat RAT_*.m8 > reads.m8
cat *read2classification.txt > read2classification.txt
```

15. Run **get_abundance_and_taxonomy.pl** to quantify peptides in copies/cell and get taxonomic classifications for each peptide based on the LCA of all reads that map to them. The RAT_otu_table output is useful if you want to calculate alpha-diversity metrics for each peptide afterwards.

```
perl get_abundance_and_taxonomy_v2.pl --smf [sample_name_smf] --r2c [read2classification.txt] --m8 [reads.m8] --otu_table [RAT_otu_table] --max_evalue 100 --min_pident 100 --min_qcov 0 --min_alnlen 7 --min_tcov 1 > RAT_abundance_and_tax
```
