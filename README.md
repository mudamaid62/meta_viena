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

## Get ORFs from the metaspades and penguin assemblies, and tranlate them into proteins using Prodigal

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


