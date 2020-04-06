<h1>FuSe</h1>
<h2>Functional grouping of transcripts for RNA-Seq analyses</h2>

After cloning the repository, first run the following commands

On Linux machine
  *cat bi_do.part0\* > ~/bi_do.rar*
  *unrar e bi_do.ra*

On Windows, unrar using Winrar
Click on any .rar file  (bi_do.part0*) and click extract here


0.	**Environment setup**: FuSe requires Perl (v5.26.1 or higher) and certain Perl modules. When using FuSe for the first time, install the Perl modules by running the following command from command prompt.

```cpan File::Find::Rule Storable List::Util List::MoreUtils Array::Utils Data::Dumper Getopt::Long```

Another prerequisite for FuSe is a **relationship file** for Ensembl gene, transcript and protein ids; which can be obtained from Biomart; refer *FuSe/data/sample_biomart.txt*.


1.	**Protein pairs**: All the protein pair combinations are created and both confidence scores (KS and DS) are calculated for them. If you have used the Ensmebl Homo_sapiens.GRCh38 for alignment, then the precomputed data object ‘bi_do’ can be used for analyses.

```perl /path/to/script/cal_pp_conf.pl --rel /path/to/file/Biomart_rel.txt --bi_do /path/to/data_object/bi_do.data --ss /path/to/ss/scoring_scheme.txt --out_path /path/to/outfile/ --pp_do prot_pairs.data```


2.	**SFPGs**: The overlapping protein pairs which are over the given CSC are used to create SFPGs. The SFPG confidence is calculated by averaging the protein pair scores. One SFPG is formed for each protein coding transcript that has other similar protein coding transcripts.

```perl /path/to/script/make_sfpgs.pl --pp_do /path/to/file/prot_pairs.data --score_type KS --csc 95 --out_path /path/to/outfile/ --name sfpg sfpg.data```


3.	**SFPGs expression**: Using the FPKM values for the samples and SFPG data, expression is calculated for all SFPGs.

```perl /path/to/script/recal_expression.pl --input /path/to/file/exp_file.txt --sfpg /path/to/file/sfpg.data --out_path /path/to/outfile/ --recal recal_exp.txt```

<h3>Creating your own BLAST Interpro data object (bi_do)</h3>


4.	**Data preparation**: Sequence alignment, protein domain, motifs and family information are required to make the bi_do. Protein sequences can be obtained from Ensembl. BLAST+ [16] and Interpro [17] were run as explained in their user manuals. For BLAST+, first a protein reference database was created using all sequences and then aligned to all sequences. The alignment results were then obtained in out format 7.

```-outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"```

The results from Interpro were obtained in tsv format using the default values. To use in future available updates to Ensembl, the data should be generated as explained in this step.


5.	**Preprocessing**: The data generated from BLAST+ and Interpro were put together to create the bi_do.

```perl /path/to/script/preprocessing_data.pl --interpro path/to/file/interpro.tsv --blast /path/to/file/blast.txt --out_path /path/to/outfile/ --bi_do bi_do.txt```
