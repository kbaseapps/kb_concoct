# kb_concoct release notes
=========================================

ToDo
-----
* Add kb_parallels support - see kbase-metabat example
* Add different ways to extract coverage (e.g. jgi_summarize_bam_contig_depths vs coverM, vs ?)
* Update documentation and inline comments

1.3.4
-----
* samtools and bbmap now installed via docker
* bbmap now has multiple modes of operation (default, fast, and very-sensitive)

1.3.3
-----
* Seed parameters set for all read mappers except bbmap.sh
* Removed BWA-MEM read mapper because no seed parameter available
* Changed to use a new deinterlace reads script

1.3.2
-----
* Changed to not display binned contig output widget

1.3.1
-----
* Added HISAT2 read aligner as mapping option
* Large code de-linting effort
* Added tests for different read mapping tools
* Added bowtie2 very sensitive mode

1.3.0
-----
* Refactored most code for kb_concoctTest, streamlined variable and function usage
* Added new method (stats.sh) to populate summary files, exposed many other statistics for genome bin summary files
* Removed bash script dependency
* Added support for multiple read file inputs
* Added auto-detection of read file type for read mapping - single vs paired-end interleaved file

1.2.0
-----
* Added multiple read mapping software (bowtie2, minimap2, and bwa), in addition to existing bbmap as selectable options

1.1.0
-----
* Added parameter options for: no_cov_normalization, no_total_coverage, and total_percentage_pca

1.0.0
-----
* Initial release

0.0.0
-----
* Module created by kb-sdk init
