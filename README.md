# cleanCall
Correcting for DNA contamination in genotype calling

## Quick Start
   	git clone https://github.com/hyunminkang/cleancall.git
   	cd cleancall
   	./configure --prefix [/path/to/install]
   	make
   	make install
   	[/path/to/install]/cctools pileup --loci [in.sites.vcf.gz] \\
	     --index [index.file.with.sample.id.and.bam.file.txt] \\
	     --out [output.pileup] --ref [ref.fa] --run [num.parallel.jobs]
   	[/path/to/install]/cctools verify --index [output.pileup.ped] \\
	     --out [output.verify] --vcf [in.sites.vcf.gz] \\
	     --run [num.parallel.jobs]
   	[/path/to/install]/cctools genotype --invcf [in.sites.vcf.gz] \\
	     --ped [output.verify.ped] --out [output.genotype.vcf.gz] \\
	     --ref [ref.fa] --run [num.parallel.jobs]

## Introduction

cleanCall is a software package for detecting, estimating, correcting for contamination in DNA sequence data. It contains one C++ program, cleanCall, as a core binary file that runs statistical algorithms. It contains three wrapper scripts, cctools-pileup for producing pileups from sequence data, cctools-verify for detecting and estimating contamination, and cctools-genotype for correcting for contaminations. These wrapper scripts allows to run cleanCall across a large number of sequenced genomes or exomes in a highly parallel manner.

## Detailed Usage

The basic usage can be obtained by simply typing the command without any argument (e.g. `cctools-genotype`). More detailed usage can be obtained with -man option (e.g. `cctools-genotype -man`). If you have questions about cleanCall, you may send me an email at hmkang@umich.edu

## Citing cleanCall

* Flickinger M, Jun G, Abecasis GR, Boehnke M, Kang HM. Am J Hum Genet. (2015) Correcting for Sample Contamination in Genotype Calling of DNA Sequence Data. *Am J Hum Genet* **97(2)**:284-90. doi: 10.1016/j.ajhg.2015.07.002.

* Jun G, Flickinger M, Hetrick KN, Romm JM, Doheny KF, Abecasis GR, Boehnke M, Kang HM (2012) Detecting and estimating contamination of human DNA samples in sequencing and array-based genotype data, *Am J Hum Genet* **91 (5)**, 839-848

## Frequently Asked Questions (FAQs)

#### What is the purpose of the VCF file given to the -loci argument to cctools-pileup? What sort of data should this contain?

The variant site VCF file should contain VCF representing known common variant sites, such as HapMap, with allele frequency information available. For example, You can use HapMap3 variant sites available in our [GotCloud Repository](https://github.com/statgen/gotcloud/blob/master/test/chr20Ref/hapmap_3.3.b37.sites.chr20.vcf.gz)

#### What is the format of `index.file.with.sample.id.and.bam.file.txt`?

The format of index file is simply tab-delimited files with the first column as sample ID and second column as full path of the BAM file corresponding to the sample.
