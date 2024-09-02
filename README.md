# PrimerDesigner
PrimerDesigner designs nested primers for either NGS or various PCR applications. It requires one of two inputs: a vcf or a csv with raw sequence 
(and can handle both simultaneously, if you so choose), a reference genome, a chromsome lengths file. These primer designs are mostly predicated 
on circulating tumor DNA lengths. The bedtools slop function can be adjusted for more design space, but will likely require tweaking of the Primer3 parameters.

For the -d flag, you designate the design logic for either:
1 ) NGS; This option tiles primers around your ROI, intending to amplify your region of choice.
2 ) PCR; This option flanks your ROI while forcing either the forward or reverse primer for the nested pair of primers to have some degree of overlap with your ROI. The intention is your PCR reaction will not occur if your breakpoint is present (i.e., only amplification of WT). This is intended for certain RT-PCR or ddPCR strategies.
3 ) Internal Probe PCR; This option designs flanking primers around your ROI will dropping an internal probe over your target region. This is primarily for RT-PCR or certain ddPCR applications. 

------------------

The primer designer program relies on the following dependencies. Download these programs from the provided URL. 

Primer3 (https://github.com/primer3-org/primer3/releases)
MFEPrimer3.1 (https://github.com/quwubin/MFEprimer-3.0/releases)
BLASTn2.2.18 (https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/)
bedtools (https://bedtools.readthedocs.io/en/latest/)

All of these programs must be installed into your $PATH variable to be called within the script. Here are two examples of how to install a program into your PATH variable:

------------------

Example 1:
```
cd src # Navigate to directory containing the executables
sudo cp primer3_core /usr/local/bin # Allow privileges and copy to your path
sudo cp -r primer3_config /opt/ # Allow privileges and copy to your path recursively, if needed
```
------------------

Example 2:
```
cd src
sudo cp mfeprimer /usr/local/bin
```
------------------

You may use any reference genome for the blast query. Hg38 is recommended for the current purposes of this script. This can be installed using this command:

`wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz`

Additionally, you need a chromosome length file which can be acquired using this command:

`wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes `

Variant call format, or vcf, are the typical output of an NGS pipeline. The intention of the NGS portion of this script is to create more streamlined panels following a investigative sequencing event. Here is a link to a representative vcf file:

https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/

If you do not have wget, download this program from here: (https://ftp.gnu.org/gnu/wget/?C=M;O=D), or (better) using Homebrew (https://brew.sh/).

------------------

Future additions: Add further PCR option granularity in regards to targeting mutations.

------------------

Please direct any questions, bugs, or comments to: hummer.blake@gmail.com
