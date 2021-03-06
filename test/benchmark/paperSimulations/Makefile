# Makefile to run simulation for paper

SHELL := /bin/bash

HGCHR1=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz
STAPH=http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=embl&id=FN433596&format=fasta&style=raw
FACAT=../../../libs/seq_file/bin/facat
SAMTOOLS=samtools

# PhiX.1.fq.gz needs to be downloaded
PROFILE_READS=../PhiX.1.fq.gz

NAMES=HsOneErrorReads100 StaphOneErrorReads100 \
			HsOneErrorReads250 StaphOneErrorReads250 \
			HsOneNoErrorReads100 StaphOneNoErrorReads100 \
			HsOneNoErrorReads250 StaphOneNoErrorReads250 \
			HsManyErrorReads100 StaphManyErrorReads100 \
			HsManyErrorReads250 StaphManyErrorReads250

SCRIPTS=$(addsuffix /run-sim.sh, $(NAMES))

REFS=chr1/chr1.fa.gz chr1/chr1.fa chr1/chr1.fa.fai \
     staph/staph.fa.gz staph/staph.fa staph/staph.fa.fai

# Number of events
HsOneSNPs=225000
HsOneIndels=25000
HsManySNPs=798242
HsManyIndels=88694

StaphOneSNPs=27000
StaphOneIndels=3000
StaphManySNPs=202081
StaphManyIndels=22453

all: $(PROFILE_READS) $(REFS) $(SCRIPTS)

clean:
	@echo rm -rf staph chr1
	@echo rm -rf $(NAMES)

test:
	echo $$@

$(PROFILE_READS):
	@echo Please fetch $(PROFILE_READS) from Illumina or pass other with PROFILE_READS=
	exit 1;

chr1/chr1.orig.fa.gz:
	mkdir -p chr1
	wget -O chr1/chr1.orig.fa.gz '$(HGCHR1)'

staph/staph.orig.fa.gz:
	mkdir -p staph
	wget -O staph/staph.orig.fa '$(STAPH)'
	gzip staph/staph.orig.fa

%.fa: %.orig.fa.gz
	zcat $< | tr -d 'N' | $(FACAT) -w 50 > $@

%.fa.gz: %.fa
	gzip -c $< > $@

%.fa.fai: %.fa
	$(SAMTOOLS) faidx $<

#$* is {Hs,Staph}{One,Many}{Error,NoError}Reads{100,250}
$(NAMES): $(SCRIPTS)
%/run-sim.sh: Makefile
	mkdir -p $*
	@set -e; \
	species=`echo $* | grep -oE '(Hs|Staph)'`; \
	quant=`echo $* | grep -oE '(One|Many)'`; \
	err=`echo $* | grep -oE '(NoError|Error)'`; \
	readlen=`echo $* | grep -oE 'Reads[0-9]*' | sed 's/Reads//g'`; \
	errprof=; \
	if [ $$species == 'Hs' ]; then \
		memwidth=26; memheight=26; \
		ref=../chr1/chr1.fa; ploidy=2; allelecovg=30; subrate=0.001; \
		if [ $$quant == 'Many' ]; then num=10; snps=$(HsManySNPs); indels=$(HsManyIndels); \
		else num=1; snps=$(HsOneSNPs); indels=$(HsOneIndels); fi; \
	else \
		memwidth=21; memheight=21; \
		ref=../staph/staph.fa; ploidy=1; allelecovg=100; subrate=0.01; \
		if [ $$quant == 'Many' ]; then num=1000; snps=$(StaphManySNPs); indels=$(StaphManyIndels); \
		else num=1; snps=$(StaphOneSNPs); indels=$(StaphOneIndels); fi \
	fi; \
	if [ $$err == 'Error' ]; then errprof="ERRPROF=../$(PROFILE_READS)"; fi; \
	( echo "make -f ../../calling-comparison.mk" \
	  echo "SEQ=$$ref NUM_INDIVS=$$num PLOIDY=$$ploidy KMER=31" \
	  echo "SNPS=$$snps INDELS=$$indels INV=0 INVLEN=100" \
	  echo "READLEN=$$readlen MPSIZE=250 ALLELECOVG=$$allelecovg $$errprof" \
	  echo "MEMWIDTH=$$memwidth MEMHEIGHT=$$memheight" \
	  echo "MAPARGS='--substitutionrate=$$subrate' "'$$@' ) > $*/run-sim.sh
	chmod +x $*/run-sim.sh

.PHONY: all clean test
