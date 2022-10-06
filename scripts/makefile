NOW=`date --utc +%y%m%d.%H%M%SZ`

include makefile.vars

# pre-processing ---
# optional

# sort/index bams, rename 
SRC=utils
get-sorted-bams:
	$(MKCMDBASH) $(SRC)/sortnrename $(SAMPLELIST) $(BAMLOC_0) $(BAMLOC_1) $(TMPDIR) > $(LOGS)/$@.$(NOW).log 2>&1 &

# remove duplicates, re-index
get-dedups:
	$(MKCMDBASH) $(SRC)/markdups $(SAMPLELIST) $(BAMLOC_1) $(TMPDIR) > $(LOGS)/$@.$(NOW).log 2>&1 &

# 26.09.22 re-header BAM files
get-headers:
	$(MKCMDBASH) $(SRC)/reheader $(BAMLOC_1) > $(LOGS)/$@.$(NOW).log 2>&1 &


# run GRAND-SLAM ---
# optional - sbatch/bash only

# create GRAND-SLAM index
SRC=GRAND-SLAM
get-gs-index:
	sbatch $(SRC)/getindex $(GEDI) $(INDEXLOC_GS) $(INDEXLOC) $(FASTA) $(GTF) $(INDEX_GS) > $(LOGS)/$@.$(NOW).log 2>&1 &
# link BAM files to GRAND-SLAM results directory
get-bam-list:
	bash $(SRC)/linkbams $(SAMPLELIST) $(BAMLOC_1) $(BAMLOC_GS) $(OUTPUT_GS) > $(LOGS)/$@.$(NOW).log 2>&1 &
# submit job
run-gs:
	sbatch $(SRC)/rungs $(GEDI) $(BAMLOC_GS) $(INDEX_GS) $(OUTPUT_GS) $(NUM_CPUS) $(GSOPTS) > $(LOGS)/$@.$(NOW).log 2>&1 &


# run BFCtools ---
# optional

# run BCFtools
SRC=utils
run-bcftools:
	$(MKCMDBASH) callsnps $(SAMPLELIST) $(BAMLOC_1) $(BCFLOC) $(OUTPUT_BCF) $(INDEXLOC) $(FASTA) $(NUM_CPUS) > $(LOGS)/$@.$(NOW).log 2>&1 &


# workflow (prep BAM files) ---

# run splbam for all samples in config
run-workflow:
	$(MKCMDPYTHON) run-all-splbam $(YML) $(RUN_OPTS) $(ALL_OPTS) --log-file $(LOGS)/$@.$(NOW).out &

# plot mismatch details
SRC=R/analysis
plot-mm:
	$(MKCMDR) $(SRC)/plot_mismatches.R $(RESLOC) > $(LOGS)/$@.$(NOW).out 2>&1 &


# abundance estimation (prep count tables) ---

# read counting - featureCounts
# make count-fc
fc=$(BAMLOC_1)/*.bam
RES=$(RESLOC)/tables/featureCounts/$*$(OUTPUT_FC)
count-%:
	$(MKCMDBASH) featureCounts $(FC_OPTS) -T $(NUM_CPUS) --tmpDir $(TMPDIR) -a $(INDEXLOC)/$(GTF) -o $(RES) $($*) > $(LOGS)/$@.$(NOW).log 2>&1 &

# read counting - Salmon	
SRC=Salmon
SALMON_PATH := $(or $(SALMON_PATH),$(shell which salmon))
GFFREAD_PATH := $(or $(GFFREAD_PATH),$(shell which gffread))

# Salmon - generate index (if using mapping-based mode)
get-salmon-index:
	$(MKCMDBASH_HMEM) $(SRC)/getindex $(SALMON_PATH) $(GFFREAD_PATH) $(INDEXLOC) $(FASTA) $(GTF) $(INDEX_SALMON) $(NUM_CPUS)

# Salmon - convert BAM files to FASTQ
run-bam2fastq:
	$(MKCMDBASH) $(SRC)/runbam2fastq $(RESLOC) $(NUM_CPUS)

# TODO: read counting - Salmon in alignment-based mode

# read counting - Salmon in (quasi-)mapping-based mode 
run-salmon-qm:
	$(MKCMDBASH) $(SRC)/runsalmonqm $(SALMON_PATH) $(RESLOC) $(INDEX_SALMON) $(SALMON_OPTS) $(NUM_CPUS)
	
	
# pulseR ---

SRC=R/pulser
# prepare pulseR input from featureCounts tables
prep-fc:
	${MKCMDR} $(SRC)/prep_fc.R $(RESLOC) $(RESLOC_PULSER) $(OUTPUT_FC) > $(LOGS)/$@.$(NOW).log 2>&1 &
# prepare pulseR input from Salmon tables
# prep-salmon:	
	
fc=$(RESLOC_PULSER)/featureCounts
salmon=$(RESLOC_PULSER)/salmon
# all time points
# make pulse-all-fc or pulse-all-salmon 
pulse-all-%:
	${MKCMDRLONG} $(SRC)/fit_rtc.R $(SRC) $($*) > $(LOGS)/$@.$(NOW).out 2>&1 &
# subset of time points
# make pulse-set-fc or pulse-set-salmon 
pulse-set-%:
	${MKCMDRLONG} $(SRC)/fit_rtc.R $(SRC) $($*) 1 > $(LOGS)/$@.$(NOW).out 2>&1 &	
	


