NOW=`date --utc +%y%m%d.%H%M%SZ`

include makefile.vars
SRC_UTILS=utils
SRC_GS=GRAND-SLAM
SRC_R=R/analysis
SRC_SALMON=Salmon
SRC_PULSER=R/pulser
SRC_ANALYSIS=R/analysis

# pre-processing ---
# optional

# sort/index bams, rename 
get-sorted-bams:
	$(MKCMDBASH) $(SRC_UTILS)/sortnrename $(SAMPLELIST) $(BAMLOC_0) $(BAMLOC_1) $(TMPDIR) > $(LOGS)/$@.$(NOW).log 2>&1 &

# remove duplicates, re-index
get-dedups:
	$(MKCMDBASH) $(SRC_UTILS)/markdups $(SAMPLELIST) $(BAMLOC_1) $(TMPDIR) > $(LOGS)/$@.$(NOW).log 2>&1 &

# 26.09.22 re-header BAM files
get-headers:
	$(MKCMDBASH) $(SRC_UTILS)/reheader $(BAMLOC_1) > $(LOGS)/$@.$(NOW).log 2>&1 &


# run GRAND-SLAM ---
# optional - sbatch/bash only

# create GRAND-SLAM index
get-gs-index:
	sbatch $(SRC_GS)/getindex $(GEDI) $(INDEXLOC_GS) $(INDEXLOC) $(FASTA) $(GTF) $(INDEX_GS) > $(LOGS)/$@.$(NOW).log 2>&1 &
# link BAM files to GRAND-SLAM results directory
get-bam-list:
	bash $(SRC_GS)/linkbams $(SAMPLELIST) $(BAMLOC_1) $(BAMLOC_GS) $(OUTPUT_GS) > $(LOGS)/$@.$(NOW).log 2>&1 &
# submit job
run-gs:
	sbatch $(SRC_GS)/rungs $(GEDI) $(BAMLOC_GS) $(INDEX_GS) $(OUTPUT_GS) $(NUM_CPUS) $(GSOPTS) > $(LOGS)/$@.$(NOW).log 2>&1 &


# run BFCtools ---
# optional

# run BCFtools
run-bcftools:
	$(MKCMDBASH) $(SRC_UTILS)/callsnps $(SAMPLELIST) $(BAMLOC_1) $(BCFLOC) $(OUTPUT_BCF) $(INDEXLOC) $(FASTA) $(NUM_CPUS) > $(LOGS)/$@.$(NOW).log 2>&1 &


# workflow (prep BAM files) ---

# run splbam for all samples in config
run-workflow:
	$(MKCMDPYTHON) run-all-splbam $(YML) $(RUN_OPTS) $(ALL_OPTS) --log-file $(LOGS)/$@.$(NOW).out &

# plot mismatch details
plot-mm:
	$(MKCMDR) $(SRC_R)/plot_mismatches.R $(RESLOC) > $(LOGS)/$@.$(NOW).out 2>&1 &


# abundance estimation (prep count tables) ---

# read counting - featureCounts
# make count-fc
fc=$(RESLOC)/mapping/*.bam
RES=$(RESLOC)/tables/featureCounts/$*$(OUTPUT_FC)
count-%:
	$(MKCMDBASH) featureCounts $(FC_OPTS) -T $(NUM_CPUS) --tmpDir $(TMPDIR) -a $(INDEXLOC)/$(GTF) -o $(RES) $($*) > $(LOGS)/$@.$(NOW).log 2>&1 &

# read counting - Salmon	
SALMON_PATH := $(or $(SALMON_PATH),$(shell which salmon))
GFFREAD_PATH := $(or $(GFFREAD_PATH),$(shell which gffread))
# Salmon - generate index (if using mapping-based mode)
get-salmon-index:
	$(MKCMDBASH_HMEM) $(SRC_SALMON)/getindex $(SALMON_PATH) $(GFFREAD_PATH) $(INDEXLOC) $(FASTA) $(GTF) $(INDEX_SALMON) $(NUM_CPUS) > $(LOGS)/$@.$(NOW).log 2>&1 &

# Salmon - convert BAM files to FASTQ
run-bam2fastq:
	$(MKCMDBASH) $(SRC_SALMON)/runbam2fastq $(RESLOC) $(NUM_CPUS) > $(LOGS)/$@.$(NOW).log 2>&1 &

# TODO: read counting - Salmon in alignment-based mode

# read counting - Salmon in (quasi-)mapping-based mode 
run-salmon-qm:
	$(MKCMDBASH) $(SRC_SALMON)/runsalmonqm $(SALMON_PATH) $(RESLOC) $(INDEX_SALMON) $(SALMON_OPTS) $(NUM_CPUS) > $(LOGS)/$@.$(NOW).log 2>&1 &
	
	
# pulseR ---

# prepare pulseR input from featureCounts tables
prep-fc:
	${MKCMDR} $(SRC_PULSER)/prep_fc.R $(RESLOC) $(RESLOC_PULSER) $(OUTPUT_FC) > $(LOGS)/$@.$(NOW).log 2>&1 &
# prepare pulseR input from Salmon tables
# prep-salmon:	
	
fc-counts=$(RESLOC_PULSER)/featureCounts
salmon-counts=$(RESLOC_PULSER)/salmon
# all time points
# make pulse-all-fc or pulse-all-salmon 
pulse-all-%:
	${MKCMDRLONG} $(SRC_PULSER)/fit_rtc.R $(SRC_PULSER) $($*) > $(LOGS)/$@.$(NOW).out 2>&1 &
# subset of time points
# make pulse-set-fc or pulse-set-salmon 
pulse-set-%:
	${MKCMDRLONG} $(SRC_PULSER)/fit_rtc.R $(SRC_PULSER) $($*) 1 > $(LOGS)/$@.$(NOW).out 2>&1 &	
	

# downstream analyses ---

# get all estimates into tables and match pulseR and GRAND-SLAM results
# ONLY for gene quantification, i.e. featureCounts
match-gs:
	${MKCMDR} $(SRC_ANALYSIS)/get_estimates.R $(RESLOC_PULSER) $(BAMLOC_GS) $(OUTPUT_GS) $(SRC_PULSER) > $(LOGS)/$@.$(NOW).out 2>&1 &

# annotate pulseR results into a new table using biomaRt
# ONLY for gene quantification, i.e. featureCounts
ann-pulse:
	${MKCMDR} $(SRC_ANALYSIS)/annotate_pulser.R $(RESLOC_PULSER) $(SRC_PULSER) > $(LOGS)/$@.$(NOW).out 2>&1 &

