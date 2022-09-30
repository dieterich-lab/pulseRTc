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
	sbatch $(SRC)/getindex $(GEDI) $(INDEXLOC_GS) $(INDEXLOC) $(FASTA) $(GTF) $(INDEX_GS)
# link BAM files to GRAND-SLAM results directory
get-bam-list:
	bash $(SRC)/linkbams $(SAMPLELIST) $(BAMLOC_1) $(BAMLOC_GS) $(OUTPUT_GS)
# submit job
run-gs:
	sbatch $(SRC)/rungs $(GEDI) $(BAMLOC_GS) $(INDEX_GS) $(OUTPUT_GS) $(NUM_CPUS) $(GSOPTS)
	
	
# TODO run bcftools
	
	
# workflow ---
# TODO dealwith RUN_OPTS to config
CFG=run
run-workflow:
	$(MKCMDPYTHON) run $(CFG)/config.yaml $(RUN_OPTS) --trim5p 5 --trim3p 5 $(ALL_OPTS) --log-file $(LOGS)/$@.$(NOW).out &
	
# plot mismatch details 
plot-mm:
	$(MKCMDR) workflow/mismatches/plot_mismatches.R default > $(LOGS)/$@.$(NOW).out 2>&1 &
	
# read counting - featureCounts
# make count-fc
fc=$(BAMLOC_1)/*.bam
RES=$(RESLOC)/tables/$*$(OUTPUT_NAME)
count-%:
	$(MKCMDBASH) featureCounts -t exon -s 2 -p -B -C -M --fraction --ignoreDup -T $(NUM_CPUS) --tmpDir $(TMPDIR) -a $(INDEXLOC)/$(GTF) -o $(RES) $($*) > $(LOGS)/$@.$(NOW).log 2>&1 &

	
	



	
	
	
# pulseR ---

PSRC=pulser
PLOG=pulser/logs
prep:
	${MKCMDR} $(PSRC)/prep.R > $(PLOG)/$@.$(NOW).log 2>&1 &
	
# all time points
pulse-all:
	${MKCMDRLONG} $(PSRC)/fit_rtc.R all > $(PLOG)/$@.$(NOW).out 2>&1 &
# subset of time points - 3
pulse-all-sets-1:
	${MKCMDRLONG} $(PSRC)/fit_rtc.R all 1 > $(PLOG)/$@.$(NOW).out 2>&1 &
# subset of time points - 4
pulse-all-sets-2:
	${MKCMDRLONG} $(PSRC)/fit_rtc.R all 2 > $(PLOG)/$@.$(NOW).out 2>&1 &	

