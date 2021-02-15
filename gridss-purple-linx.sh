#!/bin/bash
#
# Stand-alone GRIDSS-PURPLE-LINX pipeline
#
# Example: ./gridss-purple-linx.sh -n /data/COLO829R_dedup.realigned.bam -t /data/COLO829T_dedup.realigned.bam -v /data/colo829snv.vcf.gz -s colo829 -v /data/COLO829v003T.somatic_caller_post_processed.vcf.gz
# docker run  gridss/gridss-purple-linx
EX_USAGE=64
EX_NOINPUT=66
EX_CANTCREAT=73
EX_CONFIG=78
set -o errexit -o pipefail -o noclobber -o nounset
! getopt --test > /dev/null 
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
	echo '`getopt --test` failed in this environment.'
	exit $EX_CONFIG
fi

run_dir=/data
ref_dir=/refdata
install_dir=/opt/
tumour_bam=""
normal_bam=""
snvvcf=""
threads=$(nproc)
sample=""
normal_sample=""
tumour_sample=""
jvmheap="25g"
ref_genome_version="HG37"
validation_stringency="STRICT"
purple_args=""
amber_args=""
cobalt_args=""
linx_args=""
gripss_args=""
gridss_args=""
tumour_only="false"

usage_msg="Usage: gridss-purple-linx.sh

Required command line arguments:
	--tumour_bam: tumour BAM file
	--normal_bam: matched normal BAM file
	--sample: sample name
Optional parameters:
	--output_dir: output directory. (/data)
	--ref_genome_version: reference genome. HG37 or HG38 ($ref_genome_version)
	--ref_dir: path to decompressed Hartwig reference data package. ($ref_dir)
	--snvvcf: A somatic SNV VCF with the AD genotype field populated.
	--nosnvvcf: Indicates a somatic SNV VCF will not be supplied. This will reduce the accuracy of PURPLE ASCN.
	--threads: number of threads to use. ($threads)
	--install_dir: root directory of gridss-purple-linx release ($install_dir)
	--normal_sample: sample name of matched normal ({sample}_N) 
	--tumour_sample: sample name of tumour. Must match the somatic \$snvvcf sample name. ({sample}_T) 
	--jvmheap: maximum java heap size for high-memory steps ($jvmheap)
	--gridss_args: additional arguments to GRIDSS
	--gripss_args: additional arguments to GRIPSS
	--amber_args: additional arguments to AMBER
	--cobalt_args: additional arguments to COBALT
	--purple_args: additional arguments to PURPLE
	--linx_args: additional arguments to LINX
	--validation_stringency: htsjdk validation_stringency ($validation_stringency)
	--help: print this message and exit
"
usage() {
	echo "$usage_msg" 1>&2
	exit $EX_USAGE
}

OPTIONS=v:o:t:n:s:r:b:h
LONGOPTS=snvvcf:,nosnvvcf,output_dir:,tumour_bam:,normal_bam:,sample:,threads:,jvmheap:,ref_dir:,ref_genome_version:,normal_sample:,tumour_sample:,rundir:,validation_stringency:,install_dir:,help,gridss_args:,gripss_args:,amber_args:,cobalt_args:,purple_args:,linx_args:
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
	# e.g. return value is 1
	#  then getopt has complained about wrong arguments to stdout
	exit 2
fi
eval set -- "$PARSED"
while true; do
	case "$1" in
		--rundir)
			run_dir="$2"
			shift 2
			;;
		-v|--snvvcf)
			snvvcf="$2"
			shift 2
			;;
		--nosnvvcf)
			snvvcf="nosnvvcf"
			shift 1
			;;
		--ref_genome_version)
			ref_genome_version="$2"
			shift 2
			;;
		-n|--normal_bam)
			normal_bam="$2"
			shift 2
			;;
		-o|--output_dir)
			run_dir="$2"
			shift 2
			;;
		-t|--tumour_bam)
			tumour_bam="$2"
			shift 2
			;;
		-s|--sample)
			sample="$2"
			shift 2
			;;
		--normal_sample)
			normal_sample="$2"
			shift 2
			;;
		--tumour_sample)
			tumour_sample="$2"
			shift 2
			;;
		--threads)
			printf -v threads '%d\n' "$2" 2>/dev/null
			printf -v threads '%d' "$2" 2>/dev/null
			shift 2
			;;
		--jvmheap)
			jvmheap="$2"
			shift 2
			;;
		--install_dir)
			install_dir="$2"
			shift 2
			;;
		--ref_dir)
			ref_dir="$2"
			shift 2
			;;
		-h|--help)
			usage
			exit 1
			;;
		--validation_stringency)
			validation_stringency="$2"
			shift 2
			;;
		--gridss_args)
			gridss_args="$2"
			shift 2
			;;
		--gripss_args)
			gripss_args="$2"
			shift 2
			;;
		--amber_args)
			amber_args="$2"
			shift 2
			;;
		--cobalt_args)
			cobalt_args="$2"
			shift 2
			;;
		--purple_args)
			purple_args="$2"
			shift 2
			;;
		--linx_args)
			linx_args="$2"
			shift 2
			;;
		--)
			shift
			break
			;;
		*)
			echo "Command line parsing error ($1)"
			echo "$@"
			exit 3
			;;
	esac
done
write_status() { # Before logging initialised
	echo "$(date): $1" 1>&2
}
# $1: variable containing filename
# $2: command line argument name
assert_file_exists() {
	if [[ ! -f "$1" ]] ; then
		write_status "File $1 not found. Specify using the command line argument --$2"
		exit $EX_NOINPUT
	fi
}
assert_directory_exists() {
	if [[ ! -d "$1" ]] ; then
		write_status "Directory $1 not found. Specify using the command line argument --$2"
		exit $EX_NOINPUT
	fi
}
assert_directory_exists $install_dir/gridss "install_dir"
assert_directory_exists $install_dir/hmftools "install_dir"
assert_directory_exists $install_dir/gridss-purple-linx "install_dir"
assert_file_exists $install_dir/gridss/gridss.sh "install_dir"
assert_file_exists $install_dir/gridss/libgridss.R "install_dir"

viralreference=refgenomes/human_virus/human_virus.fa
blacklist=dbs/gridss/ENCFF001TDO.bed
repeatmasker=dbs/repeatmasker/rm.fa.out.bed
bafsnps=dbs/germline_het_pon/GermlineHetPon.vcf.gz
gcprofile=dbs/gc/GC_profile.1000bp.cnp
gridss_properties=dbs/gridss/gridss.properties
breakpoint_hotspot=dbs/knowledgebases/KnownFusionPairs.bedpe
breakend_pon=dbs/gridss_pon/gridss_pon_single_breakend.bed
breakpoint_pon=dbs/gridss_pon/gridss_pon_breakpoint.bedpe
viral_hosts_csv=dbs/knowledgebases/viral_host_ref.csv
known_fusion_csv=dbs/knowledgebases/known_fusion_data.csv
fragile_sites=dbs/knowledgebases/fragile_sites_hmf.csv
line_elements=dbs/knowledgebases/line_elements.csv
replication_origins=dbs/knowledgebases/heli_rep_origins.bed
ensembl_data_dir=dbs/ensembl_data_cache
driver_gene_panel=dbs/knowledgebases/DriverGenePanel.tsv
known_hotspots_vcf=dbs/knowledgebases/KnownHotspots.vcf.gz
rlib=rlib/
case "$ref_genome_version" in
	"HG37")
		ref_genome=refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta
		;;
	"HG38")
		ref_genome=refgenomes/Homo_sapiens.GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
		;;
	*)
		write_status "Invalid reference genome version: $ref_genome_version"
		exit $EX_CONFIG
		;;
esac
write_status "Running reference genome version: $ref_genome_version"

rlib=$ref_dir/$rlib
ref_genome=$ref_dir/$ref_genome
viralreference=$ref_dir/$viralreference
blacklist=$ref_dir/$blacklist
repeatmasker=$ref_dir/$repeatmasker
gridss_properties=$ref_dir/$gridss_properties
bafsnps=$ref_dir/$bafsnps
gcprofile=$ref_dir/$gcprofile
breakpoint_hotspot=$ref_dir/$breakpoint_hotspot
breakend_pon=$ref_dir/$breakend_pon
breakpoint_pon=$ref_dir/$breakpoint_pon
viral_hosts_csv=$ref_dir/$viral_hosts_csv
known_fusion_csv=$ref_dir/$known_fusion_csv
fragile_sites=$ref_dir/$fragile_sites
line_elements=$ref_dir/$line_elements
replication_origins=$ref_dir/$replication_origins
ensembl_data_dir=$ref_dir/$ensembl_data_dir
driver_gene_panel=$ref_dir/$driver_gene_panel
known_hotspots_vcf=$ref_dir/${known_hotspots_vcf}

if [[ "$normal_bam" == "" ]] && [[ "$normal_sample" == "" ]] ; then
	tumour_only=true
fi

if [[ "$snvvcf" == "nosnvvcf" ]] ; then
	write_status "No somatic SNV VCF supplied."
elif [[ ! -f "$snvvcf" ]] ; then
	write_status "Missing somatic SNV VCF. A SNV VCF with the AD genotype field populated is required."
	write_status "Use the script for generating this VCF with strelka if you have not already generated a compatible VCF."
	exit $EX_NOINPUT
fi
if [[ ! -f "$tumour_bam" ]] ; then
	write_status "Missing tumour BAM"
	exit $EX_NOINPUT
fi
if [[ ! $tumour_only == true ]] ; then
	if [[ ! -f "$normal_bam" ]] ; then
		write_status "Missing normal BAM"
		exit $EX_NOINPUT
	fi
fi
mkdir -p "$run_dir"
if [[ ! -d "$run_dir" ]] ; then
	write_status "Unable to create $run_dir"
	exit $EX_CANTCREAT
fi
if [[ ! -d "$ref_dir" ]] ; then
	write_status "Could not find reference data directory $ref_dir"
	exit $EX_NOINPUT
fi
if [[ ! -f "$ref_genome" ]] ; then
	write_status "Missing reference genome $ref_genome - specify with -r "
	exit $EX_NOINPUT
fi
if [[ -z "$sample" ]] ; then
	sample=$(basename $tumour_bam .bam)
fi
if [[ "$threads" -lt 1 ]] ; then
	write_status "Illegal thread count: $threads"
	exit $EX_CONFIG
fi
joint_sample_name=$sample
if [[ $tumour_only == true ]] ; then
	if [[ -z "$tumour_sample" ]] ; then
		tumour_sample=${sample}
	fi
else
	if [[ -z "$normal_sample" ]] ; then
		normal_sample=${sample}R
	fi
	if [[ -z "$tumour_sample" ]] ; then
		tumour_sample=${sample}T
	fi
fi
export R_LIBS="$rlib:${R_LIBS:-}"
base_path=$(dirname $(readlink $0 || echo $0))

### Find the jars
find_jar() {
	env_name=$1
	if [[ -f "${!env_name:-}" ]] ; then
		echo "${!env_name}"
	else
		write_status "Unable to find $2 jar. Specify using the environment variant $env_name"
		exit $EX_CONFIG
	fi
}

gridss_jar=$(find_jar GRIDSS_JAR gridss)
gripss_jar=$(find_jar GRIPSS_JAR gripss)
amber_jar=$(find_jar AMBER_JAR amber)
cobalt_jar=$(find_jar COBALT_JAR cobalt)
purple_jar=$(find_jar PURPLE_JAR purple)
linx_jar=$(find_jar LINX_JAR sv-linx)

for program in bwa samtools circos Rscript java ; do
	if ! which $program > /dev/null ; then
		write_status "Missing required dependency $program. $program must be on PATH"
		exit $EX_CONFIG
	fi
done
for rpackage in tidyverse devtools assertthat testthat NMF stringdist stringr argparser R.cache "copynumber" StructuralVariantAnnotation "VariantAnnotation" "rtracklayer" "BSgenome" "org.Hs.eg.db" ; do
	if ! Rscript -e "installed.packages()" | grep $rpackage > /dev/null ; then
		write_status "Missing R package $rpackage"
		exit $EX_CONFIG
	fi
done

if ! java -Xms$jvmheap -cp $gridss_jar gridss.Echo ; then
	write_status "Failure invoking java with --jvmheap parameter of \"$jvmheap\". Specify a JVM heap size (e.g. \"20g\") that is valid for this machine."
	exit $EX_CONFIG
fi

if [[ ! -s $ref_genome.bwt ]] ; then
	write_status "Missing bwa index for $ref_genome. Creating (this is a once-off initialisation step)"
	bwa index $ref_genome
fi

if [[ ! -s $ref_genome.bwt ]] ; then
	write_status "bwa index for $ref_genome not found."
	write_status "If you are running in a docker container, make sure refdata has been mounted read-write."
	exit $EX_NOINPUT
fi

mkdir -p $run_dir/logs $run_dir/gridss $run_dir/gripss $run_dir/amber $run_dir/purple
log_prefix=$run_dir/logs/$(date +%Y%m%d_%H%M%S).$HOSTNAME.$$

jvm_args=" \
	-Dsamjdk.reference_fasta=$ref_genome \
	-Dsamjdk.use_async_io_read_samtools=true \
	-Dsamjdk.use_async_io_write_samtools=true \
	-Dsamjdk.use_async_io_write_tribble=true \
	-Dsamjdk.buffer_size=4194304 \
	-Dsamjdk.async_io_read_threads=$threads"

write_status "run_dir=$run_dir"
write_status "ref_dir=$ref_dir"
write_status "install_dir=$install_dir"
write_status "tumour_bam=$tumour_bam"
write_status "normal_bam=$normal_bam"
write_status "snvvcf=$snvvcf"
write_status "threads=$threads"
write_status "sample=$sample"
write_status "normal_sample=$normal_sample"
write_status "tumour_sample=$tumour_sample"
write_status "jvmheap=$jvmheap"
write_status "validation_stringency=$validation_stringency"
write_status "ref_genome_version=$ref_genome_version"
write_status "rlib=$rlib"
write_status "ref_genome=$ref_genome"

gridss_dir=$run_dir/gridss
assembly_bam=$gridss_dir/$joint_sample_name.assembly.bam
gridss_driver_vcf=$gridss_dir/${tumour_sample}.gridss.driver.vcf.gz
gridss_unfiltered_vcf=$gridss_dir/${tumour_sample}.gridss.unfiltered.vcf.gz

if [[ ! -f $gridss_driver_vcf ]] ; then
	write_status "Running GRIDSS"
	if [[ "$normal_sample" != "" ]] ; then
		gridss_args=" --labels $normal_sample,$tumour_sample $gridss_args"
	else
		gridss_args=" --labels $tumour_sample $gridss_args"
	fi
	$install_dir/gridss/gridss.sh \
		-o ${gridss_driver_vcf} \
		-a $assembly_bam \
		-w ${gridss_dir} \
		-r ${ref_genome} \
		-j ${gridss_jar} \
		-t $threads \
		-b ${blacklist} \
		-c ${gridss_properties} \
		--repeatmaskerbed ${repeatmasker} \
		--jvmheap $jvmheap \
		$gridss_args \
		${normal_bam} ${tumour_bam}
	if [[ ! -f $gridss_driver_vcf ]] ; then
		write_status  "Error creating $gridss_driver_vcf. Aborting" 1>&2
		exit 1
	fi
	# GRIDSS #382 work-around
	if [[ ! -f ${gridss_driver_vcf}.tbi ]] ; then
		mv $(dirname ${gridss_driver_vcf})/gridss.tmp.$(basename ${gridss_driver_vcf}).tbi ${gridss_driver_vcf}.tbi
	fi
else
	write_status "Skipping GRIDSS - ${gridss_driver_vcf} exists"
fi

if [[ ! -f $gridss_unfiltered_vcf ]] ; then
	write_status "Running GRIDSS Annotations"
	java -Xmx8G -Dsamjdk.create_index=true \
		-Dsamjdk.use_async_io_read_samtools=true \
		-Dsamjdk.use_async_io_write_samtools=true \
		-Dsamjdk.use_async_io_write_tribble=true \
		-Dsamjdk.buffer_size=4194304 \
		-cp ${gridss_jar} gridss.AnnotateInsertedSequence \
		REFERENCE_SEQUENCE=${viralreference} \
		INPUT=${gridss_driver_vcf} \
		OUTPUT=${gridss_unfiltered_vcf} \
		ALIGNMENT=APPEND WORKER_THREADS=${threads} \
		 2>&1 | tee $log_prefix.gridss.AnnotateInsertedSequence.log
	if [[ ! -f $gridss_unfiltered_vcf ]] ; then
		write_status "Error creating $gridss_unfiltered_vcf. Aborting" 1>&2
		exit 1
	fi
else
	write_status "Skipping GRIDSS Annotations = ${gridss_unfiltered_vcf} exists"
fi


gripss_dir=$run_dir/gripss
gripss_somatic_vcf=$gripss_dir/${tumour_sample}.gripss.somatic.vcf.gz
gripss_somatic_filtered_vcf=$gripss_dir/${tumour_sample}.gripss.somatic.filtered.vcf.gz
if [[ ! -f $gripss_somatic_vcf ]] ; then
	write_status "Running GRIPSS"
      if [[ $tumour_only != true ]] ; then
            gripss_args="-reference $normal_sample $gripss_args"
      fi
	java -Xmx24G -cp ${gripss_jar} com.hartwig.hmftools.gripss.GripssApplicationKt \
		-ref_genome ${ref_genome} \
		-breakpoint_hotspot ${breakpoint_hotspot} \
		-breakend_pon ${breakend_pon} \
		-breakpoint_pon ${breakpoint_pon} \
		-input_vcf ${gridss_unfiltered_vcf} \
		-output_vcf ${gripss_somatic_vcf} \
		-tumor ${tumour_sample} \
		$gripss_args 2>&1 | tee $log_prefix.gripss.log
	if [[ ! -f $gripss_somatic_vcf ]] ; then
		write_status "Error creating $gripss_somatic_vcf. Aborting"
		exit 1
	fi
	java -Xmx24G -cp ${gripss_jar} com.hartwig.hmftools.gripss.GripssHardFilterApplicationKt \
		-input_vcf ${gripss_somatic_vcf} \
		-output_vcf ${gripss_somatic_filtered_vcf} \

	if [[ ! -f ${gripss_somatic_filtered_vcf} ]] ; then
		write_status "Error creating ${gripss_somatic_filtered_vcf} - aborting"
		exit 1
	fi
else
	write_status "Skipping GRIPSS - ${gripss_somatic_vcf} exists"
fi

mkdir -p $run_dir/amber
amber_vcf=$run_dir/amber/$tumour_sample.amber.baf.vcf.gz
if [[ ! -f ${amber_vcf} ]] ; then
	write_status "Running AMBER"
	if [[ $tumour_only != true ]] ; then
		amber_args="-reference $normal_sample -reference_bam $normal_bam $amber_args"
	else
		amber_args="-tumor_only $amber_args"
	fi
	java -Xmx10G $jvm_args \
		-jar $amber_jar \
		-threads $threads \
		-tumor $tumour_sample \
		-tumor_bam $tumour_bam \
		-loci $bafsnps \
		-ref_genome $ref_genome \
		-validation_stringency $validation_stringency \
		-output_dir $run_dir/amber \
		$amber_args 2>&1 | tee $log_prefix.amber.log
	if [[ ! -f ${amber_vcf} ]] ; then
		write_status "Error running AMBER - aborting"
		exit 1
	fi
else
	write_status "Skipping AMBER - ${amber_vcf} exists"
fi


mkdir -p $run_dir/cobalt
cobalt_file=$run_dir/cobalt/$tumour_sample.cobalt.ratio.pcf
if [[ ! -f ${cobalt_file} ]] ; then
	write_status "Running COBALT"
	if [[ $tumour_only != true ]] ; then
		cobalt_args="-reference $normal_sample -reference_bam $normal_bam $cobalt_args"
	else
		cobalt_args="-tumor_only $cobalt_args"
	fi
	java -Xmx10G $jvm_args \
		-cp ${cobalt_jar} com.hartwig.hmftools.cobalt.CountBamLinesApplication \
		-threads ${threads} \
		-tumor ${tumour_sample} \
		-tumor_bam ${tumour_bam} \
		-ref_genome $ref_genome \
		-output_dir ${run_dir}/cobalt \
		-gc_profile ${gcprofile} \
		-threads ${threads} \
		$cobalt_args \
		2>&1 | tee $log_prefix.cobalt.log
	if [[ ! -f ${cobalt_file} ]] ; then
		write_status "Error running COBALT - aborting"
		exit 1
	fi
else
	write_status "Skipping COBALT - ${cobalt_file} exists"
fi

mkdir -p $run_dir/purple

# circos requires /home/$LOGNAME to exist
if [[ -z "${LOGNAME:-}" ]] ; then
	export LOGNAME=$(whoami)
	mkdir -p /home/$LOGNAME
fi

purple_vcf=$run_dir/purple/$tumour_sample.purple.sv.vcf.gz
if [[ ! -f ${purple_vcf} ]] ; then
	write_status "Running PURPLE"
	if [[ -f "$snvvcf" ]] ; then
		purple_args="-somatic_vcf $snvvcf $purple_args"
	fi
	if [[ $tumour_only != true ]] ; then
		purple_args="-reference $normal_sample $purple_args"
	else
		purple_args="-reference DIPLOID -tumor_only $purple_args"
	fi
	java -Dorg.jooq.no-logo=true -Xmx10G $jvm_args \
		-jar ${purple_jar} \
		-output_dir $run_dir/purple \
		-tumor $tumour_sample \
		-amber $run_dir/amber \
		-cobalt $run_dir/cobalt \
		-gc_profile $gcprofile \
		-ref_genome $ref_genome \
		-structural_vcf ${gripss_somatic_filtered_vcf} \
		-sv_recovery_vcf ${gripss_somatic_vcf} \
		-driver_catalog -hotspots ${known_hotspots_vcf} \
		-driver_gene_panel ${driver_gene_panel} \
		-circos circos \
		-threads ${threads} \
		$purple_args \
		2>&1 | tee $log_prefix.purple.log
	if [[ ! -f ${purple_vcf} ]] ; then
		write_status "Error running PURPLE - aborting"
		exit 1
	fi
else
	write_status "Skipping PURPLE - ${purple_vcf} exists"
fi

mkdir -p $run_dir/linx

write_status "Running LINX"
java -Xmx8G -Xms4G -jar ${linx_jar} \
	-ref_genome ${ref_genome} \
	-ref_genome_version ${ref_genome_version} \
	-sample ${tumour_sample} \
	-purple_dir $run_dir/purple \
	-sv_vcf ${purple_vcf} \
	-output_dir $run_dir/linx \
	-fragile_site_file ${fragile_sites} \
	-line_element_file ${line_elements} \
	-replication_origins_file ${replication_origins} \
	-viral_hosts_file ${viral_hosts_csv} \
	-gene_transcripts_dir ${ensembl_data_dir} \
	-check_fusions \
	-known_fusion_file ${known_fusion_csv} \
	-driver_gene_panel ${driver_gene_panel} \
	-check_drivers \
	-write_vis_data \
	$linx_args \
	2>&1 | tee $log_prefix.linx.log
write_status "Generating LINE visualisations"
java -cp ${linx_jar} com.hartwig.hmftools.linx.visualiser.SvVisualiser \
	-sample ${tumour_sample} \
	-plot_out $run_dir/linx/plot/ \
	-data_out $run_dir/linx/circos/ \
	-vis_file_dir $run_dir/linx/ \
	-circos circos \
	-threads $threads

write_status "GRIDSS - Purple - Linx script complete"
trap - EXIT
exit 0 # success!

