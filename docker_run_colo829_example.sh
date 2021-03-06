#!/bin/bash

docker run --ulimit nofile=100000:100000 \
	-v d:/refdata/hg19:/refdata \
	-v d:/colo829:/data/ \
	gridss/gridss-purple-linx:latest \
	-n /data/COLO829R_dedup.realigned.bam \
	-t /data/COLO829T_dedup.realigned.bam \
	-v /data/gridss-purple-linx-latest.vcf \
	-s colo829 \
	--snvvcf /data/COLO829v003T.somatic_caller_post_processed.vcf \
	--jvmheap 15g
