#! /usr/bin/env nextflow

log.info """
	bpnet Pipeline
	===================================
	samplesheet: ${params.samplesheet}
	results_dir: ${params.results_dir}
	"""
	.stripIndent()

/* input files:	
 * samplesheet
 * merged, sorted and indexed bam files for each ChIP sample
 * Optional negative control (should be IgG or similar, not input)
 * peaks in narrowPeak format with summits
 * genome fasta, chromosome sizes file, and list of chromosomes to keep
 * chromosome splits
 * json file for model architecture
 * additional training parameters
 */

/*
 * convert bam files to strand-specific bigwigs
 */

process stranded_bigwig {
	tag "$meta.sample"
	publishDir "${params.results_dir}/bigwigs/", mode: 'copy'
	container = "vivekramalingam/tf-atlas:gcp-modeling_v2.1.0-rc.1"
	
	input:
	tuple val(meta), path(bam), path(bam_index), path(control_bam), path(control_index)
	path(chrom_sizes)
	
	output:
	tuple val(meta), path("${meta.sample}_plus.bw"), path("${meta.sample}_minus.bw"), path("${meta.sample}_control_plus.bw"), path("${meta.sample}_control_minus.bw")
	
	script:
    """
    # get coverage of 5’ positions of the plus strand
	samtools view -b $bam \$(cut -f 1 ${chrom_sizes}) | \
		bedtools genomecov -5 -bg -strand + -ibam stdin | \
		sort -k1,1 -k2,2n > plus.bedGraph

	# get coverage of 5’ positions of the minus strand
	samtools view -b $bam \$(cut -f 1 ${chrom_sizes}) | \
        bedtools genomecov -5 -bg -strand - -ibam stdin | \
        sort -k1,1 -k2,2n > minus.bedGraph

	# Convert bedGraph files to bigWig files
	bedGraphToBigWig plus.bedGraph ${chrom_sizes} ${meta.sample}_plus.bw
	bedGraphToBigWig minus.bedGraph ${chrom_sizes} ${meta.sample}_minus.bw
	
	# get coverage of 5’ positions of the control plus strand
	samtools view -b $control_bam \$(cut -f 1 ${chrom_sizes}) | \
    	bedtools genomecov -5 -bg -strand + -ibam stdin | \
    	sort -k1,1 -k2,2n > control_plus.bedGraph

	# get coverage of 5' positions of the control minus strand
	samtools view -b $control_bam \$(cut -f 1 ${chrom_sizes}) | \
			bedtools genomecov -5 -bg -strand - -ibam stdin | \
			sort -k1,1 -k2,2n > control_minus.bedGraph
	
	# Convert bedGraph files to bigWig files
	bedGraphToBigWig control_plus.bedGraph ${chrom_sizes} ${meta.sample}_control_plus.bw
	bedGraphToBigWig control_minus.bedGraph ${chrom_sizes} ${meta.sample}_control_minus.bw
    """
    
    stub:
    """
    touch ${meta.sample}_plus.bw
    touch ${meta.sample}_minus.bw
    touch ${meta.sample}_control_plus.bw
    touch ${meta.sample}_control_minus.bw
    """
}


/*
 * remove outlier peaks? Is this necessary?
 */

process remove_outliers {
	tag "$meta.sample"
	publishDir "${params.results_dir}/filtered_peaks/", mode: 'copy'
	container = "vivekramalingam/tf-atlas:gcp-modeling_v2.1.0-rc.1"
	
	input:
	tuple val(meta), path(peaks), path(plus_bw), path(minus_bw), path(control_plus_bw), path(control_minus_bw)
	path(chrom_sizes)
	path(keep_chroms)
	path(blacklist)
	
	output:
	tuple val(meta), path("${meta.sample}_filtered_peaks.bed")
	
	script:
    """
    touch input.json
	echo "{" >> input.json
	echo "\\"0\\": {" >> input.json
	echo "\\"signal\\": {" >> input.json
	echo "\\"source\\": [\\"${plus_bw}\\"," >> input.json
	echo "\\"${minus_bw}\\"]" >> input.json
	echo "}," >> input.json
	echo "\\"loci\\": {" >> input.json
	echo "\\"source\\": [\\"${peaks}\\"]" >> input.json
	echo "}," >> input.json
	echo "\\"bias\\": {" >> input.json
	echo "\\"source\\": [\\"${control_plus_bw}\\"," >> input.json
	echo "\\"${control_minus_bw}\\"]," >> input.json
	echo "\\"smoothing\\": [null, null]" >> input.json
	echo "}" >> input.json
	echo "}" >> input.json
	echo "}" >> input.json
    
    bpnet-outliers \
    --input-data input.json  \
    --quantile 0.99 \
    --quantile-value-scale-factor 1.2 \
    --task 0 \
    --chrom-sizes ${chrom_sizes} \
    --chroms \$(paste -s -d ' ' ${keep_chroms}) \
    --sequence-len 1000 \
    --blacklist $blacklist \
    --global-sample-weight 1.0 \
    --output-bed ${meta.sample}_filtered_peaks.bed
    """
    stub:
    """
    touch input.json
	echo "{" >> input.json
	echo "\\"0\\": {" >> input.json
	echo "\\"signal\\": {" >> input.json
	echo "\\"source\\": [\\"${plus_bw}\\"," >> input.json
	echo "\\"${minus_bw}\\"]" >> input.json
	echo "}," >> input.json
	echo "\\"loci\\": {" >> input.json
	echo "\\"source\\": [\\"${peaks}\\"]" >> input.json
	echo "}," >> input.json
	echo "\\"bias\\": {" >> input.json
	echo "\\"source\\": [\\"${control_plus_bw}\\"," >> input.json
	echo "\\"${control_minus_bw}\\"]," >> input.json
	echo "\\"smoothing\\": [null, null]" >> input.json
	echo "}" >> input.json
	echo "}" >> input.json
	echo "}" >> input.json
    
    touch ${meta.sample}_filtered_peaks.bed
    """
}

/*
 * Generate GC matched background regions
 */

process gc_matched_negatives {
	tag "$meta.sample"
	publishDir "${params.results_dir}/modeling/background_regions/", mode: 'copy'
	container = "vivekramalingam/tf-atlas:gcp-modeling_v2.1.0-rc.1"
	
	input:
	tuple val(meta), path(peaks)
	path(fasta)
	path(chrom_sizes)
	
	output:
	tuple val(meta), path("${meta.sample}_gc_negatives.bed")
	
	script:
    """
	bpnet-gc-reference \
			--ref_fasta $fasta \
			--chrom_sizes $chrom_sizes \
			--out_prefix genomewide_gc.bed \
			--inputlen 2114 \
			--stride 1000
        
    
	bpnet-gc-background \
			--peaks_bed $peaks \
			--out_dir . \
			--ref_gc_bed genomewide_gc.bed \
			--out_prefix ${meta.sample}_gc_negatives.bed \
			--flank_size 1057 \
			--neg_to_pos_ratio_train 4       
	"""
    
    stub:
    """
    touch ${meta.sample}_gc_negatives.bed
    """	
}


/*
 * model training
 */

process train_bpnet {
	tag "$meta.sample"
	publishDir "${params.results_dir}/bpnet_models/${meta.sample}/", mode: 'copy'
	container = "vivekramalingam/tf-atlas:gcp-modeling_v2.1.0-rc.1"
	
	input:
	tuple val(meta), path(peaks), path(plus_bw), path(minus_bw), path(control_plus_bw), path(control_minus_bw) //, path(gc_background)
	path(fasta)
	path(chrom_sizes)
	path(keep_chroms)
	path(chrom_splits)
	path(bpnet_params)
	
	output:
	tuple val(meta), path("model_split000"), path("${meta.sample}_input.json")
	
	script:
	// NOTE: how to handle input json files?
    """
    touch ${meta.sample}_input.json
    echo "{" >> ${meta.sample}_input.json
	echo "    \\"0\\": {" >> ${meta.sample}_input.json
	echo "        \\"signal\\": {" >> ${meta.sample}_input.json
	echo "            \\"source\\": [\\"${plus_bw}\\", " >> ${meta.sample}_input.json
	echo "                     \\"${minus_bw}\\"]" >> ${meta.sample}_input.json
	echo "        }," >> ${meta.sample}_input.json
	echo "        \\"loci\\": {" >> ${meta.sample}_input.json
	echo "            \\"source\\": [\\"${peaks}\\"]" >> ${meta.sample}_input.json
	echo "        }," >> ${meta.sample}_input.json
	echo "        \\"bias\\": {" >> ${meta.sample}_input.json
	echo "            \\"source\\": [\\"${control_plus_bw}\\"," >> ${meta.sample}_input.json
	echo "                     \\"${control_minus_bw}\\"]," >> ${meta.sample}_input.json
	echo "            \\"smoothing\\": [null, null]" >> ${meta.sample}_input.json
	echo "        }" >> ${meta.sample}_input.json
	echo "    }" >> ${meta.sample}_input.json
	echo "}" >> ${meta.sample}_input.json
    
	bpnet-train \
        --input-data ${meta.sample}_input.json \
        --output-dir . \
        --reference-genome $fasta \
        --chroms \$(paste -s -d ' ' ${keep_chroms}) \
        --chrom-sizes $chrom_sizes \
        --splits $chrom_splits \
        --model-arch-name BPNet \
        --model-arch-params-json $bpnet_params \
        --sequence-generator-name BPNet \
        --model-output-filename model \
        --input-seq-len 2114 \
        --output-len 1000 \
        --shuffle \
        --threads 1 \
        --epochs 100 \
	   --batch-size 64 \
		--reverse-complement-augmentation \
		--early-stopping-patience 10 \
		--reduce-lr-on-plateau-patience 5 \
        --learning-rate 0.001
	"""
	
	stub:
	"""
     touch ${meta.sample}_input.json
    echo "{" >> ${meta.sample}_input.json
	echo "    \\"0\\": {" >> ${meta.sample}_input.json
	echo "        \\"signal\\": {" >> ${meta.sample}_input.json
	echo "            \\"source\\": [\\"${plus_bw}\\", " >> ${meta.sample}_input.json
	echo "                     \\"${minus_bw}\\"]" >> ${meta.sample}_input.json
	echo "        }," >> ${meta.sample}_input.json
	echo "        \\"loci\\": {" >> ${meta.sample}_input.json
	echo "            \\"source\\": [\\"${peaks}\\"]" >> ${meta.sample}_input.json
	echo "        }," >> ${meta.sample}_input.json
	echo "        \\"bias\\": {" >> ${meta.sample}_input.json
	echo "            \\"source\\": [\\"${control_plus_bw}\\"," >> ${meta.sample}_input.json
	echo "                     \\"${control_minus_bw}\\"]," >> ${meta.sample}_input.json
	echo "            \\"smoothing\\": [null, null]" >> ${meta.sample}_input.json
	echo "        }" >> ${meta.sample}_input.json
	echo "    }" >> ${meta.sample}_input.json
	echo "}" >> ${meta.sample}_input.json

	mkdir model_split000
	touch model_split000/saved_model.pb
	"""
}

/*
 * generate predicted signal bigwigs on test regions
 */

process predicted_bw_test {
	tag "$meta.sample"
	publishDir "${params.results_dir}/predicted_test/", mode: 'copy'
	container = "vivekramalingam/tf-atlas:gcp-modeling_v2.1.0-rc.1"
	
	input:
	tuple val(meta), path(bpnet_model), path(input_json), path(peaks), path(plus_bw), path(minus_bw), path(control_plus_bw), path(control_minus_bw)
	path(fasta)
	path(chrom_sizes)
	val(test_chroms)
	
	output:
	tuple val(meta), path("${meta.sample}_pred")
	
	script:
	def test_chroms_str = test_chroms.join(" ")
    """
    mkdir ${meta.sample}_pred
    
    
	bpnet-predict \
        --model $bpnet_model\
        --chrom-sizes $chrom_sizes \
        --chroms $test_chroms_str \
        --test-indices-file None \
        --reference-genome $fasta \
        --output-dir ${meta.sample}_pred \
        --input-data $input_json \
        --sequence-generator-name BPNet \
        --input-seq-len 2114 \
        --output-len 1000 \
        --output-window-size 1000 \
        --batch-size 64 \
        --reverse-complement-average \
        --threads 1 \
        --generate-predicted-profile-bigWigs
	"""
	stub:
	def test_chroms_str = test_chroms.join(" ")
	"""
	ls ${bpnet_model}/saved_model.pb
	echo $test_chroms_str
	mkdir ${meta.sample}_pred
	
	"""
}

/*
 * generate predicted signal bigwigs on all regions
 */

process predicted_bw_all {
	tag "$meta.sample"
	publishDir "${params.results_dir}/predicted_all/", mode: 'copy'
	container = "vivekramalingam/tf-atlas:gcp-modeling_v2.1.0-rc.1"
	
	input:
	tuple val(meta), path(bpnet_model), path(input_json), path(peaks), path(plus_bw), path(minus_bw), path(control_plus_bw), path(control_minus_bw)
	path(fasta)
	path(chrom_sizes)
	val(all_chroms)
	
	output:
	tuple val(meta), path("${meta.sample}_pred")
	
	script:
	def all_chroms_str = all_chroms.join(" ")
    """
	 mkdir ${meta.sample}_pred
	
	bpnet-predict \
        --model $bpnet_model\
        --chrom-sizes $chrom_sizes \
        --chroms $all_chroms_str \
        --test-indices-file None \
        --reference-genome $fasta \
        --output-dir ${meta.sample}_pred \
        --input-data $input_json \
        --sequence-generator-name BPNet \
        --input-seq-len 2114 \
        --output-len 1000 \
        --output-window-size 1000 \
        --batch-size 64 \
        --reverse-complement-average \
        --threads 1 \
        --generate-predicted-profile-bigWigs
	"""
	stub:
	def all_chroms_str = all_chroms.join(" ")
	"""
	ls ${bpnet_model}/saved_model.pb
	echo $all_chroms_str
	mkdir ${meta.sample}_pred
	
	"""

}

/*
 * compute importance scores
 */

process compute_importance {
	tag "$meta.sample"
	publishDir "${params.results_dir}/contribution_all/", mode: 'copy'
	container = "vivekramalingam/tf-atlas:gcp-modeling_v2.1.0-rc.1"
	
	input:
	tuple val(meta), path(bpnet_model), path(input_json), path(peaks), path(plus_bw), path(minus_bw), path(control_plus_bw), path(control_minus_bw)
	path(fasta)
	path(chrom_sizes)
	val(all_chroms)
	
	output:
	tuple val(meta), path("${meta.sample}_shap")
	
	script:
	def all_chroms_str = all_chroms.join(" ")
    """
	 mkdir ${meta.sample}_shap
	
	bpnet-shap \
        --reference-genome $fasta \
        --model $bpnet_model  \
        --bed-file $peaks \
        --chroms $all_chroms_str \
        --output-dir ${meta.sample}_shap \
        --input-seq-len 2114 \
        --control-len 1000 \
        --task-id 0 \
        --input-data $input_json
	"""
	
	stub:
	def all_chroms_str = all_chroms.join(" ")
	"""
	ls ${bpnet_model}/saved_model.pb
	echo $all_chroms_str
	mkdir ${meta.sample}_shap
	
	"""

}

/*
 * motif discovery
 */

/*
 * identification of motif instances
 */


/*
 * Run workflow
 */
 
workflow {
	
	bam_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)
	| splitCsv( header:true )
    | map { row ->
        bam_meta = row.subMap('sample')
        [
        	bam_meta, 
        	file(row.reads, checkIfExists: true),
            file(row.read_index, checkIfExists: true),
            file(row.control_reads, checkIfExists: true),
            file(row.control_index, checkIfExists: true)]
    }
	
	bw_ch = stranded_bigwig(
	bam_ch,
	file("${launchDir}/${params.chrom_sizes}", checkIfExists: true)
	)
	
	peak_bw_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)
	| splitCsv( header:true )
    | map { row ->
        bam_meta = row.subMap('sample')
        [
        	bam_meta, 
        	file(row.peaks, checkIfExists: true)]
    }
    | combine(bw_ch, by: 0)
	
	filtered_peak_ch = remove_outliers(
	peak_bw_ch,
	file("${launchDir}/${params.chrom_sizes}", checkIfExists: true),
	file("${launchDir}/${params.keep_chroms}", checkIfExists: true),
	file("${launchDir}/${params.blacklist}", checkIfExists: true)
	)
	
// 	gc_background_ch = gc_matched_negatives(
// 	filtered_peak_ch,
// 	file("${launchDir}/${params.fasta}", checkIfExists: true),
// 	file("${launchDir}/${params.chrom_sizes}", checkIfExists: true)
// 	)
	
	train_ch = filtered_peak_ch.combine(bw_ch, by: 0)
// 	| combine(gc_background_ch, by: 0)
	
	model_ch = train_bpnet(
	train_ch,
	file("${launchDir}/${params.fasta}", checkIfExists: true),
	file("${launchDir}/${params.chrom_sizes}", checkIfExists: true),
	file("${launchDir}/${params.keep_chroms}", checkIfExists: true),
	file("${launchDir}/${params.chrom_splits}", checkIfExists: true),
	file(params.bpnet_params, checkIfExists: true)	
	)
	| combine(train_ch, by: 0)
// 	| view
	
	test_chroms_ch = Channel.fromPath(params.chrom_splits, checkIfExists: true)
	| splitJson(path: '0.test')
	| collect
// 	| view
	
	val_chroms_ch = Channel.fromPath(params.chrom_splits, checkIfExists: true)
	| splitJson(path: '0.val')
	| collect
// 	| view
	
	train_chroms_ch = Channel.fromPath(params.chrom_splits, checkIfExists: true)
	| splitJson(path: '0.train')
	| collect
// 	| view
	
	all_chroms_ch = Channel.fromPath(params.keep_chroms, checkIfExists: true)
	| splitCsv( header:false)
	| collect
// 	| view
	
	predicted_bw_test(
	model_ch,
	file("${launchDir}/${params.fasta}", checkIfExists: true),
	file("${launchDir}/${params.chrom_sizes}", checkIfExists: true),
	test_chroms_ch
	)
	
	predicted_bw_all(
	model_ch,
	file("${launchDir}/${params.fasta}", checkIfExists: true),
	file("${launchDir}/${params.chrom_sizes}", checkIfExists: true),
	all_chroms_ch
	)

	compute_importance(
	model_ch,
	file("${launchDir}/${params.fasta}", checkIfExists: true),
	file("${launchDir}/${params.chrom_sizes}", checkIfExists: true),
	all_chroms_ch
	)

		
}

