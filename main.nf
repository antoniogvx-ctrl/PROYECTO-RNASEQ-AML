nextflow.enable.dsl=2

// Parámetros iniciales
params.samplesheet = 'samplesheet.csv'
params.outdir      = 'results'

//Funciones de los procesos 

// Descarga de los FASTQ
process DOWNLOAD_SRA {
    container 'quay.io/biocontainers/sra-tools:3.0.3--h87f3376_0'
    publishDir "${params.outdir}/fastq", mode: 'copy'
    cpus 4

    input:
    tuple val(sample), val(group), val(srr)

    output:
    tuple val(sample), val(group), path("${sample}_R*.fastq.gz"), emit: reads

    script:
    """
    fasterq-dump --split-files --threads ${task.cpus} ${srr}
    gzip *.fastq
    mv ${srr}_1.fastq.gz ${sample}_R1.fastq.gz
    mv ${srr}_2.fastq.gz ${sample}_R2.fastq.gz
    """
}

// Análisis de calidad de los FASTQC
process FASTQC {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    cpus 2

    input:
    tuple val(sample), val(group), path(reads)

    output:
    path("*_fastqc.{zip,html}"), emit: qc_results

    script:
    """
    fastqc -q -t ${task.cpus} ${reads}
    """
}

// Alineamiento con el genoma mediante HISAT2
process HISAT2_ALIGN {
    container 'quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2388ff67fc407dad75774291ca5038f40cac4be0-0'
    publishDir "${params.outdir}/bam", mode: 'copy'
    cpus 4

    input:
    tuple val(sample), val(group), path(reads)
    path index_dir

    output:
    path("${sample}.bam"), emit: bam

    script:
    """
    hisat2 -p ${task.cpus} \\
           -x ${index_dir}/genome \\
           -1 ${reads[0]} \\
           -2 ${reads[1]} | \\
    samtools view -bS - | samtools sort -@ ${task.cpus} -o ${sample}.bam
    """
}

// Conteo de genes.Toma todos los BAMs a la vez para generar una matriz única
process FEATURECOUNTS {
    container 'quay.io/biocontainers/subread:2.0.1--hed695b0_0'
    publishDir "${params.outdir}/counts", mode: 'copy'
    cpus 4

    input:
    path bams
    path gtf

    output:
    path "gene_counts_matrix.txt", emit: counts
    path "gene_counts_matrix.txt.summary", emit: summary

    script:
    """
    featureCounts -p -T ${task.cpus} -a ${gtf} -o gene_counts_matrix.txt ${bams}
    """
}

// Análisis de calidad múltiple de los alineamientos 
process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}


// Flujo de trabajo 
workflow {
    // Definimos dónde están las referencias 
    ch_hisat2_index = file("ref/hg38")
    ch_gtf          = file("ref/Homo_sapiens.GRCh38.110.gtf")

    // Leemos el .csv con metadatos
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row -> tuple(row.sample, row.group, row.srr) }
        .set { ch_samples }

    // Descarga de SRA y análisis de calidad FASTQC
    DOWNLOAD_SRA(ch_samples)
    FASTQC(DOWNLOAD_SRA.out.reads)
    
    // Alineamiento HISAT2
    HISAT2_ALIGN(DOWNLOAD_SRA.out.reads, ch_hisat2_index)
    
    // Agrupamos todos los BAM terminados y los pasamos a FeatureCounts. (.collect() para esperar a que terminen todos los alineamientos primero)
    ch_all_bams = HISAT2_ALIGN.out.bam.collect()
    FEATURECOUNTS(ch_all_bams, ch_gtf)

    // Agrupamos los resultados de FastQC y Featurecounts para hacer el multiqc
    ch_multiqc_input = FASTQC.out.qc_results.mix(FEATURECOUNTS.out.summary).collect()
    MULTIQC(ch_multiqc_input)
}
