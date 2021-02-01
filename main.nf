#!/usr/bin/env nextflow
/*
========================================================================================
                         mpozuelo/MGI_demux
========================================================================================
mpozuelo/MGI_demux Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/mpozuelo/MGI_demux
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info mpozueloHeader()
    log.info """

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run mpozuelo/MGI_demux --input '*.txt' -profile docker

    Mandatory arguments:
      --input [file]                Samplesheet with indexes and samples information
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity.
      --cluster_path                Cluster path to store data and to search for input data (Default: /datos/ngs/dato-activo)

    Optional:
      --single_index                In case the sequencing has been done single index although there were two indexes (i5 and i7). Only run with i7 and avoid the second index2 removal

    Demultiplexing parameters:
      --max_errors                  Maximum error rate accepted. For 8bp adapters we need 0.15 to allow 1bp of error. (Default: 0.15)
      --save_untrimmed              Saves untrimmed reads when demultiplexing (Default: FALSE)

    QC:
      --skipQC                      Skip all QC steps apart from MultiQC
      --skipFastQC                  Skip FastQC

    Other options
      --outdir                      The output directory where the results will be saved
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    """.stripIndent()
}


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


/*
 * SET UP CONFIGURATION VARIABLES
 */


 // Has the run name been specified by the user?
 //  this has the bonus effect of catching both -name and --name
 custom_runName = params.name
 if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
   custom_runName = workflow.runName
 }
 else{
   workflow.runName = params.user + " " + params.timestamp
   custom_runName = workflow.runName
 }


// Define regular variables so that they can be overwritten
max_errors = params.max_errors

// Validate inputs

if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Input samplesheet file not specified!" }

if (!params.outdir) {
  params.outdir = params.run
}

cluster_path = params.cluster_path


// Header log info
log.info mpozueloHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name'] = custom_runName ?: workflow.runName
summary['Input'] = params.input
summary['Single index'] = params.single_index
summary['Demultiplexing max error rate'] = max_errors
summary['Max Resources'] = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['User'] = workflow.userName

summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"


// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'mpozuelo-MGI_demux-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'mpozuelo/MGI_demux Workflow Summary'
    section_href: 'https://github.com/mpozuelo/MGI_demux'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}



/*
 * LOAD SAMPLESHEET and assign get the columns we will use for demultiplexing

 /*
  * LOAD SAMPLESHEET and assign get the columns we will use for demultiplexing
  It contains the following columns:
 1- Name given to the sample
 2- Index
 3- Index2
 4- Barcode
 5- Run ID
 6- Lane
 7- Sequencing date
 8- Protocol
 9- Platform
 10- Sample Source (place where the samples come from)
 11- Genome
 12- User
 13- Coverage
 */



Channel
  .from( ch_input )
  .splitCsv(header:false, sep:',')
  .map { it = ["${it[0]}", "${it[1]}", "${it[2]}", "${it[3]}", "${it[4]}", "${it[5]}", "${it[7]}", "${it[8]}", "${it[11]}",
  [file("${cluster_path}/02_rfastq/${it[8]}/${it[4]}/${it[5]}/${it[4]}_${it[5]}_read_1.fq.gz", checkIfExists: true),
  file("${cluster_path}/02_rfastq/${it[8]}/${it[4]}/${it[5]}/${it[4]}_${it[5]}_read_2.fq.gz", checkIfExists: true)]]}
  .set { ch_demux }


/*
 * STEP 1 - Demultiplex - Index1
 */
//Detect index in the end of read2

process demux_index {
  tag "$sample"
  label 'process_high'
  publishDir "${cluster_path}/03_intermediate/${platform}/${run_id}/${lane}/Index-removal/", mode: 'copy',
  saveAs: { filename ->
    filename.endsWith(".log") ? "logs/$filename" : filename
  }

  input:
  set val(sample), val(index), val(index2), val(barcode), val(run_id), val(lane), val(protocol), val(platform), val(user), path(reads) from ch_demux

  output:
  set val(sample), path("*.fq.gz"), val(index), val(index2), val(barcode), val(run_id), val(lane), val(protocol), val(platform), val(user) into ch_demux_index2
  path("*.{fq.gz,log}")

  script:
  //discard = params.save_untrimmed ? '' : '--discard-untrimmed'
  read1 = reads[0]
  read2 = reads[1]
  read1_index = "${sample}_${run_id}_${lane}_${index}_R1.fq.gz"
  read2_index = "${sample}_${run_id}_${lane}_${index}_R2.fq.gz"
  errors = index.length() > 6 ? "-e 0.15" : "-e 0.2"


  if (index == "NNNNNNNN") {
    """
    length=($(echo -e `zcat $read1 | head -2 | tail -1 | awk '{print length(\$0)}'`))
    cutadapt -l ${length} -o $read1_index $read1 -j 0 > "${sample}_${run_id}_${lane}_${index}_R1.log"
    cutadapt -l ${length} -o $read2_index $read2 -j 0 > "${sample}_${run_id}_${lane}_${index}_R2.log"
    """
  } else {
    """
    cutadapt \
    $errors \
    --no-indels \
    -a $sample=\"$index\$\" \
    -o $read2_index -p $read1_index \
    $read2 $read1 \
    -j 0 \
    --discard-untrimmed > ${sample}_${run_id}_${lane}_${index}.log
    """
  }

}



/*
 * STEP 2 - Demultiplex - Index2
 */
 //After removing index, remove index2 from the end (again) of read2

 if (!params.single_index) {
   process demux_index2 {
     tag "$sample"
     label 'process_high'
     publishDir "${cluster_path}/03_intermediate/${platform}/${run_id}/${lane}/Index2-removal/", mode: 'copy',
     saveAs: { filename ->
       filename.endsWith(".log") ? "logs/$filename" : filename
     }

     input:
     set val(sample), path(reads), val(index), val(index2), val(barcode), val(run_id), val(lane), val(protocol), val(platform), val(user) from ch_demux_index2

     output:
     set val(sample), path("*.fq.gz"), val(index), val(index2), val(barcode), val(run_id), val(lane), val(protocol), val(platform), val(user) into ch_demux_BC
     path("*.{fq.gz,log}")

     script:
     //discard = params.save_untrimmed ? '' : '--discard-untrimmed'
     read1 = reads[0]
     read2 = reads[1]
     read1_index2 = "${sample}_${run_id}_${lane}_${index}_${index2}_R1.fq.gz"
     read2_index2 = "${sample}_${run_id}_${lane}_${index}_${index2}_R2.fq.gz"
     errors = index2.length() > 6 ? "-e 0.15" : "-e 0.2"

     if (index2 == "NNNNNNNN") {
       if (index == "NNNNNNNN") {
         """
         mv $read1 $read1_index2
         mv $read2 $read2_index2
         """
       } else {
         """
         length=($(echo -e `zcat $read1 | head -2 | tail -1 | awk '{print length(\$0)}'`))
         cutadapt -l ${length} -o $read1_index2 $read1 -j 0 > "${sample}_${run_id}_${lane}_${index}_${index2}_R1.log"
         cutadapt -l ${length} -o $read2_index2 $read2 -j 0 > "${sample}_${run_id}_${lane}_${index}_${index2}_R2.log"
         """
       }
     } else {
       """
       cutadapt \
       $errors \
       --no-indels \
       -a $sample=\"$index2\$\" \
       -o $read2_index2 -p $read1_index2 \
       $read2 $read1 \
       -j 0 \
       --discard-untrimmed > "${sample}_${run_id}_${lane}_${index}_${index2}.log"
       """
     }
   }
 } else {
   ch_demux_index2
    .set { ch_demux_BC }
 }




/*
 * STEP 3 - Demultiplex - BC
 */
//Detect barcode sequence in the beginning or read1

process demux_BC {
  tag "$sample"
  label 'process_high'
  publishDir "${cluster_path}/04_pfastq/${platform}/${run_id}/${lane}/${user}/demux_fastq/", mode: 'copy',
  saveAs: { filename ->
    filename.endsWith(".log") ? "logs/$filename" : filename
  }

  input:
  set val(sample), path(reads), val(index), val(index2), val(barcode), val(run_id), val(lane), val(protocol), val(platform), val(user) from ch_demux_BC

  output:
  set val(sample), path("*.fq.gz"), val(run_id), val(lane), val(platform), val(user) into ch_fastqc
  set val(sample), path("*.fq.gz"), val(index), val(barcode), val(run_id), val(lane), val(protocol), val(platform), val(user) into ch_single_cell_header
  path("*.{fq.gz,log}")

  script:
  //discard = params.save_untrimmed ? '' : '--discard-untrimmed'
  read1 = reads[0]
  read2 = reads[1]
  read1_BC = "${sample}_${run_id}_${lane}_R1.fq.gz"
  read2_BC = "${sample}_${run_id}_${lane}_R2.fq.gz"
  errors = barcode.length() > 6 ? "-e 0.15" : "-e 0.2"

  if (barcode == "NNNNNNNN" | barcode == "NNNNNN") {
    """
    mv $read1 $read1_BC
    mv $read2 $read2_BC
    """
  } else {
    """
    cutadapt \
    $errors \
    --no-indels \
    -g $sample=\"^$barcode\" \
    -o $read1_BC -p $read2_BC \
    $read1 $read2 \
    -j 0 \
    --discard-untrimmed > "${sample}_${run_id}_${lane}.log"
    """
  }
}


process single_cell_fastq {
  tag "$sample"
  label 'process_high'
  publishDir "${cluster_path}/04_pfastq/${platform}/${run_id}/${lane}/${user}/single_cell_header/", mode: 'copy',

  when:
  protocol == "scRNAseq"

  input:
  set val(sample), path(reads), val(index), val(barcode), val(run_id), val(lane), val(protocol), val(platform), val(user) from ch_single_cell_header

  output:
  path("*.fastq.gz")

  script:
  fqheader1 = reads[0].minus(".fq.gz") + "_BC.fq.gz"
  fqheader2 = reads[1].minus(".fq.gz") + "_BC.fq.gz"

  // Re-write reads header to Illumina format, taking info from MGI headers
  // For this step, BC sequence is collected from header (BC was incuded in the header in previous step)

  """
  zcat ${reads[0]} | awk -v var="$index" '{if (NR%4 == 1){print \$1"_"var} else{print \$1}}' | gzip > $fqheader1 &
  zcat ${reads[1]} | awk -v var="$index" '{if (NR%4 == 1){print \$1"_"var} else{print \$1}}' | gzip > $fqheader2
  File_ID_new=$(echo "${sample}" | rev | cut -c 3- | rev)
  File_ID_number=$(echo "${sample}" | rev | cut -c 1 | rev)
  Lane_ID_number=$(echo "${lane}" | rev | cut -c 1 | rev)
  python convertHeaders.py -i $fqheader1 -o ${File_ID_new}_S1_L00${Lane_ID_number}_R1_00${File_ID_number}.fastq.gz &
  python convertHeaders.py -i $fqheader2 -o ${File_ID_new}_S1_L00${Lane_ID_number}_R2_00${File_ID_number}.fastq.gz
  """
}


/*
 * STEP 4 - FastQC
 */
 if (!params.skipQC || !params.skipFastQC) {

   process fastqc {
     tag "$sample"
     label 'process_medium'
     publishDir "${cluster_path}/04_pfastq/${platform}/${run_id}/${lane}/${user}/fastqc/${sample}", mode: 'copy',
     saveAs: { filename ->
       filename.endsWith(".zip") ? "zips/$filename" : filename
     }

     input:
     set val(sample), path(reads), val(run_id), val(lane), val(platform), val(user) from ch_fastqc

     output:
     set path("*_fastqc.{zip,html}"), val(run_id), val(lane), val(platform), val(user) into fastqc_results

     script:
     """
     fastqc --quiet --threads $task.cpus $reads
     """
   }
 } else {
   fastqc_results = Channel.empty()
 }



// Generate the fastq files for the samples

 /*
  * Parse software version numbers
  */
 process get_software_versions {
     publishDir "${params.outdir}/pipeline_info", mode: 'copy',
         saveAs: { filename ->
             if (filename.indexOf(".csv") > 0) filename
             else null
         }

     output:
     file 'software_versions_mqc.yaml' into software_versions_yaml
     file "software_versions.csv"

     script:
     """
     fastqc --version &> v_fastqc.txt
     cutadapt --version &> v_cutadapt.txt
     scrape_software_versions.py &> software_versions_mqc.yaml
     """
 }




/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[mpozuelo/MGI_demux] Successful: $workflow.runName"

    if (!workflow.success) {
      subject = "[mpozuelo/MGI_demux] FAILED: $workflow.runName"
    }



    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";


    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
        log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
        log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "${c_purple}[mpozuelo/MGI_demux]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[mpozuelo/MGI_demux]${c_red} Pipeline completed with errors${c_reset}"
    }

}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def mpozueloHeader() {
  // Log colors ANSI codes
  c_blue = params.monochrome_logs ? '' : "\033[0;34m";
  c_dim = params.monochrome_logs ? '' : "\033[2m";
  c_white = params.monochrome_logs ? '' : "\033[0;37m";
  c_reset = params.monochrome_logs ? '' : "\033[0m";


  return """    -${c_dim}--------------------------------------------------${c_reset}-
  ${c_blue}  __  __  __   __  ___         ${c_reset}
  ${c_blue}  | \\/ | |__| |  |  /  |  |     ${c_reset}
  ${c_blue}  |    | |    |__| /__ |__|         ${c_reset}
  ${c_white}  mpozuelo/MGI_demux v${workflow.manifest.version}${c_reset}
  -${c_dim}--------------------------------------------------${c_reset}-
  """.stripIndent()
}


def checkHostname() {
  def c_reset = params.monochrome_logs ? '' : "\033[0m"
  def c_white = params.monochrome_logs ? '' : "\033[0;37m"
  def c_red = params.monochrome_logs ? '' : "\033[1;91m"
  def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
  if (params.hostnames) {
    def hostname = "hostname".execute().text.trim()
    params.hostnames.each { prof, hnames ->
      hnames.each { hname ->
        if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
          log.error "====================================================\n" +
          "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
          "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
          "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
          "============================================================"
        }
      }
    }
  }
}
