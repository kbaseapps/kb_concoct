import errno
import json
import os
import subprocess
import sys
import time
import uuid
import zipfile
import copy

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from installed_clients.ReadsUtilsClient import ReadsUtils
# from installed_clients.KBParallelClient import KBParallel


from .fasta import HeaderRenaming

from random import seed
from random import randint
# seed random number generator
seed(1)


# for future expansion
# from kb_concoct.BinningUtilities import BinningUtil as bu


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class ConcoctUtil:
    CONCOCT_BASE_PATH = '/kb/deployment/bin/CONCOCT'
    BINNER_RESULT_DIRECTORY = 'concoct_output_dir'
    BINNER_BIN_RESULT_DIR = 'final_bins'
    MAPPING_THREADS = 16
    BBMAP_MEM = '30g'
    MAX_NODES = 6  # for KBParallels

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.shock_url = config['shock-url']
        self.ws_url = config['workspace-url']
        self.dfu = DataFileUtil(self.callback_url)
        self.ru = ReadsUtils(self.callback_url)
        self.au = AssemblyUtil(self.callback_url)
        self.mgu = MetagenomeUtils(self.callback_url)
        self.fhr = HeaderRenaming()
        # self.parallel_runner = KBParallel(self.callback_url)

    # def set_up_parallel_tasks(self, task_params):

    #     # lets create the tasks that will be run in parallel
    #     tasks = []
    #     reads_list = task_params['reads_list']

    #     for fastq in reads_list:
    #         print("ADDING FASTQ FILE".format(fastq))
    #         task_params = copy.deepcopy(task_params)
    #         task_params['fastq'] = fastq
    #         module_input = 'kb_concoct'
    #         function_for_parallelizing = 'generate_alignment_bams'
    #         tasks.append( {'module_name': module_input,
    #                         'function_name': function_for_parallelizing,
    #                         'version': 'dev',
    #                         'parameters': { task_params, assembly_clean }
    #                       } )

    #     # calculate how many nodes we need (max=5)
    #     num_nsjw_nodes_required=0
    #     num_fastq = len(reads_list)
    #     if num_fastq < self.MAX_NODES:
    #         num_nsjw_nodes_required = num_fastq
    #     else:
    #         num_nsjw_nodes_required = self.MAX_NODES
    #         num_nsjw_nodes_required -= 1

    #     # kbparallel requires these params
    #     batch_run_params = {'tasks': tasks,
    #                         'runner': 'parallel',
    #                         'concurrent_local_tasks': 1,
    #                         'concurrent_njsw_tasks': 0, #num_nsjw_nodes_required,
    #                         'max_retries': 2}

    #     # We'll run the alignments in parallel here
    #     log('--->\nlaunching parallel jobs\n')
    #     kbparallel_results = self.parallel_runner.run_batch(batch_run_params)

    #     # print('KBPARALLLEL RESULTS')
    #     # pprint(kbparallel_results)
    #     # sys.exit()

    #     # ======================== #
    #     # back on the master node.
    #     # write depth.txt files to json files before trying to
    #     # merge them into one file for debugging reasons
    #     # -------------------------------------------------
    #     # check for errors returned by njsw remote jobs
    #     for fd in kbparallel_results['results']:
    #         if fd['result_package']['error'] is not None:
    #             log('kb_parallel failed to complete without throwing an error on at least one of the kbparallel-njsw nodes.')
    #             sys.exit(1)

    def _validate_run_concoct_params(self, task_params):
        """
        _validate_run_concoct_params:
                validates params passed to run_concoct method
        """
        log('Start validating run_concoct params')

        # check for required parameters
        for p in ['assembly_ref', 'binned_contig_name', 'workspace_name', 'reads_list', 'read_mapping_tool']:
            if p not in task_params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def _run_command(self, command):
        """
        _run_command: run command and print result
        """
        os.chdir(self.scratch)
        log('Start executing command:\n{}'.format(command))
        log('Command is running from:\n{}'.format(self.scratch))
        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        output, stderr = pipe.communicate()
        exitCode = pipe.returncode

        if (exitCode == 0):
            log('Executed command:\n{}\n'.format(command) +
                'Exit Code: {}\n'.format(exitCode))
        else:
            error_msg = 'Error running command:\n{}\n'.format(command)
            error_msg += 'Exit Code: {}\nOutput:\n{}\nStderr:\n{}'.format(exitCode, output, stderr)
            raise ValueError(error_msg)
            sys.exit(1)
        return (output, stderr)

    # this function has been customized to return read_type variable (interleaved vs single-end library)
    def stage_reads_list_file(self, reads_list):
        """
        stage_reads_list_file: download fastq file associated to reads to scratch area
                          and return result_file_path
        """

        log('Processing reads object list: {}'.format(reads_list))

        result_file_path = []
        read_type = []

        # getting from workspace and writing to scratch. The 'reads' dictionary now has file paths to scratch.
        reads = self.ru.download_reads({'read_libraries': reads_list, 'interleaved': None})['files']

        # reads_list is the list of file paths on workspace? (i.e. 12804/1/1).
        # "reads" is the hash of hashes where key is "12804/1/1" or in this case, read_obj and
        # "files" is the secondary key. The tertiary keys are "fwd" and "rev", as well as others.
        for read_obj in reads_list:
            files = reads[read_obj]['files']    # 'files' is dictionary where 'fwd' is key of file path on scratch.
            result_file_path.append(files['fwd'])
            read_type.append(files['type'])
            if 'rev' in files and files['rev'] is not None:
                result_file_path.append(files['rev'])

        return result_file_path, read_type

    def _get_contig_file(self, assembly_ref):
        """
        _get_contig_file: get contig file from GenomeAssembly object
        """
        contig_file = self.au.get_assembly_as_fasta({'ref': assembly_ref}).get('path')

        sys.stdout.flush()
        contig_file = self.dfu.unpack_file({'file_path': contig_file})['file_path']

        return contig_file

    def retrieve_and_clean_assembly(self, task_params):
        if os.path.exists(task_params['contig_file_path']):
            assembly = task_params['contig_file_path']
            print("FOUND ASSEMBLY ON LOCAL SCRATCH")
        else:
            # we are on njsw so lets copy it over to scratch
            assembly = self._get_contig_file(task_params['assembly_ref'])

        # remove spaces from fasta headers because that breaks bedtools
        assembly_clean = os.path.abspath(assembly).split('.fa')[0] + "_clean.fa"

        # file to dump the mapping from new headers to old
        header_mapping_filepath = os.path.join(self.scratch, 'header_mapping_' + str(uuid.uuid4()))

        self.fhr.rename_fasta_headers(assembly, fasta_filepath_new=assembly_clean, cache_filepath=header_mapping_filepath)

        """
        command = '/bin/bash reformat.sh in={} out={} addunderscore overwrite=true'.format(assembly, assembly_clean)

        log('running reformat command: {}'.format(command))
        out, err = self._run_command(command)
        """
        return assembly_clean


    def generate_stats_for_genome_bins(self, task_params, genome_bin_fna_file, bbstats_output_file):
        """
        generate_command: bbtools stats.sh command
        """
        log("running generate_stats_for_genome_bins on {}".format(genome_bin_fna_file))
        genome_bin_fna_file = os.path.join(self.scratch, self.BINNER_RESULT_DIRECTORY, genome_bin_fna_file)
        command = '/bin/bash stats.sh in={} format=3 > {}'.format(genome_bin_fna_file, bbstats_output_file)
        self._run_command(command)
        bbstats_output = open(bbstats_output_file, 'r').readlines()[1]
        n_scaffolds = bbstats_output.split('\t')[0]
        n_contigs = bbstats_output.split('\t')[1]
        scaf_bp = bbstats_output.split('\t')[2]
        contig_bp = bbstats_output.split('\t')[3]
        gap_pct = bbstats_output.split('\t')[4]
        scaf_N50 = bbstats_output.split('\t')[5]
        scaf_L50 = bbstats_output.split('\t')[6]
        ctg_N50 = bbstats_output.split('\t')[7]
        ctg_L50 = bbstats_output.split('\t')[8]
        scaf_N90 = bbstats_output.split('\t')[9]
        scaf_L90 = bbstats_output.split('\t')[10]
        ctg_N90 = bbstats_output.split('\t')[11]
        ctg_L90 = bbstats_output.split('\t')[12]
        scaf_max = bbstats_output.split('\t')[13]
        ctg_max = bbstats_output.split('\t')[14]
        scaf_n_gt50K = bbstats_output.split('\t')[15]
        scaf_pct_gt50K = bbstats_output.split('\t')[16]
        gc_avg = float(bbstats_output.split('\t')[17]) * 100  # need to figure out if correct
        gc_std = float(bbstats_output.split('\t')[18]) * 100  # need to figure out if correct

        log('Generated generate_stats_for_genome_bins command: {}'.format(command))

        return {'n_scaffolds': n_scaffolds,
                'n_contigs': n_contigs,
                'scaf_bp': scaf_bp,
                'contig_bp': contig_bp,
                'gap_pct': gap_pct,
                'scaf_N50': scaf_N50,
                'scaf_L50': scaf_L50,
                'ctg_N50': ctg_N50,
                'ctg_L50': ctg_L50,
                'scaf_N90': scaf_N90,
                'scaf_L90': scaf_L90,
                'ctg_N90': ctg_N90,
                'ctg_L90': ctg_L90,
                'scaf_max': scaf_max,
                'ctg_max': ctg_max,
                'scaf_n_gt50K': scaf_n_gt50K,
                'scaf_pct_gt50K': scaf_pct_gt50K,
                'gc_avg': gc_avg,
                'gc_std': gc_std
                }

    def deinterlace_raw_reads(self, fastq):
        fastq_forward = fastq.split('.fastq')[0] + "_forward.fastq"
        fastq_reverse = fastq.split('.fastq')[0] + "_reverse.fastq"
        command = 'deinterleave_fastq.sh < {} {} {}'.format(fastq, fastq_forward, fastq_reverse)
        try:
            self._run_command(command)
        except:
            raise Exception("Cannot deinterlace fastq file!")
        return (fastq_forward, fastq_reverse)

    def run_read_mapping_interleaved_pairs_mode(self, task_params, assembly_clean, fastq, sam):
        read_mapping_tool = task_params['read_mapping_tool']
        log("running {} mapping in interleaved mode.".format(read_mapping_tool))
        random_seed_int = randint(0, 999999999)
        log("randomly selected seed (integer) used for read mapping is: {}".format(random_seed_int))
        if task_params['read_mapping_tool'] == 'bbmap_fast':
            log("Warning: bbmap does not support setting random seeds, so results are not reproducible.")
            command = 'bbmap.sh -Xmx{} '.format(self.BBMAP_MEM)
            command += 'threads={} '.format(self.MAPPING_THREADS)
            command += 'ref={} '.format(assembly_clean)
            command += 'in={} '.format(fastq)
            command += 'out={} '.format(sam)
            command += 'fast interleaved=true mappedonly nodisk overwrite'
        elif task_params['read_mapping_tool'] == 'bbmap_default':
            log("Warning: bbmap does not support setting random seeds, so results are not reproducible.")
            command = 'bbmap.sh -Xmx{} '.format(self.BBMAP_MEM)
            command += 'threads={} '.format(self.MAPPING_THREADS)
            command += 'ref={} '.format(assembly_clean)
            command += 'in={} '.format(fastq)
            command += 'out={} '.format(sam)
            command += 'interleaved=true mappedonly nodisk overwrite'
        elif task_params['read_mapping_tool'] == 'bbmap_very_sensitive':
            log("Warning: bbmap does not support setting random seeds, so results are not reproducible.")
            command = 'bbmap.sh -Xmx{} '.format(self.BBMAP_MEM)
            command += 'threads={} '.format(self.MAPPING_THREADS)
            command += 'ref={} '.format(assembly_clean)
            command += 'in={} '.format(fastq)
            command += 'out={} '.format(sam)
            command += 'vslow=true '
            command += 'interleaved=true mappedonly nodisk overwrite'
        elif task_params['read_mapping_tool'] == 'bowtie2_default':
            (fastq_forward, fastq_reverse) = self.deinterlace_raw_reads(fastq)
            bt2index = os.path.basename(assembly_clean) + '.bt2'
            command = 'bowtie2-build -f {} '.format(assembly_clean)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '--seed {} '.format(random_seed_int)
            command += '{} && '.format(bt2index)
            command += 'bowtie2 -x {} '.format(bt2index)
            command += '-1 {} '.format(fastq_forward)
            command += '-2 {} '.format(fastq_reverse)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '-S {}'.format(sam)
        elif task_params['read_mapping_tool'] == 'bowtie2_very_sensitive':
            (fastq_forward, fastq_reverse) = self.deinterlace_raw_reads(fastq)
            bt2index = os.path.basename(assembly_clean) + '.bt2'
            command = 'bowtie2-build -f {} '.format(assembly_clean)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '--seed {} '.format(random_seed_int)
            command += '{} && '.format(bt2index)
            command += 'bowtie2 --very-sensitive -x {} '.format(bt2index)
            command += '-1 {} '.format(fastq_forward)
            command += '-2 {} '.format(fastq_reverse)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '-S {}'.format(sam)
        elif task_params['read_mapping_tool'] == 'minimap2':
            (fastq_forward, fastq_reverse) = self.deinterlace_raw_reads(fastq)
            command = 'minimap2 -ax sr -t {} '.format(self.MAPPING_THREADS)
            command += '--seed {} '.format(random_seed_int)
            command += '{} '.format(assembly_clean)
            command += '{} '.format(fastq_forward)
            command += '{} > '.format(fastq_reverse)
            command += '{}'.format(sam)
        elif task_params['read_mapping_tool'] == 'hisat2':
            (fastq_forward, fastq_reverse) = self.deinterlace_raw_reads(fastq)
            ht2index = os.path.basename(assembly_clean) + '.ht2'
            command = 'hisat2-build {} '.format(assembly_clean)
            command += '{} && '.format(ht2index)
            command += 'hisat2 -x {} '.format(ht2index)
            command += '-1 {} '.format(fastq_forward)
            command += '-2 {} '.format(fastq_reverse)
            command += '-S {} '.format(sam)
            command += '--seed {} '.format(random_seed_int)
            command += '--threads {}'.format(self.MAPPING_THREADS)
        log('running alignment command: {}'.format(command))
        out, err = self._run_command(command)

    def run_read_mapping_unpaired_mode(self, task_params, assembly_clean, fastq, sam):
        read_mapping_tool = task_params['read_mapping_tool']
        log("running {} mapping in single-end (unpaired) mode.".format(read_mapping_tool))
        random_seed_int = randint(0, 999999999)
        log("randomly selected seed (integer) used for read mapping is: {}".format(random_seed_int))
        if task_params['read_mapping_tool'] == 'bbmap_fast':
            log("Warning: bbmap does not support setting random seeds, so results are not reproducible.")
            command = 'bbmap.sh -Xmx{} '.format(self.BBMAP_MEM)
            command += 'threads={} '.format(self.MAPPING_THREADS)
            command += 'ref={} '.format(assembly_clean)
            command += 'in={} '.format(fastq)
            command += 'out={} '.format(sam)
            command += 'fast interleaved=false mappedonly nodisk overwrite'
        elif task_params['read_mapping_tool'] == 'bbmap_default':
            log("Warning: bbmap does not support setting random seeds, so results are not reproducible.")
            command = 'bbmap.sh -Xmx{} '.format(self.BBMAP_MEM)
            command += 'threads={} '.format(self.MAPPING_THREADS)
            command += 'ref={} '.format(assembly_clean)
            command += 'in={} '.format(fastq)
            command += 'out={} '.format(sam)
            command += 'interleaved=false mappedonly nodisk overwrite'
        elif task_params['read_mapping_tool'] == 'bbmap_very_sensitive':
            log("Warning: bbmap does not support setting random seeds, so results are not reproducible.")
            command = 'bbmap.sh -Xmx{} '.format(self.BBMAP_MEM)
            command += 'threads={} '.format(self.MAPPING_THREADS)
            command += 'ref={} '.format(assembly_clean)
            command += 'in={} '.format(fastq)
            command += 'out={} '.format(sam)
            command += 'vslow=true '
            command += 'interleaved=false mappedonly nodisk overwrite'
        elif task_params['read_mapping_tool'] == 'bowtie2_default':
            bt2index = os.path.basename(assembly_clean) + '.bt2'
            command = 'bowtie2-build -f {} '.format(assembly_clean)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '--seed {} '.format(random_seed_int)
            command += '{} && '.format(bt2index)
            command += 'bowtie2 -x {} '.format(bt2index)
            command += '-U {} '.format(fastq)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '-S {}'.format(sam)
        elif task_params['read_mapping_tool'] == 'bowtie2_very_sensitive':
            bt2index = os.path.basename(assembly_clean) + '.bt2'
            command = 'bowtie2-build -f {} '.format(assembly_clean)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '--seed {} '.format(random_seed_int)
            command += '{} && '.format(bt2index)
            command += 'bowtie2 --very-sensitive -x {} '.format(bt2index)
            command += '-U {} '.format(fastq)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '-S {}'.format(sam)
        elif task_params['read_mapping_tool'] == 'minimap2':
            command = 'minimap2 -ax sr -t {} '.format(self.MAPPING_THREADS)
            command += '--seed {} '.format(random_seed_int)
            command += '{} '.format(assembly_clean)
            command += '{} > '.format(fastq)
            command += '{}'.format(sam)
        elif task_params['read_mapping_tool'] == 'hisat2':
            ht2index = os.path.basename(assembly_clean) + '.ht2'
            command = 'hisat2-build {} '.format(assembly_clean)
            command += '{} && '.format(ht2index)
            command += 'hisat2 -x {} '.format(ht2index)
            command += '-U {} '.format(fastq)
            command += '-S {} '.format(sam)
            command += '--seed {} '.format(random_seed_int)
            command += '--threads {}'.format(self.MAPPING_THREADS)
        log('running alignment command: {}'.format(command))
        out, err = self._run_command(command)

    def convert_sam_to_sorted_and_indexed_bam(self, sam):
        # create bam files from sam files
        sorted_bam = os.path.abspath(sam).split('.sam')[0] + "_sorted.bam"

        command = 'samtools view -F 0x04 -uS {} | '.format(sam)
        command += 'samtools sort - -o {}'.format(sorted_bam)

        log('running samtools command to generate sorted bam: {}'.format(command))
        self._run_command(command)

        # verify we got bams
        if not os.path.exists(sorted_bam):
            log('Failed to find bam file\n{}'.format(sorted_bam))
            sys.exit(1)
        elif(os.stat(sorted_bam).st_size == 0):
            log('Bam file is empty\n{}'.format(sorted_bam))
            sys.exit(1)

        # index the bam file
        command = 'samtools index {}'.format(sorted_bam)

        log('running samtools command to index sorted bam: {}'.format(command))
        self._run_command(command)

        return sorted_bam

    def generate_alignment_bams(self, task_params, assembly_clean):
        """
            This function runs the selected read mapper and creates the
            sorted and indexed bam files from sam files using samtools.
        """

        reads_list = task_params['reads_list']

        (read_scratch_path, read_type) = self.stage_reads_list_file(reads_list)

        sorted_bam_file_list = []

        # list of reads files, can be 1 or more. assuming reads are either type unpaired or interleaved
        # will not handle unpaired forward and reverse reads input as seperate (non-interleaved) files

        for i in range(len(read_scratch_path)):
            fastq = read_scratch_path[i]
            fastq_type = read_type[i]

            sam = os.path.basename(fastq).split('.fastq')[0] + ".sam"
            sam = os.path.join(self.BINNER_RESULT_DIRECTORY, sam)

            if fastq_type == 'interleaved':  # make sure working - needs tests
                log("Running interleaved read mapping mode")
                self.run_read_mapping_interleaved_pairs_mode(task_params, assembly_clean, fastq, sam)
            else:  # running read mapping in single-end mode
                log("Running unpaired read mapping mode")
                self.run_read_mapping_unpaired_mode(task_params, assembly_clean, fastq, sam)

            sorted_bam = self.convert_sam_to_sorted_and_indexed_bam(sam)

            sorted_bam_file_list.append(sorted_bam)

        return sorted_bam_file_list

    def generate_make_coverage_table_command(self, task_params, sorted_bam_file_list):
        # create the depth file for this bam
        #
        min_contig_length = task_params['min_contig_length']
        sorted_bam = task_params['sorted_bam']

        depth_file_path = os.path.join(self.scratch, str('concoct_depth.txt'))
        command = '/kb/module/lib/kb_concoct/bin/jgi_summarize_bam_contig_depths '
        command += '--outputDepth {} '.format(depth_file_path)
        command += '--minContigLength {} '.format(min_contig_length)
        command += '--minContigDepth 1 {}'.format(sorted_bam)

        log('running summarize_bam_contig_depths command: {}'.format(command))
        self._run_command(command)

        return depth_file_path

    def fix_generate_concoct_command_ui_bug(self, task_params):
        # needed to get checkbox for UI to work with string objects, for some
        # reason strings are converted to numerics when running inside KBase UI.
        parameter_no_total_coverage = task_params['no_total_coverage']
        parameter_no_cov_normalization = task_params['no_cov_normalization']

        if task_params['no_total_coverage'] is 1:
            parameter_no_total_coverage = '--no_total_coverage'
        elif task_params['no_total_coverage'] is 0:
            parameter_no_total_coverage = ' '

        if task_params['no_cov_normalization'] is 1:
            parameter_no_cov_normalization = '--no_cov_normalization'
        elif task_params['no_cov_normalization'] is 0:
            parameter_no_cov_normalization = ' '

        return (parameter_no_total_coverage, parameter_no_cov_normalization)

    def generate_concoct_cut_up_fasta_command(self, task_params):
        """
        generate_command: concoct cut_up_fasta
        """
        contig_file_path = task_params['contig_file_path']
        contig_split_size = task_params['contig_split_size']
        contig_split_overlap = task_params['contig_split_overlap']

        log("\n\nRunning generate_concoct_cut_up_fasta_command")

        command = 'python {}/scripts/cut_up_fasta.py '.format(self.CONCOCT_BASE_PATH)
        command += '{} '.format(contig_file_path)
        command += '-c {} '.format(contig_split_size)
        command += '-o {} '.format(contig_split_overlap)
        command += '--merge_last -b temp.bed > {}/split_contigs.fa'.format(self.BINNER_RESULT_DIRECTORY)
        log('Generated concoct_cut_up_fasta command: {}'.format(command))

        self._run_command(command)

    def generate_concoct_coverage_table_from_bam(self, task_params):
        """
        generate_command: concoct generate coverage table
        """
        log("\n\nRunning generate_concoct_coverage_table_from_bam")
        command = 'python {}/scripts/concoct_coverage_table.py temp.bed '.format(self.CONCOCT_BASE_PATH)
        command += '{}/*_sorted.bam > '.format(self.BINNER_RESULT_DIRECTORY)
        command += '{}/coverage_table.tsv'.format(self.BINNER_RESULT_DIRECTORY)
        log('Generated concoct generate coverage table from bam command: {}'.format(command))

        self._run_command(command)

    def generate_concoct_command(self, task_params):
        """
        generate_command: concoct
        """

        min_contig_length = task_params['min_contig_length']
        kmer_size = task_params['kmer_size']
        max_clusters_for_vgmm = task_params['max_clusters_for_vgmm']
        max_iterations_for_vgmm = task_params['max_iterations_for_vgmm']
        total_percentage_pca = task_params['total_percentage_pca']
        parameter_no_total_coverage, parameter_no_cov_normalization = \
            self.fix_generate_concoct_command_ui_bug(task_params)

        log("\n\nRunning generate_concoct_command")
        command = 'python {}/bin/concoct '.format(self.CONCOCT_BASE_PATH)
        command += '--composition_file {}/split_contigs.fa '.format(self.BINNER_RESULT_DIRECTORY)
        command += '-l {} '.format(min_contig_length)
        command += '-b {} '.format(self.BINNER_RESULT_DIRECTORY)
        command += '--coverage_file {}/coverage_table.tsv '.format(self.BINNER_RESULT_DIRECTORY)
        command += '-t {} '.format(self.MAPPING_THREADS)
        command += '-k {} '.format(kmer_size)
        command += '-c {} '.format(max_clusters_for_vgmm)
        command += '-i {} '.format(max_iterations_for_vgmm)
        command += '--total_percentage_pca {} '.format(total_percentage_pca)
        command += '{} '.format(parameter_no_cov_normalization)
        command += '{}'.format(parameter_no_total_coverage)
        log('Generated concoct command: {}'.format(command))

        self._run_command(command)

    def generate_concoct_post_clustering_merging_command(self, task_params):
        """
        generate_command: concoct post cluster merging
        """
        min_contig_length = task_params['min_contig_length']
        log("\n\nRunning generate_concoct_post_clustering_merging_command")

        command = 'python {}/scripts/merge_cutup_clustering.py '.format(self.CONCOCT_BASE_PATH)
        command += '{}/clustering_gt{}.csv > '.format(self.BINNER_RESULT_DIRECTORY, min_contig_length)
        command += '{}/clustering_merged.csv'.format(self.BINNER_RESULT_DIRECTORY)
        log('Generated generate_concoct_post_clustering_merging command: {}'.format(command))

        self._run_command(command)

    def generate_concoct_extract_fasta_bins_command(self, task_params):
        """
        generate_command: concoct extract_fasta_bins
        """
        log("\n\nRunning generate_concoct_extract_fasta_bins_command")

        contig_file_path = task_params['contig_file_path']

        bin_result_directory = self.BINNER_RESULT_DIRECTORY + '/' + self.BINNER_BIN_RESULT_DIR
        self._mkdir_p(bin_result_directory)
        command = 'python {}/scripts/extract_fasta_bins.py '.format(self.CONCOCT_BASE_PATH)
        command += '{} '.format(contig_file_path)
        command += '{}/clustering_merged.csv '.format(self.BINNER_RESULT_DIRECTORY)
        command += '--output_path {}/{}'.format(self.BINNER_RESULT_DIRECTORY, self.BINNER_BIN_RESULT_DIR)
        log('Generated generate_concoct_extract_fasta_bins_command command: {}'.format(command))

        self._run_command(command)

    def rename_and_standardize_bin_names(self, task_params):
        """
        generate_command: generate renamed bins
        """
        log("\n\nRunning rename_and_standardize_bin_names")
        path_to_concoct_result_bins = os.path.abspath(self.BINNER_RESULT_DIRECTORY) + \
            '/' + self.BINNER_BIN_RESULT_DIR + '/'
        for dirname, subdirs, files in os.walk(path_to_concoct_result_bins):
            for file in files:
                if file.endswith('.fa'):
                    os.rename(os.path.abspath(path_to_concoct_result_bins) + '/' +
                              file, os.path.abspath(path_to_concoct_result_bins) + '/bin.' +
                              file.split('.fa')[0].zfill(3) + '.fasta')  # need to change to 4 digits

    def revert_fasta_headers(self, task_params):
        path_to_concoct_result_bins = os.path.join(
            self.scratch,
            self.BINNER_RESULT_DIRECTORY,
            self.BINNER_BIN_RESULT_DIR
            )
        for filepath in [os.path.join(path_to_concoct_result_bins, filename) for filename in os.listdir(path_to_concoct_result_bins)]:
            self.fhr.revert_fasta_headers(filepath)


    def make_binned_contig_summary_file_for_binning_apps(self, task_params):
        """
        generate_command: generate binned contig summary command
        """
        log("\n\nRunning make_binned_contig_summary_file_for_binning_apps")
        path_to_concoct_result = os.path.abspath(self.BINNER_RESULT_DIRECTORY)
        path_to_concoct_result_bins = '{}/{}/'.format(path_to_concoct_result, self.BINNER_BIN_RESULT_DIR)
        path_to_summary_file = path_to_concoct_result_bins + 'binned_contig.summary'
        with open(path_to_summary_file, 'w+') as f:
            f.write("Bin name\tCompleteness\tGenome size\tGC content\n")
            for dirname, subdirs, files in os.walk(path_to_concoct_result_bins):
                for file in files:
                    if file.endswith('.fasta'):
                        genome_bin_fna_file = os.path.join(self.BINNER_BIN_RESULT_DIR, file)
                        bbstats_output_file = os.path.join(self.scratch, self.BINNER_RESULT_DIRECTORY,
                                                           genome_bin_fna_file).split('.fasta')[0] + ".bbstatsout"
                        bbstats_output = self.generate_stats_for_genome_bins(task_params,
                                                                             genome_bin_fna_file,
                                                                             bbstats_output_file)
                        f.write('{}\t0\t{}\t{}\n'.format(genome_bin_fna_file.split("/")[-1],
                                                         bbstats_output['contig_bp'],
                                                         bbstats_output['gc_avg']))
        f.close()
        log('Finished make_binned_contig_summary_file_for_binning_apps function')

    def generate_output_file_list(self, result_directory):
        """
        generate_output_file_list: zip result files and generate file_links for report
        """
        log('Start packing result files')
        output_files = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file = os.path.join(output_directory, 'concoct_result.zip')

        with zipfile.ZipFile(result_file, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:

            for dirname, subdirs, files in os.walk(result_directory):
                for file in files:
                    if (file.endswith('.sam') or
                        file.endswith('.bam') or
                        file.endswith('.bai') or
                       file.endswith('.summary')):
                            continue
                    if (dirname.endswith(self.BINNER_BIN_RESULT_DIR)):
                            continue
                    zip_file.write(os.path.join(dirname, file), file)
                if (dirname.endswith(self.BINNER_BIN_RESULT_DIR)):
                    baseDir = os.path.basename(dirname)
                    for file in files:
                        full = os.path.join(dirname, file)
                        zip_file.write(full, os.path.join(baseDir, file))

        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'Files generated by CONCOCT App'})

        return output_files

    def generate_html_report(self, result_directory, assembly_ref, binned_contig_obj_ref):
        """
        generate_html_report: generate html summary report
        """

        log('Start generating html report')
        html_report = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'report.html')

        # get summary data from existing assembly object and bins_objects
        Summary_Table_Content = ''
        Overview_Content = ''
        (binned_contig_count, input_contig_count, total_bins_count) = \
            self.generate_overview_info(assembly_ref, binned_contig_obj_ref, result_directory)

        Overview_Content += '<p>Binned contigs: {}</p>'.format(binned_contig_count)
        Overview_Content += '<p>Input contigs: {}</p>'.format(input_contig_count)
        Overview_Content += '<p>Number of bins: {}</p>'.format(total_bins_count)

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'report_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Overview_Content</p>',
                                                          Overview_Content)
                report_template = report_template.replace('Summary_Table_Content',
                                                          Summary_Table_Content)
                result_file.write(report_template)

        html_report.append({'path': result_file_path,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for kb_concoct App'})
        return html_report

    def generate_overview_info(self, assembly_ref, binned_contig_obj_ref, result_directory):
        """
        _generate_overview_info: generate overview information from assembly and binnedcontig
        """

        # get assembly and binned_contig objects that already have some data populated in them
        assembly = self.dfu.get_objects({'object_refs': [assembly_ref]})['data'][0]
        binned_contig = self.dfu.get_objects({'object_refs': [binned_contig_obj_ref]})['data'][0]

        input_contig_count = assembly.get('data').get('num_contigs')
        binned_contig_count = 0
        total_bins_count = 0
        total_bins = binned_contig.get('data').get('bins')
        total_bins_count = len(total_bins)
        for bin in total_bins:
            binned_contig_count += len(bin.get('contigs'))

        return (binned_contig_count, input_contig_count, total_bins_count)

    def generate_report(self, binned_contig_obj_ref, task_params):
        """
        generate_report: generate summary report
        """
        log('Generating report')

        result_directory = os.path.join(self.scratch, "concoct_output_dir")

        task_params['result_directory'] = result_directory

        output_files = self.generate_output_file_list(task_params['result_directory'])

        output_html_files = self.generate_html_report(task_params['result_directory'],
                                                      task_params['assembly_ref'],
                                                      binned_contig_obj_ref)

        report_params = {
            'message': '',
            'workspace_name': task_params['workspace_name'],
            'file_links': output_files,
            'html_links': output_html_files,
            'direct_html_link_index': 0,
            'html_window_height': 266,
            'report_object_name': 'kb_concoct_report_' + str(uuid.uuid4())
        }

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def create_dict_from_depth_file(self, depth_file_path):
        # keep contig order (required by metabat2)
        depth_file_dict = {}
        with open(depth_file_path, 'r') as f:
            header = f.readline().rstrip().split("\t")
            # print('HEADER1 {}'.format(header))
            # map(str.strip, header)
            for line in f:
                # deal with cases were fastq name has spaces.Assume first
                # non white space word is unique and use this as ID.
                # line = line.rstrip()
                vals = line.rstrip().split("\t")
                if ' ' in vals[0]:
                    ID = vals[0].split()[0]
                else:
                    ID = vals[0]
                depth_file_dict[ID] = vals[1:]
            depth_file_dict['header'] = header
        return depth_file_dict

    def run_concoct(self, task_params):
        """
        run_concoct: concoct app

        required params:
            assembly_ref: Metagenome assembly object reference
            binned_contig_name: BinnedContig object name and output file header
            workspace_name: the name of the workspace it gets saved to.
            reads_list: list of reads object (PairedEndLibrary/SingleEndLibrary)
            upon which CONCOCT will be run

        optional params:
            min_contig_length: minimum contig length; default 1000

            ref: https://github.com/BinPro/CONCOCT/blob/develop/README.md
        """
        log('--->\nrunning ConcoctUtil.run_concoct\n' +
            'task_params:\n{}'.format(json.dumps(task_params, indent=1)))

        self._validate_run_concoct_params(task_params)

        # get assembly
        contig_file = self._get_contig_file(task_params['assembly_ref'])
        task_params['contig_file_path'] = contig_file

        # clean the assembly file so that there are no spaces in the fasta headers
        assembly_clean = self.retrieve_and_clean_assembly(task_params)
        task_params['contig_file_path'] = assembly_clean

        # get reads
        (reads_list_file, read_type) = self.stage_reads_list_file(task_params['reads_list'])
        task_params['read_type'] = read_type
        task_params['reads_list_file'] = reads_list_file

        # prep result directory
        result_directory = os.path.join(self.scratch, self.BINNER_RESULT_DIRECTORY)
        self._mkdir_p(result_directory)

        cwd = os.getcwd()
        log('changing working dir to {}'.format(result_directory))
        os.chdir(result_directory)

        # set up tasks for kbparallel to run alignments
        # this also submits run_alignments function in parallel
        # self.set_up_parallel_tasks(params)

        # run alignments, and update input contigs to use the clean file
        # this function has an internal loop to generate a sorted bam file for each input read file
        #
        # self.set_up_parallel_tasks(task_params)

        self.generate_alignment_bams(task_params, assembly_clean)

        # not used right now
        # depth_file_path = self.generate_make_coverage_table_command(task_params, sorted_bam_file_list)
        # depth_dict = self.create_dict_from_depth_file(depth_file_path)

        # run concoct prep, cut up fasta input
        self.generate_concoct_cut_up_fasta_command(task_params)

        # run concoct make coverage table from bam
        self.generate_concoct_coverage_table_from_bam(task_params)

        # run concoct prep and concoct
        self.generate_concoct_command(task_params)

        # run concoct post cluster merging command
        self.generate_concoct_post_clustering_merging_command(task_params)

        # run extract bins command
        self.generate_concoct_extract_fasta_bins_command(task_params)

        # run fasta renaming
        self.rename_and_standardize_bin_names(task_params)

        # revert fasta headers in bins
        self.revert_fasta_headers(task_params)

        self.make_binned_contig_summary_file_for_binning_apps(task_params)

        # file handling and management
        os.chdir(cwd)
        log('changing working dir to {}'.format(cwd))

        log('Saved result files to: {}'.format(result_directory))
        log('Generated files:\n{}'.format('\n'.join(os.listdir(result_directory))))

        # make new BinnedContig object and upload to KBase
        generate_binned_contig_param = {
            'file_directory': os.path.join(result_directory, self.BINNER_BIN_RESULT_DIR),
            'assembly_ref': task_params['assembly_ref'],
            'binned_contig_name': task_params['binned_contig_name'],
            'workspace_name': task_params['workspace_name']
        }

        binned_contig_obj_ref = \
            self.mgu.file_to_binned_contigs(generate_binned_contig_param).get('binned_contig_obj_ref')

        # generate report
        reportVal = self.generate_report(binned_contig_obj_ref, task_params)
        returnVal = {
            'result_directory': result_directory,
            'binned_contig_obj_ref': binned_contig_obj_ref
        }
        returnVal.update(reportVal)

        return returnVal
