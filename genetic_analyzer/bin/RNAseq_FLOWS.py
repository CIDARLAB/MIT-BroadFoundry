"""
Created on Mar 23, 2016

@author: dbg
         A Cristofaro
"""

import json
import argparse
import csv
from collections import OrderedDict
import sys

##TODO Based on input, how many steps to complete (Mapping, Counts + Normalization, DE, Part Profiling)

BIN_PATH = '/btl/foundry/software/rnaSeq/'
MASTER_SCRIPT = '/btl/foundry/software/flows_tunnel.sh'


def main():
    # This block parses command-line arguments
    parser = argparse.ArgumentParser(description="Generate Job JSON")
    arggroup = parser.add_argument_group("Arguments")
    arggroup.add_argument('--settings', '-s', action="store", default=None, dest="settingsFile",
                          help='Full path to settings file')
    arggroup.add_argument('--groups', '-g', action="store", default=None, dest="groupsFile",
                          help='Full path to groups file')
    arggroup.add_argument('--format', '-f', action="store", default='json', dest="format",
                          help='Output format, json or shell')
    options = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    settingsFile = options.settingsFile
    groupsFile = options.groupsFile
    format = options.format

    rnaseqTask = createRNAseqJob(settingsFile, groupsFile)

    # Prints JSON
    if format == 'json':
        print(rnaseqTask.asJSON())
    # Prints SHELL COMMANDS
    elif format == 'shell':
        rnaseqTask.asSHELL()


def createRNAseqJob(settingsFile, group_file):
    """
    Differential expression sample numbers correspond to sorted settings sample names
    Changing the order of the settings.txt file has no effect on sample numbering
    """
    # Read the settings file
    settings = load_settings(settingsFile)

    # Instantiate the Task.  It will "hold" all jobs.
    theTask = Task("RNA-seq Task")

    # Settings is sorted by keys due to the numbering method utilized in the DE-analysis method
    # Sorted ordering is mimicked in genetic_analyzer.de_analysis during group comparison set-up
    sampleNames = [x for x in sorted(settings)]
    sampleNames.remove('None')
    required_files = []
    for sample in settings:
        if sample != None:
            required_files.append(settings[sample]['R1_fastq_file'])
            required_files.append(settings[sample]['R2_fastq_file'])
            required_files.append(settings[sample]['bed_file'])
            required_files.append(settings[sample]['fasta_file'])
            required_files.append(settings[sample]['gff_file'])
    for file in set(required_files):
        if file != "":
            theTask.addRequiredFile(file)
    project_reference_name = settings[sample]['fasta_file'].split('/')[-1].split('.')[0]
    promoter_part_expression_datafile = '/btl/foundry/data/rna-seq/' + project_reference_name + '/promoter_reu.txt'
    terminator_part_expression_datafile = '/btl/foundry/data/rna-seq/' + project_reference_name + '/terminator_reu.txt'
    theTask.addRequiredFile(promoter_part_expression_datafile)
    theTask.addRequiredFile(terminator_part_expression_datafile)
    # Set up sub directories in the working directory
    theTask.requiredSubDirectories = ["log_out", "log_err", "results", "alignment"]
    theTask.requiredSubDirectories.extend(["alignment/" + sample for sample in sampleNames])
    theTask.requiredSubDirectories.extend(["results/" + sample for sample in sampleNames])

    # Jobs are constructed to call a shell script that runs qsub with appropriate
    # parameters which executes a script file that sources dotkits then executes
    # the cmd

    ###############################################################################
    # Step 01: Map Reads
    ###############################################################################

    python_script = BIN_PATH + 'map_reads.py'

    for sampleName in sampleNames:
        jobName = "MapReads_%s" % (sampleName)
        outputFile = 'log_out/01_map_reads_%s.out.log' % (sampleName)
        errorFile = 'log_err/01_map_reads_%s.err.log' % (sampleName)

        cmd = 'python %s -settings %s -samples %s' % (python_script, settingsFile, sampleName)

        job = Job(jobName, cmd=cmd)
        job.outputFile = outputFile
        job.errorFile = errorFile
        job.dependencies = []  # Doesn't depend on any other jobs
        job.memory = 8000

        theTask.addJob(job)

    ###############################################################################
    # Step 02: Count Reads
    ###############################################################################

    python_script = BIN_PATH + 'count_reads.py'

    for sampleName in sampleNames:
        jobName = "CountReads_%s" % (sampleName)
        outputFile = 'log_out/02_count_reads_%s.out.log' % (sampleName)
        errorFile = 'log_err/02_count_reads_%s.err.log' % (sampleName)
        cmd = 'python %s -settings %s -samples %s -feature gene -attribute Name -strand_opt reverse' % (
        python_script, settingsFile, sampleName)

        job = Job(jobName, cmd=cmd)
        job.outputFile = outputFile
        job.errorFile = errorFile
        job.dependencies = ['MapReads_%s' % (sampleName)]  # Only depends on one job (in this case)
        features = ['gene', 'promoter', 'terminator', 'ribozyme', 'promoter_unit', 'transcript']
        counts = []
        for f in features:
            counts.append('results/' + sampleName + '/' + sampleName + '.{feat}'.format(feat=f) + '.counts.txt')
        genelengths = ['results/' + sampleName + '/' + sampleName + '.gene.lengths.txt']
        mappedreads = ['results/' + sampleName + '/' + sampleName + '.mapped.reads.txt']
        job.requiredResultFiles.extend(counts)
        job.requiredResultFiles.extend(genelengths)
        job.requiredResultFiles.extend(mappedreads)
        job.memory = 16000

        theTask.addJob(job)

    ###############################################################################
    # Step 03: Fragment Distributions
    ###############################################################################

    python_script = BIN_PATH + 'fragment_distributions.py'

    for sampleName in sampleNames:
        jobName = "FragDist_%s" % (sampleName)
        outputFile = 'log_out/03_fragment_distributions_%s.out.log' % (sampleName)
        errorFile = 'log_err/03_fragment_distributions_%s.err.log' % (sampleName)
        cmd = 'python %s -settings %s -samples %s' % (python_script, settingsFile, sampleName)

        job = Job(jobName, cmd=cmd)
        job.outputFile = outputFile
        job.errorFile = errorFile
        job.dependencies = ['MapReads_%s' % (sampleName)]
        fragmentdistribution = ['results/' + sampleName + '/' + sampleName + '.fragment.distribution.txt']
        job.requiredResultFiles.extend(fragmentdistribution)

        theTask.addJob(job)

    ###############################################################################
    # Step 04: Transcription_Profile
    ###############################################################################

    python_script = BIN_PATH + 'transcription_profile.py'

    for sampleName in sampleNames:
        jobName = "TranProf_%s" % (sampleName)
        outputFile = 'log_out/04_transcription_profile_%s.out.log' % (sampleName)
        errorFile = 'log_err/04_transcription_profile_%s.err.log' % (sampleName)
        cmd = 'python %s -settings %s -samples %s' % (python_script, settingsFile, sampleName)

        job = Job(jobName, cmd=cmd)
        job.outputFile = outputFile
        job.errorFile = errorFile
        job.dependencies = ['MapReads_%s' % (sampleName)]
        fwdprofile = ['results/' + sampleName + '/' + sampleName + '.fwd.profiles.txt']
        revprofile = ['results/' + sampleName + '/' + sampleName + '.rev.profiles.txt']
        job.requiredResultFiles.extend(fwdprofile)
        job.requiredResultFiles.extend(revprofile)

        theTask.addJob(job)

    ###############################################################################
    # Step 05: Profile Normalization
    ###############################################################################

    python_script = BIN_PATH + 'profile_normalization.py'

    for sampleName in sampleNames:
        jobName = "ProfNorm_%s" % (sampleName)
        outputFile = 'log_out/05_profile_normalization_%s.out.log' % (sampleName)
        errorFile = 'log_err/05_profile_normalization_%s.err.log' % (sampleName)
        cmd = 'python %s -settings %s -samples %s' % (python_script, settingsFile, sampleName)

        job = Job(jobName, cmd=cmd)
        job.outputFile = outputFile
        job.errorFile = errorFile
        job.dependencies = ['ReadAnalysis', 'TranProf_%s' % (sampleName)]
        fwdprofnorm = ['results/' + sampleName + '/' + sampleName + '.fwd.norm.profiles.txt']
        revprofnorm = ['results/' + sampleName + '/' + sampleName + '.rev.norm.profiles.txt']
        job.requiredResultFiles.extend(fwdprofnorm)
        job.requiredResultFiles.extend(revprofnorm)

        theTask.addJob(job)

    ###############################################################################
    # Step 06: Read Analysis
    ###############################################################################

    python_script = BIN_PATH + 'read_analysis.py'

    jobName = "ReadAnalysis"
    outputFile = 'log_out/06_read_analysis.out.log'
    errorFile = 'log_err/06_read_analysis.err.log'
    cmd = 'python %s -settings %s -bin_path %s' % (python_script, settingsFile, BIN_PATH)

    job = Job(jobName, cmd=cmd)
    job.outputFile = outputFile
    job.errorFile = errorFile
    job.dependencies = ['CountReads_%s' % (sample) for sample in sampleNames]
    matrices = ['results/counts.matrix.txt', 'results/mapped.reads.matrix.txt', 'results/gene.lengths.matrix.txt',
                'results/norm.factors.matrix.txt', 'results/fpkm.normed.matrix.txt']
    job.requiredResultFiles.extend(matrices)

    theTask.addJob(job)

    ###############################################################################
    # Step 07: Differential Expression
    ###############################################################################

    python_script = BIN_PATH + 'de_analysis.py'

    de_groups = load_groups(group_file)
    for c in range(len(de_groups[0])):
        jobName = "DEanalysis_%s" % (c)
        outputFile = 'log_out/07_de_analysis_%s.out.log' % (c)
        errorFile = 'log_err/07_de_analysis_%s.err.log' % (c)
        cmd = 'python %s -settings %s -group1 %s -group2 %s -output_prefix group1_vs_group2_%s -bin_path %s' % (
        python_script, settingsFile, de_groups[0][c], de_groups[1][c], c, BIN_PATH)

        job = Job(jobName, cmd=cmd)
        job.outputFile = outputFile
        job.errorFile = errorFile
        job.dependencies = ['ReadAnalysis']
        job.requiredResultFiles.append('results/group1_vs_group2_%s.de.analysis.txt' % c)

        theTask.addJob(job)

    ###############################################################################
    # Step 08: Profile Analysis
    ###############################################################################

    python_script = BIN_PATH + 'part_profile_analysis_circuit.py'

    jobName = "PartProf"
    outputFile = 'log_out/08_part_profile_analysis_circuit.out.log'
    errorFile = 'log_err/08_part_profile_analysis_circuit.err.log'
    cmd = 'python %s -settings %s' % (python_script, settingsFile)

    job = Job(jobName, cmd=cmd)
    job.outputFile = outputFile
    job.errorFile = errorFile
    dependencies = ['FragDist_%s' % (sample) for sample in sampleNames] + ['ProfNorm_%s' % (sample) for sample in
                                                                           sampleNames]
    for j in dependencies:
        job.dependencies.append(j)
    profileperf = ['results/promoter.profile.perf.txt', 'results/ribozyme.profile.perf.txt',
                   'results/terminator.profile.perf.txt']
    sampleprofileperf = []
    for sample in sampleNames:
        sampleprofileperf.append('results/' + sample + '/' + sample + '.promoter.profile.perf.txt')
        sampleprofileperf.append('results/' + sample + '/' + sample + '.ribozyme.profile.perf.txt')
        sampleprofileperf.append('results/' + sample + '/' + sample + '.terminator.profile.perf.txt')
    job.requiredResultFiles.extend(profileperf)
    job.requiredResultFiles.extend(sampleprofileperf)

    theTask.addJob(job)

    ###############################################################################
    # Step 09: Part Performance Comparison Graphing
    ###############################################################################

    python_script = BIN_PATH + 'comparison_plots.py'

    jobName = "CompGraph"
    outputFile = 'log_out/09_comparison_graphs.out.log'
    errorFile = 'log_err/09_comparison_graphs.err.log'
    cmd = 'python %s -settings %s' % (python_script, settingsFile)

    job = Job(jobName, cmd=cmd)
    job.outputFile = outputFile
    job.errorFile = errorFile
    job.dependencies.append('PartProf')
    job.requiredResultFiles.extend(
        ['results/promoter_comparisons_by_part.png', 'results/promoter_comparisons_by_sample.png',
         'results/terminator_comparisons_by_part.png', 'results/terminator_comparisons_by_sample.png'])

    theTask.addJob(job)

    ###############################################################################
    # Step 10: Variant Analysis
    ###############################################################################

    python_script = BIN_PATH + 'variant_analysis.py'

    for sampleName in sampleNames:
        jobName = "VarCall_%s" % (sampleName)
        outputFile = 'log_out/10_variant_analysis_%s.out.log' % (sampleName)
        errorFile = 'log_err/10_variant_analysis_%s.err.log' % (sampleName)
        cmd = 'python %s -settings %s -samples %s' % (python_script, settingsFile, sampleName)

        job = Job(jobName, cmd=cmd)
        job.outputFile = outputFile
        job.errorFile = errorFile
        job.dependencies = ['MapReads_%s' % (sampleName)]
        marked_dupes = ['alignment/' + sampleName + '/' + sampleName + '.mdup.bam']
        rg_bams = ['alignment/' + sampleName + '/' + sampleName + '.mdup.rg.bam']
        rg_bais = ['alignment/' + sampleName + '/' + sampleName + '.mdup.rg.bam.bai']
        variants = ['results/' + sampleName + '/' + sampleName + '.vcf']
        job.memory = 8000
        job.requiredResultFiles.extend(marked_dupes)
        job.requiredResultFiles.extend(rg_bams)
        job.requiredResultFiles.extend(rg_bais)
        job.requiredResultFiles.extend(variants)

        theTask.addJob(job)

    ###############################################################################
    # Step 11: Context Analysis
    ###############################################################################

    python_script = BIN_PATH + 'context_effects.py'

    jobName = "Context"
    outputFile = 'log_out/10_context_analysis.out.log'
    errorFile = 'log_err/10_context_analysis.err.log'
    cmd = 'python %s -settings %s' % (python_script, settingsFile)

    job = Job(jobName, cmd=cmd)
    job.outputFile = outputFile
    job.errorFile = errorFile
    dependencies = ['VarCall_%s' % sample for sample in sampleNames]
    for c in range(len(de_groups[0])):
        dependencies += ["DEanalysis_%s" % c]
    dependencies += ['PartProf']
    job.dependencies = dependencies
    job.requiredResultFiles.extend(['results/context_data.txt'])

    theTask.addJob(job)

    return theTask


def load_settings(filename):
    settings = {}
    data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
    # Ignore header
    header = next(data_reader)
    """ Settings dict specification
    'flask_1': {   'R1_fastq_file':
        '/btl/projects/Foundry/Yongjin/Tom/circuit_fastq/H9FULADXX.ATTATGTT.unmapped.1_Circuit_sample_T1_2.fastq',
                   'R2_fastq_file':
        '/btl/projects/Foundry/Yongjin/Tom/circuit_fastq/H9FULADXX.ATTATGTT.unmapped.2_Circuit_sample_T1_2.fastq',
                   'bed_file': './data/bed/0x58v50.bed',
                   'fasta_file': './data/fasta/0x58v50.fasta',
                   'gff_file': './data/gff/0x58v50.gff',
                   'output_path': './results/flask_1/',
                   'temp_path': '/broad/hptmp/acristo/circuit_flask_1/'},
    """
    for row in data_reader:
        if len(row) == len(header):
            sample = row[0]
            sample_data = {}
            for el_idx, el in enumerate(header[1:]):
                sample_data[el] = row[el_idx + 1]
            settings[sample] = sample_data
    return settings


def load_groups(group_file):
    with open(group_file) as data:
        data = data.readlines()
        # Ignore header
        data = [x.rstrip() for x in data[1:]]
        de_groups = [[], []]
        # Append group 1 to list at index 0 and group 2 to second list at index 1 of de_groups
        try:
            for comparison in data:
                group1, group2 = comparison.split()
                de_groups[0].append(group1)
                de_groups[1].append(group2)
        except:
            raise ValueError("Need 2 groups for each comparison, group1 and group2")
    return de_groups


class Task(object):
    def __init__(self, name=None):
        self.requiredFiles = []
        self.name = name
        self.jobs = []
        self.environment = {}
        self.abortAllIfFailed = False
        self.requiredSubDirectories = []

    def asObject(self):
        a = OrderedDict([('taskName', self.name),
                         ('requiredSubDirectories', self.requiredSubDirectories),
                         ('requiredFiles', self.requiredFiles),
                         ('environment', self.environment),
                         ('abortAllIfFailed', self.abortAllIfFailed),
                         ('jobs', [x.asObject() for x in self.jobs]),
                         ('jobDependencies', {x.name: x.dependencies for x in self.jobs})])
        return a

    def addJob(self, job):
        self.jobs.append(job)

    def addRequiredFile(self, fname):
        if type(fname) == type([]):
            self.requiredFiles.extend(fname)
        else:
            self.requiredFiles.append(fname)

    def asJSON(self):
        return json.dumps(self.asObject(), indent=4)

    def asSHELL(self):
        for job in self.jobs:
            job.fixargs()
            if job.errorFile:
                cmd = job.cmd + ' ' + ' '.join(job.args) + " > " + job.outputFile + " 2> " + job.errorFile
            else:
                cmd = job.cmd + ' ' + ' '.join(job.args) + " > " + job.outputFile
            print(cmd)


class Job(object):
    def __init__(self, name, cmd=None, outputFile=None):
        self.name = name
        self.cmd = cmd
        self.args = []
        self.dependencies = []
        self.requiredResultFiles = []
        self.abortAllIfFailed = True
        self.outputFile = outputFile
        self.errorFile = None
        self.softRunDurationLimit = 0
        self.hardRunDurationLimit = 0
        self.fixed = False
        self.memory = None

    def fixargs(self):
        if self.fixed:
            return
        if not self.cmd:
            raise "CmdNotSpecifiedError"
        if not self.outputFile:
            raise "OutPutFileNotSpecifiedError"

        toks = self.cmd.split()
        self.cmd = MASTER_SCRIPT
        self.args = toks
        self.fixed = True

    def asObject(self):
        self.fixargs()
        a = OrderedDict([('jobName', self.name),
                         ('remoteCommand', self.cmd),
                         ('commandArgs', self.args),
                         ('outputFile', self.outputFile),
                         ('errorFile', self.errorFile),
                         ('softRunDurationLimit', self.softRunDurationLimit),
                         ('hardRunDurationLimit', self.hardRunDurationLimit),
                         ('memoryRequired', self.memory),
                         ('requiredResultFiles', self.requiredResultFiles)])
        return a

    def asJSON(self):
        return json.dumps(self.asObject(), indent=4, sort_keys=False)


if __name__ == '__main__':
    main()
