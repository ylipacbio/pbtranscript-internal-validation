from pbtranscript.Utils import execute, realpath, mkdir, rmpath
import os.path as op
import os

def bam2fasta(i_bam, o_fasta):
    cmd = 'bam2fasta {0} -o {1} && gunzip {1}.fasta.gz && mv {1}.fasta {2}'.\
        format(i_bam, o_fasta + '.tmp', o_fasta)
    execute(cmd)

class SMRTLinkIsoSeq3Files(object):
    def __init__(self, root_dir):
        self.root_dir = realpath(root_dir)
        def f(*args):
            return op.join(self.root_dir, 'tasks', *args)
        self.polished_bam = f('pbcoretools.tasks.gather_transcripts-1', 'file.transcriptset.xml')
        self.unpolished_bam = f('isoseq3.tasks.cluster-0', 'unpolished.bam')
        self.isoseq_flnc_bam = f('isoseq3.tasks.cluster-0', 'unpolished.flnc.bam')
        self.hq_isoforms_transcriptset = f('pbcoretools.tasks.gather_transcripts-1', 'file.transcriptset.xml')
        self.hq_isoforms_fa = f('pbcoretools.tasks.bam2fasta_transcripts-0', 'hq_transcripts.fasta')
        self.lq_isoforms_fa = f('pbcoretools.tasks.bam2fasta_transcripts-0', 'lq_transcripts.fasta')
        self.hq_isoforms_fq = f('pbcoretools.tasks.bam2fastq_transcripts-0', 'hq_transcripts.fastq')
        self.lq_isoforms_fq = f('pbcoretools.tasks.bam2fastq_transcripts-0', 'lq_transcripts.fastq')
        self.pbscala_sh = f('../pbscala-job.sh')
        self.subreads_xml = [fn for fn in os.listdir(f('../entry-points')) if op.isfile(op.join(f('../entry-points/'), fn)) and fn.endswith('subreadset.xml')][0]
        self.ccs_xml = f('pbcoretools.tasks.gather_ccsset-1', 'file.consensusreadset.xml')

    def export_isoseq_flnc_fa(self, isoseq_flnc_fa):
        bam2fasta(self.isoseq_flnc_bam,  isoseq_flnc_fa)

    def export_unpolished_fa(self, unpolished_fa):
        bam2fasta(self.unpolished_bam,  unpolished_fa)
