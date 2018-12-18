import os.path as op
from isocollapse.independent.system import execute, realpath, mkdir, rmpath
from ..Utils import get_subread_xml_from_job_path, _bam2fastx


class SMRTLinkIsoSeq3Files(object):
    def __init__(self, root_dir):
        self.root_dir = realpath(root_dir)

        def f(*args):
            return op.join(self.root_dir, 'tasks', *args)
        self.polished_bam = f('pbcoretools.tasks.gather_transcripts-1', 'file.transcriptset.xml')
        self.unpolished_bam = f('isoseq3.tasks.cluster-0', 'unpolished.bam')
        self.hq_isoforms_transcriptset = f('pbcoretools.tasks.gather_transcripts-1', 'file.transcriptset.xml')
        self.hq_isoforms_fa = f('pbcoretools.tasks.bam2fasta_transcripts-0', 'hq_transcripts.fasta')
        self.lq_isoforms_fa = f('pbcoretools.tasks.bam2fasta_transcripts-0', 'lq_transcripts.fasta')
        self.hq_isoforms_fq = f('pbcoretools.tasks.bam2fastq_transcripts-0', 'hq_transcripts.fastq')
        self.lq_isoforms_fq = f('pbcoretools.tasks.bam2fastq_transcripts-0', 'lq_transcripts.fastq')
        self.pbscala_sh = f('../pbscala-job.sh')
        try:
            self.subreads_xml = get_subread_xml_from_job_path(root_dir)
        except Exception:
            self.subreads_xml = None
        self.ccs_xml = f('pbcoretools.tasks.gather_ccsset-1', 'file.consensusreadset.xml')
        self.ccs_report_json = f('pbreports.tasks.ccs_report-0', 'ccs_report.json')
        self.hq_transcript_ds = f('pbcoretools.tasks.consolidate_transcripts-0',
                                  'combined.hq.transcriptset.xml')  # HQ transcriptset
        self.transcript_ds = f('pbcoretools.tasks.gather_transcripts-1',
                               'file.transcriptset.xml')  # all transcriptset

    @property
    def flnc_bam(self):
        f0 = op.join(self.root_dir, 'tasks', 'isoseq3.tasks.cluster-0', 'unpolished.flnc.bam')
        f1 = op.join(self.root_dir, 'tasks', 'isoseq3.tasks.refine-0', 'flnc.bam')
        return f0 if op.exists(f0) else f1 if op.exists(f1) else None

    def export_unpolished_fa(self, unpolished_fa):
        _bam2fastx('bam2fasta', self.unpolished_bam,  unpolished_fa)
