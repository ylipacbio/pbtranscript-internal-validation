#!/usr/bin/env python

import os.path as op
from pbtranscript2.independent.workspace import WorkSpace
from pbtranscript2.independent.utils import realpath

class SMRTLinkIsoSeq2Files(object):
    def __init__(self, root_dir):
        self.root_dir = realpath(root_dir)
        def f(*args):
            return op.join(self.root_dir, *args)

        self.tasks_dir = f('tasks')
        self.ws_json = f('tasks', 'pbtranscript2tools.tasks.create_workspace-0', 'isoseq2_workspace.json')
        self.ws_obj = WorkSpace.from_json(self.ws_json)
        self.ws_obj.root_dir = f(op.dirname(self.ws_json), 'workspace')
        self.subreads_xml = self.ws_obj.subreads_xml
        self.isoseq_flnc_ds = f('tasks', 'pbcoretools.tasks.gather_contigset-3', 'file.contigset.xml')
        self.isoseq_nfl_ds = f('tasks', 'pbcoretools.tasks.gather_contigset-2', 'file.contigset.xml')
        self.isoseq_flnc_fa = f(self.ws_obj.root_dir, self.ws_obj.flnc_fa)
        self.isoseq_nfl_fa = f(self.ws_obj.root_dir, self.ws_obj.nfl_fa)
        self.ccs_ds = self.ws_obj.ccs_xml
        self.ccs_fa = f('tasks', 'pbcoretools.tasks.bam2fasta_ccs-0', 'ccs.fasta')
        self.ccs_fa_gz = f('tasks', 'pbcoretools.tasks.bam2fasta_ccs-0', 'ccs.tar.gz')
        self.hq_isoforms_fa = f(self.ws_obj.root_dir, self.ws_obj.hq_fa)
        self.hq_isoforms_fq = f(self.ws_obj.root_dir, self.ws_obj.hq_fq)
        self.lq_isoforms_fa = f(self.ws_obj.root_dir, self.ws_obj.lq_fa)
        self.lq_isoforms_fq = f(self.ws_obj.root_dir, self.ws_obj.lq_fq)
        self.consensus_isoforms_fa = f(self.ws_obj.root_dir, self.ws_obj.consensus_isoforms_fa)
        self.ccs_report_json = f('tasks', 'pbreports.tasks.ccs_report-0', 'ccs_report.json')
        self.classify_report_json = f('tasks', "pbreports.tasks.isoseq_classify-0", "isoseq_classify_report.json")
        self.cluster_report_json = f('tasks', "pbreports.tasks.isoseq_cluster-0", "isoseq_cluster_report.json")
        self.cluster_report_csv = f(self.ws_obj.root_dir, 'cluster_report.FL_nonFL.csv')
        self.sample_prefix_to_uc_pickle_json = f(self.ws_obj.root_dir, 'sample_to_uc_pickle.json')
        # pickle was used for pbtranscript while json was used for pbtranscript2
        self.hq_lq_prefix_pickle = self.sample_prefix_to_uc_pickle_json
