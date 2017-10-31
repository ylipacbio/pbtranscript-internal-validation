#!/usr/bin/env python

import sys
import os.path as op
from .io import SMRTLinkIsoSeq2Files
from .validate_smrtlink_isoseq_rc0 import run, get_parser
from pbtranscript2.mains.post_mapping_to_genome import post_mapping_to_genome_runner


def main(args=sys.argv[1:]):
    """main"""
    def post_mapping_to_genome_runner_f(in_isoforms, in_sam, in_pickle,
            out_isoforms, out_gff, out_abundance, out_group, out_read_stat,
            min_aln_coverage, min_aln_identity, min_flnc_coverage,
            max_fuzzy_junction, allow_extra_5exon, min_count):
        return post_mapping_to_genome_runner(in_isoforms=in_isoforms, in_sam=in_sam,
                in_json=in_pickle, out_isoforms=out_isoforms, out_gff=out_gff,
                out_abundance=out_abundance, out_group=out_group, out_read_stat=None,
                min_aln_coverage=min_aln_coverage, min_aln_identity=min_aln_identity,
                max_fuzzy_junction=max_fuzzy_junction, allow_extra_5exon=allow_extra_5exon,
                min_count=min_count)

    run(SMRTLinkIsoSeq2Files, post_mapping_to_genome_runner_f, get_parser().parse_args(args))

if __name__ == "__main__":
    main()
