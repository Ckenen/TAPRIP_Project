#!/usr/bin/env python
import sys
from collections import defaultdict
import pysam


inbam, outbam = sys.argv[1:]

with pysam.AlignmentFile(inbam) as f:
    with pysam.AlignmentFile(outbam, "wb", f) as fw:
        data = defaultdict(list)
        for s in f:
            if s.is_proper_pair:
                continue
            if s.is_secondary:
                continue
            if not s.mate_is_unmapped and s.next_reference_id != s.reference_id:
                continue
            data[(s.reference_name, s.query_name, s.is_read1)].append(s)
        ss = []
        for k, v in data.items():
            if len(v) == 2:
                s1, s2 = v
                if s1.is_reverse == s2.is_reverse:
                    ss.append(s1)
                    ss.append(s2)
        for s in ss:
            fw.write(s)