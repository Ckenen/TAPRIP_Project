#!/usr/bin/env python
import sys
# import os
# import shutil
import multiprocessing as mp
import numpy as np
import pyBigWig


def worker(bwfile, chrom, outfile=None):
    width = 1000
    num = 10000
    bulk = width * num
    fw = None
    rows = None
    if outfile is None:
        rows = []
    else:
        fw = open(outfile, "w+")
        
    with pyBigWig.open(bwfile) as f:
        size = f.chroms()[chrom]
        for start in range(0, size, bulk):
            end = min(start + bulk, size)
            covs = np.nan_to_num(f.values(chrom, start, end))
            for i1 in range(0, len(covs), width):
                i2 = min(i1 + width, len(covs))
                v = np.mean(covs[i1:i2])
                if v == 0:
                    continue
                row = [chrom, start + i1, start + i2, v]
                if fw is None:
                    rows.append(row)
                else:
                    fw.write("\t".join(map(str, row)) + "\n")
    if rows is None:
        return outfile
    else:
        return rows


def main():
    infile, threads, outfile = sys.argv[1:]
    threads = int(threads)
    
    # tempdir = outfile + ".TMP"
    # if not os.path.exists(tempdir):
    #     os.mkdir(tempdir)
    
    results = []
    pool = mp.Pool(threads)
    with pyBigWig.open(infile) as f:
        for chrom in f.chroms():
            # outfile2 = tempdir + "/%s.txt" % chrom
            outfile2 = None
            r = pool.apply_async(worker, (infile, chrom, outfile2))
            results.append(r)
    pool.close()
    pool.join()
    
    with open(outfile, "w+") as fw:
        for r in results:
            for row in r.get():
                line = "\t".join(map(str, row))
                fw.write(line + "\n")
                
    # shutil.rmtree(tempdir)
    
    
if __name__ == "__main__":
    main()