import os, glob

##### Function for loading files
def get_all_fqgz(wildcards):
    """ Return all samples containing read 1 and read 2 under current folder """
    r1 = glob.glob('{sample}/*R1*.fastq.gz'.format(sample=wildcards.all_samples))
    r2 = glob.glob('{sample}/*R2*.fastq.gz'.format(sample=wildcards.all_samples))
    assert len(r1) == 1, '{sample}'.format(wildcards.all_samples)
    assert len(r2) == 1
    return [r1[0], r2[0]]

# ##### Functions for getting file names and unzipping files
def get_all_files(d):
    return [d+"/"+f for f in os.listdir(d) if os.path.isfile(os.path.join(d, f))]

def unzip_fastq(f):
    cmd = "gunzip " + f
    p = subprocess.Popen(cmd, shell=True)
    return None

## Functions for merging tables
def merge_htseq_tables(matches, outfile):
# Print directory containing each file
    print(matches)
    print(outfile)
    print("%s files" %(len(matches)))
# Get gene list from first column of first file (assumes that gene lists are the same for all files)
    genes = []
    with open(matches[0]) as f:
        for line in f:
            gene = line.rstrip().split("\t")[0]
            genes.append(gene)
# Get counts from each file
    samples = []
    counts = []
    for match in matches:
        if match.split("/")[-2] == 'htseq_output':
            sample = match.split("/")[-3]
        else:
            sample = match.split("/")[-2]
        with open(match) as f:
            my_counts = []
            for line in f:
                count = line.rstrip().split("\t")[1]
                my_counts.append(count)
        samples.append(sample)
        counts.append(my_counts)
# Print output
    with open(outfile, 'w') as out:
        header = ["symbol"] + samples # header
        out.write("\t".join(header) + "\n")
        for i in range(len(genes)):
            my_counts = [counts[x][i] for x in range(len(counts))]
            line = "\t".join([genes[i]] + my_counts)
            out.write(line + "\n")

def merge_star_tables(matches, outfile):
# Print directory containing each file
    print(matches)
    print(outfile)
    print ("%s files" %(len(matches)))
# Parse statistics from each file
    samples = []
    num_input_reads = []
    num_input_read_length = []
    num_uniquely_mapped_reads = []
    num_reads_multiple_loci = []
    num_reads_too_many_loci = []
    percent_reads_unmapped_too_many_mismatches = []
    percent_reads_unmapped_too_short = []
    percent_reads_unmapped_other = []
    percent_reads_unmapped_too_many_loci = []
    for match in matches:
        if match.split("/")[-2] == 'STAR_output':
            sample = match.split("/")[-3]
        else:
            sample = match.split("/")[-2]
        samples.append(sample)
        with open(match) as f:
            for line in f:
                if "Number of input reads" in line:
                    x = line.rstrip().split()[-1]
                    num_input_reads.append(x)
                elif "Average input read length" in line:
                    x = line.rstrip().split()[-1]
                    num_input_read_length.append(x)
                elif "Uniquely mapped reads number" in line:
                    x = line.rstrip().split()[-1]
                    num_uniquely_mapped_reads.append(x)
                elif "Number of reads mapped to multiple loci" in line:
                    x = line.rstrip().split()[-1]
                    num_reads_multiple_loci.append(x)
                elif "Number of reads mapped to too many loci" in line:
                    x = line.rstrip().split()[-1]
                    num_reads_too_many_loci.append(x)
                elif "% of reads unmapped: too many mismatches" in line:
                    x = line.rstrip().split()[-1]
                    percent_reads_unmapped_too_many_mismatches.append(x)
                elif "% of reads unmapped: too short" in line:
                    x = line.rstrip().split()[-1]
                    percent_reads_unmapped_too_short.append(x)
                elif "% of reads unmapped: other" in line:
                    x = line.rstrip().split()[-1]
                    percent_reads_unmapped_other.append(x)
                elif "% of reads mapped to too many loci" in line:
                    x = line.rstrip().split()[-1]
                    percent_reads_unmapped_too_many_loci.append(x)
# Write output
    features = ["input",
                "read_length",
                "uniquely_mapped",
                "multiple_loci",
                "too_many_loci",
                "percent_unmapped_too_many_mismatches",
                "percent_unmapped_too_short",
                "percent_unmapped_other",
                "percent_unmapped_too_many_loci"]
    stats = [num_input_reads, num_input_read_length, num_uniquely_mapped_reads,
             num_reads_multiple_loci, num_reads_too_many_loci,
             percent_reads_unmapped_too_many_mismatches,
             percent_reads_unmapped_too_short, percent_reads_unmapped_other,
             percent_reads_unmapped_too_many_loci]
    with open(outfile, 'w') as out:
        header = ["name"] + features # header
        out.write("\t".join(header) + "\n")
        for i in range(len(samples)):
            sample = samples[i]
            my_stats = [stats[x][i] for x in range(len(stats))]
            line = "\t".join([sample] + my_stats)
            out.write(line + "\n")

## Experimental: separate inputs into N chunks
def chunks(l, n=40):
    out = []
    chunkid = []
    cid = 1
    for i in range(0, len(l), n):
        out.append(l[i:i+n])
        chunkid.append('group_'+str(cid)+'.in')
        cid += 1
    return (out, chunkid)

def write_chunks(wildcards):
    for g, c in zip(wildcards.groups, wildcards.chunkid):
        with open(os.path.join(wildcards.root_dir, c), 'w') as f:
            for sample in g:
                f.write('%s\n' %sample)

def load_chunks(wildcards):
    pass
