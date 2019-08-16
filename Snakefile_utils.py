import os, glob
import pandas as pd
import numpy as np
import logging

##### Function for loading all R1 and R2 fastqs
def get_all_fqgz(wildcards):
    '''
    Return all samples containing read 1 and read 2 under current folder
    '''
    r1 = glob.glob('{sample}/*R1*.fastq.gz'.format(sample=wildcards.all_samples))
    r2 = glob.glob('{sample}/*R2*.fastq.gz'.format(sample=wildcards.all_samples))
    assert len(r1) == 1, '{sample}'.format(wildcards.all_samples)
    assert len(r2) == 1
    return [r1[0], r2[0]]

def get_all_unmapped(wildcards):
    '''
    Return all unmapped samples containing read 1 and read 2 under current
    folder
    '''
    r1 = glob.glob('{sample}/*mate1'.format(sample=wildcards.all_samples))
    r2 = glob.glob('{sample}/*mate2'.format(sample=wildcards.all_samples))
    assert len(r1) == 1, '{sample}'.format(wildcards.all_samples)
    assert len(r2) == 1
    return [r1[0], r2[0]]

##### Function for loading all bams
def get_all_bams(wildcards):
    '''
    Return all samples containing bams under current folder
    '''
    bams = glob.glob('{sample}/*/*.bam'.format(sample=wildcards.all_samples))
    return [bams]

##### Functions for getting file names and unzipping files
def get_all_files(d):
    '''
    return all files in a directory
    '''
    return [d+"/"+f for f in os.listdir(d) if os.path.isfile(os.path.join(d, f))]

def unzip_fastq(f):
    '''
    function to gunzip .gz
    '''
    cmd = "gunzip " + f
    p = subprocess.Popen(cmd, shell=True)
    return None

## Functions for merging tables
def merge_htseq_tables(matches, outfile):
    '''
    merge a list of htseq tables defined by the list matches. Saves the merged
    output as outfile.
    '''
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
            if len(my_counts) > 0:
                samples.append(sample)
                counts.append(my_counts)
            else:
                logging.warning('{:} is empty'.format(match))

# Print output
    with open(outfile, 'w') as out:
        header = ["symbol"] + samples # header
        out.write("\t".join(header) + "\n")
        for i in range(len(genes)):
            my_counts = [counts[x][i] for x in range(len(counts))]
            line = "\t".join([genes[i]] + my_counts)
            out.write(line + "\n")

def merge_star_tables(matches, outfile):
    '''
    merge a list of star tables defined by the list matches. Saves the merged
    output as outfile.
    '''
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

## Preprocessing of expression matrix

def load_dataframes(input_list, axis=1, remove_last_nrows=0,
                    strip_header=True, sample_prefix=[], transpose=False):
    '''
    Load expression matrix as pandas dataframe. Typically we remove the last
    five rows from htseq output as they contain counting parameters.
    '''
    df = pd.DataFrame()
    if type(input_list) is str:
        input_list = [input_list]
    if len(sample_prefix) > 0:
        assert len(sample_prefix) == len(input_list), 'Number of sample prefix\
                must match the number of provided input dirs'

    for ind, f in enumerate(input_list):
        if ind == 0:
            df = pd.read_csv(f, index_col=0, sep='\t')
            if transpose == True:
                df = df.T
            if len(sample_prefix) > 0:
                df.columns = [sample_prefix[ind] + '_' + x for x \
                              in df.columns.values]
        else:
            tmp = pd.read_csv(f, index_col=0, sep='\t')
            if transpose == True:
                tmp = tmp.T
            if len(sample_prefix) > 0:
                tmp.columns = [sample_prefix[ind] + '_' + x for x \
                              in tmp.columns.values]
            df = pd.concat([df, tmp], axis=axis)
    if remove_last_nrows > 0:
        df = df.iloc[:-remove_last_nrows, :]
    if strip_header == True:
        df.index = [x.replace('transcript:', '').replace('gene:', '') \
                    if 'transgene:' not in x \
                    else x for x in df.index.values]
    return df

def convert_transcript_to_gene(infile, transcript_annotation, outfile,
                       remove_last_nrows=5):
    '''
    Convert transcript expression matrix to gene expression matrix
    '''
    X = load_dataframes(infile, remove_last_nrows=remove_last_nrows)
    ann = pd.read_csv(transcript_annotation, sep='\t', index_col=0)
    Y = convert_id_to_gene(X.T, ann, start='transcript', end='gene').T
    Y.to_csv(outfile, sep='\t', compression='gzip')


def convert_id_to_gene(X, gene_df, start='transcript', end='gene'):
    '''
    Convert feature_id expression matrix to transcript expression matrix
    '''
    dictionary = {}
    for key, value in zip(gene_df[start].values, gene_df[end].values):
        if type(key) is str:
            if value == 'Unknown':
                dictionary[key] = key
            else:
                dictionary[key] = value
    out = X.copy()
    new_names = [dictionary[x] if x in dictionary.keys() else x \
                 for x in out.columns.values]
    out.columns = new_names
    out = group_sum(out, rows=False)
    return out

## Experimental: separate inputs into N chunks
def chunks(l, n=40, prefix=""):
    out = []
    chunkid = []
    cid = 1
    for i in range(0, len(l), n):
        out.append(l[i:i+n])
        chunkid.append(os.path.join(prefix, 'group_'+str(cid)))
        cid += 1
    return (out, chunkid)

def write_chunks(d, postfix=''):
    for key in d:
        chunk_file = str(key) + postfix
        if os.path.exists(chunk_file):
            os.remove(chunk_file)
        with open(chunk_file, 'w') as f:
            for pairs in d[chunk_file]:
                f.write("{}\t".format(pairs[0]))
                f.write("{}\n".format(pairs[1]))

def load_chunks(infile, postfix=''):
    load = []
    chunk_file = infile + postfix
    with open(chunk_file, 'r') as f:
        for line in f:
            load.append(line.strip().split('\t'))
    return load

## Scripts to make dataframe mapping exon/transcript to gene and products
def transcript_annotation_name(infile):
    '''
    Figure out transcript annotation name
    '''
    filename = infile.split('/')[-1].split('.')[0]
    dirname = '/'.join(infile.split('/')[:-1])
    return os.path.join(dirname, filename)

def make_annotation_dataframe(infile, outfile):
    '''
    Create a dataframe that stores exon/transript names with their
    corresponding transcript/gene/product names.
    '''
    if outfile.split('.')[-1] != 'gz':
        outfile += '.gz'
    outpath = os.path.abspath(outfile)
    outpath = '/'.join(outpath.split('/')[:-1])
    if '.' in outfile.split('/')[-1]:
        outname_exon = outpath+'/'+outfile.split('/')[-1].split('.')[0] + \
               '_exon_annotation.tsv.gz'
        outname_transcript = outpath+'/'+outfile.split('/')[-1].split('.')[0]\
                +'_transcript_annotation.tsv.gz'
    else:
        outname_exon = outpath+'/exon_annotation.tsv.gz'
        outname_transcript = outpath+'/transcript_annotation.tsv.gz'
    print('Saving to {:}'.format(outname_exon))
    print('Saving to {:}'.format(outname_transcript))
    t_interest_type = ['transcript', 'mRNA', 'tRNA', 'ncRNA', 'rRNA',
                       'pseudogenic_transcript']
    g_interest_type = ['gene', 'pseudogene', 'ncRNA_gene']
    special_transcripts = ['ERCC-', 'transgene:']
    property_keywords = ['description=', 'Parent=', 'ID=']
    null_word = 'Unknown'

    def check_if_transcript(X, feature_type):
        '''
        Check if the current feature is a transcript type. The feature is
        identified as transcript if feature name contains "transcript:", or
        "ENMUST" (mouse ensembl transcript prefix), or if the labeled feature
        type is in our list of transcript feature types.
        '''
        is_transcript = False
        if X.split(':')[0] == 'transcript':
            is_transcript = True
        elif 'ENSMUST' in X:
            is_transcript = True
        else:
            for u in t_interest_type:
                if feature_type == u:
                    is_transcript = True
        return is_transcript

    def check_if_gene(X, feature_type):
        '''
        Check if the current feature is a special type. The feature is
        identified as transcript if feature name contains "gene:", or
        "ENMUSG" (mouse ensembl gene prefix), or if the labeled feature
        type is in our list of gene feature types.
        '''
        is_gene = False
        if X.split(':')[0] == 'gene':
            is_gene = True
        elif 'ENSMUSG' in X:
            is_gene = True
        else:
            for u in g_interest_type:
                if feature_type == u:
                    is_gene = True
        return is_gene

    def check_if_special(X):
        '''
        Check if the feature is a special type. If it is, then its feature,
        transcript, gene, and product names are all the same.
        '''
        is_special = False
        for u in special_transcripts:
            if u in X:
                is_special = True
        return is_special

    id_length = {}
    id_to_transcript = {}
    transcript_to_gene = {}
    all_t_to_gene = {}
    id_desc = {}
    parent_type = {}
    all_id = []
    all_transcript = []
    with open(infile, 'r') as f:
        for line in f:
            if line[0] != '#':
                split = line.rstrip().split('\t')
                if len(split) > 8:
                    feature_type = split[2]
                    length = int(split[4]) - int(split[3])
                    strand = split[6]
                    description, parent, ID = null_word, null_word, null_word
                    for u in split[8].split(';'):
                        for i, keyword in enumerate(property_keywords):
                            if keyword in u and i == 0:
                                description = u.replace(keyword, '').\
                                        split('[Source:')[0].replace('%2C ', ' ')
                            elif keyword in u and i == 1:
                                parent = u.replace(keyword, '')
                            elif keyword in u and i == 2:
                                ID = u.replace(keyword, '')
                    if ID != null_word and feature_type != 'CDS':
                        if ID not in parent_type.keys():
                            parent_type[ID] = feature_type
                        if ID not in id_length.keys():
                            id_length[ID] = length
                        if check_if_gene(ID, feature_type) == True:
                            iam = 'gene'
                        elif check_if_transcript(ID, feature_type) == True:
                            iam = 'transcript'
                        else:
                            iam = 'id'
                        if parent != null_word:
                            if check_if_gene(parent, parent_type[parent]) == True:
                                mom = 'gene'
                            else:
                                mom = 'transcript'
                        if iam == 'id' or iam == 'transcript':
                            is_special = check_if_special(ID)
                            if iam == 'id':
                                all_id.append(ID)
                                if is_special is True:
                                    all_transcript.append(ID)
                            elif iam == 'transcript':
                                all_transcript.append(ID)
                            if mom == 'transcript':
                                if ID not in id_to_transcript.keys():
                                    id_to_transcript[ID] = [parent]
                                elif parent not in id_to_transcript[ID]:
                                    id_to_transcript[ID].append(parent)
                            elif mom == 'gene':
                                if ID not in transcript_to_gene.keys():
                                    transcript_to_gene[ID] = [parent]
                                elif parent not in transcript_to_gene[ID]:
                                    transcript_to_gene[ID].append(parent)
                                if iam == 'transcript' and \
                                   ID not in all_t_to_gene.keys():
                                    all_t_to_gene[ID] = [parent]
                        else:
                            if ID not in id_desc.keys():
                                if description == null_word:
                                    id_desc[ID] = 'hypothetical protein'
                                else:
                                    id_desc[ID] = description
    df = pd.DataFrame(all_id, columns=['id'], dtype='object')
    df['transcript'] = [id_to_transcript[x][0] for x in all_id]
    all_t = [x for x in df['transcript'].values]
    df['gene'] = [transcript_to_gene[x][0] if x in transcript_to_gene.keys() else \
                  'Unknown' for x in all_t]
    all_g = [x for x in df['gene'].values]
    df['product'] = [id_desc[x] if x in id_desc.keys() else \
                     'Unknown' for x in all_g]
    df['length'] = [id_length[x] for x in df['id'].values]
    df['gene'] = [x if 'ERCC-' not in y and 'transgene:' not in y else y\
                  for x, y in zip(df['gene'].values, df['id'].values)]
    df['product'] = [x if 'ERCC-' not in y and 'transgene:' not in y else y\
                  for x, y in zip(df['product'].values, df['id'].values)]
#### clean up the transcript / gene names
    df['transcript']= [x.replace('transcript:','') for x in \
                        df['transcript'].values]
    df['gene']= [x.replace('gene:','') for x in \
                        df['gene'].values]
    df.to_csv(outname_exon, sep='\t', compression='gzip')
#### we make a dictionary just for transcripts. Turns out there are many
#### transcripts that all map to the same exons. So it's better to do feature
#### counting in transcript space and convert transcript to genes.
    df = pd.DataFrame(all_transcript, columns=['transcript'], dtype='object')
    df['gene'] = [all_t_to_gene[x][0] if x in all_t_to_gene.keys() else \
                  'Unknown' for x in all_transcript]
    df['product'] = [id_desc[x] if x in id_desc.keys() else \
                     'Unknown' for x in df['gene'].values]
    df['transcript']= [x.replace('transcript:','') for x in \
                        df['transcript'].values]
    df['gene']= [x.replace('gene:','') for x in df['gene'].values]
    df['gene'] = [x if 'ERCC-' not in y and 'transgene:' not in y else y\
                  for x, y in zip(df['gene'].values, df['transcript'].values)]
    df['product'] = [x if 'ERCC-' not in y and 'transgene:' not in y else y\
                  for x, y in zip(df['product'].values, df['transcript'].values)]
    df.to_csv(outname_transcript, sep='\t', compression='gzip')

def group_sum(input_df, rows=False):
    '''
    Pandas groupby sum function is way too slow. So I came up with an
    alternative.
    '''
    if rows:
        input_df = input_df.copy().T
    names, keys = np.unique(input_df.columns.values.astype(str),
                            return_inverse=True)
    n_keys = max(keys) + 1
    result = np.zeros((np.shape(input_df)[0], n_keys))
    gene_names = input_df.columns.values
    new_cols = []
    for i, new_name in zip(range(n_keys), names):
        result[:, i] = np.sum(input_df.loc[:, keys == i], axis=1)
        new_cols.append(new_name)
    result = pd.DataFrame(result, columns=new_cols, index=input_df.index)
    if rows:
        result = result.T
    return result

