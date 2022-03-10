from Bio import SeqIO
import tqdm
import os
import numpy as np

from common import generator
from common import utilsBlast

##########################################################################

def get_protein_info(fasta_file):
    '''
    Function that get the protein length of all proteins
    and put it in a dictionnary

    :param fasta_file: path the the main fasta file used to blast all vs all
    :type: list of string
    :return: a dataframe with the information about the proteins and their cluster, a dictionnary with the protein and their length
    :rtype: pandas.DataFrame, dict 
    '''


    print('\n#######################')
    print('# Parsing cluster files')
    print('#######################\n')

    dict_prot = {}

    numline = generator.buf_count_prot_gen(fasta_file)
    
    parser = SeqIO.parse(fasta_file, 'fasta')

    for protein in tqdm.tqdm(parser, total=numline):
        dict_prot[protein.id] = len(protein.seq)

    return dict_prot

##########################################################################

def get_cluster_info(fasta_clusters):
    '''
    Function that parse the fasta cluster file and create a dataframe 
    that associate protein and cluster. 

    :param fasta_cluster: list of path to fasta files
    :type: list of string
    :return: a dictionnary with the information about the proteins and their cluster
    :rtype: dict
    '''


    print('\n#######################')
    print('# Parsing cluster files')
    print('#######################\n')

    dict_fam = {}

    for fasta in tqdm.tqdm(fasta_clusters):
        # infer family fron cluser file
        family = os.path.split(fasta)[-1].split('.')[0]

        parser = SeqIO.parse(fasta, 'fasta')

        for protein in parser:
            dict_fam[protein.id] = family

    return dict_fam

##########################################################################


def create_table_threshold(blast_out, families, protein_dict, output, output_removed, length_treshold=100, option_cov='mean', option_pid='mean'):
    '''
    Function that take the blast all vs all file and create a dataframe 
    that will summarized the information about the different pairs of hits

    :param blast_out: Path to the blast all vs all file
    :type: string
    :param families: Dictionnary that contains all the informations about which cluster the sequence is in (key: protein_id, value:family)
    :type: dict
    :param protein_dict: Dictionnary with the protein sequence in key and lenght in value
    :type: dict
    :param output: path to a file that summarise all the informations about the sequence from the blast, family (columns: 'protein1', 'protein2', 'pident', 'evalue', 'coverage', 'fam')
    :type: string
    :param output_removed: path to a file a file that summarise if by mistake protein are prensent in the blast anad not in the fasta
    :type: string
    :param length_threshold: Minimum length to accept hsp (default: 100)
    :type: int
    :param option_cov: Option to calculate the coverage: ["mean","subject", "query", "shortest", "longest"] (default: 'mean')
    :type: string
    :param option_id: Option to calculate the p[ercentage of identity: ["mean","subject", "query", "shortest", "longest", "HSP"] (default: 'mean')
    :type: string    
    :return: Nothing
    '''

    print('\n##########################')
    print('# Creating table threshold')
    print('##########################\n')

    # Opening blast_out and preparation
    blast_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                'sstart', 'send', 'evalue', 'bitscore']

    # Get the types of te columns for multiple HSPs dataframe
    blast_dtypes = [('qseqid','S100'),
                    ('sseqid','S100'),
                    ('pident',np.float64),
                    ('length',np.int32),
                    ('mismatch',np.int32),
                    ('gapopen',np.int32),
                    ('qstart',np.int32),
                    ('qend',np.int32),
                    ('sstart',np.int32),
                    ('send',np.int32),
                    ('evalue',np.float64),
                    ('bitscore',np.float64),
                ]


    # Header of the output
    final_header = ['protein1', 'protein2', 'pident', 'evalue', 'coverage', 'fam']

    not_in_fasta = []

    with open(output, 'wt') as w_file:
        # Write the header in two times because format string need that
        header = '\t'.join(final_header)
        w_file.write(f"{header}\n")

        numline = utilsBlast.count_hsps(blast_out=blast_out)

        # Read the blast hsp by hsp
        for sub_blast in tqdm.tqdm(utilsBlast.iterrator_on_blast_hsp(blast_out=blast_out), total=numline) :
            # Get the number of hsps
            num_HSPs = len(sub_blast)

            qseqid = sub_blast[0][0]
            sseqid = sub_blast[0][1]

            if num_HSPs == 1:
                pident_blast, coverage_blast, evalue_blast, score = utilsBlast.summarize_hit_only(split_line = sub_blast[0], 
                                                                                                length_query = protein_dict[qseqid], 
                                                                                                length_subject = protein_dict[sseqid],
                                                                                                option_cov = option_cov,
                                                                                                option_pid = option_pid)
            else:
                df_hsps = utilsBlast.prepare_df_hsps(list_hsps = sub_blast,
                                                    blast_dtypes = blast_dtypes, 
                                                    HSPMIN = length_treshold)

                if df_hsps.shape[0] == 1:
                    pident_blast, coverage_blast, evalue_blast, score = utilsBlast.summarize_hit_only(split_line = df_hsps[0], 
                                                                                                    length_query = protein_dict[qseqid], 
                                                                                                    length_subject = protein_dict[sseqid],
                                                                                                    option_cov = option_cov,
                                                                                                    option_pid = option_pid)
                else:
                    delta_lg, coverage_blast, pident_blast, evalue_blast, score = utilsBlast.summarize_hits(df_hsps = df_hsps, 
                                                                                                            length_query = protein_dict[qseqid], 
                                                                                                            length_subject = protein_dict[sseqid],
                                                                                                            option_cov = option_cov, 
                                                                                                            option_pid = option_pid)

            # Look if both proteins are in the family
            if qseqid in families and sseqid in families:
                # If exist put in the table because both are in the family
                line2write = f'{qseqid}\t{sseqid}\t{pident_blast}\t{evalue_blast}\t{coverage_blast}\tin_family_{families[qseqid]}\n'
                w_file.write(line2write)

            # Look if one protein is in the family
            elif qseqid in families or sseqid in families:
                # Get which is in the family
                focus_id = qseqid if qseqid in families else sseqid

                # If by mistake fasta is present in blast output but not in fasta file
                if (qseqid not in protein_dict) and (qseqid not in not_in_fasta):
                    not_in_fasta.append(qseqid)
                elif (sseqid not in protein_dict) and (sseqid not in not_in_fasta):
                    not_in_fasta.append(sseqid)
                elif (qseqid in protein_dict) and (sseqid in protein_dict):
                    # If exist put in the table because one of them is in the family
                    line2write = f'{qseqid}\t{sseqid}\t{pident_blast}\t{evalue_blast}\t{coverage_blast}\tout_family_{families[focus_id]}\n'
                    w_file.write(line2write)

    if not_in_fasta:
        with open(output_removed, 'wt') as w_file:
            for id in not_in_fasta:
                w_file.write(f'{id}\n')
                        
    return 

##########################################################################