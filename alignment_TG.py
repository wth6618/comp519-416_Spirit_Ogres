import mappy as mp
import timeit
import subprocess
import os
from My_Class import *
from util import *
import numpy as np
from Bio import SeqIO
from Pre_process import preprocessing
import distance


output_dict = dict()
"""
Inputs: reference sequence
Outputs: Aligner
"""
def construct_aligner(input_seq):
    a = mp.Aligner(seq=input_seq, preset='map-ont')
    # a = mp.Aligner("tryB_transcripts.fasta",preset='map-ont')  # load or build index
    if not a: raise Exception("ERROR: failed to load/build index")

    return a

def get_children(id,DB,type):

    return DB.children(id, featuretype=type)


"""
Inputs: source database; target database; matching window, transcripts within the window
Output: classified calls and category

Candidates with accuracy above the first threshold will be considered as potential 
unique_transcript matches. 
 In the mRNA-level, check the number of transcript and marginal position in the matching window 
 In the exon-level, check the number of exons and the marginal position of the first and last exon.
"""
def check_features(name_S,f_arrary,low_boundary, upper_boundary,FeatureDB_T,FeatureDB_S,count,length_S):

    print(len(f_arrary))
    if len(f_arrary) == 1:

        if f_arrary[0].start> low_boundary and f_arrary[0].end < upper_boundary:

            print(f_arrary[0].id)
            #parent = FeatureDB.parents(f_arrary[0].id,featuretype = 'gene')
            return 'unique_transcript', 'exact_match',f_arrary[0].id
        children_T = list(get_children(f_arrary[0].id, FeatureDB_T, 'exon'))
        children_S = list(get_children(name_S, FeatureDB_S, 'exon'))
        # check if they have same number of exons
        if len(list(children_T)) == len(list(children_S)) :
            if (children_T[0].end <low_boundary) or (children_T[-1].start >upper_boundary):
                if count > 0:
                    return '', '', ''
                else:
                    parent = list(FeatureDB_T.parents(f_arrary[0].id, featuretype='gene'))
                    return 'absent_transcript', 'new_exons', parent[0].id
            else:
                print(f_arrary[0].id)

                return 'unique_transcript', 'all_jxn_match', f_arrary[0].id

        # check exon boundary
        else:

            if len(list(children_T)) > len(list(children_S)):
                window = upper_boundary - low_boundary + 0.10* length_S

                # if the number of exons in the match windows equals to the target exons
                if length_S <= window:
                    print(f_arrary[0].id)

                    return 'unique_transcript', 'source_contained', f_arrary[0].id
                else:
                    if count > 0:
                        return '', '', ''
                    else:
                        parent = list(FeatureDB_T.parents(f_arrary[0].id, featuretype='gene'))
                        return 'absent_transcript', '', parent[0].id

            else:
                region = (f_arrary[0].seqid, low_boundary, upper_boundary)

                exon_features = list(FeatureDB_T.region(region, featuretype='exon'))


                # if the match window contains the entire source transcript
                if len(exon_features) == len(children_T):

                    print(f_arrary[0].id)
                    return 'unique_transcript', 'target_contained', f_arrary[0].id
                else:
                    if count > 0:
                        return '', '', ''
                    else:
                        parent = list(FeatureDB_T.parents(f_arrary[0].id, featuretype='gene'))
                        return 'absent_transcript', '', parent[0].id
    else:
        gf_set = []
        window = upper_boundary - low_boundary +30

        if count == 0 and window >= length_S * 0.25:
        #if count == 0:
            same_gene = True
            for f in f_arrary:
                for gene in FeatureDB_T.parents(f.id,featuretype='gene'):

                    if len(gf_set) == 0:
                        gf_set.append(gene.id)
                        continue

                    # if empty means no same parent
                    elif gene.id not in gf_set:
                        gf_set.append(gene.id)
                        same_gene = False

            if same_gene:
                print('same gene')
                parent = list(FeatureDB_T.parents(f_arrary[0].id, featuretype='gene'))
                return 'absent_transcript', 'changed_exons', parent[0].id
            else:
                target_genes = ';'.join(gf_set)
                print(target_genes)
                return 'gene_fusion','new_exons',target_genes
        else:
            parent = list(FeatureDB_T.parents(f_arrary[0].id, featuretype='gene'))
            return 'absent_transcript', '', parent[0].id


"""
Inputs:
all the exons within the matching window;matching window boundaries; 

Output:
return the type of absent and corresponding catergory 

"""
def check_absent(exon_arrary,low_boundary, upper_boundary,target_gene):

    print('check_absent')


    # no exon present
    if len(exon_arrary) == 0 :

        return 'absent_transcript','novel_retained_intron',target_gene

        # check exon boundary
    else:
        if low_boundary < exon_arrary[0].start and upper_boundary > exon_arrary[-1].end:

            return 'absent_transcript', 'changed_exons', target_gene
        else:
            return 'absent_transcript', '', target_gene

"""
Inputs: Match information; database; the output file handle
Output: None

base on the match information and write result to the output file.
the calls are written with priority: 
unique_transcript > multiple_transcript>gene_fusion>absent_transcript
> absent_gene > absent_genome

"""
def generate_output(matches,source_id,gene,sourceA,sourceB,DB_T,file):
    u_count = 0
    g_count = 0
    g_info = None
    u_info = []
    abt, abg,abgenome = False, False,False
    if len(gene) == 0:
        print('{} not in database'.format(source_id))
        output = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sourceA, sourceB,
                                                               source_id,
                                                               '', 'absent_gene', 0,
                                                               source_id, '', 'match_not_same_gene')
        file.write(output)
        return
    for M in matches:
        if M[0] == 'unique_transcript':
            u_count += 1
            u_info.append(M)
        if M[0] == 'gene_fusion':
            g_count += 1
            g_info = M
        if M[0] == 'absent_transcript':
            abt = True
            abt_info = M
        if M[0] == 'absent_gene':
            abg = True
            abg_info = M



    if u_count == 0:
        if g_count > 0:
            # write gene_fusion
            output = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sourceA, sourceB,
                                                                   source_id,
                                                                   '', 'gene_fusion', 0,
                                                                   gene[0].id, g_info[2], g_info[1])
            file.write(output)
            return
        if abt:
            output = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sourceA, sourceB,
                                                                   source_id,
                                                                   '', 'absent_transcript', 0,
                                                                   gene[0].id,abt_info[2] , abt_info[1])
            file.write(output)
            return
        if abg:

            output = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sourceA, sourceB,
                                                                   source_id,
                                                                   '', 'absent_gene', 0,
                                                                   gene[0].id, '',abg_info[1])
            file.write(output)
            return
        else:
            output = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sourceA, sourceB,
                                                                   source_id,
                                                                   '', 'absent_genome', 0,
                                                                   gene[0].id, '', 'unmapped')
            file.write(output)
            return

    else:
        if u_count == 1:
            call, category, target_id = u_info[0]
            geneT = list(DB_T.parents(target_id, featuretype='gene'))
            output = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sourceA, sourceB,
                                                                   source_id,
                                                                   target_id, 'unique_transcript', 0,
                                                                   gene[0].id, geneT[0].id, category)
            file.write(output)
            return

        # multiple transcript
        else:
            for u in u_info:
                call, category, target_id = u
                #print(target_id)
                geneT = list(DB_T.parents(target_id, featuretype='gene'))
                #if len(geneT) == 0:
                #    print('error')
                #if len(gene) == 0:
                #    print('gene_error')
                output = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sourceA, sourceB,
                                                                       source_id,target_id,
                                                                       'multiple_transcript', 0,
                                                                       gene[0].id,
                                                                       geneT[0].id, '')


                file.write(output)
            return



"""
Input: actual sequence of the gene in the matching window
Output: boolean whether it can be classified as antisense or not
"""
def check_antisense(gene,trans,start):

    print('antisense')

    h_dist = distance.hamming(gene[start:len(trans)+1],trans)
    threshold = len(trans)/10
    if h_dist < threshold:
        return True
    else:
        return False


# source B target A
"""
Inputs: source transcripts; target genome; output file handle; source/target database

output: None

Align every source transcript to the target genome and calculate accuracy from each hit;
use two accuracy thresholds to distinguish hit/match for further classification
"""
def Alignment_TG(FeatureDB_S, FeatureDB_T,SourceS,SourceT,write_file,m_file = None, source_file = None,target_file = None):

    #match_dict = np.load(m_file,allow_pickle=True).item()
    match_dict = dict()
    #target

    # source
    s_filehandle = open(source_file)
    trans_file = SeqIO.parse(s_filehandle, 'fasta')
    for trans in trans_file:

        trans_id, trans_seq = trans.id, str(trans.seq)
        print('mapping {} to'.format(trans_id))
        match_dict[trans_id] = []
        count = 0
        t_filehandle = open(target_file)
        chromo_file = SeqIO.parse(t_filehandle, 'fasta')
        for genome in chromo_file:
            seqid, genome_seq = genome.id, str(genome.seq)
            #print(seqid)
            aligner = construct_aligner(genome_seq)
            for hit in aligner.map(trans_seq):
                length = len(trans_seq)
                accuracy = 1.0 - float(hit.NM / length)
                start, end = hit.r_st - 20, hit.r_en + 20
                #print(accuracy)
                if accuracy > 0.90:
                    region = (seqid,start,end)
                    interv_features= list(FeatureDB_T.region(region,featuretype='mRNA'))
                    if len(interv_features) != 0:
                        call, category,target_t = check_features(trans_id,interv_features,start,end,FeatureDB_T,FeatureDB_S,count,len(trans_seq))
                    else:
                        if count == 0:
                            gene_features = list(FeatureDB_T.region(region, featuretype='gene'))

                            if len(gene_features) == 0:

                                call, category, target_t = 'absent_gene', 'match_not_same_gene', ''
                            else:
                                '''gene = gene_features[0]
                                gen_seq = gene.sequence(target_file, use_strand=True)
                                if check_antisense(gen_seq,trans_seq,start-gene.start):

                                    call, category, target_t = 'absent_gene', 'antisense', ''
                                else:
                                    call, category, target_t = 'absent_transcript', '', '' '''
                                call, category, target_t = 'absent_transcript', '', ''

                        else:
                            call, category, target_t = '', '', ''


                    if call == 'unique_transcript' or call == 'gene_fusion':
                        count += 1

                    if call != '':
                        match_dict[trans_id].append((call,category,target_t))

                # not reaching the accuracy threshold
                elif accuracy > 0.30:
                    if count == 0:

                        region = (seqid, start, end)

                        interv_features = list(FeatureDB_T.region(region, featuretype='mRNA'))
                        if len(interv_features) != 0:
                            parent = list(FeatureDB_T.parents(interv_features[0].id, featuretype='gene'))

                            exon_features = list(FeatureDB_T.region(region, featuretype='exon'))
                            call, category,target_t = check_absent(exon_features,start,end,parent[0].id)
                        else:
                            gene_features = list(FeatureDB_T.region(region, featuretype='gene'))
                            if len(gene_features) == 0:

                                call, category, target_t = 'absent_gene', 'match_not_same_gene', ''
                            else:

                                call, category, target_t = 'absent_transcript', '', ''
                    else:
                        call, category, target_t = '', '', ''

                    if call != '':
                        match_dict[trans_id].append((call,category,target_t))
                else:
                    call, category, target_t = 'absent_genome', 'unmmapped', ''
                    match_dict[trans_id].append((call, category, target_t))

        gene = list(FeatureDB_S.parents(trans_id,featuretype='gene'))
        generate_output(match_dict[trans_id],trans_id,gene,SourceS,SourceT,FeatureDB_T,write_file)

        t_filehandle.close()


    return match_dict



"""
Testing Session:

start = timeit.default_timer()
FeatureDB_A,FeatureDB_B,FeatureDB_C = preprocessing()
print('finish preprocess')
file_handle = open('Outputs/output_BA.txt','w+')
file_handle.write('#SourceA	SourceB	SourceA_Transcript_ID	SourceB_Transcript_ID	Call	Score	SourceA_Gene	SourceB_Gene	Category\n')
output = Alignment_TG(FeatureDB_B,FeatureDB_A,'B','A',file_handle,m_file='saved_data/BA.npy',source_file='tryB_transcripts.fasta',
                                                      target_file='Challenge_9934185_A.chromosomes/A.chromosomes.fasta')
np.save('saved_data/BA.npy',output,allow_pickle=True)

output = Alignment_TG(FeatureDB_B,FeatureDB_A,'B','A',file_handle,m_file='saved_data/BA.npy',source_file='tryB_transcripts.fasta',
                                                      target_file='Challenge_9934185_A.chromosomes/A.chromosomes.fasta')
np.save('saved_data/AB.npy',output,allow_pickle=True)
file_handle.close()
print('done')

stop = timeit.default_timer()
print('Time: ', stop - start)"""