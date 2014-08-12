from random import random
from xlrd import *
from common_functions import *
from tfbs_selector import *

def maine() :
    tf_sub(uas, strength, seq, parameterD)
    polyAT_sub(uas, strength, seq, parameterD)

def tf_sub(uas, strength, seq, parameterD) :

    sites = parameterD['tf_sites'][uas]

    count = 0

    tf_sub_recordD = {}
    
    while count < sites :
        
        # This ensures the first site is always inserted with a TF site.  For UAS2, the other
        # three are inserted with the probability in the ProGenie_Parameters.xlsx.  For UAS1,
        # there is only one other site that has the probability of insertion.
        if count < 1 :

            # Choose a transcription factor site to substitute
            tf_seq_list = tf_choose_sequence(uas, strength, parameterD)
            
            # Determine placement of the TF in the sequence
            # Returns list of the form ['left sequence cutoff','right sequence cutoff]
            tf_loc_list = tf_choose_location(uas, count, tf_seq_list[2], parameterD)
            
            # Write the modified sequence with the TF site inserted.
            seq_L = seq[:tf_loc_list[0]]
            tf_seq = tf_seq_list[2]
            seq_R = seq[tf_loc_list[1]:]
            
            seq = seq_L+tf_seq+seq_R
        
            tf_sub_recordD[count+1] = {'name' : tf_seq_list,
                                       'location' : tf_loc_list}

        if count >= 1:
            
            # This random generator decides whether to insert a TFBS sequence at subsequent sites, based on
            # the probability defined in ProGenie_Parameters.xlsx for each strength level.            
            tf_insert = random()
            tf_yes = parameterD['tf_add'][uas][strength]

            if tf_insert <= tf_yes :
                
                # Choose a transcription factor site to substitute
                tf_seq_list = tf_choose_sequence(uas, strength, parameterD)
                
                # Determine placement of the TF in the sequence
                # Returns list of the form ['left sequence cutoff','right sequence cutoff]
                tf_loc_list = tf_choose_location(uas, count, tf_seq_list[2], parameterD)

                # Write the modified sequence with the TF site inserted.
                seq_L = seq[:tf_loc_list[0]]
                tf_seq = tf_seq_list[2]
                seq_R = seq[tf_loc_list[1]:]

                seq = seq_L+tf_seq+seq_R
                    
                tf_sub_recordD[count+1] = {'name' : tf_seq_list,
                                           'location' : tf_loc_list}
      
        count += 1

    return [seq, tf_sub_recordD]

def tf_choose_sequence(uas, strength, parameterD) :
 
    # Generate a random number to choose which TF to substitute, likelihood
    #of choosing a particular transcription factor is determined in the tf_probs list.
    #(0): REB1, (0-1): RAP1, (1-2): GCR1, (2-3): ABF1,(3-4): MCM1, (4): RSC3.
    # In UAS2, to limit the number of patterns, as well as eliminate those TFs that act
    # close to the transcription start site, only allow REB1, RAP1, and GCR1 sites.
    tf_choose = random()
            
    tf_probs = parameterD['tf_choice_prob'][uas]

    if tf_choose <= tf_probs[0] :
        tf_list = reb1_site_chooser(strength, parameterD)
        
    if tf_probs[0] < tf_choose <= tf_probs[1] :
        tf_list = rap1_site_chooser(strength, parameterD)
        
    if tf_probs[1] < tf_choose <= tf_probs[2] :
        tf_list = gcr1_site_chooser(strength, parameterD)
        
    if tf_probs[2] < tf_choose <= tf_probs[3] :
        tf_list = abf1_site_chooser(strength, parameterD)
        
    if tf_probs[3] < tf_choose <= tf_probs[4] :
        tf_list = mcm1_site_chooser(strength, parameterD)
        
    if tf_probs[4] < tf_choose <= 1 :
        tf_list = rsc3_site(strength, parameterD)

    return tf_list

def tf_choose_location(uas, count, tf, parameterD) :
    
    location = parameterD['tf_loci'][uas][count]

    # In UAS2, the polyAT will be before the TFBS.  In UAS1, the TFBS will be before the polyAT site.
    # The purpose of this is to wrap the TFBS region in nucleosome-disfavoring sequence and also to
    # separate the TFBS region further from the core TSS, in the manner of native promoters.
    tfsliceD = {1 : [location-len(tf),
                     location],
                
                2 : [location,
                     location+len(tf)]}

    location_list = tfsliceD[uas]
    
    return location_list

def polyAT_sub(uas, strength, seq, parameterD):

    sites = parameterD['pAT_sites'][uas]
    
    count = 0

    pAT_sub_recordD = {}

    while count < sites :

        # This random generator decides whether to insert a
        # poly dA:dT sequence, based on the probability defined
        # in the polyatD dictionary for each strength level.
        pAT_insert = random()

        pAT_yes = parameterD['pAT_add'][uas][strength]
        
        if pAT_insert <= pAT_yes :

            # Choose a sequence to insert
            pAT_seq_list = pAT_choose_sequence(uas, strength, parameterD)

            # Determine placement of the poly dA:dT
            pAT_loc_list = pAT_choose_location(uas, count, parameterD)
  
            # Write the modified sequence with the poly dA:dT inserted:
            seq_L = seq[:pAT_loc_list[0]]
            pAT_seq = pAT_seq_list[1]
            seq_R = seq[pAT_loc_list[1]:]
            
            seq = seq_L+pAT_seq+seq_R

            pAT_sub_recordD[count+1] = {'name': pAT_seq_list,
                                        'location': pAT_loc_list}

        count += 1
    
    return [seq, pAT_sub_recordD]

def pAT_choose_sequence(uas, strength, parameterD) :
    
        # This random generator decides which of the polydA:dT sequences to insert.  It will choose either
        # (0): polyT, (1): polyA, (2): random mix.  As polyA will likely decrease strength of the promoter
        # Or at least increase strength less than T or a mix), make polyA less likely as the strength increases.
        pAT_choose = random()

        pAT_probs = parameterD['pAT_choice_prob'][strength]

        index = 0
            
        if pAT_choose <= pAT_probs[0] :
            index = 1
            
        if pAT_probs[0] < pAT_choose <= pAT_probs[1] :
            index = 2

        pAT_list = parameterD['pAT_seqs'][index]

        return pAT_list

def pAT_choose_location(uas, count, parameterD) :
    
        # I want to insert the poly dA:dT before the probable TFBS location, so the slice location chosen is
        # the endslice for the polyAT.  It will be the startslice for the TFBS.
        location = parameterD['pAT_loci'][uas][count]
        
        pAT_length = 13
            
        psliceD = {1 : [location, location+pAT_length],
                       
                   2 : [location-pAT_length, location]}

        location_list = psliceD[uas]

        return location_list
        
if __name__ == "__maine__" :
    maine()
