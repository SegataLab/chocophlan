#!/usr/bin/env python3

__author__ = ('Fabio Cumbo (fabio.cumbo@unitn.it)')
__version__ = '0.01'
__date__ = 'Jan 21, 2019'

import sys, getopt, os, time

def load_association_file( filepath, retrieveOnly='', fieldSep='\t' ):
    associations = set()
    with open( filepath ) as association_file:
        for line in association_file:
            if line.strip():
                if line.startswith( retrieveOnly ):
                    associations.add( tuple( line.strip().split( fieldSep )[1:] ) ) # skip th reference genome label
    return associations

# index 0 = SGB id
# index 1 = GGB id
# index 2 = FGB id
def get_id( associations, searchFor, searchIdx, retrieveIdx ):
    for association in associations:
        if association[ searchIdx ] == searchFor:
            return association[ retrieveIdx ]
    return None

def build_description_line( fields, fieldSep=';', newLine=True ):
    description_line = fieldSep.join( field for field in fields )
    if newLine:
        return description_line + '\n'
    return description_line

def create_description_file( descriptions_filepath, associations, output_filepath, sgb_idx=0, levelEstTax_idx=4, estTax_idx=5, closestTax_idx=7, fieldSep=";" ):
    out_file = open( output_filepath, 'w+' )
    with open( descriptions_filepath ) as descriptions_file:
        skipLine = True
        for line in descriptions_file:
            if line.strip():
                if skipLine: # header of description file
                    skipLine = False
                    out_file.write( line )
                else:
                    fields = line.strip().split( fieldSep ) # split description line on field separator
                    sgb_id = fields[ sgb_idx ] # SGB ID
                    estimated_taxonomy = fields[ estTax_idx ] # estimated taxonomy
                    estimated_taxonomic_level = fields[ levelEstTax_idx ] # level of the estimated taxonomy
                    if estimated_taxonomic_level == 'Species':
                        # if the estimated taxonomy level is 'Species', then add a new taxonomic level 't__SGBid'
                        estimated_taxonomy += '|t__SGB' + sgb_id
                        # update the estimated taxonomic level
                        #estimated_taxonomic_level = '' # TODO
                    elif estimated_taxonomic_level == 'Genus':
                        # if the estimated taxonomy level is 'Genus', then add two new taxonomic levels 's__GENUS_SGBid|t__SGBid'
                        s_level = '|' + estimated_taxonomy.split( '|' )[-1].replace( 'g__', 's__' ) + '_SGB' + sgb_id
                        t_level = '|t__SGB' + sgb_id
                        estimated_taxonomy += s_level + t_level
                        # update the estimated taxonomic level
                        #estimated_taxonomic_level = '' # TODO
                    elif estimated_taxonomic_level == 'Family':
                        # if the estimated taxonomy level is 'Family', then add three new taxonomic levels 'g__GGBid|s__GGBid_SGBid|t__SGBid'
                        ggb_id = get_id( associations, sgb_id, 0, 1 )
                        if ggb_id is not None:
                            g_level = '|g__GGB' + ggb_id 
                            s_level = '|s__GGB' + ggb_id + '_SGB' + sgb_id
                            t_level = '|t__SGB' + sgb_id
                            estimated_taxonomy += g_level + s_level + t_level
                            # update the estimated taxonomic level
                            #estimated_taxonomic_level = '' # TODO
                    elif estimated_taxonomic_level == 'Other':
                        # if the estimated taxonomy level is 'Other', build the new taxonomic levels 
                        # with k__ and p__ retrieved from the 'Full Taxonomic Label of the Closest Genome'
                        # and 'c__CFBGid|o__OFBGid|f__FBGid|g__GGBid|s__GGBid_SGBid|t__SGBid'
                        full_tax_label_of_closest_genome = fields[ closestTax_idx ].strip().split( '|' )
                        k_p_closest = full_tax_label_of_closest_genome[0] + '|' + full_tax_label_of_closest_genome[1]
                        fbg_id = get_id( associations, sgb_id, 0, 2 )
                        ggb_id = get_id( associations, sgb_id, 0, 1 )
                        if fbg_id is not None and ggb_id is not None:
                            c_level = '|c__CFBG' + fbg_id
                            o_level = '|o__OFBG' + fbg_id
                            f_level = '|f__FBG' + fbg_id
                            g_level = '|g__GGB' + ggb_id
                            s_level = '|s__GGB' + ggb_id + '_SGB' + sgb_id
                            t_level = '|t__SGB' + sgb_id
                            estimated_taxonomy = k_p_closest + c_level + o_level + f_level + g_level + s_level + t_level
                            # update the estimated taxonomic level
                            #estimated_taxonomic_level = '' # TODO
                    # update the estimated taxonomy
                    fields[ estTax_idx ] = estimated_taxonomy
                    # update the level of the estimated taxonomy
                    fields[ levelEstTax_idx ] = estimated_taxonomic_level
                    out_file.write( build_description_line( fields ) )
    out_file.close()

if __name__ == '__main__':
    argv = sys.argv[1:]
    sgbs_description_filepath = ''
    all_005_015_030_average_clusters_filepath = ''
    new_sgbs_description_filepath = ''

    try:
        opts, args = getopt.getopt( argv, 'hs:c:o:', [ 'sgbs=', 'associations=', 'out=' ] )
    except getopt.GetoptError:
        print( 'fix_taxonomy.py -s <sgbsfile> -c <associationsfile> -o <outputfile>' )

    for opt, arg, in opts:
        if opt == '-h':
            print( 'fix_taxonomy.py -s <sgbsfile> -c <associationsfile>' )
            sys.exit()
        elif opt in ( '-s', '--sgbs' ):
            sgbs_description_filepath = arg
        elif opt in ( '-c', '--associations' ):
            all_005_015_030_average_clusters_filepath = arg
        elif opt in ( '-o', '--output' ):
            new_sgbs_description_filepath = arg
    
    if os.path.isfile( sgbs_description_filepath ) and os.path.isfile( all_005_015_030_average_clusters_filepath ):
        if not os.path.exists( new_sgbs_description_filepath ):
            # assume --clusters formatted as a tsv
            # retrieveOnly='GCA_' : retrieve the reference genomes only
            # 'GCA_' is the prefix of the reference genome labels
            #associations = load_association_file( all_005_015_030_average_clusters_filepath, retrieveOnly='GCA_' )
            associations = load_association_file( all_005_015_030_average_clusters_filepath )
            # assume --sgbs formatted as a csv
            t0 = time.time()
            create_description_file( sgbs_description_filepath, associations, new_sgbs_description_filepath )
            t1 = time.time()
            print( 'New description file created in {}s'.format( int( t1 - t0 ) ) )
        else:
            print( 'The specified output file already exists' )
    else:
        print( 'One or more input file paths are not valid' )
