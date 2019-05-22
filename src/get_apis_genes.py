#!/usr/bin/env python3

import logging
import requests


#############
# FUNCTIONS #
#############

def get_og_annot(og, annot_url):
    my_params = {'id': og}
    annot_resp = requests.get(url=annot_url, params=my_params)
    annot_data = annot_resp.json()['data']
    annot_name = annot_data['name']
    return(annot_name)


def get_og_genes(og, ortho_url, ortho_params):
    my_params = ortho_params.copy()
    my_params['id'] = og
    og_resp = requests.get(url=ortho_url, params=my_params)
    og_data = og_resp.json()['data']
    my_gene_ids = []
    for match in og_data:
        genes = match['genes']
        for gene in genes:
            my_gene_id = gene['gene_id']['id']
            my_gene_ids.append(my_gene_id)
    return(sorted(set(my_gene_ids)))


def search_for_orthogroup(fbid, search_url, search_params):
    my_params = search_params.copy()
    my_params['query'] = fbid
    search_resp = requests.get(url=search_url, params=my_params)
    search_data = search_resp.json()
    orthogroup_matches = search_data['data']
    return(orthogroup_matches)


def tidy_dict(messy_dict):
    my_dict = messy_dict.copy()
    for key in my_dict.keys():
        my_dict[key] = sorted(set(my_dict[key]))
    return(my_dict)


###########
# GLOBALS #
###########

fbgn_list = snakemake.input[0]
outfile = snakemake.output[0]
# fbgn_list = 'data/dros_eye_genes.txt'
# outfile = 'test.csv'

search_url = 'http://www.orthodb.org/search'
ortho_url = 'http://www.orthodb.org/orthologs'
annot_url = 'http://www.orthodb.org/group'
search_params = {
    'level': '33392'    # holometabolous insects
}
ortho_params = {
    'species': '7460'   # apis mellifera
}

########
# MAIN #
########

def main():
    # set up log
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=snakemake.log[0],
        level=logging.INFO)

    # read the list of droso genes
    with open(fbgn_list, 'rt') as f:
        fbids = sorted(set(x.rstrip() for x in f.readlines()))

    # search for orthogroups
    fbid_to_og = {}
    for fbid in fbids:
        logging.info(f'Searching for {fbid}')
        fbid_to_og[fbid] = search_for_orthogroup(
            fbid,
            search_url,
            search_params)

    # search for apis genes
    og_to_loc = {}
    og_to_annot = {}
    for fbid in fbid_to_og:
        my_ogs = fbid_to_og[fbid]
        for og in my_ogs:
            logging.info(f'Getting genes for {og}')
            my_locs = get_og_genes(
                og,
                ortho_url,
                ortho_params)
            if og in og_to_loc:
                for loc in my_locs:
                    og_to_loc[og].append(loc)
            else:
                og_to_loc[og] = my_locs
            logging.info(f'Getting annotation for {og}')
            my_annot = get_og_annot(og, annot_url)
            if og in og_to_annot:
                og_to_annot[og].append(my_annot)
            else:
                og_to_annot[og] = [my_annot]

    # tidy up dicts
    fbid_to_og = tidy_dict(fbid_to_og)
    og_to_annot = tidy_dict(og_to_annot)
    og_to_loc = tidy_dict(og_to_loc)

    # write a table
    logging.info(f'Generating result table')
    all_oglines = []
    ogidx = []
    for fbid in fbid_to_og:
        my_ogs = fbid_to_og[fbid]
        for og in my_ogs:
            my_oglines = list(f'{fbid},{loc},{og}' for loc in og_to_loc[og])
            for ogline in my_oglines:
                all_oglines.append(ogline)
                ogidx.append(og)

    outlines = ['fbid,loc,orthodb_og,orthodb_og_annotation']
    for i in range(0, len(all_oglines)):
        for annot in og_to_annot[ogidx[i]]:
            my_annotline = f'{all_oglines[i]},{annot}'
            outlines.append(my_annotline)

    logging.info(f'Writing results to {outfile}')
    with open(outfile, 'wt') as f:
        f.write('\n'.join(outlines))
        f.write('\n')


if __name__ == '__main__':
    main()