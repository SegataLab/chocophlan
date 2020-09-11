    #!/usr/bin/env python3


__author__ = ('Francesco Beghini (francesco.beghini@unitn.it)'
              'Nicolai Karcher (karchern@gmail.com),'
              'Francesco Asnicar (f.asnicar@unitn.it),'
              'Lauren McIver (lauren.j.mciver@gmail.com),'
              'Nicola Segata (nicola.segata@unitn.it)')

__date__ = '31 Aug 2020'

from chocophlan import *
import chocophlan.utils as utils
from chocophlan._version import __UniRef_version__

def init_parse(terminating_):
    global terminating
    terminating = terminating_

NCBI_URL="https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"
GCA_COLUMN=0
BSA_COLUMN=2
TAX_COLUMN=5
FTP_COLUMN=19
GFF_EXTENSION="_genomic.gff.gz"
FNA_EXTENSION="_genomic.fna.gz"
logger = utils.setup_logger('.','CHOCOPhlAn_download_{}'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M')))

def do_download(inputs, verbose, logger):
    if not terminating.is_set():
        try:
            ftp_base, full_link, output_path = inputs
            
            ftp = ftplib.FTP(ftp_base)  # Login to ftp server
            ftp.login()

            dir_path = os.path.dirname(output_path)
            os.makedirs(dir_path, exist_ok=True)

            if not os.path.exists(output_path):
                with open(output_path, "wb") as fileout:
                    if verbose:
                        logger.info("Downloading {}".format(full_link))
                    ftp.retrbinary("RETR " + full_link, fileout.write)
                    ftp.quit()
            else:
                ftp.quit()
                if verbose:
                    logger.info("File {} already present".format(output_path))
        except Exception as e:
            terminating.set()
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)
            logger.error(str(e))
            raise
    else:
        terminating.set()


def download(config, logger, verbose=False):
    os.makedirs(config['download_base_dir']+os.path.split(config['relpath_uniref100'])[0], exist_ok=True)
    os.makedirs(config['download_base_dir']+os.path.split(config['relpath_taxdump'])[0], exist_ok=True)
    os.makedirs(config['download_base_dir']+os.path.split(config['relpath_taxonomic_catalogue'])[0], exist_ok=True)
    os.makedirs(config['download_base_dir']+os.path.split(config['relpath_uniprot_trembl'])[0], exist_ok=True)
    os.makedirs(config['download_base_dir']+os.path.split(config['relpath_uniparc'])[0], exist_ok=True)

    ### Download UniProt version
    terminating = mp.Event()
    chunksize = config['nproc']
    with mp.Pool(initializer=init_parse, initargs=(terminating,), processes=1) as pool:
        try:
            do_download_partial = partial(do_download,verbose=config['verbose'], logger=logger)

            [_ for _ in pool.imap_unordered(do_download_partial, [(config['uniprot_ftp_base'],
                  config['relnotes'],
                  config['download_base_dir'] + "/" + config['relpath_relnotes'])],
                                  chunksize=chunksize)]
        except Exception as e:
            logger.error(str(e))
            logger.error('Download failed', exit=True)
            raise

    ### uniprot XML ###
    argument_list = [(config['uniprot_ftp_base'],
                    config['uniprot_uniref100'],
                    config['download_base_dir'] + config['relpath_uniref100'])]
    
    argument_list.append((config['uniprot_ftp_base'],
                          os.path.split(config['uniprot_uniref100'])[0] + "/RELEASE.metalink",
                          config['download_base_dir'] + os.path.split(config['relpath_uniref100'])[0] + "/RELEASE_100.metalink"))

    argument_list.append((config['uniprot_ftp_base'],
                          config['uniprot_uniref90'],
                          config['download_base_dir'] + config['relpath_uniref90']))

    argument_list.append((config['uniprot_ftp_base'],
                          os.path.split(config['uniprot_uniref90'])[0] + "/RELEASE.metalink",
                          config['download_base_dir'] + os.path.split(config['relpath_uniref90'])[0] + "/RELEASE_90.metalink"))
    
    argument_list.append((config['uniprot_ftp_base'],
                          config['uniprot_uniref50'],
                          config['download_base_dir'] + config['relpath_uniref50']))
    
    argument_list.append((config['uniprot_ftp_base'],
                          os.path.split(config['uniprot_uniref50'])[0] + "/RELEASE.metalink",
                          config['download_base_dir'] + os.path.split(config['relpath_uniref50'])[0] + "/RELEASE_50.metalink"))

    argument_list.append((config['uniprot_ftp_base'],
                          config['uniprot_sprot'],
                          config['download_base_dir'] + config['relpath_uniprot_sprot']))

    argument_list.append((config['uniprot_ftp_base'],
                          config['uniprot_trembl'],
                          config['download_base_dir'] + config['relpath_uniprot_trembl']))
    
    argument_list.append((config['uniprot_ftp_base'],
                          os.path.split(config['uniprot_sprot'])[0] + "/RELEASE.metalink",
                          config['download_base_dir'] + os.path.split(config['relpath_uniprot_sprot'])[0] + "/RELEASE.metalink"))

    argument_list.append((config['uniprot_ftp_base'],
                          config['uniparc'],
                          config['download_base_dir'] + config['relpath_uniparc']))
    
    argument_list.append((config['uniprot_ftp_base'],
                          os.path.split(config['uniparc'])[0] + "/RELEASE.metalink",
                          config['download_base_dir'] + os.path.split(config['relpath_uniparc'])[0] + "/RELEASE.metalink"))

    argument_list.append((config['ebi_ftp_base'],
                          config['uniprot_protomes'],
                          config['download_base_dir'] + config['relpath_proteomes_xml']))
    
    ### idmapping file ###
    argument_list.append((config['uniprot_ftp_base'],
                          config['uniprot_idmapping'],
                          config['download_base_dir'] + config['relpath_idmapping']))
    ### Refseq taxdump ###
    argument_list.append((config['refseq_ftp_base'],
                          config['refseq_taxdump'],
                          config['download_base_dir'] + config['relpath_taxdump']))

    terminating = mp.Event()
    chunksize = config['nproc']
    with mp.Pool(initializer=init_parse, initargs=(terminating,), processes= len(argument_list)//2 ) as pool:
        try:
            if verbose:
                logger.info("Starting parallel download")
            
            # Need to define a partial function because function passed via imap require only one argument
            do_download_partial = partial(do_download,verbose=config['verbose'], logger=logger)

            [_ for _ in pool.imap_unordered(do_download_partial, argument_list,
                                  chunksize=1)]
        except Exception as e:
            logger.error(str(e))
            logger.error('Download failed', exit=True)
            raise

    if verbose:
        logger.info('Download succesful')

    files_to_check = [(config['download_base_dir'] + config['relpath_uniref100'], config['download_base_dir'] + os.path.split(config['relpath_uniref100'])[0] + "/RELEASE_100.metalink"),
    (config['download_base_dir'] + config['relpath_uniref90'], config['download_base_dir'] + os.path.split(config['relpath_uniref90'])[0] + "/RELEASE_90.metalink"),
    (config['download_base_dir'] + config['relpath_uniref50'], config['download_base_dir'] + os.path.split(config['relpath_uniref50'])[0] + "/RELEASE_50.metalink"),
    (config['download_base_dir'] + config['relpath_uniprot_sprot'],  config['download_base_dir'] + os.path.split(config['relpath_uniprot_sprot'])[0] + "/RELEASE.metalink"),
    (config['download_base_dir'] + config['relpath_uniparc'],  config['download_base_dir'] + os.path.split(config['relpath_uniparc'])[0] + "/RELEASE.metalink"),
    (config['download_base_dir'] + config['relpath_uniprot_trembl'],  config['download_base_dir'] + os.path.split(config['relpath_uniprot_trembl'])[0] + "/RELEASE.metalink")]

    find_file = etree.ETXPath("//{http://www.metalinker.org/}file")
    find_veri = etree.ETXPath("//{http://www.metalinker.org/}verification")

    with mp.Pool(processes= len(files_to_check) ) as pool:
        check_hash_partial = partial(check_hash, verbose=verbose, logger=logger)

        ret = {x[0]:x[1] for x in pool.imap_unordered(check_hash_partial, files_to_check,
                                chunksize=1)}

    for d, m in files_to_check:
        calc_hash = ret[d]
        down_hash = ""
        meta = etree.parse(m)
        entry = [x for x in find_file(meta) if x.get('name') == os.path.split(d)[1]][0]
        down_hash = entry.getchildren()[1].getchildren()[0].text


        if (calc_hash is None) or (down_hash is None):
            sys.exit("MD5 checksums not found, something went wrong!")
        # compare checksums
        if calc_hash != down_hash:
            sys.exit("MD5 checksums do not correspond! Delete previously downloaded files so they are re-downloaded")
        else:
            logger.info("{} checksum correspond".format(d))

    with open('{}/DONE'.format( config['download_base_dir']), 'w') as d_status:
        d_status.write('{}\n'.format(__UniRef_version__))

def check_hash(d, verbose, logger):
    if verbose:
        logger.info("Checking {} md5 hash...".format(d[0]))
    md5hash = hashlib.md5()
    with open(d[0],'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5hash.update(chunk)
    return (d[0], md5hash.hexdigest()[:32])

def download_file(url,file):
    try:
        response = requests.get(url, stream=True)
        with open(file, "wb") as fileout:
            for chunk in response.iter_content(chunk_size=512):
                if chunk:
                    fileout.write(chunk)
        status = 0
    except:
        status = 1
    return status

def  get_ncbi_assembly_info():
    # create a tempfile for the download
    file_handle, new_file = tempfile.mkstemp(prefix="ncbi_download")

    logger.info("Downloading the assembly data info from NCBI")
    download_file(NCBI_URL,new_file)
    
    # read in the file, ignoring headers
    logger.info("Loading assembly data...")
    data = []
    for line in open(new_file):
        if not line.startswith("#"):
            data.append(line.rstrip().split("\t"))

    # remove the temp file
    os.remove(new_file)
    return data

def get_gca_ftp(name,data=None):
    """ Get the full ftp path to the basename for all downloads for a given gca """

    if not data:
        data = get_ncbi_assembly_info()

    # allow for different versions of the gca
    name=name.split(".")[0]+"."

    #u tils.info("Searching for ftp: {}\n".format(name))
    ftp=None
    for line in data:
        if line[GCA_COLUMN].startswith(name):
            ftp=line[FTP_COLUMN]+"/"+os.path.basename(line[FTP_COLUMN])
            break
    return ftp


def download_gff_fna(gca_name, folder, data=None, delete=None, total_size=0):
    """ Download the gff and fna for the given gca to the folder specified """

    # get the url for the gff and fna
    status = 0
    base_ftp=get_gca_ftp(gca_name, data=data)
    if not base_ftp:
        return status
    gff_ftp = base_ftp.replace('ftp:','https:') + GFF_EXTENSION
    fna_ftp = base_ftp.replace('ftp:','https:') + FNA_EXTENSION

    # download the files
    gff_file=os.path.join(folder,os.path.basename(gff_ftp))
    fna_file=os.path.join(folder,os.path.basename(fna_ftp))

    for url,file in zip([gff_ftp,fna_ftp],[gff_file,fna_file]):
        status += download_file(url,file)
        if delete:
            os.remove(file)
            logger.info("Deleted: " + file)
    return status

def process(item, data, config):
    if not terminating.is_set():
        k,v = item
        ncbi_ids = dict(v['ncbi_ids'])
        id_ = ''
        ret = ()
        if 'GCSetAcc' in ncbi_ids:
            id_ = ncbi_ids['GCSetAcc']
        elif 'Biosample' in ncbi_ids:
            id_ = [line[GCA_COLUMN] for line in data if line[BSA_COLUMN] == ncbi_ids['Biosample'] and int(line[TAX_COLUMN]) == v['tax_id']]
            if len(id_): id_=id_[0]
            ret = (k, id_)
        if len(id_):
            ncbi_ids['GCSetAcc'] = id_
            try:
                download_folder = '{}/{}/{}/{}'.format(config['download_base_dir'], config['relpath_genomes'], id_.split('_')[1][0:6],id_.split('_')[1][6:9])
            except:
                logger.info('{}'.format(id_))
            os.makedirs(download_folder, exist_ok=True)
            try:
                status = download_gff_fna(id_,download_folder,data)
            except Exception as e:
                logger.error('Failed to download {}'.format(id_))
                terminating.set()

            if status:
                logger.info("{} was not found!".format(id_))
                return id_
            elif ret:
                return ret
    else:
        terminating.set()

def process_from_file(item, data, basefolder):
    if not terminating.is_set():
        try:
            download_folder = '{}/{}/{}'.format(basefolder, item.split('_')[1][0:6],item.split('_')[1][6:9])
        except:
            logger.info('{}'.format(item))
        os.makedirs(download_folder, exist_ok=True)
        try:
            status = download_gff_fna(item,download_folder,data)
        except Exception as e:
            logger.error('Failed to download {}'.format(item))
            terminating.set()

        if status:
            logger.info("{} was not found!".format(item))
            return item
    else:
        terminating.set()

def download_ncbi_from_txt(input_file, basefolder):
    # download the assembly data once
    data = get_ncbi_assembly_info()
    terminating = mpdummy.Event()
    print("Loading input file with GCA ids to download...")
    assembly_ids = [x.strip() for x in open(input_file).readlines()]
    partial_process = partial(process_from_file, data=data, basefolder=basefolder)
    with mpdummy.Pool(initializer=init_parse, initargs=(terminating, ), processes=20) as pool:
        failed = [f for f in pool.imap_unordered(partial_process, assembly_ids, chunksize=10)]

    failed = [x.split('.')[0] for x in filter(None,failed)]

    with open(os.path.join(basefolder,'failed_GCA.txt'),'w') as f:
        f.write('\n'.join(failed))

def download_ncbi_from_proteome_pickle(config, proteomes=None, taxontree=None):
    # download the assembly data once
    data = get_ncbi_assembly_info()
    terminating = mpdummy.Event()
    if not proteomes:
        logger.info("Loading proteomes data...")
        proteomes = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_proteomes']), 'rb'))
    if not taxontree:
        logger.info("Loading NCBI taxontree...")
        taxontree = pickle.load(open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_taxontree']), 'rb'))

    partial_process = partial(process, data=data, config=config)
    with mpdummy.Pool(initializer=init_parse, initargs=(terminating,), processes=config['nproc']) as pool:
        failed = [f for f in pool.imap_unordered(partial_process, ((k,v) for k, v in proteomes.items() if 'ncbi_ids' in v), chunksize=config['nproc'])]

    proteomes_to_update = filter(lambda x:type(x) == tuple, failed)
    failed = filter(lambda x:type(x) == str, failed)

    to_update = False
    for proteome, gca in proteomes_to_update:
        to_update = True
        temp_acc = list(proteomes[proteome]['ncbi_ids'])
        temp_acc.append(('GCSetAcc',gca))
        proteomes[proteome]['ncbi_ids'] = tuple(temp_acc)

    if to_update:
        pickle.dump(proteomes, open('{}/{}'.format(config['download_base_dir'], config['relpath_pickle_proteomes']), 'wb'))

    with open('failed_GCA.txt','w') as f:
        f.write('\n'.join(list(failed)))
    config['relpath_gca2taxa'] = config['relpath_gca2taxa'].replace('DATE', __UniRef_version__)
    with open('{}{}'.format(config['export_dir'],config['relpath_gca2taxa']), 'w') as f:
        f.write('GCA_accession\tUProteome\tNCBI_taxid\ttaxstr\ttaxidstr\n')
        for upid, gca, taxid in ((p, dict(proteomes[p]['ncbi_ids']).get('GCSetAcc','None'), proteomes[p]['tax_id']) for p in proteomes if 'ncbi_ids' in proteomes[p]):
            taxstr = taxontree.print_full_taxonomy(taxid)
            f.write('{}\t{}\t{}\t{}\n'.format(gca.split('.')[0], upid, taxid, '\t'.join(taxstr)))

if __name__ == '__main__':
    t0 = time.time()

    args = utils.read_params()
    utils.check_params(args, verbose=args.verbose)

    config = utils.read_configs(args.config_file, verbose=args.verbose)
    config = utils.check_configs(config)
    logger = utils.setup_logger('.','CHOCOPhlAn_download_{}'.format(datetime.datetime.today().strftime('%Y%m%d_%H%M')))
    
    download(config['download'], logger, verbose=config['download']['verbose'])

    t1 = time.time()

    logger.info('Total elapsed time {}s'.format(int(t1 - t0)))
    sys.exit(0)