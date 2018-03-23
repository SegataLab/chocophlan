#!/usr/bin/env python

""" Given an gca and an output folder, download the corresponding gff and fna files

To run:
$ ./download_gff_fna.py id_file.txt downloads_folder

Will download the gff and fna files for all ids in the input file to a
folder named "downloads_folder" in the current working directory.
"""

import os
import sys
import tempfile
import time
import datetime

# try to import urllib.request.urlretrieve for python3
try:
    from urllib.request import urlretrieve
    from urllib.request import urlcleanup
except ImportError:
    from urllib import urlretrieve
    from urllib import urlcleanup

NCBI_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"
GCA_COLUMN=0
FTP_COLUMN=19
GFF_EXTENSION="_genomic.gff.gz"
FNA_EXTENSION="_genomic.fna.gz"

def download_file(url,file):
    """ Download the url to the file location """

    # try to download the file
    try:
        print("Downloading: " + url)
        file, headers = urlretrieve(url,file)
    except EnvironmentError:
        sys.exit("ERROR: Unable to download "+url)

    # clean the cache, fixes issue with downloading
    # two ftp files one after the next in python2 (some versions)
    urlcleanup()

def get_ncbi_assembly_info():
    """ Download the assembly data from NCBI """

    # create a tempfile for the download
    file_handle, new_file = tempfile.mkstemp(prefix="ncbi_download")

    # try to download the file
    print("Downloading the assembly data info from NCBI")
    download_file(NCBI_URL,new_file)
    
    # read in the file, ignoring headers
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

    print("Searching for ftp: "+name)
    ftp=None
    for line in data:
        if line[GCA_COLUMN].startswith(name):
            ftp=line[FTP_COLUMN]+"/"+os.path.basename(line[FTP_COLUMN])
            break

    if not ftp or not ftp.startswith("ftp"):
        ftp=None
        print("ERROR: Unable to find ftp: " + name)

    return ftp

def file_size(file):
    """ Return the size of the file in MB """

    try:
        size = os.path.getsize(file) / (1024.0**2)
    except OSError:
        size = 0

    return size

def download_gff_fna(gca_name, folder, data=None, delete=None, total_size=0):
    """ Download the gff and fna for the given gca to the folder specified """

    # get the url for the gff and fna
    base_ftp=get_gca_ftp(gca_name, data=data)
    if not base_ftp:
        return total_size
    gff_ftp=base_ftp+GFF_EXTENSION
    fna_ftp=base_ftp+FNA_EXTENSION

    # download the files
    gff_file=os.path.join(folder,os.path.basename(gff_ftp))
    fna_file=os.path.join(folder,os.path.basename(fna_ftp))

    for url,file in zip([gff_ftp,fna_ftp],[gff_file,fna_file]):
        download_file(url,file)
        print("Downloaded file: " + file)
        file_size_mb = file_size(file)
        print(file+ " size : " + str(file_size_mb) + " MB")
        if delete:
            os.remove(file)
            print("Deleted: " + file)
        total_size+=file_size_mb

    return total_size

def current_time():
    print datetime.datetime.fromtimestamp(time.time()).strftime('%m-%d-%Y %H:%M:%S')

def main():
    # get the two options from the command line
    gca_file=sys.argv[1]
    download_folder=sys.argv[2]

    # create download folder if it does not exist
    if not os.path.isdir(download_folder):
        os.makedirs(download_folder)

    # download the assembly data once
    data = get_ncbi_assembly_info()

    # download the files
    current_time()
    total_size = 0
    for line in open(gca_file):
        id=line.rstrip()
        total_size=download_gff_fna(id,download_folder,data,True,total_size)
        current_time()
    print("Total downloads: " + str(total_size) + " MB")
    current_time()

if __name__ == "__main__":
    main()

