import os

__CHOCOPhlAn_version__ = '0.2'
__MetaPhlAn2_db_version__ = '23'
__UniRef_version__ = open('data/relnotes.txt').readline().strip().split()[-1].replace('_','') if os.path.exists('data/relnotes.txt') else ''