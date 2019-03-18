import os

__CHOCOPhlAn_version__ = '0.2'
__MetaPhlAn2_db_version__ = '25'
__UniRef_version__ = open('data/_relnotes.txt').readline().strip().split()[-1].replace('_','') if os.path.exists('data/_relnotes.txt') else ''
