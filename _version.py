import os

__CHOCOPhlAn_version__ = '201901'
__MetaPhlAn2_db_version__ = '25'
__UniRef_version__ = open('data_201901/_relnotes.txt').readline().strip().split()[-1].replace('_','') if os.path.exists('data_201901/_relnotes.txt') else ''
