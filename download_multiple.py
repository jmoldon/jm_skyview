import numpy as np
import os

infile = 'source_list.txt'
mypath = './images/'
sources = np.atleast_2d(np.genfromtxt(infile, delimiter=',', dtype="|U", autostrip=True))
surveys = ['NVSS', 'VLA FIRST (1.4 GHz)', 'TGSS ADR1', 'SDSSr', 'DSS']
for source in sources:
    for survey in surveys:
        command = f'python -W ignore jm_skyview.py -s "{survey}" -p {mypath} {source[0]}'
#        print(command)
        os.system(command)
    

