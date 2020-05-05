from os import listdir
import os
from os.path import isfile, join

mypath = './test/'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
print(onlyfiles)

for x in onlyfiles:
    if('cpp' not in x and 'pr' not in x):
        os.system('./run.sh ./inputs/'+x.split('.txt')[0])
