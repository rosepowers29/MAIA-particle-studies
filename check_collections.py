from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from optparse import OptionParser
import os

parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
(options, args) = parser.parse_args()


to_process = []

if os.path.isdir(options.inFile):
    for r, d, f in os.walk(options.inFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(options.inFile)

filenum=0
nevts_proc = 0
for file in to_process:
    # create a reader and open an LCIO file
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    try:
        reader.open(file)
    except Exception:
        #let it skip the bad files without breaking
        print("skipping file", filenum)
        filenum = filenum+1
        continue
    filenum=filenum+1
    # loop over all events in the file
    for ievt, event in enumerate(reader):
        if filenum % 10 == 0:
            print(" ")
            print("File "+str(filenum))
            print("Processing event " + str(ievt))

        print(event.getCollectionNames())

    reader.close()
