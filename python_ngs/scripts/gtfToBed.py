#!~/anaconda3/bin/python
import sys, getopt
import pandas as pd
import numpy as np
import os.path
import sys

def main(argv):
    gtfname = ''
    feature = ''
    outfile = None
    try:
        opts, args = getopt.getopt(argv,"hi:f:o:",["ifile=","feat=","out="])
    except getopt.GetoptError:
        print ('usage: gtfToBed.py -i <gtf-file to convert> -f <features> -o <outfile (optional)>')
        sys.exit(2)
    if len(opts) > 1:
        for opt, arg in opts:
            if opt == 'gi':
                print ('usage: gtfToBed.py -i <gtf-file to convert> -f <features> -o <outfile (optional)>')
                sys.exit()
            elif opt in ("-i","--ifile"):
                gtfname = arg
            elif opt in ("-f","--feat"):
                feature = arg
            elif opt in ("-o","--out"):
                outfile = arg
        print("Inputfile: "+gtfname)
        print("Features:  "+feature)
        if outfile==None:
            outfile = gtfname.replace('.gtf','.' + feature + '.bed')
        print("Outfile:   "+outfile)
    else:
        print ('usage: gtfToBed.py -i <gtf-file to convert> -f <features>  -o <outfile (optional)>')
        sys.exit()

    while not os.path.isfile(gtfname) :
        print("> '" + gtfname + "' does not exist..")
        gtfname = input(">> Please enter GTF file to convert to BED or exit ('q'): ")
        if gtfname == 'q':
            sys.exit()

    print("> " + gtfname + " will be converted to BED")

    # read gtf file
    print("> Reading gtf file '" + gtfname + "'..")
    gtf = pd.read_csv(gtfname,sep='\t',header=None,comment='#')

    if feature=='all':
        bed = pd.DataFrame({'a':gtf[0],'b':gtf[3],'c':gtf[4],'d':gtf[8].apply(lambda x: x.replace('\"','\'')),'e':gtf[5],'f':gtf[6]})
        print(bed.head())
        bed.to_csv(outfile,sep="\t",header=None,index=False)
        print("> Printed to: " + outfile)
    else:
        # function to get feature from each row
        def getFeature(annotation,feature):
            """
            getIndex: get index of the feature in current annotation list
            """
            idx=0
            annot.ind=-1
            for i in annotation:
                if i.strip().find(feature)==0:
                    annot.ind=idx
                idx+=1
            if annot.ind > -1:
                return annotation[annot.ind].strip().split(" ")[1].replace("\"","")
            else:
                return None
        annot = gtf[8].apply(lambda x: x.split(';'))
        a = pd.DataFrame(annot.apply(lambda x: getFeature(x,feature)))

        df1 = pd.DataFrame(gtf[[0,3,4]])
        df2 = pd.DataFrame(a)
        df3 = pd.DataFrame(gtf[[5,6]])

        bed = pd.concat([df1,df2,df3],axis=1)
        bed.to_csv(outfile,sep="\t",header=None,index=False)
        print("> Printed to: " + outfile)





if __name__ == "__main__":
    main(sys.argv[1:])
