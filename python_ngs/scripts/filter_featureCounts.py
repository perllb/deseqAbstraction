filename = input("Enter filename: ")
filt = input("How many reads on average is the cutoff?: ")
pattern1 = input("Enter header string to remove: ")
pattern2 = input("Enter another header string to remove: ")

rep = ""


## function to get average reads of samples 
def getAverage(lineList):
    ncol = len(lineList)
    count = 0 # variable to count the total number of reads 
    for samp in lineList[6:ncol]:
        count = count + float(samp)
    return count/(ncol-6)


## function to filter file
def filter_fCounts(file,filter=5):
    newfile = open(file+".filtered","w")
    with open(file,"r") as fin:
        i = 1
        lines=fin.readlines()
        lines[1] = lines[1].replace(pattern1,rep)
        lines[1] = lines[1].replace(pattern2,rep)
        print(lines[1])
        for line in lines:
            # if
            if i > 1:
                if i == 2:
                    newfile.write(line)
                else:
                    # if more than 5 reads on average, write
                    linelist = line.split()
                    if getAverage(linelist)>float(filter):
                        newfile.write(line)
            i = i + 1

filter_fCounts(filename,filt)

print("New, filtered count file created, named:",filename+".filtered")
