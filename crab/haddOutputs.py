#A script to combine the outputs of crab jobs


import os
import argparse
import subprocess

#-------------------------------------------------------------------------------------------------------------------------------------------------------#

def parseArgs():
    argparser = argparse.ArgumentParser(description="hadd together output root files from crab jobs")
    argparser.add_argument("-d", "--dir", required = True, action = "store", help="CRAB submission directory relative to ./CrabSubmits/")
    argparser.add_argument("-o", "--out", required = False, action = "store", help="Output directory relative to /store/user/bbarton/TaustarToTauTauZ/")
    argparser.add_argument("-n", "--noHadd", action = "store_true", help="If specified, will not try to hadd files and will only produce dir list")
    argparser.add_argument("-f", "--force", required = False, action="store_true", help="If specified, will add the -f flags to hadd and xrdcp commands")
    args = argparser.parse_args()

    return args

#-------------------------------------------------------------------------------------------------------------------------------------------------------#

#Gets all of the necessary information from the crab submission logs and status to construct the EOS directory names containing the output root files
def buildDirList(args):
    print("\nBuilding directory list..." )
    
    eosDirList = []
    dirBase = "/store/user/bbarton/TaustarToTauTauZ/Crab/"
    dirBase += args.dir

    dirPath = "CrabSubmits/" + args.dir + "/"
    for subDir in os.listdir(dirPath):
        command = 'crab status -d ' + dirPath + subDir + ' | grep "Task name:" ' 
        stdout, stderr  = subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        splitOut = stdout.split(":")
        subTime = splitOut[1].strip()
        jobName = splitOut[2][8:].strip()
        #print("subTime=" + subTime)
        #print("jobName=" + jobName)

        command = 'less ' + dirPath + subDir + '/crab.log | grep "config.Data.inputDataset ="'
        stdout, stderr  = subprocess.Popen(command, universal_newlines=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if len(stdout) == 0:
            print("ERROR: No dataset found")
            exit(2)
        dataset = stdout.split("'")[-2].split("/")[1]
        #print("dataset="+dataset)

        eosDirPath = dirBase + "/" + dataset + "/" + jobName + "/" + subTime + "/0000/"
        #print(eosDirPath)
        eosDirList.append(eosDirPath)

    print("... done building directory list\n")
    return eosDirList

#-------------------------------------------------------------------------------------------------------------------------------------------------------#

def haddFiles(dirList, args):
    print("\nhadd'ing files together...")

    if args.force:
        force = "-f "
    else:
        force = ""
    
    fileList = []
    for dirPath in dirList:
        targetName = dirPath.split("/")[-4][5:] + ".root"
        #print(targetName)
        fileList.append(targetName)
        
        os.system("hadd " + force + targetName + " `xrdfsls -u " + dirPath + "`")

    return fileList
    print("...done hadd'ing files together")

#-------------------------------------------------------------------------------------------------------------------------------------------------------#

def copyFiles(fileList, args):
    print("\nCopying files...")

    outdir = os.environ["XRDURL"] + "/store/user/bbarton/TaustarToTauTauZ/" + args.out
    if outdir[-1] != "/":
        outdir += "/"

    if args.force:
        force = "-f "
    else:
        force = ""
    
    for fil in fileList:
        os.system("xrdcp " + force + fil + " " + outdir)
        
    print("...done copying files")

#-------------------------------------------------------------------------------------------------------------------------------------------------------#

if __name__ == "__main__":
    args = parseArgs()
    dirList = buildDirList(args)
    if not args.noHadd:
        fileList = haddFiles(dirList, args)
        if args.out:
            copyFiles(fileList, args)
    else:
        print("dirList = " + str(dirList))

        
    
