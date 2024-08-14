#Script to generate crab configurations and/or submit crab jobs

import argparse
import datetime
import os 

## ---------------------------------------------------------------------------------------------------------------------------- ##

def parseArgs():
    argparser = argparse.ArgumentParser(description="Tool to generate and/or submit CRAB configuration files")
    argparser.add_argument("-a", "--action", action="store", required=True, choices=["GEN", "GEN_SUB", "SUB"], help="Whether to generate, generate and submit, or only submit configs")
    argparser.add_argument("-t", "--type", action="store", choices=["SIG", "MC", "DATA"], help="Whether to run on signal MC, background MC, or data. Required if action includes generating configs")
    argparser.add_argument("-y", "--year", action="append", choices=["RUN2", "RUN3", "2016", "2016post", "2017", "2018", "2022", "2022post", "2023", "2023post"], help="A year or run period to process")
    argparser.add_argument("-d", "--dataset", action="append", choices=["SIG", "M250", "M500", "M750", "M1000", "M1250", "M1500", "M1750", "M2000", "M2500", "M3000", "M3500", "M4000", "M4500", "M5000"], help="A specific dataset to process. Must also provide one or more years. Dataset takes precedent of --inFiles")
    argparser.add_argument("-e", "--executable", action="store", default="./crab_script.py", help="The executable script to use")
    argparser.add_argument("-i", "--inFiles", action="append", help="A file to use as input. Multiple files can be specified. If action=GEN/GEN_SUB, these should be the input root files. If action=SUB, this should be directory/file to submit")
    #argparser.add_argument("-n", "--nEvents", action="store", type=int, help="The maximum number of events to process")
    #argparser.add_argument("-o", "--outDir", action="store", default="./", help="The desired output directory")
    argparser.add_argument("-k", "--keepAndDrop", action="store", default="./keep_and_drop.txt", help="A path to a text file listing which branches to keep/drop")
    argparser.add_argument("--echo", action="store_true", help="Echo the arguments used")
    argparser.add_argument("--log", action="store", choices=["TRUE", "FALSE"], default="TRUE", help="Whether to log this configuration creation/submission in crabLog.txt")

    raw_args = argparser.parse_args()  
    args = {}

    args["GEN"] = raw_args.action == "GEN" or raw_args.action == "GEN_SUB"
    args["SUB"] = raw_args.action == "SUB" or raw_args.action == "GEN_SUB"

    if args["GEN"] and not (raw_args.inFiles or raw_args.year):
        print("Must provide either input file(s) via -i or year(s) via -y to generate configs")
        exit(1)
    if raw_args.dataset and not raw_args.year:
        print("Must provide a year(s) as well as a dataset via -y")
        exit(1)

    if raw_args.type == "SIG":
        args['TYPE'] = "SignalMC"
    elif raw_args.type == "SIG":
        args['TYPE'] = "BkgdMC"
    else:
        args['TYPE'] = "Data"

    args["ERA"] = 0
    args["YEARS"] = []
    if raw_args.year:
        
        for year in raw_args.year:
            #if year == "ALL":
            #    args["YEARS"].extend(["2015", "2016", "2017", "2018", "2022", "2023"])
            if year == "RUN2":
                args["YEARS"].extend(["2016", "2016post", "2017", "2018"])
                args["ERA"] = 2
            elif year == "RUN3":
                args["YEARS"].extend(["2022", "2022post", "2023", "2023post"])
                args["ERA"] = 3
            else:
                args["YEARS"].append(year)

        if args["ERA"] == 0:
            if "2022" in raw_args.year or "2022post" in raw_args.year or  "2023" in raw_args.year or "2023post":
                args["ERA"] = 3
            else:
                args["ERA"] = 2
        args["YEARS"] = set(args["YEARS"])

    args["DATASETS"] = []
    if raw_args.type == "SIG" and raw_args.dataset:
        if args["ERA"] == 3: 
            if "SIG" in raw_args.dataset or "ALL" in raw_args.dataset:
                masses = ["M250","M500", "M750", "M1000", "M1250", "M1500", "M1750", "M2000", "M2500", "M3000", "M3500", "M4000", "M4500", "M5000"]
            else:
                masses = raw_args.dataset
            
            for year in args["YEARS"]:
                for mass in masses:
                    if year == "2022":
                        args["DATASETS"].append("/TaustarToTauZ_"+mass.lower()+"_TuneCP5_13p6TeV_pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM")
                    elif year == "2022post":
                        args["DATASETS"].append("/TaustarToTauZ_"+mass.lower()+"_TuneCP5_13p6TeV_pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM")
                    elif year == "2023":
                        args["DATASETS"].append("/TaustarToTauZ_"+mass.lower()+"_TuneCP5_13p6TeV_pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM")
                    elif year == "2023post":
                        args["DATASETS"].append("/TaustarToTauZ_"+mass.lower()+"_TuneCP5_13p6TeV_pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v2/NANOAODSIM")
        else: 
            args["DATASETS"] = []
    elif raw_args.dataset:
        if raw_args.type == "MC":
            args["DATASETS"] = []
        elif raw_args.type == "DATA":
            args["DATASETS"] = []

    
    args["FILES"] = []
    if raw_args.inFiles and len(args["DATASETS"]) == 0: #DATASETS takes precedent over files
        for fil in raw_args.inFiles:
            if os.path.isfile(fil):
                args["FILES"].append(fil)
            elif os.path.isdir(fil):
                dirContents = os.listdir(fil)
                for f in dirContents:
                    if os.path.isfile(fil + "/" + f):
                        args["FILES"].append(fil + "/" + f)
            else:
                print("WARNING: could not file or directory: " + fil)
    elif raw_args.type == "SIG" and raw_args.dataset and args["ERA"] == 2: #Run2 MC was private production so not a DAS dataset
        if "SIG" in raw_args.dataset or "ALL" in raw_args.dataset:
            masses = ["M250"," M500", "M750", "M1000", "M1250", "M1500", "M1750", "M2000", "M2500", "M3000", "M3500", "M4000", "M4500", "M5000"]
        else:
            masses = raw_args.dataset
        
        for year in args["YEARS"]:
                if year == "2016":
                    year = "2015"
                elif year == "2016post":
                    year = "2016"
                for mass in masses:
                    args["FILES"].append("/store/user/bbarton/TaustarToTauTauZ/SignalMC/TauZ/taustarToTauZ_"+mass.lower()+"_"+year+".root")
        

    if args["GEN"]:
        if os.path.isfile(raw_args.executable) and os.path.isfile(raw_args.executable[:-2] + "sh"):
            args["EXE"] = raw_args.executable
        else:
            print("Could not locate desired executable file and/or associated .sh file")
            exit(1)

    
    if args["TYPE"] == "Data":
        args["LUMI"] = {"2016" : "", "2016post": "", "2017": "", "2018":"", "2022":"", "2022post":"", "2023":"", "2023post":""}
    else:
        args["LUMI"] = {}

    #args["N"] = -1
    #if raw_args.nEvents:
    #    args["N"] = raw_args.nEvents
    
    if args["GEN"]:
        if os.path.isfile(raw_args.keepAndDrop):
            args["KEEP_DROP"] = raw_args.keepAndDrop
        else:
            print("ERROR: Could not find keep/drop file")
            exit(1)

    args["LOG"] = raw_args.log == "TRUE"

    if raw_args.echo:
        for key in args.keys():
            print(key + "\t:\t" + str(args[key]))

    return args

## ---------------------------------------------------------------------------------------------------------------------------- ##

def main(args):
    if args["GEN"]:
        filelist = makeConfigs(args)
        if args["SUB"]:
            submit(args, filelist)
    elif args["SUB"]:
        submit(args, filelist=[])
        

## ---------------------------------------------------------------------------------------------------------------------------- ##

#Create the CRAB configuration files
def makeConfigs(args):
    date = datetime.datetime.now()
    typeStr = args["TYPE"] + "/"
    if not os.path.isdir("./CrabConfigs/" + args["TYPE"] + "/"):
        os.system("mkdir ./CrabConfigs/" + args["TYPE"] + "/")
    dirName = "./CrabConfigs/" + typeStr + date.strftime("%d%b%Y_%H-%M") + "/"
    if not os.path.isdir(dirName):
        os.system("mkdir " + dirName)
    print("\n\nWill write config files to " + dirName)

    #requestName = args["TYPE"] + "_run" + str(args["ERA"]) + "_" + date.strftime("%d%b%Y_%H%M")
    outDir = "/store/user/bbarton/TaustarToTauTauZ/Crab/" + typeStr + date.strftime("%d%b%Y_%H-%M") + "/"
    print("Please execute: eosmkdir " + outDir + "\n\n")
    createdFiles = []

    
    inputs = []
    dataset = False
    if len(args["DATASETS"]) > 0:
        inputs = args["DATASETS"]
        dataset = True
    else:
        inputs = args["FILES"]

    for inputName in inputs:
        filename = ""
        if dataset:
            filename, year =  datasetToName(inputName)
            filename = "crabConf_" + filename + ".py"
        else:
            filename = "crabConf_" + inputName.split("/")[-1].split(".")[0] + ".py"
            year = ""
            yrOptions = ["2016post", "2016", "2017", "2018", "2022post", "2022",  "2023post", "2023"]
            for yrOp in yrOptions: 
                if inputName.find(yrOp) >= 0: #Note, this relies on the "post" periods appearing first in the yrOptions list
                    year = yrOp
                    break
        
        with open(dirName + filename, "w+") as confFile:
            confFile.write("from WMCore.Configuration import Configuration\n")
            confFile.write("from CRABClient.UserUtilities import config\n\n")

            confFile.write("config = config()\n\n")

            confFile.write("config.section_('General')\n")
            confFile.write("config.General.requestName = '" + filename[9:-3] + "'\n")
            confFile.write("config.General.workArea = 'CrabSubmits/" + typeStr + date.strftime("%d%b%Y_%H-%M") + "'\n")
            confFile.write("config.General.transferOutputs = True\n")
            confFile.write("config.General.transferLogs = True\n\n")

            confFile.write("config.section_('JobType')\n")
            confFile.write("config.JobType.pluginName = 'Analysis'\n")
            confFile.write("config.JobType.psetName = 'PSet.py'\n")
            confFile.write("config.JobType.scriptExe = '" + args["EXE"][:-3] +".sh'\n")
            scriptArgs = str(args["ERA"])
            confFile.write("config.JobType.scriptArgs = ['arg1=%s']\n" % scriptArgs)
            confFile.write("config.JobType.inputFiles = ['" + args["KEEP_DROP"] + "', '" + args["EXE"] + "', '../scripts/haddnano.py']\n")
            confFile.write("config.JobType.allowUndistributedCMSSW = True\n\n")

            confFile.write("config.section_('Data')\n")
            if dataset:
                confFile.write("config.Data.inputDataset = '%s'\n" % inputName )
                if args["TYPE"] == "Data":
                    confFile.write("config.Data.lumiMask='%s'\n" % args["LUMI"][year])
            confFile.write("config.Data.inputDBS = 'global'\n")
            confFile.write("config.Data.splitting = 'FileBased'\n")
            confFile.write("config.Data.unitsPerJob = 1\n")
            confFile.write("config.Data.outLFNDirBase = '"+ outDir + "/'\n")
            confFile.write("config.Data.publication = False\n\n")

            confFile.write("config.section_('Site')\n")
            confFile.write("config.Site.storageSite = 'T3_US_FNALLPC'\n")

        createdFiles.append(dirName + filename)
    return createdFiles
## ---------------------------------------------------------------------------------------------------------------------------- ##

#Convert DAS dataset name to a version friendly for naming config files, histograms, etc
def datasetToName(dataset):
    year = ""
    yrOptions = {"2016":"", "2016post":"", "2017":"2017", "2018":"2018", "22EE":"2022post", "2022":"2022", "23BPix":"2023post", "2023":"2023"} #TODO 2016 handling
    for yrOp in yrOptions.keys(): 
        if dataset.find(yrOp) >= 0:
            year = yrOptions[yrOp]
            break

    if year == "":
        print("WARNING: No year found!")
        
        
    if dataset.startswith("/TaustarToTauZ"):
        name = "taustarToTauZ_" + dataset.split("_")[1] + "_" + year
        return name, year
    else:
        print("WARNING: Unrecognized dataset!")
        return dataset

## ---------------------------------------------------------------------------------------------------------------------------- ##


#Submit CRAB jobs for the files in filelist
# filelist: A list of files/directories containing files to submit
def submit(args, filelist=[]):
    print("\nSubmitting CRAB jobs:")
    logStr = "\nSubmitting CRAB jobs:\n"
    logStr += "Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M") + "\n"
    logStr += "Configs:\n"

    if len(filelist) == 0:
        filelist = args["FILES"]
    for f in filelist:
        print("Submitting: " + f)
        os.system("crab submit -c " + f)
        if args["LOG"]:
            logStr += f + "\n"
    
    if args["LOG"]:
        with open("crabLog.txt", "a+") as logFile:
            logFile.write(logStr)


## ---------------------------------------------------------------------------------------------------------------------------- ##

if __name__ == "__main__":
    args = parseArgs()
    main(args)
