import glob
import os
import sys
import importlib

myDir = os.getcwd()
inputDirName = "/eos/cms/store/group/phys_egamma/tnpTuples/rverma/crab_CMSHLT_3313/EGamma0/crab_CMSHLT_3313/240823_124805/0000/"
jobsDirName = myDir+"/Jobs/"
#triggerName = sys.argv[1]
cfgFileName = "makeMini_cfg.py"
shName = myDir+"/RunCondor.sh"
fileList = filter(os.path.isfile, glob.glob(inputDirName + "*.root"))

os.system("rm Jobs/ -rf")

try:
    os.mkdir(jobsDirName)
except OSError as error:
    print(error)


nJob = 0    
for f in fileList: #Loop over the rootfiles in the selected directory

    print("Preparing job number %s"%str(nJob+1)+"\n")
    jobSubDir = jobsDirName+"Job_%s/"%str(nJob)
    os.mkdir(jobSubDir)
    os.system("cp "+shName+" "+jobSubDir) #Copy the bash script that runs cmsRun on condor to each Jobs directory

    
    cfgInputFile = open(myDir+"/"+cfgFileName,"r") #Open the configuration file and save its lines as a list    
    lines = cfgInputFile.readlines()
    print("The input file to be used is:",f)
    
    with open(jobSubDir+cfgFileName, 'a+') as cfgHandle: #Create a new configuration file in each Jobs directory

        #[line == "'file:"+f for line in lines if line.strip().startswith("'file:/eos/user/p/pellicci/")]
        line_index = -1
        for line in lines:
            line_index +=1
            if line.strip().startswith("fileNames = cms.untracked.vstring"):
                line = "fileNames = cms.untracked.vstring('file:"+f+"'),\n"
                lines[line_index] = line
        cfgHandle.writelines(lines) #Copy the content of the initial config file in there, but changing the input file name

    cfgInputFile.close()

    nJob += 1


#Prepare the content of the condor submission file
condor_str  = "Executable = $(filename)\n"
condor_str += "Universe = vanilla\n"
condor_str += "Should_Transfer_Files = YES\n"
condor_str += "WhenToTransferOutput = ON_EXIT\n"
condor_str += "Transfer_Input_Files = $(filename), $Fp(filename)"+cfgFileName+"\n"
#condor_str += "x509userproxy = $ENV(X509_USER_PROXY)\n"
condor_str += "getenv      = True\n"
condor_str += '+JobFlavour = "workday"\n'
condor_str += "request_memory = 4000\n"
condor_str += "request_cpus = 8\n"
condor_str += "Output = $Fp(filename)log_$(Cluster)_$(Process).stdout\n"
condor_str += "Error  = $Fp(filename)log_$(Cluster)_$(Process).stderr\n"
condor_str += "Log  = $Fp(filename)log_$(Cluster)_$(Process).log\n"
condor_str += "Arguments = $(Cluster) $(Process) "+cfgFileName+"\n"
condor_str += "Queue filename matching ("+myDir+"/Jobs/Job_*/*.sh)"

#Write condor submission file (equivalent to jdl)
condor_name = myDir+"/condor_cluster.sub"
condor_file = open(condor_name, "w")
condor_file.write(condor_str)

#Write script that will call the condor submission file created above
submit_all = open("submit_all.jobb","w")
submit_all.write("condor_submit %s\n"%condor_name)
os.system("chmod +x submit_all.jobb")
