#!/bin/bash

OUTPUT=InputFile_
GREP=root
MERGEFILE=$2 #L1ITMbltMerged.root
DIRECTORY=/nfs/dust/cms/user/pooja/scratch/hadron_outer/emulation_HO/CMSSW_7_2_2_patch2/src/L1Trigger/L1IntegratedMuonTrigger/test/crab3/modifiedBird_v1/
FINALFILE=MergedInputTxtFiles.txt
#EOSPATH="srm://dcache-se-cms.desy.de:8443/pnfs/desy.de/cms/tier2/store/user/psaxena/L1Trigger/HOUpgrade/Generation/SingleMuonGun/GEN-SIM-DIGI/SingleMuMinus_Winter15_FlatPt-0to200_MCRUN2_72_V3_GEN_SIM_DIGI//150316_134931"
#EOSPATH="srm://dcache-se-cms.desy.de:8443/pnfs/desy.de/cms/tier2/store/user/psaxena/privatelyGeneratedSamples/SingleMuonGun/GEN-SIM-DIGI-L1/SingleMuMinus_Summer15_FlatPt-3to200_MCRUN2_72_V3_GEN_SIM_DIGI_L1/150609_151046"
#EOSPATH="srm://dcache-se-cms.desy.de:8443/pnfs/desy.de/cms/tier2/store/user/psaxena/analyzedData/DttfFromCombinedPrimitive/GEN-SIM-DIGI-L1/crab_PrivateProduction_DttfFromCombinedPrimitives/150716_075954/"
EOSPATH="srm://dcache-se-cms.desy.de:8443/pnfs/desy.de/cms/tier2/store/user/psaxena/analyzedData/DttfFromCombinedPrimitive/HOQualityTrue/GEN-SIM-DIGI-L1/crab_PrivateProduction_DttfFromCombinedPrimitives/150721_083613/"

FILEPATH[1]=$EOSPATH/0000/failed
FILEPATH[2]=$EOSPATH/0001/failed
#FILEPATH[3]=$EOSPATH/0002
#FILEPATH[4]=$EOSPATH/0003
#FILEPATH[5]=$EOSPATH/0004
#FILEPATH[6]=$EOSPATH/0005
#FILEPATH[7]=$EOSPATH/0006
#FILEPATH[8]=$EOSPATH/0007
#FILEPATH[9]=$EOSPATH/0008
#FILEPATH[10]=$EOSPATH/0009
#FILEPATH[11]=$EOSPATH/0010



usage () {
    echo "---------------------------------------"
    echo " Copy the Text Files from /pnfs/eos to txt-file"
    echo "---------------------------------------"
    echo "      ==Possible arguments=="
    echo "  cp      - to copy the txt file name from dCache/sub-directory"
    echo "  cpfailed- to copy the txt file name from dCache"
    echo "  cpfile  - to copy the txt file name locally"
    echo "  rmtf    - to remove the txt files already created"
    echo "  rmlog   - to remove the Log files"
    echo "  rmcfg   - to remove the config files"
    echo "  rmall   - to remove the txt and root files already created"
    echo "  rmextra - to remove extra files, config files BUT not Input_txt files"
    echo "  status  - status of the bird jobs "
    echo "  del     - delete the bird jobs"
    echo "  hadd    - hadd the root files"
    echo "  help    - shows this help"
}

pause(){
    echo "Press [ENTER] to continue.."
    read -p "$*"
}


### copy the FileName from eos/subdirectory to $3
CopyFileNameFromCastor()  {
    for FileNameIndx in "${FILEPATH[@]}"
    do
	if [[ ! -e "dest_path/$FileNameIndx" ]]; then
	    SubIndex=$(echo $FileNameIndx | awk '{split($FileNameIndx, a, "/00"); print "0"a[2]}')
	    Index=$(echo $SubIndex | awk '{split($SubIndex, a, "/failed"); print a[1]}')
	    echo 'Index modified is : '$Index
	    echo "Copying fileName \"$FileNameIndx  | grep root\" to $OUTPUT$Index"
	    #srmls $FileNameIndx --count 99999 --offset 2 | grep "$GREP" | awk -F'tier2' '{print string path $GREP}' string="" path=""  > $OUTPUT$Index
	    srmls $FileNameIndx --count 99999 --offset 2 | grep "$GREP" | awk '{split($FileNameIndx, a, "tier2"); print a[2]}' > $OUTPUT$Index
	    FINALFILE=$OUTPUT$Index
	    
	    echo "@@@@@@@@@@@@ Progressing ... Please be patient..."
	    
            ### split the files into small files for job submittion
	    awk '
            NR==1 {outfile = sprintf("%s_%01d.txt", FILENAME, 0)
            }
            {print > outfile
            }
            (NR % 5) == 0 {
            close(outfile)
            outfile = sprintf("%s_%01d.txt", FILENAME, int(NR/5))
            }'  $FINALFILE 
	fi
    done
    rm -f ${OUTPUT}*[0-9][0-9][0-9]
}


copyFileFromDir() {
    if [[ ! -e "dest_path" ]]; then
	  #  Index=$(echo $FileNameIndx | awk '{split($FileNameIndx, a, "/00"); print "0"a[2]}')
        echo "Copying fileName \"$DIRECTORY  | grep root\" to $OUTPUT$Index"
	ls $DIRECTORY | grep "$GREP"  > $FINALFILE
    fi
}


### copy the FileName from eos to $3, not being used now! but it does work in copyign the file from single directory
CopyFileNameFromCastor_()  {
    if [[ ! -e "dest_path/$EOSPATH" ]]; then
	SubIndex=$(echo $EOSPATH | awk '{split($EOSPATH, a, "/00"); print "0"a[2]}')
	Index=$(echo $SubIndex | awk '{split($SubIndex, a, "/failed"); print a[1]}')
	echo "Copying fileName \"$EOSPATH  | grep root\" to $OUTPUT$Index"
	srmls $EOSPATH --count 99999 --offset 2 | grep "$GREP" | awk '{split($EOSPATH, a, "tier2"); print a[2]}' > $OUTPUT$Index
	FINALFILE=$OUTPUT$Index
	
	echo "@@@@@@@@@@@@ Progressing ... Please be patient..."
	
        ### split the files into small files for job submittion
	awk '
            NR==1 {outfile = sprintf("%s_%01d.txt", FILENAME, 0)
            }
            {print > outfile
            }
            (NR % 100) == 0 {
            close(outfile)
            outfile = sprintf("%s_%01d.txt", FILENAME, int(NR/100))
            }'  $FINALFILE 
    fi
    rm -f ${OUTPUT}*[0-9][0-9][0-9]
}


### copy the FileName from eos/subdirectory to $3
CopyFileNameFromCastorSubDir()  {
    for FileNameIndx in "${FILEPATH[@]}"
    do
	if [[ ! -e "dest_path/$FileNameIndx" ]]; then
	    Index=$(echo $FileNameIndx | awk '{split($FileNameIndx, a, "/00"); print "0"a[2]}')
	    echo 'Index is :' $Index
	    echo "Copying fileName \"$FileNameIndx  | grep root\" to $OUTPUT$Index"
	    #srmls $FileNameIndx --count 99999 --offset 2 | grep "$GREP" | awk -F'tier2' '{print string path $GREP}' string="" path=""  > $OUTPUT$Index
	    srmls $FileNameIndx --count 99999 --offset 2 | grep "$GREP" | awk '{split($FileNameIndx, a, "tier2"); print a[2]}' > $OUTPUT$Index
	    FINALFILE=$OUTPUT$Index
	    
	    echo "@@@@@@@@@@@@ Progressing ... Please be patient..."
	    
            ### split the files into small files for job submittion
	    awk '
            NR==1 {outfile = sprintf("%s_%01d.txt", FILENAME, 0)
            }
            {print > outfile
            }
            (NR % 100) == 0 {
            close(outfile)
            outfile = sprintf("%s_%01d.txt", FILENAME, int(NR/100))
            }'  $FINALFILE 
	fi
    done
    rm -f ${OUTPUT}*[0-9][0-9][0-9]
}


copyFileFromDir() {
    if [[ ! -e "dest_path" ]]; then
	  #  Index=$(echo $FileNameIndx | awk '{split($FileNameIndx, a, "/00"); print "0"a[2]}')
        echo "Copying fileName \"$DIRECTORY  | grep root\" to $OUTPUT$Index"
	ls $DIRECTORY | grep "$GREP"  > $FINALFILE
    fi
}

# Remove any earlier versions of the split output files.
RemoveInputFile() {
    rm -f ${OUTPUT}*[0-9][0-9][0-9]*.txt
    rm -f ${OUTPUT}*[0-9][0-9][0-9]
}

RemoveAllInputFile() {
#    rm -f ${OUTPUT}*[0-9][0-9][0-9]*.txt
#    rm -f ${OUTPUT}*[0-9][0-9][0-9]
    rm -f runall.zsh*
    rm allfiles
    rm *00*root
    rm run2_*_00*py
}

RemoveLogFiles() {
    rm -f runall.zsh*
    rm allfiles
}

RemoveConfigFiles() {
    rm run2_*_00*py
}



#hadd the root files                                                             
AddRootFiles() {
    hadd $MERGEFILE *root
#    if [ -s $MERGEFILE ] && (( $(stat -c %s $MERGEFILE) > 104857 )); then
#        rm *root
#    fi
    echo "@@@@@@ Success HaddED to outputfile:" $MERGEFILE
}


if [ $# -lt 1 ]; then
    usage
    exit 1
fi

if [ $1 = "help" ]; then
    usage
    exit 0

elif [ "$1" = "rmtf" ]; then
    RemoveInputFile

elif [ "$1" = "rmall" ]; then
    RemoveAllInputFile

elif [ "$1" = "rmlog" ]; then
    RemoveLogFiles

elif [ "$1" = "rmcfg" ]; then
    RemoveConfigFiles

elif [ "$1" = "cpfile" ]; then
    copyFileFromDir

elif [ "$1" = "cp" ]; then
    CopyFileNameFromCastorSubDir 

elif [ "$1" = "cpfailed" ]; then
    CopyFileNameFromCastor
    
elif [ "$1" = "hadd" ]; then
    if [ $# -ne 2 ]; then
	echo 'Error: OutputFile.root Missing '
	exit 2
    fi
    AddRootFiles

elif [ "$1" = "status" ]; then
    qstat -u 'pooja'

elif [ "$1" = "del" ]; then
    qdel -f -u 'pooja'
fi


#srm://dcache-se-cms.desy.de:8443/pnfs/desy.de/cms/tier2//store/user/psaxena/L1Trigger/HOUpgrade/Analysis/PrivatelyGeneratedMuonGun/SingleMuonGun/crab_PrivateProduction_MuonGun_Pt-2to200_gun/150206_135953/

#407202849 /abcd/directory/sub_directory/Samsung/150206_135953/0000/L1ITMBLT_54.root
#     infile=['/abcd/directory/sub_directory/Samsung/150206_135953/0000/L1ITMBLT_54.root']
#     infile.append('/abcd/directory/sub_directory/Samsung/150206_135953/0000/L1ITMBLT_93.root')
#awk -vs="'" -vs1="infile=[" -vs2="]" -vs3="infile.append(" -vs4=")" '(NR==1){print s1 s $2 s s2} (NR>1){print s3 s $2 s s4}'  Input_file






