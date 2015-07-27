#!/bin/sh
source /usr/sge/default/common/settings.sh
#CONFIGFILE=run2_L1ITMuBarrelRPCHOAlgoAnalyzer_cfg.py
CONFIGFILE=$1
ls -1 *txt > allfiles

### Checking the config files
if [ $# -ne 1 ]; then
    echo 'Error: Input Config file Missing'
    exit 2 
fi

cat > runall.zsh <<EOF
#!/bin/zsh
# use current dir and current environment 
#$ -cwd
#$ -V                                      
#$ -t 1-`cat allfiles | wc -l`
#$ -l os=sld6

CONFIGFILE=$CONFIGFILE

FILE=\$(awk "NR==\$SGE_TASK_ID" allfiles)
echo 'File is: ' $FILE

ScriptName=\$(echo \$FILE | awk '{split(\$FILE, a, ".txt"); print a[1]}')
echo 'scriptName is :' $ScriptName

Index=\$(echo \$FILE | awk '{match(\$0,/\_.*\_/);print substr(\$0,RSTART,RLENGTH+2)}')
echo 'Index is :' $Index
NewFile=\$(echo \$CONFIGFILE | awk '{split(\$CONFIGFILE, a, "_"); print a[1]"_"a[2]}')
NewCONFIG=\$NewFile\$Index"_cfg.py"

cp \$CONFIGFILE \$NewCONFIG
sed -i "s/InputFile[_0-9]*.txt/\${ScriptName}.txt/g" \$NewCONFIG
sed -i s/L1ITMuonBarrelResolutionComparePlots.root/L1ITMuonBarrelResolutionComparePlots\$Index.root/g \$NewCONFIG
cmsRun \$NewCONFIG
EOF

chmod u+x runall.zsh
qsub runall.zsh

