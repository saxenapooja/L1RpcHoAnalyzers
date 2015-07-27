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
#$ -l h_vmem=5G

CONFIGFILE=$CONFIGFILE

FILE=\$(awk "NR==\$SGE_TASK_ID" allfiles)
ScriptName=\$(echo \$FILE | awk '{split(\$FILE, a, ".txt"); print a[1]}')

SubIndex=\$(echo \$FILE | awk '{split(\$File, a, ".txt"); print a[1]}')
Index=\$(echo \$SubIndex | awk '{split(\$SubIndex, a, "File"); print a[2]}')

NewFile=\$(echo \$CONFIGFILE | awk '{split(\$CONFIGFILE, a, "_"); print a[1]"_"a[2]}')
NewCONFIG=\$NewFile\$Index"_cfg.py"

cp \$CONFIGFILE \$NewCONFIG
sed -i "s/InputFile[_0-9]*.txt/\${ScriptName}.txt/g" \$NewCONFIG
sed -i s/L1ITMuonBarrelResolutionComparePlots.root/L1ITMuonBarrelResolutionComparePlots\$Index.root/g \$NewCONFIG
cmsRun \$NewCONFIG
EOF

chmod u+x runall.zsh
qsub runall.zsh

