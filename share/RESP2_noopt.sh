# The same as RESP2.sh, but without automatic geometry optimization
# Written by Tian Lu (sobereva@sina.com)
# Last update: 2022-Aug-16
# Examples:
# RESP2(0.5) for singlet neutral molecule with water solvent: ./RESP2_noopt.sh maki.pdb
# RESP2(0.5) for triplet neutral molecule with water solvent: ./RESP2_noopt.sh nozomi.xyz 0 3
# RESP2(0.5) for singlet anion with ethanol solvent: ./RESP2_noopt.sh nico.mol -1 1 ethanol

#!/bin/bash
delta=0.5
level_SP="B3LYP/def2TZVP"
Gaussian=g09

export inname=$1
filename=${inname%.*}
suffix=${inname##*.}

if [ $2 ];then
	echo "Net charge = $2"
	chg=$2
else
	echo "Net charge was not defined. Default to 0"
	chg=0
fi

if [ $3 ];then
	echo "Spin multiplicity = $3"
	multi=$3
else
	echo "Spin multiplicity was not defined. Default to 1"
	multi=1
fi

if [ $4 ];then
	echo Solvent is $4
	solvent="scrf(solvent="$4")"
else
	solvent="scrf(solvent=water)"
	echo "Solvent name was not defined. Default to water"
fi

echo delta parameter is $delta

keyword_SP_gas="# "$level_SP" pop=MK IOp(6/33=2,6/42=6)"
keyword_SP_solv="# "$level_SP" "$solvent" pop=MK IOp(6/33=2,6/42=6)"

#### Convert input file to .xyz file
Multiwfn $1 > /dev/null << EOF
100
2
2
tmp.xyz
0
q
EOF


#### Single point in gas
cat << EOF > gau.gjf
%chk=gau.chk
$keyword_SP_gas

test

$chg $multi
EOF
awk '{if (NR>2) print }' tmp.xyz >> gau.gjf
cat << EOF >> gau.gjf


EOF
rm tmp.xyz

echo
echo Running single point task in gas phase via Gaussian...
$Gaussian < gau.gjf > gau.out

if grep -Fq "Normal termination" gau.out
then
	echo Done!
else
	echo The task has failed! Exit the script...
	exit 1
fi

echo Running formchk...
formchk gau.chk > /dev/null

echo Running Multiwfn...
Multiwfn gau.fchk -ispecial 1 > /dev/null << EOF
7
18
8
1
gau.out
y
0
0
q
EOF

mv gau.chg gas.chg
echo RESP charge in gas phase has been outputted to gas.chg

#### Single point in solvent
cat << EOF > gau.gjf
%chk=gau.chk
$keyword_SP_solv geom=allcheck guess=read


EOF

echo
echo Running single point task in solvent phase via Gaussian...
$Gaussian < gau.gjf > gau.out

if grep -Fq "Normal termination" gau.out
then
	echo Done!
else
	echo The task has failed! Exit the script...
	exit 1
fi

echo Running formchk...
formchk gau.chk > /dev/null

echo Running Multiwfn...
Multiwfn gau.fchk -ispecial 1 > /dev/null << EOF
7
18
8
1
gau.out
y
0
0
q
EOF

mv gau.chg solv.chg
echo RESP charge in solvent phase has been outputted to solv.chg

#### Calculate RESP2
chgname=${1//$suffix/chg}

paste gas.chg solv.chg |awk '{printf ("%-3s %12.6f %12.6f %12.6f %15.10f\n",$1,$2,$3,$4,(1-d)*$5+d*$10)}' d=$delta > $chgname
rm gau.gjf gau.fchk gau.chk gau.out

echo
echo Finished! The inputted atomic coordinates with RESP2 charges \(the last column\) have been exported to $chgname in current folder
echo Please properly cite Multiwfn in your publication according to \"How to cite Multiwfn.pdf\" in Multiwfn package
