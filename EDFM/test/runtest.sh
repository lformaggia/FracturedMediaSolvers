#!/bin/bash
function run
{
    numtest=$1
    echo "******    RUNNING TEST$numtest   ******"
    outdir=output$1
    outtestdir=output_test$1
    logfile="test$1".log
    datafile="data$1".pot
    status=0;
    /bin/rm -r -f ${outdir} 
    /bin/rm -f ${logfile}
    mkdir ${outdir}
    ./main_test InputFile=${datafile} > ${logfile}
    echo "******    CHECKING FILES   ******"
    for f in `ls ${outdir}`;  do
	filename=`basename $f`
	gunzip ${outtestdir}/${filename}.gz
	diff -q ${outdir}/$filename ${outtestdir}/${filename}
	if test ! $?
	then
	    status=1
	fi
	gzip -9 ${outtestdir}/${filename}
    done
    if test ! ${status}
    then
	echo "TEST FAILED"    
    else
	echo "TEST SUCCESFULL"
	/bin/rm -r -f ${outdir} 
	/bin/rm -f ${logfile}
    fi
    return $status
}
run 1
run 2

