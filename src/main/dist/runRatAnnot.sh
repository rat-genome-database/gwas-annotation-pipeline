#!/usr/bin/env bash
#
# gwas-annotation-pipeline
#
. /etc/profile
APPNAME=gwas-annotation-pipeline
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

#EMAILLIST=mtutaj@mcw.edu,llamers@mcw.edu
EMAILLIST="mtutaj@mcw.edu llamers@mcw.edu"

APPDIR=/home/rgddata/pipelines/$APPNAME
cd $APPDIR

java -Dspring.config=$APPDIR/../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$APPDIR/properties/log4j2.xml \
    -Xmx20g -jar lib/$APPNAME.jar -ratAnnotRun "$@" > runRat.log 2>&1

mailx -s "[$SERVER] GWAS Annotation Pipeline Run" $EMAILLIST < $APPDIR/logs/ratSummary.log
