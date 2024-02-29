#!/usr/bin/env bash
#
# gwas-annotation-pipeline
#
. /etc/profile
APPNAME=gwas-annotation-pipeline
SERVER=`hostname -s | tr '[a-z]' '[A-Z]'`

APPDIR=/home/rgddata/pipelines/$APPNAME
cd $APPDIR

java -Dspring.config=$APPDIR/../properties/default_db2.xml \
    -Dlog4j.configurationFile=file://$APPDIR/properties/log4j2.xml \
    -jar lib/$APPNAME.jar -checkDbDnp "$@" > snpRun.log 2>&1

mailx -s "[$SERVER] GWAS check with DBSNP table run" mtutaj@mcw.edu,llamers@mcw.edu < $APPDIR/logs/snpSummary.log
