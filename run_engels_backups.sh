#!/bin/sh
#Copies over files onto the department shared backup system.
#
# To use this script you need a folder called ~/backups. In this folder _link_
# to the files/directories you want backed up with ln -s path/to/file_or/dir.
# Make sure you don't link the ~/backups directory!
sdate=`date`

rm -f ~/homebackup.log
rm -f ~/monkbackup.log
rm -f ~/homebackup.err
rm -f ~/monkbackup.err

#First, home: 
#
# The way our directory is set up we use -L to copy the linked files, but -a
# implies -l, which we don't want, so turn it off with --no-l.
touch ~/homebackup.log
rsync -avzL --no-l ~/backups/home/* /d/engels2/eigenbrot/home/ > ~/homebackup.log 2> ~/homebackup.err
rsync -az ~/homebackup.log /d/engels2/eigenbrot/

#monk:
touch ~/monkbackup.log
rsync -avzL --no-l ~/backups/monk /d/engels2/eigenbrot/ > ~/monkbackup.log 2> ~/monkbackup.err
rsync -az ~/monkbackup.log /d/engels2/eigenbrot/

touch ~/engels_backuplog.txt
date=`date`
echo "Backup started on" ${sdate} >> ~/engels_backuplog.txt
echo "    Completed on" ${date} >> ~/engels_backuplog.txt

#Copy this file and the log:
rsync -az ~/snakes/run_engels_backups.sh /d/engels2/eigenbrot/
rsync -az ~/engels_backuplog.txt /d/engels2/eigenbrot/
