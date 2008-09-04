#!/bin/sh
## Copyright (c) Mark Olesen
## All rights reserved.
##
##     This software is distributed WITHOUT ANY WARRANTY; without even
##     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##     PURPOSE.
#
# 20080809 vtkOpenFOAMReader installation preparation script
# 20080814 Added VTK/Hybrid/vtkPOpenFOAMReader.{cxx,h}

set -e

Script=${0##*/}
srcdir=${0%/*}

# use OpenFOAM ParaView_INST_DIR variable if possible
# or try to install into the ParaView3 directory
# or try with the current working directory

unset pvdir
for i in $ParaView_INST_DIR ParaView3 .
do
   if [ -d "$i/VTK/IO" ]
   then
      pvdir=$i
      break
   fi
done


if [ -d "$pvdir" ]
then
   echo "$Script: using paraview directory '$pvdir'"
else
   echo "$Script: no suitable paraview3 directory found:"
   echo "    ParaView_INST_DIR=$ParaView_INST_DIR"
   exit 1
fi

noRebuild=1  # assume we don't need to rebuild

for i in vtkOpenFOAMReader.cxx vtkOpenFOAMReader.h
do
   if cmp $srcdir/$i $pvdir/VTK/IO/$i >/dev/null 2>&1
   then
      echo "already up-to-date:  VTK/IO/$i "
   else
      unset noRebuild
      echo "replacing:  VTK/IO/$i "
      cp -f $srcdir/$i $pvdir/VTK/IO
   fi
done

tmpFile=${TMPDIR:-/tmp}/vtkOpenFOAM.$$
# Clean up on termination and on Ctrl-C
trap 'rm -f $tmpFile 2>/dev/null; exit 0' EXIT TERM INT

readers=Servers/ServerManager/Resources/readers.xml

# replace the FOAM reader section of
#    Servers/ServerManager/Resources/readers.xml
# with our content.
perl -x $0 $srcdir/vtkOpenFOAMReader.xml $pvdir/$readers > $tmpFile
rc=$?

if [ "$rc" -eq 0 -a -f $tmpFile ]
then
   if cmp $tmpFile $pvdir/$readers > /dev/null 2>&1
   then
      echo "already up-to-date:  $readers"
      rc=0
   else
      unset noRebuild
      echo "updating:  $readers"
      cp -f $pvdir/$readers $pvdir/$readers-bak
      cp -f $tmpFile $pvdir/$readers
      rc=$?
   fi

   echo
   echo "preparation done ..."
   if [ "x$pvdir" == "x$ParaView_INST_DIR" ]
   then
      echo "If you have already built ParaView, you can accelerate the rebuild:"
      echo "    buildParaView3.3-cvs -fast"
      echo "otherwise build ParaView:"
      echo "    buildParaView3.3-cvs"
   else
      echo "   No ParaView_INST_DIR set, build ParaView3 as usual."
   fi
   if [ -n "$noRebuild" ]
   then
      echo
      echo "NOTE: all files were already up-to-date"
      echo
   fi
else
   echo "$Script: problem altering readers.xml"
fi

exit $rc
# end-of-shell

# ----------------------------------------------------------------------
# perl script for adjusting readers.xml
# we could also use sed or awk for the same thing
#
#!/usr/bin/perl -w
use strict;
my $name = qr{(?:Open)?(?i:FOAM)\s+Reader};

@ARGV == 2 or die "needs two arguments\n";

my ($newXML) = shift;
while (<>) {
    if ( /<!--\s*Beginning\s+of\s+$name\s*-->/ ...
        /<!--\s*End\s+of\s+$name\s*-->/ )
    {
        if (/<!--\s*End\s+of\s+$name\s*-->/) {
            open XML, $newXML or die "cannot read file: $newXML\n";
            print while <XML>;
        }
    }
    else {
        print;    ## pass thru
    }
}

# ---------------------------------------------------------------- end-of-file
