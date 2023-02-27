#!/usr/bin/sh

echo $1
outputfile=`echo $1 | sed "s/.dat/.wfn/g" `
echo $outputfile

if [ -f "$outputfile" ]; then
  echo file exists
  if [ -s "$outputfile" ]; then
    echo file not empty
  else
    echo file empty
    fi
else 
  beg_line="----- TOP OF INPUT FILE FOR BADER'S AIMPAC PROGRAM -----"
  end_line="----- END OF INPUT FILE FOR BADER'S AIMPAC PROGRAM -----"
  sed  "/$beg_line/,/$end_line/!d;//d" $1 > $outputfile
  fi
