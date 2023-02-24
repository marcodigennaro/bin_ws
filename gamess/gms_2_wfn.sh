#!/usr/bin/sh

echo $1
outputfile=`echo $1 | sed "s/.dat/.wfn/g" `
echo $outputfile

beg_line="----- TOP OF INPUT FILE FOR BADER'S AIMPAC PROGRAM -----"
end_line="----- END OF INPUT FILE FOR BADER'S AIMPAC PROGRAM -----"

if [ -f "$outputfile" ]; then
  echo file exists $outputfile
  if [ -s "$outputfile" ]; then
    echo file not empty
  else
    echo file empty 
    fi
else 
  grep "TOP OF INPUT FILE FOR BADER'S AIMPAC PROGRAM" $1
  if [ $? = 0 ]; then
    sed  "/$beg_line/,/$end_line/!d;//d" $1 > $outputfile
  else
    echo missing BADER information
  fi
fi
