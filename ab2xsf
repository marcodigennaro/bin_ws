#!/usr/bin/perl
#
#    Copyright (C) 2010 ABINIT group
#
#    Written by Gian-Marco Rignanese in perl
#    This is free software, and you are welcome to redistribute it
#    under certain conditions (GNU General Public License,
#    see ~abinit/COPYING or http://www.gnu.org/copyleft/gpl.txt).
#
#    ABINIT is a project of the Universite Catholique de Louvain,
#    Corning Inc. and other collaborators, see ~abinit/doc/developers/contributors.txt.
#    Please read ~abinit/doc/users/acknowledgments.html for suggested
#    acknowledgments of the ABINIT effort.
#
#    For more information, see http://www.abinit.org .
#

$r0=0.5291772083;
$rprim1[0] = 1.0; $rprim1[1] = 0.0; $rprim1[2] = 0.0;
$rprim2[0] = 0.0; $rprim2[1] = 1.0; $rprim2[2] = 0.0;
$rprim3[0] = 0.0; $rprim3[1] = 0.0; $rprim3[2] = 1.0;

$ianimstep = 1;

for($ifil=0;$ifil<=$#ARGV;$ifil++){
   @file=`cat $ARGV[$ifil]`;

   while($line=shift(@file)){
#      print "### ", $line;
      if($line =~ / natom =/){
#         print "line with natom =\n";
#         print $line;
         @parts=split(' ',$line);
#         print "parts0 = $parts[0] \n";
         while ($parts[0] ne "natom") {shift @parts;}
         $natom=$parts[2];
	 #print "natom = ", $natom;
      }
      if($line =~ /     acell[0-9]*/){
         @acell = split(' ',$line);
         shift (@acell);
         @acell[0] *= $r0;
         @acell[1] *= $r0;
         @acell[2] *= $r0;
	 #print "acell = ", @acell;
      }
      if($line =~ /  acell=/){
         @acell = split(' ',$line);
         shift (@acell);
         @acell[0] *= $r0;
         @acell[1] *= $r0;
         @acell[2] *= $r0;
	 #print "acell = ", @acell;
      }
      if($line =~ /     rprim[0-9]*/){
         @tmp=split(' ',$line);
         $rprim1[0] = $tmp[1];
         $rprim1[1] = $tmp[2];
         $rprim1[2] = $tmp[3];
         $line=shift(@file);
         @tmp=split(' ',$line);
         $rprim2[0] = $tmp[0];
         $rprim2[1] = $tmp[1];
         $rprim2[2] = $tmp[2];
         $line=shift(@file);
         @tmp=split(' ',$line);
         $rprim3[0] = $tmp[0];
         $rprim3[1] = $tmp[1];
         $rprim3[2] = $tmp[2];
	 #print "rprim = ", @rprim1,@rprim2,@rprim3;
      }
      if($line =~ /  rprim=/){
         @tmp=split(' ',$line);
	 $rprim1[0] = $tmp[1];
	 $rprim1[1] = $tmp[2];
	 $rprim1[2] = $tmp[3];
         $line=shift(@file);
         @tmp=split(' ',$line);
         $rprim2[0] = $tmp[0];
         $rprim2[1] = $tmp[1];
         $rprim2[2] = $tmp[2];
         $line=shift(@file);
         @tmp=split(' ',$line);
         $rprim3[0] = $tmp[0];
         $rprim3[1] = $tmp[1];
         $rprim3[2] = $tmp[2];
         #print "rprim = ", @rprim1,@rprim2,@rprim3;
      }

      if($line =~ /     typat[0-9]*/){
#         for($ii=1;$ii<50;$ii++){
#            if($natom>=$ii*20){
#               $line=$line.shift(@file);
##	       print "now line = ", $line
#            }
#         }
         @typ=split(' ',$line);
         shift (@typ);
         #print "typat 0 = ", @typ, "\n";
      }
      if($line =~ /      type[0-9]*/){
         for($ii=1;$ii<50;$ii++){
            if($natom>=$ii*20){
               $line=$line.shift(@file);
            }
         }
         @typ=split(' ',$line);
         shift (@typ);
         #print "type 0 = ", $typ;
      }
      #print $line,$line =~ /    xangst /;
      if($line =~ /    xangst[0-9]* /){
         @pos = split (' ', $line);
#	 print "line ", $line, "\n";
#	 print "pos ", $pos[0], $pos[1], $pos[2], $pos[3], "\n";
#   already in angstroems
         $xcoord[0] = $pos[1];
         $ycoord[0] = $pos[2];
         $zcoord[0] = $pos[3];
         $iat=1;
         while(!($line =~ /    xcart[0-9]* /)){
#            print "in xangst clause ", $line;
            $line=shift(@file);
            @pos=split(' ',$line);
#   already in angstroems
            $xcoord[$iat] = $pos[0];
            $ycoord[$iat] = $pos[1];
            $zcoord[$iat] = $pos[2];
            $iat++;
         }
      }
      if($line =~ /Cartesian coordinates \(bohr\)/){
         $line=shift(@file);
#         print "zatnum = ", $zatnum;
         print $acell;
         print $rprim1,$rprim2,$rprim3;
         $iat=0;
         while(!($line =~ /Cartesian force/)){
            @pos=split(' ',$line);
            $xcoord[$iat] = $pos[0]*$r0;
            $ycoord[$iat] = $pos[1]*$r0;
            $zcoord[$iat] = $pos[2]*$r0;
            $iat++;
            $line=shift(@file);
         }
         print "DIM-GROUP\n";
         print "3   1\n";
         print "PRIMVEC\n";
         printf ("%20.10f %20.10f %20.10f  \n",
                $rprim1[0]*$acell[0], $rprim1[1]*$acell[0], $rprim1[2]*$acell[0]); 
         printf ("%20.10f %20.10f %20.10f  \n",
                $rprim2[0]*$acell[1], $rprim2[1]*$acell[1], $rprim2[2]*$acell[1]); 
         printf ("%20.10f %20.10f %20.10f  \n",
                $rprim3[0]*$acell[2], $rprim3[1]*$acell[2], $rprim3[2]*$acell[2]); 

         print "PRIMCOORD\n";
	 $ianimstep++;
         printf ("$natom   1 \n");
         for ($iatom=0; $iatom<$natom; $iatom++)
          {
          printf ("%d  %20.10f %20.10f %20.10f\n", $zatnum[$typ[$iatom]],
              $xcoord[$iatom], $ycoord[$iatom], $zcoord[$iatom]);
          }
#         print "ATOMS\n";
#         for ($iatom=0; $iatom<$natom; $iatom++)
#          {
#          printf ("%d  %20.10f %20.10f %20.10f\n", $zatnum[$typ[$iatom]],
#              $xcoord[$iatom], $ycoord[$iatom], $zcoord[$iatom]);
#          }
      }
      if($line =~ /    zatnum   / || $line =~ /    znucl   /){
         @zatnum=split(' ',$line);
#         print "zatnum 0 = ", $zatnum, $line;
#         print "xcoord ", $xcoord[0], $xcoord[1], $xcoord[2];
#         print "ycoord ", $ycoord[0], $ycoord[1], $ycoord[2];
#         print "zcoord ", $zcoord[0], $zcoord[1], $zcoord[2];
#
         print "DIM-GROUP\n";
         print "3   1\n";
         print "PRIMVEC\n";
         printf ("%20.10f %20.10f %20.10f  \n",
                $rprim1[0]*$acell[0], $rprim1[1]*$acell[0], $rprim1[2]*$acell[0]); 
         printf ("%20.10f %20.10f %20.10f  \n",
                $rprim2[0]*$acell[1], $rprim2[1]*$acell[1], $rprim2[2]*$acell[1]); 
         printf ("%20.10f %20.10f %20.10f  \n",
                $rprim3[0]*$acell[2], $rprim3[1]*$acell[2], $rprim3[2]*$acell[2]); 

         print "PRIMCOORD\n";
	 $ianimstep++;
         printf ("$natom   1 \n");
         for ($iatom=0; $iatom<$natom; $iatom++)
          {
          printf ("%d  %20.10f %20.10f %20.10f\n", $zatnum[$typ[$iatom]],
              $xcoord[$iatom], $ycoord[$iatom], $zcoord[$iatom]);
          }
#         print "ATOMS\n";
#         for ($iatom=0; $iatom<$natom; $iatom++)
#          {
#          printf ("%d  %20.10f %20.10f %20.10f\n", $zatnum[$typ[$iatom]],
#              $xcoord[$iatom], $ycoord[$iatom], $zcoord[$iatom]);
#          }
#
      }
   }
}

#***************************************************8
