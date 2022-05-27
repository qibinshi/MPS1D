#!/usr/bin/perl
#
# modified from Zhu Lupei's cap.pl

$usage = 
" =====  Multiple  point  sources  inversion  =====

  Data preparation:
  
   1) 3-component waveforms: *.HN[ZNE] 
   2) set origin time as 0, and set latitude,longitude in header
   3) weight file 

  Greens function library:

   So far, we use FK to compute the library

  Time window determination:
    
   1) If data has been marked with t1 t2 t3 t4, will use them to cut P and S
   2) If data is not marked but apparent velocities are given, they will be 
      used to estimate t1 t2 t3 t4
   3) Otherwise, Using the tp, ts in Greens functions

  Usage: mpsmp.pl [-D<w1/p1/p2>] [-Ggreen/dep_int] [-Hdt] [-L<tau>] [-N<n>] 
  [-P[<Yscale[/Xscale[/k]]]>] [-Qnof] [-S<s1/s2[/tie]>] [-T<m1/m2>] [-V<vp/vl/vr>] 
  [-Wi] [-Xn] [-Zstring] event_dirs

    -A  taper option (default is off)
    -B  vertical weight
    -D	weight for Pnl (w1) and distance scaling powers for Pnl (p1) and 
        surface waves (p2). (2/1/0.5)
    -G  Green's function directory and depth interal of green's function lib
    -H  dt (0.2).
    -R  search range for strike, dip and rake
    -L  length of the Markov chain
    -P	generate waveform-fit plot with plotting scale.
    	Yscale: inch per cm/s at 100 km. (100000)
	Xscale: seconds per inch. (40)
	append k if one wants to keep those waveforms.
    -Q  number of freedom per sample (1)
    -S	max. time shifts in sec for Pnl and surface waves (1/5) and
    -T	max. time window lengths for Pnl and surface waves (35/70)
    -V	apparent velocities for Pnl, Love, and Rayleigh waves (off)
    -R  the range of longitude and latitude.
    -W  use displacement for inversion; 1=> data in velocity; 2=> data in disp
    -X  specify the file of multiple point sources, with format of: lon lat
        dep_src dist_in_fault vr t_rupture slip str dip rake rise mu
    -Z  specify a different weight file name (weight.dat)
";

@ARGV > 1 || die $usage;
#================defaults==================
$green="/home/Qibin/greenlib/premkyushu";
$dir='';
$disp=0;
$mltp=0;
$io_taper=0;
$weight="weight.dat";

$plot = 0;
$amplify = 100000;
$sec_per_inch = 40;
$keep = 0;
$dura = 0;
$rise = 0.25;

($m1, $m2) = (35,70);

$max_shft1=1;
$max_shft2=5;

$weight_of_pnl=2;
$power_of_body=1;
$power_of_surf=0.5;

($vp, $love, $rayleigh) = (7.8, 3.5, 3.1);

$nof = 0.1;
$dt  = 0.1;

$itv10 = 0.05;
$itv11 = 0.05;

foreach (grep(/^-/,@ARGV)) {
   $opt = substr($_,1,1);
   @value = split(/\//,substr($_,2));
   if ($opt eq "A") {
     $io_taper = 1;
   } elsif ($opt eq "D") {
     ($weight_of_pnl,$power_of_body,$power_of_surf)=@value;
   } elsif ($opt eq "G") {
     ($green) = substr($_,2);
   } elsif ($opt eq "H") {
     ($dt,$ddd) = @value;
   } elsif ($opt eq "P") {
     $plot = 1;
     $amplify = $value[0] if $#value >= 0;
     $sec_per_inch = $value[1] if $#value > 0;
     $keep = 1 if $#value > 1;
   } elsif ($opt eq "Q") {
     $nof = $value[0];
   } elsif ($opt eq "S") {
     ($max_shft1, $max_shft2) = @value;
   } elsif ($opt eq "R") {
     ($lo1, $la1, $lo2, $la2, $min_dp, $max_dp, $slop0, $slop1) = @value;
   } elsif ($opt eq "I") {
     ($itv1, $itv2, $itv3, $itv4, $itv5, $itv6, $itv7, $itv8, $itv9) = @value;
   } elsif ($opt eq "C") {
     ($itv10, $itv11) = @value;
   } elsif ($opt eq "T") {
     ($m1, $m2) = @value;
   } elsif ($opt eq "L") {
     $nchain = $value[0];
   } elsif ($opt eq "B") {
     $weight_vertical = $value[0];
   } elsif ($opt eq "E") {
     ($energy, $expo) = @value;
   } elsif ($opt eq "N") {
     $no_proc = $value[0];
   } elsif ($opt eq "V") {
     ($vp, $love, $rayleigh) = @value;
   } elsif ($opt eq "W") {
     $disp = $value[0];
   } elsif ($opt eq "X") {
     $mltpl_src = $value[0];
   } elsif ($opt eq "Z") {
     $weight = $value[0];
   } else {
     printf STDERR $usage;
     exit(0);
   }
}
@event = grep(!/^-/,@ARGV);

foreach $eve (@event) {

  next unless -d $eve;

  open(WEI, "$eve/$weight") || die "could open $weight\n";
  @wwf=<WEI>;
  close(WEI);
  
  open(SSS,"$eve/$mltpl_src") || die "can't open $mltpl_src!\n";
  @psrcs = <SSS>;
  close(SSS);
  $nosrc=$#psrcs + 1;

  $cap = "/Users/qb/codes/MPS_MP/mpsmp";

  open(SRC, "| $cap $eve $green") || die "can not run $cap\n";
  print "\nEvent : $eve\nGreens: $green\n";
  print "Location Range: (lon) $lo1 $lo2  (lat) $la2 $la1 (dep) $min_dp $max_dp \n\n";
  print SRC "$itv1 $itv2 $itv3 $itv4 $itv5 $itv6 $itv7 $itv8 $itv9 $itv10 $itv11\n";
  print SRC "$power_of_body $power_of_surf $weight_of_pnl $weight_vertical $nof\n";
  print SRC "$lo1 $la1 $lo2 $la2 $min_dp $max_dp $slop0 $slop1\n";
  print SRC "$vp $love $rayleigh $dt\n";
  print SRC "$m1 $m2 $max_shft1 $max_shft2 $energy $expo\n";
  print SRC "$nchain $plot $io_taper $disp $no_proc\n";
  print SRC "$green\n";
  printf SRC "%d %f\n",$#psrcs + 1, $ddd;
  print SRC @psrcs;
  printf SRC "%d\n",$#wwf + 1;
  print SRC @wwf;
  close(SRC);
  print STDERR "inversion done\n";

  require "/Users/qb/codes/MPS_MP/mps_plt.pl";

  $model = "model_mps";
  plot:
  if ( $plot > 0 && ($? >> 8) == 0 ) {
     chdir($eve);
     &plot($model, $m1, $m2, $amplify, 6, $sec_per_inch,$f1_pnl,$f2_pnl,$f1_sw,$f2_sw,$dura,$lo1,$lo2,$la2,$la1,${nosrc});
     unlink(<${model}*.[1-2][0-6]>) unless $keep;
     chdir("../");
  }

}
exit(0);
