# plot MPS solution

$pssac = "pssac2";
sub plot {
  system("gmtset MEASURE_UNIT cm");
  system("gmtset ANNOT_OFFSET_PRIMARY 0.1");
  system("gmtset ANNOT_FONT_SIZE_PRIMARY 8p");
  system("gmtset HEADER_FONT_SIZE 10p");
  system("gmtset LABEL_FONT_SIZE 10p");
  system("gmtset LABEL_OFFSET 0.0");
  system("gmtset PLOT_DEGREE_FORMAT ddd.x");
  system("gmtset ANNOT_OFFSET_PRIMARY 0.1c");

  local($model,$t1,$t2,$am,$num_com,$sec_per_inch,$c1,$c2,$c3,$c4,$dur,$lo1,$lo2,$la1,$la2,$nosrc) = @_;
  local($nn,$tt,$plt1,$plt2,$plt3,$plt4,$i,$nam,$com1,$com2,$j,$x,$y,@aa,$rslt,@name,@aztk,@dist,@src_idx);
  local $keepBad = 0;
   
  open(TTT,"> aztk.txt");  
  
  @trace = ("1/255/255/255","3/0/0/0");
  #  @name = ("P V","P R","P T","S V","S R","S T");
  @name = ("P Z","P N","P E","S Z","S N","S E");
  @src_idx = ("A","B","C","D","E","F","G","H","I","J");
  
  ($nn,$hight) = (12,26.5); ## Plot (12-2) stations on each page; the height of each page is 26.5 cm.

## the scaling numbers 
  $dla = (${la2}-${la1}) / 2;
  $dlo = (${lo2}-${lo1}) / 2;
  $scal = 3.0 * $dlo / $dla;
  $ddla = $dla / 5;
  $ddlo = $dlo / 5;
  $sepa = 0.1*$sec_per_inch;
  $tt = 3*$t1 + 3*$t2 + 5*$sepa;
  $width = 0.1*int(10*$tt/$sec_per_inch+0.5);
  @x0 = ($t1+$sepa, $t1+$sepa, $t1+$sepa, $t2+$sepa, $t2+$sepa, $t2);
  $ttt  = $t1/3;
  $tttp = $t1/$sec_per_inch;
  $ttts = $t2/$sec_per_inch;  
  $tshf = ($x0[0] * 3 + $x0[3] * 2) / $sec_per_inch;
  $dy = $nosrc - 5.6;
## interface to GMT4
  $plt1  = "| $pssac -JX$width/$hight -R0/$tt/0/$nn -Y0.4 -Ent-2 -W0.5p -M$am/0 -K -P >> $model.ps";
  $plt2  = "| pstext -JX -R -O -K -N >> $model.ps";
  $plt3  = "| psmeca -JX1/1 -R-1/1/-1/1 -Sa5  -X-2 -Y29.5 -O -K -G100 -N>> $model.ps";
  $plt4  = "| psxy -JPa1 -R0/360/0/1 -Sx0.1 -W1/255/0/0 -Y-1 -K -O -N >> $model.ps";
  $plt51 = "psbasemap -JX$tttp/0.5 -R0/$t1/0/0.5 -X-13.2 -Y-23.6 -K -O -Ba${ttt}S>>$model.ps";
  $plt52 = "psbasemap -JX$ttts/0.5 -R0/$t2/0/0.5 -X$tshf -O -Ba${ttt}S>>$model.ps"; # last plot, no -K
  $plt61 = "pscoast -R${lo1}/${la1}/${lo2}/${la2}r -JC$lo1/$la1/$scal -Df -Ba${dlo}f${ddlo}/a${dla}f${ddla}WSen -X15.5 -Y$dy -S240/255/255 -Gwhite -W0.1p,grey -K -O >> $model.ps";
  $plt62 = "psxy Active_fault_0517.gmt -R -J -M -W1.0p,grey60 -K -O >> $model.ps";
  $plt63 = "| pstext -R -J -K -O -N -G100 >> $model.ps";

  open(FFF,"$model.out");
  @rslt = <FFF>;
  close(FFF);
  @meca = split('\s+',shift(@rslt));
  @others = grep(/^#/,@rslt); @rslt=grep(!/^#/,@rslt);
  
  unlink("$model.ps") if -e "$model.ps";
  while (@rslt) {
    open(PLT, $plt1);
    $i = 0; @aztk=();
    @aaaa = splice(@rslt,0,$nn-2);
    ## plot waveforms for all stations
    foreach (@aaaa) {
      @aa = split;
      $nam = "${model}_$aa[0].";
      $nams= "$aa[0].";
      $x=0;
      $com1=$num_com+10;
      $com2=$com1+10;
      for($j=0;$j<$num_com;$j++) {
	if ($aa[5*$j+2]>0) {
	   printf PLT "%s %f %f 3/0/0/0\n",$nam.$com1,$x,$nn-$i-2; ## plot the data.
           printf PLT "%s %f %f 3/255/0/0\n",$nam.$com2,$x,$nn-$i-2; ## plot the synthetics.
	} elsif ($keepBad) {
           printf PLT "%s %f %f 3/0/255/0\n",$nam.$com1,$x,$nn-$i-2;
           printf PLT "%s %f %f 3/255/0/0\n",$nam.$com2,$x,$nn-$i-2;
        }
        $x = $x + $x0[$j];
        $com1-=1;
        $com2-=1;
      }
      $aztk[$i] = `saclst az user1 f ${nam}11`;
#$dist[$i] = `saclst gcarc f ${nams}HHE`;
      $dist[$i] = `saclst dist f ${nams}HHE`;
      $i++;
    }
    close(PLT);
   
    open(PLT, $plt2);
    $y = $nn-2;
    $i = 0;
    ## plot info for all wave fits
    foreach (@aaaa) {
      @bb = split(/\s+/,$aztk[$i]);
      $bb[1]=substr($bb[1],0,5);
      @cc = split(/\s+/,$dist[$i]);
      $cc[1]=substr($cc[1],0,4);
      @aa = split;
      $x = 0;

      printf PLT "%f %f 10 0 1 1 $aa[0]\n",$x-2*$sec_per_inch,$y; ## station name
      printf PLT "%f %f 10 0 0 1 $bb[1]\n",$x-2*$sec_per_inch,$y-0.25; ## azimuth
      printf PLT "%f %f 10 0 0 1 $cc[1]\n",$x-2*$sec_per_inch,$y+0.25; ## distance
      for($j=0;$j<$num_com;$j++) {
        printf PLT "%f %f 10 0 0 1 $aa[5*$j+6]\n",$x,$y-0.4; ## timeshift
        printf PLT "%f %f 10 0 0 1 $aa[5*$j+4]\n",$x,$y-0.6; ## cross-correlation coefficient
        $x = $x + $x0[$j];
      }
      $y--;
      $i++;
    }

    $x=0-0.8*$sec_per_inch;
    $y=$nn+1.1;
    ## plot inversion results
    for( $a = 0; $a < $nosrc; $a = $a + 1 ){
      my @meca_conjugate=`gawk -f /opt/CAP/src/subhir/to_conjugate.awk $meca[5+12*$a] $meca[6+12*$a] $meca[7+12*$a]`;
      @meca_conj=split('\s+',@meca_conjugate[0]);
      printf PLT "%f %f 11 0 0 0 $src_idx[$a]   @meca[5+12*$a]/@meca[6+12*$a]/@meca[7+12*$a] (%d/%d/%d)     dp: @meca[9+12*$a]    Mw@meca[11+12*$a]    delay: @meca[13+12*$a]    dura: @meca[15+12*$a]\n", $x, $y-0.45*$a, int(0.5+@meca_conj[0]),int(0.5+@meca_conj[1]),int(0.5+@meca_conj[2]);
    }
    ## plot misfit error
    printf PLT "%f %f 10 0 2 1  Err : @meca[5+12*$nosrc]\n",$x,$y+0.3;
    $x = 0;
    ## plot component names
    for($i=0;$i<$num_com;$i++) {
      printf PLT "%f %f 9 0 1 1 $name[$i]\n",$x,$nn-1.60;
      $x = $x+$x0[$i];
    }
    close(PLT);

    open(PLT, $plt3);
      for( $a = 0; $a < $nosrc; $a = $a + 1 ){
        printf PLT "0 %f  0 @meca[5+12*$a, 6+12*$a, 7+12*$a] 1\n", -2-2*$a;
      }
    close(PLT);

    for( $a = 0; $a < $nosrc; $a = $a + 1 ){
      open(PLT, $plt4); 
      foreach (@aztk) {
        @aa = split;
        @bb = split(/\_/,$aa[0]);
        if ($aa[2]>90.) {$aa[1] += 180; $aa[2]=180-$aa[2];}	
        printf PLT "$aa[1] %f\n",sin($aa[2]*3.14159/360)/cos($aa[2]*3.14159/360);  
      }
      close(PLT);
    }

    system("gmtset TICK_LENGTH -0.1c");
    system("$plt61");
    # system("$plt62"); 
    open(PLT, $plt63);
      for( $a = 0; $a < $nosrc; $a = $a + 1 ){
        printf PLT "@meca[8+12*$nosrc+2*$a] @meca[7+12*$nosrc+2*$a] 10 0 1 MC $src_idx[$a]\n";
      }
    close(PLT); 
    system("gmtset TICK_LENGTH 0.2c");
    system("$plt51");
    system("$plt52");
    
    print TTT @aztk;
  }
  close(TTT);
}
1;
