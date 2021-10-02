use v5.12;
use warnings;
use strict;
use Time::HiRes qw(time);
use PDL;
use PDL::NiceSlice;

my @testcase = ( # radius, index_l, dispersion?, distance?
	[ 27.05, 1.5137, 63.6, 0.52 ],
	[ -16.68, 1, 1, 0.138 ], # dispersion factor for vacuum is changed to 1
	[ -16.68, 1.6164, 36.7, 0.38 ],
	[ -78.1, 1, 1, 0 ]
);

# Process the number of iterations argument, if one is supplied.
my $niter = @ARGV?shift @ARGV:1000;	    # Iteration counter
usage() if $niter !~ m/^\d+$/ || $niter < 1;

my $nik = $niter / 1000;

my $max_size=10000; # max size of ndarray
my $size=$niter>$max_size?$max_size:$niter;
my $loops=ceil $niter/$size;
$niter=$loops*$size; # Actual number of iterations.

my $all_lines=pdl(7621.0, 6869.955, 6562.816, 5895.944, 5269.557, 4861.344, 4340.477, 3968.494);
my $cdf_lines=$all_lines->index(pdl(2,3,5));
my $d_line=$all_lines->index(3);
my $clear_aperture = pdl(4);
my $height=zeroes($size)->xlinvals(0,$clear_aperture/2);
my ($aberr_lspher, $aberr_osc, $aberr_lchrom, $max_lspher,
    $max_osc, $max_lchrom);

my ($distance_paraxial, $slope_paraxial, $distance_general, $slope_general);

print << "EOD";
Ready to begin John Walker's floating point accuracy
and performance benchmark.  $niter iterations will be made.

Measured run time in seconds should be divided by $nik
to normalise for reporting results.  For archival results,
adjust iteration count so the benchmark runs about five minutes.

EOD
#    print("Press return to begin benchmark: ");
#    <>;

my $start=time;
my $count=0;
while(++$count<=$loops){
    ($distance_paraxial, $slope_paraxial)
	=trace_line_paraxial($d_line->dummy(1,$size)->copy, $height(*1)->copy);
    ($distance_general, $slope_general)
	=trace_line_general($cdf_lines->dummy(1,$size)->copy, $height(*3)->copy);
    $aberr_lspher
	=$distance_paraxial((0))-$distance_general((1));
    $aberr_osc = 1
	- ($distance_paraxial((0))*$slope_paraxial((0)))
	/ ($distance_general((1))*sin($slope_general((1))));
    $aberr_lchrom = $distance_general((2)) - $distance_general((0));
    $max_lspher = sin($slope_general((1)));
    # D light
    $max_lspher = 0.0000926 / ($max_lspher * $max_lspher);
    $max_osc = 0.0025;
    $max_lchrom = $max_lspher;
}
#    printf("Stop the timer:\a ");
#    <>;

my $elapsed=time - $start;

# Now evaluate the accuracy of the results from the last ray trace

my @outarr;

$outarr[0] = sprintf "%15s   %21.11f  %14.11f",
   "Marginal ray", $distance_general((1),(-1)), $slope_general((1),(-1));
$outarr[1] = sprintf "%15s   %21.11f  %14.11f",
   "Paraxial ray", $distance_paraxial((0),(-1)), $slope_paraxial((0),(-1));
$outarr[2] = sprintf "Longitudinal spherical aberration:      %16.11f", $aberr_lspher((-1));
$outarr[3] = sprintf "    (Maximum permissible):              %16.11f",
    $max_lspher((-1));
$outarr[4] = sprintf
   "Offense against sine condition (coma):  %16.11f",
   $aberr_osc((-1));
$outarr[5] = sprintf
   "    (Maximum permissible):              %16.11f",
   $max_osc;
$outarr[6] = sprintf
   "Axial chromatic aberration:             %16.11f",
   $aberr_lchrom((-1));
$outarr[7] = sprintf
   "    (Maximum permissible):              %16.11f",
   $max_lchrom((-1));

# Now compare the edited results with the master values from
# reference executions of this program.

my @refarr = ( # Reference results.  These happen to
	       # be derived from a run on Microsoft
	       # Quick BASIC on the IBM PC/AT.
	'   Marginal ray          47.09479120920   0.04178472683',
	'   Paraxial ray          47.08372160249   0.04177864821',
	'Longitudinal spherical aberration:        -0.01106960671',
	'    (Maximum permissible):                 0.05306749907',
	'Offense against sine condition (coma):     0.00008954761',
	'    (Maximum permissible):                 0.00250000000',
	'Axial chromatic aberration:                0.00448229032',
	'    (Maximum permissible):                 0.05306749907'
);


my $errors = 0;
for (my $i = 0; $i < 8; $i++) {
   if ($outarr[$i] ne $refarr[$i]) {
      printf("\nError in results on line %d...\n", $i + 1);
      print("Expected:  $refarr[$i]\n");
      print("Received:  $outarr[$i]\n");
      print("(Errors)    ");
      my $k = length($refarr[$i]);
      for (my $j = 0; $j < $k; $j++) {
	 print(substr($refarr[$i], $j, 1) eq substr($outarr[$i], $j, 1) ? ' ' : '^');
	 if (substr($refarr[$i], $j, 1) ne substr($outarr[$i], $j, 1)) {
	    $errors++;
	}
      }
      print("\n");
   }
}
if ($errors > 0) {
   printf("\n%d error%s in results.  This is VERY SERIOUS.\n",
      $errors, $errors > 1 ? "s" : "");
} else {
   printf("\nNo errors in results.\n");
}

say "Iterations: $niter, Elapsed time: $elapsed";


sub trace_line_paraxial {
    my ($lines, $height)=@_;
    my $object_distance=$lines->zeroes;
    my $from_index=1; # start
    my $axis_slope_angle=$lines->zeroes;
    my  $iang=$lines->zeroes;	    	# Incidence angle
    my $rang=$lines->zeroes;	    	# Refraction angle
    my $iang_sin=$lines->zeroes;	# Incidence angle sin
    my $rang_sin=$lines->zeroes;	# Refraction angle sin
    my $old_axis_slope_angle;
    for (@testcase){
	my ($radius_of_curvature, $to_index, $dispersion_factor, $position)=@$_;
	$to_index +=
	    (($all_lines(3)-$lines) /
	     ($all_lines(2) - $all_lines(5)))
	    * (($to_index - 1) / $dispersion_factor);
	if ($radius_of_curvature != 0) {
	    my ($asa_slice0, $is_slice0, $height_slice0)
		=whereND $axis_slope_angle, $iang_sin, $height, $object_distance==0;
	    $asa_slice0.=0;
	    $is_slice0.=$height_slice0/$radius_of_curvature;
	    my ($is_slice1, $od_slice1, $asa_slice1, $rh_slice1)=whereND $iang_sin,
		 $object_distance, $axis_slope_angle, $height, $object_distance!=0;
	    $is_slice1.= (($od_slice1-$radius_of_curvature)/$radius_of_curvature)
		*$asa_slice1;
	    $rang_sin = ($from_index/$to_index)*$iang_sin;
	    $rh_slice1.=$od_slice1 * $asa_slice1;
	    $axis_slope_angle += $iang_sin-$rang_sin;
	    $object_distance = $height / $axis_slope_angle;
	} else {
	    $object_distance = $object_distance * ($to_index / $from_index);
	    $axis_slope_angle = $axis_slope_angle * ($from_index / $to_index);
	}
	$from_index=$to_index;
	$object_distance -= $position;
    }
    return($object_distance, $axis_slope_angle);
}

sub trace_line_general {
    my ($lines, $height)=@_;
    my $object_distance=$lines->zeroes;
    my $from_index=1; # start
    my $axis_slope_angle=$lines->zeroes;
    #my  $iang=$lines->zeroes;	    	# Incidence angle
    my $rang=$lines->zeroes;	    	# Refraction angle
    my $iang_sin=$lines->zeroes;	# Incidence angle sin
    my $rang_sin=$lines->zeroes;	# Refraction angle sin
    my $old_axis_slope_angle;
    for (@testcase){
	my ($radius_of_curvature, $to_index, $dispersion_factor, $position)=@$_;
	$to_index +=
	    (($all_lines(3)-$lines) /
	     ($all_lines(2) - $all_lines(5)))
	    * (($to_index - 1) / $dispersion_factor);
	if ($radius_of_curvature != 0) {
	    my ($asa_slice0, $is_slice0, $height_slice0)
		=whereND $axis_slope_angle, $iang_sin, $height, $object_distance==0;
	    $asa_slice0.=0;
	    $is_slice0.=$height_slice0/$radius_of_curvature;
	    my ($is_slice1, $od_slice1, $asa_slice1)=whereND $iang_sin,
		 $object_distance, $axis_slope_angle, $object_distance!=0;
	    $is_slice1.= (($od_slice1-$radius_of_curvature)/$radius_of_curvature)
		*$asa_slice1->sin;
	    my $iang=asin($iang_sin);
	    $rang_sin = ($from_index/$to_index)*$iang_sin;
	    $old_axis_slope_angle = $axis_slope_angle->copy;
	    $axis_slope_angle += $iang-asin($rang_sin);
	    my $sagitta=2*$radius_of_curvature*sin(($old_axis_slope_angle+$iang)/2)**2;
	    $object_distance =
		$radius_of_curvature*sin($old_axis_slope_angle+$iang)/$axis_slope_angle->tan
		+ $sagitta;
	} else {
	    $rang=-asin(($from_index/$to_index)*sin($axis_slope_angle));
	    $object_distance *= $to_index * cos(-$rang)/ ($from_index*cos($axis_slope_angle));
	    $axis_slope_angle = -$rang;
	}
	$from_index=$to_index;
	$object_distance -= $position;
    }
    return($object_distance, $axis_slope_angle);
}




sub usage {
    print <<'EOD';

This is John Walker's floating point accuracy and
performance benchmark program.  You call it with

    perl fbench.pl <itercount>

where <itercount> is the number of iterations
to be executed.  Archival timings should be made
with the iteration count set so that roughly five
minutes of execution is timed.

EOD
	  exit 1;
}
