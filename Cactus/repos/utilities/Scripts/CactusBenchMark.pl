#!/usr/local/bin/perl -s

if (($help || $h) || ($ARGV[0] =~ /\-?\-?help\b/i)) {
print <<EOC;
--> CactusBenchMark.pl <--

This program will take a parameter file, and change the local_nx/local_ny/local_nz
or global_nx/global_ny/global_nz variables, run the simulation, and return the wallclock
time.

So, you specify a minimum (grid size) value, a maximum, and a step value.  It will 
then loop through from min to max, at whatever step, and produce the time values.

Usage:
\t\$./CactusBenchMark.pl <min> <max> <step> <executable> <param file>
  
Example:
\t\$./CactusBenchMark.pl 10 50 5 Cactus/exe/cactus_bench BenchTest.par

Output: 

Size    Time
10      0.16560100
15      0.41357700
20      1.03835800
25      2.94272200
30      4.59066600
35      5.61678100
40      7.57280400
45      17.92507300
50      17.70047400
EOC

exit 1;
}

($minimum, $maximum, $step, $executable, $parameterfile) = @ARGV;

$scope = "";

# read in the original param file
open (IN, $parameterfile) || die "\nCannot open parameter file ($parameterfile): $!";
while (<IN>) {
   if ((/driver::local\_n/) && ($scope eq "")) {
      $scope = "local";
   }
   $originalcontents .= $_;
}
close (IN);

# assign a global scope, as we did not find a local one 
if ($scope eq "") {
   $scope = "global";
}

print "\nSize\tTime";

# loop through the different test grid sizes
for (my $gridsize=$minimum; $gridsize <= $maximum; $gridsize += $step) 
{
   my $newparam = $originalcontents;

   # replace the variables local_nx..y..z or global_nx..y..z 
   foreach my $cord (x..z) {
      $newparam =~ s/(driver::${scope}\_n${cord}\s*?=\s*?)\d+/${1}${gridsize}/;
   }

   $outfile = "BenchTest_${gridsize}.par";

   open (TMP, ">$outfile") || die "\nCannot open output param file $outfile: $!";
   print TMP $newparam;
   close (TMP);

   $wallclocktime = "unknown";
   $execute = "$executable $outfile"; 

   open (PRG, "$execute |") || die "\nCannot open the cactus simulation $execute: $!";
   while (<PRG>) {
      if (/Total time for simulation.*?(\d+\.\d+)/) {
         $wallclocktime = $1;
      }
   }
   close (PRG);
 
   print "\n$gridsize\t$wallclocktime";
}

print STDERR "\n\nFinished.\n";
