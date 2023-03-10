use warnings;
use strict;

my $chr = "LG17_CM035489.1";
#my $start = 3143416;
my $start = 3123500;
#my $end = 3154107;
my $end = 3160000;

open(OUT, ">./Spn88Ea_LOC123757200.sync");
my $file = "../pooled_7pops.indelfiltered.sync";
open(IN, $file);
while(<IN>){
chomp;
	my @li = split(/\s+/, $_);
#	print "$li[0]\n";
	last if($li[0] =~ /^LG18/);
	next if($li[0] ne $chr);
#	next if($_!~/LG17_CM035489.1/);
#	if($li[0]  /^LG17/){
	next if($li[1] < $start);
	last if($li[1] > $end);
	print OUT "$_\n";
#	}
}
