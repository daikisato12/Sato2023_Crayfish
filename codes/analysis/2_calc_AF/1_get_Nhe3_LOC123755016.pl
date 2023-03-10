use warnings;
use strict;

my $chr = "LG11_CM035483.1";
#my $start = 11934700;
my $start = 11880000;
#my $end = 11965843;
my $end = 12010000;

open(OUT, ">./Nhe3_LOC123755016.sync");
my $file = "../pooled_7pops.indelfiltered.sync";
open(IN, $file);
while(<IN>){
chomp;
	my @li = split(/\s+/, $_);
#	print "$li[0]\n";
	last if($li[0] =~ /^LG12/);
	next if($li[0] ne $chr);
#	next if($_!~/LG17_CM035489.1/);
#	if($li[0]  /^LG17/){
	next if($li[1] < $start);
	last if($li[1] > $end);
	print OUT "$_\n";
#	}
}
