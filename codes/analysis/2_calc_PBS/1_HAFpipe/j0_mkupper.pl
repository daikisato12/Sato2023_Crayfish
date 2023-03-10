use warnings;
use strict;

open(IN, $ARGV[0]);
open(OUT, ">$ARGV[1]");
while(<IN>){
	chomp;
	if($_=~/^>/){
		print OUT "$_\n";
	}else{
		my $str = uc($_);
		print OUT "$str\n";
	}
}
