use warnings;
use strict;
# perl code1_extract_EvolvingGenes.pl

my $infile1 = "Gamma_branch_probabilities.tab";
my $infile2 = "Gamma_change.tab";

my %spename;
my %hash;
open(IN, $infile1);
while(<IN>){
	chomp;
	my @li = split(/\s+/, $_);
	if($_=~/^#/){
		for (my $a=1; $a<$#li; $a++){
			$spename{$a} = $li[$a];
		}
	}else{
		for (my $a=1; $a<$#li; $a++){
			my $spename = $spename{$a};
			$hash{$spename}{$li[0]} = 1 if($li[$a] < 0.05);
		}
	}
}
close(IN);

my %incr;
my %decr;
open(IN, $infile2);
while(<IN>){
	chomp;
	my @li = split(/\s+/, $_);
	if($_=~/FamilyID/){
		for (my $a=1; $a<$#li; $a++){
			$spename{$a} = $li[$a];
		}
	}else{
		OUTER:for (my $a=1; $a<$#li; $a++){
			my $spename = $spename{$a};
			next OUTER unless(defined $hash{$spename}{$li[0]});
			if($li[$a] =~ /\+/){
				$incr{$spename}{$li[0]} = 1;
			}else{
				$decr{$spename}{$li[0]} = 1;
			}
		}
	}
}
close(IN);

open(OUT, ">./analysis/List_EvolvingGenes.txt");
print OUT "Node\tIncrease\tDecrease\n";
foreach my $key (keys %hash){
	my $incr = keys %{$incr{$key}};
	my $decr = keys %{$decr{$key}};
	print OUT "$key\t$incr\t$decr\n";
	open(OUT2, ">./analysis/List_EvolvingGenes_$key\_incr.txt");
	open(OUT3, ">./analysis/List_EvolvingGenes_$key\_decr.txt");
	foreach my $gene (keys %{$incr{$key}}){
		print OUT2 "$gene\n";
	}
	close(OUT2);
	foreach my $gene (keys %{$decr{$key}}){
		print OUT3 "$gene\n";
	}
	close(OUT3);
}
close(OUT);
