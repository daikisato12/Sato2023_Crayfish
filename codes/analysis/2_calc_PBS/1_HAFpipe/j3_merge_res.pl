use warnings;
use strict;

my @pop = ("Aomori", "atchafalaya2", "kamakura2", "NewOrleans", "okinawa2", "tonden2", "triunfo2");
my %chr;
my $t = 0;
open(IN, "/codes/analysis/2_calc_PBS/1_HAFpipe/id_only.txt");
while(<IN>){
	chomp;
	next if($_!~/^LG/);
	$t++;
	$chr{$t} = $_;
}
close(IN);

my %hash;
foreach my $key (sort {$a <=> $b} keys %chr){
	for (my $a=0; $a<=$#pop; $a++){
		open(IN, "${data}/${mapped_reads_pool}/mpileup/HAFpipe/$pop[$a]/$pop[$a].sort.filtered.bam.$chr{$key}.afSite");
		while(<IN>){
			chomp;
			next if($_ =~ /^pos/);
			my @li = split(/\,/, $_);
			my $pos = $li[0];
			$hash{$key}{$pos}{$a} = $li[1];
		}
		close(IN);
	}
}

open(OUT, ">${data}/${mapped_reads_pool}/mpileup/HAFpipe/merge_result.tsv");
print OUT "Chr\tPos";
for my $pop (@pop){
	print OUT "\t$pop";
}
print OUT "\n";

foreach my $key (sort {$a <=> $b} keys %hash){
	foreach my $pos (sort {$a <=> $b} keys %{$hash{$key}}){
		print OUT "$chr{$key}\t$pos";
		foreach my $a (sort {$a <=> $b} keys %{$hash{$key}{$pos}}){
			print OUT "\t$hash{$key}{$pos}{$a}"
		}
		print OUT "\n";
	}
}
close(OUT);
