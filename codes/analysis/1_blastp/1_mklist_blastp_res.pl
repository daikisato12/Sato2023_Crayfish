use warnings;
use strict;

my $infile = '/data/analyzed_data/blastp/GCF_020424385.1_ASM2042438v2_translated_cds_longest_dros.blastp';
my $outfile = '/data/analyzed_data/blastp/GCF_020424385.1_ASM2042438v2_translated_cds_longest_dros.best_hits_gn.txt';
my $drosfile = '${data}/${ref_drosophila}/Drosophila_melanogaster.BDGP6.32.pep.all.fa';

my %drosid;
open(IN, $drosfile);
while(<IN>){
chomp;
	if($_=~/^>(\S+) pep primary_assembly(\S+) gene:(\S+) transcript/){
		print "$_\n";
		$drosid{$1} = $3;
	}
}

my %hash;
my %score;
open(IN, $infile);
while(<IN>){
chomp;
	next if($_=~/^#/);
	my @li = split(/\t/,$_);
	my $gene = $li[0];
	if(defined $hash{$gene}){
		if($li[-1] == $score{$gene}){
			push @{$hash{$gene}}, $li[1];
		}
	}else{
		push @{$hash{$gene}}, $li[1];
		$score{$gene} = $li[-1];
	}
}

my $tmp1 = "aa";
my $tmp2 = "aa";
open(OUT, ">$outfile");
foreach my $gene (sort {$a cmp $b} keys %hash){
	foreach	my $drospep (@{$hash{$gene}}){
		if(defined $drosid{$drospep}){
			my $drosgene = $drosid{$drospep};
			print OUT "$gene\t$drosgene\t$drospep\t$score{$gene}\n";
		}else{
			print OUT "$gene\tNA\t$drospep\t$score{$gene}\n";
		}
	}
}
