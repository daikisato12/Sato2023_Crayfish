use warnings;
use strict;

my %code = (0 => 'A', 1 => 'T', 2 => 'C', 3 => 'G', 4 => 'N', 5 => 'Del');
my $file = $ARGV[0]; #Nhe3_LOC123755016.sync
my $outfile = $file;
substr($outfile, -5, 5,'_AF.txt');
open(OUT, ">./$outfile");
print OUT "Chr\tCoord\tRef\tAomori\tAtchafalaya\tKamakura\tNewOrleans\tOkinawa\tTonden\tTriunfo\n";
open(IN, $file);
OUTER: while(<IN>){
chomp;
	my @li = split(/\s+/, $_);
	my $loc = $li[1];
	my %hash;
	my %hash2;
	my %al;
	my $af;
	for (my $a = 3; $a<=$#li; $a++){
		my $cov = 0;
		my $cov_ref;
		my @al = split(/\:/, $li[$a]);
		for (my $b=0; $b<=$#al; $b++){
			my $nuc = $code{$b};
			if($al[$b] > 0){
				$hash2{$a} = $al[$b] if($nuc eq uc($li[2]));
				$hash{$a}{$nuc} = $al[$b];
				$al{$nuc} = 1;
			}
			$cov += $al[$b];
		}
		next OUTER if($cov < 4);
		if(defined $hash2{$a}){
			$hash2{$a} = $hash2{$a} / $cov;
		}else{
			$hash2{$a} = 0;
		}
	}
	next if(keys(%al) < 2);
	print OUT "$li[0]\t$li[1]\t$li[2]";
	for my $key (sort {$a <=> $b} keys %hash){
		my $t = 0;
		for my $key2 (sort {$a cmp $b} keys %{$hash{$key}}){
			if($t == 0){
				print OUT "\t$hash2{$key} \($key2\:$hash{$key}{$key2}";
			}else{
				print OUT ", $key2\:$hash{$key}{$key2}";
			}
			$t++;
		}
		print OUT ")";
	}
	print OUT "\n";
}
