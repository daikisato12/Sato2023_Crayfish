#!/usr/bin/perl

use 5.010;
use warnings;
use strict;

=head2 Fst calculation from popoolation2

By default Fst is calculated from the allele-frequencies (not from the allele-counts) using the standard equation as shown in "Hartl and Clark (2007): Principles of Population Genetics"
 
 Fst = (Pi_total - Pi_within) / Pi_total
 Pi_within = (Pi_population1 + Pi_population2)/ 2
 Pi: 1 - fA ^ 2 - fT ^ 2 -fC ^ 2 - fG ^ 2
 fN:  frequency of nucleotide N
 Pi_total: for the total Pi the allele frequencies of the two
     populations are averaged and Pi is calculated as shown above

if the alternative method C<--karlsson-fst> is used, Fst is calculated  from the allele-counts according to:
Karlsson et al. (2007): Efficient mapping of mendelian traits in dogs through genome-wide association, Nature genetics
=head1 AUTHORS

Robert Kofler

Christian Schloetterer

=cut

#perl 2_calc_fst_sliding_window.pl merge_result_pi.tsv ./tmp/HAF_fst_5kbp_step1kbp.tsv 5000 1000
my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $window_size = $ARGV[2];
my $step_size = $ARGV[3];

my %comb = ("Aomori-Kamakura" => 2,
            "Aomori-Okinawa" => 4,
            "Aomori-Sapporo" => 6,
            "Kamakura-Okinawa" => 8,
            "Kamakura-Sapporo" => 10,
            "Okinawa-Sapporo" => 12);

my %last;
my %pi_within;
my %pi_total;
open(IN, $infile);
while(<IN>){
    chomp;
    my @li = split(/\t/, $_);
    next if($_ =~ /^Chr/);
    foreach my $comb (keys %comb){
        my $index1 = $comb{$comb};
        my $index2 = $index1 + 1;
        my $chr = $li[0];
        my $coord = $li[1];
        $pi_within{$chr}{$coord}{$comb} = $li[$index1];
        $pi_total{$chr}{$coord}{$comb} = $li[$index2];
        if(defined $last{$chr}){
            $last{$chr} = $coord if($coord > $last{$chr});
        }else{
            $last{$chr} = $coord
        }
    }
}
close(IN);

open(OUT, ">./$outfile");
print OUT "Chr\tCoord\tNSNPs";
foreach my $comb (sort {$a cmp $b} keys %comb){
    print OUT "\t$comb"
}
print OUT "\n";
foreach my $chr (sort {$a cmp $b} keys %last){
    for (my $start = 0; $start < $last{$chr}; $start += $step_size){
        my $end = $start + $window_size;
        my $mid = ( $start + $end ) / 2;
        my $nsnps = 0;
        my %fst_denom;
        my %fst_numer;
        foreach my $comb (sort {$a cmp $b} keys %comb){
            $fst_denom{$comb} = 0;
            $fst_numer{$comb} = 0;
        }
#        print "$pi_within{$chr}\n";
        foreach my $key (sort {$a <=> $b} keys %{$pi_within{$chr}}){
#            print "$key\t$start\n";
            next if($key < $start);
#            print "$key\t$end\n";
            last if($key > $end);
#            print "$chr\t$mid\t$key\n";
            $nsnps++;
            foreach my $comb (sort {$a cmp $b} keys %comb){
                $fst_denom{$comb} += $pi_total{$chr}{$key}{$comb};
                $fst_numer{$comb} += ($pi_total{$chr}{$key}{$comb} - $pi_within{$chr}{$key}{$comb});
            }
        }
        next if($nsnps == 0);
        print OUT "$chr\t$mid\t$nsnps";
        print "$chr\t$mid\t$nsnps";
        foreach my $comb (sort {$a cmp $b} keys %comb){
#            print "total: $fst_denom{$comb}\twithin: $fst_numer{$comb}\n";
            my $fst;
            if($fst_denom{$comb} == 0){
                $fst = "NA";
            }else{
                $fst = $fst_numer{$comb} / $fst_denom{$comb};
            }
            print OUT "\t$fst";
            print "\t$fst";
        }
        print OUT "\n";
        print "\n";
    }
}
