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

my @arr0 = (0, 0, 0, 0, 0, 0, 0);
my @arr1 = (1, 1, 1, 1, 1, 1, 1);
my @index = (2, 4, 6, 7);
my %pop;
open(OUT, ">$ARGV[1]"); #${data}/mapped_reads/mpileup/hafpipe/merge_result_pi.tsv
open(IN, "$ARGV[0]"); #${data}/mapped_reads/mpileup/hafpipe/merge_result.tsv
while(<IN>){
    chomp;
    my @li = split(/\t/, $_);
    if($_ =~ /^Chr/){
        for my $t (@index){
            print "$t\t$li[$t]\n";
            $pop{$t} = $li[$t];
        }
        print OUT "$li[0]\t$li[1]";
        for (my $a=0; $a<$#index; $a++){
            my $ind1 = $index[$a];
            for (my $b=$a+1; $b<=$#index; $b++){
                my $ind2 = $index[$b];
                print OUT "\tpi_within.$pop{$ind1}-$pop{$ind2}\tpi_total.$pop{$ind1}-$pop{$ind2}";
            }
        }
        print OUT "\n";
    }else{
        next if(@li[2..8] ~~ @arr0 or @li[2..8] ~~ @arr1);
        print OUT "$li[0]\t$li[1]";
        for (my $a=0; $a<$#index; $a++){
            my $ind1 = $index[$a];
            my $pop1 = $pop{$ind1};
            my $pi_within1 = 1 - $li[$ind1] ** 2 - (1 - $li[$ind1]) ** 2;
            for (my $b=$a+1; $b<=$#index; $b++){
                my $ind2 = $index[$b];
                my $pop2 = $pop{$ind2};
                my $pi_within2 = 1 - $li[$ind2] ** 2 - (1 - $li[$ind2]) ** 2;
                my $pi_within = ( $pi_within1 + $pi_within2) / 2;
                my $af_ave = ( $li[$ind1] + $li[$ind2] ) / 2;
                my $pi_total = 1 - $af_ave ** 2 - (1 - $af_ave) ** 2;
                print OUT "\t$pi_within\t$pi_total";
            }
        }
        print OUT "\n";
    }
}
