#!/usr/bin/perl

use 5.010;
use warnings;
use strict;

#perl 3_calc_PBS.pl HAF_fst_5kbp_step1kbp.tsv HAF_PBS_5kbp_step1kbp.tsv
my $infile = $ARGV[0]; #HAF_fst_5kbp_step1kbp.tsv
my $outfile = $ARGV[1]; #HAF_PBS_5kbp_step1kbp.tsv
my %pop = ("Aomori-Kamakura-Sapporo" => [3, 5, 7],
            "Aomori-Okinawa-Sapporo" => [4, 5, 8],
            "Kamakura-Okinawa-Sapporo" => [6, 7, 8]);

my %comb;
my %maxfst;
open(IN, $infile);
while(<IN>){
    chomp;
    my @li = split(/\t/, $_);
    next if($_ =~ /NA/);
    if($_ =~ /^Chr/){
        for (my $a=3; $a<=$#li; $a++){
            $comb{$a} = $li[$a];
            my @pop = split(/\-/, $li[$a]);
            $maxfst{$pop[0]}{$pop[1]} = 0;
        }
    }else{
        for (my $a = 3; $a <= $#li; $a++){
            my $comb = $comb{$a};
            my @pop = split(/\-/, $comb);
            $maxfst{$pop[0]}{$pop[1]} = $li[$a] if($li[$a] < 1 and $li[$a] > $maxfst{$pop[0]}{$pop[1]});
        }
    }
}

open(OUT, ">$outfile");
print OUT "Chr\tCoord\tNSNPs";
foreach my $key (keys %pop){
    my @pop = split(/\-/, $key);
    for my $pop (@pop){
        print OUT "\t$pop.$key"
    }
}
print OUT "\n";

open(IN, $infile);
while(<IN>){
    chomp;
    my @li = split(/\t/, $_);
    next if($_ =~ /NA/);
    next if($_ =~ /^Chr/);
    print OUT "$li[0]\t$li[1]\t$li[2]";
    foreach my $key (keys %pop){
        my @pop = split(/\-/, $key);
        my $pop1 = $pop[0];
        my $pop2 = $pop[1];
        my $pop3 = $pop[2];
        my $id1 = $pop{$key}[0];
        my $id2 = $pop{$key}[1];
        my $id3 = $pop{$key}[2];
        my %T; my %PBS;
        if($li[$id1]==1){
            $T{$pop1}{$pop2}=-log(1-$maxfst{$pop1}{$pop2});
        }else{
            $T{$pop1}{$pop2}=-log(1-$li[$id1]);
            $T{$pop1}{$pop2}=0 if($T{$pop1}{$pop2} eq "-0");
        }
        if($li[$id2]==1){
            $T{$pop1}{$pop3}=-log(1-$maxfst{$pop1}{$pop3});
        }else{
            $T{$pop1}{$pop3}=-log(1-$li[$id2]);
            $T{$pop1}{$pop3}=0 if($T{$pop1}{$pop3} eq "-0");
        }
        if($li[$id3]==1){
            $T{$pop2}{$pop3}=-log(1-$maxfst{$pop2}{$pop3});
        }else{
            $T{$pop2}{$pop3}=-log(1-$li[$id3]);
            $T{$pop2}{$pop3}=0 if($T{$pop2}{$pop3} eq "-0");
        }
        $PBS{$pop1}=($T{$pop1}{$pop2}+$T{$pop1}{$pop3}-$T{$pop2}{$pop3})/2;
        $PBS{$pop2}=($T{$pop1}{$pop2}+$T{$pop2}{$pop3}-$T{$pop1}{$pop3})/2;
        $PBS{$pop3}=($T{$pop1}{$pop3}+$T{$pop2}{$pop3}-$T{$pop1}{$pop2})/2;
        print OUT "\t$PBS{$pop1}\t$PBS{$pop2}\t$PBS{$pop3}";
    }
    print OUT "\n";
}
