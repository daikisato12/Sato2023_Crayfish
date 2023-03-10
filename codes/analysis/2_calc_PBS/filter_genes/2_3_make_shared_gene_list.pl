 use warnings;
use strict;

#my @bp = ("1bp", "100bp", "500bp", "1kbp", "10kbp");
#my @bp = ("1kbp", "10kbp");
#my @pop = ("Aomori_kamakura2_tonden2", "Aomori_okinawa2_tonden2", "kamakura2_okinawa2_tonden2");
#my %pop = ("Aomori_kamakura2_tonden2" => "AKT", "Aomori_okinawa2_tonden2" => "AOT", "kamakura2_okinawa2_tonden2" => "KOT");
#my @threshold = ("1", "0.1");
#my @threshold = ("0.1");
#my $cov = "mincov4_maxcov16";
my $dir = "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/manuscript_tmp/data/HAF/PBS/analyzed_data";

my $resblast = '/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/data/blastp/GCF_020424385.1_ASM2042438v2_translated_cds_longest_dros.best_hits_gn.txt';
my $dros_gtf = '/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/Data/Ref/Ref_Drosophila/Drosophila_melanogaster.BDGP6.32.104.gtf';
my $deg = '/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/data/RNA-seq/unique_mapping/DEGs/RNAseq_DEG2_CrayfishID.txt';
my $ctrg = '/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/figures/tmp/WGCNA_Module.tsv';

my %drospp_length;
open(IN, $dros_gtf);
while(<IN>){
chomp;
	next if($_ =~ /^#/);
	my @li = split(/\t/,$_);
	next if($li[2] ne "CDS");
	if($li[8] =~ /gene_id \"(\S+)\"\;/){
		my $flyid = $1;
		my $length = $li[4] - $li[3] + 1;
		if(defined $drospp_length{$flyid}){
			$drospp_length{$flyid} = $drospp_length{$flyid} + $length;
		}else{
			$drospp_length{$flyid} = $length;
		}
	}
}
close(IN);

my %hash;
open(IN, $resblast);
while(<IN>){
chomp;
	my @li = split(/\t/,$_);
	my $gene = $li[0];
	my $flyid = $li[1];
	if(defined $hash{$gene}){
		my $flyid_old = $hash{$gene};
		if($drospp_length{$flyid} > $drospp_length{$flyid_old}){
			$hash{$gene} = $flyid;
		}
	}else{
		$hash{$gene} = $flyid;
	}
}
close(IN);

my %deg;
open(IN, $deg);
while(<IN>){
	chomp;
	$deg{$_} = 1;
}

my %ctrg;
open(IN, $ctrg);
while(<IN>){
	chomp;
    my @li = split(/\t/,$_);
    next if($_ =~ /^GeneID/);
	$ctrg{$li[0]} = 1;
}

#for my $bp (@bp){
#    for my $threshold (@threshold){
#        for my $pop (@pop){
my $bp = "5kbp_step1kbp";
my $threshold = "0.1";
my $infile = "$dir\/HAF_PBS_$bp\-top$threshold\_2combsharedloci_gene.tsv";
my $outfile = "$dir\/HAF_PBS_$bp\-top$threshold\_2combsharedloci_gene_info.tsv";
my %deg3;
my %ctrg3;
open(OUT, ">$outfile");
open(IN, $infile);
while(<IN>){
    chomp;
    my @li = split(/\t/, $_);
    if($_=~/^Chr/){
        print OUT "$li[0]\t$li[6]\t$li[1]\t$li[2]\t$li[7]\tGenes_FlyID\tDEGs\tCTRGs\t$li[4]\t$li[5]\t$li[3]\n";
    }else{
        print OUT "$li[0]\t$li[6]\t$li[1]\t$li[2]\t$li[7]\t";
        my @genes = split(/\, /, $li[7]);
        my %deg2;
        my %ctrg2;
        my $t = 0;
        for my $gene (@genes){
            my $gene_tmp = $gene;
            $gene_tmp = "$gene \($hash{$gene}\)" if(defined $hash{$gene});
            if(defined $deg{$gene}){
                $deg2{$gene} = "$gene_tmp";
                $deg3{$gene} = "$gene_tmp";
            }
            if(defined $ctrg{$gene}){
                $ctrg2{$gene} = "$gene_tmp";
                $ctrg3{$gene} = "$gene_tmp";
            }
            if($t == 0){
                print OUT "$gene_tmp";
            }else{
                print OUT ", $gene_tmp";
            }
            $t++;
        }
        $t = 0;
        print OUT "\t";
        foreach my $key (keys %deg2){
            if($t == 0){
                print OUT "$deg2{$key}";
            }else{
                print OUT ", $deg2{$key}";
            }
            $t++;
        }
        $t = 0;
        print OUT "\t";
        foreach my $key (keys %ctrg2){
            if($t == 0){
                print OUT "$ctrg2{$key}";
            }else{
                print OUT ", $ctrg2{$key}";
            }
            $t++;
        }
        print OUT "\t$li[4]\t$li[5]\t$li[3]\n";
    }
}
close(OUT);

open(OUT2, ">$dir\/HAF_PBS_$bp\-top$threshold\_2combsharedloci_DEGs_genelist.txt");
foreach my $key (keys %deg3){
    print OUT2 "$key\n";
}
open(OUT3, ">$dir\/HAF_PBS_$bp\-top$threshold\_2combsharedloci_CTRGs_genelist.txt");
foreach my $key (keys %ctrg3){
    print OUT3 "$key\n";
}
