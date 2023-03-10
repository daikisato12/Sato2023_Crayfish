use warnings;
use strict;
use File::Basename 'basename', 'dirname';
use File::Find::Rule;

my $resblast = '/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/data/blastp/GCF_020424385.1_ASM2042438v2_translated_cds_longest_dros.best_hits_gn.txt';
my $dros_gtf = '/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/Data/Ref/Ref_Drosophila/Drosophila_melanogaster.BDGP6.32.104.gtf';
my $gtf = '/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/Data/Ref/Xu2021/v2/GCF_020424385.1_ASM2042438v2_genomic.gtf';
my $fa = '/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/Data/Ref/Xu2021/v2/GCF_020424385.1_ASM2042438v2_genomic.fna';

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

my %fly;
open(IN, $resblast);
while(<IN>){
chomp;
	my @li = split(/\t/,$_);
	my $gene = $li[0];
	my $flyid = $li[1];
	if(defined $fly{$gene}){
		my $flyid_old = $fly{$gene};
		if($drospp_length{$flyid} > $drospp_length{$flyid_old}){
			$fly{$gene} = $flyid;
		}
	}else{
		$fly{$gene} = $flyid;
	}
}
close(IN);


my %index;
my %index2;
open(IN, $fa);
while(<IN>){
chomp;
    if($_=~/^>/){
        if($_=~/>(\S+) Procambarus clarkii strain Jiangsu (\S+), ASM2042438v2, whole genome shotgun sequence/ or $_=~/>(\S+) Procambarus clarkii strain Jiangsu linkage group (\S+), ASM2042438v2, whole genome shotgun sequence/){
            $index{$1} = $2;
            $index2{$2} = $1;
#            print "$1\t$2\n";
        }
    }
}

#print(@files);

my %gene;
open(IN, $gtf);
while(<IN>){
    chomp;
    next if($_=~/^#/);
    my @li=split(/\t/, $_);
#    print "$li[8]\n";
    my $chr = $index{$li[0]} if(defined $index{$li[0]});
    next unless(defined $chr);
    if($li[2] eq "gene"){
        if($li[8] =~ /gene_id \"(\S+)\"\; transcript_id/){
            my $tmp = $1;
#           print "$tmp\n";
#           $tmp = "undefined" unless(defined $tmp);
            @{$gene{$chr}{$tmp}}=($li[3],$li[4]); #start, end
        }
    }
}
close(IN);

my $dir = "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/manuscript_tmp/data/HAF/PBS/analyzed_data";
my @files = File::Find::Rule->in($dir);
#print "@files\n";

#=pod
foreach my $file (@files){ #$file: PBS_Aomori_kamakura2_tonden2_Tonden_annotatedCHR_mincov4_maxcov16_1kbp-top0.1.tsv
    next if($file !~ /top0.1_2combsharedloci.tsv/);
    print "$file\n";
    my %genes;
    my $outdir = dirname $file;
    my $outfile = basename $file;
    substr($outfile, -4, 0, '_gene');
    open(OUT, ">$outdir/$outfile");

    my %hash;
    my %genes_all;
    open(IN, $file);
    while(<IN>){
        chomp;
        print OUT "$_\tCHR_NEW\tGenes\n" if($_=~/^Chr/);
        my @li = split(/\t/, $_);
        my @li2 = split(/\_/, $li[0]);
        my $chr = $li2[0];
        my $t = 0;
        my $genes;
        my $chr_new = "$chr\_$index2{$chr}" if(defined $index2{$chr});
        next unless(defined $chr_new);
#        $hash{$t}=($li[1],$li[2]); #chromosome, coordinate
        next unless(defined $gene{$chr});
        foreach my $key (keys %{$gene{$chr}}){
#            print "$key\n";
            next if(int($gene{$chr}{$key}[0])-50000 > $li[1] or int($gene{$chr}{$key}[1])+50000 < $li[1]); #search genes within 50kb up/downstream
            $genes="$genes, $key" if($t>0);
            $genes="$key" if($t==0);
            $t++;
#            print "$fly{$key}\n";
            $genes_all{$key} = $fly{$key};
        }
        if($t>0){
            print OUT "$_\t$chr_new\t$genes\n";
        }else{
            print OUT "$_\t$chr_new\tNA\n";
        }
    }
    close(IN);
    close(OUT);

    $outfile = basename $file;
    substr($outfile, -4, 4, '_genelist.txt');
    open(OUT, ">$outdir/$outfile");
    my $outfile2 = basename $file;
    substr($outfile2, -4, 4, '_genelist_FlyID.tsv');
    open(OUT2, ">$outdir/$outfile2");
    print OUT2 "CrayfishID\tFlyID\n";
    foreach my $key (sort {$a cmp $b} keys %genes_all){
        print OUT "$key\n";
        if(defined $genes_all{$key}){
            print OUT2 "$key\t$genes_all{$key}\n";
        }else{
            print OUT2 "$key\t\n";
        }
    }
}
#=cut