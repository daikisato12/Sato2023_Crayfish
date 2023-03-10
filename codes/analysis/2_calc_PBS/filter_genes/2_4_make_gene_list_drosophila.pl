use warnings;
use strict;
use File::Find::Rule;

# This code converts crayfish geneIDs into Drosphila. 
# If a given crayfish gene was blasted to multiple Drosophila genes with the same scores, genes The longest CDS sequences were selected.

#my $dir = '/home/daikisato177/Projects/Crayfish/Analysis/Result/220519';
#my $all = '/home/daikisato177/Projects/Crayfish/Data_2022/RNA-seq/RNAseq_analyzed_CrayfishID.txt';
my $resblast = '/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/data/blastp/GCF_020424385.1_ASM2042438v2_translated_cds_longest_dros.best_hits_gn.txt';
my $gtf = '/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/Data/Ref/Xu2021/v2/GCF_020424385.1_ASM2042438v2_genomic.gtf';
my $rnaseq = '/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/data/RNA-seq/rawdata/crayfish_hisat2_stringtie_tpm20220519_LOC.tsv';
my $dros_gtf = '/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/Data/Ref/Ref_Drosophila/Drosophila_melanogaster.BDGP6.32.104.gtf';
my $lglistfile =  '/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/data/blastp/LG_index.txt'; #'/home/daikisato177/Data/crayfish/Ref_Xu2021/v2/LG_index.txt';

my %lglist;
open(IN, $lglistfile);
while(<IN>){
chomp;
	$lglist{$_} = 1;
}

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
#	my $fly_ppid = $li[2];
	if(defined $hash{$gene}){
		my $flyid_old = $hash{$gene};
		if($drospp_length{$flyid} > $drospp_length{$flyid_old}){
			$hash{$gene} = $flyid;
		}else{
			print "$flyid_old\: $drospp_length{$flyid_old}\t$flyid\: $drospp_length{$flyid}\n";
		}
	}else{
		$hash{$gene} = $flyid;
	}
#	push @{$hash{$gene}}, $flyid;
}
close(IN);

my %flyid_all;
open(OUT, ">/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/data/blastp/Allgenes_LG_crayfish_newID.txt");
open(IN, $gtf);
while(<IN>){
chomp;
	my @li = split(/\t/, $_);
	next unless(defined $lglist{$li[0]});
	next unless($li[2] eq "gene");
	if($li[8] =~ /gene_id \"(\S+)\"\; transcript_id/){
		my $gene = $1;
		print OUT "$gene\n";
		if(defined $hash{$gene}){
			my $flyid = $hash{$gene};
			$flyid_all{$flyid} = 1;
		}
#		for my $key (@{$hash{$gene}}){
#			$flyid_all{$key} = 1;
#		}
	}
}
close(IN);
close(OUT);

my %rnaseq;
my %rnaseq_fly;
open(IN, $rnaseq);
while(<IN>){
chomp;
	my @li = split(/\t/, $_);
	next if($_=~/^ID/);
	$rnaseq{$li[1]} = 1;
	if(defined $hash{$li[1]}){
		my $flygene = $hash{$li[1]};
		$rnaseq_fly{$flygene} = 1;
	}
}
close(IN);

open(OUT, ">/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/data/RNA-seq/rawdata/RNAseq_analyzed_CrayfishID.txt");
foreach my $key (keys %rnaseq){
	print OUT "$key\n";
}
close(OUT);
open(OUT2, ">/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/data/RNA-seq/rawdata/RNAseq_analyzed_CrayfishID_flyID.txt");
foreach my $key (keys %rnaseq_fly){
	print OUT2 "$key\n";
}
close(OUT2);

open(OUT2, ">Allgenes_LG_crayfish_newID_flyID.txt");
foreach my $key (sort {$a cmp $b} keys %flyid_all){
	print OUT2 "$key\n";
}
close(OUT2);

my $dir = "/Users/dsato/Dropbox/研究室/投稿論文/2020/ザリガニ論文/ザリガニ論文_Sato_Makino/data/PBS/";
my @files = File::Find::Rule->in($dir);
foreach my $file (@files) {
	next if($file !~ /sharedgene.txt/);
	my $outfile = $file;
	substr($outfile, -4, 4,'_flyID.txt');
	open(OUT, ">$outfile");

	my %hash2;
	open(IN, $file);
	while(<IN>){
		chomp;
		my $gene = $_;
		next unless(defined $hash{$gene});
		my $flyid = $hash{$gene};
		$hash2{$flyid} = 1;
#		for my $key (@{$hash{$gene}}){
#			print "$key\n";
#			$hash2{$key} = 1;
#		}
	}
	close(IN);

	foreach my $key (sort {$a cmp $b} keys %hash2){
#		print OUT "$key\n";
		print OUT "$key\n";
#		print OUT "$flyid{$key}\n";
	}	
}
