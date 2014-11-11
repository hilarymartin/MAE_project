use warnings;
use strict;

my $ped= $ARGV[0];#ped file with original genotype data used for phasing
my $map= $ARGV[1];#map file with original genotype data used for phasing
my $errors = $ARGV[2];#generr file from duoHMM
my $newped = $ARGV[3];#new ped file for output
my $individual_count = $ARGV[4]; #output file for count of errors per individual
my $snps_count=$ARGV[5]; #output file for count of errors per SNP
my %positions;
open(MAP,$map);
my $count=0;
while(<MAP>){
    chomp $_;
    $count++;
    my @d=split("\t",$_);

    $positions{$d[3]} = $count;
}
close(MAP);

my %err;
my %snp_err;
open(ERRORS,$errors);
<ERRORS>;
#ID father mother RSID1 POS PROB_ERROR
while(<ERRORS>){
    chomp $_;
    my ($ID,$father,$mother,$chr,$pos,$prob_error)=split("\t",$_);
#most SNPs in these files have prob_error >0.9, so will just set all of them to missing
    $err{$ID}{$pos}++;
    if($father!~/^0/){
	$err{$father}{$pos}++;
    }
    if($mother!~/^0/){
	$err{$mother}{$pos}++;
    }
    $snp_err{$pos}++;
}
close(ERRORS);

open(PED,$ped);
open(NEWPED,">".$newped);
while(<PED>){
    chomp $_;
    my @d=split("\t",$_);
    if(exists $err{$d[1]}){
	foreach my $pos(keys %{$err{$d[1]}}){
#	    print $d[1],",\t,",$snp,",\n";
	    if(not exists $positions{$pos}){
		print $pos,"\n";
	    }
	    my $index=$positions{$pos}+5;
	    $d[$index]="0 0";
	}
    }
    print NEWPED join("\t",@d),"\n";
}
close(PED);

open(SNPERR,">".$snps_count);
foreach my $snp(keys %snp_err){
    print SNPERR $snp,"\t",$snp_err{$snp},"\n";
}
close(SNPERR);

open(INDIVERR,">".$individual_count);
foreach my $sample(keys %err){
    my $n=scalar keys %{$err{$sample}};
    print INDIVERR $sample,"\t",$n,"\n";
}
close(INDIVERR);


    
