#!/usr/bin/perl -w
#
# Clusterize annotated results ".transitions.annotated" (with score and probabilities) with mcl
#
use strict;
#use Getopt::Long;

my $usage="$0 cvs_file transition_file\n";

my $CUT_OFF=1e-3; #min probability
my $Inf=2.5; #inflation parameter for mcl (eg. 1.4, 2, 4,6) - 2=default
my $skip_negatives = 1;

my($csv,$transitions)=(shift,shift);
my $org=$1 if ($transitions=~/(\S+).transitions.*annotated/) or die ("not a .transitions.annotated file");
open T, "<$transitions" or die $!;

open MCL, ">$org.mcl_input" or die $!;

# write mcl input (cluster1 cluster2 score)
my %Score;
foreach (<T>){
	my @A=split /\t/;
	next if($A[8]!~/\d/ or $A[8]<0 and $skip_negatives);
	if($A[10]=~/\d/ and $A[10]<$CUT_OFF){
		$A[10]=1e-300 if ($A[10]<1e-300);
		printf MCL "$A[0]\t$A[1]\t%s\n", -(log($A[10])/log(10)); #negative base10 log of the expectation value (NLE)
		$Score{$A[0]}{$A[1]}=$Score{$A[1]}{$A[0]}=$A[8];
	}
}
close MCL; close T;

#MCL clustering 
system ("mcl $org.mcl_input -I $Inf --abc -o $org.clusters");
printf STDERR "mcl input in $org.mcl_input\n";
printf STDERR "written $org.clusters\n";


#annotate with descriptions
open CSV, "<$csv" or die $!;
my (%D,%G);
foreach (<CSV>){
  my @A=split /\t/;
    my $genes=grep /\S/, @A[2..$#A]; # number of genes; 
    $G{$A[0]}=$genes;
    $D{$A[0]}=$A[1]; #description
}
close CSV;

open CLU, "<$org.clusters" or die $!;
open CLU_ANN, ">$org.clusters.annotated" or die $!;
my $c_count=0;
foreach (<CLU>){
  chomp;
  my @B=split;
  next if ($#B<1); #skip mono clusters
  $c_count++;
  my $t_score=0;
  for (my $i=0;defined $B[$i+1];$i++) { 
     for (my $t=$i+1;defined $B[$t];$t++) { 
		if (defined $Score{$B[$i]}{$B[$t]}){ $t_score+=$Score{$B[$i]}{$B[$t]} } else {$t_score+=0.1} 
	}
  }
  printf CLU_ANN "%.5f %s \(%s\): ", $t_score/( (scalar @B * (scalar @B -1))/2 ),$G{$B[0]}, scalar @B; 
  foreach (@B){printf CLU_ANN "%s:%s ### ", $_,$D{$_}  } 
  printf CLU_ANN "\n";
}
close CLU; close CLU_ANN;
printf STDERR "written $c_count modules in $org.clusters.annotated\n";
