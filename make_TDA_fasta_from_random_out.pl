#!/usr/bin/perl
use warnings;
use strict;

my $random_csv = shift(@ARGV);
my $TDA_out_prefix = shift(@ARGV);
my $nr_target_out_prefix = shift(@ARGV);
my $nr_decoy_out_prefix = shift(@ARGV);
my $max_number_per_split = shift(@ARGV); #30000000

open(TABLE,"$random_csv");
my $i = 1;
my $t_size = 0;
my $orfm_targets = 0;
my $sixgill_targets = 0;
my $d_size = 0;
my $final_splits_number = 0;
while(my $x = <TABLE>){
	chomp($x);
	my $splits_number = 0;
	if($x =~ m/seq_name,peptide,random_peptide/){
		next;
	}
	my ($name,$seq,$rand) = split(/\,/,$x);
	my $sixgill = "sixgill";
	$name =~ s/[\:\;\, \t\n]+//g;
	if($name =~ m/forward/ or $name =~ m/reverse/ or $name =~ m/rverse/){
		$sixgill = "orfM";
		$orfm_targets++;
	}else{
		$sixgill_targets++;
	}
	my $new_name = "$i\_$sixgill\_$name";
	until((($splits_number * $max_number_per_split)/$i) >= 1){
		$splits_number++;
	}
	if($splits_number > $final_splits_number){
		$final_splits_number = $splits_number;
	}
	open my $tda, ">>", "$TDA_out_prefix\_$splits_number\.faa";
        open my $target, ">>", "$nr_target_out_prefix\_$splits_number\.faa";
        open my $decoy, ">>", "$nr_decoy_out_prefix\_$splits_number\.faa";
        if(length($rand) > 1){
                print $tda ">target\_$new_name\n$seq\n>decoy\_$new_name\n$rand\n";
                print $target ">target\_$new_name\n$seq\n";
                print $decoy ">decoy\_$new_name\n$rand\n";
                $t_size++;
                $d_size++;
        }else{
                print $tda ">target\_$new_name\n$seq\n";
                print $target ">target\_$new_name\n$seq\n";
                $t_size++;
        }
        close $tda;
        close $target;
        close $decoy;
	if(($i % 1000000) == 0){
                print STDERR "Printed $i sequences\n";
        }
	$i++;
}
close TABLE;
	
print STDERR "Targets --> $t_size\norfM targets --> $orfm_targets\nSixgill targets --> $sixgill_targets\nDecoys --> $d_size\nSplits --> $final_splits_number\n";
