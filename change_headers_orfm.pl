#!/usr/bin/perl
use warnings;
use strict;

my $forward_orfs = shift(@ARGV);
my $reverse_orfs = shift(@ARGV);

my $n = 0;
open(FORWARD,"$forward_orfs");
while(my $x = <FORWARD>){
	chomp($x);
	if($x =~ m/^>/){
		my $pre = substr($x,1);
		my $new = ">$n\_forward\_$pre";
		print "$new\n";
		$n++;
	}else{
		print "$x\n";
	}
}
close FORWARD;
open(REVERSE,"$reverse_orfs");
while(my $x = <REVERSE>){
        chomp($x);
        if($x =~ m/^>/){
                my $pre = substr($x,1);
                my $new = ">$n\_reverse\_$pre";
                print "$new\n";
                $n++;
        }else{
                print "$x\n";
        }
}
close REVERSE;
