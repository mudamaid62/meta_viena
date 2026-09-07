#!/usr/bin/perl
use warnings;
use strict;

my $table_1 = shift(@ARGV);
my $table_2 = shift(@ARGV);
my $out_1 = shift(@ARGV);
my $out_2 = shift(@ARGV);
my $out_s = shift(@ARGV);
open(TABLE_1,"$table_1");

#query,target,pident,qcov,tcov,evalue,bits,qlen,tlen,alnlen,qseq

my %seqs_1;
my %seqs_2;
my %queries_1;
my %queries_2;
my %all_queries;

while(my $x = <TABLE_1>){
	chomp($x);
	my @array = split(/\t/,$x);
	my $name = shift(@array);
	my $query = $name;
	$name =~ s/\:/_/g;
	my $seq = pop(@array);
	$seqs_1{$name} = $seq;
	$queries_1{$name} = $query;
	$all_queries{$name} = $query;
}
close TABLE_1;
open(TABLE_2,"$table_2");

while(my $x = <TABLE_2>){
        chomp($x);
        my @array = split(/\t/,$x);
        my $name = shift(@array);
        my $query = $name;
        $name =~ s/\:/_/g;
        my $seq = pop(@array);
        $seqs_2{$name} = $seq;
        $queries_2{$name} = $query;
        $all_queries{$name} = $query;
}
close TABLE_2;
my $n1 = 0;
my $n2 = 0;
my $ns = 0;
my $pairs = 0;
open(OUT_1, ">$out_1");
open(OUT_2, ">$out_2");
open(OUT_S, ">$out_s");
foreach my $q(sort keys %all_queries){
	if(exists($queries_1{$q}) and exists($queries_2{$q})){
		my $len_1 = length($seqs_1{$q});
		my $len_2 = length($seqs_2{$q});
		my $q1 = "=" x $len_1;
		my $q2 = "=" x $len_2;
		print OUT_1 "\@$queries_1{$q}\n$seqs_1{$q}\n+\n$q1\n";
		print OUT_2 "\@$queries_2{$q}\n$seqs_2{$q}\n+\n$q2\n";
		$n1++;
		$n2++;
		$pairs++;
	}elsif(exists($queries_1{$q})){
		my $len_1 = length($seqs_1{$q});
		my $q1 = "=" x $len_1;
		print OUT_S "\@$queries_1{$q}\n$seqs_1{$q}\n+\n$q1\n";
		$ns++;
		$n1++;
	}elsif(exists($queries_2{$q})){
                my $len_2 = length($seqs_2{$q});
                my $q2 = "=" x $len_2;
                print OUT_S "\@$queries_2{$q}\n$seqs_2{$q}\n+\n$q2\n";
                $ns++;
		$n2++;
	}
}
my $head = "-" x 20;
print STDERR "$head\nPrinted $n1 forward reads\nPrinted $n2 reverse reads\nPrinted $pairs proper pairs\nPrinted $ns single reads\n$head\n";
