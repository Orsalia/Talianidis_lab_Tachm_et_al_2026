#!/bin/perl
#use warnings;
use strict;

my $inputfile=$ARGV[0];
my $inputfile2=$ARGV[1];
my $inputfile3=$ARGV[2];
open (INPUT_FILE,$inputfile);

open (file2,">$inputfile3.csv") or die "$!\n";

my %matrix1;
my %matrix2;
my %matrix3;

my @matrix1;
my @matrix2;
my @matrix3;
my @matrix4;
my @matrix5;
my @matrix6;
my @matrix7;
my @matrix8;
my @matrix9;
my @matrix10;
my @matrix11;
my @matrix12;
my @matrix13;
my @keep=0;
my @temp2=0;

my $temp1=0;
my $temp2=0;
my $temp3=0;

my $var1=0;
my $var2=0;
my $var3=0;my $var3b=0;my $var3c=0;my $var3d=0;
my $var4=0;
my $var5=0;my $var5b=0;my $var5c=0;my $var5d=0;
my $var6=0;
my $var7=1;
my $var8=0;my $var8b=0; my $var8c=0;my $var8d=0;
my $var9=0;
my $var10=0;
my $signal1=0;my $diff=0;my $keepvar=0;my $keephold=0;my $signal2=0;

my $test1=0;
my $test2=0;

my $varkeep1=0; my $varkeep2=0;

while(my $line1=<INPUT_FILE>)  
 {   
     @{$matrix1{$temp1}}=split(/\,/,$line1);
     $var1=${$matrix1{$temp1}}[1]+20000;
     if ((${$matrix1{$temp1}}[5]=~/"Differential"/)||(${$matrix1{$temp1}}[5]=~/"Shifted"/))
    { print file2 "$inputfile2\t${$matrix1{$temp1}}[1]\t$var1\t${$matrix1{$temp1}}[2]\t${$matrix1{$temp1}}[3]\t${$matrix1{$temp1}}[4]\t${$matrix1{$temp1}}[5]\t${$matrix1{$temp1}}[6]\t${$matrix1{$temp1}}[7]\n";
    }
  }
 


