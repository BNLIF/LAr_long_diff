#!/usr/bin/perl

open(infile,"data.txt");
open(outfile,">output.txt");
while(<infile>){
    my $current = $_;
    chomp($current);
    @temp = split(/\s+/,$current);
    if (-e "./processed/$temp[0]"){
    }else{
	system("mkdir ./processed/$temp[0]");
    }
    for (my $i = 2; $i!=@temp;$i++){
	if (-e "./processed/$temp[0]/$temp[$i]"){
	}else{
	    system("mkdir ./processed/$temp[0]/$temp[$i]");
	}
	print outfile "$temp[0] $temp[$i] $temp[1]\n";
    }
}
close(infile);
