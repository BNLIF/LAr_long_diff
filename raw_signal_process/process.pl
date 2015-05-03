#!/usr/bin/perl

open(infile,"output.txt");

while(<infile>){
    my $current = $_;
    chomp($current);
    ($dir1, $dir2, $dis) = split(/\s+/,$current);
    my @array;
    my $flag = 0;
    open(infile1,"./raw/$dir1\/$dir2\/Data\ Spec.txt");
    open(outfile,">./processed/$dir1\/$dir2\/list.dat");
    #print "./processed/$dir1\/$dir2\/list.dat\n";
    while(<infile1>){
	my $t1 = $_;
	chomp($t1);
	my @temp = split(/\s+/,$t1);
	if ($temp[1]=~m/^\d/){
	    if ($flag == 0){
		@array = @temp;
		$flag = 1;
	    }else{
		for (my $i=1;$i!=@temp;$i++){
		    print outfile "$temp[$i] $temp[0] $array[$i] $dis\n";
		}
	    }
	}
    }
}
