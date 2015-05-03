#!/usr/bin/perl
open(infile,"output.txt");
while(<infile>){
    my $current = $_;
#my $current = <infile>;
    chomp($current);
    ($dir1, $dir2, $dis) = split(/\s+/,$current);
    open(infile1,"./processed/$dir1\/$dir2\/list.dat");
    while(<infile1>){
	my $curr = $_;
#my $curr = <infile1>;
	chomp($curr);
	($num,$time,$hv,$dis) = split(/\s+/,$curr);
	my $filename;
	if ($num<10){
	    $filename = "tek000$num\.csv";
	}elsif ($num<100){
	    $filename = "tek00$num\.csv";
	}elsif ($num<1000){
	    $filename = "tek0$num\.csv";
	}elsif ($num<10000){
	    $filename = "tek$num\.csv";
	}
	my $flag= 0;
	open(infile2,"./raw/$dir1\/$dir2\/$filename");
	open(outfile,">./processed/$dir1\/$dir2\/$num\.dat");
	print outfile "$time $hv $dis\n";
#print "./processed/$dir1\/$dir2\/$num\/.dat $time $hv $dis\n";
	while(<infile2>){
	    $cu = $_;
	    chomp($cu);
	    @temp = split(/,/,$cu);
	    if ($flag==1){
		print outfile "$temp[0] $temp[1] $temp[2] $temp[3]\n";
	    }
	    if ($temp[0] =~ m/TIME/){
		$flag = 1;
	    }
	}
close(outfile);
close(infile2);
    }
}
