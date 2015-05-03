#!/usr/bin/perl
open(infile,"thomas.data");
while(<infile>){
    ($dir,$filename,$dis) = split(/\s+/,$_);
    print "$dir, $filename, $dis";
    chomp($dis);
    open(infile1,"./raw/$dir\/$filename");
    open(outfile,">./processed/thomas/$filename");
    @temp = split(/\s+/,<infile1>);
    my $num = scalar(@temp) -1;
    print outfile "$num $dis\n";
    print outfile "@temp\n";
    while(<infile1>){
	print outfile "$_";
    }
    close(outfile);
}
