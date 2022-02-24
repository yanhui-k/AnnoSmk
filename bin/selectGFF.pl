#!/usr/bin/perl -w
open LIST, $ARGV[0] or die;
open IN, $ARGV[1] or die;

while(<LIST>)
{
	chomp;
    s/\s.*//;
	$name{$_} = 1;
}

while(<IN>)
{
    my $name = "";
    my $name2 = "";
	if(/ID=(.*?);Name/)
    {
        $name=$1;
    }
	if(/Parent=(.*)/)
    {
        $name2=$1;
    }
    #@e=split;
    #if(exists $name{$e[0]})
	#{
	#	print $_;
	#}
	if(exists $name{$name} or exists $name{$name2})
	{
		print $_;
	}
}
