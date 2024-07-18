#! /usr/bin/perl

# Author: Laurent Facq 5/10/2021

# find . \( -name *.cpp -o -name *.hpp \) -type f -exec bash -c "cat {} | Scripts/perl/clean-cxx-head.pl `basename {}` > /tmp/.t ; diff {} /tmp/.t " \;

use File::Basename;

$comment=0;

$header= ".licence-header";

$head="";
{ local $/; undef $/; open(IN,$header) or die "cannot open header $header" ; $head=<IN>; close(IN); }

$filename = basename($ARGV[0]);
$head =~ s/%F%/$filename/;
shift @ARGV;

print $head;

while ($l=<>)
{
    if ($l =~ m/\@brief(:|\s)+$/) { next; } # empty brief ignored
    #if ($l =~ m/\@brief/) { print "// $l" ; next; } # keep brief 
    
    if ($l =~ m/\*\//) { $comment=0; next; } # comment / end
    if ($l =~ m/^\s*$/) { print ""; next; } # empty
    if ($l =~ m/^\s*\/\*/) { $comment=1; next; } # comment / start
    
    if ($comment) { next; }
    
    print "$l"; last; # header finished 
}

# print remaining content
while ($l=<>) { print $l; }
