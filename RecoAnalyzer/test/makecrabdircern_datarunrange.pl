#!/usr/local/bin/perl

#use strict;
#use warnings;

($#ARGV==6) || die"usage: [ DataSetName   TagName trigProcName runMin-runMax crab_scheduler (condor for cmslpc or glite in general) lumiPerJob cerncastorDir]\n";


my $fname = $ARGV[0]; 
my $dataset = $ARGV[0];
my $tagname = $ARGV[1];
my $hltname = $ARGV[2];
my $runMinToMax = $ARGV[3];
my $crabscheduler = $ARGV[4];
my $lumiperjob = $ARGV[5];
my $cerncastordir = $ARGV[6];


$fname =~ s/\///g; 

print $dataset , "\n"; 

print $fname, "dd\n"; 
    
my $newdir = "crab_analysis/${fname}" ;

if( -d $newdir) {
    if( -d "$newdir/$runMinToMax"){
	print "directory already createdd.. $newdir/$runMinToMax \n";
	exit;
    }else{
	print "make new directory $newdir/$runMinToMax\n";
	system("mkdir $newdir/$runMinToMax");
    }

}else{
    print "making directory $newdir ..\n"; 
    system("mkdir $newdir");
    if( -d "$newdir/$runMinToMax"){
	print "directory already createdd.. $newdir/$runMinToMax \n";
	exit;
    }else{
	print "make new directoy $newdir/$runMinToMax\n";
	system("mkdir $newdir/$runMinToMax");
    }
}

my $newdir = "crab_analysis/${fname}/$runMinToMax" ;

system("cp crab_temp/analyzertemp_runranges.py $newdir/test.py");
system("cp crab_temp/crab.cfg $newdir/crab.cfg");

system("perl ReplaceString.pl CRAB_Scheduler $crabscheduler $newdir/crab.cfg");
system("perl ReplaceString.pl runMinToMax $runMinToMax $newdir/crab.cfg");
system("perl ReplaceString.pl maindataname $fname $newdir/crab.cfg");
system("perl ReplaceString.pl DataSetName $dataset $newdir/crab.cfg");
system("perl ReplaceString.pl cern_castor_dir $cerncastordir $newdir/crab.cfg");
system("perl ReplaceString.pl LumiPerJob $lumiperjob $newdir/crab.cfg");

system("perl ReplaceString.pl TagName $tagname $newdir/test.py");
system("perl ReplaceString.pl maindataname $fname $newdir/test.py");
system("perl ReplaceString.pl hltProcName $hltname $newdir/test.py");
system("perl ReplaceString.pl runMinToMax $runMinToMax $newdir/test.py");


print "crab directory ready...\n"; 

#system("cd $newdir; pwd; source /uscmst1/prod/grid/gLite_SL5_CRAB_27x.sh; eval `scramv1 runtime -sh`; which cmsRun; source /uscmst1/prod/grid/CRAB/crab.sh; which crab; crab -create -submit");
#print "crab job submitted ..\n"; 

