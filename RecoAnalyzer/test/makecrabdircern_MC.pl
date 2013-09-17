#!/usr/local/bin/perl

#use strict;
#use warnings;

($#ARGV==5) || die"usage: [ DataSetName TagName trigProcName crab_scheduler (condor for cmslpc or glite in general) eventPerJob cerncastorDir]\n";


my $fname = $ARGV[0]; 
my $dataset = $ARGV[0];
my $tagname = $ARGV[1];
my $hltname = $ARGV[2];
my $crabscheduler = $ARGV[3];
my $eventperjob = $ARGV[4];
my $cerncastordir = $ARGV[5];

$fname =~ s/\///g; 

print $dataset , "\n"; 

print $fname, "dd\n"; 
    
my $newdir = "crab_analysis/${fname}" ;

if( -d $newdir) {
    print "directory already createdd.. $newdir \n";
    exit;
}else{
    print "making directory $newdir ..\n"; 
    system("mkdir $newdir");
}


system("cp crab_temp/analyzertemp.py $newdir/test.py");
system("cp crab_temp/crab.cfgMC $newdir/crab.cfg");

system("perl ReplaceString.pl CRAB_Scheduler $crabscheduler $newdir/crab.cfg");
system("perl ReplaceString.pl maindataname $fname $newdir/crab.cfg");
system("perl ReplaceString.pl DataSetName $dataset $newdir/crab.cfg");
system("perl ReplaceString.pl cern_castor_dir $cerncastordir $newdir/crab.cfg");
system("perl ReplaceString.pl EvtPerJOB $eventperjob $newdir/crab.cfg");

system("perl ReplaceString.pl TagName $tagname $newdir/test.py");
system("perl ReplaceString.pl maindataname $fname $newdir/test.py");
system("perl ReplaceString.pl hltProcName $hltname $newdir/test.py");


print "crab directory ready...\n"; 

#system("cd $newdir; pwd; source /uscmst1/prod/grid/gLite_SL5_CRAB_27x.sh; eval `scramv1 runtime -sh`; which cmsRun; source /uscmst1/prod/grid/CRAB/crab.sh; which crab; crab -create -submit");
#print "crab job submitted ..\n"; 

