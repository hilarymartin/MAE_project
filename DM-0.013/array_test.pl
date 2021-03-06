#!perl
use Test::More tests => 3;
use Test::Files;
use DM 0.001003;

my $test_dir = 't/50-job_array.t.tmp';

my $dm = DM->new(
    "dryRun"   => 1,
    "numJobs"  => 1,
    outputFile => "$test_dir/output.log",
    );

my @prereqs = ( "$test_dir/prereq1", "$test_dir/prereq2", "$test_dir/prereq3" );
system( "touch " . join( ' ', @prereqs ) );

my @targets = ( "$test_dir/target1", "$test_dir/target2" );

my $jobArrayObject = $dm->startJobArray(
    globalTmpDir => $test_dir,
    target       => "$test_dir/target_array.flag"
    );
$dm->addJobArrayRule(
    target  => $targets[0],
    prereqs => $prereqs[0],
    command => "echo 'hi world 1' > $targets[0]",
    engineName => 'SGE'
    );
$dm->addJobArrayRule(
    target  => $targets[1],
    prereqs => [ @prereqs[ 1 .. 2 ] ],
    command => "echo 'hi world 2' > $targets[1]",
    engineName => 'SGE'
    );

$dm->endJobArray();

compare_ok( $jobArrayObject->targetsFile,
	    "$test_dir/targets.expected", "targets file was created correctly" );
compare_ok( $jobArrayObject->prereqsFile,
	    "$test_dir/prereqs.expected", "prereqs file was created correctly" );
compare_ok( $jobArrayObject->commandsFile,
	    "$test_dir/commands.expected", "commands file was created correctly" );




unlink @prereqs;
