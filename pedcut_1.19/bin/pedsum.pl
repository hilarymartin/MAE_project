#! /usr/bin/perl 
# use strict;
# this script summarize a given set of pedigrees

&pedsum(@ARGV);
sub pedsum {
  die 
"Died because of misuse of the subroutine pedsum.
Usage: pedsum (inped,output);
In total of 2 arguments should be provided.
The pedsum sub recognized these arguments: @_\n" if (@_<2);
  my $inped = $_[0]; my $output = $_[1]; 
  my (@pids);
  my (%pid_N,%pid_f,%pid_nf,%pid_m,%pid_fm,%pid_aff,%pid_unaff,%pid_g,%pid_bit);
  my $g=0;
  open (IN, "$inped")||die"can not open input pedigree file $inped,$!\n";
  while (my $line = <IN>){
    chomp $line;
    my @line = split/\s+/,$line;
    if (!$line[0]){shift (@line);}
    if(@line==1){@line=split/,/,$line[0]}
    if (@line>4){
      push (@pids,$line[0]);
      $pid_N{$line[0]}++;
      if($line[2]==0 and $line[3]==0){$pid_f{$line[0]}++;}
      else{$pid_nf{$line[0]}++;}
      if ($line[4]==1){$pid_m{$line[0]}++;}
      elsif ($line[4]==2){$pid_fm{$line[0]}++;}
      else {
print "Input error: sex of pedi $line[0] indiv $line[1] is $line[4]\n";
      }
      if ($line[5]==2){$pid_aff{$line[0]}++;}
      elsif ($line[5]==1){$pid_unaff{$line[0]}++;}
      else {}
      if ($line[5] != 0){$pid_g{$line[0]}++;}
      if ($line[6] != 0){$pid_d{$line[0]}++;}
    }
  }@pids = sort{$a<=>$b}keys %{{map{$_,1} @pids}};
  foreach my $p(@pids){
    if (!$pid_N{$p}){$pid_N{$p}=0;}
    if (!$pid_aff{$p}){$pid_aff{$p}=0;}
    if (!$pid_g{$p}){$pid_g{$p}=0;}
    if (!$pid_unaff{$p}){$pid_unaff{$p}=0;}
    if (!$pid_bit{$p}){$pid_bit{$p}=0;}
    if (!$pid_m{$p}){$pid_m{$p}=0;}
    if (!$pid_f{$p}){$pid_f{$p}=0;}
    if (!$pid_d{$p}){$pid_d{$p}=0;}
  }
  open (OUT,">$output")||die"$!\n";select (OUT);
  print "pid,N,n_SOI,n_aff,n_unaff,bits,n_male,n_founder,n_duplicate\n";
  foreach  my $p(@pids){
    $pid_bit{$p}=$pid_nf{$p}*2-$pid_f{$p};
print"$p,$pid_N{$p},$pid_g{$p},$pid_aff{$p},$pid_unaff{$p},$pid_bit{$p},$pid_m{$p},$pid_f{$p},$pid_d{$p}\n";
  }close (OUT); select(STDOUT);
  print "Created Output $output\n";
}
