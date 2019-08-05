#!/usr/bin/perl
use strict;

my $dir="/well/mccarthy/production/chip-seq/data/HNF4A_HNF1A_Singapore";
my $sample_file="$dir/sample_key.txt";

open(IN,'<',$sample_file) or die;
my $fline = <IN>;
while (my $line = <IN>){
	chomp $line;
	my @samples=split("\t",$line);
	if ($samples[1] =~ m/(.+)_input/){
		print $samples[1],"\t", $1,"\n";
		my $chip_id = `grep $1 $sample_file | grep -v input | awk '{printf \$1}'`;
		my $chip_sample = "$chip_id.fastq.gz";

		my $input_sample = "$samples[0].fastq.gz";
		#print $1,"\t",$chip_id,"\t", $chip_sample,"\t",$input_sample,"\n";
		my $aquas_command = "python /well/got2d/agata/bin/TF_chipseq_pipeline/chipseq.py TF idr -species hg19 ";
		$aquas_command .= "-peak_caller macs2 -fastq1 $dir/concat_raw_fastq/$chip_sample ";
		$aquas_command .= "-ctl_fastq1 $dir/concat_raw_fastq/$input_sample --out-dir $dir/aquas/$1 --title $1\n";
		print $aquas_command,"\n";
		system $aquas_command;	
	}
}
close IN;

