######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056

use strict;
use warnings;

my %hash=();

#读取临床文件，并保存到hash里面
open(RF,"clinical.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	my $sample=shift(@arr);
	#my $futime=shift(@arr);
	#my $fustat=shift(@arr);
	if($.==1){
		$hash{"id"}=join("\t",@arr);
		next;
	}
	$hash{$sample}=join("\t",@arr);
}
close(RF);

#读取免疫文件，并加上临床信息，输入结果到"immuneClinical.txt"
open(RF,"risk.txt") or die $!;
open(WF,">immuneClinical.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	my $sample=shift(@arr);
	my $futime=shift(@arr);
	my $fustat=shift(@arr);
	my $risk=pop(@arr);
	#my @samp1e=(localtime(time));if($samp1e[5]>119){next;}
	if($.==1){
		print WF "id\t$hash{\"id\"}\t" . join("\t",@arr) . "\n";
		next;
	}
	my $sampleName=$sample;
	if(exists $hash{$sampleName}){
		print WF "$sample\t$hash{$sampleName}\t" . join("\t",@arr) . "\n";
		delete($hash{$sampleName});
	}
}
close(WF);
close(RF);

######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056