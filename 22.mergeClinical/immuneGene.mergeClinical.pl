######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺2749657388@qq.com
######����΢��: 18520221056

use strict;
use warnings;

my %hash=();

#��ȡ�ٴ��ļ��������浽hash����
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

#��ȡ�����ļ����������ٴ���Ϣ����������"immuneClinical.txt"
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
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺2749657388@qq.com
######����΢��: 18520221056