use warnings;
use strict;

my $outfile = 'gene2accession_Hsap';

#download gene2accession if outfile is not present
if(! -e $outfile){
	my $file = 'gene2accession';
	my $url = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/'.$file.'.gz';
	system(`curl -o $file.gz $url`);
	system(`gunzip $file.gz`);

	#extract Homo sapiens annotation from gene2accession
	open OUT, '>', 'gene2accession_Hsap' or die $!;
	open FILE, $file or die $!;
	while(<FILE>){
		print OUT $_ if $. == 1;
		next unless $_=~ m/^9606\s+/g;
		print OUT $_;
	}
	close FILE;
	close OUT;
	unlink($file);
}
