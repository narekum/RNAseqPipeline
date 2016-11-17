#! /usr/bin/perl

####################################################
##  Narendra Kumar, PhD                            #
##  Epigenetics Unit                               #
##  Institute of Cancer Sciences                   #
##  University of Glasgow, UK                      #
##  narekum@gmail.com,narendra.kumar@glasgow.ac.uk #
####################################################

BEGIN {
        use Cwd qw(realpath cwd);
        use File::Basename;
        our ($fn, $dir) = fileparse(realpath($0));
}
$dir=~s/\/$//;

#use strict ;
use Getopt::Long;
use File::Copy 'move';
use Parallel::ForkManager;

print '
   _____  _   _
  |  __ \| \ | |   /\
  | |__) |  \| |  /  \ ______ ___  ___  __ _
  |  _  /| . ` | / /\ \______/ __|/ _ \/ _` |
  | | \ \| |\  |/ ____ \     \__ \  __/ (_| |
  |_|  \_\_| \_/_/    \_\    |___/\___|\__, |
  Tophat2, fastqc, picard, trim galore    | |
  Samtools, bedtools, ucsc tools          |_|
   ____  _            _ _
  |  _ \(_)_ __   ___| (_)_ __   ___
  | |_) | | ._ \ / _ \ | | ._ \ / _ \
  |  __/| | |_) |  __/ | | | | |  __/
  |_|   |_| .__/ \___|_|_|_| |_|\___|
          |_|

';

my $cmd=$0." ".join(" ",@ARGV); ### command line copy

my $mapdir;
my $species;
my $quals;
my $threads=2;
my $dontzip="y";
my $extra;
my $gtf;
#my $insertsize=150;

my $start_time ;
my $help;
my $time_tag=$start_time=time;

GetOptions ('h|help'=>\$help,                     # --help         : print this help
            "d|dir=s" => \$mapdir,                # --dir          : directory where sequences are kept  (required)
            "g|gtf=s" => \$gtf,                   # --gtf          : gtf file (required)
            "s|species=s" => \$species,           # --species      : species e.g. hg19  (required)
            "q|quals=s" => \$quals,               # --quals        : quals (required)
            "t|threads=s" => \$threads,           # --threads      : number of processors to be used (default 2)
            "e|extra=s" => \$extra,               # --extra        : extra bowtie options enclosed in quote ""
            "z|zip=s" => \$dontzip,               # --zip          : zip the fastq files after the run y=yes; n=no (default y)
           )
or &PRINTHELP($fn);


if(defined $help || !$mapdir || !$species || !$quals || !$gtf) {
        &PRINTHELP($fn);
        exit;
}

my @options=( "$cmd",
           "--dir            $mapdir",
           "--gtf            $gtf",
           "--species        $species",
           "--quals          $quals",
           "--threads        $threads",
           "--extra          \"$extra\"",
           "--zip            $dontzip",
);

print "** Running $fn with the following options\n\n" ;
print "      $_\n" foreach @options ;
print "\n";

###############################################################################
###Set the values of following variables to the their respective locations of 
#your system
my $home="/home/naren";
my $tophatpath="$home/local/builds/tophat-2.0.11.Linux_x86_64/tophat2";
my $fastqc="$home/local/builds/FastQC/fastqc";
my $bedGraphToBigWigPath="$home/local/ucscTools/bedGraphToBigWig";
my $markduplicatespath="$home/local/builds/picard-tools-1.98/MarkDuplicates.jar";
my $trimgalore="$home/local/bin/trim_galore";
my $genomeCoverageBed="$home/local/src/bedtools-2.17.0/bin/genomeCoverageBed";
my $bedToBigBed="$home/local/ucscTools/bedToBigBed";
my $contaminant_list="$home/local/builds/FastQC/Contaminants/contaminant_list.txt";
my $fetchChromSizes="fetchChromSizes"; 
my $wigToBigWig="$home/local/ucscTools/wigToBigWig";
my $compressor="pigz";
my $BOWTIE2_INDEXES="/mnt/WDMyBook2TB/library/indexes/genomes/bowtie2-indexes/$species";
my $BOWTIE_INDEXES="/mnt/WDMyBook2TB/library/indexes/genomes/bowtie-indexes/$species";
###############################################################################

my $indextouse=$BOWTIE2_INDEXES;
my $transcriptomesuffix="";
if ( $extratophatopts =~ /--bowtie/  ){
	my $transcriptomesuffix=".bowtie1";
	my $indextouse=$BOWTIE_INDEXES;
	
}

################################################################################
## COPY SEQUENCES INTO SEQUENCE FOLDER
## CREATE ONE IF NOT EXIST
if ( -d "$mapdir/sequences" ) {
        print  "$mapdir/sequences data folder already exists - Not moving to sequences folder\n";
} else {
        mkdir "$mapdir/sequences" ;
        opendir(my $dh, $mapdir) || die "Can't opendir $mapdir: $!";
        my @fastqfiles = grep { /fastq.gz$|fastq$/ && -f "$mapdir/$_" } readdir($dh);
        closedir $dh;
        foreach (@fastqfiles){
                print "moving $mapdir/$_ to $mapdir/sequences\n";
                move ("$mapdir/$_", "$mapdir/sequences/") or die "move failed: $!";
        }
}

################################################################################
## QUALITY CHECK BY FASTQC 
## if the quality folder already exists we can skip, otherwise call FastQC

if ( -d "$mapdir/quality" ) {
        print  "$mapdir/quality data folder already exists - Not checking quality\n";
} else {
        mkdir "$mapdir/quality" ;
        system ("$fastqc --nogroup --threads $threads --contaminants $contaminant_list -f fastq $mapdir/sequences/*.fastq $mapdir/sequences/*.fastq.gz");
        system ("mv $mapdir/sequences/*fastqc* $mapdir/quality/");
}

################################################################################
## UNZIPPING sequences if Zipped.
	my @readsequences=READDIR("$mapdir/sequences",'fastq.gz$');
	print "\nUnzipping fastq sequences \n" if grep {defined($_)} @readsequences;
	foreach (@readsequences){
		print "Unzipping $mapdir/sequences/$_ ...";
		system ("gunzip $mapdir/sequences/$_");
		print "done\n";
	}

################################################################################
## TRIM sequences

if ( -d "$mapdir/trimmed-seq" ) {
	print "Trimmed sequence exists - Not trimming\n";
} else {
	print "Trimming sequences....\n";
	mkdir "$mapdir/trimmed-seq";
	my $qual;
	#ADD definition for other quals later.
	$qual="phred64" if $quals eq "phred64-quals";
	$qual="phred33" if $quals eq "phred33-quals";
	$qual="phred64" if $quals eq "solexa1.3-quals";
	my @trim_commands=();
	my @readsequences=READDIR("$mapdir/sequences",'_1.fastq$');
	foreach (@readsequences){
		(my $seqfile=$_) =~ s/_1.fastq$// ;
		push @trim_commands, "$trimgalore --quality 20 --$qual --paired --suppress_warn --retain_unpaired $mapdir/sequences/${seqfile}_1.fastq $mapdir/sequences/${seqfile}_2.fastq --output_dir $mapdir/trimmed-seq\n";
	}

	&PARALLELRUN(\@trim_commands,$threads);
	
	if ( -d "$mapdir/trimmed-seq/unpaired" ) {
		print "Not making unpaired folder\n";
	} else {
		mkdir "$mapdir/trimmed-seq/unpaired";
	}

	system ("mv $mapdir/trimmed-seq/*unpaired*fq* $mapdir/trimmed-seq/unpaired/");
	system ("rename \'s/_val_[12].fq/.fastq/\' $mapdir/trimmed-seq/*fq*");
	system ("rename \'s/_trimmed.fq/.fastq/\' $mapdir/trimmed-seq/*fq*");

	if ( -d "$mapdir/logs" ) {
		print "Logs folder exists\n";
	} else {
		mkdir "$mapdir/logs";
	}

	system ("mv $mapdir/trimmed-seq/*trimming_report.txt $mapdir/logs/");
	
}

###############################################################################
# do the alignment

if ( -d "$mapdir/aligned"){
	print "Alignment already done - Not aligning\n";
} else {
	mkdir "$mapdir/aligned" ;
	my @fastseq=READDIR("$mapdir/trimmed-seq",'_1.fastq$');
	foreach (@fastseq){
		(my $seqfile=$_) =~ s/_1.fastq$// ;
		my $tophatlog="$mapdir/trimmed-seq/$seqfile.tophat.log";
		print "$_\n";
		system ("date > $tophatlog");
		system ("echo $tophatlog");
		system ("$tophatpath --version >> $tophatlog");

		my $qual;	
		if ($quals eq "phred64-quals") {$qual="solexa1.3-quals"} 
		elsif ($quals eq "phred33-quals") {$qual="solexa-quals"} 
		else { $qual=$quals}

		print "Using quality: $qual\n";

		if ( $gtf ne "none") {
			if ( -d "$gtf.transcriptome-index$transcriptomesuffix" ) {
				print "No need to make transcriptome\n";
			} else {
				print "Need to make transcriptome $gtf.transcriptome-index$transcriptomesuffix\n";
				mkdir "$gtf.transcriptome-index$transcriptomesuffix";
			}
		}
		my $tophatgtf="";
		if ( $gtf eq "none" ){
			$tophatgtf="";
		} else {
			$tophatgtf="-M -G $gtf --transcriptome-index=$gtf.transcriptome-index$transcriptomesuffix";
			
		}
		my $outfolder="$mapdir/aligned/$seqfile";
		my $tophatcmd="$tophatpath --b2-sensitive --no-coverage-search -z pigz $extratophatopts --$qual -p $threads $tophatgtf -o $outfolder $indextouse $mapdir/trimmed-seq/${seqfile}_1.fastq $mapdir/trimmed-seq/${seqfile}_2.fastq";
		system ("echo $tophatcmd >> $tophatlog");
		print "***\n$tophatcmd\n***\n";
                system ("$tophatcmd >> $tophatlog 2>&1");
	}

        if ( -d "$mapdir/logs" ) {
                print "Logs folder exists\n";
        } else {
                mkdir "$mapdir/logs";
        }
	
	system ("mv $mapdir/trimmed-seq/*log $mapdir/logs/");

}

###############################################################################
## Remove duplicates
if ( -e "$mapdir/.dupsremoved" ){
	print "Dups already removed\n";
} else {
	print "Making bam files with duplicates removed\n";
	my @rmdupcommands=();
	opendir(my $dh, "$mapdir/aligned") || die "Can't opendir $dir: $!";
	my @accepted_hits = grep { ! /\.\.|\./ } readdir($dh); closedir $dh;
	foreach (@accepted_hits){
		push @rmdupcommands, "java -Xmx4g -jar $markduplicatespath REMOVE_DUPLICATES=true INPUT=$mapdir/aligned/$_/accepted_hits.bam OUTPUT=$mapdir/aligned/$_/accepted_hits.rmdup.bam METRICS_FILE=$mapdir/aligned/$_/dup.txt\n";
	}

	&PARALLELRUN(\@rmdupcommands,$threads);		
	system ("touch $mapdir/.dupsremoved");

}

###############################################################################
## INDEXING Bam files

if ( -e "$mapdir/.bamindexed" ){
        print "BAMS already indexed\n";
} else {
	print "Indexing BAMs\n";
	my @indexbamcommands=();
	opendir(my $dh, "$mapdir/aligned") || die "Can't opendir $dir: $!";
	my @samples = grep { ! /\.\.|\./ } readdir($dh); closedir $dh;
	foreach (@samples){
		print "$_\n";
		push @indexbamcommands, "samtools index $mapdir/aligned/$_/accepted_hits.bam\n";
		push @indexbamcommands, "samtools index $mapdir/aligned/$_/accepted_hits.rmdup.bam\n";
	}

	&PARALLELRUN(\@indexbamcommands,$threads);
	system ("touch $mapdir/.bamindexed");

}

###############################################################################
## Checking if Chromosome size file exists

if ( ! -f "$dir/../chrmSizes/chrmSizes.$species" ){
	print "$dir/../chrmSizes/chrmSizes.$species not found\n--fetching now\n";
        system ("$fetchChromSizes $species > $dir/../chrmSizes/chrmSizes.$species");
}

###############################################################################
## Create bigWig for genome browser

if ( -e "$mapdir/.wigscreated" ) {
        print "Wigs already created\n";
} else {
        print "Creating wigs\n";
	my @makewigcommands=();
        opendir(my $dh, "$mapdir/aligned") || die "Can't opendir $dir: $!";
        my @samples = grep { ! /\.\.|\./ } readdir($dh); closedir $dh;
        foreach (@samples){
                print "$_\n";
                push @makewigcommands, "$dir/../bin/createwig.sh $mapdir/aligned/$_/accepted_hits.bam $genomeCoverageBed $dir/../chrmSizes/chrmSizes.$species $wigToBigWig $species";
                push @makewigcommands, "$dir/../bin/createwig.sh $mapdir/aligned/$_/accepted_hits.rmdup.bam $genomeCoverageBed $dir/../chrmSizes/chrmSizes.$species $wigToBigWig $species";
        }

#	print "$_\n" foreach @makewigcommands ;
	&PARALLELRUN(\@makewigcommands,$threads);
	system ("touch $mapdir/.wigscreated");
}

if ( -e "$mapdir/.junctionscreated" ){
	print "Junctions already created\n";
} else {
	print "Creating junctions\n";
	my @junctioncommands=();
	opendir(my $dh, "$mapdir/aligned") || die "Can't opendir $dir: $!";
	my @samples = grep { ! /\.\.|\./ } readdir($dh); closedir $dh;
	foreach (@samples){
		system ("cat $mapdir/aligned/$_/junctions.bed | awk 'FNR > 1' | sort -k1,1 -k2,2n > $mapdir/aligned/$_/junctions.bed.trim");
		push @junctioncommands, "$bedToBigBed $mapdir/aligned/$_/junctions.bed.trim $dir/../chrmSizes/chrmSizes.$species $mapdir/aligned/$_/junctions.bigBed";
	}

	&PARALLELRUN(\@junctioncommands,$threads);
	system ("touch $mapdir/.junctionscreated");
}

###############################################################################
## Zip the sequence files in the end to same space if $dontzip eq y

if  ($dontzip eq "y" ) {
	print "Compressing the sequence files\n";
	system ("$compressor -r $mapdir/sequences/*.fastq");
	print "Compressing the trimmed sequence files\n";
        system ("$compressor -r $mapdir/trimmed-seq/*.fastq");
	print "Compressing the unpaired files\n";
        system ("$compressor -r $mapdir/trimmed-seq/unpaired/*.fq");
}


sub READDIR {
	my ($dir,$pattern) = @_;
	opendir(my $dh, $dir) || die "Can't opendir $dir: $!";
        my @files = grep { /$pattern/ && -f "$dir/$_" } readdir($dh);
        closedir $dh;
	return (@files);

}

sub PARALLELRUN {
	my ($commands,$nproc)=@_;
        my $pm = new Parallel::ForkManager($nproc);
        $pm->run_on_finish(     # callback function to retrieve and assemble matrix from each chromosome matrices.
        sub{
                my($pid, $exit_code, $ident, $exit_signal, $core_dump, $dr) = @_;
                        print "** $ident finished, pid: $pid\n";
           }
        );

        $pm->run_on_start( sub {
                my ($pid,$ident)=@_;
                print "** $ident started, pid: $pid\n";
        });

        foreach my $keys ( @$commands ) {
                my $pid = $pm->start($keys) and next;
                system("$keys");
                $pm->finish(0);

        }

        $pm->wait_all_children;
}

sub PRINTHELP {
        my ($fn)=@_;

        print <<DES;
Usage: $fn < --dir MAPDIRECTORY --gft GTFFILE --species SPECIES --quals QUALS > [Options]

Options:
  -h, --help       Show this help message.
  -d, --dir        directory where sequences are kept  (required)
  -g, --gtf        gtf file (required)
  -s, --species    species e.g. hg19  (required)
  -q, --quals  	   quals (required)
                   available quals: phred33-quals
                                    phred64-quals
                                    solexa1.3-quals
  -t, --threads    number of processors to be used (default 2)
  -e, --extra      extra bowtie options enclosed in quote ""
  -z, --zip        zip the fastq files after the run y=yes; n=no (default y)

Examples:
\$ perl $fn --dir mapdirectory --gft gtffile.gtf --species hg19 --quals phred33-quals
\\or\\
\$ perl $fn  -d mapdirectory -g gtf gtffile.gtf -s hg19 -q phred33-quals -t 8 -z n

DES
exit;
}
