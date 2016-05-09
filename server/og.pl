#!/usr/bin/perl -w
use strict;

# written by Bas E. Dutilh 2011

my $pid = $$;
my $start_time = time ();

my $usage = "\nUsage:\t$0 protein_fasta_files_directory blast_all_vs_all_directory\n\nNote:\tyou only need to search A-vs-A, A-vs-B and B-vs-B; B-vs-A can be skipped as I will assume it is a mirror of A-vs-B\n\n";

if (scalar @ARGV != 2) {
	die $usage; }

&print_time ();
print STDERR "starting up: $0 @ARGV\n";

my $fasta_directory = $ARGV[0];
my $blast_directory = $ARGV[1];
my $output_inps = "$pid.inparalogous_groups.txt";
my $output_OGs = "$pid.OGs.txt";

if ((! (-e $blast_directory)) || (! (-e $fasta_directory))) {
	&print_time ();
	die "$blast_directory or $fasta_directory do not exist\n"; }
if ((-e $output_inps) || (-e $output_OGs)) {
	&print_time ();
	die "$output_inps or $output_OGs exists: NOT overwriting\n"; }
if (! (open (INP, ">$output_inps"))) {
	&print_time ();
	die "Can't write to $output_inps: $!\n"; }
if (! (open (OG, ">$output_OGs"))) {
	&print_time ();
	die "Can't write to $output_OGs: $!\n"; }

# determine protein lengths
my %lengths = ();
if (! (opendir (DIR, $fasta_directory))) {
	&print_time ();
	die "Can't open directory $fasta_directory: $!\n"; }
my @files = readdir (DIR);
closedir (DIR);
&print_time ();
print STDERR "reading protein lengths from $fasta_directory\n";
FILE: foreach my $file (@files) {
	if ($file !~ /\w/o) {
		&print_time ();
		print STDERR "skipping file $fasta_directory/$file\n";
		next FILE; }
	if (! (open (LEN, "<$fasta_directory/$file"))) {
		&print_time ();
		die "can't open $fasta_directory/$file: $!\n"; }
	my $g = $file;
	my $i = "";
	my $q = "";
	while (my $line = <LEN>) {
		chomp ($line);
		$line =~s/\r+//og;
		if ($line =~ s/^>//o) {
			if ($q) {
				$q =~ s/\W+//og;
				$lengths{$g}{$i} = length ($q);
				$q = ""; }
			$i = $line; }
		else {
			$q .= $line; } }
	$q =~ s/\W+//og;
	$lengths{$g}{$i} = length ($q); }
close (LEN);

# read all_vs_all similarities
if (! (opendir (DIR, $blast_directory))) {
	&print_time ();
	die "Can't open directory $blast_directory: $!\n"; }
@files = readdir (DIR);
closedir (DIR);
my %sim = ();
my %proteins = ();
# 0      1      2      3      4        5       6      7    8      9    10     11
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
FILE: foreach my $file (sort @files) {
	if ($file !~ /\w/o) {
		next FILE; }
	my $tmp = $file;
	my $genome_from = "";
	FROM: foreach my $this_genome (sort { length($b) <=> length ($a); } keys %lengths) {
		if ($tmp =~ s/$this_genome//) {
			$genome_from = $this_genome;
			last FROM; } }
	my $genome_to = "";
	TO: foreach my $this_genome (sort { length($b) <=> length ($a); } keys %lengths) {
		if ($tmp =~ s/$this_genome//) {
			$genome_to = $this_genome;
			last TO; } }
	if (($genome_from eq "") || ($genome_to eq "")) {
		&print_time ();
		print STDERR "skipping file $blast_directory/$file\n";
		next FILE; }
	&print_time ();
	print STDERR "reading file $blast_directory/$file\n";
	my %top_scores = ();
	if (! (exists ($lengths{$genome_from}))) {
		&print_time ();
		die "no length info for $genome_from\n"; }
	if (! (exists ($lengths{$genome_to}))) {
		&print_time ();
		die "no length info for $genome_to\n"; }
	if (! (open (IN, "<$blast_directory/$file"))) {
		die "Can't open $blast_directory/$file: $!\n"; }
	LINE: while (my $line = <IN>) {
		chomp ($line);
		$line =~s/\r+//og;
		my @a = split /\t/o, $line;
		my $query = "";
		if (exists ($lengths{$genome_from}{$a[0]})) {
			 $query = $a[0]; }
		elsif (exists ($lengths{$genome_from}{$a[1]})) {
			$query = $a[1]; }
		my $hit = "";
		if (exists ($lengths{$genome_to}{$a[0]})) {
			$hit = $a[0]; }
		elsif (exists ($lengths{$genome_to}{$a[1]})) {
			$hit = $a[1]; }
		if (($query eq "") || ($hit eq "")) {
			&print_time ();
			die "no length info for $a[0] or $a[1]\n"; }
		if (($a[3] < $lengths{$genome_from}{$query} * .5) || ($a[3] < $lengths{$genome_to}{$hit} * .5)) {
			next LINE; }
		if (! (exists ($top_scores{$query}))) {
			$top_scores{$query} = $a[11]; }
		elsif ($top_scores{$query} < $a[11]) {
			$top_scores{$query} = $a[11];
			if ($genome_from ne $genome_to) {
				delete ($sim{$genome_from}{$query}{$genome_to}); } }
		if (! (exists ($top_scores{$hit}))) {
			$top_scores{$hit} = $a[11]; }
		elsif ($top_scores{$hit} < $a[11]) {
			$top_scores{$hit} = $a[11];
			if ($genome_from ne $genome_to) {
				delete ($sim{$genome_to}{$hit}{$genome_from}); } }
		++$proteins{$genome_from}{$query};
		++$proteins{$genome_to}{$hit};
		if ((($top_scores{$query} == $a[11]) || ($genome_from eq $genome_to))
			&& ((! (exists ($sim{$genome_from}{$query}{$genome_to}{$hit}))) || ($sim{$genome_from}{$query}{$genome_to}{$hit} < $a[11]))) {
			$sim{$genome_from}{$query}{$genome_to}{$hit} = $a[11]; }
		if ((($top_scores{$hit} == $a[11]) || ($genome_from eq $genome_to))
			&& ((! (exists ($sim{$genome_to}{$hit}{$genome_from}{$query}))) || ($sim{$genome_to}{$hit}{$genome_from}{$query} < $a[11]))) {
			$sim{$genome_to}{$hit}{$genome_from}{$query} = $a[11]; } }
	close (IN); }

if (scalar keys %sim == 0) {
	&print_time ();
	die "no similarity scores in input files\n"; }

# find best hits outside genome for each protein
&print_time ();
print STDERR "finding best one-directional hits for each protein\n";
my %best_hits = ();
foreach my $genome_from (keys %sim) {
	foreach my $protein_from (keys %{$sim{$genome_from}}) {
		TO: foreach my $genome_to (keys %{$sim{$genome_from}{$protein_from}}) {
			# outside own genome only, otherwise you will find the same gene
			if ($genome_to eq $genome_from) {
				next TO; }
			$best_hits{$genome_from}{$protein_from}{$genome_to}{"score"} = 0;
			foreach my $protein_to (keys %{$sim{$genome_from}{$protein_from}{$genome_to}}) {
				if ($best_hits{$genome_from}{$protein_from}{$genome_to}{"score"} < $sim{$genome_from}{$protein_from}{$genome_to}{$protein_to}) {
					delete ($best_hits{$genome_from}{$protein_from}{$genome_to});
					$best_hits{$genome_from}{$protein_from}{$genome_to}{"score"} = $sim{$genome_from}{$protein_from}{$genome_to}{$protein_to}; }
				# record all equally high soring top hits
				if ($best_hits{$genome_from}{$protein_from}{$genome_to}{"score"} == $sim{$genome_from}{$protein_from}{$genome_to}{$protein_to}) {
					++$best_hits{$genome_from}{$protein_from}{$genome_to}{"hits"}{$protein_to}; } } } } }

# compose in-paralogous groups for all genomes
my %inparalogs = ();
my %prot2inp = ();
foreach my $genome (sort keys %sim) {
	&print_time ();
	print STDERR "making greedy in-paralogous groups for $genome\n";
	my %this_inpar = ();
	foreach my $protein_from (keys %{$sim{$genome}}) {
		# find in-paralogs
		foreach my $protein_to (keys %{$sim{$genome}{$protein_from}{$genome}}) {
			if ($sim{$genome}{$protein_from}{$genome}{$protein_to} >= $sim{$genome}{$protein_from}{$genome}{$protein_from}) {
				++$this_inpar{$genome}{$protein_from}{$protein_to};
				++$this_inpar{$genome}{$protein_to}{$protein_from}; } } } # reciprocal = greedy!
	# make in-paralogous groups
	%{$inparalogs{$genome}} = ();
	my $i = 0;
	P: foreach my $protein_from (keys %{$sim{$genome}}) {
		foreach my $inp_g (keys %{$inparalogs{$genome}}) {
			if (exists ($inparalogs{$genome}{$inp_g}{$protein_from})) {
				foreach my $protein_to (keys %{$this_inpar{$genome}{$protein_from}}) {
					++$inparalogs{$genome}{$inp_g}{$protein_to}; }
				next P; } }
		++$inparalogs{$genome}{"INP_G.$i.$genome"}{$protein_from};
		foreach my $protein_to (keys %{$this_inpar{$genome}{$protein_from}}) {
			++$inparalogs{$genome}{"INP_G.$i.$genome"}{$protein_to}; }
		++$i; }
	foreach my $inp_g (keys %{$inparalogs{$genome}}) {
		foreach my $protein_from (keys %{$inparalogs{$genome}{$inp_g}}) {
			if (exists ($prot2inp{$genome}{$protein_from})) {
				foreach my $protein_to (keys %{$inparalogs{$genome}{$prot2inp{$genome}{$protein_from}}}) {
					$prot2inp{$genome}{$protein_to} = $inp_g;
					++$inparalogs{$genome}{$inp_g}{$protein_to}; }
				delete ($inparalogs{$genome}{$prot2inp{$genome}{$protein_from}}); }
			$prot2inp{$genome}{$protein_from} = $inp_g; } }
	foreach my $prot (keys %{$proteins{$genome}}) {
		if (! (exists ($prot2inp{$genome}{$prot}))) {
			++$inparalogs{$genome}{"INP_G.$i.$genome"}{$prot};
			$prot2inp{$genome}{$prot} = "INP_G.$i.$genome";
			++$i; } } }
%sim = ();			

# output inparalogous groups
&print_time ();
print STDERR "output inparalogous groups\n";
foreach my $genome (sort keys %inparalogs) {
	foreach my $inp (keys %{$inparalogs{$genome}}) {
		foreach my $protein (keys %{$inparalogs{$genome}{$inp}}) {
			print INP "$genome\t$inp\t$protein\n"; } } }
#			print INP "$genome\t$inp\t$genome.peg.$protein\n"; } } }
close (INP);

# find BBH of inparalogous groups between genomes
my %bbh = ();
foreach my $genome_from (sort keys %inparalogs) {
	&print_time ();
	print STDERR "finding BBHs for $genome_from ...\n";
	foreach my $inp_from (keys %{$inparalogs{$genome_from}}) {
		my %this_hit_inps = ();
		foreach my $protein_from (keys %{$inparalogs{$genome_from}{$inp_from}}) {
			foreach my $genome_to (keys %{$best_hits{$genome_from}{$protein_from}}) {
				foreach my $protein_to (keys %{$best_hits{$genome_from}{$protein_from}{$genome_to}{"hits"}}) {
					if (! (exists ($prot2inp{$genome_to}{$protein_to}))) {
						&print_time ();
						print STDERR "no prot2inp{$genome_to}{$protein_to}\n"; }
					++$this_hit_inps{$genome_to}{$prot2inp{$genome_to}{$protein_to}}; } } }
		foreach my $genome_to (keys %this_hit_inps) {
			foreach my $inp_to (keys %{$this_hit_inps{$genome_to}}) {
				foreach my $protein_to (keys %{$inparalogs{$genome_to}{$inp_to}}) {
					if (exists ($best_hits{$genome_to}{$protein_to}{$genome_from})) {
						foreach my $protein_from (keys %{$best_hits{$genome_to}{$protein_to}{$genome_from}{"hits"}}) {
							if (! (exists ($prot2inp{$genome_from}{$protein_from}))) {
								&print_time ();
								die "No prot2inp{$genome_from}{$protein_from}\n"; }
							if ($prot2inp{$genome_from}{$protein_from} eq $inp_from) {
								++$bbh{$genome_from}{$inp_from}{$genome_to}{$inp_to}; } } } } } } } }
%prot2inp = ();

# define OGs
my %inp2OGs = ();
my %OGs2inp = ();
my $this_OG = 0;
foreach my $genome_from (sort keys %bbh) {
	&print_time ();
	print STDERR "defining OGs for $genome_from\n";
	foreach my $inp_from (keys %{$bbh{$genome_from}}) {
		if (exists ($inp2OGs{$inp_from})) {
			my $old_OG = $inp2OGs{$inp_from};
			foreach my $inp_to (keys %{$OGs2inp{$old_OG}}) {
				$inp2OGs{$inp_to} = $this_OG;
				++$OGs2inp{$this_OG}{$inp_to}; }
			delete ($OGs2inp{$old_OG}); }
		$inp2OGs{$inp_from} = $this_OG;
		++$OGs2inp{$this_OG}{$inp_from};
		foreach my $genome_to (keys %{$bbh{$genome_from}{$inp_from}}) {
			foreach my $inp_to (keys %{$bbh{$genome_from}{$inp_from}{$genome_to}}) {
				if (exists ($inp2OGs{$inp_to})) {
					if ($inp2OGs{$inp_to} ne $this_OG) {
						my $old_OG = $inp2OGs{$inp_to};
						foreach my $inp_to_to (keys %{$OGs2inp{$old_OG}}) {
							$inp2OGs{$inp_to_to} = $this_OG;
							++$OGs2inp{$this_OG}{$inp_to_to}; }
						delete ($OGs2inp{$old_OG});
						$inp2OGs{$inp_to} = $this_OG;
						++$OGs2inp{$this_OG}{$inp_to}; } } } }
		++$this_OG; } }
%bbh = ();

# output OGs
&print_time ();
print STDERR "output OGs\n";
$this_OG = 0;
foreach my $og (keys %OGs2inp) {
	++$this_OG;
	INP: foreach my $inp (keys %{$OGs2inp{$og}}) {
		my $genome = $inp;
		$genome =~ s/^INP_G\.\d+\.//o;
		if (! (exists ($inparalogs{$genome}{$inp}))) {
			&print_time ();
			print STDERR "no inparalogs{$genome}{$inp}\n";
			next INP; }
		foreach my $protein (keys %{$inparalogs{$genome}{$inp}}) {
			print OG "$protein\t$inp\tVibrioOG$this_OG\n"; } } }
#			print OG "$genome.peg.$protein\t$inp\tVibrioOG$this_OG\n"; } } }
close (OG);

sub print_time () {
        my $time = time () - $start_time;
        my $min = 0;
        while ($time - 60 >= 0) {
                ++$min;
                $time -= 60; }
        if ($time < 10) {
                $time = "0$time"; }
        print STDERR "$min:$time\t$pid\t"; }

