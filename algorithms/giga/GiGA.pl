#!/usr/bin/perl -w
#use strict;

printhelp() unless @ARGV;

my $starttime = time();

my ($annfile, $namefile, $genenamefile, $sensitivity, $infile, $outfile, $random, $networkfile, $network_description, $max_number);

my $upregulated = 1;
my $textoutput = 0;

foreach (@ARGV) {
    $genenamefile = substr($_,2) if (/^-G/i);
    $sensitivity = substr($_,2) if (/^-T/i);
    $infile = substr($_,2) if (/^-I/i);
    $networkfile = substr($_,2) if (/^-N/i);
    $network_description = substr($_,2) if (/^-X/i);	
    $outfile = substr($_,2) if (/^-O/i);
    $max_number = substr($_,2) if (/^-M/i);	
	$upregulated = 0 if (/^-D/i);
	$textoutput = 1 if (/^-Ftxt/i);
    $random = 1 if (/^-R/i);
    printhelp() if (/^-H/i);
}

$sensitivity = 1 unless ($sensitivity);
$max_number = 20 unless ($max_number);

printhelp("N") unless (($networkfile) or ($network_description));

my (%neighbors, %node_ranks, %node_real_ranks, %node_values, %edge_labels);

if ($networkfile) {
    open NETWORKFILE, "<$networkfile" or die "Can't open $networkfile: $!";
    while (<NETWORKFILE>) {
	next if /^#/;
	chomp;
	my ($a, $b, $label) = split /\t/;    
        
	$a =~ s/^"*//;
        $a =~ s/"*$//;
	$b =~ s/^"*//;
        $b =~ s/"*$//;
  
	if ($a gt $b) {($b, $a) = ($a, $b)};
	push @{ $edge_labels{$a."|".$b} }, $label;
	push @{ $neighbors{$a} }, $b;
	push @{ $neighbors{$b} }, $a;
    }
    close NETWORKFILE;
    foreach (keys %neighbors) {
	my %n;
	@n{@{ $neighbors{$_} }} = ();
	@{$neighbors{$_}} = keys %n;
    }

} else {
	open NETWORK, "<$network_description" or die "can't open $network_description: $!";
	my %labels;	
	while (<NETWORK>) {
		next if /^#/;
		chomp;
		my ($gene, @labels) = split /\t|\|/;
		$gene =~ s/^"*//;
		$gene =~ s/"*$//;
		foreach my $l (@labels) {
			push @{$labels{$l}}, $gene;
		}		
	}
	
	foreach my $l (keys %labels) {
		if (@{$labels{$l}} > 150) {print STDERR "Unspecific annotation $l ignored.\n"; next};
		my @sort_genes = sort {$b cmp $a} @{$labels{$l}};
		foreach my $n1 (0..$#sort_genes-1) {
			foreach my $n2 ($n1+1..$#sort_genes) {
				my $a = $sort_genes[$n1];
				my $b = $sort_genes[$n2];
				if ($a gt $b) {($b, $a) = ($a, $b)};
				push @{ $edge_labels{$a."|".$b} }, $l;
				push @{ $neighbors{$a} }, $b;
				push @{ $neighbors{$b} }, $a;
			}	 
		}
	}
	close NETWORK;
	foreach (keys %neighbors) {
		my %n;
		@n{@{ $neighbors{$_} }} = ();
		@{$neighbors{$_}} = keys %n;
	}
    
	
}


my %genenames;

if ($genenamefile) {
    open GENENAMEFILE, "<$genenamefile" or die "Can't open $genenamefile: $!";
    while (<GENENAMEFILE>) {
	next if /^#/;
	chomp;
	my ($g, $n) = split /\t/;
	($g) = ($g =~ /"*\s*(.*[^"]+)/);
        ($n) = ($n =~ /"*\s*(.*[^"]+)/);
	$genenames{$g} = $n;
        $genenames{$g} = "" unless ($genenames{$g});
    }
    close GENENAMEFILE;
} else {$genenamefile = ""}


if ($infile) {
    open STDIN, "<$infile" or die "Can't open $infile: $!";
} else {$infile = "STDIN"}

my @datafile = ();

while (<STDIN>) {
    push @datafile, $_ unless (/^#/);
} 

if ($random) {
    my @random;
    for my $i (0..$#datafile) {
	$random[$i] = rand();
    }
    my @randomindex = sort {$random[$a] <=> $random[$b]} (0..$#datafile);
    @datafile[@randomindex] = @datafile;
}

if ($outfile) {
    open OUTFILE, ">$outfile" or die "Can't open $outfile: $!";
} else {
	$outfile = "STDOUT"; 
    open OUTFILE, ">-" or die "Can't open $outfile: $!";

}


my $r = 0;
my $real_r = 0; 

foreach (@datafile) {
    next if (/^#/ or /^"#/ or /^\W*$/);
   	     
    
    chomp;
    my ($g, $exp) = split /\t/;
    
    unless ($exp) {$exp = 0};
    
    $g =~ s/^"*//;
    $g =~ s/"*$//;

    $real_r++;	
    $node_real_ranks{$g} = $real_r;	

    next unless ($neighbors{$g});

    $r++;

    $node_ranks{$g} = $r;
    $node_values{$g} = $exp;
	     
}

my $genenumber = $r;

#prune network to remove nodes not represented in experiment

my %old_neighbors;
%old_neighbors = %neighbors;
foreach my $g (keys %neighbors) {
    my %n;
    @n{@{$neighbors{$g}}} = ();
    foreach my $n (keys %n) {
        delete $n{$n} unless ($node_ranks{$n});
        $genenames{$n} = "" unless ($genenames{$n});
    }
    @{$neighbors{$g}} = keys %n;
	delete $neighbors{$g} unless (@{$neighbors{$g}});
}

# identify local minima;

my @loc_min;

foreach my $n1 (keys %node_ranks) {
    my $min = 1;
    foreach my $n2 (@{$neighbors{$n1}}) {
        next unless $node_ranks{$n2};
        if ($node_ranks{$n1} > $node_ranks{$n2}) {$min = 0; last} 
    }
    push @loc_min, $n1 if ($min);
}

# delimit regulated neighborhood around each local minimum

my %p_min;
my $threshold = $sensitivity/$genenumber;
my %min_groups;
my %max_ranks;
my %current;


foreach my $l (@loc_min) {

    my @min = $l;
    $max_ranks{$l} = $node_ranks{$l};
	@{$min_groups{$l}} = @min;

    $p_min{$l} = p_value (1, $genenumber, $max_ranks{$l} );

    while (1==1) {
         my ($new_max_rank, @new_min);
         ($new_max_rank, @new_min) = expand ($l, @min);
    
         last if ($new_max_rank == -1);
         last if (@new_min >= $max_number);

         @min = @new_min;
         my $new_p = p_value ( scalar @new_min, $genenumber, $new_max_rank);
         if ($new_p < $p_min{$l}) {
	     $p_min{$l} = $new_p;
	     @{$min_groups{$l}} = @min;
	     $max_ranks{$l} = $new_max_rank;
	 }


    }

}

# output

my %already_done;

if ($textoutput) {

# Output in plain text format

foreach my $l (sort {$p_min{$a} <=> $p_min{$b}} @loc_min) {

    if ($p_min{$l} < $threshold) {
        unless (exists $already_done{$l}) {
           print OUTFILE "$l\t$genenames{$l}\t$p_min{$l}\t$max_ranks{$l}\n";

	    my $cc = 0;
            foreach (sort {$node_ranks{$a} <=> $node_ranks{$b}} @{$min_groups{$l}}) {
		$cc++;
                print OUTFILE "-$cc-\t$_\t$genenames{$_}\t$node_ranks{$_}\t$node_real_ranks{$_}\n";
            }
	    @already_done{@{$min_groups{$l}}} = ();
	}
    } 
}

print OUTFILE "total genes measured in network: $genenumber.\n";

} else {

# Graph Description Language Output

my $first = 1;

print OUTFILE "graph: {\nlayout_algorithm: forcedir\n\ninfoname1: \"locus name\" infoname2: \"rank\" infoname3: \"expression value\"\n";
print OUTFILE "edge.arrowstyle: none\n"; # original contained: display_edge_labels: yes\n
print OUTFILE "node.fontname: \"helvR12\"\n";

print OUTFILE 'node: {title: "VIRTUAL CENTER" label: "" scaling: 0.0 focus}'; print "\n"; 
print OUTFILE "colorentry 42: 230 230 230\n";
print OUTFILE "colorentry 43: 240 240 210\n";
print OUTFILE "colorentry 44: 255 0 0\n";
print OUTFILE "colorentry 45: 255 70 70\n";
print OUTFILE "colorentry 46: 255 140 140\n";
print OUTFILE "colorentry 47: 255 210 210\n";
print OUTFILE "colorentry 48: 0 255 0\n";
print OUTFILE "colorentry 49: 70 255 70\n";
print OUTFILE "colorentry 50: 140 255 140\n";
print OUTFILE "colorentry 51: 210 255 210\n";

my @colors = (44,45,46,47,48,49,50,51);

my $back_color;
my $color;

my $class_counter = 0;

foreach my $l (sort {$p_min{$a} <=> $p_min{$b}} @loc_min) {

    if ($p_min{$l} < $threshold) {
        unless (exists $already_done{$l}) {
           print STDERR "$l\t$genenames{$l}\t$p_min{$l}\t$max_ranks{$l}\n";
			
		print OUTFILE "edge: {source: \"VIRTUAL CENTER\" target:\"SUBGRAPH $l\" label: \"\" linestyle: invisible}\n";
		$class_counter++;
		$back_color = ($back_color == 42) ? 43 : 42;
		my @title_parts = split /\s/, $genenames{$l};
		if (@title_parts > 5) {@title_parts = @title_parts[0..4]; $title_parts[4] .= "..."; }
		my $top_title = join '\n', @title_parts;
		my $p_string = sprintf "%3.2e", $p_min{$l};
		print OUTFILE "graph: {title: \"SUBGRAPH $l\" status: folded color: $back_color label: \"";
		print OUTFILE scalar(@{$min_groups{$l}}), " genes\\n$top_title et al.\\np=$p_string\"\n";

		if ($first) {$first = 0; print "borderstyle: triple bordercolor: red\n"};


		my %all_labels;
		my %label_connect;
		$first = 0;

		my %loc;
		@loc{@{$min_groups{$l}}} = ();

		foreach my $g (sort keys %loc) {
			my @title_parts = split /\s/, $genenames{$g};
			if (@title_parts > 5) {@title_parts = @title_parts[0..4]; $title_parts[4] .= "..."; }
			if ($upregulated) {
				if ($node_ranks{$g} < $genenumber/5) {$color = $colors[3];}	
				if ($node_ranks{$g} < $genenumber/10) {$color = $colors[2];}
				if ($node_ranks{$g} < $genenumber/20) {$color = $colors[1];}
				if ($node_ranks{$g} < $genenumber/100) {$color = $colors[0];}
			
				if ($node_ranks{$g} > $genenumber - $genenumber/5) {$color = $colors[7];}
				if ($node_ranks{$g} > $genenumber - $genenumber/10) {$color = $colors[6];}
				if ($node_ranks{$g} > $genenumber - $genenumber/20) {$color = $colors[5];}
				if ($node_ranks{$g} > $genenumber - $genenumber/100) {$color = $colors[4];}
			} else {
				if ($node_ranks{$g} < $genenumber/5) {$color = $colors[7];}	
				if ($node_ranks{$g} < $genenumber/10) {$color = $colors[6];}
				if ($node_ranks{$g} < $genenumber/20) {$color = $colors[5];}
				if ($node_ranks{$g} < $genenumber/100) {$color = $colors[4];}
			
				if ($node_ranks{$g} > $genenumber - $genenumber/5) {$color = $colors[3];}
				if ($node_ranks{$g} > $genenumber - $genenumber/10) {$color = $colors[2];}
				if ($node_ranks{$g} > $genenumber - $genenumber/20) {$color = $colors[1];}
				if ($node_ranks{$g} > $genenumber - $genenumber/100) {$color = $colors[0];}
			}

			my $node_title = join '\n', @title_parts;
			print OUTFILE "node: {title: \"$g\" label: \"$node_title\" color: $color info1: \"$g\"";
			print OUTFILE "info2: \"$node_ranks{$g}\" info3: \"$node_values{$g}\"}\n";
			NEIGHBOR: foreach my $n (sort @{$neighbors{$g}}) {
			
				next NEIGHBOR unless (exists $loc{$n});
				my ($a, $b);
				if ($g gt $n) {($b, $a) = ($g, $n)} else {($a, $b) = ($g, $n)};
				
				next NEIGHBOR unless (exists $edge_labels{$a."|".$b} );
				my @edges = @{$edge_labels{$a."|".$b}};

				LABEL: foreach my $label (@edges) {
					next LABEL if (exists $label_connect{$n."|".$label}); 
					$all_labels{$label} = ();
					$label_connect{$n."|".$label} = ();
					print OUTFILE "edge: {source:\"$n\" target: \"$label\"}\n";
				}

			}
				

		}
		foreach my $label (keys %all_labels) {
			print OUTFILE "node: {title: \"$label\" scaling: 0.5}\n";
		}
		
		print OUTFILE "}\n";

	    my $cc = 0;
            foreach (sort {$node_ranks{$a} <=> $node_ranks{$b}} @{$min_groups{$l}}) {
		$cc++;
                print STDERR "-$cc-\t$_\t$genenames{$_}\t$node_ranks{$_}\n";
            }
	    @already_done{@{$min_groups{$l}}} = ();
	}
    } 
}


print OUTFILE "}\n";

print STDERR "total genes measured in network: $genenumber.\n";
}


my $elapsed = time()-$starttime;
print STDERR "$elapsed sec. elapsed.\n";

exit;

sub expand {
	my ($l, @array) = @_;
    my @m = sort {$node_ranks{$b} <=> $node_ranks{$a}} @array;
    my $max = $node_ranks{$m[0]};

    %current = ();
    @current{@m} = ();

    my %neigh;
    foreach my $g (keys %current) {
	@neigh{@{$neighbors{$g}}} = ();
    }

    foreach my $g (keys %current) {
	delete $neigh{$g};
    }
    

    my @neigh = sort {$node_ranks{$a} <=> $node_ranks{$b}} keys %neigh;

    return (-1, undef) unless (@neigh); #no expansion possible

    my $new_max = $node_ranks{$neigh[0]};

    foreach my $c (0..$#neigh) {
	last if ($node_ranks{$neigh[$c]} > $new_max);
	$current{$neigh[$c]} =();
    }

    while (1) {
	my %new_neigh;
	foreach my $g (keys %current) {
	    @new_neigh{@{$neighbors{$g}}} = ();
	}

	foreach my $g (keys %current) {
	    delete $new_neigh{$g};
	}
    

	my @new_neigh = sort {$node_ranks{$a} <=> $node_ranks{$b}} keys %new_neigh;

	last unless (@new_neigh); #no further expansion possible

	last if ($node_ranks{$new_neigh[0]} > $new_max);

	foreach my $c (0..$#new_neigh) {
	    last if ($node_ranks{$new_neigh[$c]} > $new_max);
	    $current{$new_neigh[$c]} = ();
	
	}	
    }
    

    my @new_group = keys %current;

	my @sort_gr = sort {$node_ranks{$b} <=> $node_ranks{$a}} @new_group;
	print STDERR "$new_max - $node_ranks{$sort_gr[0]}\n" if ($new_max != $node_ranks{$sort_gr[0]});


    return ($new_max, @new_group);

}


sub p_value {
    my ($group_genes, $total_genes, $group_max_rank) = @_;
    
    return ($group_max_rank/$total_genes) if ($group_genes == 1);

    my $p = $group_max_rank/$total_genes;
    for my $c (1..($group_genes-1)) {
	$p *= ($group_max_rank-$c)/($total_genes-$c);
    }
    return $p;

}


sub follow {
    my ($g, $max) = @_;
    return if ($node_ranks{$g} > $max);

    $current{$g} = ();

    foreach my $n (@{$neighbors{$g}}) {
	next if ($node_ranks{$n} > $max);
	next if (exists $current{$n});	

	follow($n, $max);
    }

}


sub printhelp {

    print STDERR "Network specification file (option -n or -x) is required!\n\n" if (defined $_[0] && $_[0] eq "N");

    print STDERR <<ENDHELP;
Parameters of $0:
-i[INPUTFILE]: list of genes sorted by some parameter of differential expression.
               The first column must contain the gene identifiers as used in the
               annotation file. If no input file is specified, reads from STDIN.
               The (optional) second column can contain an expression value, e.g.
               the log fold-change.

-o[OUTPUTFILE]: Output file that will contain the results of GroupAnalysis. If no 
                output file is specified writes to STDOUT.

-n[NETWORKFILE]: Three tab-delimited columns, describing the network structure
                 to be analysed. The first two columns specify two nodes that
                 are connected by an edge, the third column the label of this 
                 edge.

-x[ALTERNATIVE NETWORKFILE]: Tab-delimited text file, the first column containing
		             a gene identifier, the following columns the annotations 
			     (labels) associated with this gene. This file will only
                             be used if no network file (-n) is specified. It provides
                             a more economical description of the graph. 

-g[GENENAMESFILE]: Gene names associated with gene identifiers. Column 1 contains
                   the gene identifiers as used in input and annotation files, 
                   column 2 contains the gene name or description.

-t[Number]: Threshold sensitivity of the analysis. Smaller numbers yield a more restricted
            list of groups in the result. Defaults to 1.

-m[Number]: Specifies the maximal number of genes in each cluster. Default = 20.

-f[TXT|GDL]: Output format plain text (TXT) or graph description language (GDL). Default is GDL.

-r : Randomize the order of genes in the input file.

-h : prints this help and exits.


ENDHELP

    exit;
}
