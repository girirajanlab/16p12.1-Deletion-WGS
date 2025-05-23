# script to caluculate the proportion of a region that is unique.

# alterations must have a header row

# columns: chr bps bpe name (additional columsn are optional)

# SegDups file must be non-redundent (mergeBED)



if (@ARGV != 3) {

    print STDOUT "Usage: Script [RegionList] [SDList] [OutputFile]\n";

    print STDOUT "Region file format: <Chr><BPS><BPE><Name>...\n";

    print STDOUT "SD file format: <Chr><S><E>\n";

    exit;

}



open(CLONEFILE,"<$ARGV[0]");

open(SDFILE,"<$ARGV[1]");

open(OUTFILE,">$ARGV[2]");



%sdTable = ();

%cloneTable = ();



print STDOUT "Processing Control CNV file...\n";



# put SD data in hash table

while ($line = <SDFILE>) {

    chomp $line;

    ($Cchr,$CS,$CE) = split(/\t/,$line);

    ${$sdTable{$Cchr}}{$line} = "$CS\t$CE\t$Cstate)";

}

close(SDFILE);



print STDOUT "Processing Regions Calls...\n";



# parse your list of calls and find the cnvs



$clone = <CLONEFILE>;

chomp $clone;

#print OUTFILE "$clone\tUniqueBP\n";



while ($clone = <CLONEFILE>) {

    chomp $clone;

    ($chr,$bps,$bpe,$sid,@rest) = split(/\t/,$clone);

     $cnvs = 0;

     foreach $entry (sort keys %{$sdTable{$chr}}) {

	($CS,$CE) = split(/\t/,${$sdTable{$chr}}{$entry});

	# dup region is inside SD

	if (($CS <= $bps) && ($CE >= $bpe)) {

		$cnvs = $cnvs + ($bpe - $bps);    

	}

	# dup region is larger than SD

	elsif (($CS > $bps) && ($CE < $bpe)) {

		$cnvs = $cnvs + ($CE - $CS);

	}

	# first portion of SD overlap region 

	elsif (($CS >= $bps) && ($CE >= $bpe) && ($CS <= $bpe)) {

	    	$cnvs = $cnvs + ($bpe - $CS); 

	}

	# second half of SD overlaps region

	elsif (($CS <= $bps) && ($CE >= $bps) && ($CE <= $bpe)) {

	    $cnvs = $cnvs + ($CE - $bps); 

	}

    }

    print STDOUT "Done: $sid $chr $bps $bpe\n";

    $cnvs = ($bpe - $bps) - $cnvs;

    $columns = join("\t",@rest);

    print OUTFILE "$chr\t$bps\t$bpe\t$sid\t$cnvs\t$columns\n";

}



print STDOUT "Done\n";



close(CLONEFILE);

close(OUTFILE);



