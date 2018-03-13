#!/usr/local/bin/perl -w

# stage_rankmethod.pl

# AUTHOR: Adam Reid
# Copyright (C) 2017, 2018 Genome Research Ltd.
# This program is distributed under the terms of the GNU General Public License

use strict;

# Query dataset parameters
my $min_genes = 100;

# Reference dataset parameters (These can be adjusted for each new reference dataset depending on how variable genes are between conditions)
my $min_sd = 3;
my $min_mean = 30;
my $max_mean = 70;
#my $min_sd = 2;
#my $min_mean = 20;
#my $max_mean = 80;

# Command line arguments
my($query, $ref, $min_val, $ref_data) = @ARGV;

if(!defined $query)
{
	print STDERR "USAGE: <normalised query data> <ref data> <query data expression cutoff e.g. 10 [default] for FPKM> <ref meta>\n";

	exit;
}

if(!defined $min_val)
{
	$min_val = 10;
}

# Get stage meta data if present
my %meta = ();

if(defined $ref_data)
{
	open(R, "<$ref_data") or die "$!";

	while(<R>)
	{
		chomp;

		my($id, $name) = split /\t/;

		$meta{$id} = $name;
	}
	close R;
}

# Read in query data
open(Q, "<$query") or die "$!";

my $junk;
my @head = ();
my %query = ();

while(<Q>)
{
        chomp;

        if(/^id/)
        {
                ($junk, @head) = split /\t/;
        }
        else
        {
                my($id, @vals) = split /\t/;


                for(my $i=0;$i<scalar @head;$i++)
                {
                        $query{$id}->{$head[$i]} = $vals[$i];
			#print "$id\t$head[$i]\t$vals[$i]\n";
                }
        }
}
close Q;


# Make profiles of genes for each cell including only those genes expressed above a certain level
my %profiles = ();


foreach my $id (sort keys %query)
{
        foreach my $head (@head)
        {
                if($query{$id}->{$head} >= $min_val)
                {
                        $profiles{$head}->{$id} = $query{$id}->{$head};

                        #$used_genes{$head}++;
                }
        }
}

# Read in reference data
open(R, "<$ref") or die "$!";

my %ref = ();

my $total_ref_genes = ();

while(<R>)
{
        chomp;

        if(/^id/)
        {
                ($junk, @head) = split /\t/;
        }
        else
        {
                my($id, @vals) = split /\t/;

		# Determine whether this gene falls within acceptable parameters for usefullness
		my $sd = sd(\@vals);
		my $mean = mean(\@vals);

		if($sd < $min_sd || $mean < $min_mean || $mean > $max_mean)
		{
			next;
		}
		else
		{
                	for(my $i=0;$i<scalar @head;$i++)
                	{
                	        $ref{$head[$i]}->{$id} = $vals[$i];
                	}

			$total_ref_genes++;
		}
        }
}
close R;

print STDERR "$total_ref_genes reference genes passed the filter\n";

my %data = ();
my %results = ();

my %used_genes = ();

# Compare each query ranking to each reference ranking
foreach my $q_cond (sort keys %profiles)
{

        foreach my $t_cond (sort keys %ref)
        {
		$used_genes{$q_cond} = 0;

                my @q_vals = ();
                my @t_vals = ();

                my $t_total = 0;

                foreach my $id (sort keys %{$profiles{$q_cond}})
                {
                        if(exists $ref{$t_cond}->{$id})
                        {
                                push @q_vals, $profiles{$q_cond}->{$id};
                                push @t_vals, $ref{$t_cond}->{$id};

                                $t_total += $ref{$t_cond}->{$id};

				$used_genes{$q_cond}++;
                        }
                }

                if($t_total == 0)
                {
                        print STDERR "No positive values in reference for $t_cond and $q_cond\n";
                        next;
                }

		if($used_genes{$q_cond} < $min_genes)
		{
			#print STDERR "Not enough genes for $q_cond\n";
		
			$data{$q_cond}->{$t_cond} = 0;

			$results{$q_cond}->{'stage'} = 'none';
			$results{$q_cond}->{'r'} = 0;

			next;
		}

		my $r = 0;

		eval
		{
                	$r = spearman(\@q_vals, \@t_vals);
		};
		if($@)
		{
			print "@q_vals\n@t_vals\n";
			exit;
		}

                #print "$q_cond\t@q_vals\n$t_cond\t@t_vals\n\n";
                #print "$q_cond\t$t_cond\t$r\n";

                $data{$q_cond}->{$t_cond} = $r;

                unless(exists $results{$q_cond}->{'stage'} && $results{$q_cond}->{'r'} > $r)
                {
                        $results{$q_cond}->{'stage'} = $t_cond;
                        $results{$q_cond}->{'r'} = $r;
                }
        }
}

print "Sample\tPrediction\tGenes_used\tSpearmans_r\tstdev\tRes\n";

foreach my $cell (sort keys %results)
{
        my @data = ();

	my @cor_vals = ();

        foreach my $t_cond (@head)
        {
                #print "$t_cond\t$data{$cell}->{$t_cond}\n";
                if(!exists $data{$cell}->{$t_cond})
                {
                        $data{$cell}->{$t_cond} = 0;
                }

		my $t_cond_meta = $t_cond;

		if(exists $meta{$t_cond})
		{
			$t_cond_meta = $t_cond.'('.$meta{$t_cond}.')';
		}

                push @data, "$t_cond_meta:$data{$cell}->{$t_cond}";

		push @cor_vals, $data{$cell}->{$t_cond};
        }

	my $mean_cor_vals = mean(\@cor_vals);

	my $cov = 0;
	my $sd = 0;

	if($mean_cor_vals > 0)
	{
		$sd = sd(\@cor_vals);
		$cov = $sd / $mean_cor_vals;
	}

        my $data = join " ", @data;

	# Convert prediction using reference metadata if present
	if(exists $meta{$results{$cell}->{'stage'}})
	{
		$results{$cell}->{'stage'} = $meta{$results{$cell}->{'stage'}};
	}

        print "$cell\t$results{$cell}->{'stage'}\t$used_genes{$cell}\t$results{$cell}->{'r'}\t$sd\t$data\n";
}

# Pearson correlation
sub pearson
{
	my($a, $b) = @_;

	my $a_mean = mean($a);
	my $b_mean = mean($b);

	my $AoB = 0;
	my $var_a = 0;
	my $var_b = 0;

	for(my $i=0;$i<scalar @$a;$i++)
	{
		my $a_norm = $a->[$i] - $a_mean;
                my $b_norm = $b->[$i] - $b_mean;

                $AoB += (($a_norm) * ($b_norm));

                $var_a += ($a_norm * $a_norm);
                $var_b += ($b_norm * $b_norm);
	}

	if(sqrt($var_a) * sqrt($var_b) == 0)
	{
		print STDERR "Zero value sqrt variance with $a and $b\n";

		return 0;

	}
	else
	{
		my $r = ($AoB / (sqrt($var_a) * sqrt($var_b)));

		return $r;
	}
}

# Spearman correlation
sub spearman
{
        my($rand_profa, $rand_profb) = @_;

        my @ranka = ();
        my %ranka = ();
        my @rankb = ();
        my %rankb = ();

        my @sorta = sort {$a<=>$b} @$rand_profa;
        my @sortb = sort {$a<=>$b} @$rand_profb;

        for(my $i=0;$i<scalar @sorta;$i++)
        {
                push @{$ranka{$sorta[$i]}}, $i+1;
                push @{$rankb{$sortb[$i]}}, $i+1;
        }

        foreach(@{$rand_profa})
        {
                if(scalar @{$ranka{$_}} == 1)
                {
                        push @ranka, $ranka{$_}->[0];
                }
                else
                {
                        push @ranka, mean($ranka{$_});
                }
        }
        foreach(@{$rand_profb})
        {
                if(scalar @{$rankb{$_}} == 1)
                {
                        push @rankb, $rankb{$_}->[0];
                }
                else
                {
                        push @rankb, mean($rankb{$_});
                }
        }

        my $a_mean = mean($rand_profa);
        my $b_mean = mean($rand_profb);

        my $Ex = 0;
        my $Ey = 0;
        my $Ey2 = 0;
        my $Ex2 = 0;
        my $Exy = 0;

        for(my $i=0;$i<scalar @ranka;$i++)
        {
                $Ex += $ranka[$i];
                $Ey += $rankb[$i];
                $Ex2 += $ranka[$i]*$ranka[$i];
                $Ey2 += $rankb[$i]*$rankb[$i];
                $Exy += $ranka[$i]*$rankb[$i];
        }

        my $n = scalar @ranka;

        my $r = (($n*$Exy) - ($Ex*$Ey)) / sqrt((($n*$Ex2) - ($Ex*$Ex)) * (($n*$Ey2) - ($Ey*$Ey)));

        return $r;
}

# Calulate mean
sub mean
{
	my($vals) = @_;

	my $sum = 0;
	my $n = scalar @$vals;

	foreach (@$vals)
	{
		$sum += $_;
	}

	return($sum / $n);
}

# Randomly reorder an array
sub fy1 
{
	my $deck = shift;  # $deck is a reference to an array

	my $i = @$deck;

	while (--$i) 
	{
		my $j = int rand ($i+1);

		@$deck[$i,$j] = @$deck[$j,$i];
	}
}

# Perform multiple hypothesis testing using R
sub multihypo_test
{
        my($pvalues) = @_;

        my $p_string = join ",", @$pvalues;

        open(R, ">r.in") or die "$!";
        print R "p.adjust(c($p_string), method=c(\"BH\"))\n";
        close R;

        chomp(my $result = `R --no-save < r.in`);

        my @p_adjust;

        while ($result =~ /\[\d+\] (.*)/g)
        {
                push @p_adjust, split /\s/, $1;
        }

        return \@p_adjust;
}

# Calculate standard deviation
sub sd
{
        my($vals) = @_;

        my $mean = mean($vals);

        my $sum_sq_dev = 0;

        foreach my $val (@$vals)
        {
                $sum_sq_dev += (($val - $mean)**2);
        }

        my $var = $sum_sq_dev / scalar @$vals;

        return sqrt($var);
}
