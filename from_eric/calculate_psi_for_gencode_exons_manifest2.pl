use warnings;
use strict;

my %triplet2exons;
my %exons_for_position;
my %all_features;
my %enst2type;
my %enst2ensg;
my %ensg2name;
my %ensg2type;

my %a5ss2exons;
my %a3ss2exons;
my %a3ss_hash;
my %a5ss_hash;
my $hashing_value = 100000;
#defaults to hg19                                                                                                                                                                                    
my $species = "hg19";
if (exists $ARGV[1]) {
    $species = $ARGV[1];
}


my $gencode_gtf_file = "/projects/ps-yeolab/genomes/hg19/gencode_v19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf";
my $gencode_tablebrowser_file = "/projects/ps-yeolab/genomes/hg19/gencode_v19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat";
#my $gencode_tablebrowser_file = "/projects/ps-yeolab/genomes/hg19/gencode_v19/gencodev19_comprehensive";                                                                                            

if ($species eq "mm9") {
    $gencode_gtf_file = "/projects/ps-yeolab/genomes/mm9/gencode.vM1.annotation.gtf";
    $gencode_tablebrowser_file = "/projects/ps-yeolab/genomes/mm9/gencode.vM1.annotation.gtf.parsed_ucsc_tableformat";
} elsif ($species eq "hg38") {
    $gencode_gtf_file = "/projects/ps-yeolab/genomes/GRCh38/gencode/v24/gencode.v24.chr_patch_hapl_scaff.annotation.gtf";
    $gencode_tablebrowser_file = "/projects/ps-yeolab/genomes/GRCh38/gencode/v24/gencode.v24.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat";
}

my %exon_triplets;
print STDERR "using $species $gencode_tablebrowser_file\n";
&read_gencode_gtf($gencode_gtf_file);
&read_gencode($gencode_tablebrowser_file);

my %dataset_counter;
#my $temp_folder = "/home/elvannostrand/data/ENCODE/RNAseq/scripts/exon_junction_counts/encode_rnaseq_psis_20161208/";
my $temp_folder = "/home/elvannostrand/data/ENCODE/RNAseq/scripts/exon_junction_counts/encode_rnaseq_psis_20170321/";

my $bam_counter=0;
#my $manifest_fi = "/home/elvannostrand/data/ENCODE/RNAseq/Brent_RNASEQlist_fromDCC.tsv.file_accessions_K562.hg19_V19.20161207.tsv";
#my $manifest_fi = "/home/elvannostrand/data/ENCODE/RNAseq/Brent_RNASEQlist_fromDCC.tsv.file_accessions_HepG2.hg19_V19.20161207.tsv";
#my $manifest_fi = "/projects/ps-yeolab3/encode/k562_brenton-graveley_ambiguous_bams_for_integrated_analysis.txt";
my $manifest_fi = "/projects/ps-yeolab3/encode/hepg2_brenton-graveley_ambiguous_bams_for_integrated_analysis.txt";

my %jxncounts_hash;
my %alreadydone;
open(M,$manifest_fi);
for my $line (<M>) {
    chomp($line);
    my @tmp = split(/\t/,$line);
    next if ($tmp[1] eq "control_rep1");
    my $encacc = $tmp[0];
    my $bam_rep1 = $tmp[5];
    my $bam_rep2 = $tmp[7];
    my $control_bam1 = $tmp[1];
    my $control_bam2 = $tmp[3];

    for my $bamfi ($bam_rep1,$bam_rep2) {
        $bamfi =~ s/\.bam$//;
	next if (exists $alreadydone{$bamfi});
	my $sorted_bam_fi = $temp_folder.$bamfi.".bam.namesort.bam";
        my $output = $sorted_bam_fi.".jxn_counts";

	unless (-e $output) {
	    print STDERR "missing output fi $output\n";
	    next;
	}

	&read_file($output,$bamfi);
	$alreadydone{$bamfi} = "expt";
    }

    for my $bamfi ($control_bam1,$control_bam2) {
        $bamfi =~ s/\.bam$//;
	next if (exists $alreadydone{$bamfi});
        my $sorted_bam_fi = $temp_folder.$bamfi.".bam.namesort.bam";
	my $output = $sorted_bam_fi.".jxn_counts";

	unless (-e $output) {
            print STDERR "missing output fi $output\n";
            next;
        }

        &read_file($output,$bamfi);
        $alreadydone{$bamfi} = "control";
    }
    $bam_counter++;
#    last if ($bam_counter > 4);
}
close(M);

my %filelist;
for my $bamfi (keys %alreadydone) {
    my $label = $alreadydone{$bamfi};
    push @{$filelist{$label}},$bamfi;
}
my $n_controls = scalar(@{$filelist{"control"}});
my $n_expts = scalar(@{$filelist{"expt"}});

print "#N controls = ".$n_controls."\n";
print "#N expts = ".$n_expts."\n";


my %flags;
my %toprint;
for my $chr (keys %exon_triplets) {
    for my $str ("+","-") {
        for my $triplet (keys %{$exon_triplets{$chr}{$str}}) {
    #	    my ($jxns,$id) = split(/\:/,$triplet);
            my $jxns = $triplet;
            my $id = join("|",@{$exon_triplets{$chr}{$str}{$triplet}});
            my ($se_jxn_up,$se_jxn_dn,$se_jxn_ex) = split(/\|/,$jxns);

            my %counts;
            for my $datatype ("control","expt") {
                for my $label (@{$filelist{$datatype}}) {

                    my ($inclusion_n,$exclusion_n) = (0,0);
                    for my $incl_jxn ($se_jxn_up,$se_jxn_dn) {
                        if (exists $jxncounts_hash{$chr.":".$str.":".$incl_jxn}{$label}) {
                            $inclusion_n += $jxncounts_hash{$chr.":".$str.":".$incl_jxn}{$label};
                        }
                    }
                    if (exists $jxncounts_hash{$chr.":".$str.":".$se_jxn_ex}{$label}) {
                        $exclusion_n = $jxncounts_hash{$chr.":".$str.":".$se_jxn_ex}{$label};
                    }

                    if ($exclusion_n + $inclusion_n >= 30) {
                        my $psi = sprintf("%.5f",($inclusion_n) / ($inclusion_n + 2 * $exclusion_n));

                        if ($exclusion_n eq "0") {
                            $counts{"strict_CE"}{$datatype}++;
                            $dataset_counter{$label.":".$datatype}++;

                        } elsif ($psi > 0.95) {
                        } elsif ($psi < 0.05) {
                        } elsif ($psi > 2/3) {
                            $counts{"incl"}{$datatype}++;
                        } elsif ($psi < 1/3) {
                            $counts{"excl"}{$datatype}++;
                        } elsif ($psi <= 2/3 && $psi >= 1/3) {
                            $counts{"center"}{$datatype}++;
                        } elsif ($psi < 0.9 && $psi > 0.1) {
                            $counts{"other"}{$datatype}++;
                        } else {
                        }
                        if ($psi >= 0.05 && $psi <= 0.95) {
                            $counts{"SE_all"}{$datatype}++;
                        }

                    }
                }
            }

            for my $type ("strict_CE","center","incl","excl","SE_all") {
                for my $datatype ("control","expt") {
                    $counts{$type}{$datatype} = 0 unless (exists $counts{$type}{$datatype});
                }
            }

    #	    for (my $j=0;$j<10;$j++) {
    #		if ($counts{"strict_CE"}{"control"} == $n_controls - $j) {
                    # exon is a strict CE in all controls                                                                           #
    #		    $flags{"strict_CE_".$j}++;
    #		}
    #	    }

    #  push @{$exons_for_position{stop}{$up_stop}},$up_start."-".$up_stop;
    #  push @{$exons_for_position{start}{$dn_start}},$dn_start."-".$dn_stop;
    # push @{$triplet2exons{$chr.":".$str.":".$se_jxn_up."|".$se_jxn_dn."|".$se_jxn_ex}},$enst."_".$i."|".$up_start."-".$up_stop."|".$ex_start."-".$ex_stop."|".$dn_start."-".$dn_stop;

            my ($current_upex,$current_uplen);
            my ($current_dnex,$current_dnlen);
            my $current_seex;
            for my $exon_triplet (@{$triplet2exons{$chr.":".$str.":".$triplet}}) {
                my ($enstinfo,$up_ex,$se_ex,$dn_ex) = split(/\|/,$exon_triplet);

                my ($upex_a,$upex_b) = split(/\-/,$up_ex);
                my ($dnex_a,$dnex_b) = split(/\-/,$dn_ex);
                my $new_uplen = $upex_b - $upex_a;
                my $new_dnlen = $dnex_b - $dnex_a;

                if ($current_upex) {
                    if ($new_uplen > $current_uplen) {
                    $current_upex = $up_ex;
                    $current_uplen = $new_uplen;
                    }

                    if ($new_dnlen > $current_dnlen) {
                    $current_dnex = $dn_ex;
                    $current_dnlen = $new_dnlen;
                    }

                    unless ($se_ex eq $current_seex) {
                    print STDERR "this shouldn't happen - different SE exons? $se_ex $current_seex $triplet\n";
                    }
                } else {
                    $current_upex = $up_ex;
                    $current_uplen = $new_uplen;
                    $current_dnex = $dn_ex;
                    $current_dnlen = $new_dnlen;
                    $current_seex = $se_ex;
                }
            }


            if ($counts{"strict_CE"}{"control"} == $n_controls) {
            # exon is a strict CE in all controls
            my $label = "strict_CE_all";
                    $flags{$label}++;
                    push @{$toprint{$label}},$chr."|".$str."|".$triplet."\t".$current_upex."\t".$current_seex."\t".$current_dnex;
                }

            if ($counts{"center"}{"control"} >= 0.50 * ($n_controls)) {
            my $label = "nSEcenter_0.5";
                    $flags{$label}++;
                    push @{$toprint{$label}},$chr."|".$str."|".$triplet."\t".$current_upex."\t".$current_seex."\t".$current_dnex;
            }

            if ($counts{"excl"}{"control"} >= 0.50 * ($n_controls)) {
            my $label = "nSEexcl_0.5";
                    $flags{$label}++;
                    push @{$toprint{$label}},$chr."|".$str."|".$triplet."\t".$current_upex."\t".$current_seex."\t".$current_dnex;
                }

            if ($counts{"SE_all"}{"control"} >= 0.50 * ($n_controls)) {
                    my $label = "nSEall_0.5";
                    $flags{$label}++;
                    push @{$toprint{$label}},$chr."|".$str."|".$triplet."\t".$current_upex."\t".$current_seex."\t".$current_dnex;
                }

            if ($counts{"incl"}{"control"} >= 0.50 * ($n_controls)) {
                    my $label = "nSEincl_0.5";
                    $flags{$label}++;
                    push @{$toprint{$label}},$chr."|".$str."|".$triplet."\t".$current_upex."\t".$current_seex."\t".$current_dnex;
                }

            if ($counts{"center"}{"control"} + $counts{"center"}{"expt"} >= 0.1 * ($n_controls + $n_expts)) {
            my $label = "aSEcenter_0.1";
                    $flags{$label}++;
                    push @{$toprint{$label}},$chr."|".$str."|".$triplet."\t".$current_upex."\t".$current_seex."\t".$current_dnex;
                }
	}


	for my $a3ss_event (keys %{$a3ss_hash{$chr}{$str}}) {
            my %counts;


	    my ($unchanged_pos,$alt3ss_basic,$alt3ss_extension) = split(/\|/,$a3ss_event);
	    my ($basic_jxn,$extension_jxn);
	    if ($str eq "+") {
		$basic_jxn = $chr.":".$str.":".$unchanged_pos."-".$alt3ss_basic;
		$extension_jxn = $chr.":".$str.":".$unchanged_pos."-".$alt3ss_extension;
	    } elsif ($str eq "-") {
		$basic_jxn = $chr.":".$str.":".$alt3ss_basic."-".$unchanged_pos;
		$extension_jxn = $chr.":".$str.":".$alt3ss_extension."-".$unchanged_pos;
	    } else {
		print STDERR "strand error xx $str\n";
	    }

	    for my $datatype ("control","expt") {
                for my $label (@{$filelist{$datatype}}) {

		    my ($basic_n,$extension_n) = (0,0);		    
		    $basic_n = $jxncounts_hash{$basic_jxn}{$label} if (exists $jxncounts_hash{$basic_jxn}{$label});
		    $extension_n = $jxncounts_hash{$extension_jxn}{$label} if (exists $jxncounts_hash{$extension_jxn}{$label});
		    
		    if ($basic_n + $extension_n >= 30) {
                        my $psi = sprintf("%.5f",($extension_n) / ($extension_n + $basic_n));
                        
			if ($psi > 0.95) {
#                            $counts{"a3ss_extension"}{$datatype}++;
                        } elsif ($psi < 0.05) {
#                            $counts{"a3ss_basic"}{$datatype}++;
			} elsif ($psi < 1/3) {
                            $counts{"a3ss_basic"}{$datatype}++; 
                        } elsif ($psi > 2/3) {
                            $counts{"a3ss_extension"}{$datatype}++; 
                        } elsif ($psi < 0.9 && $psi > 0.1) {
                            $counts{"a3ss_center"}{$datatype}++;
                        } else {
                            $counts{"a3ss_other"}{$datatype}++;
                        }
			if ($psi <= 0.95 && $psi >= 0.05) {
			    $counts{"a3ss_all"}{$datatype}++;
			}
                    }

		}
	    }
	
	    for my $type ("a3ss_center","a3ss_basic","a3ss_extension","a3ss_all") {
                for my $datatype ("control","expt") {
                    $counts{$type}{$datatype} = 0 unless (exists $counts{$type}{$datatype});
                }
            }

	    my ($up_stable,$dn_extended,$dn_basic);
	    my ($current_uplen,$current_dnlenbasic);
	    for my $exon_pairs (@{$a3ss2exons{$chr.":".$str.":".$a3ss_event}}) {
		my ($new_upstable,$new_dnextended,$new_dnbasic) = split(/\|/,$exon_pairs);
		my ($new_upa,$new_upb) = split(/\-/,$new_upstable);
		my $new_uplen = $new_upb - $new_upa;
		my ($new_dnbasica,$new_dnbasicb) = split(/\-/,$new_dnbasic);
		my $new_dnlenbasic = $new_dnbasicb - $new_dnbasica;

		if ($up_stable) {
		    if ($new_uplen > $current_uplen) {
			$up_stable = $new_upstable;
			$current_uplen = $new_uplen;
		    }
		    if ($new_dnlenbasic > $current_dnlenbasic) {
			$dn_extended = $new_dnextended;
			$dn_basic = $new_dnbasic;
			$current_dnlenbasic = $new_dnlenbasic;
		    }

		} else {
		    $up_stable = $new_upstable;
		    $dn_extended = $new_dnextended;
		    $dn_basic = $new_dnbasic;
		    $current_uplen = $new_uplen;
		    $current_dnlenbasic = $new_dnlenbasic;
		}
	    }

	    if ($counts{"a3ss_center"}{"control"} > 0.50 * ($n_controls)) {
		my $label = "nA3SScenter_0.5";
                $flags{$label}++;
                push @{$toprint{$label}},$chr."|".$str."|".$a3ss_event."\t".$up_stable."\t".$dn_extended."\t".$dn_basic;
            }

	    if ($counts{"a3ss_basic"}{"control"} > 0.50 * ($n_controls)) {
		my $label = "nA3SSbasic_0.5";
                $flags{$label}++;
                push @{$toprint{$label}},$chr."|".$str."|".$a3ss_event."\t".$up_stable."\t".$dn_extended."\t".$dn_basic;
            }
            if ($counts{"a3ss_extension"}{"control"} > 0.50 * ($n_controls)) {
		my $label = "nA3SSextension_0.5";
                $flags{$label}++;
                push @{$toprint{$label}},$chr."|".$str."|".$a3ss_event."\t".$up_stable."\t".$dn_extended."\t".$dn_basic;
            }
	    if ($counts{"a3ss_all"}{"control"} > 0.50 * ($n_controls)) {
                my $label = "nA3SSall_0.5";
                $flags{$label}++;
                push @{$toprint{$label}},$chr."|".$str."|".$a3ss_event."\t".$up_stable."\t".$dn_extended."\t".$dn_basic;
            }

	    if ($counts{"a3ss_center"}{"control"} + $counts{"a3ss_center"}{"expt"} > 0.10 * ($n_controls + $n_expts)) {
                my $label = "aA3SScenter_0.1";
                $flags{$label}++;
                push @{$toprint{$label}},$chr."|".$str."|".$a3ss_event."\t".$up_stable."\t".$dn_extended."\t".$dn_basic;
            }

	}

	
	    #a5ss2exons{$chr:$str:$event},up_basic|up_extended|dn_stable
	    #a3ss2exons{$chr:$str:$event},up_stable|dn_extended|dn_basic

	for my $a5ss_event (keys %{$a5ss_hash{$chr}{$str}}) {
            my %counts;

	    my ($alt5ss_basic,$alt5ss_extension,$unchanged_pos) = split(/\|/,$a5ss_event);
#            my ($unchanged_pos,$alt3ss_basic,$alt3ss_extension) = split(/\|/,$a5ss_event);
            my ($basic_jxn,$extension_jxn);
            if ($str eq "+") {
                $basic_jxn = $chr.":".$str.":".$alt5ss_basic."-".$unchanged_pos;
                $extension_jxn = $chr.":".$str.":".$alt5ss_extension."-".$unchanged_pos;
            } elsif ($str eq "-") {
                $basic_jxn = $chr.":".$str.":".$unchanged_pos."-".$alt5ss_basic;
                $extension_jxn = $chr.":".$str.":".$unchanged_pos."-".$alt5ss_extension;
            } else {
                print STDERR "strand error xx $str\n";
            }

            for my $datatype ("control","expt") {
                for my $label (@{$filelist{$datatype}}) {

                    my ($basic_n,$extension_n) = (0,0);
                    $basic_n = $jxncounts_hash{$basic_jxn}{$label} if (exists $jxncounts_hash{$basic_jxn}{$label});
                    $extension_n = $jxncounts_hash{$extension_jxn}{$label} if (exists $jxncounts_hash{$extension_jxn}{$label});

                    if ($basic_n + $extension_n >= 30) {
                        my $psi = sprintf("%.5f",($extension_n) / ($extension_n + $basic_n));

                        if ($psi > 0.95) {
			    #extension is always used - possible annotation error / very low abundance
#                            $counts{"a5ss_extension"}{$datatype}++;
                        } elsif ($psi < 0.05) {
			    #basic is always used - possible annotation error / very low abundance
#                            $counts{"a5ss_basic"}{$datatype}++;
			} elsif ($psi < 1/3) {
			    $counts{"a5ss_basic"}{$datatype}++; 
			} elsif ($psi > 2/3) {
			    $counts{"a5ss_extension"}{$datatype}++;   
                        } elsif ($psi < 0.95 && $psi > 0.1) {
                            $counts{"a5ss_center"}{$datatype}++;
                        } else {
                            $counts{"a5ss_other"}{$datatype}++;
                        }
			if ($psi >= 0.05 && $psi <= 0.95) {
			    $counts{"a5ss_all"}{$datatype}++;
			}
                    }

                }
            }

            for my $type ("a5ss_center","a5ss_basic","a5ss_extension","a5ss_all") {
                for my $datatype ("control","expt") {
                    $counts{$type}{$datatype} = 0 unless (exists $counts{$type}{$datatype});
                }
            }

            my ($dn_stable,$up_extended,$up_basic);
            my ($current_dnlen,$current_uplenbasic);
            for my $exon_pairs (@{$a5ss2exons{$chr.":".$str.":".$a5ss_event}}) {
                my ($new_upbasic,$new_upextended,$new_dnstable) = split(/\|/,$exon_pairs);
                my ($new_dna,$new_dnb) = split(/\-/,$new_dnstable);
                my $new_dnlen = $new_dnb - $new_dna;

                my ($new_upbasica,$new_upbasicb) = split(/\-/,$new_upbasic);
                my $new_uplenbasic = $new_upbasicb - $new_upbasica;

                if ($dn_stable) {
                    if ($new_dnlen > $current_dnlen) {
                        $dn_stable = $new_dnstable;
                        $current_dnlen = $new_dnlen;
                    }
                    if ($new_uplenbasic > $current_uplenbasic) {
                        $up_extended = $new_upextended;
                        $up_basic = $new_upbasic;
                        $current_uplenbasic = $new_uplenbasic;
                    }

                } else {
                    $dn_stable = $new_dnstable;
                    $up_extended = $new_upextended;
                    $up_basic = $new_upbasic;
                    $current_dnlen = $new_dnlen;
                    $current_uplenbasic = $new_uplenbasic;
                }
            }

            if ($counts{"a5ss_center"}{"control"} > 0.50 * ($n_controls)) {
		my $label = "nA5SScenter_0.5";
                $flags{$label}++;
                push @{$toprint{$label}},$chr."|".$str."|".$a5ss_event."\t".$up_basic."\t".$up_extended."\t".$dn_stable;
            }

	    if ($counts{"a5ss_basic"}{"control"} > 0.50 * ($n_controls)) {
		my $label = "nA5SSbasic_0.5";
                $flags{$label}++;
                push @{$toprint{$label}},$chr."|".$str."|".$a5ss_event."\t".$up_basic."\t".$up_extended."\t".$dn_stable;
            }

	    if ($counts{"a5ss_extension"}{"control"} > 0.50 * ($n_controls)) {
		my $label = "nA5SSextension_0.5";
                $flags{$label}++;
                push @{$toprint{$label}},$chr."|".$str."|".$a5ss_event."\t".$up_basic."\t".$up_extended."\t".$dn_stable;
            }

	    if ($counts{"a5ss_all"}{"control"} > 0.50 * ($n_controls)) {
                my $label = "nA5SSall_0.5";
                $flags{$label}++;
                push @{$toprint{$label}},$chr."|".$str."|".$a5ss_event."\t".$up_basic."\t".$up_extended."\t".$dn_stable;
            }
	    
	    if ($counts{"a5ss_center"}{"control"} + $counts{"a5ss_center"}{"expt"} > 0.10 * ($n_controls + $n_expts)) {
                my $label = "aA5SScenter_0.1";
                $flags{$label}++;
                push @{$toprint{$label}},$chr."|".$str."|".$a5ss_event."\t".$up_basic."\t".$up_extended."\t".$dn_stable;
            }


	    

        }


#	push @{$a3ss_hash{$chr}{$str}{$dn_start."|".$alt3ss_basic."|".$alt3ss_extension}},$ensg;
#	push @{$a5ss_hash{$chr}{$str}{$alt5ss_basic."|".$alt5ss_extension."|".$up_stop}},$ensg;

    }
}

for my $k (keys %flags) {
    print STDERR "".$k."\t".$flags{$k}."\n";
}

#for my $k2 (keys %dataset_counter) {
#    print STDERR "DATA\t$k2\t$dataset_counter{$k2}\n";
#}

for my $k (keys %toprint) {
    my $outfi = $manifest_fi.".".$k."_20170401";
    open(OUT,">$outfi");
    for my $triplet (@{$toprint{$k}}) {
	print OUT "".$triplet."\n";
    }
    close(OUT);
}

sub read_gencode_gtf {

    my $file = shift;
#    my $file = "/projects/ps-yeolab/genomes/hg19/gencode_v19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf";                                                                                      
    print STDERR "Reading in $file\n";
    open(F,$file);
    for my $line (<F>) {
        chomp($line);
        next if ($line =~ /^\#/);
        my @tmp = split(/\t/,$line);

        my $stuff = $tmp[8];
        my @stufff = split(/\;/,$stuff);
        my ($ensg_id,$gene_type,$gene_name,$enst_id,$transcript_type);

        for my $s (@stufff) {
            $s =~ s/^\s//g;
            $s =~ s/\s$//g;

            if ($s =~ /gene_id \"(.+?)\"/) {
		if ($ensg_id) {
                    print STDERR "two ensg ids? $line\n";
                }
                $ensg_id = $1;
            }
            if ($s =~ /transcript_id \"(.+?)\"/) {
                if ($enst_id) {
                    print STDERR "two enst ids? $line\n";
                }
                $enst_id = $1;
            }
            if ($s =~ /gene_type \"(.+?)\"/) {
                if ($gene_type) {
                    print STDERR "two gene types $line\n";
                }
                $gene_type = $1;

            }

            if ($s =~ /transcript_type \"(.+?)\"/) {
                $transcript_type = $1;
            }
            if ($s =~ /gene_name \"(.+?)\"/) {
                $gene_name = $1;
            }
        }



        $ensg2name{$ensg_id}{$gene_name}=1;
        $ensg2type{$ensg_id}{$gene_type}=1;

	if ($enst_id) {
            if (exists $enst2ensg{$enst_id} && $ensg_id ne $enst2ensg{$enst_id}) {
		print STDERR "error two ensgs for enst $enst_id $ensg_id $enst2ensg{$enst_id}\n";
            }

            $enst2ensg{$enst_id} = $ensg_id;
            $enst2type{$enst_id} = $transcript_type;
	}
    }
    close(F);

}


sub read_gencode {
    ## eric note: this has been tested for off-by-1 issues with ucsc brower table output!                                                                                                           \
                                                                                                                                                                                                     
    my $fi = shift;
#    my $fi = "/projects/ps-yeolab/genomes/hg19/gencode_v19/gencodev19_comprehensive";                                                                                               
    my %exon_pairs_tmp;
    my %exon_pairs_tmp_hashed;
                
    print STDERR "reading in $fi\n";
    open(F,$fi);
    while (<F>) {
        chomp($_);
        my @tmp = split(/\t/,$_);
	my $ensg = $tmp[0];
        my $enst = $tmp[1];
        next if ($enst eq "name");
        my $chr = $tmp[2];
        my $str = $tmp[3];
        my $txstart = $tmp[4];
        my $txstop = $tmp[5];
        my $cdsstart = $tmp[6];
        my $cdsstop = $tmp[7];

        my @starts = split(/\,/,$tmp[9]);
        my @stops = split(/\,/,$tmp[10]);

        my @tmp_features;


	
        my $transcript_type = $enst2type{$enst};
        unless ($transcript_type) {
            print STDERR "error transcript_type $transcript_type $enst\n";
        }
        if ($transcript_type eq "protein_coding") {

	    if (scalar(@starts) >= 3) {
		for (my $i=1;$i<scalar(@starts)-1;$i++) {
		    my ($up_start,$up_stop) = ($starts[$i-1],$stops[$i-1]);
		    my ($ex_start,$ex_stop) = ($starts[$i],$stops[$i]);
		    my ($dn_start,$dn_stop) = ($starts[$i+1],$stops[$i+1]);

		    my $se_jxn_up = $up_stop."-".$ex_start;
		    my $se_jxn_dn = $ex_stop."-".$dn_start;
		    my $se_jxn_ex = $up_stop."-".$dn_start;

		    push @{$exon_triplets{$chr}{$str}{$se_jxn_up."|".$se_jxn_dn."|".$se_jxn_ex}},$enst."_".$i;

		    
		    push @{$exons_for_position{stop}{$up_stop}},$up_start."-".$up_stop;
		    push @{$exons_for_position{start}{$dn_start}},$dn_start."-".$dn_stop;
		    
		    push @{$triplet2exons{$chr.":".$str.":".$se_jxn_up."|".$se_jxn_dn."|".$se_jxn_ex}},$enst."_".$i."|".$up_start."-".$up_stop."|".$ex_start."-".$ex_stop."|".$dn_start."-".$dn_stop;
		}
	    }

	    if (scalar(@starts) >= 2) {
		for (my $i=1;$i<scalar(@starts);$i++) {
		    my ($up_start,$up_stop) = ($starts[$i-1],$stops[$i-1]);
		    my ($dn_start,$dn_stop) = ($starts[$i],$stops[$i]);

		    $exon_pairs_tmp{$chr}{$str}{$up_start."-".$up_stop."|".$dn_start."-".$dn_stop."|".$ensg} = 1;

		    my $xx = int($up_start / $hashing_value);
		    my $yy = int($dn_stop / $hashing_value);
		    for my $jj ($xx..$yy) {
			push @{$exon_pairs_tmp_hashed{$chr}{$str}{$jj}},$up_start."-".$up_stop."|".$dn_start."-".$dn_stop."|".$ensg;
		    }
		}

	    }

#            for (my $i=0;$i<@starts;$i++) {
#            }
            for (my $i=0;$i<scalar(@starts)-1;$i++) {
                push @tmp_features,$enst."|intron|".($stops[$i]+1)."-".$starts[$i+1];
            }
        } else {

            for (my $i=0;$i<@starts;$i++) {
                push @tmp_features,$enst."|noncoding_exon|".($starts[$i]+1)."-".$stops[$i];
            }
            for (my $i=0;$i<scalar(@starts)-1;$i++) {
                push @tmp_features,$enst."|noncoding_intron|".($stops[$i]+1)."-".$starts[$i+1];
            }
        }



    }
    close(F);

    for my $chr (keys %exon_pairs_tmp) {
	for my $str ("+","-") {
	    my %alreadydone;

	    for my $exon_pair (keys %{$exon_pairs_tmp{$chr}{$str}}) {
		my ($up_ex,$dn_ex,$enst) = split(/\|/,$exon_pair);
		my ($up_start,$up_stop) = split(/\-/,$up_ex);
		my ($dn_start,$dn_stop) = split(/\-/,$dn_ex);
		my $ensg = $exon_pairs_tmp{$chr}{$str}{$exon_pair};

		my $xx = int($up_start / $hashing_value);
		my $yy = int($dn_stop / $hashing_value);
		for my $jj ($xx..$yy) {
		    for my $exon_pair2 (@{$exon_pairs_tmp_hashed{$chr}{$str}{$jj}}) {
			next if ($exon_pair eq $exon_pair2);
			my $ensg2 = $exon_pairs_tmp{$chr}{$str}{$exon_pair2};
			next unless ($ensg eq $ensg2);

			my ($up_ex2,$dn_ex2,$enst2) = split(/\|/,$exon_pair2);
			next unless ($enst2 eq $enst);
			my ($up_start2,$up_stop2) = split(/\-/,$up_ex2);
			my ($dn_start2,$dn_stop2) = split(/\-/,$dn_ex2);

			if ($up_start == $up_start2 && $dn_stop == $dn_stop2) {
			    # flanking_introns_are_same
			    ## for these - "basic" is the larger intron version, "extension" is the 5' or 3' extension version (shorter intron)
			    if ($str eq "+") {
				if ($dn_start == $dn_start2) {
				    #alt 5ss
				    my $alt5ss_basic = &min($up_stop,$up_stop2);
                                    my $alt5ss_extension = &max($up_stop,$up_stop2);
				    my $event_key = $alt5ss_basic."|".$alt5ss_extension."|".$dn_start;

#				    print "event $event_key alt5ssbasic $alt5ss_basic alt5ss_ext $alt5ss_extension up $up_ex up2 $up_ex2 dn $dn_ex\n";

				    if ($up_stop < $up_stop2) {
					push @{$a5ss2exons{$chr.":".$str.":".$event_key}},$up_ex."|".$up_ex2."|".$dn_ex;
				    } else {
                                        push @{$a5ss2exons{$chr.":".$str.":".$event_key}},$up_ex2."|".$up_ex."|".$dn_ex;
				    }
				    push @{$a5ss_hash{$chr}{$str}{$event_key}},$ensg;
				}

				if ($up_stop == $up_stop2) {
				    #alt 3ss
				    my $alt3ss_basic = &max($dn_start,$dn_start2);
                                    my $alt3ss_extension = &min($dn_start,$dn_start2);
				    my $event_key = $up_stop."|".$alt3ss_basic."|".$alt3ss_extension;
				    if ($dn_start < $dn_start2) {
					push @{$a3ss2exons{$chr.":".$str.":".$event_key}},$up_ex."|".$dn_ex."|".$dn_ex2;
                                    } else {
                                        push @{$a3ss2exons{$chr.":".$str.":".$event_key}},$up_ex."|".$dn_ex2."|".$dn_ex;
                                    }

                                    push @{$a3ss_hash{$chr}{$str}{$event_key}},$ensg;
				}
			    } elsif ($str eq "-") {
				if ($dn_start == $dn_start2) {
                                    #alt 3ss on rev strand
                                    my $alt3ss_basic = &min($up_stop,$up_stop2);
                                    my $alt3ss_extension = &max($up_stop,$up_stop2);
				    my $event_key = $dn_start."|".$alt3ss_basic."|".$alt3ss_extension;
				    if ($up_stop < $up_stop2) {
					push @{$a3ss2exons{$chr.":".$str.":".$event_key}},$dn_ex."|".$up_ex2."|".$up_ex;
				    } else {
                                        push @{$a3ss2exons{$chr.":".$str.":".$event_key}},$dn_ex."|".$up_ex."|".$up_ex2;
                                    }

                                    push @{$a3ss_hash{$chr}{$str}{$event_key}},$ensg;
				}
                                if ($up_stop == $up_stop2) {
                                    #alt 5ss on rev strand                                                                                                                           
                                    my $alt5ss_basic = &max($dn_start,$dn_start2);
                                    my $alt5ss_extension = &min($dn_start,$dn_start2);
				    my $event_key = $alt5ss_basic."|".$alt5ss_extension."|".$up_stop;

				    if ($dn_start < $dn_start2) {
                                        push @{$a5ss2exons{$chr.":".$str.":".$event_key}},$dn_ex2."|".$dn_ex."|".$up_ex;
                                    } else {
                                        push @{$a5ss2exons{$chr.":".$str.":".$event_key}},$dn_ex."|".$dn_ex2."|".$up_ex;
                                    }
                                    push @{$a5ss_hash{$chr}{$str}{$event_key}},$ensg;
				}
			    } else {
				print STDERR "error - strand err $str $exon_pair\n";
			    }
			}
		    }
		}
	    }
	}
    }

}

sub max {
    my $x = shift;
    my $y = shift;
    if ($x > $y) {
	return($x);
    } else {
	return($y);
    }
}

sub min {
    my $x = shift;
    my $y = shift;
    if ($x < $y) {
	return($x);
    } else {
	return($y);
    }
}




sub read_file {
    my $jxncounts_fi = shift;
    my $label = shift;

    print STDERR "doing $label\n";
    open(J,$jxncounts_fi);
    for my $line (<J>) {
	chomp($line);
	my ($chr,$str,$pos,$count) = split(/\t/,$line);
	$jxncounts_hash{$pos}{$label} = $count;
    }
    close(J);
}
