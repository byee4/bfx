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
# print STDERR "using $species $gencode_tablebrowser_file\n";
&read_gencode_gtf($gencode_gtf_file);
&read_gencode($gencode_tablebrowser_file);

my %dataset_counter;
#my $temp_folder = "/home/elvannostrand/data/ENCODE/RNAseq/scripts/exon_junction_counts/encode_rnaseq_psis_20161208/";
my $temp_folder = "/projects/ps-yeolab3/bay001/tmp/eric_jxc/";

my $bam_counter=0;
#my $manifest_fi = "/home/elvannostrand/data/ENCODE/RNAseq/Brent_RNASEQlist_fromDCC.tsv.file_accessions_K562.hg19_V19.20161207.tsv";
#my $manifest_fi = "/home/elvannostrand/data/ENCODE/RNAseq/Brent_RNASEQlist_fromDCC.tsv.file_accessions_HepG2.hg19_V19.20161207.tsv";
#my $manifest_fi = "/projects/ps-yeolab3/encode/k562_brenton-graveley_ambiguous_bams_for_integrated_analysis.txt";
# my $manifest_fi = "/projects/ps-yeolab3/encode/hepg2_brenton-graveley_ambiguous_bams_for_integrated_analysis.txt";
my $manifest_fi = "/home/bay001/projects/codebase/bfx/from_eric/rna_seq_manifest.txt.2";

my %jxncounts_hash;

my %alreadydone;
# print "Parsing manifest";
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

# print "#N controls = ".$n_controls."\n";
# print "#N expts = ".$n_expts."\n";


my %flags;
my %toprint;

# print "Calculating psi values";

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
                        print "${psi}\t${chr}:${str}:${label}:${exon_triplets{$chr}{$str}{$triplet}}\n";
                    }
                }
            }
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
        }
    }
}

sub read_gencode_gtf {

    my $file = shift;
#    my $file = "/projects/ps-yeolab/genomes/hg19/gencode_v19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf";                                                                                      
    # print STDERR "Reading in $file\n";
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
                
    # print STDERR "reading in $fi\n";
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

    # print STDERR "doing $label\n";
    open(J,$jxncounts_fi);
    for my $line (<J>) {
	chomp($line);
	my ($chr,$str,$pos,$count) = split(/\t/,$line);
	$jxncounts_hash{$pos}{$label} = $count;
    }
    close(J);
}
