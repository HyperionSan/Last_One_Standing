# GCOV tool
# (1) Compile cactus with gcc using the options
#     -fprofile-arcs and -ftest-coverage
# (2) Run the test suite
# (3) from the Cactus dir, run GcovAnalysis.pl
# usage:
#   perl GcovAnalysis.pl configuration_name

use strict;
use FileHandle;
use Data::Dumper;
use Cwd;

my $cwd = Cwd::getcwd()."/";

# for debugging only
my $count = 0;
my $done = 0;

# $ARGV[0] should be a configuration name
my $dir = "configs/$ARGV[0]";

# Die if we can't find the config.
die $dir unless(-d $dir);

my %executed = ();
my $lines = {};

# search the directories and call gcov
do_dir($dir);

# do post gcov analysis
finish();

sub do_dir {
    my $fd = new FileHandle;
    my $f = shift;
    if($f =~ /\.(gcno|gcda)$/) {
        my $pre = $`;
        if(-r "$pre.gcda" and -r "$pre.gcno" and -r "$pre.o") {
            # We can only call gcov if we have both the
            # gcda and the gcno file
            gcov("$pre.o");
            #$done = 1 if($count++ == 20);
            last if($done);
        } else {
            # We are missing one of the files needed by
            # gcov. We mark line 0 as not being executed
            # in order to make the src file show up on
            # index.html. Since there is no line 0, no
            # red mark shows up when displaying the file.
            my $unused = "$cwd$pre";

            # There are no capital F's in the config dir
            $unused =~ s/\.F$/.f/;
            $unused =~ s/\.F90$/.f90/;
            $unused =~ s/\.F77$/.f77/;

            $lines->{$unused}->{0} = 0;
        }
    }
    if(-d $f) {
        opendir($fd,$f) or die $f;
        while(my $fn = readdir($fd)) {
            next if($fn eq ".");
            next if($fn eq "..");
            last if($done);
            my $full = "$f/$fn";
            # ignore automatically generated code as well as ExternalLibraries
            unless ($full =~ m'.*/cctk_Bindings$' or 
                    $full =~ m'.*/bindings$' or
                    $full =~ m'.*/scratch$' or
                    $full =~ m'.*/build/CactusBindings$') {
                do_dir($full);
            }
        }
        closedir($fd);
    }
}

sub finish {
    # Create the index.html file
    my $colors = new FileHandle;
    open($colors,">index.html");
    my $xdir = undef;
    my $table = {};
    my $fileids = {};
    my $fileid = 1;
    my $yes = 0;
    my $no = 0;

    # Count up lines with zero count (no's), or
    # more than zero (yes's). Keep track by
    # directory, etc.
    for my $src (sort keys %$lines) {
        my $cn = 0;
        my $tot = 0;
        my $lo = undef;
        my $hi = undef;
        my $lsrc = $src;
        my $pre = "";
        my $range = "";
        unless($lsrc =~ s{^${cwd}}{}) {
            #print "SKIPPING $lsrc\n";
            next;
        }
        $lsrc =~ s{^${dir}/}{};

        # Find the directory for this file.
        # The arrangement dir is special.
        next unless($lsrc =~ m{arrangements/[^/]*/[^/]*/|.*/});

        my $ndir = $&;
        my $nfile = $';

        # Assign a unique id to each file
        $fileids->{$src} = $fileid++ unless(defined($fileids->{$src}));

        for my $line (sort { $a <=> $b } keys %{$lines->{$src}}) {
            if($lines->{$src}->{$line} == 0) {
                $no++;
                $table->{$ndir}->{dir}->{no}++;
                $table->{$ndir}->{file}->{$nfile}->{no}++;
            } else {
                $yes++;
                $table->{$ndir}->{dir}->{yes}++;
                $table->{$ndir}->{file}->{$nfile}->{yes}++;
            }
        }
        my $sfd = new FileHandle;
        my $wfd = new FileHandle;
        my $outfile = sprintf("file%05d.html",$fileids->{$src});

        # Skip files that we can't locate
        open($sfd,$src) or open($sfd,"/dev/null");
        open($wfd,">$outfile") or die $outfile;
        my $lineno = 1;
        print $wfd "<html><head><style type='text/css'>\n";
        print $wfd "td { font-family: Courier; }\n";
        print $wfd ".code { background: white; }\n";
        print $wfd ".gcode { background: #eeffee; }\n";
        print $wfd ".bcode { background: #ffeeee; }\n";
        print $wfd ".gline { background: #00ff00; }\n";
        print $wfd ".bline { background: #ff0000; }\n";
        print $wfd "</style></head></body>\n";
        print $wfd "<table border=0 cellpadding=1 cellspacing=0>\n";
        print $wfd "<tr><td><b>File:</b></td><td><b>$src</b></td></tr>\n";
        while(my $line = <$sfd>) {
            my $col = "code";
            my $lcol = "code";
            if(defined($lines->{$src}->{$lineno})) {
                if($lines->{$src}->{$lineno} == 0) {
                    $col = "bline";
                    $lcol = "bcode"
                } else {
                    $col = "gline";
                    $lcol = "gcode";
                }
            }

            # Some characters are treated specially in html
            my %e = ("<" => "&lt;",">" => "&gt;","&" => "&amp;"," " => "&nbsp;","\t" => "&nbsp;");

            $line =~ s/[<>& \t]/$e{$&}/g;
            printf $wfd "<tr><td class='$col'>%5d</td><td class='$lcol'>:%s</td></tr>\n",$lineno,$line;
            $lineno++;
        }
        print $wfd "</table></body></html>\n";
        close($wfd);
    }
    my $col = get_color($yes,$no);
    my $per = percent($yes,$no);
    print $colors "<html><head><style type='text/css'>\n";
    print $colors "td { font-family: Courier; color: white; }\n";
    print $colors ".totals { background: #FFFF99; color: black; }\n";
    print $colors ".dirtotals { background: #FFFFCC; color: black; }\n";
    print $colors ".entry { background: #FFFFFF; color: black; }\n";
    print $colors "</style></head></body>\n";
    print $colors "<table border=1 cellpadding=5 cellspacing=0>\n";
    print $colors "<tr><td class='totals'><b>TOTALS</b></td><td style='background: #$col'>$per</td></tr>\n";
    for my $ndir (sort keys %$table) {
        my $lno = $table->{$ndir}->{dir}->{no};
        my $lyes = $table->{$ndir}->{dir}->{yes};
        my $col = get_color($lyes,$lno);
        my $per = percent($lyes,$lno);
        my $type = "DIR";
        if($ndir =~ m{^arrangements/}) {
            $type = "THORN";
        }
        print $colors "<tr><td class='dirtotals'><b>$type: $ndir</b></td><td style='background: #$col'>$per</td></tr>\n";
        for my $nfile (sort keys %{$table->{$ndir}->{file}}) {
            my $src = "$ndir$nfile";

            # Find the file on disk
            for my $f (("$cwd$dir$src","$cwd$src","$cwd$dir/$src")) {
                if(-r $f) {
                    $src = $f;
                }
            }

            my $outfile = sprintf("file%05d.html",$fileids->{$src});
            my $lno = $table->{$ndir}->{file}->{$nfile}->{no};
            my $lyes = $table->{$ndir}->{file}->{$nfile}->{yes};
            my $col = get_color($lyes,$lno);
            my $per = percent($lyes,$lno);
            if(defined($fileids->{$src})) {
                print $colors "<tr><td class='entry'><a href='$outfile'>$nfile</td></td><td style='background: #$col;color: #FFFFFF'>$per</td></tr>\n";
            } else {
                print $colors "<tr><td class='entry'>$nfile</td><td style='background: #$col'>$per</td></tr>\n";
            }
        }
    }
    print $colors "</table></body></html>\n";
    close($colors);
    exit(0);
}

# Call gcov on the file. Only do it once.
sub gcov {
    my $file = shift;
    return if(defined($executed{$file}));
    print "gcov $file\n";
    $executed{$file}++;
    system("gcov $file > /dev/null");
    for my $gcov (<*.gcov>) {
        analyze($gcov);
        unlink($gcov);
    }
}

# Populate the lines data structure
# with execution counts.
sub analyze {
    my $file = shift;
    my $fd = new FileHandle;
    open($fd,$file) or die $file;
    my $src = undef;
    while(<$fd>) {
        if(/Source:(.*)/) {
            $src = $1;
        } elsif(/^\s*(\d+):\s*(\d+)/) {
            $lines->{$src}->{$2} += 1*$1;
        } elsif(/^\s*#####:\s*(\d+)/) {
            $lines->{$src}->{$1} += 0;
        }
    }
}

sub get_color {
    my $yes = shift;
    my $no = shift;
    return "0000FF" if($yes == 0 and $no == 0);
    my $green = (255.0*$yes)/($yes+$no);
    my $red = 255.0-$green;
    return sprintf("%02x%02x00",$red,$green);
}

sub percent {
    my $yes = shift;
    my $no = shift;
    return sprintf("<b>%2.1f%% of %d</b>",(100.0*$yes)/($yes+$no),$yes+$no);
}
