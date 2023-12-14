#!/usr/bin/perl -w
use strict;
use Switch;
# author: Frank Löffler

my %headings = ("t" => "Thorn",
                "a" => "Author(s)",
                "m" => "Maintainer(s)",
                "l" => "Licence",
                "p" => "Purpose");

# get column specification from stdin or take default
my @columns = split("", "tap");
if (defined($ARGV[0])) {
  @columns = split("", $ARGV[0]);
}

# get a list of all README files
my @readmes = `find -L -maxdepth 3 -mindepth 3 -name README`;

print "<!-- Do not edit this file directly. It is automatically generated.\n".
      "     All changes will be overwritten. Resistance is futile. -->\n";

# save the table headers for later use
my $table_header = "<table border=1>\n";
foreach my $c (@columns) {
  $table_header.="<th>".$headings{$c}."</th>";
}
$table_header.="\n";

my $old_arrangement = "start";

sub name_list() {
  my ($name_str) = @_;
  my $ret = "<td>";
  foreach my $name (split("\n", $name_str)) {
    chomp($name);
    $name =~ s/^\s+//;
    $name =~ s/\s+$//;
    $name =~ s/ /&nbsp;/g;
    $ret.= "$name<br>\n";
  }
  $ret.="  </td>\n";
  return $ret;
}

# loop over all README files
foreach my $readme_filename (sort @readmes) {
  # parse for the arrangement and thorn name
  my $arrangement;
  my $thornname = "unknown";
  if ($readme_filename =~ /\.\/([-_a-zA-Z0-9]+)\/([-_a-zA-Z0-9]+)\/README/ ) {
    $arrangement = $1;
    $thornname = $2;
  }
  if ($thornname eq "unknown") { die("Could not parse thorn name"); }
  # close and reopen table when we arrive at a new arrangement
  if ($old_arrangement ne $arrangement) {
    if ($old_arrangement ne "start") { print "</table>\n"; }
    print "<h3>$arrangement</h3>\n$table_header";
  }
  # finally read the README file
  open README, $readme_filename or die ("Could not open $readme_filename");
  my $readme_str = join("", <README>);
  close README;
  # parse it with a single regexp: could be prettier
  if ($readme_str =~ /^Cactus Code Thorn $thornname\nAuthor\(s\)    :((.+\n)+?)Maintainer\(s\):((.+\n)+?)Licence      : (.*)\n-{74}\n\n1. Purpose\n\n((.+\n)+?)((\n\n)|$)/m) {
    my $author_str     = $1;
    my $maintainer_str = $3;
    my $licence_str    = $5;
    my $purpose_str    = $6;
    print " <tr>\n";
    foreach my $c (@columns) {
      switch ($c) {
        case "t" { print "<td>$thornname</td>"; }
        case "a" { print &name_list($author_str); }
        case "m" { print &name_list($maintainer_str); }
        case "l" { print "<td>$licence_str</td>\n"; }
        case "p" { print "<td>$purpose_str</td>\n"; }
      }
    }
    print " </tr>\n";
  }
  else {
    # die("\nError: $arrangement/$thornname not parsable");
    print " <tr>\n";
    foreach my $c (@columns) {
      switch ($c) {
        case "t" { print "<td>$thornname</td>"; }
        else     { print "<td>not parseable</td>"; }
      }
    }
    print " </tr>\n";
  }
  $old_arrangement = $arrangement;
}
print "</table>\n";
