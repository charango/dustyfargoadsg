#!/usr/bin/perl
open VAR, "var.c";
@var = <VAR>;
close VAR;
open PARAM, ">param.h";
open PARAMNOEX, ">param_noex.h";
print PARAM <<EOF;
/********************************/
/*                              */
/* This file is created         */
/* automatically during         */
/* compilation. Do not edit.    */
/* See perl script              */
/* "varparser.pl" for details   */
/*                              */
/********************************/
EOF
print PARAMNOEX <<EOF;
/********************************/
/*                              */
/* This file is created         */
/* automatically during         */
/* compilation. Do not edit.    */
/* See perl script              */
/* "varparser.pl" for details   */
/*                              */
/********************************/
EOF
print "Running variable parser\n";
foreach (@var) {
  if (/^\s*var\("/) {
    /.*"\w+",\s(\S+?),\s(\S+?),/;
    $type = $2;
    $v = $1;
    if (($type ne "STRING") && ($v !~ /&/)) {
      die ("Non-string $v should be dereferenced. Add a '&'\n");
    }
    if (($type eq "STRING") && ($v =~ /&/)) {
      die ("String $v should not be dereferenced. Remove '&'\n");
    }
    if ($v =~ /&(\w*)$/) {$v = $1};
    $tf = "real";
    $trail = "";
    if ($type eq "INT") { $tf = "int";}
    if ($type eq "STRING") { $tf = "char"; $trail="\[512\]";}
    $line = $tf."\t".$v.$trail.";\n";
    print PARAM "extern\t".$line;
    print PARAMNOEX $line;
  }
}
close PARAM;
close PARAMNOEX;
open GLOBAL, "global.h";
@global = <GLOBAL>;
close GLOBAL;
open GLOBALEX, ">global_ex.h";
print GLOBALEX <<EOF;
/********************************/
/*                              */
/* This file is created         */
/* automatically during         */
/* compilation. Do not edit.    */
/* See perl script              */
/* "varparser.pl" for details   */
/*                              */
/********************************/
EOF
foreach (@global) {
  if (/^\s*[a-zA-Z]+/) {
    $line = $_;
    $line =~ s/\{.*?\}//g;
    $line =~ s/(\s*=.*?)(,|;)/$2/g;
    print GLOBALEX "extern ".$line;
  }
}
close GLOBALEX;
#
open SG, "sg.h";
@global = <SG>;
close SG;
open SGEX, ">sg_ex.h";
print SGEX <<EOF;
/********************************/
/*                              */
/* This file is created         */
/* automatically during         */
/* compilation. Do not edit.    */
/* See perl script              */
/* "varparser.pl" for details   */
/*                              */
/********************************/
EOF
foreach (@global) {
  if (/^\s*[a-zA-Z]+/) {
    $line = $_;
    $line =~ s/\{.*?\}//g;
    $line =~ s/(\s*=.*?)(,|;)/$2/g;
    print SGEX "extern ".$line;
  }
}
close SGEX;
