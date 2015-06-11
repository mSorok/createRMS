/*
 * SSCAN PROGRAM Version 1.00
 * BioRetroSynth Lab - Institute of Systems and Synthetic Biology 
 * Nov 2012
*/

Stereo Signature Molecular Descriptor

Usage: sscan <file.mol> <output-type>
       <file> is mol file 
       <output-type> 
output types are (x is the diameter, any integer between 0 and 100000):
	scanx: signature DAG
	sscanx: stereo signature DAG (CIP rules)
	fscanx: fast stereo signature mode
with scan and sscan output modes:
	if x is omitted then the code will canonicalize the structure
	and print the maximum signature (according to the invariant computed from the CIP rules).
with sig and ssig output modes:
	if x is omitted, x takes the value 0
	if x is greater than 12, x takes the value 12.

