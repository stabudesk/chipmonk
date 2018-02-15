# chipmonk (previous name, bedtack)

Chipmonk is the new name for this project. It is a more general name indicating that ChIP analysis is the final aim. Obviously it would have been cuter to call it chipmunk, but that is also very obvious and probably that name has been used several times (though I didn't check). Also "monk" lends an authoritative air to the program in the way that perlmonks website achieves. Yes, it's fuzzy, but this allows functioanl flexibility.

What follows is a discussion on the previous name, obsolete.

The name was chosen principally because it "tackles" BED files and the name BEDtools and BEDutils imply a much bigger platform of tools while this project is rather small.

The name also has a fuzzy evocation of "BED attack" or "BED tuck" which all imply a certain handling of the BED file format which actually varies quite a bit, principally on the number of columns it has.

If you do not have a BED file, but rather have a GFF file, there are tools out there to do the conversion, so it's not part of the current scope.

## 0 indexing
GFF uses 1-indexing. The range is closed. Start and end are both part of the interval. However, internally, the program converts these
to zero indexing, which is the bed convention. Be careful when eyeballing gffs and comparing to program output.
BED files use 0 indexing (semi-open ranges) , as does (by default) the genome size files.
This means semi-open intervals. The first number is inclusive of the range, while the second is one above the final included position.

When ranges appear on the ID line of outputted fasta files, they use the bed convention.

## Sometimes chromosome order 4, 9, 5
Lexicographic ordering means 4,9,5 because Roman numerals are used for chromosome names

## Other times normal
this is because samtools depth is changes to numeric ordering and IX follows VIII not IV

# introducing gfmatchup
Attempting to match up gff's to see which features are in which.

# what's difference between chipmonk and chipmonk2
Well, the key difference is that chipmonk2 is the more modern, and so it's recommended.

# First run and test files

These are the different yeast reference genomes.
Particularly the unannotated W303LYZE
so in as363's gff file we have five entries for chr18 (although with subnumbers 1-5). Only chr 1-17 are true chromosomes really.
The reference not only has 18, but also 19, 20, and 21. QUite likely these are contigs that could
not not be packaged into the other chromosomes.
What we shall do is take the 5 annotations in chr18, and render their subnumbers as 1.

# data files
Beware W303\_LYZE\_ioptrf\_v2.gff, this is an output of chip monk
there isn't really a way of inputting it,
use the raw  W303\_LYZE\_IGVoptimised\_rf.gff
instead

# Activities
uf.depth, large 16M depth positions across S288 of some sort .. had to delete it

# matching gf22 and blast output
hurry and preempting led me to match gf22 against target ... which of course in many occurences is fine. However, it's lunacy. These are blast results, you need to match gf22 ID strings to the query's ID. This will ofcourse give multiple matches, and then you use eval pcti to select which you print.

# gffsimp
This was coded for more complicated gff files, namely the salmo salar gff file. The Genbank gene name unusually includes a .1 in its name, they all do. This creates problems
as they need to be stripped off. Also this code by default focusses only on mRNA entries in the GFF.

# Worries
chipmonk2 has got too big? Even executing it is difficult. I.e.
> Famous errors: presenting -ystr to the executable ... it gives a plain segfault no info. obvious it shoudl be -y!
> presenting the same gff file to both -a and -y option ... not looking properly at args .. another blank segfault.

# Making the blast output hashtable
Funnily enough only the first entries were stored anyway. Kind of unusual
I need to go over why the other were rejected.
It's actually something that I wanted to do, but I need to discover why it's being done automatically.
this is easy enough to discover the specifically rejects all except first.

# the issue of strand-sense
when matching of a fasta file, that only has the one forward sense. If you have feature which is Crick sense, you decide to translate or not. My convention is to translate, and to accompany
the annotation somehow with the "-" reverse strand symbol. IS there scope for confusion? Of course, bu the idea and the convention I think, is that a fasta should have the literal sequence. 

# feature name encoding
some of this has turned out to have percetn-econding: %20 for spaces, %28 for ( and %29 for ) 

# Is there a need for a blast output plus gff processor?
I would have thought that this is obvious ... it is. Why, you want to know the feature associated with the hit.
Well, I suppose you need to google to see who has dal with it ... surely alot of people?
* https://www.researchgate.net/post/How_to_retrieve_the_full_annotation_in_blast_report-it_only_gives_accession_number
* https://github.com/nterhoeven/blast2bed
