#!perl 
use Getopt::Long;
my %opts = ();
GetOptions(\%opts,'c=s', 'i=s', 'p=i', 'r=s', 'o=s', 'g=s', 'n=s', 'h', 'v', 'f=i');
if(($opts{'h'})||(scalar(keys(%opts))==0)||(!$opts{'i'})){
my $os=$^O;
print"Running in: $os\n";
print <<HELP;
NAME
        version2discriminante.pl

PARAMETTERS 

	-i Grow analysis system output, a text file with x0D as line 
	   separator.(INPUT)      *obligatory

	-h Display this.

	-n Cut point proportion n*(mean(positive control ,negative control))
	   default 1.

	-v verbose

	-f Format to export plots, default pdf.
	   0 for pdf
	   1 for bmp
	   2 for jpeg
	   3 for png
	   4 for tiff


AUTHOR
        Roberto Galindo Ramirez
        LCG 4a Generacion
        Instituto de Fisiologia Celular
        Ciudad Universitaria UNAM

DESCRIPTION
        Takes the results from grow analysis system and obtains yeastgenome's 
        description.
        This program call to version2graficante.R in order to discriminate and 
        make graphics about strains which grow at experiment to measure.
        
        
INPUT
       Strains should be at this order into the plates:
       
       
            PROBE 1         PROBE 2    CONTROL
cols   A|C|E|G B|D|F|H A|C|E|G B|D|F|H
       01  11  21  31  41  51  61  71  81  91
       *   *   *   *   *   *   *   *   
       *   *   *   *   *   *   *   *   +
       *       *       *       *       *
       *       *       *       *       *
       *       *       *       *       -
       *       *       *       *
       *       *       *       *
       *       *       *       *
       *       *       *       *
       *       *       *       *

       
       
EXAMPLE
  perl version2discriminante.pl -c colecciona.txt -i Libro1.txt -p 1 -r A -o bambi.out -g graph
  perl version2discriminante.pl -c directory/colecciona.txt -i placa2a1b12rep.txt -n 0.75 
	
        
FORMAT 
        Output graphical files with the graphical view of grow.
	If there is almost one strain with meanly grow write into a plain text
	file strain names as a list.

HELP
}
else{
	$fname=$opts{'i'};
	print"\nEscribe los nombres de los mutantes por columna,\nen caso de no haber mutante solo da enter\n\n";
	$tname=substr($fname,0,length($fname)-4);
	open(FI, $fname);
	@fi=<FI>;
	close(FI);
	@headers=split(",",$fi[0]);
	foreach$pozo(@headers){
		$pozo=~s/Well (\d)(\d)/_/;
	}
	$mutantes[0]="Tiempo";
	$k=0;
	for($j=1; $j<scalar(@headers); $j++){
		if(($j-1)%10==0){
			$col=(($j-1)/10)+1;
			print"Cual es la mutante de la columna $col ?\n";
			$help=<STDIN>;
			chop($help);
		}
		if($headers[$j] eq "_0"){
			$headers[$j]="_10";
		}
		$mutantes[$j]=$help.$headers[$j];
		$mutantes[$j]=~/(.+)_(\d+)/;
		if($1 ne ""){
			#print"$1\t";
			next;
		}
		else{
			$empty[$k]=$j;
			$k+=1;
		}
	#print"$mutantes[$j]\n";	
	}
	print"$k\n";
	for($j=1;$j<scalar(@fi);$j++){
		@slot=split(",",$fi[$j]);
		@trueslot=();
		for($l=1;$l<scalar(@slot);$l++){
			$flag=0;
			for($i=0;$i<$k;$i++){
				if($l==$empty[$i]){
					$flag=1;
					last;				
				}
				else{
					next;
				}
			}
			if($flag==0){
				$ayuda[0]=$slot[$l];
				push(@trueslot,@ayuda);
			}
		}
		$fi[$j]=join(",",@trueslot);
	}
#	print"@mutantes\n";
	$mut[0]=join(",",@mutantes);
	$mut[0]=~s/,_(\d+)//g;			#la g sirve para hacerlo global
	$mut[0]=~s/Tiempo,//;
	$mut[0]=~s/\n//;
	$x=scalar(@mutantes)-1;	
	print"@mut \n Longitud de mutantes : $x\n";			
	shift(@fi);
	$y=scalar(split(",",$fi[2]));
	print"Longitud de \$fi[2]  despues del shift = $y\n"; 
	unshift(@fi,@mut);
#	print"$fi[0]\n";
	
	$impr=join("\n",@fi);		###Aqui no ponemos el ultimo enter del archivo, causara problemas??   NOOO
	open(TMP, '>'.$tname.'_tmp.csv');
	print TMP "$impr\n";
	close(TMP);
	$R='R --slave --args '.$tname.'_tmp.csv < "mutanteColumna.R"';
	system($R);
	#system('rm ');
	
}
