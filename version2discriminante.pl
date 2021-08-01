#!perl -w
use Getopt::Long;
my %opts = ();
GetOptions(\%opts,'c=s', 'i=s', 'p=i', 'r=s', 'o=s', 'g=s', 'n=s', 'h');
if(($opts{'h'})||(scalar(keys(%opts))==0)||(!$opts{'c'})||(!$opts{'i'})){
my $os=$^O;
print"Running in: $os\n";
print <<HELP;
NAME
        version2discriminante.pl

PARAMETTERS 

	-c Colection file name.   *obligatory

	-i Grow analysis system output, a text file with x0D as line 
	   separator.(INPUT)      *obligatory

	-p Plate to analize.

	-r Plate row to start analysis.

	-h Display this.

	-n Cut point proportion n*(mean(positive control ,negative control))
	   default 1.

	-o Output name.(OUTPUT)  -obsolet(see previous version)

	-g Graphic name(OUTPUT)


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
        Output plots.pdf with the graphical view of grow.


HELP
}
else{
	$file_coleccion=$opts{'c'};
	$file_creci=$opts{'i'};
	if(!$opts{'p'}||!$opts{'r'}){
		$file_creci=~/placa(\d)(\w)(.+)/;
		$plate=$1;
		$row=uc($2);
	}
	else{
		$plate=$opts{'p'};
		$row=$opts{'r'};
	}
	if($opts{'g'}){
		$graphic=$opts{'g'};
	}
	else{
		$graphic=substr($file_creci,0,length($file_creci)-4);
	}
	if($opts{'o'}){
		$output=$opts{'o'};
	}
	else{
		$output=substr($file_creci,0,length($file_creci)-3)."out";
	}
	if($opts{'n'}){
		$nivel=$opts{'n'};
	}
	else{
		$nivel=1;
	}


#####Programa principal
	
	$mod_creci="mod_creci.txt";
	$mod_cole="mod_cole.txt";
	if($file_creci=~/.csv/){
		csv_to_txt($file_creci, '>'.$mod_creci);
	}
	else{
		cambialineas($file_creci, '>'.$mod_creci);
	}
		cambialineas($file_coleccion, '>'.$mod_cole);
	open(GROW, $mod_creci);
	@grow=<GROW>;
	close(GROW);
	open(COLE, $mod_cole);
	@col=<COLE>;
	close(COLE);
	$samples=$grow[0];
	@pozo=split("\t",$samples);
	$k=0;
	for($j=0; $j<scalar(@col);$ j++){
		@mutante=split("\t",$col[$j]);
		if($mutante[4] == $plate){
			for($i=1; $i<scalar(@pozo); $i++){
			$helpy=$pozo[$i];
				for($h=0; $h<10; $h++){
				#print"$h  ";
					$patron1="Well 1".$h;
					$patron2="Well 2".$h;
					if($h%2==0){
						$cambio="\t-\t".$h."\t";
					}
					else{
						$cambio="\t-\t".$h."\t+";
					}
					if(($row eq "A")||($row eq "B")){
						$cambio=~tr/0123456789/AABBAABBZZ/;
					}
					elsif(($row eq "C")||($row eq "D")){
						$cambio=~tr/0123456789/CCDDCCDDZZ/;
					}
					elsif(($row eq "E")||($row eq "F")){
						$cambio=~tr/0123456789/EEFFEEFFZZ/;
					}
					elsif(($row eq "G")||($row eq "H")){
						$cambio=~tr/0123456789/GGHHGGHHZZ/;
					}
					else{
						die"\nError, grong row.\n";
					}
					$cambio=~tr/+/1/;
					$cambio=~s/-/$plate/;
					if($helpy=~/$patron1/){
						$helpy=~s/$patron1/$cambio/;
						@verificador=split("\t",$helpy);						
					}
					elsif($helpy=~/$patron2/){
						$helpy=~s/$patron2/$cambio/;
						@verificador=split("\t",$helpy);
					}
					if(($col[$j]=~/$helpy/)&&($verificador[3]==$mutante[6])){
						print"pozo[$i]= $pozo[$i]\n";
						$pozo[$i]=$mutante[1];
						$indices[$k]=$i;
						$k+=1;
						last;
					}
				}
			}
		}
	}
	$to_r="creci_pl_r.txt";
	open(PARA_R, '>'.$to_r);
	for($i=0; $i<$k; $i++){
		print PARA_R "$pozo[$indices[$i]]\t";
	}
	print PARA_R "\n";
	for($i=1; $i<scalar(@grow); $i++){
		@do=split("\t", $grow[$i]);
		for($j=0; $j<$k; $j++){
			print PARA_R "$do[$indices[$j]]\t";
		}
		print PARA_R "\n";
	}
	close(PARA_R);

	for($i=1;$i<scalar(@grow);$i++){
		for($j=0;$j<scalar(@pozo);$j++){
			if($pozo[$j]=~/85/){
				@ayuda=split("\t",$grow[$i]);
				$wt[$i-1]=$ayuda[$j];
				$wt_pt[$i-1]=$ayuda[$j-3];
			}
		}	
	}
	open(WT,">wild.txt");
	open(WTP,">wildp.txt");
	foreach($i=0;$i<scalar(@wt);$i++){
		print WT $wt[$i];
		print WT "\n";
		print WTP $wt_pt[$i];
		print WTP "\n";
	}
	close(WT);
	close(WTP);



	$graphic1=$graphic."_lineplot";
	$graphic2=$graphic."_barplot";
	$R='R --save --args '.$to_r.' wild.txt wildp.txt '.$graphic1.' '.$graphic2.' '.$output.' '.$nivel.' < "version2graficante.R"';
	system($R);
	system('rm wild.txt');
	system('rm wildp.txt');
	system('rm creci_pl_r.txt');
	system('rm mod_creci.txt');
	system('rm mod_cole.txt');
}
###FUNCIONES
sub cambialineas{
	open(ARCHI, $_[0]);
	open(MODI, $_[1]);
	@archi=<ARCHI>;
	close(ARCHI);
	foreach $linea(@archi){
		@espli=split("\\x0D", $linea);
		foreach $linea_v(@espli){
			print MODI "$linea_v\n";
		}
	}
	close(MODI);
}
###
sub csv_to_txt{
	open(ARCHI, $_[0]);
	open(MODI, $_[1]);
	@archi=<ARCHI>;
	close(ARCHI);
	foreach $linea(@archi){
		$linea=~tr/,/\t/;
		print MODI $linea;
	}
	close(MODI);
}


