#!usr/bin/perl
#questo programma divide un file contenente più pdb
#(ogni PDB inizia con "TITLE...." e finisce con "ENDMDL")
#in file contenenti un singolo pdb ognuno



open (INPUT,'models.pdb');
while ($riga=<INPUT>)
        {chomp $riga;

	if ($riga=~/^(TITLE).*([0-9]\.)/)

		{open (OUTPUT,">>".$2."pdb");
		print OUTPUT "$riga\n";
		}
        else {if ($riga=~/^(ENDMDL)/)
		{print OUTPUT "$riga\n";
		close OUTPUT;
		}
		else {print OUTPUT "$riga\n" }
	     }	


	}




