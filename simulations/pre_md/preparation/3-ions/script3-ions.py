# remove waters
cmd.remove("resn SOL")

s=cmd.select("sel1","resn Na")
cmd.select("sel2","resn Cl")

if (s==0):
	t = cmd.select("sel2 around 4")
        if(t==0):
                print "				We have a winner! Structure has been selected."
                os._exit(1)
        else:
                print "				This is bad...not far enough!!"
                os._exit(0)
else:
	t = cmd.select("sel1 around 4")
	if(t==0):
		print "				We have a winner! Structure has been selected."
		os._exit(1)
	else:
		print "				This is bad...not enough distant!"
		os._exit(0)
