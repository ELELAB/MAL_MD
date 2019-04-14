#set dir "/data/user/shared_projects/chalmers/mal/simulations/dimer/analyses/9-md_MD2_2microsec_done/caver_2.0/pdbA"

mol load pdb ../data/frame500.pdb.pdb

after idle { 
  mol representation NewCartoon 
  mol delrep 0 top
  mol addrep top
  mol modcolor 0 top "ColorID" 8
} 

