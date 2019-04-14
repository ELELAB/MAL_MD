from pymol import cmd

if not os.path.exists("../data/clusters_timeless"):
	cmd.cd("/data/user/shared_projects/chalmers/mal/simulations/dimer/analyses/9-md_MD2_2microsec_done/caver_2.0/resultsA/pymol")

cmd.cd("modules")
import caver
cmd.cd("..")

execfile('./modules/rgb.py')

color = 1
list = os.listdir("../data/clusters_timeless")
list.sort()
name = ''
for fn in list:
	old_name = name
	name = fn

	if color < 1000 and caver.new_cluster(old_name, name):
		color += 1

	cmd.load('../data/clusters_timeless/' + fn, name)
	cmd.color('caver' + str(color), name)
	cmd.alter(name, 'vdw=b')
cmd.do('set all_states,1')

cmd.set('two_sided_lighting', 'on')
cmd.set('transparency', '0.2')

cmd.load('../data/origins.pdb', 'origins')
cmd.load('../data/v_origins.pdb', 'v_origins')
cmd.show('nb_spheres', 'origins')
cmd.show('nb_spheres', 'v_origins')

cmd.load('../data/frame500.pdb.pdb', 'structure')
cmd.hide('lines', 'structure')
cmd.show('cartoon', 'structure')
cmd.color('gray', 'structure')

