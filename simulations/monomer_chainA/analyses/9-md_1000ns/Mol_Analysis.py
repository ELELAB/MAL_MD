#    Copyright (C) 2018 Matteo Lambrughi, Matteo Tiberti, Maria Francesca Allega, Valentina Sora, Mads Nygaard, Agota Toth, Juan Salamanca Viloria, Emmanuelle Bignon, Elena Papaleo <elenap@cancer.dk>, Computaitonal Biology Laboratory, Danish Cancer Society Research Center, 2100, Copenhagen, Denmark
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#!/usr/bin/env python
__author__ = "Mads Nygaard, Agota Toth"
__date__ = "20171211"

### Global dependencies
import os
import matplotlib as mpl
if os.environ.get('DISPLAY') is None:
    mpl.use("Agg")

import matplotlib.pyplot as plt
import gromacs
import numpy as np
from copy import deepcopy

### Monkey patch for mpi
if False:
    def rename_attribute(object_, old_attribute_name, new_attribute_name):
        setattr(object_, new_attribute_name, getattr(object_,
                old_attribute_name))

    rename_attribute(gromacs.tools, "Make_ndx_mpi", "Make_ndx")
    rename_attribute(gromacs.tools, "Trjconv_mpi", "Trjconv")
    rename_attribute(gromacs.tools, "Gmxcheck_mpi", "Gmxcheck")
    rename_attribute(gromacs.tools, "G_mindist_mpi", "G_mindist")
    rename_attribute(gromacs.tools, "Editconf_mpi", "Editconf")
    rename_attribute(gromacs.tools, "G_rms_mpi", "G_rms")
    rename_attribute(gromacs.tools, "G_gyrate_mpi", "G_gyrate")
    rename_attribute(gromacs.tools, "G_rmsf_mpi", "G_rmsf")


### Classes

### Extended dict with variables for diff variables

class fnamesclass(dict):
    def set_group(self, group):
        self.group = group

    def get_group(self):
        return self.group

    def set_folder(self, folder):
        self.folder = folder

    def get_folder(self):
        return self.folder

    def set_chains(self, chains):
        if type(chains) is not list:
            self.chains = [chains]
        else:
            self.chains = chains[:]

    def get_chains(self):
        return self.chains[:]

    def merge(self, fnamesobj):
        pass

    def __repr__(self):
        return "Group: %s , Folder: %s " % (self.group, self.folder) + \
            super(fnamesclass, self).__repr__()

#    def __getitem__(self, y):
#        item = super(fnamesclass, self).__getitem__(y)
#        temp = fnamesclass({y: item})
#        temp.set_folder(self.get_folder())
#        temp.set_group(self.get_group())
#        return temp


class gTools_runner:
    """Runs gromacs wrapper, keeps track of input/output filenames etc."""
    def __init__(self, wdir, odir="Mol_An", splitlen=[10], f=None, tpr=None, ndxfile=None,
                 outf=None, outs=None, chains=None, center=None,
                 ndxparams=None, dryrun=False, skip=None, verbose=False,
                 dt=None, group=0):

        # Globally sets verbose mode
        gromacs.environment.flags['capture_output'] = not verbose
        # Store different variables of self
        self.splitlen = splitlen
        self.dryrun = dryrun
        if not os.path.isdir(wdir):
            raise IOError("Folder:%s is not found" % wdir)
        else:
            self.wdir = os.path.abspath(wdir) + "/"

        realodir = wdir+"/"+odir

        # Create output directory if not present
        if not os.path.isdir(realodir):
            os.makedirs(realodir)

        for i in splitlen:
            if not os.path.isdir(realodir+"/rmsf"+str(i)):
                os.makedirs(realodir+"/rmsf"+str(i))

        if (os.path.isfile(wdir+"/residuetypes.dat")) and \
                (not os.path.isfile(os.getcwd()+"/residuetypes.dat")):
            import shutil
            print("Copy residuetypes.dat to current dir. Can be removed after"
                  " script is completed.")
            shutil.copyfile(wdir+"/residuetypes.dat", os.getcwd() +
                            "/residuetypes.dat")

        self.odir = os.path.abspath(realodir) + "/"
        self.ndxfile = self.odir + ndxfile
        self.outf = self.odir + outf
        self.outs = self.odir + outs
        self.f = self.wdir + f
        self.tpr = self.wdir + tpr
        self.logfile = self.odir + "file.log"
        self.errfile = self.odir + "errfile.log"
        self.fnames = fnamesclass()
        self.fnames.set_group(group)
        self.fnames.set_folder(self.odir)
        
        if center is None:
            self.center = "Protein"
        else:
            self.center = center
        self.skip = skip
        self.dt = dt
        if chains is None:
            self.fnames.set_chains("Protein")
        else:
            self.fnames.set_chains(chains)
        # Try to guess the ndx parameters if missing from config
        # (high chance of faliure here)
        self.ndxparams = self._guessndxparam(ndxparams)

    def _checkpath(self, fname, wdir):
        """Check if path is absolute, not used..."""
        if fname[0] == "/":
            return os.path.abspath(fname)
        else:
            return os.path.abspath(wdir + fname)

    def _guessndxparam(self, ndxparams):
        """Very crude script to guess parameters for creation of ndx file.
        Returns list of params. Self.chains, self.center needs to be
        specified"""
        
        params = []
        if ndxparams is not None: 
            if type(ndxparams) is not list:
                ndxparams = [ndxparams]

            for param in ndxparams:
                params.append(param.lower())
                params.append(param.lower() + " & 3")

        for chain in self.fnames.get_chains():
            if not any(n in chain.lower() for n in ["ch", "protein"]):
                print "Something fishy is going on with the chains"
            elif "protein" in chain.lower():
                params.append("1 & 3")
            else:
                params.append("chain " + chain.lower().split("ch ")[-1])
                params.append("chain " + chain.lower().split("ch ")[-1] + " & 3")
        params.append(self.center.replace("_", " ").replace("C-alpha", "3"))
        return list(set(params)) 
        
    def _timesplit(self, newdt, stop, start=0):
        """Returns list of (start,stop) points to split a trajectory
        based on the newdt value."""
        splits = []
        if stop is None:
            with open(self.odir + "r_gyrate.xvg", "r") as f:
                lastline = f.readlines()[-1]
                stop = int(float(lastline.split()[0]))
                self.simlen = stop
        if newdt < self.dt:
            return ((None, None), )
        else:
            splits.append((start, (start+newdt)))
            for n in range(start+newdt, stop, newdt):
#		print self.trjdt
                splits.append((n+self.trjdt, (n+newdt)))
            return splits

    def _errorcatch(self, excode, stdo, erro, funname):
        """Not used"""
        print erro
        print stdo
        print "Something went wrong, %s exited with errorcode %s"
        raw_input("Press any key to continiue")

    def _addfilename(self, key, fname):
        """Add filename to the filenamedict"""
        if fname not in self.fnames.get(key, []):  # Add filename if unique
            self.fnames.setdefault(key, []).append(fname)

    def getfilenames(self):
        """Get filenamedict"""
        return deepcopy(self.fnames)

    def getodir(self):
        """Get output dir location(str)"""
        return self.odir

    def getinfo(self, fname):
        """Run gromacs check to get info from traj.
        Returns (nframes,dt, simlen(fs))"""
        g_check = gromacs.tools.Gmxcheck(
                f=fname,
                stdout=False,
                stderr=False
                )
        print "Getting info from traj"
        excode, out, err = g_check.run()
        for line in err.split("\n"):
            sline = line.split()
            if all(k in sline for k in ("don't", "match")):
                print " ".join(sline)
                import re
                if len(re.findall(r"\((\d+)[ ,]", " ".join(sline))) > 0:
                    self.trjdt = int(re.findall(r"\((\d+)[ ,]", " ".join(sline))[0])
            if len(sline) > 1:
                if "Time" == sline[0]:
                    comments = map(int, sline[1:])
                    break
        self.frames = comments[0]
        if len(comments) > 1:
            self.trjdt = comments[1]
            self.simlen = (comments[0]-1)*comments[1]
        else:
            self.simlen = None
        return (self.frames, self.trjdt, self.simlen)

    def mkndx(self, commands, ffile=None, **kwargs):
        """Runs gromacs make_ndx with commands(list)"""
        if commands[-1].lower() is not "q":
            commands = list(commands)
            commands.append("q")
        if ffile is None:
            ffile = self.tpr
        g_mkndx = gromacs.tools.Make_ndx(
                f=ffile,
                o=self.ndxfile,
                input=commands
                )
        if not self.dryrun:
            print "Making ndx file"
            g_mkndx.run(**kwargs)

    def mindist(self, **kwargs):
        """Runs gromacs mindist on the new, centered traj."""
        outfile = self.odir + "min_pbc_dist.xvg"
        self._addfilename("mindist", outfile)
        g_mindist = gromacs.tools.G_mindist(
                f=self.outf,
                s=self.tpr,
                pi=True,
                od=outfile,
                input="Protein",
                stdout=False,
                stderr=False
                )

        if not self.dryrun:
            print "Running mindist"
            excode, out, err = g_mindist.run(**kwargs)
        if not self.dryrun:
            for line in out.split("\n"):
                sline = line.split()
                if len(sline) > 2:
                    if ["The", "shortest", "periodic"] == sline[0:3]:
                        self.minpdist = float(sline[5])
                        self.mintime = float(sline[9])
                        break

    def editconf(self, **kwargs):
        """Generates a new conf file"""
        g_editconf = gromacs.tools.Editconf(
                f=self.tpr,
                n=self.ndxfile,
                o=self.outs,
                input="Protein",
                )
        if not self.dryrun:
            print "Updating conf"
            excode, out, err = g_editconf.run()

    def trjconv(self, **kwargs):
        """Generates a new centered trjconv file"""
        center = self.center
        g_trjconv = gromacs.tools.Trjconv(
                f=self.f,
                s=self.tpr,
                o=self.outf,
                n=self.ndxfile,
                skip=self.skip,
                dt=self.dt,
                pbc="cluster",
                ur="compact",
                center=True,
                input=("Protein", center, "Protein"),
                )  # trjconv -f traj.trr -s npt.tpr -dump time -o confout.dump.gro -sep -pbc mol -ur compact
        if not self.dryrun:
            print "Centering trajectory"
            excode, out, err = g_trjconv.run()
        self.getinfo(self.outf)

    def rmsd(self, **kwargs):
        """Calcs the rmsd for the two chains"""
        chains = self.fnames.get_chains()
        for chain in chains:
            outfile = self.odir + "rmsd_%s.xvg" % (chain)
            self._addfilename("rmsd", outfile)
            g_rmsd = gromacs.tools.G_rms(
                    f=self.outf,
                    s=self.outs,
                    o=outfile,
                    n=self.ndxfile,
                    input=(chain, chain),
                    )
            if not self.dryrun:
                print "Calculating RMSD"
                excode, out, err = g_rmsd.run()

    def rg(self, **kwargs):
        """Calcs the rg for the complete protein"""
        outfile = self.odir + "r_gyrate.xvg"
        self._addfilename("rg", outfile)
        g_rg = gromacs.tools.G_gyrate(
                f=self.outf,
                s=self.outs,
                o=outfile,
                input="protein",
                )
        if not self.dryrun:
            print "Calculating Rg"
            excode, out, err = g_rg.run()

    def rmsf(self, splitlen, **kwargs):
        """Calcs the rmsf broken into chunks of splitlen(ns) size"""
        #chains = [x for x in self.fnames.get_chains() if "C-alpha" in x]
        chains = self.fnames.get_chains()
        splitlenfs = splitlen*1000
        try:
            splits = self._timesplit(splitlenfs, self.simlen)
        except AttributeError:
            self.getinfo(self.outf)
            splits = self._timesplit(splitlenfs, self.simlen)
        for chain in chains:
            chain = chain + "_&_C-alpha"
            chain_fname = chain.replace("&", "n")  # Fix for not using & in filename
            for begin, end in splits:
                outfile = self.odir + "rmsf"+str(splitlen) + "/rmsf_%s_%s.xvg" % (chain_fname, end)
                self._addfilename("rmsf" + str(splitlen) + "_" + chain_fname, outfile)
                g_rmsf = gromacs.tools.G_rmsf(
                        f=self.outf,
                        s=self.outs,
                        o=outfile,
                        od=self.odir + "rmsf"+str(splitlen) + "/rmsf_od_%s_%s.xvg" % (chain_fname, end),
                        b=begin,
                        e=end,
                        n=self.ndxfile,
                        res=True,
                        input=(chain, chain),
                        )
                if not self.dryrun:
                    print "Calculating RMSF t %s to %s" % (begin, end)
                    excode, out, err = g_rmsf.run()

    def pdbout(self, pdbname="pdbmovie.pdb", **kwargs):
        """Outputs a pdb trajectory, protein fitted of the centered traj"""
        g_pdbout = gromacs.tools.Trjconv(
                f=self.outf,
                s=self.outs,
                o=self.odir + pdbname,
                fit="rot+trans",
                input=["Protein", "Protein"]
                )
        if not self.dryrun:
            print "Printing .pdb movie"
            excode, out, err = g_pdbout.run(**kwargs)

    def runall(self, **kwargs):
        self.torun = {"mkndx": True, "mindist": True, "trjconv": True,
                      "editconf": True, "rmsd": True, "rg": True, "rmsf": True,
                      "pdbout": True}
        self.torun.update(kwargs)

        if self.torun["mkndx"] is True:   # To create the .gro file
            self.mkndx(commands=self.ndxparams, ffile=self.tpr)
            #print self.ndxparams
        if self.torun["editconf"] is True:
            self.editconf()
#        if self.torun["mkndx"] is True:          #Consisting of only protein
#            self.mkndx(commands = self.ndxparams, ffile = self.outs)
        if self.torun["trjconv"] is True:
            self.trjconv()
        else:
            self.getinfo(self.outf)

        if self.torun["mindist"] is True:
            self.mindist()

        if self.torun["rmsd"] is True:
            self.rmsd()

        if self.torun["rg"] is True:
            self.rg()

        if self.torun["rmsf"] is True:
            for i in self.splitlen:                
                self.rmsf(i)

        if self.torun["pdbout"] is True:
            self.pdbout(skip=20)


class gTools_plotter:
    """Plotting tool for xvg files made with Gromacs.
    Uses the gromacs wrapper"""
    class xvg_read(gromacs.fileformats.XVG):
        """Updated class of the XVG reader to handle current axes better"""
        def plot(self, ax, **kwargs):
            self.fig.sca(ax)
            gromacs.fileformats.XVG.plot(self, **kwargs)
            return self.fig.gca()

    def __init__(self, fnames, odir, begin=0, end=None, splitlen=[10], chains=["Protein"], figsize=(8.27, 11.69)):
        """Makes the figure etc"""
        from cycler import cycler
        self.fnames = fnames
        self.odir = odir
        self.begin = begin
        self.end = end
        self.splitlen = splitlen
        self.fig = plt.figure(figsize=figsize)
        plt.rc('axes', prop_cycle=(cycler('color', ['k', 'r', 'g', 'b', 'm', 'y', 'c'])))
        self.axlist = []
        print self.odir
    
    def _filter_rmsfiles(self, flist, begin, end):
        # Removing rmsf times if begin is set
        if end is not None:
            shortend_list = [f for f in flist if (int(
                            re.search("(\d+)\.xvg", f).group(1)) >= begin) and 
                            (int(re.search("(\d+)\.xvg", f).group(1)) <= end)]
        else:
            shortend_list = [f for f in flist if (int(
                            re.search("(\d+)\.xvg", f).group(1)) >= begin)]

        if len(shortend_list) <= 0:
            print("Less than 0 rmsf files after filtering, using original list")
            shortend_list = flist
        return shortend_list 

    def resplotstd(self, flist, label, ax=None, std_fill=True, std_bar=True,
                   shift=0, **kwargs):
        """Plots residues vs val. Calculates the stddev for the filelist"""
        if ax is None:
            ax = plt.gca()
        filter_rmsflist = self._filter_rmsfiles(flist=flist, begin=self.begin, end=self.end)
        data, xaxis, yaxis = self.stddata(filter_rmsflist)
        upperstd = data[1]+data[2]
        lowerstd = data[1]-data[2]
        data[0] = data[0]+shift
        lines = ax.plot(data[0], data[1], label=label, **kwargs)
        plotcol = lines[-1].get_color()
        if std_fill:
            ax.fill_between(data[0],
                            upperstd,
                            lowerstd,
                            alpha=0.1,
                            facecolor=plotcol,
                            lw=0)
        if std_bar:
            ax.bar(data[0], data[2], color=plotcol)
        ax.set_xlabel(xaxis)
        ax.set_ylabel(yaxis)
        #ax.set_ylim((0, upperstd[2:-3].max()))
        ax.set_xlim((data[0].min(), data[0].max()))
        ax.tick_params(direction='out') 
        ax.yaxis.tick_left()
        ax.xaxis.tick_bottom()

        # ax.set_ylim((0, 1.2))
#        if "chB" in label:
#            ylims = [0, 0.3]
#        elif "chA" in label:
#            ylims = [0, 0.45]
#        ax.set_ylim(ylims)

        # self.axlist[-1].xaxis.set_minor_locato(AutoMinorLocator(1))

    def stddata(self, filelist):
        """Calculates the stddev for the files in the filelist.
        Returns (np.array(resnr|data|stddev),xlabel,ylabel)"""
        datalist = []
        for fname in filelist:
            if not os.path.isfile(fname):
                fname = fname.replace("_n_", "_&_")
            data = self.xvg_read(fname)
            datalist.append(data.array[1])
        xaxis, yaxis = data.xaxis, data.yaxis
        dataarray = np.array(datalist)
        arrayout = np.array((data.array[0], dataarray.mean(axis=0),
                             dataarray.std(axis=0)))
        return (arrayout, xaxis, yaxis)

    def tplot(self, flist, label, winsize=20, **kwargs):
        """Plots time vs value for all files in flist with a smoothed version
        on top."""
        colours = ["k", "r", "g", "b", "m", "y", "c"]
        for idx, fname in enumerate(flist):
            color = colours[idx % len(colours)]
            boxlabel = None
            if "rmsd" in label:
                boxlabel = fname[-11:-4]
                if (fname[-11] == "r") and ("ch" not in boxlabel):
                    boxlabel = fname[-9:-4]
                if "ch" in boxlabel:
                    boxlabel = fname[-7:-4]
            data = self.xvg_read(fname)
            data.parse()
            if "Time (ps)" in data.xaxis:
                xlabel = "Time (ns)"
                xvals = data.array[0]*0.001
                ylabel = data.yaxis
                smoothnum = winsize//(xvals.max()/len(xvals))
                if not smoothnum % 2:
                    smoothnum -= 1
                if smoothnum < 2:
                    smoothnum = 3
            self.axlist[-1].plot(xvals, data.array[1], label=boxlabel, markersize=15,
                                 c=color, alpha=0.5, **kwargs)
            if winsize < xvals.max():
                self.axlist[-1].plot(xvals, savitzky_golay(data.array[1],
                                     smoothnum, 1), c=color, **kwargs)
            self.axlist[-1].set_xlabel(xlabel)
            self.axlist[-1].set_ylabel(ylabel)
            if "rmsd" in label:
                self.axlist[-1].legend(loc="lower right", fontsize="x-small")
            if self.end is not None:
                self.axlist[-1].set_xlim(xmin=self.begin/1000, xmax=self.end/1000)
            else:
                self.axlist[-1].set_xlim(xmin=self.begin/1000)

            self.axlist[-1].tick_params(direction='out') 
            self.axlist[-1].yaxis.tick_left()
            self.axlist[-1].xaxis.tick_bottom()
#           ylims = None
#            if label == "rmsd":
#                ylims= [0, 1.0]
#            if ylims is not None:
#                self.axlist[-1].set_ylim(ylims)

    def showfig(self):
        """Shows the fig in a window"""
        self.fig.show()
        raw_input("Press key to close window")

    def savefig(self, filename, local=False):
        """Saves the fig as filename"""
        if local is True:
            if not os.path.isdir("./Mol_Analysis"):
                os.makedirs("./Mol_Analysis")
            self.fig.savefig("./Mol_Analysis/" + filename, format="pdf")
        else:
            self.fig.savefig(self.odir + filename, format="pdf")

    def addplot(self, ncols, nrows, count):
        self.axlist.append(self.fig.add_subplot(ncols, nrows, count))

    def plotall(self, local=True, ncols=3, nrows=2, **kwargs):
        """Plots all files in self.fnames in one figure"""
        self.torun = {"mindist": True, "trjconv": True, "editconf": True,
                      "rmsd": True, "rg": True, "rmsf": True, "pdbout": True}

        self.torun.update(**kwargs)
        nplots = len(self.fnames.keys())
        if nrows==2 and nplots>6:
            sys.exit("More than 6 plot! Please use -large")
        print "Plotting %s plots" % (nplots)
        self.fig.suptitle(self.odir)
        self.fig.add_subplot(ncols, nrows, 1)
        count = 0
        for key, flist in sorted(self.fnames.iteritems()):
            count += 1
            self.addplot(ncols, nrows, count)
            if "rmsf" in key:
                self.resplotstd(flist=flist, label=key, ax=self.axlist[-1])
            elif "group" in key:
                pass
            else:
                self.tplot(flist=flist, label=key)
            self.axlist[-1].set_title(key.replace("_", " "))

        self.fig.tight_layout()
        self.fig.subplots_adjust(top=0.92)
        if local is True:
            localpath = os.getcwd().split("/")
            tags = [x for x in self.odir.split("/") if x not in localpath]
            fname = "plot_"  + "_".join(tags) + ".pdf"
            self.savefig(fname, local=True)
        else:
            fname = "plot.pdf"
            self.savefig(fname)

### Functions


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """

    # From http://scipy.github.io/old-wiki/pages/Cookbook/SavitzkyGolay
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in
               range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window+1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')

### Main


if __name__ == "__main__":
    import argparse
    import ConfigParser
    import tempfile
    import re
    # Classes in Main

    class ListConfigParser(ConfigParser.SafeConfigParser):
        """Modified ConfigParser to read configfile with lists"""
        def listread(self, fname):
            """Reads a config with a list under [Folders] tag"""
            self.dirs = []  # For the list of folders with labels
            tfile = tempfile.TemporaryFile()
            with open(fname, "r") as f:
                for line in f.readlines():
                    if line.strip() == "[Folders]":
                        folderlist = True
                        continue
                    elif "[" in line.strip():
                        folderlist = False
                        tfile.write(line)
                    # Sort out blank lines and comments
                    elif folderlist is True and line.strip() is not "" and not line.startswith("#"):  
                        if len(line.strip().split()) == 1:  # If no label, put in group 0
                            self.dirs.append(tuple(["0", line.strip()]))
                        elif len(line.strip().split()) == 2:
                            self.dirs.append(tuple(line.strip().split()))
                    elif folderlist is False:  # Parse the rest to a tempfile for standard configparser
                        tfile.write(line)
            tfile.seek(0)  # Read tempfile from line 0
            self.readfp(tfile)
            return self.getlist()

        def listwrite(self, f):
            """Writes the config to the file f including the list"""
            if len(self.dirs) > 0:
                f.write("[Folders]\n")
                f.writelines(self.forprint(self.dirs))
                f.write("")
            self.write(f)

        def forprint(self, alist):
            newlist = []
            for n in alist:
                newlist.append(" ".join(n)+"\n")
            return newlist

        def appendfolder(self, dirname):
            """Append to the list"""
            self.dirs.append(dirname+"\n")

        def getlist(self):
            """Gets the list cleaned"""
            return self.dirs

    def multithreader(l_processes, n_threads):  # A default mulithreader
        import Queue
        import threading
        import time
        exitFlag = 0
        if len(l_processes) < n_threads:
            n_threads = len(l_processes)

        class updatedThread(threading.Thread):
            def __init__(self, threadID, name, q):
                threading.Thread.__init__(self)
                self.threadID = threadID
                self.name = name
                self.q = q

            def run(self):
                print "Starting T" + self.name
                process_data(self.name, self.q)
                print "Exiting T" + self.name

        def process_data(threadName, q):
            while not exitFlag:
                queueLock.acquire()
                if not workQueue.empty():
                    data = q.get()
                    queueLock.release()
                    time.sleep(float(threadName)*0.01)
                    print "T%s running:%s" % (threadName, data[1])
                    runandplot(data[1], group=data[0])
                else:
                    queueLock.release()

        threadList = range(n_threads)
        nameList = l_processes
        queueLock = threading.Lock()
        workQueue = Queue.Queue(0)
        threads = []
        threadID = 1
        tnow = time.time()

        # Create new threads
        for tName in threadList:
            thread = updatedThread(threadID, tName, workQueue)
            thread.start()
            threads.append(thread)
            threadID += 1

        # Fill the queue
        queueLock.acquire()
        for word in nameList:
            workQueue.put(word)
        queueLock.release()

        # Wait for queue to empty
        while not workQueue.empty():
            pass

        # Notify threads it's time to exit
        exitFlag = 1

        # Wait for all threads to complete
        for t in threads:
            t.join()
        print "Exiting Main Thread after %d sec" % (time.time()-tnow)

    ### Argparser

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=("This is a script that can "
                                                  "assist in plotting data of "
                                                  "MD simulations with "
                                                  "matplotlib."))
    parser.add_argument("-f", metavar="config.cfg", default="config.cfg",
                        help="Specify the input config file")
    parser.add_argument("-n", action="store_true",
                        help="Dryrun, no gromacs analysis are run")
    parser.add_argument("-nomindist", action="store_false",
                        help="No gromacs mindist analysis are run")
    parser.add_argument("-nopdb", action="store_false",
                        help="No .pdb movie are output")
    parser.add_argument("-v", action="store_true", help="Verbose mode")
    parser.add_argument("-nt", metavar="N", type=int, default=1,
                        help="Number of threads")
    parser.add_argument("-local", action="store_true",
                        help="Locate plots in folder of script")
    parser.add_argument("-nogroup", action="store_true", help="No groupplot")
    parser.add_argument("-begin", type=int, default=0, help="Begin from time X (fs)")
    parser.add_argument("-end", type=int, default=None, help="End at time X (fs)")
    parser.add_argument("-large", action="store_true", help="Can plot 9 plots")
    parser.add_argument("-splitlen", metavar="slen", default=10, help="Length of timewindow for RMSF calculation (ns). If multiple than separated with comma")
    args = parser.parse_args()
    timescales = []
    try:
        timescales = args.splitlen.split(",")
    except AttributeError:
        timescales.append(args.splitlen)
    timescales = [int(i) for i in timescales]
#    print timescales

    # Comment following lines out
    #args.f = "./config_1us.cfg"
    #args.n = True
    #args.local = True
    #args.nogroup = True
    #args.nomindist = False
    #args.nopdb = False
    #args.begin = 10000

    # Start by reading config file
    config = ListConfigParser()
    config.listread(args.f)
    groupfolderlist = config.getlist()
    # Storing defaults here:
    settings = {"odir": "Mol_An", "f": None, "tpr": None,
                "outf": "center_traj.xtc", "outs": "updated.gro",
                "ndxfile": "index.ndx", "ndxparams": None,
                "chains": ["Protein"], "center": "Protein", "skip": None,
                "dt": None}

    # Sorting out the config file, switching out the defaults:
    for key in sorted(settings.keys()):
        try:
            option = config.get("Settings", key)
            if "," in option:
                option = map(lambda s: s.strip(), option.split(",")) 
            elif "none" in option.lower():
                option = None
                continue 
            elif "chains" in key:
                option = [option]
            settings[key] = option
        except ConfigParser.NoOptionError:
            if key == "f" or key == "tpr":
                raise LookupError("Please specify name for "'"%s"'" in "
                                  "the config file." % (key))
            print("No setting for '"'%s'"' found, "
                  "will try default: %s" % (key, settings[key]))
    missing_folder = False 

    # Checking if folders/files exist

    for single_folder in zip(*groupfolderlist)[1]:
        if os.path.isdir(single_folder):
            if args.n is False:
                for test_file in ["f", "tpr"]:
                    if not os.path.isfile(single_folder+"/"+settings[test_file]):
                        print "Cannot find file   %s" % (single_folder+"/"+settings[test_file])
                        missing_folder = True
        else:
            print "Cannot find folder %s" % (single_folder) 
            missing_folder = True

    if missing_folder is True:
        raise IOError("Some folders/files not found")

    ngroups = len(set(zip(*groupfolderlist)[0]))
    nmodels = len(groupfolderlist)
    
    print "Found %s models in %s groups" % (nmodels, ngroups)
    print "\n".join(map(lambda x: " ".join(str(y) for y in x), groupfolderlist)) 

    dummyrunner = gTools_runner(groupfolderlist[0][1], splitlen=timescales, dryrun=args.n, verbose=args.v,
                                group=groupfolderlist[0][0], **settings)
    
    print "Making ndx with '%s'" % dummyrunner.ndxparams
    print "Centering with  '%s'" % dummyrunner.center
    print "Chain[s] set to '%s'" % dummyrunner.fnames.get_chains()

    # List for grabbing data from the multithreader...
    files_to_plot = []
    
    # Function to run in the multithreader
    def runandplot(folder, group=0):
        runner = gTools_runner(wdir=folder, splitlen=timescales, dryrun=args.n, verbose=args.v,
                               group=group, **settings)
        runner.runall(mindist=args.nomindist, pdbout=args.nopdb)
        files_to_plot.append(runner.getfilenames())
    
    # Adding it all to the multithreader
    multithreader(groupfolderlist, args.nt)
    
    for filenames in files_to_plot:
        if args.large:
            plotning = gTools_plotter(filenames, filenames.get_folder(), begin=args.begin, splitlen=timescales, figsize=(11.69, 11.69))
            plotning.plotall(local=args.local, ncols=3, nrows=3)
        else:
            plotning = gTools_plotter(filenames, filenames.get_folder(), begin=args.begin, splitlen=timescales)
            plotning.plotall(local=args.local)
        plt.close()

    if args.nogroup is False:
        groups = sorted(list(set(map(lambda x: x.get_group(), files_to_plot))))  # Get unique groups
        group_sorted = [[y for y in files_to_plot if y.get_group() == x] for x in groups]  # Sort groups

        # Printing rmsf together
        allchains = settings["chains"][:]
        # Remove chains from plot here:
        # allchains.remove("chA")
        for i in timescales:
            groupplotter = gTools_plotter(None, None, begin=args.begin, end=args.end, splitlen=timescales)  # Get new plotter (only for its functions)
            for idx, chain in enumerate(allchains):  # settings["chains"]):
                groupplotter.addplot(len(allchains), 1, idx+1)  # settings["chains"]), 1, idx+1)
                for (label_g, model_g) in zip(groups, group_sorted):
                    groupedfnames = []
                    shift = 0
                    # Align chains here:
                    #if ("p62" in label_g) and (chain == "chB"):
                    #    shift = 709
                    for model in model_g:                    
                        try:
                            groupedfnames.extend(model["rmsf"+str(i)+"_"+chain+"_n_C-alpha"])
                        except IOError:
                            groupedfnames.extend(model["rmsf"+str(i)+"_"+chain+"_&_C-alpha"])
                    label_cleaned = label_g.replace("_", " ")
                    groupplotter.resplotstd(groupedfnames, label_cleaned, std_bar=False, shift=shift)  # , std_fill=False)
                    groupplotter.axlist[-1].legend()
            groupplotter.savefig("grouped_plot_rmsf"+str(i)+".pdf", local=True)            

