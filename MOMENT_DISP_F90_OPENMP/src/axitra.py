# Python wrappers to reflectivity axitra code using discrete wavenumber
#
#   
#   axitra.Axitra() to create an instance of Axitra class that contains all parameters
#   methods associated to the class:
#
#   Axitra.read(): fill the class from existing configuration file
#   Axitra.print(): print some member values
#   Axitra.clean(): clean the class and the associated files on disk
#
# Usefull functions are for moment tensor source:
#   axitra.moment.green() to compute Green's functions
#   axitra.moment.conv()  to convolve Green's functions with a source time function
#
# And for unidirectionnal forces:
#   axitra.force.green()
#   axitra.force.conv()

import numpy as np



class Axitra:
    '''
    Class defining the parameters needed to compute Greens's functions using axitra reflectivity code.
    A class instance is created for each run of Green's function computation.

    Contains:
    - id = a unique id per computation
    - sid = a string of form 'axi_id'
    - nlayer = number of layer
    - nstation = number of stations
    - nsource = number of sources
    - model =    (nlayer x 5) ndarray defining the velocity model
    - stations = (nStations x 4) ndarray defining the station coordinates and index
    - sources =  (nstation x 4) ndarray defining the source coordinate and index
    - fmax = maximum frequency computed
    - npt = number of time points
    - duration = duration in sec
    - xl = cylindrical periodicity (default=0. optional, computed if not supplied)
    - latlon = are coordinates geographic (True) or cartesian (False, default)
    - freesurface = with (True, default) or without (False) free surface at Z=0
    - ikmax = max number of iteration (default=10000)
    - axpath = path to axitra binaries (default='.')
    '''
    def __init__(self, model, stations, sources, fmax, duration, xl=0., latlon=False, freesurface=True, ikmax=100000, axpath='.', id=None,aw=None):

        import sys
        import random

        self.model = model
        self.stations = stations
        self.sources = sources
        self.fmax = fmax
        self.duration = duration
        self.latlon = latlon
        self.nfreq = int(fmax*duration)
        self.freesurface = freesurface
        self.ikmax = ikmax
        self.axpath = axpath
        self.aw = aw

        # check model size
        #try:
        #    self.nlayer = model.shape[1]
        #except(IndexError):
        #    print('model array must be of size [6 x nlayer]')
        #    sys.exit(1)
        if model.ndim != 2 or model.shape[1] != 6:
            print('model array must be of size [nlayer x 6]')
            print('top_depth/thickness, vp, vs, rho, qp, qs')
            sys.exit(1)
        else:
            self.nlayer = model.shape[0]

        # check station list
        if stations.ndim != 2 or stations.shape[1] != 4:
            print('stations array must be of size [nstation x 4]')
            print('index, x, y, z')
            sys.exit(1)
        else:
            self.nstation = stations.shape[0]

        # check source list
        if sources.ndim != 2 or sources.shape[1] != 4:
            print('sources array must be of size [nsource x 4]')
            print('index, x, y, z')
            sys.exit(1)
        else:
            self.nsource = sources.shape[0]

        # estimate periodicity if not supplied
        if xl == 0.:
            # vp = model[:,2]
            # vs = model[:,3]
            vmax = max(model[:,2].max(), model[:,3].max())
            self.xl = 1.1 * duration * vmax
        else:
            self.xl = xl

        # compute the number of time point
        # same formula as convm.f90
        xmm = np.log(np.double(self.nfreq)) / np.log(2.)
        mm = int(xmm) + 1
        if (xmm - mm + 1 > 0):
            mm = mm + 1
        self.npt = 2 ** mm

        # id
        if id:
            self.id = id
            self.sid = 'axi_' + str(self.id)
        else:
            self.id = random.randint(1,1000)
            self.sid = 'axi_'+str(self.id)

    def clean(self):
        '''
        Clean all files written to disk relative to the current class instance of Axitra:
        <header>.data, <header>.res, <header>.head , .....
        '''

        import subprocess
        try:
            subprocess.run(["rm "+self.sid+".*"], shell=True)
        except:
            print('ERROR: could not clean axitra file on disk')
            return

    def print(self):
        '''
        print Axitra class parameters
        :return: None
        '''
        print('xl=',self.xl,' duration=',self.duration, ' nfreq=',self.nfreq,' fmax=',self.fmax,' path_to_binary=',self.axpath)
        print('nsource=', self.nsource, ' nstat=', self.nstation,' id=',self.id)

    @classmethod
    def read(self, suffix=None, axpath='.'):
        '''
        read an existing axi_suffix.data file and the associate "source" and "station" files
        and return an Axitra class instance.

        Note: this is done internally by calling a fortran programme <axitra2py> that reads all
        the files and serialized them into a temporary file <axi.datapy> that we read here.
        The axitra2py forytran executable must be compiled before by "make python"

        :param suffix of axi_suffix.data_file, or read axi.Data if no suffix supplied
        :return: an Axitra class instance
        '''

        import subprocess

        if suffix:
            try:
                ier = subprocess.run([axpath + "/axitra2py", suffix])
            except (FileNotFoundError):
                print('ERROR: axitra2py was not found at path: ' + axpath + '/axitra2py')
                return None
        else:
            try:
                ier = subprocess.run([axpath + "/axitra2py"])
            except (FileNotFoundError):
                print('ERROR: axitra2py was not found at path: ' + axpath + '/axitra2py')
                return None

        if ier.returncode == 0:
            # everything ran file
            with open('axi.datapy','r') as f:
                line = f.readline()   # read suffix line
                suffix = line.strip() # remove '\n'
                if suffix == 'no_suffix':
                    suffix = None
                line  = f.readline()  # read namelist parameter
                line  = line.strip()
                col   = line.split()
                nfreq = int(col[0])
                tl    = float(col[1])
                aw    = float(col[2])
                xl    = float(col[3])
                ikmax = int(col[4])
                latlon =  True if col[5] == 'T' else False
                freesurface = True if col[6] == 'T' else False
                sourcefile = f.readline().strip()
                statfile = f.readline().strip()
                # read velocity model
                nc = int(f.readline().strip())
                model = np.empty((nc,6))
                for i in range(0,nc):
                    col = f.readline().strip().split()
                    for r in range(0,6):
                        model[i,r] = float(col[r])
                # read source list
                ns = int(f.readline().strip())
                source = np.empty((ns, 4))
                for i in range(0, ns):
                    col = f.readline().strip().split()
                    for r in range(0, 4):
                        source[i, r] = float(col[r])
                # read station list
                nr = int(f.readline().strip())
                station = np.empty((nr, 4))
                for i in range(0, nr):
                    col = f.readline().strip().split()
                    for r in range(0, 4):
                        station[i, r] = float(col[r])
                # read updated xl and tl
                col = f.readline().strip().split()
                if xl == 0.:
                    xl = float(col[0])
                if tl == 0.:
                    tl = float(col[1])


        else:
            print('ERROR: could not run axitra2py successfully')
            return None

        ap = Axitra(model, station, source, nfreq/tl, tl, xl, latlon, freesurface, ikmax, axpath, id=suffix , aw=aw)
        return ap
# -----------------------------------------------------------------------------------------------
#
#                                    MOMENT CLASS
#
# -----------------------------------------------------------------------------------------------
class moment:
    @classmethod
    def green(self, axitra_param):
        '''
        Compute the Green's functions for the set of parameters supplied and for moment tensor sources
        The (6 x nfreq x nsources x nstations) Green's functions are stored in a file on disk, the next step is to
        convolve them with the source function(s) using moment_conv() in order to obtain the
        seismograms.

        This function write several input files on disk and call the fortran program "axitra" according
        to the path 'axpath'. The input/output files are of the form "axi_???.suffix" where axi_??? can be obtain
        from the returned Axitra class instance by class.sid

        Input parameters:
        - an instance of Axitra class obtained either by class constructor or Axitra.rea()

        Return:
        - ap = the instance of Axitra class
        '''

        # save parameters in standard axitra input files
        self.save_files(axitra_param)

        # call green's function computation
        self.run_green_fortran(axitra_param)

        return axitra_param

    @classmethod
    def conv(self, ap, hist, source_type, t0, t1=0., unit=1, sfunc=None):
        '''
        Compute the convolution of Green's function obtained by a previous call to moment_green()
        by source time functions.

        Input:
        ap = Axitra class instance from moment_green
        hist = source history array (nsource x 6)
               index, moment(Nm), strike, dip, rake, 0., 0., delay
               index, slip, strike, dip, rake, width, height, delay
        source_type =
            0 : Dirac
            1 : Ricker
            2 : smooth step (not causal)
            3 : source time function stored in file <header>.sou
            4 : integral of triangle step
            5 : ramp step
            6 : not used....
            7 : True step (watch high frequencies cutoff!!)
            8 : integral of trapezoid step (~Haskel model)
        add +10 to source_type to get only source time function but no convolution
        (useful for checking)

        t0 = rise time
        t1 = optional time for some sources
        unit = unit on output (1=disp, 2=vel, 3=acc)
        sfunc = source function signal given as numpy array

        Return:
            A time 1D ndarray and three 2D ndarrays of size (nstations x ntimes) for each component
            or None in case of error
        '''
        try:
            import convmPy
        except:
            import sys
            sys.path.append(ap.axpath)
            import convmPy
        import os

        #check if axi_??.data files exists
        try:
            size = os.path.getsize(ap.sid+'.data')
            if size == 0:
                print('Error: No Greens function found, '+ap.sid+'.data is empty')
                return None
        except(FileNotFoundError):
            print('Error: No Greens function found, '+ap.sid+'.data does not exist')
            return None


        # check history file consistency
        if hist.ndim != 2 or hist.shape[1] != 8:
            print('source_history_array must be of size [nsource x 6]')
            print('index, moment(Nm), strike, dip, rake, 0., 0., delay')
            print('or')
            print('index, slip, strike, dip, rake, width, height, delay')
            return None
        else:
            np.savetxt(ap.sid + '.hist', hist, fmt=['%d', '%.3g', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f'], delimiter=" ")

        # allocate memory for results
        sismox = np.zeros((ap.nstation, ap.npt), dtype='float64')
        sismoy = np.zeros((ap.nstation, ap.npt), dtype='float64')
        sismoz = np.zeros((ap.nstation, ap.npt), dtype='float64')

        # check whether source time function is given or not
        if source_type == 3 or source_type == 13:
            try:
                if sfunc.ndim == 1 and sfunc.size == ap.npt:
                    np.savetxt(ap.sid + '.sou', sfunc, fmt='%.3g')
                else:
                    print('wrong array dimension of source time function ('+str(sfunc.shape))
                    print('must be a numpy vector of length '+str(ap.npt))
                    return None
            except:
                print('source_type = [1]3 requires the source time function to be given as sfunc argument')
                return None

        # run convolution
        convmPy.moment_conv(ap.id, source_type, t0, t1, unit, sismox, sismoy, sismoz)

        # create time vector
        dt = ap.duration/ap.npt
        time = np.arange(0, ap.duration, dt,)

        return time, sismox, sismoy, sismoz

    @classmethod
    def save_files(self, ap):
        '''
        Write standard axitra input files: <axi_???.data>, <axi_???.station>, <axi_???.source>
        :param ap:
        :return:
        '''

        from io import StringIO

        # axi.data input file
        with open(ap.sid+'.data','w') as f:
            f.writelines(['&input','\n','nfreq='+str(ap.nfreq),',aw=2.'])
            f.writelines([',tl='+str(ap.duration), ',xl='+str(ap.xl),',ikmax='+str(ap.ikmax),',latlon='+l2f(ap.latlon)])
            f.writelines([',freesurface='+l2f(ap.freesurface),',sourcefile="'+ap.sid+'.source"'])
            f.writelines([',statfile="'+ap.sid+'.stat"','\n','//','\n'])
            # write layer model in a memory file
            fm = StringIO()
            np.savetxt(fm, ap.model, fmt='%.3f', delimiter=" ")
            fm.seek(0)
            fsm = fm.read()
            fm.close()
            # and append its content to axi.data
            f.write(fsm)


        # source coordinates input file
        np.savetxt(ap.sid + '.source', ap.sources, fmt=['%d','%.3f','%.3f','%.3f'], delimiter=" ")

        # station coordinate input file
        np.savetxt(ap.sid + '.stat', ap.stations, fmt=['%d','%.3f','%.3f','%.3f'], delimiter=" ")

    @classmethod
    def run_green_fortran(self,ap):
        '''
        Call the fortran code that compute the Green's function
        ap: axitra class instance that contains required parameters
        :return: ier = return code 0 if ok, >=1 if error
        ier = 1 => executable not found
        ier = 2 => max iteration number reached
        ier = 3 => other error
        '''
        import subprocess
        try:
            ier=subprocess.run([ap.axpath+"/axitra", str(ap.id)])
        except (FileNotFoundError):
            print('ERROR: axitra was not found at path: '+ap.axpath+'/axitra')
            return 1

        if ier.returncode == 0:
            print(ap.axpath+'/axitra ran sucessfully')
            return 0
        elif ier.returncode == 1:
            print('ERROR: '+ap.axpath+'/axitra returned error code 1')
            print('most likely the max iteration number was reached')
            print(' run axitra from outside python in the current directory to see more details')
            return 2
        else:
            print('ERROR: ' + ap.axpath + '/axitra returned a non zero error code')
            print(' run axitra from outside python in the current directory to see more details')
            return 3




# -----------------------------------------------------------------------------------------------
#
#                                    FORCE CLASS
#
# -----------------------------------------------------------------------------------------------
class force:
    @classmethod
    def green(self, axitra_param):
        '''
        Compute the Green's functions for the set of parameters supplied and for force sources
        The (3 x nfreq x nsources x nstations) Green's functions are stored in a file on disk, the next step is to
        convolve them with the source function(s) using moment_conv() in order to obtain the
        seismograms.

        This function write several input files on disk and call the fortran program "axitra" according
        to the path 'axpath'. The input/output files are of the form "axi_???.suffix" where axi_??? can be obtain
        from the returned Axitra class instance by class.sid

        Input parameters:
        - axitra_param = an instance of Axitra class

        Return:
        - axitra_param = a copy of the instance of Axitra class
        '''

        # save parameters in standard axitra input files
        self.save_files(axitra_param)

        # call green's function computation
        self.run_green_fortran(axitra_param)

        return axitra_param

    @classmethod
    def conv(self, ap, hist, source_type, t0, t1=0., unit=1, sfunc=None):
        '''
        Compute the convolution of Green's function obtained by a previous call to moment_green()
        by source time functions.

        Input:
        ap = Axitra class instance from moment_green
        hist = source history array (nsource x 6)
               index, moment(Nm), strike, dip, rake, 0., 0., delay
               index, slip, strike, dip, rake, width, height, delay
        source_type =
            0 : Dirac
            1 : Ricker
            2 : step
            3 : source time function stored in file <header>.sou
            4 : triangle
            5 : ramp
            6 : not used....
            7 : True step (watch high frequencies cutoff!!)
            8 : Trapezoid

        t0 = rise time
        t1 = optional time for some sources
        unit = unit on output (1=disp, 2=vel, 3=acc)
        source_func = source function signal given as numpy array

        Return:
            A time 1D ndarray and three 2D ndarrays of size (nstations x ntimes) for each component
        '''

        try:
            import convfPy
        except:
            import sys
            sys.path.append(ap.axpath)
            import convfPy
        import os

        #check if axi_??.data files exists
        try:
            size = os.path.getsize(ap.sid+'.data')
            if size == 0:
                print('Error: No Greens function found, '+ap.sid+'.data is empty')
                return None
        except(FileNotFoundError):
            print('Error: No Greens function found, '+ap.sid+'.data does not exist')
            return None

        # check history file consistency
        if hist.ndim != 2 or hist.shape[1] != 6:
            print('source_history_array must be of size [nsource x 6]')
            print('index, fx, fy, fz, Amp (N), delay')
            return None
        else:
            np.savetxt(ap.sid + '.hist', hist, fmt=['%d', '%.3g', '%.3f', '%.3f', '%.3f', '%.3f' ], delimiter=" ")

        # allocate memory for results
        sismox = np.zeros((ap.nstation, ap.npt), dtype='float64')
        sismoy = np.zeros((ap.nstation, ap.npt), dtype='float64')
        sismoz = np.zeros((ap.nstation, ap.npt), dtype='float64')

        # check whether source time function is given or not
        if source_type == 3:
            try:
                if sfunc.ndim == 1 and sfunc.size == ap.npt:
                    np.savetxt(ap.sid + '.sou', sfunc, fmt='%.3g')
                else:
                    print('wrong array dimension of source time function ('+str(sfunc.shape))
                    print('must be a numpy vector of length '+str(ap.npt))
                    return
            except:
                print('source_type = 3 requires the source time function to be given as sfunc argument')
                return


        # run convolution
        convfPy.force_conv(ap.id, source_type, t0, t1, unit, sismox, sismoy, sismoz)

        # create time vector
        dt = ap.duration/ap.npt
        time = np.arange(0, ap.duration, dt,)

        return time, sismox, sismoy, sismoz

    @classmethod
    def save_files(self, ap):
        '''
        write standard axitra input files: <axi_???.data>, <station>, <source>
        :param ap:
        :return:
        '''

        from io import StringIO

        # axi.data input file
        with open(ap.sid+'.data','w') as f:
            f.writelines(['&input','\n','nfreq='+str(ap.nfreq),',aw=2.'])
            f.writelines([',tl='+str(ap.duration), ',xl='+str(ap.xl),',ikmax='+str(ap.ikmax),',latlon='+l2f(ap.latlon)])
            f.writelines([',freesurface='+l2f(ap.freesurface),',sourcefile="'+ap.sid+'.source"'])
            f.writelines([',statfile="'+ap.sid+'.stat"','\n','//','\n'])
            # write layer model in a memory file
            fm = StringIO()
            np.savetxt(fm, ap.model, fmt='%.3f', delimiter=" ")
            fm.seek(0)
            fsm = fm.read()
            fm.close()
            # and append its content to axi.data
            f.write(fsm)


        # source coordinates input file
        np.savetxt(ap.sid + '.source', ap.sources, fmt=['%d','%.3f','%.3f','%.3f'], delimiter=" ")

        # station coordinate input file
        np.savetxt(ap.sid + '.stat', ap.stations, fmt=['%d','%.3f','%.3f','%.3f'], delimiter=" ")

    @classmethod
    def run_green_fortran(self, ap):
        '''
        Call the fortran code that compute the Green's function
        ap: axitra class instance that contains required parameters
        :return: ier = return code 0 if ok, >=1 if error
        ier = 1 => executable not found
        ier = 2 => max iteration number reached
        ier = 3 => other error
        '''
        import subprocess
        try:
            ier=subprocess.run([ap.axpath+"/axitra", str(ap.id)])
        except (FileNotFoundError):
            print('ERROR: axitra was not found at path: '+ap.axpath+'/axitra')
            return 1

        if ier.returncode == 0:
            print(ap.axpath+'/axitra ran sucessfully')
            return 0
        elif ier.returncode == 1:
            print('ERROR: '+ap.axpath+'/axitra returned error code 1')
            print('most likely the max iteration number was reached')
            print(' run axitra from outside python in the current directory to see more details')
            return 2
        else:
            print('ERROR: ' + ap.axpath + '/axitra returned a non zero error code')
            print(' run axitra from outside python in the current directory to see more details')
            return 3


def l2f(var):
    if var:
        return '.true.'
    else:
        return '.false.'
