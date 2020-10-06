# Python wrappers to reflectivity axitra code using discrete wavenumber
#
# Usefull functions are for unidirectionnal forces:
#   axitra.force_green()
#   axitra.force_conv()

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
    def __init__(self, model, stations, sources, fmax, duration, xl, latlon, freesurface, ikmax, axpath):

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

def force_green(model, stations, sources, duration, fmax, xl=0., latlon=False, freesurface=True, ikmax=10000, axpath='.'):
    '''
    Compute the Green's functions for the set of parameters supplied and for force sources
    The (3 x nfreq x nsources x nstations) Green's functions are stored in a file on disk, the next step is to
    convolve them with the source function(s) using moment_conv() in order to obtain the
    seismograms.

    This function write several input files on disk and call the fortran program "axitra" according
    to the path 'axpath'. The input/output files are of the form "axi_???.suffix" where axi_??? can be obtain
    from the returned Axitra class instance by class.sid

    Input parameters:
    - model = (nlayer x 5) ndarray defining the velocity model
            top_depth/thickness, vp, vs, rho, qp, qs
    - stations = (nstations x 4) ndarray giving the station index and coordinates
    - sources = (nsources x 4) ndarray giving the source index and coordinates
    - duration = time length of seismogram in sec
    - fmax = maximum frequency computed by axitra
    - xl = (optionnal) cylindrical periodicity
    - latlon = (optionnal) coordinates are geographical (True) or cartesian (False)
    - freesurface = (optionnal) compute with a free surface at Z = 0.
    - ikmax = max iteration number for the discrete wavenumber summation

    Return:
    - ap = an instance of Axitra class
    '''

    ## fill in Axitra class with all parameters
    axitra_param = Axitra(model, stations, sources, fmax, duration, xl, latlon, freesurface, ikmax,axpath)

    # save parameters in standard axitra input files
    force_save_files(axitra_param)

    # call green's function computation
    __call_green_fort__(axitra_param)

    return axitra_param

def __call_green_fort__(ap):
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


def force_conv(ap, hist, source_type, t0, t1=0., unit=1, sfunc=None):
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

    import convfPy

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


def l2f(var):
    if var:
        return '.true.'
    else:
        return '.false.'

def force_save_files(ap):
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
