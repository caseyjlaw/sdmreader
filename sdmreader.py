# Copyright 2014 Casey Law <caseyjlaw@gmail.com> and collaborators.
# Licensed under GNU GPL v2.

""" sdmreader -- functions for reading data and metadata from SDM format files

Functions:

read_bdf -- reads data from binary data format and returns numpy array.
calc_uvw -- parses metadata to calculate uvw coordinates for given scan (or time/direction). returns (u,v,w) tuple. Requires CASA libraries.
read_metadata -- parses metadata of SDM file (xml format) to return tuple with two dictionaries (scaninfo, sourceinfo). 
call_qatime -- helper function for calculating times. Requires CASA libraries.

BDFData class does the heavy lifting to parse binary data format and return numpy array. Does not yet parse flag table.

Note: baseline order used in the bdf is a bit unusual and different from what is assumed when working with a measurement set.
Order of uvw and axis=1 of data array Pythonically would be [i*nants+j for j in range(nants) for i in range(j)], so [ (1,2), (1,3), (2,3), (1,4), ...].
"""

# set up CASA tools
try:
    # try casapy-free casa
    import casautil
    me = casautil.tools.measures()
    qa = casautil.tools.quanta()
    print 'Accessing CASA libraries with casautil.'
except ImportError:
    try:
        from casa import quanta as qa
        from casa import measures as me
        print 'Accessing CASA libraries with casapy.'
    except ImportError:
        print 'No CASA library access. '

import numpy as np
import os, glob, mmap
import xml.etree.ElementTree as et
from email.feedparser import FeedParser
from email.message import Message


def read_bdf(sdmpath, scan, nskip=0, readints=0):
    """ Reads given range of integrations from sdm of given scan.
    Uses BDFData object to read.
    readints=0 will read all of bdf (skipping nskip).
    """

    assert os.path.exists(sdmpath)
    scans, sources = read_metadata(sdmpath)
    bdffile = glob.glob(sdmpath + '/ASDMBinary/*' + str(scans[scan]['bdfnum']))[0]

    fp = open(bdffile)
    bdf = BDFData(fp).parse()
    if readints == 0:
        readints = bdf.n_integrations - nskip

    print 'Reading %d ints starting at int %d' % (readints, nskip)
    data = np.empty( (readints, bdf.n_baselines, bdf.n_channels, bdf.n_basebands), dtype='complex64', order='C')
    for i in xrange(readints):
        data[i] = bdf.get_data ('crossData.bin', i+nskip)
#        flag[i] = bdf.get_data ('flags.bin', i+nskip)  # need to get auto+cross parsing right to implement this
    fp.close()
    return data


def calc_uvw(sdmpath, scan=0, datetime=0, radec=()):
    """ Calculates and returns uvw in meters for a given SDM, time, and pointing direction.
    sdmpath is path to sdm directory that includes "Station.xml" file.
    scan is scan number defined by observatory (first scan == 1).
    datetime is time (as string) to calculate uvw (format: '2014/09/03/08:33:04.20')
    radec is (ra,dec) as tuple in units of degrees (format: (180., +45.))
    """

    assert os.path.exists(sdmpath+'/Station.xml')

    # get scan info
    scans, sources = read_metadata(sdmpath)

    # default is to use scan info
    if (datetime == 0) and (len(radec) == 0):
        assert scan != 0   # default scan value not valid

        print 'Calculating uvw for first integration of scan %d of source %s' % (scan, scans[scan]['source'])
        datetime = scans[scan]['start']
        sourcenum = [kk for kk in sources.keys() if sources[kk]['source'] == scans[scan]['source']][0]
        direction = me.direction('J2000', str(np.degrees(sources[sourcenum]['ra']))+'deg', str(np.degrees(sources[sourcenum]['dec']))+'deg')

    # secondary case is when datetime is also given
    elif (datetime != 0) and (len(radec) == 0):
        assert scan != 0   # default scan value not valid
        assert '/' in datetime   # datetime should be in be in yyyy/mm/dd/hh:mm:ss.sss

        print 'Calculating uvw at %s for scan %d of source %s' % (datetime, scan, scans[scan]['source'])
        sourcenum = [kk for kk in sources.keys() if sources[kk]['source'] == scans[scan]['source']][0]
        direction = me.direction('J2000', str(np.degrees(sources[sourcenum]['ra']))+'deg', str(np.degrees(sources[sourcenum]['dec']))+'deg')

    else:
        assert '/' in datetime   # datetime should be in be in yyyy/mm/dd/hh:mm:ss.sss
        assert len(radec) == 2  # radec is (ra,dec) tuple in units of degrees

        print 'Calculating uvw at %s in direction %s' % (datetime, direction)
        ra = radec[0]; dec = radec[1]
        direction = me.direction('J2000', str(ra)+'deg', str(dec)+'deg')

    # define metadata "frame" for uvw calculation
    root = et.parse(sdmpath + '/ExecBlock.xml').getroot()
    telescopename = root.iter('row').next().find('telescopeName').text

    me.doframe(me.observatory(telescopename))
    me.doframe(me.epoch('utc', datetime))
    me.doframe(direction)

    # read antpos
    root = et.parse(sdmpath + '/Station.xml').getroot()
    positions = [rr.find('position').text.split(' ') for rr in root.iter('row') if 'ANTENNA' in rr.find('type').text]
    x = [float(positions[i][2]) for i in range(len(positions))]
    y = [float(positions[i][3]) for i in range(len(positions))]
    z = [float(positions[i][4]) for i in range(len(positions))]
    ants = me.position('itrf', qa.quantity(x, 'm'), qa.quantity(y, 'm'), qa.quantity(z, 'm'))

    # calc bl
    bls = me.asbaseline(ants)
    uvwlist = me.expand(me.touvw(bls)[0])[1]['value']

    # define new bl order to match sdm binary file bl order
    u = np.empty(len(uvwlist)/3); v = np.empty(len(uvwlist)/3); w = np.empty(len(uvwlist)/3)
    nants = len(ants['m0']['value'])
    ord1 = [i*nants+j for i in range(nants) for j in range(i+1,nants)]
    ord2 = [i*nants+j for j in range(nants) for i in range(j)]
    key=[]
    for new in ord2:
        key.append(ord1.index(new))
    for i in range(len(key)):
        u[i] = uvwlist[3*key[i]]
        v[i] = uvwlist[3*key[i]+1]
        w[i] = uvwlist[3*key[i]+2]

    return u, v, w


def read_metadata(sdmfile):
    if (os.path.exists(sdmfile) == False):
        print "Could not find the SDM file = ", sdmfile
        return([],[])
    if (os.path.exists(sdmfile+'/Scan.xml') == False):
        print "Could not find the Scan.xml file.  Are you sure this is an SDM?"
        return([],[])
        
    try:
        from xml.dom import minidom
    except ImportError, e:
        print "failed to load xml.dom.minidom:\n", e
        exit(1)

    # read Scan.xml into dictionary also and make a list
    xmlscans = minidom.parse(sdmfile+'/Scan.xml')
    scandict = {}
    rowlist = xmlscans.getElementsByTagName("row")
    for rownode in rowlist:
        rowfid = rownode.getElementsByTagName("scanNumber")
        fid = int(rowfid[0].childNodes[0].nodeValue)
        # number of subscans
        try:
            # ALMA
            rowsubs = rownode.getElementsByTagName("numSubScan")
            nsubs = int(rowsubs[0].childNodes[0].nodeValue)
        except:
            # EVLA
            rowsubs = rownode.getElementsByTagName("numSubscan")
            nsubs = int(rowsubs[0].childNodes[0].nodeValue)
        # intents
        rownint = rownode.getElementsByTagName("numIntent")
        nint = int(rownint[0].childNodes[0].nodeValue)

        rowintents = rownode.getElementsByTagName("scanIntent")
        sint = str(rowintents[0].childNodes[0].nodeValue)
        sints = sint.split()
        rint = ''
        for r in range(nint):
            intent = sints[2+r]
            if rint=='':
                rint = intent
            else:
                rint += ' '+intent

        # start and end times in mjd ns
        rowstart = rownode.getElementsByTagName("startTime")
        start = int(rowstart[0].childNodes[0].nodeValue)
        startmjd = float(start)*1.0E-9/86400.0
        try:
            t = qa.quantity(startmjd,'d')
            starttime = call_qatime(t,form="ymd",prec=8)
        except:
            pass
        rowend = rownode.getElementsByTagName("endTime")
        end = int(rowend[0].childNodes[0].nodeValue)
        endmjd = float(end)*1.0E-9/86400.0
        try:
            t = qa.quantity(endmjd,'d')
            endtime = call_qatime(t,form="ymd",prec=8)
        except:
            pass

        # source name
        rowsrc = rownode.getElementsByTagName("sourceName")
        if (len(rowsrc) < 1):
            print "Scan %d appears to be corrupt." % (len(scandict)+1)
        else:
            src = str(rowsrc[0].childNodes[0].nodeValue)
            # to find out what all is available,
#            print rownode.getElementsByTagName("*")
            scandict[fid] = {}
            try:
                scandict[fid]['start'] = starttime
            except:
                pass
            scandict[fid]['startmjd'] = startmjd
            try:
                scandict[fid]['end'] = endtime
            except:
                pass
            scandict[fid]['endmjd'] = endmjd
#            print "starttime = ", starttime
#            print "endtime = ", endtime
            try:
                timestr = starttime+'~'+endtime
                scandict[fid]['timerange'] = timestr
            except:
                pass
            scandict[fid]['source'] = src
            scandict[fid]['intent'] = rint
            scandict[fid]['nsubs'] = nsubs
            scandict[fid]['duration'] = endmjd-startmjd
#    print '  Found ',rowlist.length,' scans in Scan.xml'

    xmlmain = minidom.parse(sdmfile+'/Main.xml')
    # iteration option 1
    rowlist = xmlmain.getElementsByTagName("row")
    for rownode in rowlist:
        rowfid = rownode.getElementsByTagName("scanNumber")
        fid = int(rowfid[0].childNodes[0].nodeValue)
        bdfnumstr = rownode.getElementsByTagName('EntityRef')[0].getAttribute('entityId').split('/')[-1]
        try:
            bdfnum = int(bdfnumstr)   # bad scans (killed) seem to have bdfnumstr='X1'
        except ValueError:
            bdfnum = 0
        scandict[fid]['bdfnum'] = bdfnum

    # read Source.xml into dictionary also and make a list
    xmlsources = minidom.parse(sdmfile+'/Source.xml')
    sourcedict = {}
    sourcelist = []
    sourceId = []
    rowlist = xmlsources.getElementsByTagName("row")
    for rownode in rowlist:
        rowfid = rownode.getElementsByTagName("sourceId")
        fid = int(rowfid[0].childNodes[0].nodeValue)

        # source name
        rowsrc = rownode.getElementsByTagName("sourceName")
        src = str(rowsrc[0].childNodes[0].nodeValue)
        try:
            rowsrc = rownode.getElementsByTagName("directionCode")
            directionCode = str(rowsrc[0].childNodes[0].nodeValue)
        except:
            directionCode = ''
        rowsrc = rownode.getElementsByTagName("direction")
        (ra,dec) = rowsrc[0].childNodes[0].nodeValue.split()[2:4]
        ra = float(ra)
        dec = float(dec)
        if (src not in sourcelist):
            sourcelist.append(src)
            sourceId.append(fid)
            sourcedict[fid] = {}
#            sourcedict[fid]['sourceName'] = src
            sourcedict[fid]['source'] = src
            sourcedict[fid]['directionCode'] = directionCode
            sourcedict[fid]['ra'] = ra
            sourcedict[fid]['dec'] = dec
#            print "Loading source %s to index %d" % (src,fid)
        else:
            ai = sourceId[sourcelist.index(src)]
#            print "Source %s is already at index %d = ID:%d" % (src,sourcelist.index(src),ai)
            if (ra != sourcedict[ai]['ra'] or dec != sourcedict[ai]['dec']):
                print "WARNING: Multiple directions found for source %d = %s" % (fid,src)
                ras = (ra - sourcedict[ai]['ra'])*180*3600*math.cos(dec)/math.pi
                decs = (dec - sourcedict[ai]['dec'])*180*3600/math.pi
                print "The difference is (%f,%f) arcseconds." % (ras,decs)
#    for src in range(len(sourcedict)):
#        print "%s direction = %f, %f" % (sourcedict[src]['sourceName'],
#                                         sourcedict[src]['ra'],
#                                         sourcedict[src]['dec'])
        
    # return the dictionary for later use
    return [scandict, sourcedict]


def call_qatime(arg, form='', prec=0):
    """
    This is a wrapper for qa.time(), which in casa 4.0 returns a list of 
    strings instead of just a scalar string.  In this case, return the first 
    value in the list.
    - Todd Hunter
    """

    result = qa.time(arg, form=form, prec=prec)
    if (type(result) == list or type(result)==np.ndarray):
        return(result[0])
    else:
        return(result)


"""
Class for reading Binary Data File, the raw format for SDM data.

Abuse the email.feedparser class to read the structure of a BDF binary file.
We remember where the binary blobs are so we can mmap them into numpy arrays
later.

f = open ('/path/to/bdf_file')
bdf = BDFData (f).parse ()
for i in xrange (bdf.n_integrations):
    c = bdf.get_data ('crossData.bin', i)
    nbl, nchan, npol = c.shape
    ....

The following arrays can be retrieved with the get_data() function:

crossData.bin: the basic cross-correlation data
  dtype: complex64
  shape: (nbaselines, nchans, npol)

autoData.bin: the basic auto-correlation data
  dtype: complex64
  shape: (nantennas, nchans, nfeeds=2)

flags.bin: flags on the auto and cross correlations
  dtype: uint32
  shape: (nbaselines + nantennas, nchans, npol)
  FIXME: I don't know how the baseline/antenna dimension should be indexed!
         Consult the BDF spec. Also, this will break if the BDF does not
         contain auto+cross data.

Shortcomings:

We hardcode the array axis orderings and which axis have non-unity size. This
could theoretically all change under us.

The BDF spec makes it sound like the different binary blobs are allowed to
have differing sizes -- the "size" attribute in the header is a maximum. If
this ever happens in practice, we're kind of dicked -- I don't see how we
can guess the datachunk size without just reading it in. Let's hope that
never happens.

BDF is little-endian as are x86 processors, so we ignore endianness issues.
"""

__all__ = ['BDFData']

_datatypes = {
    'autoData.bin': np.complex64,
    'crossData.bin': np.complex64,
    'flags.bin': np.uint32,
}

nanttag = 'numAntenna'
basebandtag = 'baseband'

class BDFData (object):
    def __init__ (self, fp):
        """fp is an open, seekable filestream."""
        self.fp = fp
        self.mmdata = mmap.mmap (fp.fileno (), 0, mmap.MAP_PRIVATE, mmap.PROT_READ)


    def parse (self):
        """Parse the BDF mime structure and record the locations of the binary
        blobs. Sets up various data fields in the BDFData object."""

        feedparser = FeedParser (Message)
        binarychunks = {}
        sizeinfo = None
        self.fp.seek (0, 0)

        while True:
            data = self.fp.readline ()
            if not data:
                break

            feedparser.feed (data)

            skip = (data == '\n' and
                    len (feedparser._msgstack) == 3 and
                    feedparser._msgstack[-1].get_content_type () in ('application/octet-stream',
                                                                     'binary/octet-stream'))

            if skip:
                # We just finished reading the headers for a huge binary blob.
                # Time to remember where the data chunk is and pretend it doesn't
                # exist.
                msg = feedparser._msgstack[-1]
                ident = msg['Content-Location']
                assert ident.endswith ('.bin'), 'confusion #1 in hacky MIME parsing!'
                binarychunks[ident] = self.fp.tell ()
                if sizeinfo is None:
                    headxml, sizeinfo, tagpfx = _extract_size_info (feedparser)
                kind = ident.split ('/')[-1]
                assert kind in sizeinfo, 'no size info for binary chunk kind %s in MIME!' % kind
                self.fp.seek (sizeinfo[kind] + 1, 1) # skip ahead by data chunk size
                sample = self.fp.read (16)
                assert sample.startswith ('--MIME'), 'crap, unexpected chunk size in MIME parsing: %r' % sample
                self.fp.seek (-16, 1) # go back

        if headxml is None:
            raise RuntimeError ('never found any binary data')

        self.mimemsg = feedparser.close ()
        self.headxml = headxml
        self.sizeinfo = sizeinfo
        self.binarychunks = binarychunks

        # Compute some miscellaneous parameters that we'll need.
        self.n_integrations = len (self.mimemsg.get_payload ()) - 1
        self.n_antennas = int (headxml.find (tagpfx + nanttag).text)
        self.n_baselines = (self.n_antennas * (self.n_antennas - 1)) // 2

        ds = headxml.find (tagpfx + dstag)
        nbb = 0
        nspw = 0
        nchan = 0
        crosspolstr = None

        for bb in ds.findall (tagpfx + basebandtag):
            nbb += 1

            for spw in bb.getchildren ():
                nspw += 1
                nchan += int (spw.get ('numSpectralPoint'))

                if crosspolstr is None:
                    crosspolstr = spw.get ('crossPolProducts')
                elif spw.get ('crossPolProducts') != crosspolstr:
                    raise Exception ('can only handle spectral windows with identical cross pol products')

        self.n_basebands = nbb
        self.n_spws = nspw
        self.n_channels = nchan
        self.crosspols = crosspolstr.split ()
        return self # convenience


    def get_data (self, datakind, integnum):
        """Given an integration number (0 <= integnum < self.n_integrations) and a
        data kind ('crossData.bin', 'autoData.bin'), memory-map the corresponding data
        and return a wrapping numpy array."""

        if integnum < 0 or integnum >= self.n_integrations:
            raise ValueError ('illegal integration number %d' % integnum)

        size = self.sizeinfo.get (datakind)
        if size is None:
            raise ValueError ('unrecognized data kind "%s"' % datakind)

        dtype = _datatypes[datakind]
        key = '/%d/%s' % (integnum + 1, datakind) # numbers are 1-based here

        for ident, offset in self.binarychunks.iteritems ():
            if ident.endswith (key):
                break
        else:
            # Gets executed if we don't break out of the loop.
            raise ValueError ('can\'t find integration #%d of kind %s in BDF'
                              % (integnum, datakind))

        dslice = self.mmdata[offset:offset+size]
        data = np.fromstring (dslice, dtype=dtype)

        if datakind == 'crossData.bin':
            data = data.reshape ((self.n_baselines, self.n_channels, len (self.crosspols)))
        elif datakind == 'autoData.bin':
            data = data.reshape ((self.n_antennas, self.n_channels, 2))
        elif datakind == 'flags.bin':
            data = data.reshape ((self.n_baselines + self.n_antennas, self.n_channels,
                                  len (self.crosspols)))

        return data


tagprefixes = ['{http://Alma/XASDM/sdmbin}', '']
dstag = 'dataStruct'
cdtag = 'crossData'
adtag = 'autoData'
fgtag = 'flags'

def _extract_size_info (feedparser):
    # This parses the XML of the header section

    text = feedparser._msgstack[0].get_payload ()[0].get_payload ()
    headxml = et.fromstring (text)

    # The XML may or may not have an xmlns attribute which manifests itself
    # as a prefix to the tags we need to use.

    sizeinfo = {}

    for tp in tagprefixes:
        dselem = headxml.find (tp + dstag)
        if dselem is not None:
            break
    else:
        raise RuntimeError ('cannot find dataStruct item in any known XML namespace')

    # Now pull out the dataStruct bit and chunk size info. We return
    # sizes in bytes, not data elements.

    e = dselem.find (tp + cdtag)
    if e is not None:
        sizeinfo['crossData.bin'] = 4 * int (e.attrib['size'])

    e = dselem.find (tp + adtag)
    if e is not None:
        sizeinfo['autoData.bin'] = 4 * int (e.attrib['size'])

    e = dselem.find (tp + fgtag)
    if e is not None:
        sizeinfo['flags.bin'] = 4 * int (e.attrib['size'])

    # ... could fill in more if needed ...

    return headxml, sizeinfo, tp

