sdmreader
=========

Python toolkit for reading the radio astronomy data format known as the Science Data Model (SDM). This data format is the primary format used by the National Radio Astronomy Observatory archive for telescopes such as the VLA and ALMA.

Uses CASA libraries either via casautil (from pwkit) or casapy session. Default is pwkit, but script can be modified for use in casapy session by commenting out definition of "qa" and "me" objects.

Requirements:
	-- Python 2.5 or higher
	-- pwkit 0.3.0 (casautil for casapy-free CASA in python) or CASA
	-- numpy

Example syntax (gets scan info, bdf data, and uvw for one first integration:
> import sdmreader
> (scandict, sourcedict) = sdmreader.read_metadata(sdmfile)
> data = sdmreader.read_bdf(sdmfile, scan)
> (u, v, w) = sdmreader.calc_uvw(sdmfile, scan)

Contributors:
* Casey Law, @caseyjlaw
* Peter Williams, @pkgw
* Adam Deller
* Steve Myers
* Todd Hunter
