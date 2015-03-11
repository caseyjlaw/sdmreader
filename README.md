sdmreader
=========

Python toolkit for reading the radio astronomy data format known as the Science Data Model (SDM). This data format is the primary format used by the National Radio Astronomy Observatory archive for telescopes such as the VLA and ALMA.

Provides basic data/metadata parsing functions. Metadata for scan and source lists comes by parsing xml files (e.g., "Main.xml"), while data comes by parsing BDF (binary data format) as MIME. Simple uvw calculator also provided, though this requires CASA library access (either through casautil (in pwkit) or by running in a casapy session.

Requirements:
---------
* Python 2.5 or higher
* numpy
* Optional (for uvw calculation): pwkit 0.3.0 (casapy-free CASA) or run in casapy

Usage:
------
Example syntax (gets scan info, bdf data, and uvw for one first integration:  
`> import sdmreader`  
`> (scandict, sourcedict) = sdmreader.read_metadata(sdmfile)`  
`> data = sdmreader.read_bdf(sdmfile, scan)`  
`> (u, v, w) = sdmreader.calc_uvw(sdmfile, scan)`  

Contributors:
* Casey Law, @caseyjlaw
* Peter Williams, @pkgw
* Adam Deller
* Steve Myers
* Todd Hunter
