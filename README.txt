Mehul Ghosal 02/08/2023
	in collaboration with Dr. Rob Jedicke and Dr. Bryce Bolin

DATA
	Each night of data for an asteroid is separated into folders labeled with designation and date of observation. Each folder contains the fits images from CFHT MegaCAM, preprocessed by Elixir.
	There are two sets of images. 
		The first received set can be identified with the 'o' before the chip number in the filename and the '.flt' extension, eg. '2016_GE1_2016_04_04_UTC/1917066o13.flt'. These in general contain the asteroid streak
		The second set was received as a compressed folder of multi-extension FITS files. The source is hidden away on a massive drive you can probably ask for. I uncompressed these and extracted the chips to separate FITS files, with 'on' before the chip number in the filename and the '.fits' extension '2016_GE1_2016_04_04_UTC/1917066o14.fits'


METHOD
	1. Extracting sources: using SExtractor 