.SUFFIXES : .fig .mps .eps .mp .0 .pdf

EPSFIGS = \

#  AmbientFields.eps \
#  Profiles.eps \
#	ScatteredProfile1.eps \
#	map.eps

FIGS = ${EPSFIGS:.eps=.fig}
PDFFIGS = ${EPSFIGS:.eps=.pdf}
PDFSVG = ${SVGFIGS:.svg=.pdf}

all: report.pdf

report.pdf:	report.tex ${PDFFIGS}
	pdflatex report.tex
	cat report.pdf > reportview.pdf
	osascript "$(HOME)/school/utils/PreviewReload.scpt" reportview.pdf

.eps.pdf:
  # Note: epstopdf handles the bounding box.. --nogs means I want to 
  # handle the gs step myself
	epstopdf --nocompress --nogs  $< | gs -q -sDEVICE=pdfwrite -dColorImageFilter=/FlateEncode -dUseFlateCompression=false -dAutoRotatePages=/None -dAutoFilterColorImages=false -sOutputFile="$@" - -c quit
# this needs the -dAutoFilterColorImages=false flag to work properly
#	epstopdf --outfile=$@ --nocompress --debug $<

.svg.pdf:
	/Applications/Inkscape.app/Contents/Resources/bin/inkscape -A $@ $<
