LON1 = 109.43
LON2 = 109.62

DZDT = 0.2

ALL = naga1.east.270
ALL+= naga1.west.90
ALL+= naga2.east.270
ALL+= naga2.west.90
ALL+= naga1.east.CC
ALL+= naga1.west.CC
ALL+= naga2.east.CC
ALL+= naga2.west.CC
ALL+= naga1.east.toLine
ALL+= naga1.west.toLine
ALL+= naga2.east.toLine
ALL+= naga2.west.toLine

TARGETS = $(foreach TPW, $(ALL), $(TPW).$(DZDT).nc)

.phony: all

all: $(TARGETS)

naga1.east.270.$(DZDT).nc: naga1.nc flight.py
	./flight.py --lonEast=$(LON1) --heading=270 --qEastward --dzdt=$(DZDT) $< $@

naga2.east.270.$(DZDT).nc: naga2.nc flight.py
	./flight.py --lonEast=$(LON2) --heading=270 --qEastward --dzdt=$(DZDT) $< $@

naga1.west.90.$(DZDT).nc: naga1.nc flight.py
	./flight.py --lonEast=$(LON1) --heading=90 --qWestward --dzdt=$(DZDT) $< $@

naga2.west.90.$(DZDT).nc: naga2.nc flight.py
	./flight.py --lonEast=$(LON2) --heading=90 --qWestward --dzdt=$(DZDT) $< $@

naga1.east.CC.$(DZDT).nc: naga1.nc flight.py
	./flight.py --lonEast=$(LON1) --current --qEastward --dzdt=$(DZDT) $< $@

naga2.east.CC.$(DZDT).nc: naga2.nc flight.py
	./flight.py --lonEast=$(LON2) --current --qEastward --dzdt=$(DZDT) $< $@

naga1.west.CC.$(DZDT).nc: naga1.nc flight.py
	./flight.py --lonEast=$(LON1) --current --qWestward --dzdt=$(DZDT) $< $@

naga2.west.CC.$(DZDT).nc: naga2.nc flight.py
	./flight.py --lonEast=$(LON2) --current --qWestward --dzdt=$(DZDT) $< $@

naga1.east.toLine.$(DZDT).nc: naga1.nc flight.py
	./flight.py --lonEast=$(LON1) --toLine --qEastward --dzdt=$(DZDT) $< $@

naga2.east.toLine.$(DZDT).nc: naga2.nc flight.py
	./flight.py --lonEast=$(LON2) --toLine --qEastward --dzdt=$(DZDT) $< $@

naga1.west.toLine.$(DZDT).nc: naga1.nc flight.py
	./flight.py --lonEast=$(LON1) --toLine --qWestward --dzdt=$(DZDT) $< $@

naga2.west.toLine.$(DZDT).nc: naga2.nc flight.py
	./flight.py --lonEast=$(LON2) --toLine --qWestward --dzdt=$(DZDT) $< $@
