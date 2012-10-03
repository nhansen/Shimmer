INC = -I printCompCounts

subdirs ::
	cd printCompCounts && make all
	mkdir -p blib/bin && cp printCompCounts/printCompCounts blib/bin/
