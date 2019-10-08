FILENAME = readme

$(FILENAME).pdf: $(FILENAME).xtx ./times/*
	xelatex $(FILENAME).xtx
	xelatex $(FILENAME).xtx

clean:
	rm $(FILENAME).aux
	rm $(FILENAME).log
	rm $(FILENAME).out
	rm $(FILENAME).toc
