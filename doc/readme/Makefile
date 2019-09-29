FILENAME = readme

$(FILENAME).pdf: $(FILENAME).xtx
	xelatex $(FILENAME).xtx
	xelatex $(FILENAME).xtx

clean:
	rm $(FILENAME).aux
	rm $(FILENAME).log
	rm $(FILENAME).out
	rm $(FILENAME).toc
