# !/bin/sh

# not yet
test:
	;

# generate doc
doc: clear_doc
	mkdir ./doc || true
	doxygen Doxyfile
	cd ./doc/latex && make
	cp ./doc/latex/refman.pdf ./doc/doc.pdf 

#clear doc folder
clear_doc:
	rm -f -r ./doc/

# run pipeline
run:
	python3 run_pipeline.py

