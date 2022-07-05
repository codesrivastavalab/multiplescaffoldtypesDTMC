.PHONY: gnu
gnu:
	make -f Makefile.gnu

.PHONY: intel
intel:
	make -f Makefile.intel

.PHONY: clean
clean:
	make -f Makefile.gnu clean
	make -f Makefile.intel clean
