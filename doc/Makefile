# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
BUILDDIR      = _build

APIBUILD      = sphinx-apidoc
APIOPTS       = -o ./
APIDIR        = ../pyphot
APICLEANFILES = pyphot.rst pyphot*.rst 

# Internal variables.
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .
# the i18n builder cannot share the environment and doctrees with the others

# User-friendly check for sphinx-build
ifeq ($(shell which $(SPHINXBUILD) >/dev/null 2>&1; echo $$?), 1)
     $(error The '$(SPHINXBUILD)' command was not found. Make sure you have Sphinx installed, then set the SPHINXBUILD environment variable to point to the full path of the '$(SPHINXBUILD)' executable. Alternatively you can add the directory with the executable to your PATH. If you don\'t have Sphinx installed, grab it from http://sphinx-doc.org/)
endif


.PHONY: help clean

# Put it first so that "make" without argument is like "make help".
help:
help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"
	@echo "  apidoc     to make/refresh the API documentation"
	@echo "  clean      to clean the build directory"

clean:
	rm -rf $(BUILDDIR)/*

apidoc:
	rm -rf $(APICLEANFILES)
	$(APIBUILD) $(APIOPTS) $(APIDIR)

html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html 
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

