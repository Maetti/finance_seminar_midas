.PHONY:
	clean
	all

# specify path<
## execution file
ex_path=execution
## source files
source_path=R
## documenatation file
# doc_path=documentation

# code to extract R code and to create Markdown Files
purl=Rscript -e "knitr::purl('$(1)', '$(2)', quiet=TRUE, documentation=0)"
doc_markdown=Rscript -e "rmarkdown::render(input = '$(1)', output_file = '$(2)', output_format='html_document')"

# set the execution files
exfiles:=$(patsubst vignettes/ex-%.Rmd,$(ex_path)/%-ex.R,$(wildcard ex-*.Rmd))

# set the documentation files
sourcefile:=$(patsubst vignettes/fun-%.Rmd,$(source_path)/%-fun.R,$(wildcard vignettes/fun-*.Rmd))

# html output documentation files
#htmlfile:=$(patsubst vignettes/FUN-%.Rmd,$(doc_path)/%-FUN.html,$(wildcard vignettes/FUN-*.Rmd))

# specify all targets
all: $(exfiles) $(sourcefile)

# run execution files
$(ex_path)/%-ex.R: vignettes/ex-%.Rmd
      $(call purl,$^,$@)

# run documentation files
$(source_path)/%-fun.R: vignettes/fun-%.Rmd
	$(call purl,$^,$@)

# run html doc files
# $(doc_path)/%-FUN.html: vignettes/FUN-%.Rmd
#	$(call doc_markdown,$^,$@)

# clean up everything
clean:
	rm -f $(ex_path)*.R
	rm -f $(source_path)*.R