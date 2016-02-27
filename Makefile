SUBDIRS = IDX_build  IDX_genEx  IDX_parse  IDX_search preProcessY

all: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

.PHONY: $(SUBDIRS)

clean:
	@for i in $(SUBDIRS); do \
	echo "Clearing in $$i..."; \
	(cd $$i; $(MAKE) clean); done

