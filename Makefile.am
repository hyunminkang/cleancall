## Helloworld C++: an example project using Automake

## Place generated object files (.o) into the same directory as their source
## files, in order to avoid collisions when non-recursive make is used.
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}
AUTOMAKE_OPTIONS = foreign
SUBDIRS = src scripts
SAMTOOLS_VERSION = 0.1.18

install-exec-hook:
	echo "Downloading and installing samtools-${SAMTOOLS_VERSION}"
	wget http://sourceforge.net/projects/samtools/files/samtools/0.1.18/samtools-0.1.18.tar.bz2/download -O samtools-0.1.18.tar.bz2
	tar xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2
	make -C samtools-${SAMTOOLS_VERSION}/ samtools
	cp samtools-${SAMTOOLS_VERSION}/samtools $(bindir)
	rm -rf samtools-${SAMTOOLS_VERSION}/
	rm samtools-${SAMTOOLS_VERSION}.tar.bz2
	echo "Finished installing samtools-${SAMTOOLS_VERSION}"

clean-local:
	rm -rf samtools-${SAMTOOLS_VERSION}/
	rm -f samtools-${SAMTOOLS_VERSION}.tar.bz2
