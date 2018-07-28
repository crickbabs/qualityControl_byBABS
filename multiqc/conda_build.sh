#!/bin/sh

module load Anaconda2/5.1.0

PKG=multiqc_plugins
BUILD_DIR=./build
OUT_DIR=./archive

conda build \
	--python 3.6 \
	--croot $BUILD_DIR \
	--output-folder $OUT_DIR \
	$PKG

# clean
rm -rfv $BUILD_DIR
mv -fv $OUT_DIR/linux-64/$PKG-*.tar.bz2 ./
rm -rfv $OUT_DIR

