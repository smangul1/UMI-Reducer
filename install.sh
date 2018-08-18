cd tools
./install-MiniConda.sh
cd MiniConda/lib
ln -s libncursesw.so.5 libtinfow.so.5



cd ..
cd bin
./conda install -c bioconda pysam
./conda install -c anaconda bzip2 
./conda install -c anaconda numpy
./conda install -c bioconda hisat2
./conda install -c r r
./conda upgrade -c conda-forge readline
./conda install -c r r-ggplot2

MiniConda="$PWD/MiniConda/bin/python"

