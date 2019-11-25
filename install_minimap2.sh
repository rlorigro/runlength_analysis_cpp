curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -
./minimap2-2.17_x64-linux/minimap2

export PATH=&(pwd)/minimap2-2.17_x64-linux:$PATH
