curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -
sudo echo -e "\n#add minimap2 to path\n" >> ~/.bashrc
sudo echo export PATH=$(pwd)/minimap2-2.17_x64-linux:$PATH >> ~/.bashrc
source ~/.bashrc
