# cd to OpenFOAM
 cd /lib/openfoam/

# Modify to rwx: sudo chmod -R 777 openfoam2012/
sudo chmod -R 777 openfoam2012/

# Install cmake, g++, wmake
sudo apt install cmake
sudo apt install g++
sudo apt install wmake

# Change to current OpenFOAM-Version
cd openfoam2012/

# Change openfoam2012/etc/config.sh/kahip - KAHIP_VERSION=kahip-3.10
cp /home/ubuntu/home-bs/scripts/kahip etc/config.sh/kahip 

# Make ThirdParty-Dir
rm -r ThirdParty
mkdir ThirdParty
cd ThirdParty

# Download Kahip to openfoam2012/ThirdParty/kahip-2.12(New Dir)
# Currently using release 2.12 (See https://github.com/KaHIP/KaHIP/issues/55)
git clone --branch v2.12 https://github.com/KaHIP/KaHIP.git
mv KaHIP kahip-2.12

# Download the content OpenFOAM ThirdParty-Common to openfoam2012/ThirdParty and extract there
git clone --branch v2012 https://develop.openfoam.com/Development/ThirdParty-common.git
mv ThirdParty-common/* .
rm -r ThirdParty-common

# Compile ./makeKAHIP from openfoam2012/ThirdParty
./makeKAHIP

# Compile OpenFOAM on 16 cores
cd ..
./Allwmake -j 16