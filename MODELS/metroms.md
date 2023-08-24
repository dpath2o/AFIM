```{=org}
#+STARTUP: indent , overview
```
# ROMS

Regional Oceanographic Model System (ROMS)

```{=org}
#+NAME: ROMS
```
Compiling a test ROMS on CentOS 7.4 -- *sealion*

Created SSH Key on NECTAR and downloaded to local machine

``` bash
mv ~/Downloads/NECTARprivate_key.pem ~/.ssh/sealion.pem
chmod 0600 ~/.ssh/sealion.pem
```

1.  Created NECTAR CentOS Instance

    NAME
    :   sealion

    ID
    :   6caeb392-4773-4b4b-86a2-ff96cff1ef49

    Description
    :   hear the lion roar

    Project ID
    :   \[d404319b1cca4f10ab265e607a129616\]

    Status
    :   Active

    Locked
    :   False

    Availability Zone
    :   QRIScloud

    Created
    :   15 Sep 2021, 7:11 p.m.

    Age
    :   7 hours, 3 minutes

    1.  Specs

        Flavor Name
        :   r3.xxlarge

        RAM
        :   128GB

        VCPUs
        :   32 VCPU

        Disk
        :   30GB

        1.  IP Addresses

            qld-data
            :   10.255.132.204

            qld
            :   203.101.231.127

    2.  Security Groups

        1.  http

            ALLOW IPv6 to
            :   ::/0

            ALLOW IPv4 443/tcp from
            :   0.0.0.0/0

            ALLOW IPv4 to
            :   0.0.0.0/0

            ALLOW IPv4 80/tcp from
            :   0.0.0.0/0

        2.  icmp

            ALLOW IPv4 to
            :   0.0.0.0/0

            ALLOW IPv4 icmp from
            :   0.0.0.0/0

            ALLOW IPv6 to
            :   ::/0

        3.  ssh

            ALLOW IPv4 22/tcp from
            :   0.0.0.0/0

            ALLOW IPv6 to
            :   ::/0

            ALLOW IPv4 to
            :   0.0.0.0/0

        4.  default

            ALLOW IPv4 from
            :   default

            ALLOW IPv6 to
            :   ::/0

            ALLOW IPv6 from
            :   default

            ALLOW IPv4 to
            :   0.0.0.0/0

        5.  test

            ALLOW IPv6 to
            :   ::/0

            ALLOW IPv4 to
            :   0.0.0.0/0

            ALLOW IPv4 22/tcp from
            :   0.0.0.0/0

    3.  Metadata

        Key Name
        :   Sealion

        Image Name
        :   NeCTAR CentOS 7 x86~64~

        Image ID
        :   0efe8a7e-5a3b-4428-853e-08ab780af042

        1.  Volumes Attached

            Volume
            :   No volumes attached.

2.  Created useful aliases and functions on local machines shell

    ``` bash
    nano .zshrc_me
    ```

    ``` bash
    SEALION_IP=[NECTAR ASSIGNED IP]
    alias sealion="ssh -i ~/.ssh/sealion.pem ec2-user@${SEALION_IP}"
    sealion_scp_from () {
                     scp -i ~/.ssh/sealion.pem ec2-user@${SEALION_IP}:$1 $2
    }
    sealion_scp_to () {
                     scp -i ~/.ssh/sealion.pem $1 ec2-user@${SEALION_IP}:$2
    }
    ```

3.  SSH to *sealion* and changed root password

    ``` bash
    sudo passwd [NECTAR DEFAULT CENTOS USERNAME]
    ```

4.  Download ROMS via git

    1.  Use `git-credential-store`{.verbatim}

        ``` text
        git config --global credential.helper 'store --file ~/.my-credentials'
        ```

    2.  Store ROMS passcode in cryptographic file

        ``` text
        cd ~
        gpg --gen-key
        ```

        Creating an *RSA* of `4096`{.verbatim} with indefinite key
        expiry and then my UTas email address and name

    3.  Created key file using *nano* editor

        ``` text
        nano .netrc
        ```

        and file contents

        ``` text
        machine www.myroms.org
        login [ROMS USNAME]
        password [ROMS PASSWORD]
        protocol https
        ```

    4.  Download the `git-credential-stor`{.verbatim}

        ``` text
        mkdir ~/bin
        curl -o ~/bin/git-credential-netrc https://raw.githubusercontent.com/git/git/master/contrib/credential/netrc/git-credential-netrc
        gpg -e -r [My UTas EMAIL] .netrc
        git config --global credential.helper "netrc -f ~/.netrc.gpg -v"
        ```

    5.  Download ROMS

        ``` text
        mkdir ~/ROMS
        git clone https://[ROMS USERNAME]@www.myroms.org/git/src ~/ROMS
        cd ~/ROMS
        git config filter.id.smudge ".git_filters/id.smudge %f"
        git config filter.id.clean ".git_filters/id.clean %f"
        rm .git/index
        git checkout HEAD -- "$(git rev-parse --show-toplevel)"
        ```

    6.  Download ROMS Documentation

        ``` text
        git clone https://github.com/kshedstrom/roms_manual
        ```

5.  Install critical software

    ``` bash
    sudo yum group install "Development Tools"
    sudo yum install texlive-*
    sudo yum install openmpi
    sudo yum install netcdf-devel
    ```

    The first command critically installs `gcc`{.verbatim},
    `gfortran`{.verbatim}, `cpp`{.verbatim}, `make`{.verbatim} and a
    number of other libraries and software tools that are useful. The
    second installs latex. The third install
    [OpenMPI](https://www.open-mpi.org) , which has not been fully
    tested and may require to be installed from source. The last
    installs `nc-config`{.verbatim} and associated netcdf binaries
    required, however, it does not compile due to the following
    dependency error.

    ``` text
      Error: Package: netcdf-4.3.3.1-5.el7.x86_64 (epel)
              Requires: hdf5(x86-64) = 1.8.12
              Available: hdf5-1.8.12-12.el7.x86_64 (epel)
                  hdf5(x86-64) = 1.8.12-12.el7
              Installing: hdf5-1.8.13-7.el7.x86_64 (centos-openstack)
                  hdf5(x86-64) = 1.8.13-7.el7
    You could try using --skip-broken to work around the problem
    You could try running: rpm -Va --nofiles --nodigest
    ```

6.  Download and build ROMS documentation

    ``` text
    cd ~/ROMS
    git clone https://github.com/kshedstrom/roms_manual
    cd roms_manual
    pdflatex roms_manual.tex
    pdflatex roms_manual.tex
    ```

    Stored locally as

    1.  ROMS Manual

        ```{=org}
        #+NAME: ROMS Manual
        ```
        [ROMS
        Manual](file:///Volumes/SEASONG/PHD/references/modelling/roms_manual.pdf)

7.  ROMS `nc-config`{.verbatim} Dependency

    As per page.6 of the [*ROMS Manual*]{.spurious-link
    target="ROMS Manual"} having this tool is nearly essential, hence,
    I\'ll now attempt to build NetCDF from source

8.  Building NetCDF from source

    Following
    [NETCDF](https://www.unidata.ucar.edu/software/netcdf/documentation/NUG/getting_and_building_netcdf.html)
    build instructions

    ``` bash
    cd ~
    mkdir source_builds
    cd source_builds
    wget https://downloads.unidata.ucar.edu/netcdf-c/4.8.1/src/netcdf-c-4.8.1.tar.gz
    wget https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-fortran-4.5.3.tar.gz
    wget http://www.zlib.net/zlib-1.2.11.tar.gz
    ```

    Downloaded `hdf5`{.verbatim} source code locally and securely copied
    to *sealion*

    1.  building `zlib`{.verbatim}

        ``` bash
        tar xvf zlib.tgz
        cd zlib
          ZDIR=/usr/local
         ./configure --prefix=${ZDIR}
         make check
        sudo make install
        cd ..
        ```

    2.  building `hdf5`{.verbatim}

        ``` bash
        tar xvf hd5.tgz
        cd hdf5
        H5DIR=/usr/local
        ./configure --with-zlib=${ZDIR} --prefix=${H5DIR} --enable-hl
        make test
        sudo make install
        cd ..
        ```

    3.  building of `netcdf-c`{.verbatim}

        ``` bash
        tar xvf netcdf-c.tgz
        cd netcdf-c
        NCDIR=/usr/local
        CPPFLAGS='-I${H5DIR}/include -I${ZDIR}/include' LDFLAGS='-L${H5DIR}/lib -L${ZDIR}/lib' ./configure --prefix=${NCDIR}
        make test
        cd ~
        ```

        This fails with logged failure in [netcdf-c build log
        failure](file:///Volumes/SEASONG/PHD/modelling/log/test-suite.log)

    4.  building `cmake`{.verbatim}

        ``` bash
        wget https://github.com/Kitware/CMake/releases/download/v3.21.2/cmake-3.21.2.tar.gz
        tar xvf cmake-3.21.2
        cd cmake-3.21.2
        ./bootstrap --prefix=/usr/local
        make -j$(nproc)
        make install
        ```

Compiling a test ROMS on Ubuntu 20.04

following [ROMS customise build
script](https://www.myroms.org/wiki/ROMS_UNSW2008#Customize_the_Build_Script)

1.  first need to ensure `EMACS` can `TRAMP` into
2.  copied
    <ssh:ubuntu@131.217.175.244:/home/ubuntu/projects/upwelling/build_roms.sh>
3.  running *upwelling* idealised test with
    `./romsG < roms_upwelling.in > my_upwelling.log`{.verbatim}
4.  Unable to *openMPI* version of ROMS compiled on *seaclub* , opened
    an IT ticket with TPAC helpdesk

ROMS Code [[bookmark]{.smallcaps}]{.tag
tag-name="bookmark"} [[code]{.smallcaps}]{.tag tag-name="code"}

[GitHub Link and Information](https://www.myroms.org/wiki/Git)

[ROMS Forum](https://www.myroms.org/forum/index.php)
[[bookmark]{.smallcaps}]{.tag
tag-name="bookmark"} [[forum]{.smallcaps}]{.tag tag-name="forum"}

[ROMS Working Group](https://www.myroms.org) [[acct]{.smallcaps}]{.tag
tag-name="acct"} [[admin]{.smallcaps}]{.tag tag-name="admin"}

# [METROMS](https://github.com/metno/metroms) [[coupled]{.smallcaps}]{.tag tag-name="coupled"} [[atmo]{.smallcaps}]{.tag tag-name="atmo"} [[ocean]{.smallcaps}]{.tag tag-name="ocean"} [[modelling]{.smallcaps}]{.tag tag-name="modelling"} {#metroms}

Coupled ROMS-CICE through MCT ; For citing of the code, see:
<https://doi.org/10.5281/zenodo.597679>

Setup Curvilinear Grid

1.  Ole Richter\'s Circumpolar Antarctic Grid

    <https://github.com/kuechenrole/antarctic_melting/blob/master/notebooks/preprocessing/make_grid.ipynb>

Curvilinear Definition [[def]{.smallcaps}]{.tag tag-name="def"}

A curvilinear grid or structured grid is a grid with the same
combinatorial structure as a regular grid, in which the cells are
quadrilaterals or (general) cuboids, rather than rectangles or
rectangular cuboids.

Bathymetry [[grid]{.smallcaps}]{.tag
tag-name="grid"} [[bathy]{.smallcaps}]{.tag tag-name="bathy"}

1.  GEBCO2020

    <https://www.gebco.net>

2.  Bedmachine

    <https://sites.uci.edu/morlighem/dataproducts/bedmachine-antarctica/>

3.  Extrapolate and Merge new GA grid

    [this](https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/135334)
    grid, from [Jodi Smith (Geosciences
    Australia)](mailto:Jodie.Smith@ga.gov.au) westward at \~2km using
    [IBCSO](https://www.scar.org/science/ibcso/home/)
